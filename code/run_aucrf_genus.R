source('code/utilities.R')
get_dependencies(c('AUCRF', 'pROC'))

get_genus_shared <- function(dataset){

	shared_file <- paste0('data/', dataset, '/', dataset, '.0.03.subsample.shared')
	shared_otus <- read.table(file=shared_file, header=T, stringsAsFactors=F, row.names=2)[,-c(1,2)]
	is_present <- apply(shared_otus, 2, sum) > 0
	shared <- shared_otus[,is_present]

	tax_file <- paste0('data/', dataset, '/', dataset, '.taxonomy')
	tax_data <- read.table(file=tax_file, header=T, stringsAsFactors=F)
	taxonomy <- tax_data$Taxonomy
	taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
	taxonomy <- gsub('"', '', taxonomy)
	taxonomy <- gsub(';', '.', taxonomy)
	names(taxonomy) <- tax_data$OTU
	taxonomy <- taxonomy[colnames(shared)]

	unique_taxa <- levels(as.factor(taxonomy))

	shared_genus <- NULL

	for(ut in unique_taxa){
		otus <- names(taxonomy[taxonomy %in% ut])
		sub_shared <- shared_otus[,colnames(shared_otus) %in% otus]

		if(is.null(dim(sub_shared))){
			shared_genus <- cbind(shared_genus, sub_shared)
		} else {
			genus_count <- apply(sub_shared, 1, sum)
			shared_genus <- cbind(shared_genus, genus_count)
		}
	}
	colnames(shared_genus) <- unique_taxa
	rownames(shared_genus) <- rownames(shared)
	return(shared_genus)
}

get_metrics <- function(predicted, reference){

	contingency_table <- table(predicted = predicted, reference=reference)

	sensitivity <- contingency_table[2,2] / sum(contingency_table[,2])
	specificity <- contingency_table[1,1] / sum(contingency_table[,1])
	accuracy <- (contingency_table[1,1] + contingency_table[2,2]) / sum(contingency_table)
	posPredValue <- contingency_table[2,2] / sum(contingency_table[2,])
	negPredValue <- contingency_table[1,1] / sum(contingency_table[1,])

	return(c(sensitivity = sensitivity, specificity = specificity,
		 			accuracy = accuracy, posPredValue = posPredValue,
					negPredValue = negPredValue))
}

run <- function(datasets){

	datasets <- unlist(strsplit(datasets, split=" "))

#	datasets <- c('baxter', 'escobar', 'turnbaugh', 'wu', 'ross', 'hmp', 'zupancic')

	test_data <- list()
	model <- list()
	roc <- list()
	model_summary <- NULL
	roc_summary <- NULL


	for(d in datasets){
		print(d)
		set.seed(1976)

		shared <- get_genus_shared(d)
		n_seqs <- sum(shared[1,])

		metadata_file <- paste0('data/', d, '/', d, '.metadata')
		metadata <- read.table(file=metadata_file, header=T)
		metadata <- metadata[metadata$sample %in% rownames(shared),]

		stopifnot(rownames(shared) == metadata$sample)

		keep_otus <- apply(shared > 0, 2, sum) > (0.1 * nrow(shared))

		rel_abund_keep <- shared[,keep_otus] / n_seqs
		test_data[[d]] <- data.frame(obese=as.factor(as.numeric(metadata$obese)), rel_abund_keep)

		model[[d]] <- AUCRF(obese ~ ., data=test_data[[d]], ntree=1000, nodesize=20)
		k_opt <- model[[d]]$Kopt
		otus <- paste(model[[d]]$Xopt, collapse=',')

		probabilities <- predict(model[[d]]$RFopt, type='prob')[, 2]
		roc[[d]] <- roc(metadata$obese~probabilities)

		model[[d]]$auc <- ci(roc[[d]])

		opt_index <- which.max(roc[[d]]$sensitivities + roc[[d]]$specificities)
		model[[d]]$opt_sensitivity <- roc[[d]]$sensitivities[opt_index]
		model[[d]]$opt_specificity <- roc[[d]]$specificities[opt_index]
		model[[d]]$opt_threshold <- roc[[d]]$thresholds[opt_index]

		predicted <- factor(probabilities >= model[[d]]$opt_threshold, levels=c("FALSE", "TRUE"))
		metrics <- get_metrics(predicted, metadata$obese)

		model[[d]]$opt_accuracy <- metrics["accuracy"]
		model[[d]]$opt_posPredValue <- metrics["posPredValue"]
		model[[d]]$opt_negPredValue <- metrics["negPredValue"]

		roc_summary <- rbind(roc_summary, cbind(d, roc[[d]]$sensitivities, roc[[d]]$specificities))
		model_summary <- rbind(model_summary, c(d, model[[d]]$auc,k_opt, otus,
					 							model[[d]]$opt_sensitivity, model[[d]]$opt_specificity,
												model[[d]]$opt_threshold))
	}

	colnames(model_summary) <- c("dataset", "auc_lci", "auc", "auc_hci", "k_opt", "otus", "sensitivity", "specificity", "threshold")
	write.table(model_summary, file="data/process/random_forest.genus.summary", quote=F, sep='\t', row.names=F)

	colnames(roc_summary) <- c("dataset", "sensitivity", "specificity")
	write.table(roc_summary, file="data/process/random_forest.genus.roc_data", quote=F, sep='\t', row.names=F)


	# need to do the round-robin where we run predict wiht each model against all
	# of the other datasets
	testing_summary <- NULL

	for(train in datasets){

		train_otus <- model[[train]]$Xopt
		forest <- model[[train]]$RFopt
		threshold <- model[[train]]$opt_threshold

		for(test in datasets){

			if(test != train){
				otus_missing <- !train_otus %in% colnames(test_data[[test]])
				if(sum(otus_missing) > 0){
					test_data[[test]][,train_otus[otus_missing]] <- 0
				}

				p <- predict(forest, test_data[[test]], type='prob')
				predicted <- factor(p[,2] >= threshold, levels=c('FALSE', 'TRUE'))

				roc_test <- roc(test_data[[test]]$obese~p[,2])

				metrics <- get_metrics(predicted, test_data[[test]]$obese)

				train_test_p <- roc.test(roc_test, roc[[train]])$p.value

				result <- c(train, test, metrics["sensitivity"], metrics["specificity"],
									metrics["accuracy"], metrics["posPredValue"],
									metrics["negPredValue"], ci(roc_test), train_test_p)
				testing_summary <- rbind(testing_summary, result)
			} else {
				result <- c(train, test, model[[train]]$opt_sensitivity,
																	model[[train]]$opt_specificity,
																	model[[train]]$opt_accuracy,
																	model[[train]]$opt_posPredValue,
																	model[[train]]$opt_negPredValue,
																	model[[train]]$auc, NA)
				testing_summary <- rbind(testing_summary, result)

			}
		}
	}

	colnames(testing_summary) <- c("train", "test", "sensitivity", "specificity",
																	"accuracy", "posPredValue", "negPredValue",
																	"auc_lci", "auc", "auc_hci", "p_value")
	rownames(testing_summary) <- NULL
	write.table(testing_summary, file="data/process/random_forest.genus.train_test", quote=F, sep='\t', row.names=F)

}

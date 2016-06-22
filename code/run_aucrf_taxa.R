source('code/utilities.R')
source('code/cross_validate.R')
get_dependencies(c('AUCRF', 'pROC'))


tax_name <- c("kingdom", "phylum", "class", "order", "family", "genus")


get_tax_substring <- function(tax, tax_level){
	substring <- unlist(strsplit(tax, ";"))[1:tax_level]
	paste(substring, collapse='.')
}

get_tax_name <- function(tax_file, tax_level){


	tax_data <- read.table(file=tax_file, header=T, stringsAsFactors=F)
	taxonomy <- tax_data$Taxonomy
	taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
	taxonomy <- gsub('"', '', taxonomy)

	tax_substring <- sapply(taxonomy, get_tax_substring, tax_level)

	names(tax_substring) <- tax_data$OTU

	tax_substring
}

get_tax_level_shared <- function(dataset, tax_level){

	shared_file <- paste0('data/', dataset, '/', dataset, '.0.03.subsample.shared')
	shared_otus <- read.table(file=shared_file, header=T, stringsAsFactors=F, row.names=2)[,-c(1,2)]
	is_present <- apply(shared_otus, 2, sum) > 0
	shared <- shared_otus[,is_present]

	tax_file <- paste0('data/', dataset, '/', dataset, '.taxonomy')
	taxonomy <- get_tax_name(tax_file, tax_level)
	taxonomy <- taxonomy[colnames(shared)]
	unique_taxa <- levels(as.factor(taxonomy))

	shared_tax_level <- NULL

	for(ut in unique_taxa){
		otus <- names(taxonomy[taxonomy %in% ut])
		sub_shared <- shared_otus[,colnames(shared_otus) %in% otus]

		if(is.null(dim(sub_shared))){
			shared_tax_level <- cbind(shared_tax_level, sub_shared)
		} else {
			tax_level_count <- apply(sub_shared, 1, sum)
			shared_tax_level <- cbind(shared_tax_level, tax_level_count)
		}
	}
	colnames(shared_tax_level) <- unique_taxa
	rownames(shared_tax_level) <- rownames(shared)
	return(shared_tax_level)
}

run <- function(datasets, tax_level){

	datasets <- unlist(strsplit(datasets, split=" "))

#	datasets <- c('baxter', 'escobar', 'turnbaugh', 'wu', 'ross', 'hmp', 'zupancic', 'zeevi', 'schubert', 'goodrich')

	test_data <- list()
	model <- list()
	roc <- list()
	model_summary <- NULL
	roc_summary <- NULL


	for(d in datasets){
		print(d)
		set.seed(1976)

		shared <- get_tax_level_shared(d, tax_level)
		n_seqs <- sum(shared[1,])

		metadata_file <- paste0('data/', d, '/', d, '.metadata')
		metadata <- read.table(file=metadata_file, header=T)
		metadata <- metadata[metadata$sample %in% rownames(shared),]

		stopifnot(rownames(shared) == metadata$sample)

		na_obesity <- is.na(metadata$obese)
		metadata <- metadata[!na_obesity,]
		shared <- shared[!na_obesity,]

		stopifnot(rownames(shared) == metadata$sample)

		keep_otus <- apply(shared > 0, 2, sum) > (0.1 * nrow(shared))

		rel_abund_keep <- shared[,keep_otus] / n_seqs
		test_data[[d]] <- data.frame(obese=as.factor(as.numeric(metadata$obese)), rel_abund_keep)

		model[[d]] <- AUCRF(obese ~ ., data=test_data[[d]], ntree=500, nodesize=10)
		k_opt <- model[[d]]$Kopt
		otus <- paste(model[[d]]$Xopt, collapse=',')

		probabilities <- predict(model[[d]]$RFopt, type='prob')[, 2]
		roc[[d]] <- roc(metadata$obese~probabilities)

		model[[d]]$auc <- ci(roc[[d]])

		model[[d]]$auc_cv <- cross_validate(model[[d]])


		opt_index <- which.max(roc[[d]]$sensitivities + roc[[d]]$specificities)
		model[[d]]$opt_sensitivity <- roc[[d]]$sensitivities[opt_index]
		model[[d]]$opt_specificity <- roc[[d]]$specificities[opt_index]
		model[[d]]$opt_threshold <- roc[[d]]$thresholds[opt_index]

		predicted <- factor(probabilities >= model[[d]]$opt_threshold, levels=c("FALSE", "TRUE"))

		accuracy <- ci.coords(roc[[d]], x=model[[d]]$opt_threshold, input = "threshold", ret="acc", progress="none")

		model[[d]]$opt_accuracy <- accuracy[2]
		model[[d]]$opt_accuracy_lci <- accuracy[1]
		model[[d]]$opt_accuracy_uci <- accuracy[3]

		max_rsq <- NA

		if(d != 'turnbaugh'){
			reg_data <- data.frame(bmi=metadata$bmi, rel_abund_keep)
			model_reg <- randomForest(bmi ~ ., data=reg_data, ntree=500, nodesize=10)

			o <- order(model_reg$importance, decreasing=T)

			limit <- ifelse(length(o) <= 30,length(o),30)
			rsq <- rep(NA, limit)

			for(i in 2:limit){
				reg_data <- data.frame(bmi=metadata$bmi, rel_abund_keep[,o[1:i]])
				rsq[i] <- randomForest(bmi ~ ., data=reg_data, ntree=500, nodesize=10)$rsq[500]
			}

			max_rsq <- max(rsq, na.rm=T)
		}

		roc_summary <- rbind(roc_summary, data.frame(dataset=d, sensitivity = roc[[d]]$sensitivities, specificity = roc[[d]]$specificities))

		model_summary <- rbind(model_summary, c(d, model[[d]]$auc,
			 									model[[d]]$auc_cv["cv_est"], k_opt, otus,
												model[[d]]$opt_sensitivity,
												model[[d]]$opt_specificity,
												model[[d]]$opt_threshold, max_rsq))
	}

	colnames(model_summary) <- c("dataset", "auc_lci", "auc", "auc_hci", "auc_cv", "k_opt", "otus", "sensitivity", "specificity", "threshold", "regression_auc")

	summary_file <- paste0("data/process/random_forest.", tax_name[tax_level], ".summary")
	write.table(model_summary, file=summary_file, quote=F, sep='\t', row.names=F)

	roc_file <- paste0("data/process/random_forest.", tax_name[tax_level], ".roc_data")
	write.table(roc_summary, file=roc_file, quote=F, sep='\t', row.names=F)


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

				accuracy <- ci.coords(roc_test, x=model[[d]]$opt_threshold, input = "threshold", ret="acc", progress="none")

				train_test_p <- roc.test(roc_test, roc[[train]])$p.value

				result <- c(train, test,
									accuracy[2], accuracy[1],
									accuracy[3], train_test_p)
				testing_summary <- rbind(testing_summary, result)
			} else {
				result <- c(train, test, model[[train]]$opt_accuracy,
																	model[[train]]$opt_accuracy_lci,
																	model[[train]]$opt_accuracy_uci, NA)
				testing_summary <- rbind(testing_summary, result)

			}
		}
	}

	colnames(testing_summary) <- c("train", "test",
																	"accuracy", "accuracy_lci", "accuracy_uci",
																	"p_value")
	rownames(testing_summary) <- NULL

	testing_file <- paste0("data/process/random_forest.", tax_name[tax_level], ".train_test")
	write.table(testing_summary, file=testing_file, quote=F, sep='\t', row.names=F)

}

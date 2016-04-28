source('code/utilities.R')
get_dependencies(c('AUCRF', 'pROC'))

run <- function(datasets){

	datasets <- unlist(strsplit(datasets, split=" "))

#	datasets <- c('baxter', 'turnbaugh', 'wu', 'escobar', 'ross', 'hmp', 'zupancic')

	model_summary <- NULL
	roc_summary <- NULL

	for(d in datasets){
		print(d)
		set.seed(1976)

		shared_file <- paste0('data/', d, '/', d, '.0.03.subsample.shared')
		shared <- read.table(file=shared_file, header=T, stringsAsFactors=F, row.names=2)[,-c(1,2)]
		n_seqs <- sum(shared[1,])

		metadata_file <- paste0('data/', d, '/', d, '.metadata')
		metadata <- read.table(file=metadata_file, header=T)
		metadata <- metadata[metadata$sample %in% rownames(shared),]

		stopifnot(rownames(shared) == metadata$sample)

		keep_otus <- apply(shared > 0, 2, sum) > (0.1 * nrow(shared))

		rel_abund_keep <- shared[,keep_otus] / n_seqs
		test_data <- data.frame(obese=as.factor(as.numeric(metadata$obese)), rel_abund_keep)

		model <- AUCRF(obese ~ ., data=test_data, ntree=500, nodesize=10)
		k_opt <- model$Kopt
		otus <- paste(model$Xopt, collapse=',')

		probabilities <- predict(model$RFopt, type='prob')[, 2]
		roc <- roc(metadata$obese~probabilities)

		auc <- ci(roc)

		auc_cv <- AUCRFcv(model, nCV=10)$cvAUC

		roc_summary <- rbind(roc_summary, cbind(d, roc$sensitivities, roc$specificities))
		model_summary <- rbind(model_summary, c(d, auc, auc_cv, k_opt, otus))

	}

	colnames(model_summary) <- c("dataset", "auc_lci", "auc", "auc_hci", "auc_cv", "k_opt",  "otus")
	write.table(model_summary, file="data/process/random_forest.otu.summary", quote=F, sep='\t', row.names=F)

	colnames(roc_summary) <- c("dataset", "sensitivity", "specificity")
	write.table(roc_summary, file="data/process/random_forest.otu.roc_data", quote=F, sep='\t', row.names=F)

}

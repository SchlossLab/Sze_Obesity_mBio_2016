source('code/utilities.R')

get_dependencies("epiR")


analyze <- function(alpha, is_obese){

	threshold <- median(alpha)
	low_high <- factor(ifelse(alpha <= threshold, "low", "high"), levels=c("low", "high"))
	obese_factor <- factor(is_obese, levels=c("TRUE", "FALSE"))

	contingency <- table(low_high, obese_factor)
	test <- epi.2by2(contingency, method="cohort.count")

	rr <- test$massoc$RR.strata.score
	p.value <- test$massoc$chisq.strata$p.value

	counts <- as.vector(contingency)
	names(counts) <- c('low_obese', 'high_obese', 'low_nonobese', 'high_obese')

	return(unlist(c(rr, p.value=p.value, counts)))
}



run <- function(datasets){

	datasets <- unlist(strsplit(datasets, split=" "))

	summary_data <- NULL

	for(d in datasets){
		alpha_file <- paste0('data/', d, '/', d, '.groups.ave-std.summary')
		alpha <- read.table(file=alpha_file, header=T, stringsAsFactors=F)
		alpha <- alpha[alpha$method == 'ave',]

		metadata_file <- paste0('data/', d, '/', d, '.metadata')
		metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=F)
		metadata <- metadata[metadata$sample %in% alpha$group,]

		stopifnot(alpha$group == metadata$sample)

		bf_relabund <- get_bacteroides_firmicutes(d)

		shannon <- analyze(alpha$shannon, metadata$obese)
		sobs <- analyze(alpha$sobs, metadata$obese)
		shannoneven <- analyze(alpha$shannoneven, metadata$obese)
		bacteroidetes <- analyze(bf_relabund[,"b"], metadata$obese)
		firmicutes <- analyze(bf_relabund[,"f"], metadata$obese)
		bf_ratio <- analyze(bf_relabund[,"bf"], metadata$obese)

		test <- rbind(shannon, sobs, shannoneven, bacteroidetes, firmicutes, bf_ratio)
		study_data <- cbind(dataset = d, metric = rownames(test), test)

		summary_data <- rbind(summary_data, study_data)
	}

	write.table(summary_data, file="data/process/relative_risk.summary", quote=F, sep='\t', row.names=F)

}

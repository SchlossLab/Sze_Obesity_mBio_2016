source('code/utilities.R')

get_dependencies(c("epiR", "metafor"))


analyze_indiv_study <- function(alpha, is_obese){

	threshold <- median(alpha)
	low_high <- factor(ifelse(alpha <= threshold, "low", "high"), levels=c("low", "high"))
	obese_factor <- factor(is_obese, levels=c("TRUE", "FALSE"))

	contingency <- table(low_high, obese_factor)
	test <- epi.2by2(contingency, method="cohort.count")

	rr <- test$massoc$RR.strata.score
	p.value <- test$massoc$chisq.strata$p.value

	counts <- as.vector(contingency)
	names(counts) <- c('low_obese', 'high_obese', 'low_nonobese', 'high_nonobese')

	return(unlist(c(rr, p.value=p.value, counts)))
}

analyze_composite <- function(subset){
	test <- rma(ai=low_obese, bi=low_nonobese, ci=high_obese, di=high_nonobese, data=subset, measure="RR", method="REML")
#	BFRes <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=raw.overall, measure="RR", slab=paste(Study, Year, Total, sep=", "), method="REML")
##tpos <- as.numeric(ShannonRRTable$tposH) #obese in low diversity group
#tneg <- as.numeric(ShannonRRTable$tnegH) #nonobese in low diversity group
#cpos <- as.numeric(ShannonRRTable$cposH) #obese in high diversity group
#cneg <- as.numeric(ShannonRRTable$cnegH) #nonobese in high diversity group

	c(exp(c(rr=test$b[[1,1]], ci_lb=test$ci.lb, ci_ub=test$ci.ub)), p_value=test$pval)
}

datasets <- "baxter escobar hmp goodrich zupancic schubert ross turnbaugh wu"

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

		shannon <- analyze_indiv_study(alpha$shannon, metadata$obese)
		sobs <- analyze_indiv_study(alpha$sobs, metadata$obese)
		shannoneven <- analyze_indiv_study(alpha$shannoneven, metadata$obese)
		bacteroidetes <- analyze_indiv_study(bf_relabund[,"b"], metadata$obese)
		firmicutes <- analyze_indiv_study(bf_relabund[,"f"], metadata$obese)
		bf_ratio <- analyze_indiv_study(bf_relabund[,"bf"], metadata$obese)

		test <- rbind(shannon, sobs, shannoneven, bacteroidetes, firmicutes, bf_ratio)
		study_data <- data.frame(dataset = d, metric = rownames(test), test)

		summary_data <- rbind(summary_data, study_data)
	}

	write.table(summary_data, file="data/process/relative_risk.summary", quote=F, sep='\t', row.names=F)


	shannon <- analyze_composite(summary_data[summary_data$metric=="shannon",])
	shannoneven <- analyze_composite(summary_data[summary_data$metric=="shannoneven",])
	sobs <- analyze_composite(summary_data[summary_data$metric=="sobs",])
	firmicutes <- analyze_composite(summary_data[summary_data$metric=="firmicutes",])
	bacteroidetes <- analyze_composite(summary_data[summary_data$metric=="bacteroidetes",])
	bf_ratio <- analyze_composite(summary_data[summary_data$metric=="bf_ratio",])

	composite <- rbind(shannon, shannoneven, sobs,
											firmicutes, bacteroidetes, bf_ratio)

	write.table(summary_data, file="data/process/relative_risk.composite", quote=F, sep='\t', row.names=F)

}

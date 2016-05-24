library("lme4")

source('code/utilities.R')

composite_analysis <- function(alpha, dataset, is_obese, pow){

	null_model <- lmer(alpha^pow ~ (1|dataset), REML=FALSE)
	obese_model <- lmer(alpha^pow ~ (1|is_obese) + (1|dataset), REML=FALSE)
	anova(null_model, obese_model)[["Pr(>Chisq)"]][2]

}

single_analysis <- function(alpha, pow, is_obese){

	p.value <- t.test(alpha^pow ~ is_obese)$p.value

	mean_overall <- mean(alpha, na.rm=T)
	sd_overall <- sd(alpha, na.rm=T)

	obese <- alpha[is_obese]
	non_obese <- alpha[!is_obese]

	mean_non <- mean(non_obese, na.rm=T)
	mean_obese <- mean(obese, na.rm=T)
	sd_non <- sd(non_obese, na.rm=T)
	sd_obese <- sd(obese, na.rm=T)
	n_non <- length(non_obese)
	n_obese <- length(obese)

	ci_obese <- unname(quantile(obese, c(0.025, 0.975)))
	ci_non <- unname(quantile(non_obese, c(0.025, 0.975)))


	return(
		c(p.value = p.value,
			mean_overall = mean_overall, sd_overall = sd_overall,
			mean_non = mean_non, mean_obese = mean_obese,
			sd_non = sd_non, sd_obese = sd_obese,
			n_non = n_non, n_obese = n_obese, low_ci_obese=ci_obese[1], high_ci_obese=ci_obese[2], low_ci_non=ci_non[1], high_ci_non=ci_non[2]
		)
	)
}

datasets <- c('baxter', 'goodrich', 'escobar', 'hmp', 'ross', 'turnbaugh', 'wu', 'zupancic')

run <- function(datasets){

	datasets <- unlist(strsplit(datasets, split=" "))

	summary_data <- NULL
	composite_data <- NULL

	for(d in datasets){
		print(d)

		alpha_file <- paste0('data/', d, '/', d, '.groups.ave-std.summary')
		alpha <- read.table(file=alpha_file, header=T, stringsAsFactors=F)
		alpha <- alpha[alpha$method == 'ave',]

		metadata_file <- paste0('data/', d, '/', d, '.metadata')
		metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=F)
		metadata <- metadata[metadata$sample %in% alpha$group,]

		stopifnot(alpha$group == metadata$sample)

		na_obesity <- is.na(metadata$obese)
		metadata <- metadata[!na_obesity,]
		alpha <- alpha[!na_obesity,]

		stopifnot(rownames(alpha) == metadata$sample)


		bf_relabund <- get_bacteroides_firmicutes(d)

		#transformations were worked out with a qqplot to check for normality
		shannon <- single_analysis(alpha$shannon, pow=2, metadata$obese)
		sobs <- single_analysis(alpha$sobs, pow=0.5, metadata$obese)
		shannoneven <- single_analysis(alpha$shannoneven, pow=4, metadata$obese)
		bacteroidetes <- single_analysis(bf_relabund[,"b"], pow=0.5, metadata$obese)
		firmicutes <- single_analysis(bf_relabund[,"f"], pow=1, metadata$obese)
		bf_ratio <- single_analysis(bf_relabund[,"bf"], pow=0.5, metadata$obese)

		test <- rbind(shannon, sobs, shannoneven,
									bacteroidetes, firmicutes, bf_ratio)

		study_data <- data.frame(dataset=d, subject=alpha$group,
			 					obese=metadata$obese,
								shannon=alpha$shannon, sobs=alpha$sobs,
								shannoneven=alpha$shannoneven, bacteroidetes=bf_relabund[,"b"],
								firmicutes=bf_relabund[,"f"], bf_ratio=bf_relabund[,"bf"])
		composite_data <- rbind(composite_data, study_data)

		study_summary <- data.frame(dataset = d, metric = rownames(test), test)
		summary_data <- rbind(summary_data, study_summary)
	}

	write.table(composite_data, file="data/process/alpha.data", quote=F, sep='\t', row.names=F)

	write.table(summary_data, file="data/process/alpha_tests.summary", quote=F, sep='\t', row.names=F)

#	composite_analysis <- function(alpha, dataset, is_obese, pow){
	composite_p <- NULL
	composite_p["shannon"] <- composite_analysis(composite_data$shannon, composite_data$dataset, composite_data$obese, 2)
	composite_p["sobs"] <- composite_analysis(composite_data$sobs, composite_data$dataset, composite_data$obese, 0.5)
	composite_p["shannoneven"] <- composite_analysis(composite_data$shannoneven, composite_data$dataset, composite_data$obese, 4)
	composite_p["bacteroidetes"] <- composite_analysis(composite_data$bacteroidetes, composite_data$dataset, composite_data$obese, 0.5)
	composite_p["firmicutes"] <- composite_analysis(composite_data$firmicutes, composite_data$dataset, composite_data$obese, 1)
	composite_p["bf_ratio"] <- composite_analysis(composite_data$bf_ratio, composite_data$dataset, composite_data$obese, 0.5)

	write.table(file="data/process/alpha_composite.summary", data.frame(composite_p), quote=F, sep='\t')
}

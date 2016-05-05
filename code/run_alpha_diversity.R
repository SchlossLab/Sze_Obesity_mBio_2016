source('code/utilities.R')

single_analysis <- function(alpha, is_obese){

	p.value <- wilcox.test(alpha ~ is_obese)$p.value

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

	return(
		c(p.value = p.value,
			mean_overall = mean_overall, sd_overall = sd_overall,
			mean_non = mean_non, mean_obese = mean_obese,
			sd_non = sd_non, sd_obese = sd_obese,
			n_non = n_non, n_obese = n_obese
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

		bf_relabund <- get_bacteroides_firmicutes(d)

		shannon <- single_analysis(alpha$shannon, metadata$obese)
		sobs <- single_analysis(alpha$sobs, metadata$obese)
		shannoneven <- single_analysis(alpha$shannoneven, metadata$obese)
		bacteroidetes <- single_analysis(bf_relabund[,"b"], metadata$obese)
		firmicutes <- single_analysis(bf_relabund[,"f"], metadata$obese)
		bf_ratio <- single_analysis(bf_relabund[,"bf"], metadata$obese)

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
}

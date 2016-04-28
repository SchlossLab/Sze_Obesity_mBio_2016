source('code/utilities.R')

datasets <- c('baxter', 'escobar', 'hmp', 'ross', 'turnbaugh', 'wu')

run <- function(datasets){

	datasets <- unlist(strsplit(datasets, split=" "))

	z_transform <- NULL
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

		shannon_z <- scale(alpha$shannon)
		bf_z <- scale(log(bf_relabund[,"bf"] + 1))
		is_obese <- metadata$obese

		z_transform <- rbind(z_transform, data.frame(study=d, shannon=shannon_z, bf=bf_z, obese=is_obese))
	}

	write.table(z_transform, file="data/process/z_transform.data", quote=F, sep='\t', row.names=F)

	shannon_test <- t.test(z_transform[,"shannon"]~z_transform[,"obese"])
	bf_test <- t.test(z_transform[,"bf"]~z_transform[,"obese"])

	test_summary <- matrix(c(shannon_test$estimate, shannon_test$p.value,
									bf_test$estimate, bf_test$p.value), nrow=2, byrow=T)
	rownames(test_summary) <- c("shannon", "bf")
	colnames(test_summary) <- c("non-obese", "obese", "p.value")

	write.table(test_summary, file="data/process/z_transform.test", quote=F, sep='\t', row.names=T)
}

source('code/utilities.R')
get_dependencies("lme4")

datasets <- c('baxter', 'goodrich', 'escobar', 'hmp', 'ross', 'turnbaugh', 'wu', 'zupancic')

re_model <- function(data, study, is_obese){
	null_model <- lmer(data ~ (1|study), REML=FALSE)
	re_model <- lmer(data ~ (1|is_obese) + (1|study), REML=FALSE)
	test <- anova(null_model, re_model)
	test$Pr[2]
}

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
		sobs_z <- scale(alpha$sobs)
		shannoneven_z <- scale(alpha$shannoneven)
		b_z <- scale(log(bf_relabund[,"b"] + 1))
		f_z <- scale(log(bf_relabund[,"f"] + 1))
		bf_z <- scale(log(bf_relabund[,"bf"] + 1))

		is_obese <- metadata$obese

		z_transform <- rbind(z_transform, data.frame(study=d, shannon=shannon_z, sobs=sobs_z, shannoneven=shannoneven_z, bacteroidetes=b_z, firmicutes=f_z, bf_ratio=bf_z, obese=is_obese))
	}

	write.table(z_transform, file="data/process/z_transform.data", quote=F, sep='\t', row.names=F)

	metrics <- c("shannon", "sobs", "shannoneven", "bacteroidetes", "firmicutes", "bf_ratio")

	p_values <- apply(z_transform[,metrics], 2, re_model, study=z_transform$study, is_obese=z_transform$obese)

	test_summary <- data.frame(metric=names(p_values), p_value=p_values)

	write.table(test_summary, file="data/process/z_transform.test", quote=F, sep='\t', row.names=F)
}

source('code/utilities.R')
get_dependencies(c('pwr','statmod'))


run_alpha <- function(datasets){

	alpha <- read.table(file="data/process/alpha_tests.summary", header=T, stringsAsFactors=F)

	datasets <- unique(alpha$dataset)
	metrics <- unique(alpha$metric)



	# get observed power for each study and metric
	observed_power <- data.frame(
											matrix(nrow=length(datasets), ncol=length(metrics)+3)
										)
	colnames(observed_power) <- c('study', metrics, 'n_non', 'n_obese')

	for(m in metrics){
		metric_table <- alpha[alpha$metric == m,]

		#observed power
		obs_cohens_d <- abs(metric_table$mean_non - metric_table$mean_obese) /
		 													metric_table$sd_overall

		observed_power[,m] <- pwr.t2n.test(d=obs_cohens_d,
																	n1=metric_table$n_non,
																	n2=metric_table$n_obese)$power

		observed_power$n_non <- metric_table$n_non
		observed_power$n_obese <- metric_table$n_obese
	}
	observed_power$study <- datasets

	write.table(observed_power, file="data/process/alpha_power.observed", quote=F, sep='\t', row.names=F)



	# get power of various studies to detect varying effect sizes w/ original
	# design and the required N to have 80% power to detect varying effect sizes
	# w/ a balanced design
	predicted_power <- NULL
	differences <- c(0.01, 0.05, 0.10, 0.15)

	predicted_power <- data.frame(matrix(vector(), nrow=0, ncol=5))
	colnames(predicted_power) <- c('metric', 'study', 'effect_size', 'power', 'balanced_n')

	for(m in metrics){
		metric_table <- alpha[alpha$metric == m,]

		for(diff in differences){
				cohens_d <- diff * metric_table$mean_overall / metric_table$sd_overall

				study_power <- pwr.t2n.test(d=cohens_d,n1=metric_table[,"n_non"], n2=metric_table[,"n_obese"])$power

				needed_n <- sapply(cohens_d,function(x){ pwr.t.test(d=x,power=0.80)$n } )

				temp <- NULL
				temp <- data.frame(metric=rep(m, length(cohens_d)),
													study=metric_table$dataset,
													effect_size=rep(diff, length(cohens_d)), power=study_power, balanced_n=needed_n)

				predicted_power <- rbind(predicted_power, temp)
		}
	}

	write.table(predicted_power, file="data/process/alpha_power.predicted", quote=F, sep='\t', row.names=F)
}




run_rr <- function(datasets){

	rr <- read.table(file="data/process/relative_risk.summary", header=T, stringsAsFactors=F)

	datasets <- unique(rr$dataset)
	metrics <- unique(rr$metric)


	# get observed power for each study and metric
	observed_power <- data.frame(
											matrix(nrow=length(datasets), ncol=length(metrics)+3)
										)
	colnames(observed_power) <- c('study', metrics, 'n_non', 'n_obese')
	observed_power$study <- datasets

	for(m in metrics){
		metric_table <- rr[rr$metric == m,]

		#observed power
		p_obese <- (metric_table$high_obese)/
								(metric_table$high_obese + metric_table$low_obese)
		p_non_obese <- (metric_table$high_nonobese)/
										(metric_table$high_nonobese + metric_table$low_nonobese)
		n_obese <- metric_table$high_obese + metric_table$low_obese
		n_non_obese <- metric_table$high_nonobese + metric_table$low_nonobese

		observed_power[, m] <- pwr.2p2n.test(ES.h(p_obese, p_non_obese),
				                      n_obese, n_non_obese,
				                      alternative="two.sided")$power

	}
	observed_power$n_non <- n_non_obese
	observed_power$n_obese <- n_obese

	write.table(observed_power, file="data/process/rr_power.observed", quote=F, sep='\t', row.names=F)



	# get power of various studies to detect varying effect sizes w/ original
	# design and the required N to have 80% power to detect varying effect sizes
	# w/ a balanced design
	differences <- c(0.01, 0.05, 0.10, 0.15)

	predicted_power <- data.frame(matrix(vector(), nrow=0, ncol=5))
	colnames(predicted_power) <- c('metric', 'study', 'effect_size', 'power', 'balanced_n')

	for(m in metrics){
		metric_table <- rr[rr$metric == m,]

		p_obese <- (metric_table$high_obese)/
								(metric_table$high_obese + metric_table$low_obese)
		n_obese <- metric_table$high_obese + metric_table$low_obese
		n_non_obese <- metric_table$high_nonobese + metric_table$low_nonobese

		for(diff in differences){

			study_power <- pwr.2p2n.test(h=ES.h(p_obese, p_obese + diff),
																	n1=n_obese, n2=n_non_obese,
																	alternative="two.sided")$power

			needed_n <- sapply(p_obese, function(x){
											pwr.2p.test(h=ES.h(x, x + diff), power=0.80)$n
										})





#			needed_n <- sapply(cohens_d,function(x){ pwr.t.test(d=x,power=0.80)$n } )

			temp <- data.frame(metric=rep(m, length(needed_n)),
												study=metric_table$dataset,
												effect_size=rep(diff, length(needed_n)), power=study_power, balanced_n=needed_n)

			predicted_power <- rbind(predicted_power, temp)
		}
	}

	write.table(predicted_power, file="data/process/rr_power.predicted", quote=F, sep='\t', row.names=F)
}

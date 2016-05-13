make_label <- function(string, p){
	p_value <- paste0("=", format(round(p, digits=3), nsmall=3L))
	p_value[p_value == "=0.000" | p_value == "=0"] <- "<0.001"
	paste(string, " (P", p_value, ")", sep="")
}


make_study_label <- function(data){
	options(scipen=10000)
	study <- ifelse(data$dataset=="hmp", "HMP", capwords(data$dataset))
	make_label(study,data$p.value)
}

make_metric_label <- function(metric){

	formatted_metrics <- c(shannon = "Shannon Diversity Index",
												shannoneven = "Shannon Evenness Index",
												sobs = "Number of Observed OTUs",
												bacteroidetes = "Relative Abundance of Bacteroidetes",
												firmicutes = "Relative Abundance of Firmicutes",
												bf_ratio = "Ratio of Bacteroidetes to Firmicutes"
											)
	formatted_metrics[metric]
}

capwords <- function(s, strict = FALSE) {
		s <- as.character(s)
		cap <- function(s) paste(toupper(substring(s, 1, 1)),
									{s <- substring(s, 2); if(strict) tolower(s) else s},
														 sep = "", collapse = " " )
		sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


diversity_plot <- function(metric, leg=NULL){

	gap <- 0.1

	alpha <- read.table(file="data/process/alpha_tests.summary", header=T)
	alpha$se_non <- alpha$sd_non/sqrt(alpha$n_non)
	alpha$se_obese <- alpha$sd_obese/sqrt(alpha$n_obese)

	composite <- read.table(file="data/process/alpha_composite.summary", header=T)

	datasets <- levels(alpha$dataset)

	subset <- alpha[alpha$metric == metric,]
	subset <- subset[order(subset$dataset),]
	stopifnot(subset$dataset == datasets)


	par(mar=c(5,8,0.5,0.5))
	plot(NA, ylim=c(0,length(datasets)), xlim=c(0,max(subset$mean_non+1.95*subset$se_non, subset$mean_obese+1.95*subset$se_obese)), axes=F, xlab=make_metric_label(metric), ylab="")

	arrows(x0=subset$mean_obese, x1=subset$mean_obese+1.95*subset$se_obese, y0=length(datasets):1+gap, y1=length(datasets):1+gap, angle=90, length=gap/2)
	arrows(x0=subset$mean_obese, x1=subset$mean_obese-1.95*subset$se_obese, y0=length(datasets):1+gap, y1=length(datasets):1+gap, angle=90, length=gap/2)

	arrows(x0=subset$mean_non, x1=subset$mean_non+1.95*subset$se_non, y0=length(datasets):1-gap, y1=length(datasets):1-gap, angle=90, length=gap/2)
	arrows(x0=subset$mean_non, x1=subset$mean_non-1.95*subset$se_non, y0=length(datasets):1-gap, y1=length(datasets):1-gap, angle=90, length=gap/2)

	points(x=subset$mean_non, y=length(datasets):1-gap, pch=21, bg="red")
	points(x=subset$mean_obese, y=length(datasets):1+gap, pch=21, bg="blue")

	points(x=c(mean(subset$mean_non), mean(subset$mean_obese)), y=c(-gap,gap), pch=21, bg=c("red", "blue"), xpd=T)

	axis(1)
	axis(2, at=length(datasets):1, label=make_study_label(subset), las=2, cex.axis=0.8)

	axis(2, at=0, label=make_label("Overall", composite[metric,1]), las=2, cex.axis=0.8)

	box()

	if(!is.null(leg)){
		legend(x=leg[1], y=leg[2], legend=c("Non-obese", "Obese"), pch=21, pt.bg=c("blue", "red"), cex=0.8)
	}

}


pdf(file="alpha.pdf")
diversity_plot('shannon', leg=c(0,9))
diversity_plot('shannoneven', leg=c(0,9))
diversity_plot('sobs', leg=c(0,9))
diversity_plot('bacteroidetes', leg=c(0,9))
diversity_plot('firmicutes', leg=c(0,9))
diversity_plot('bf_ratio', leg=c(4,9))
dev.off()

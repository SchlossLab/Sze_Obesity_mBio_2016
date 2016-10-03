make_label <- function(string, p){
	p_value <- paste0("=", format(round(p, digits=3), nsmall=3L))
	p_value[p_value == "=0.000" | p_value == "=0"] <- "<0.001"
	paste(string, " (P", p_value, ")", sep="")
}


make_study_label <- function(dataset){
#	options(scipen=10000)
	study <- ifelse(dataset=="hmp", "HMP", capwords(dataset))
#	make_label(study,data$p.value)
}

make_metric_label <- function(metric){

	formatted_metrics <- c(shannon = "Shannon Diversity Index",
												shannoneven = "Shannon\nEvenness Index",
												sobs = "Number of\nObserved OTUs",
												bacteroidetes = "Relative Abundance\nof Bacteroidetes",
												firmicutes = "Relative Abundance\nof Firmicutes",
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


diversity_plot <- function(metric){

	gap <- 0.1

	alpha <- read.table(file="data/process/alpha_tests.summary", header=T)
	alpha$se_non <- alpha$sd_non/sqrt(alpha$n_non)
	alpha$se_obese <- alpha$sd_obese/sqrt(alpha$n_obese)

	composite <- read.table(file="data/process/alpha_composite.summary", header=T)

	datasets <- levels(alpha$dataset)

	subset <- alpha[alpha$metric == metric,]
	subset <- subset[order(subset$dataset),]
	stopifnot(subset$dataset == datasets)


	par(mar=c(5,0.5,2,0.5))
	plot(NA, ylim=c(0,length(datasets)), xlim=c(0,max(subset$mean_non+1.95*subset$se_non, subset$mean_obese+1.95*subset$se_obese)), axes=F, xlab=make_metric_label(metric), ylab="")

	arrows(x0=subset$mean_obese, x1=subset$mean_obese+1.95*subset$se_obese, y0=length(datasets):1+gap, y1=length(datasets):1+gap, angle=90, length=gap/2)
	arrows(x0=subset$mean_obese, x1=subset$mean_obese-1.95*subset$se_obese, y0=length(datasets):1+gap, y1=length(datasets):1+gap, angle=90, length=gap/2)

	arrows(x0=subset$mean_non, x1=subset$mean_non+1.95*subset$se_non, y0=length(datasets):1-gap, y1=length(datasets):1-gap, angle=90, length=gap/2)
	arrows(x0=subset$mean_non, x1=subset$mean_non-1.95*subset$se_non, y0=length(datasets):1-gap, y1=length(datasets):1-gap, angle=90, length=gap/2)

	points(x=subset$mean_non, y=length(datasets):1-gap, pch=21, bg="red")
	points(x=subset$mean_obese, y=length(datasets):1+gap, pch=21, bg="blue")

	sig_star <- ifelse(subset$p.value < 0.05, "*", "")
	text(x=0, y=length(datasets):1, labels=sig_star, cex=3)

	points(x=c(mean(subset$mean_non), mean(subset$mean_obese)), y=c(-gap,gap), pch=21, bg=c("red", "blue"), xpd=T)

	if(composite[composite$metric == metric,"p"] < 0.05){
		text(x=0, y=0, labels="*", cex=3)
	}

	axis(1)
	axis(2, at=0:nrow(subset), labels=FALSE, las=2, cex.axis=0.8)

#	axis(2, at=0, label=make_label("Overall", composite[metric,1]), las=2, cex.axis=0.8)

	box()

	return(subset$dataset)
}

build_figure <- function(metrics, width=6.5, height=3.75, leg=c(8, 10)){

	n_metrics <- length(metrics)

	output_file <- paste0("results/figures/", paste(metrics, collapse="_"), ".tiff")
	tiff(file=output_file, width=width, height=height, units='in', res=300)

	axis_width <- 0.12 * n_metrics
	layout(matrix(c(n_metrics+1, seq(1:(n_metrics))), nrow=1), width=c(axis_width, rep(1, n_metrics)))

	datasets <- NULL

	for(m in 1:n_metrics){
		datasets <- diversity_plot(metrics[m])
		mtext(side=2, line=-1, at=par()$usr[4]*1.05, text=LETTERS[m], cex=1, las=2, font=2)
	}

	legend(x=leg[1], y=leg[2], legend=c("Non-obese", "Obese"), pch=21, pt.bg=c("red", "blue"), cex=1.0)

	par(mar=c(5,0.5,2,0.5))
	plot(NA, xlim=c(0,1), ylim=c(0,length(datasets)),axes=F, xlab="", ylab="")

	text(x=0.1, y=length(datasets):1, labels=make_study_label(datasets), las=2, cex=1.0, adj=0, xpd=T)
	text(x=0.1, y=0, labels="Overall", las=2, cex=1.0, adj=0, xpd=T)

	dev.off()

}

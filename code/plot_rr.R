make_label <- function(string, p){
	p_value <- paste0("=", format(round(p, digits=3), nsmall=3L))
	p_value[p_value == "=0.000"] <- "<0.001"
	paste(string, " (P", p_value, ")", sep="")
}

make_study_label <- function(datasets){
#	options(scipen=10000)
	study <- ifelse(datasets=="hmp", "HMP", capwords(datasets))

#	make_label(study, data$p.value)
}

make_metric_label <- function(metric){

	formatted_metrics <- c(shannon = "Relative Risk\nShannon Diversity Index",
												shannoneven = "Relative Risk\nShannon Evenness\nIndex",
												sobs = "Relative Risk\nNumber of\nObserved OTUs",
												bacteroidetes = "Relative Risk\nRelative Abundance\nof Bacteroidetes",
												firmicutes = "Relative Risk\nRelative Abundance\nof Firmicutes",
												bf_ratio = "Relative Risk\nRatio of Bacteroidetes to Firmicutes"
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


rr_plot <- function(metric){

	gap <- 0.15

	rr <- read.table(file="data/process/relative_risk.summary", header=T)
	datasets <- levels(rr$dataset)

	rr_composite <- read.table(file="data/process/relative_risk.composite", header=T)

	subset <- rr[rr$metric == metric,]
	subset <- subset[order(subset$dataset),]
	stopifnot(subset$dataset == datasets)

	par(mar=c(5,0.5,0.5,0.5))
	plot(NA, ylim=c(0,length(datasets)), xlim=c(min(subset$lower),max(subset$upper)), axes=F, xlab="", ylab="", log="x")

	mtext(side=1, text=make_metric_label(metric), line=4, cex=0.7)
	abline(v=1)

	arrows(x0=subset$est, x1=subset$upper, y0=length(datasets):1, y1=length(datasets):1, angle=90, length=gap/2)
	arrows(x0=subset$est, x1=subset$lower, y0=length(datasets):1, y1=length(datasets):1, angle=90, length=gap/2)

	points(x=subset$est, y=length(datasets):1, pch=19)

	sig_star <- ifelse(subset$p.value < 0.05, "*", "")
	text(x=min(subset$lower), y=length(datasets):1, labels=sig_star, cex=3)

	points(x=rr_composite[metric,"rr"], y=0, pch=19)
	arrows(x0=rr_composite[metric,"rr"], x1=rr_composite[metric,"ci_lb"], y0=0, y1=0, angle=90, length=gap/2)
	arrows(x0=rr_composite[metric,"rr"], x1=rr_composite[metric,"ci_ub"], y0=0, y1=0, angle=90, length=gap/2)

	if(rr_composite[metric,"p_value"] < 0.05){
		text(x=min(subset$lower), y=0, labels="*", cex=3)
	}


	axis(1)
	axis(2, at=0:nrow(subset), labels=FALSE)
	box()

	return(datasets)
}

build_figure <- function(metrics, width=6.5, height=3.75){
	n_metrics <- length(metrics)

	output_file <- paste0("results/figures/rr_", paste(metrics, collapse="_"), ".tiff")
	tiff(file=output_file, width=width, height=height, units='in', res=300)

	axis_width <- 0.12 * n_metrics
	layout(matrix(c(n_metrics+1, seq(1:(n_metrics))), nrow=1), width=c(axis_width, rep(1, n_metrics)))

	datasets <- NULL

	for(m in metrics){
		datasets <- rr_plot(m)
	}

	par(mar=c(5,0.5,0.5,0.5))
	plot(NA, xlim=c(0,1), ylim=c(0,length(datasets)),axes=F, xlab="", ylab="")

	text(x=0.1, y=length(datasets):1, labels=make_study_label(datasets), las=2, cex=1.0, adj=0, xpd=T)
	text(x=0.1, y=0, labels="Overall", las=2, cex=1.0, adj=0, xpd=T)

	dev.off()

}

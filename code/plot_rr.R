make_label <- function(string, p){
	p_value <- paste0("=", format(round(p, digits=3), nsmall=3L))
	p_value[p_value == "=0.000"] <- "<0.001"
	paste(string, " (P", p_value, ")", sep="")
}

make_study_label <- function(data){
	options(scipen=10000)
	study <- ifelse(data$dataset=="hmp", "HMP", capwords(data$dataset))

	make_label(study, data$p.value)
}

make_metric_label <- function(metric){

	formatted_metrics <- c(shannon = "Relative Risk\nShannon Diversity Index",
												shannoneven = "Relative Risk\nShannon Evenness Index",
												sobs = "Relative Risk\nNumber of Observed OTUs",
												bacteroidetes = "Relative Risk\nRelative Abundance of Bacteroidetes",
												firmicutes = "Relative Risk\nRelative Abundance of Firmicutes",
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

metric <- 'shannon'

rr_plot <- function(metric){

	gap <- 0.15

	rr <- read.table(file="data/process/relative_risk.summary", header=T)
	datasets <- levels(rr$dataset)

	rr_composite <- read.table(file="data/process/relative_risk.composite", header=T)

	subset <- rr[rr$metric == metric,]
	subset <- subset[order(subset$dataset),]
	stopifnot(subset$dataset == datasets)

	par(mar=c(5,8,0.5,0.5))
	plot(NA, ylim=c(0,length(datasets)), xlim=c(min(subset$lower),max(subset$upper)), axes=F, xlab=make_metric_label(metric), ylab="", log="x")

	abline(v=1)

	arrows(x0=subset$est, x1=subset$upper, y0=length(datasets):1, y1=length(datasets):1, angle=90, length=gap/2)
	arrows(x0=subset$est, x1=subset$lower, y0=length(datasets):1, y1=length(datasets):1, angle=90, length=gap/2)

	points(x=subset$est, y=length(datasets):1, pch=19)

	points(x=rr_composite[metric,"rr"], y=0, pch=19)
	arrows(x0=rr_composite[metric,"rr"], x1=rr_composite[metric,"ci_lb"], y0=0, y1=0, angle=90, length=gap/2)
	arrows(x0=rr_composite[metric,"rr"], x1=rr_composite[metric,"ci_ub"], y0=0, y1=0, angle=90, length=gap/2)

	axis(1)
	axis(2, at=length(datasets):0, label=c(make_study_label(subset), make_label("Overall", rr_composite[metric,"p_value"])), las=2, cex.axis=0.8)
	box()

}

pdf(file="rr.pdf")
rr_plot('shannon')
rr_plot('shannoneven')
rr_plot('sobs')
rr_plot('bacteroidetes')
rr_plot('firmicutes')
rr_plot('bf_ratio')
dev.off()

layout(matrix(c(1,2), nrow=1))
rr_plot('shannon')
rr_plot('bf_ratio')
layout(1)

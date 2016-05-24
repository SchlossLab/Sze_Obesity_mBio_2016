capwords <- function(s, strict = FALSE) {
		cap <- function(s) paste(toupper(substring(s, 1, 1)),
									{s <- substring(s, 2); if(strict) tolower(s) else s},
														 sep = "", collapse = " " )
		sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

plot_roc_data <- function(data_type){

	roc_curve_file <- paste0("data/process/random_forest.", data_type, ".roc_data")
	roc_curve <- read.table(file=roc_curve_file, header=T, stringsAsFactors=F)

	rf_summary_file <- paste0("data/process/random_forest.", data_type, ".summary")
	rf_summary <- read.table(file=rf_summary_file, header=T, stringsAsFactors=F)

	N <- ifelse(data_type == "otu", "OTUs", "Genera")

	rf_summary <- rf_summary[order(rf_summary$auc_cv, decreasing=T),]
	datasets <- rf_summary$dataset

	study_name <- capwords(datasets)
	study_name[study_name == "Hmp"] <- "HMP"

	basic_colors <- c("red", "lightseagreen", "dodgerblue", "black")
	colors <- rep(basic_colors, length.out=length(datasets))
	names(colors) <- datasets

	basic_lty <- 1:3
	lty <- rep(basic_lty, each=length(basic_colors), length.out=length(datasets))
	names(lty) <- datasets

	plot_roc <- function(study){
		study_data <- roc_curve[roc_curve$dataset == study,]
		points(study_data$specificity~study_data$sensitivity, type="s", col=colors[study], lty=lty[study], lwd=2)
	}


	pdf(file=paste0("results/figures/roc_curve.", data_type, ".pdf"))
	par(mar=c(4,4,0.5,0.5))
	plot(NA, xlim=c(1,0), ylim=c(0,1), axes=F, ylab="Specificity", xlab="Sensitivity")

	abline(a=1, b=-1, col="gray")

	sapply(datasets, plot_roc)

	axis(1)
	axis(2, las=2)
	box()

	legend_name <- paste(
		study_name,
		" (AUC-CV=",
		format(rf_summary$auc_cv, digits=2),
		"; ", rf_summary$k_opt, " ", N, ")", sep="")

	legend(x=0.6, y=0.35, legend=legend_name, lty=lty, col=colors, lwd=2, cex=0.9)
	dev.off()
}

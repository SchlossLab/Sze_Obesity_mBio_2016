pch <- c(baxter=21, escobar=21, goodrich=21, hmp=21, ross=21,
				schubert=21, turnbaugh=21, wu=21, zeevi=21, zupancic=21)
col <- c(baxter="black", escobar="red", goodrich="green", hmp="blue",
 				ross="orange", schubert="black", turnbaugh="red", wu="green",
				zeevi="blue", zupancic="orange")
bg <- c(baxter="white", escobar="white", goodrich="white", hmp="white",
 				ross="white", schubert="black", turnbaugh="red", wu="green",
				zeevi="blue", zupancic="orange")
lty <- c(baxter=2, escobar=2, goodrich=2, hmp=2, ross=2,
				schubert=1, turnbaugh=1, wu=1, zeevi=1, zupancic=1)

capwords <- function(s, strict = FALSE) {
		cap <- function(s) paste(toupper(substring(s, 1, 1)),
									{s <- substring(s, 2); if(strict) tolower(s) else s},
														 sep = "", collapse = " " )
		sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


plot_roc_data <- function(data_type, double=TRUE){

	roc_curve_file <- paste0("data/process/random_forest.", data_type, ".roc_data")
	roc_curve <- read.table(file=roc_curve_file, header=T, stringsAsFactors=F)

	rf_summary_file <- paste0("data/process/random_forest.", data_type, ".summary")
	rf_summary <- read.table(file=rf_summary_file, header=T, stringsAsFactors=F)

	N <- ifelse(data_type == "otu", "OTUs", "Genera")

	rf_summary <- rf_summary[order(rf_summary$auc_cv, decreasing=T),]
	datasets <- sort(rf_summary$dataset)

	study_name <- capwords(datasets)
	study_name[study_name == "Hmp"] <- "HMP"

	plot_roc <- function(study){
		study_data <- roc_curve[roc_curve$dataset == study,]
		points(study_data$specificity~study_data$sensitivity, type="s", col=col[study], lty=lty[study], lwd=2)
	}

	if(double == FALSE){
		par(mar=c(4,4,0.5,0.5))
		plot(NA, xlim=c(1,0), ylim=c(0,1), axes=F, ylab="Specificity", xlab="Sensitivity")
		axis(2, las=2)
	} else {
		par(mar=c(4,0.5,0.5,0.5))
		plot(NA, xlim=c(1,0), ylim=c(0,1), axes=F, ylab="", xlab="Sensitivity")
		axis(2, las=2, labels=NA)
	}

	abline(a=1, b=-1, col="gray")
	axis(1)
	box()

	sapply(datasets, plot_roc)



	legend_name <- paste(
		study_name,
		" (AUC=",
		format(rf_summary$auc_cv, digits=2), ")", sep="")
		#"; ", rf_summary$k_opt, " ", N, ")", sep="")


	legend('bottomright', legend=legend_name, lty=lty, col=col, lwd=2, cex=0.9)
}


build_figure <- function(){

	pdf(file=paste0("results/figures/roc_curve.pdf"), width=7.5, height=4)
	layout(matrix(c(3,1,2), nrow=1), width=c(0.15, 1,1))
	plot_roc_data("otu")
	plot_roc_data("genus")

	par(mar=c(4,0,0.5,0.5))
	plot(NA, ylim=c(0,1), xlim=c(0,1), axes=F, xlab="", ylab="")
	text(x=rep(0.8,6), y=seq(0,1,0.2), labels=format(seq(0,1,0.2), digits=2))
	text(x=0.2, y=0.5, labels="Specificity", srt=90)
	dev.off()

}

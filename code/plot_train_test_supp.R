pch <- c(baxter=21, escobar=21, goodrich=21, hmp=21, ross=21,
				schubert=21, turnbaugh=21, wu=21, zeevi=21, zupancic=21)
col <- c(baxter="black", escobar="red", goodrich="green", hmp="blue",
 				ross="orange", schubert="black", turnbaugh="red", wu="green",
				zeevi="blue", zupancic="orange")
bg <- c(baxter="white", escobar="white", goodrich="white", hmp="white",
 				ross="white", schubert="black", turnbaugh="red", wu="green",
				zeevi="blue", zupancic="orange")
names <- c(baxter="Baxter", escobar="Escobar", goodrich="Goodrich", hmp="HMP",
 				ross="Ross", schubert="Schubert", turnbaugh="Turnbaugh", wu="Wu",
				zeevi="Zeevi", zupancic="Zupancic")



plot_train_test <- function(rank){

	train_test_file <- paste0("data/process/random_forest.", tolower(rank),
													".train_test")
	z <- read.table(file=train_test_file, header=T, stringsAsFactors=F)

	o <- order(z$train, z$test)
	z <- z[o,]

	par(mar=c(0.5,1,1.5,0.5))
	plot(NA, xlim=c(0.9,10.1), ylim=c(0,1), ylab="", xlab="", axes=F)
	abline(h=0.5, col="gray", lwd=2, lty=2)
	abline(v=seq(1.5,9.5,1))
	studies <- sort(unique(z$train))
	n_studies <- length(studies)

	stagger <- seq(-0.3,0.3,length.out=10)
	for(s in 1:n_studies){
		study <- z[z$test==studies[s],]
	#	segments(x0=s+stagger, y0=study$accuracy_lci, y1=study$accuracy_uci)
		points(x=s+stagger, y=study$accuracy, col=col[studies], bg=bg[studies],
						lwd=2, pch=pch[studies])
	}
	box()
	axis(2, las=2)
	mtext(side=3, text=rank, at=0.5, adj=0, font=2)


	summary_file <- paste0("data/process/random_forest.", tolower(rank),
													".summary")
	s <- read.table(file=summary_file, header=T, stringsAsFactors=F)
	s <- s[order(s$dataset),]
	stopifnot(s$dataset==studies)
	text(x=1:10, y=rep(0,10), labels=format(s$auc_cv, digits=2))

	return(studies)
}


tiff(file="results/figures/train_test_supp.tiff", width=6.75, height=10, units='in', res=300)

layout(matrix(c(6,1,
							6,2,
							6,3,
							6,4,
							6,5,
							0,7), nrow=6, byrow=TRUE), width=c(0.1,1), height=c(1,1,1,1,1,0.2))

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")
studies <- NULL

for(r in ranks){
	studies <- plot_train_test(r)
}

mtext(1, at=1:length(studies), text=names[studies], cex=0.7, font=2, line=0.5)

plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
text(0.25, 0.5, "Accuracy", srt=90, cex=2)

par(mar=c(0,1,0,0.5))
plot(x=1:10, y=rep(0.3, 10), xlim=c(0.9,10.1), ylim=c(0,1), pch=pch[studies],
			bg=bg[studies], col=col[studies], cex=2, lwd=2, axes=F, xlab="", ylab="")

dev.off()

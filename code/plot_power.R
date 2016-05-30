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



build_plots <- function(method){

	pred <- read.table(file=paste0("data/process/", method, "_power.predicted"), header=T, stringsAsFactors=FALSE)

	metrics <- unique(pred$metric)

	for(m in metrics){

		pred_subset <- pred[pred$metric == m,]
		o <- order(pred_subset$effect_size, pred_subset$study)
		pred_subset <- pred_subset[o,]

		effects <- unique(pred_subset$effect_size)
		n_effects <- length(effects)

		studies <- unique(pred_subset$study)
		n_studies <- length(studies)

		stagger <- seq(-0.3,0.3,length.out=n_studies)

		pdf_file <- paste0("results/figures/", method, "_", m, "_power.pdf")

		pdf(file=pdf_file, width=6.0, height=5)
		layout(matrix(c(1,1,3,2,2,3,0,0,0), nrow=3, byrow=T), width=c(1,1,0.4), height=c(1,1,0.2))

		par(mar=c(0.5,5,0.5,0.5))
		plot(NA, xlim=c(0.7,4.3), ylim=c(0,1), ylab="Power to Detect Effect Size\nWith Original Sampling Effort", xlab="", axes=F)

		for(e in 1:n_effects){

			effect <- pred_subset[pred_subset$effect_size==effects[e],]
			points(x=e+stagger, y=effect$power, col=col[effect$study],
				bg=bg[effect$study], lwd=2, pch=pch[effect$study])

		}

		axis(1, at=1:n_effects, labels=FALSE)
		axis(2, las=2)
		box()
		mtext(side=2, at=1.0, line=3, text="A", las=2, font=2, cex=2)
		abline(v=c(1.5, 2.5, 3.5))


		par(mar=c(0.5,5,0.5,0.5))
		plot(NA, xlim=c(0.7,4.3), ylim=c(1,max(pred_subset$balanced_n)), ylab="Number of Samples\nNeeded per Group", xlab="", axes=F, log='y')

		for(e in 1:n_effects){

			effect <- pred_subset[pred_subset$effect_size==effects[e],]
			points(x=e+stagger, y=effect$balanced_n, col=col[effect$study],
				bg=bg[effect$study], lwd=2, pch=pch[effect$study])

		}

		axis(1, at=1:n_effects, labels=100*effects)
		axis(2, las=2)
		box()
		mtext(side=2, at=1.2*(10^par()$usr[4]), line=3, text="B", las=2, font=2, cex=2)

		if(method == 'alpha'){
			mtext(1, line=2, text = "Effect Size (%)", cex=0.7)
		} else {
			mtext(1, line=2, text = "Effect Size (Cohen's d)", cex=0.7)
		}

		abline(v=c(1.5, 2.5, 3.5))

		par(mar=c(0,0,0,0))
		plot(NA, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="")
		legend(x=0.1, y=0.66, legend=names, pch=pch, col=col, pt.bg=bg, pt.cex=1.5)
		dev.off()

	}

}

pdf(file="results/figures/train_test.pdf", width=10, height=4)

par(mar=c(4,4,0.5,0.5))
z <- read.table(file="data/process/random_forest.genus.train_test", header=T, stringsAsFactors=F)

o <- order(z$train, z$test)
z <- z[o,]

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

plot(NA, xlim=c(0.9,10.1), ylim=c(0,1), ylab="Accuracy", xlab="", axes=F)
abline(h=0.5, col="gray", lwd=2, lty=2)
abline(v=seq(1.5,9.5,1))
studies <- sort(unique(z$train))
n_studies <- length(studies)

stagger <- seq(-0.3,0.3,length.out=10)
for(s in 1:n_studies){
	study <- z[z$test==studies[s],]
#	segments(x0=s+stagger, y0=study$accuracy_lci, y1=study$accuracy_uci)
	points(x=s+stagger, y=study$accuracy, col=col[studies], bg=bg[studies], lwd=2,
		pch=pch[studies])

}

axis(1, at=1:length(studies), labels=names[studies], cex.axis=0.8)
axis(2, las=2)

points(x=1:10, y=rep(-0.23, 10), pch=pch[studies], bg=bg[studies], col=col[studies], xpd=T, cex=2, lwd=2)
box()

dev.off()

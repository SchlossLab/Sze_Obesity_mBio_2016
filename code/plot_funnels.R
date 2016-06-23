pretty_metric <- c("shannon" = "Shannon", "sobs" = "Richness",
									"shannoneven" = "Evenness",
									"bf_ratio" = "Bacteroidetes:Firmicutes",
									"bacteroidetes" = "Bacteroidetes",
									"firmicutes" = "Firmicutes")

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

my_funnel <- function(metric, rr){
	metric_frame <- rr[rr$metric == metric, ]

	metric_frame$se <- (log(metric_frame$upper) - log(metric_frame$est)) / 1.96

	overall <- mean(log(metric_frame$est))
	se <- seq(0,1,0.01)
	l <- overall - 1.96 * se
	u <- overall + 1.96 * se

	par(mar=c(2.5,2,1.5,0.5))
	plot(metric_frame$se~log(metric_frame$est), ylim=c(1,0), ylab="Standard Error", pch=pch[metric_frame$dataset], col=col[metric_frame$dataset], bg=bg[metric_frame$dataset])
	abline(v=overall)

	points(se~l, type="l")
	points(se~u, type="l")
	text(x=par()$usr[1], y=-0.1, label=pretty_metric[metric], xpd=TRUE, font=2, cex=1.3, pos=4)
}

rr <- read.table(file="data/process/relative_risk.summary", header=T)

calculators <- c("shannon", "sobs", "shannoneven", "bf_ratio", "bacteroidetes", "firmicutes")


tiff('results/figures/funnel_plot.tiff', width=6.5, height=6, unit='in', res=300)
layout(matrix(c(8,8,8,0,1,2,3,7,4,5,6,7,9,9,9,0), nrow=4), widths=c(0.1,1,1,0.4), heights=c(1,1,1,0.1))
sapply(calculators, my_funnel, rr)


par(mar=c(0,0,0,0))
plot.new()
text(x=0.5,y=0.5, label="log(Relative Risk)", cex=1.5)

par(mar=c(0,0,0,0))
plot.new()
text(x=0.5,y=0.5, label="Standard Error", srt=90, cex=1.5)

plot.new()
legend(x=0, y=0.65, legend=names, col=col, pt.bg=bg, pch=pch, cex=1.2)
dev.off()

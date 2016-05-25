## Make Custom RR Figures for Manuscript
## Obesesity and the bacterial microbiome
## Marc Sze
## March 24, 2016

library(ggplot2)
library(metafor)
library(wesanderson)

# Adaption of code from Joey ggforest wrapper

# First need to load and create the necessary rma class file
# For Shannon Diversity

ShannonRRTable <- read.csv("results/tables/denovoShannonRRTable.csv")

StudyNo <- as.numeric(ShannonRRTable$X)
Study <- as.character(ShannonRRTable$Study)
Year <- c("2015", "2015", "2014", "2014", "2012", "2011", "2012", "2009")
Total <- c("n=172", "n=63", "n=507", "n=30", "n=200", "n=256", "n=63", "n=146")
tpos <- as.numeric(ShannonRRTable$tposH) #obese in low diversity group
tneg <- as.numeric(ShannonRRTable$tnegH) #nonobese in low diversity group
cpos <- as.numeric(ShannonRRTable$cposH) #obese in high diversity group
cneg <- as.numeric(ShannonRRTable$cnegH) #nonobese in high diversity group
RR <- as.numeric(ShannonRRTable$RRH)
low <- as.numeric(ShannonRRTable$lowH)
high <- as.numeric(ShannonRRTable$highH)

numeric.data <-as.data.frame(cbind(tpos, tneg, cpos, cneg))
raw.overall <- as.data.frame(cbind(StudyNo, Study,Year, Total, numeric.data))

dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data=raw.overall, append=TRUE)
k <- length(dat$StudyNo)
dat.fm <- data.frame(study=factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(raw.overall, c(rbind(tpos, tneg, cpos, cneg)))

#Used to generate the RE Model estimate based on all the studies
ShannonRes <- rma(ai=tpos, bi=tneg, ci=cpos, 
                  di=cneg, data=raw.overall, measure="RR", 
                  slab=paste(Study, Year, Total, sep=", "), method="REML")
RRSDPvalue <- anova.rma(ShannonRes)

# One important distinction to make is that the graph now used the logFC
# values which is pretty much the log RR values (it is close but not quite
# exact)  Thus difference from 0.0 is what is more important.

test <- rbind(data.frame(
  Study = "RE Model", LogFC = ShannonRes$b, CILB=ShannonRes$ci.lb, 
  CIUB=ShannonRes$ci.ub,
  p = ShannonRes$pval, stringsAsFactors = FALSE), 
  data.frame(
    Study = ShannonRes$slab, LogFC = ShannonRes$yi, 
    CILB=ShannonRes$yi - 2*sqrt(ShannonRes$vi),
    CIUB=ShannonRes$yi + 2*sqrt(ShannonRes$vi), 
    p = ShannonRes$pval,
    stringsAsFactors = FALSE))

test2 <- transform(test, interval = CIUB - CILB)
test3R <- transform(test2, RelConf = 1/interval)
test3R$Study[1] <- "Pooled Outcome"
test3R$Study <- factor(test3R$Study, 
                      levels = c(
                        "Pooled Outcome", "HMP, 2011, n=256", "Turnbaugh, 2009, n=146", 
                        "Wu, 2012, n=63", "Zupancic, 2012, n=200", "Escobar, 2014, n=30", 
                        "Ross, 2015, n=63", "Goodrich, 2014, n=507", "Baxter, 2015, n=172"
                      ))
test3R$Region <- c("NA", "V4", "V1-V3", "V4", "V1-V3", "V1-V3", "V3-V5", "V2", "V2")
test3R$Region <- as.factor(test3R$Region)

# Make the actual graph

RRTest <- ggplot(test3R, aes(LogFC, Study, xmax=CIUB, xmin=CILB))
HRRG1 <- RRTest + coord_cartesian(xlim=c(-2.2, 2.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(aes(colour=Region), alpha=0.5, height=0) + 
  geom_errorbarh(data = subset(test3R, Study=="Pooled Outcome"), 
                 colour="black", alpha=0.5, height=0) + 
  geom_point(aes(colour=Region), size = 3) + ylab("") + 
  xlab("Log Relative Risk") + 
  geom_point(data = subset(test3R, Study=="Pooled Outcome"), 
             colour="black", size=5) + 
  scale_color_manual(values = wes_palette("FantasticFox")) + 
  scale_size(range=c(2, 5), guide=FALSE) + 
  theme(axis.line=element_line(colour="black"), 
        axis.line.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.title=element_text(size=12, face="bold"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
        text = element_text(size=15), 
        axis.text.y = element_text(face=c(
          "bold", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain"), 
          size=c(12, 10, 10, 10, 10, 10, 10, 10, 10)), 
        legend.position="none") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))

# First need to load and create the necessary rma class file
# For Bacteroidetes/Firmicutes Ratio

BFRatioRRTable <- read.csv("results/tables/denovoBFRatioRRTable.csv")

StudyNo <- as.numeric(BFRatioRRTable$X)
Study <- as.character(BFRatioRRTable$Study)
Year <- c("2015", "2015", "2014", "2014", "2012", "2011", "2012", "2009")
Total <- c("n=172", "n=63", "n=507", "n=30", "n=200", "n=256", "n=63", "n=146")
tpos <- as.numeric(BFRatioRRTable$tposBF) #obese in low diversity group
tneg <- as.numeric(BFRatioRRTable$tnegBF) #nonobese in low diversity group
cpos <- as.numeric(BFRatioRRTable$cposBF) #obese in high diversity group
cneg <- as.numeric(BFRatioRRTable$cnegBF) #nonobese in high diversity group
RR <- as.numeric(BFRatioRRTable$RRBF)
low <- as.numeric(BFRatioRRTable$lowBF)
high <- as.numeric(BFRatioRRTable$highBF)

numeric.data <-as.data.frame(cbind(tpos, tneg, cpos, cneg))
raw.overall <- as.data.frame(cbind(StudyNo, Study,Year, Total, numeric.data))

dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data=raw.overall, append=TRUE)
k <- length(dat$StudyNo)
dat.fm <- data.frame(study=factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(raw.overall, c(rbind(tpos, tneg, cpos, cneg)))

#Used to generate the RE Model estimate based on all the studies
BFRes <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=raw.overall, measure="RR", slab=paste(Study, Year, Total, sep=", "), method="REML")
RRBFPvalue <- anova.rma(BFRes)


# One important distinction to make is that the graph now used the logFC
# values which is pretty much the log RR values (it is close but not quite
# exact)  Thus difference from 0.0 is what is more important.

test <- rbind(data.frame(
  Study = "RE Model", LogFC = BFRes$b, CILB=BFRes$ci.lb, CIUB=BFRes$ci.ub,
  p = BFRes$pval, stringsAsFactors = FALSE), 
  data.frame(
    Study = BFRes$slab, LogFC = BFRes$yi, 
    CILB=BFRes$yi - 2*sqrt(BFRes$vi), CIUB=BFRes$yi + 2*sqrt(BFRes$vi), 
    p = BFRes$pval, stringsAsFactors = FALSE))

test2 <- transform(test, interval = CIUB - CILB)
test3BF <- transform(test2, RelConf = 1/interval)
test3BF$Study[1] <- "Pooled Outcome"
test3BF$Study <- factor(test3BF$Study, 
                      levels = c(
                        "Pooled Outcome", "HMP, 2011, n=256", "Turnbaugh, 2009, n=146", 
                        "Wu, 2012, n=63", "Zupancic, 2012, n=200", "Escobar, 2014, n=30", 
                        "Ross, 2015, n=63", "Goodrich, 2014, n=507", "Baxter, 2015, n=172"
                      ))
test3BF$Region <- c("NA", "V4", "V1-V3", "V4", "V1-V3", "V1-V3", "V3-V5", "V2", "V2")
test3BF$Region <- as.factor(test3BF$Region)

# Make the actual graph

BFRRTest <- ggplot(test3BF, aes(LogFC, Study, xmax=CIUB, xmin=CILB))
BFRRG1 <- BFRRTest + coord_cartesian(xlim=c(-2.2, 2.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(aes(colour=Region), alpha=0.5, height=0) + 
  geom_errorbarh(data = subset(test3BF, Study=="Pooled Outcome"), 
                 colour="black", alpha=0.5, height=0) + 
  geom_point(aes(colour=Region), size = 3) + ylab("") + 
  xlab("Log Relative Risk") + 
  geom_point(data = subset(test3BF, Study=="Pooled Outcome"), 
             colour="black", size=5) + 
  scale_color_manual(values = wes_palette("FantasticFox")) + 
  scale_size(range=c(2, 5), guide=FALSE) + 
  theme(axis.line=element_line(colour="black"), 
        axis.line.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.title=element_text(size=12, face="bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size=15), 
        axis.text.y = element_text(face=c(
          "bold", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain"), 
          size=c(12, 10, 10, 10, 10, 10, 10, 10, 10)), 
        legend.position="none") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


## Make Custom Figures for Supplemental
## Obesesity and the bacterial microbiome
## Marc Sze
## March 28, 2016

library(ggplot2)

# Loading in respective Data

OverallPValue <- read.csv("results/tables/denovoOverallPTable.csv")
PowerTable <- read.csv("results/tables/denovoPowerTable.csv")
ModelledPower1PerDiff <- read.csv("results/tables/denovoModelled1PerDiff.csv")
ModelledPower5PerDiff <- read.csv("results/tables/denovoModelled5PerDiff.csv")
ModelledPower10PerDiff <- read.csv("results/tables/denovoModelled10PerDiff.csv")
ModelledPower15PerDiff <- read.csv("results/tables/denovoModelled15PerDiff.csv")
ActualDiff <- read.csv("results/tables/denovoActualDiff.csv")
SimulatedN_B <- read.csv("results/tables/denovoSimulatedN_B.csv")
SimulatedN_F <- read.csv("results/tables/denovoSimulatedN_F.csv")
SimulatedN_S <- read.csv("results/tables/denovoSimulatedN_S.csv")
SimulatedN_J <- read.csv("results/tables/denovoSimulatedN_J.csv")
ShannonRRTable <- read.csv("results/tables/denovoShannonRRTable.csv")
BFRatioRRTable <- read.csv("results/tables/denovoBFRatioRRTable.csv")

##################### OTU RICHNESS ###########################################

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffS <- as.numeric(format(round(ActualDiff$S, 3), nsmall=3))
SPvalue <- as.numeric(format(round(OverallPValue$OTURich, 2), nsmall=2))
OneP <- as.numeric(format(round(ModelledPower1PerDiff$S, 3), nsmall=3))
FiveP <- as.numeric(format(round(ModelledPower5PerDiff$S, 3), nsmall=3))
TenP <- as.numeric(format(round(ModelledPower10PerDiff$S, 3), nsmall=3))
FifteenP <- as.numeric(format(round(ModelledPower15PerDiff$S, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(SPvalue, ActualDiffS, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", SPvalue, sep="")
DataTable1$Legend <- TableLegend

library(reshape)
test <- melt(DataTable1, id=c("ActualDiffS", "SPvalue", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
s1 <- ggplot(test, aes(EstDiff, value))
s2 <- s1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("S (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.key = element_rect(fill = "white")) + 
  geom_vline(aes(xintercept=ActualDiffS[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffS[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededS <- as.numeric(format(round(SimulatedN_S$SActualNeededN, 
                                          3), nsmall=3))
OnePNS <- as.numeric(format(round(SimulatedN_S$SOnePer, 3), nsmall=3))
FivePNS <- as.numeric(format(round(SimulatedN_S$SFivePer, 3), nsmall=3))
TenPNS <- as.numeric(format(round(SimulatedN_S$STenPer, 3), nsmall=3))
FifteenPNS <- as.numeric(format(round(SimulatedN_S$SFifteenPer, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededS, OnePNS, FivePNS, TenPNS, FifteenPNS))
DataTable2$Study <- Study
sampleNS <- melt(DataTable2, id=c("ActualNNeededS", "Study"))

sampleNS$EstDiff[sampleNS$variable=="OnePNS"] <- 1
sampleNS$EstDiff[sampleNS$variable=="FivePNS"] <- 5
sampleNS$EstDiff[sampleNS$variable=="TenPNS"] <- 10
sampleNS$EstDiff[sampleNS$variable=="FifteenPNS"] <- 15

LogActualNNeededS <- as.numeric(format(round(log10(SimulatedN_S$SActualNeededN), 
                                             3), nsmall=3))

set.seed(100)
sn1 <- ggplot(sampleNS, aes(EstDiff, log10(value)))
sn2 <- sn1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("S (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.position="none") + 
  ylim(0, 5) +  
  geom_hline(aes(yintercept=LogActualNNeededS[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededS[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


##################### Evenness ###########################################

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffJ <- as.numeric(format(round(ActualDiff$J, 3), nsmall=3))
JPvalue <- as.numeric(format(round(OverallPValue$Evenness, 2), nsmall=2))
OneP <- as.numeric(format(round(ModelledPower1PerDiff$J, 3), nsmall=3))
FiveP <- as.numeric(format(round(ModelledPower5PerDiff$J, 3), nsmall=3))
TenP <- as.numeric(format(round(ModelledPower10PerDiff$J, 3), nsmall=3))
FifteenP <- as.numeric(format(round(ModelledPower15PerDiff$J, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(JPvalue, ActualDiffJ, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", JPvalue, sep="")
DataTable1$Legend <- TableLegend

test <- melt(DataTable1, id=c("ActualDiffJ", "JPvalue", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
j1 <- ggplot(test, aes(EstDiff, value))
j2 <- j1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("J (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.key = element_rect(fill = "white")) + 
  geom_vline(aes(xintercept=ActualDiffJ[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffJ[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededJ <- as.numeric(format(round(SimulatedN_J$JActualNeededN, 
                                          3), nsmall=3))
OnePNJ <- as.numeric(format(round(SimulatedN_J$JOnePer, 3), nsmall=3))
FivePNJ <- as.numeric(format(round(SimulatedN_J$JFivePer, 3), nsmall=3))
TenPNJ <- as.numeric(format(round(SimulatedN_J$JTenPer, 3), nsmall=3))
FifteenPNJ <- as.numeric(format(round(SimulatedN_J$JFifteenPer, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededJ, OnePNJ, FivePNJ, TenPNJ, FifteenPNJ))
DataTable2$Study <- Study
sampleNJ <- melt(DataTable2, id=c("ActualNNeededJ", "Study"))

sampleNJ$EstDiff[sampleNJ$variable=="OnePNJ"] <- 1
sampleNJ$EstDiff[sampleNJ$variable=="FivePNJ"] <- 5
sampleNJ$EstDiff[sampleNJ$variable=="TenPNJ"] <- 10
sampleNJ$EstDiff[sampleNJ$variable=="FifteenPNJ"] <- 15

LogActualNNeededJ <- as.numeric(format(round(log10(SimulatedN_J$JActualNeededN), 
                                             3), nsmall=3))

set.seed(100)
jn1 <- ggplot(sampleNJ, aes(EstDiff, log10(value)))
jn2 <- jn1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                         position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("J (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.position="none") + 
  ylim(0, 5) +  
  geom_hline(aes(yintercept=LogActualNNeededJ[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededJ[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


##################### FIRMICUTES ###########################################

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffF <- as.numeric(format(round(ActualDiff$firm, 3), nsmall=3))
FPvalue <- as.numeric(format(round(OverallPValue$Firmicutes, 2), nsmall=2))
OneP <- as.numeric(format(round(ModelledPower1PerDiff$firm, 3), nsmall=3))
FiveP <- as.numeric(format(round(ModelledPower5PerDiff$firm, 3), nsmall=3))
TenP <- as.numeric(format(round(ModelledPower10PerDiff$firm, 3), nsmall=3))
FifteenP <- as.numeric(format(round(ModelledPower15PerDiff$firm, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(FPvalue, ActualDiffF, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", FPvalue, sep="")
DataTable1$Legend <- TableLegend

test <- melt(DataTable1, id=c("ActualDiffF", "FPvalue", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
f1 <- ggplot(test, aes(EstDiff, value))
f2 <- f1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("Firmicutes (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.key = element_rect(fill = "white")) + 
  geom_vline(aes(xintercept=ActualDiffF[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffF[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededF <- as.numeric(format(round(SimulatedN_F$FirmActualNeededN, 
                                          3), nsmall=3))
OnePNF <- as.numeric(format(round(SimulatedN_F$FirmOnePer, 3), nsmall=3))
FivePNF <- as.numeric(format(round(SimulatedN_F$FirmFivePer, 3), nsmall=3))
TenPNF <- as.numeric(format(round(SimulatedN_F$FirmTenPer, 3), nsmall=3))
FifteenPNF <- as.numeric(format(round(SimulatedN_F$FirmFifteenPer, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededF, OnePNF, FivePNF, TenPNF, FifteenPNF))
DataTable2$Study <- Study
sampleNF <- melt(DataTable2, id=c("ActualNNeededF", "Study"))

sampleNF$EstDiff[sampleNF$variable=="OnePNF"] <- 1
sampleNF$EstDiff[sampleNF$variable=="FivePNF"] <- 5
sampleNF$EstDiff[sampleNF$variable=="TenPNF"] <- 10
sampleNF$EstDiff[sampleNF$variable=="FifteenPNF"] <- 15

LogActualNNeededF <- as.numeric(format(round(log10(SimulatedN_F$FirmActualNeededN), 
                                             3), nsmall=3))

set.seed(100)
fn1 <- ggplot(sampleNF, aes(EstDiff, log10(value)))
fn2 <- fn1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                         position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("Firmicutes (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.position="none") + 
  ylim(0, 5) +  
  geom_hline(aes(yintercept=LogActualNNeededF[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededF[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


##################### Bacteroidetes ###########################################

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffB <- as.numeric(format(round(ActualDiff$bacter, 3), nsmall=3))
BPvalue <- as.numeric(format(round(OverallPValue$Bacteroidetes, 2), nsmall=2))
OneP <- as.numeric(format(round(ModelledPower1PerDiff$bacter, 3), nsmall=3))
FiveP <- as.numeric(format(round(ModelledPower5PerDiff$bacter, 3), nsmall=3))
TenP <- as.numeric(format(round(ModelledPower10PerDiff$bacter, 3), nsmall=3))
FifteenP <- as.numeric(format(round(ModelledPower15PerDiff$bacter, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(BPvalue, ActualDiffB, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", BPvalue, sep="")
DataTable1$Legend <- TableLegend

test <- melt(DataTable1, id=c("ActualDiffB", "BPvalue", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
b1 <- ggplot(test, aes(EstDiff, value))
b2 <- b1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("Bacteroidetes (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.key = element_rect(fill = "white")) + 
  geom_vline(aes(xintercept=ActualDiffB[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffB[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededB <- as.numeric(format(round(SimulatedN_B$BacterActualNeededN, 
                                          3), nsmall=3))
OnePNB <- as.numeric(format(round(SimulatedN_F$FirmOnePer, 3), nsmall=3))
FivePNB <- as.numeric(format(round(SimulatedN_F$FirmFivePer, 3), nsmall=3))
TenPNB <- as.numeric(format(round(SimulatedN_F$FirmTenPer, 3), nsmall=3))
FifteenPNB <- as.numeric(format(round(SimulatedN_F$FirmFifteenPer, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededB, OnePNB, FivePNB, TenPNB, FifteenPNB))
DataTable2$Study <- Study
sampleNB <- melt(DataTable2, id=c("ActualNNeededB", "Study"))

sampleNB$EstDiff[sampleNB$variable=="OnePNB"] <- 1
sampleNB$EstDiff[sampleNB$variable=="FivePNB"] <- 5
sampleNB$EstDiff[sampleNB$variable=="TenPNB"] <- 10
sampleNB$EstDiff[sampleNB$variable=="FifteenPNB"] <- 15

LogActualNNeededB <- as.numeric(format(round(log10(SimulatedN_B$BacterActualNeededN), 
                                             3), nsmall=3))

set.seed(100)
bn1 <- ggplot(sampleNB, aes(EstDiff, log10(value)))
bn2 <- bn1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                         position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("Bacteroidetes (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.position="none") + 
  ylim(0, 5) +  
  geom_hline(aes(yintercept=LogActualNNeededB[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededB[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


######################## RELATIVE RISK  ###############################################

# Load in specific Data Needed
PowerModelHRR <- read.csv("results/tables/denovoPowerModelHRR.csv")
PowerModelBFRR <- read.csv("results/tables/denovoPowerModelBFRR.csv")
RRSimNNeededH <- read.csv("results/tables/denovoRRSimNNeededH.csv")
RRSimNNeededBF <- read.csv("results/tables/denovoRRSimNNeededBF.csv")

##RR Shannon Diversity

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffRRH <- as.numeric(format(round(PowerModelHRR$ActualDiffHRR*100, 3), nsmall=3))
HRRPvalue <- as.numeric(format(round(ShannonRRTable$pValueH, 2), nsmall=2))
OneP <- as.numeric(format(round(PowerModelHRR$OnePercentHRR, 3), nsmall=3))
FiveP <- as.numeric(format(round(PowerModelHRR$FivePercentHRR, 3), nsmall=3))
TenP <- as.numeric(format(round(PowerModelHRR$TenPercentHRR, 3), nsmall=3))
FifteenP <- as.numeric(format(round(PowerModelHRR$FifteenPercentHRR, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(ActualDiffRRH, HRRPvalue, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", HRRPvalue, sep="")
DataTable1$Legend <- TableLegend

test <- melt(DataTable1, id=c("ActualDiffRRH", "HRRPvalue", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
RRH1 <- ggplot(test, aes(EstDiff, value))
RRH2 <- RRH1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Legend") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("H RR (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.key = element_rect(fill = "white")) + 
  geom_vline(aes(xintercept=ActualDiffRRH[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRH[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededRRH <- as.numeric(format(round(RRSimNNeededH$RRActSimNH, 
                                          3), nsmall=3))
OnePNHRR <- as.numeric(format(round(RRSimNNeededH$RROnePercentNH, 3), nsmall=3))
FivePNHRR <- as.numeric(format(round(RRSimNNeededH$RRFivePercentNH, 3), nsmall=3))
TenPNHRR <- as.numeric(format(round(RRSimNNeededH$RRTenPercentNH, 3), nsmall=3))
FifteenPNHRR <- as.numeric(format(round(RRSimNNeededH$RRFifteenPercentNH, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededRRH, OnePNHRR, FivePNHRR, TenPNHRR, FifteenPNHRR))
DataTable2$Study <- Study
sampleNHRR <- melt(DataTable2, id=c("ActualNNeededRRH", "Study"))

sampleNHRR$EstDiff[sampleNHRR$variable=="OnePNHRR"] <- 1
sampleNHRR$EstDiff[sampleNHRR$variable=="FivePNHRR"] <- 5
sampleNHRR$EstDiff[sampleNHRR$variable=="TenPNHRR"] <- 10
sampleNHRR$EstDiff[sampleNHRR$variable=="FifteenPNHRR"] <- 15

LogActualNNeededHRR <- as.numeric(format(round(log10(RRSimNNeededH$RRActSimNH), 
                                             3), nsmall=3))

set.seed(100)
HRRn1 <- ggplot(sampleNHRR, aes(EstDiff, log10(value)))
HRRn2 <- HRRn1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                         position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("H RR (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.position="none") + 
  ylim(0, 5) +  
  geom_hline(aes(yintercept=LogActualNNeededHRR[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededHRR[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


##RR BF Ratio

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffRRBF <- as.numeric(format(round(PowerModelBFRR$ActualDiffBFRR*100, 3), nsmall=3))
BFRRPvalue <- as.numeric(format(round(BFRatioRRTable$pValueBF, 2), nsmall=2))
OneP <- as.numeric(format(round(PowerModelBFRR$OnePercentBFRR, 3), nsmall=3))
FiveP <- as.numeric(format(round(PowerModelBFRR$FivePercentBFRR, 3), nsmall=3))
TenP <- as.numeric(format(round(PowerModelBFRR$TenPercentBFRR, 3), nsmall=3))
FifteenP <- as.numeric(format(round(PowerModelBFRR$FifteenPercentBFRR, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(ActualDiffRRBF, BFRRPvalue, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", BFRRPvalue, sep="")
DataTable1$Legend <- TableLegend

test <- melt(DataTable1, id=c("ActualDiffRRBF", "BFRRPvalue", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
RRBF1 <- ggplot(test, aes(EstDiff, value))
RRBF2 <- RRBF1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                           position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Legend") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("BF Ratio RR (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.key = element_rect(fill = "white")) + 
  geom_vline(aes(xintercept=ActualDiffRRBF[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffRRBF[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededRRBF <- as.numeric(format(round(RRSimNNeededBF$RRActSimNBF, 3), nsmall=3))
OnePNBFRR <- as.numeric(format(round(RRSimNNeededBF$RROnePercentNBF, 3), nsmall=3))
FivePNBFRR <- as.numeric(format(round(RRSimNNeededBF$RRFivePercentNBF, 3), nsmall=3))
TenPNBFRR <- as.numeric(format(round(RRSimNNeededBF$RRTenPercentNBF, 3), nsmall=3))
FifteenPNBFRR <- as.numeric(format(round(RRSimNNeededBF$RRFifteenPercentNBF, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededRRBF, OnePNBFRR, FivePNBFRR, TenPNBFRR, FifteenPNBFRR))
DataTable2$Study <- Study
sampleNBFRR <- melt(DataTable2, id=c("ActualNNeededRRBF", "Study"))

sampleNBFRR$EstDiff[sampleNBFRR$variable=="OnePNBFRR"] <- 1
sampleNBFRR$EstDiff[sampleNBFRR$variable=="FivePNBFRR"] <- 5
sampleNBFRR$EstDiff[sampleNBFRR$variable=="TenPNBFRR"] <- 10
sampleNBFRR$EstDiff[sampleNBFRR$variable=="FifteenPNBFRR"] <- 15

LogActualNNeededBFRR <- as.numeric(format(round(log10(RRSimNNeededBF$RRActSimNBF), 
                                               3), nsmall=3))

set.seed(100)
BFRRn1 <- ggplot(sampleNBFRR, aes(EstDiff, log10(value)))
BFRRn2 <- BFRRn1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                             position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("BF Ratio RR (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                     axis.line=element_line(colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(), 
                     legend.position="none") + 
  ylim(0, 5) +  
  geom_hline(aes(yintercept=LogActualNNeededBFRR[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBFRR[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


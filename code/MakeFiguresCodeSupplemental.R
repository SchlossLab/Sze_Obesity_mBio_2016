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








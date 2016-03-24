## Make Figures for Manuscript
## Obesesity and the bacterial microbiome
## Marc Sze
## March 23, 2016

library(ggplot2)


OverallPValue <- read.csv("results/tables/denovoOverallPTable.csv")
PowerTable <- read.csv("results/tables/denovoPowerTable.csv")
ModelledPower1PerDiff <- read.csv("results/tables/denovoModelled1PerDiff.csv")
ModelledPower5PerDiff <- read.csv("results/tables/denovoModelled5PerDiff.csv")
ModelledPower10PerDiff <- read.csv("results/tables/denovoModelled10PerDiff.csv")
ModelledPower15PerDiff <- read.csv("results/tables/denovoModelled15PerDiff.csv")
ActualDiff <- read.csv("results/tables/denovoActualDiff.csv")
SimulatedN_H <- read.csv("results/tables/denovoSimulatedN_H.csv")
SimulatedN_BF <- read.csv("results/tables/denovoSimulatedN_BF.csv")

##################### SHANNON DIVERSITY ###################################

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffH <- as.numeric(format(round(ActualDiff$H, 3), nsmall=3))
Shannon <- as.numeric(format(round(OverallPValue$Shannon, 2), nsmall=2))
OneP <- as.numeric(format(round(ModelledPower1PerDiff$H, 3), nsmall=3))
FiveP <- as.numeric(format(round(ModelledPower5PerDiff$H, 3), nsmall=3))
TenP <- as.numeric(format(round(ModelledPower10PerDiff$H, 3), nsmall=3))
FifteenP <- as.numeric(format(round(ModelledPower15PerDiff$H, 3), nsmall=3))

DataTable1 <- as.data.frame(
  cbind(Shannon, ActualDiffH, OneP, FiveP, TenP, FifteenP))
DataTable1$Study <- Study
TableLegend <- paste(as.character(Study), ", P-value=", Shannon, sep="")
DataTable1$Legend <- TableLegend

library(reshape)
test <- melt(DataTable1, id=c("ActualDiffH", "Shannon", "Study", "Legend"))

test$EstDiff[test$variable=="OneP"] <- 1
test$EstDiff[test$variable=="FiveP"] <- 5
test$EstDiff[test$variable=="TenP"] <- 10
test$EstDiff[test$variable=="FifteenP"] <- 15

set.seed(100)
a1 <- ggplot(test, aes(EstDiff, value))
a2 <- a1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
           axis.title=element_text(size=10, face="bold")) + 
  xlab("H (Effect Size)") + ylab("Estimated Power") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                              axis.line=element_line(colour="black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank()) + 
  geom_vline(aes(xintercept=ActualDiffH[1], colour=Legend[1]), 
                linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffH[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))

  
# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededH <- as.numeric(format(round(SimulatedN_H$HActualNeededN, 
                                          3), nsmall=3))
OnePNH <- as.numeric(format(round(SimulatedN_H$HOnePer, 3), nsmall=3))
FivePNH <- as.numeric(format(round(SimulatedN_H$HFivePer, 3), nsmall=3))
TenPNH <- as.numeric(format(round(SimulatedN_H$HTenPer, 3), nsmall=3))
FifteenPNH <- as.numeric(format(round(SimulatedN_H$HFifteenPer, 3), nsmall=3))

DataTable2 <- as.data.frame(
  cbind(ActualNNeededH, OnePNH, FivePNH, TenPNH, FifteenPNH))
DataTable2$Study <- Study
sampleNH <- melt(DataTable2, id=c("ActualNNeededH", "Study"))

sampleNH$EstDiff[sampleNH$variable=="OnePNH"] <- 1
sampleNH$EstDiff[sampleNH$variable=="FivePNH"] <- 5
sampleNH$EstDiff[sampleNH$variable=="TenPNH"] <- 10
sampleNH$EstDiff[sampleNH$variable=="FifteenPNH"] <- 15

LogActualNNeededH <- as.numeric(format(round(log10(SimulatedN_H$HActualNeededN), 
                                          3), nsmall=3))

set.seed(100)
b1 <- ggplot(sampleNH, aes(EstDiff, log10(value)))
b2 <- b1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("H (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                              axis.line=element_line(colour="black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank()) + 
  ylim(0, 5) + 
  geom_hline(aes(yintercept=LogActualNNeededH[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededH[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))



##################### Bacteroidetes/Firmicutes #############################

# A data frame that has all the necessary data for power
Study <- as.character(PowerTable$Study)
ActualDiffBF <- as.numeric(format(round(ActualDiff$BFratio, 3), nsmall=3))
BFratio <- as.numeric(format(round(OverallPValue$BFRatio, 2), nsmall=2))
OnePBF <- as.numeric(format(round(ModelledPower1PerDiff$BFratio, 
                                  3), nsmall=3))
FivePBF <- as.numeric(format(round(ModelledPower5PerDiff$BFratio, 
                                   3), nsmall=3))
TenPBF <- as.numeric(format(round(ModelledPower10PerDiff$BFratio, 
                                  3), nsmall=3))
FifteenPBF <- as.numeric(format(round(ModelledPower15PerDiff$BFratio, 
                                      3), nsmall=3))

DataTable3 <- as.data.frame(
  cbind(ActualDiffBF, OnePBF, FivePBF, TenPBF, FifteenPBF))
DataTable3$Study <- Study

TableLegend <- paste(as.character(Study), ", P-value=", BFratio, sep="")
DataTable3$Legend <- TableLegend
sampleBF <- melt(DataTable3, id=c("ActualDiffBF", "Study", "Legend"))

sampleBF$EstDiff[sampleBF$variable=="OnePBF"] <- 1
sampleBF$EstDiff[sampleBF$variable=="FivePBF"] <- 5
sampleBF$EstDiff[sampleBF$variable=="TenPBF"] <- 10
sampleBF$EstDiff[sampleBF$variable=="FifteenPBF"] <- 15
sampleBF <- rbind(sampleBF, c(NA, NA, NA, NA, NA, 30))

set.seed(25)
c1 <- ggplot(sampleBF, aes(EstDiff, value))
c2 <- c1 + geom_jitter(aes(colour = Legend), size=5, alpha=7/10, 
                       position = position_jitter(width=1.0)) + 
  scale_colour_brewer(palette="Set1", name = "Study") + 
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("BF Ratio (Effect Size)") + ylab("Estimated Power") + ylim(0,1) + 
  scale_x_discrete(breaks=c(1, 5, 10, 15, 30), 
                   labels=c("1%", "5%", "10%", "15%", "30%")) + theme(
                              axis.line=element_line(colour="black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank()) + 
  geom_vline(aes(xintercept=ActualDiffBF[1], colour=Legend[1]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[2], colour=Legend[2]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[3], colour=Legend[3]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[4], colour=Legend[4]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[5], colour=Legend[5]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[6], colour=Legend[6]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[7], colour=Legend[7]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=ActualDiffBF[8], colour=Legend[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=17.5), colour="black", linetype="solid") + 
  ggtitle('A') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


# A data frame that has all the necessary data for number sims
Study <- as.character(PowerTable$Study)
ActualNNeededBF <- as.numeric(format(round(SimulatedN_BF$BFActualNeededN, 
                                          3), nsmall=3))
OnePNBF <- as.numeric(format(round(SimulatedN_BF$BFratioOnePer, 3), nsmall=3))
FivePNBF <- as.numeric(format(round(SimulatedN_BF$BFratioFivePer, 3), 
                              nsmall=3))
TenPNBF <- as.numeric(format(round(SimulatedN_BF$BFratioTenPer, 3), 
                             nsmall=3))
FifteenPNBF <- as.numeric(format(round(SimulatedN_BF$BFratioFifteenPer, 3), 
                                 nsmall=3))

DataTable4 <- as.data.frame(
  cbind(ActualNNeededBF, OnePNBF, FivePNBF, TenPNBF, FifteenPNBF))
DataTable4$Study <- Study
sampleNBF <- melt(DataTable4, id=c("ActualNNeededBF", "Study"))

sampleNBF$EstDiff[sampleNBF$variable=="OnePNBF"] <- 1
sampleNBF$EstDiff[sampleNBF$variable=="FivePNBF"] <- 5
sampleNBF$EstDiff[sampleNBF$variable=="TenPNBF"] <- 10
sampleNBF$EstDiff[sampleNBF$variable=="FifteenPNBF"] <- 15

LogActualNNeededBF <- as.numeric(format(round(
  log10(SimulatedN_BF$BFActualNeededN), 3), nsmall=3))

set.seed(100)
d1 <- ggplot(sampleNBF, aes(EstDiff, log10(value)))
d2 <- d1 + geom_jitter(aes(colour = Study), size=5, alpha=7/10, 
                       position = position_jitter(width=1.1)) + 
  scale_colour_brewer(palette="Set1", name = "Study") +  
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10, face="bold")) + 
  xlab("BF Ratio (Effect Size)") + ylab("Log Estimated Number Needed") + 
  scale_x_discrete(breaks=c(1, 5, 10, 15), 
                   labels=c("1%", "5%", "10%", "15%")) + theme(
                              axis.line=element_line(colour="black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank()) + 
  ylim(0, 7.5) + 
  geom_hline(aes(yintercept=LogActualNNeededBF[1], colour=Study[1]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[2], colour=Study[2]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[3], colour=Study[3]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[4], colour=Study[4]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[5], colour=Study[5]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[6], colour=Study[6]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[7], colour=Study[7]), 
             linetype="dashed") + 
  geom_hline(aes(yintercept=LogActualNNeededBF[8], colour=Study[8]), 
             linetype="dashed") + 
  geom_vline(aes(xintercept=2.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=7.5), colour="black", linetype="solid") + 
  geom_vline(aes(xintercept=12.5), colour="black", linetype="solid") + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, face="bold", size=20))


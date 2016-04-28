## Power Calculations Variations: Estimating needed n
## Obesesity and the bacterial microbiome: Supplemental variables
## Bacteroidetes, Firmicutes, OTU Richness, Evenness
## Marc Sze
## March 28, 2016


library(statmod)
library(pwr)
source("code/UsedFunctions.R")

ShannonRRTable <- read.csv("results/tables/denovoShannonRRTable.csv")
BFRatioRRTable <- read.csv("results/tables/denovoBFRatioRRTable.csv")
PowerTable <- read.csv("results/tables/denovoPowerTable.csv")
SummaryStatsByObeseGroup <- read.csv("results/tables/denovoSSbyObeseGroup.csv")
SummaryStatsByStudy <- read.csv("results/tables/denovoSSbyStudy.csv")
ActualDiff <- read.csv("results/tables/denovoActualDiff.csv")

Study <- as.character(PowerTable$Study)

######### ACTUAL ###########################################################


BacterActualNeededN <- rep(0, 8)

for(i in 1:8){
  cohensD <- (SummaryStatsByObeseGroup$MeanNonObeseB[i] -  
                SummaryStatsByObeseGroup$MeanObeseB[i]) / 
    SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  BacterActualNeededN[i] <- pwr.t.test(d=cohensD, 
                                  power=0.8, sig.level=0.05, type="two.sample", 
                                  alternative="two.sided")$n
}

FirmActualNeededN <- rep(0, 8)

for(i in 1:8){
  cohensD <- abs((SummaryStatsByObeseGroup$MeanNonObeseF[i] -  
                    SummaryStatsByObeseGroup$MeanObeseF[i])) / 
    SummaryStatsByStudy$sdF[i]  
  # This needs to change based on study
  
  FirmActualNeededN[i] <- pwr.t.test(d=cohensD, 
                                   power=0.8, sig.level=0.05, type="two.sample", 
                                   alternative="two.sided")$n
}


SActualNeededN <- rep(0, 8)

for(i in 1:8){
  cohensD <- abs((SummaryStatsByObeseGroup$MeanNonObeseS[i] -  
                    SummaryStatsByObeseGroup$MeanObeseS[i])) / 
    SummaryStatsByStudy$sdS[i]  
  # This needs to change based on study
  
  SActualNeededN[i] <- pwr.t.test(d=cohensD, 
                                     power=0.8, sig.level=0.05, type="two.sample", 
                                     alternative="two.sided")$n
}


JActualNeededN <- rep(0, 8)

for(i in 1:8){
  cohensD <- abs((SummaryStatsByObeseGroup$MeanNonObeseJ[i] -  
                    SummaryStatsByObeseGroup$MeanObeseJ[i])) / 
    SummaryStatsByStudy$sdJ[i]  
  # This needs to change based on study
  
  JActualNeededN[i] <- pwr.t.test(d=cohensD, 
                                  power=0.8, sig.level=0.05, type="two.sample", 
                                  alternative="two.sided")$n
}

########## ONE PERCENT #####################################################


OnePercentDiff <- SummaryStatsByStudy$averageStudyB * 0.01
BacterOnePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  BacterOnePer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


OnePercentDiff <- SummaryStatsByStudy$averageStudyF * 0.01
FirmOnePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  FirmOnePer[i] <- pwr.t.test(d=cohensD, 
                                 power=0.8, sig.level=0.05, type="two.sample", 
                                 alternative="two.sided")$n
}


OnePercentDiff <- SummaryStatsByStudy$averageStudyS * 0.01
SOnePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  SOnePer[i] <- pwr.t.test(d=cohensD, 
                              power=0.8, sig.level=0.05, type="two.sample", 
                              alternative="two.sided")$n
}


OnePercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.01
JOnePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  JOnePer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


########## FIVE PERCENT #####################################################

FivePercentDiff <- SummaryStatsByStudy$averageStudyB * 0.05
BacterFivePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  BacterFivePer[i] <- pwr.t.test(d=cohensD, 
                                power=0.8, sig.level=0.05, type="two.sample", 
                                alternative="two.sided")$n
}


FivePercentDiff <- SummaryStatsByStudy$averageStudyF * 0.05
FirmFivePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  FirmFivePer[i] <- pwr.t.test(d=cohensD, 
                              power=0.8, sig.level=0.05, type="two.sample", 
                              alternative="two.sided")$n
}


FivePercentDiff <- SummaryStatsByStudy$averageStudyS * 0.05
SFivePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  SFivePer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


FivePercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.05
JFivePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  JFivePer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}

########## TEN PERCENT #####################################################

TenPercentDiff <- SummaryStatsByStudy$averageStudyB * 0.10
BacterTenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  BacterTenPer[i] <- pwr.t.test(d=cohensD, 
                                 power=0.8, sig.level=0.05, type="two.sample", 
                                 alternative="two.sided")$n
}


TenPercentDiff <- SummaryStatsByStudy$averageStudyF * 0.10
FirmTenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  FirmTenPer[i] <- pwr.t.test(d=cohensD, 
                               power=0.8, sig.level=0.05, type="two.sample", 
                               alternative="two.sided")$n
}


TenPercentDiff <- SummaryStatsByStudy$averageStudyS * 0.10
STenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  STenPer[i] <- pwr.t.test(d=cohensD, 
                            power=0.8, sig.level=0.05, type="two.sample", 
                            alternative="two.sided")$n
}


TenPercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.10
JTenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  JTenPer[i] <- pwr.t.test(d=cohensD, 
                            power=0.8, sig.level=0.05, type="two.sample", 
                            alternative="two.sided")$n
}


########## FIFTEEN PERCENT #####################################################

FifteenPercentDiff <- SummaryStatsByStudy$averageStudyB * 0.10
BacterFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  BacterFifteenPer[i] <- pwr.t.test(d=cohensD, 
                                power=0.8, sig.level=0.05, type="two.sample", 
                                alternative="two.sided")$n
}


FifteenPercentDiff <- SummaryStatsByStudy$averageStudyF * 0.10
FirmFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  FirmFifteenPer[i] <- pwr.t.test(d=cohensD, 
                              power=0.8, sig.level=0.05, type="two.sample", 
                              alternative="two.sided")$n
}


FifteenPercentDiff <- SummaryStatsByStudy$averageStudyS * 0.10
SFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  SFifteenPer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


FifteenPercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.10
JFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  JFifteenPer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}

#################### WRITE DATA TO A CSV FILE ##############################

SimulatedN_B <- as.data.frame(cbind(BacterActualNeededN, BacterOnePer, 
                                    BacterFivePer, BacterTenPer, 
                                    BacterFifteenPer))
rownames(SimulatedN_B) <- Study
write.csv(SimulatedN_B, "results/tables/denovoSimulatedN_B.csv")

SimulatedN_F <- as.data.frame(cbind(
  FirmActualNeededN, FirmOnePer, FirmFivePer, FirmTenPer, FirmFifteenPer))
rownames(SimulatedN_F) <- Study
write.csv(SimulatedN_F, "results/tables/denovoSimulatedN_F.csv")


SimulatedN_S <- as.data.frame(cbind(SActualNeededN, SOnePer, SFivePer, 
                                    STenPer, SFifteenPer))
rownames(SimulatedN_S) <- Study
write.csv(SimulatedN_S, "results/tables/denovoSimulatedN_S.csv")

SimulatedN_J <- as.data.frame(cbind(JActualNeededN, JOnePer, JFivePer, 
                                    JTenPer, JFifteenPer))
rownames(SimulatedN_J) <- Study
write.csv(SimulatedN_J, "results/tables/denovoSimulatedN_J.csv")






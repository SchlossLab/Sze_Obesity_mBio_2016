## Power Calculations Variations: Estimating needed n
## Obesesity and the bacterial microbiome
## Marc Sze
## March 19, 2016

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

HActualNeededN <- rep(0, 8)

for(i in 1:8){
  cohensD <- (SummaryStatsByObeseGroup$MeanNonObeseH[i] -  
                SummaryStatsByObeseGroup$MeanObeseH[i]) / 
    SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  HActualNeededN[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}

BFActualNeededN <- rep(0, 8)

for(i in 1:8){
  cohensD <- (SummaryStatsByObeseGroup$MeanNonObeseBF[i] -  
                SummaryStatsByObeseGroup$MeanObeseBF[i]) / 
    SummaryStatsByStudy$sdBF[i]  
  # This needs to change based on study
  
  BFActualNeededN[i] <- pwr.t.test(d=cohensD, 
                                 power=0.8, sig.level=0.05, type="two.sample", 
                                 alternative="two.sided")$n
}

########## ONE PERCENT #####################################################

OnePercentDiff <- SummaryStatsByStudy$averageStudyH * 0.01
HOnePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  HOnePer[i] <- pwr.t.test(d=cohensD, 
                     power=0.8, sig.level=0.05, type="two.sample", 
                     alternative="two.sided")$n
}


OnePercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.01
BFratioOnePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratioOnePer[i] <- pwr.t.test(d=cohensD, 
                      power=0.8, sig.level=0.05, type="two.sample", 
                      alternative="two.sided")$n
}


########## FIVE PERCENT #####################################################

FivePercentDiff <- SummaryStatsByStudy$averageStudyH * 0.05
HFivePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  HFivePer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


FivePercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.05
BFratioFivePer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratioFivePer[i] <- pwr.t.test(d=cohensD, 
                                 power=0.8, sig.level=0.05, type="two.sample", 
                                 alternative="two.sided")$n
}


########## TEN PERCENT #####################################################

TenPercentDiff <- SummaryStatsByStudy$averageStudyH * 0.10
HTenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  HTenPer[i] <- pwr.t.test(d=cohensD, 
                            power=0.8, sig.level=0.05, type="two.sample", 
                            alternative="two.sided")$n
}


TenPercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.10
BFratioTenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratioTenPer[i] <- pwr.t.test(d=cohensD, 
                                  power=0.8, sig.level=0.05, type="two.sample", 
                                  alternative="two.sided")$n
}


########## FIFTEEN PERCENT #####################################################

FifteenPercentDiff <- SummaryStatsByStudy$averageStudyH * 0.10
HFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  HFifteenPer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


FifteenPercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.10
BFratioFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratioFifteenPer[i] <- pwr.t.test(d=cohensD, 
                                 power=0.8, sig.level=0.05, type="two.sample", 
                                 alternative="two.sided")$n
}


SimulatedN_H <- as.data.frame(cbind(HActualNeededN, 
                                    HOnePer, HFivePer, HTenPer, HFifteenPer))
rownames(SimulatedN_H) <- Study
write.csv(SimulatedN_H, "results/tables/denovoSimulatedN_H.csv")

SimulatedN_BF <- as.data.frame(cbind(
  BFActualNeededN, 
  BFratioOnePer, BFratioFivePer, BFratioTenPer, BFratioFifteenPer))
rownames(SimulatedN_BF) <- Study
write.csv(SimulatedN_BF, "results/tables/denovoSimulatedN_BF.csv")


####################### Estimation of n for percent differences
####################### Between two proportions (ACTUAL)

# Shannon Diversity
# Calculate the Effect size (h)

effectSize <- rep(0, 8)
for(i in 1:8){
  TotalTreatment <- ShannonRRTable$tposH[i] + ShannonRRTable$tnegH[i]
  ControlTreatment <- ShannonRRTable$cposH[i] + ShannonRRTable$cnegH[i]
  
  effectSize[i] <- ES.h(ShannonRRTable$tposH[i]/TotalTreatment, 
                        ShannonRRTable$cposH[i]/ControlTreatment)
}

# Calculate Simulated N for each data set 
RRActSimNH <- rep(0, 8)
for(i in 1:8){
  try(
  RRActSimNH[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                               power=0.8, alternative = "two.sided")$n)
}

# BF Ratio
# Calculate the Effect size (h)

effectSize <- rep(0, 8)
for(i in 1:8){
  TotalTreatment <- BFRatioRRTable$tposBF[i] + BFRatioRRTable$tnegBF[i]
  ControlTreatment <- BFRatioRRTable$cposBF[i] + BFRatioRRTable$cnegBF[i]
  
  effectSize[i] <- ES.h(BFRatioRRTable$tposBF[i]/TotalTreatment, 
                        BFRatioRRTable$cposBF[i]/ControlTreatment)
}

# Calculate Simulated N for each data set 
RRActSimNBF <- rep(0, 8)
for(i in 1:8){
  try(
    RRActSimNBF[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                              power=0.8, alternative = "two.sided")$n)
}

RRSimNNeeded <- cbind(RRActSimNH, RRActSimNBF)
rownames(RRSimNNeeded) <- Study
write.csv(RRSimNNeeded, "results/tables/denovoRRSimNNeeded.csv")

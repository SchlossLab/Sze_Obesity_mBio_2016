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
  cohensD <- abs((SummaryStatsByObeseGroup$MeanNonObeseBF[i] -  
                SummaryStatsByObeseGroup$MeanObeseBF[i])) / 
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

FifteenPercentDiff <- SummaryStatsByStudy$averageStudyH * 0.15
HFifteenPer <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  HFifteenPer[i] <- pwr.t.test(d=cohensD, 
                           power=0.8, sig.level=0.05, type="two.sample", 
                           alternative="two.sided")$n
}


FifteenPercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.15
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

# First need to calculate the base proportions
TreatmentProportionH <- (ShannonRRTable$tposH)/(ShannonRRTable$tposH + 
                                                 ShannonRRTable$tnegH)

ControlProportionH <- (ShannonRRTable$cposH)/(ShannonRRTable$cposH + 
                                              ShannonRRTable$cnegH)

TreatmentProportionBF <- (BFRatioRRTable$tposBF)/(BFRatioRRTable$tposBF + 
                                                  BFRatioRRTable$tnegBF)

ControlProportionBF <- (BFRatioRRTable$cposBF)/(BFRatioRRTable$cposBF + 
                                               BFRatioRRTable$cnegBF)


# Shannon Diversity
# Calculate the Effect size (h)
# Calculate Simulated N for each data set

effectSize <- rep(0, 8)
RRActSimNH <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionH[i], 
                        TreatmentProportionH[i])
  try(
  RRActSimNH[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                               power=0.8, alternative = "two.sided")$n)
}

# BF Ratio
# Calculate the Effect size (h)
# Calculate Simulated N for each data set 

effectSize <- rep(0, 8)
RRActSimNBF <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionBF[i], 
                        TreatmentProportionBF[i])
  try(
    RRActSimNBF[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                              power=0.8, alternative = "two.sided")$n)
}


####################### Between two proportions 1% difference

# Shannon Diversity
# Calculate the Effect size (h)
# Calculate Simulated N for each data set

effectSize <- rep(0, 8)
RROnePercentNH <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionH[i], 
                        ControlProportionH[i] + 0.01)
  try(
    RROnePercentNH[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                 power=0.8, alternative = "two.sided")$n)
}

# BF Ratio
# Calculate the Effect size (h)
# Calculate Simulated N for each data set 

effectSize <- rep(0, 8)
RROnePercentNBF <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionBF[i], 
                        ControlProportionBF[i] + 0.01)
  try(
    RROnePercentNBF[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                  power=0.8, alternative = "two.sided")$n)
}




####################### Between two proportions 5% difference

# Shannon Diversity
# Calculate the Effect size (h)
# Calculate Simulated N for each data set

effectSize <- rep(0, 8)
RRFivePercentNH <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionH[i], 
                        ControlProportionH[i] + 0.05)
  try(
    RRFivePercentNH[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                 power=0.8, alternative = "two.sided")$n)
}

# BF Ratio
# Calculate the Effect size (h)
# Calculate Simulated N for each data set 

effectSize <- rep(0, 8)
RRFivePercentNBF <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionBF[i], 
                        ControlProportionBF[i] + 0.05)
  try(
    RRFivePercentNBF[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                  power=0.8, alternative = "two.sided")$n)
}


####################### Between two proportions 10% difference

# Shannon Diversity
# Calculate the Effect size (h)
# Calculate Simulated N for each data set

effectSize <- rep(0, 8)
RRTenPercentNH <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionH[i], 
                        ControlProportionH[i] + 0.10)
  try(
    RRTenPercentNH[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                      power=0.8, alternative = "two.sided")$n)
}

# BF Ratio
# Calculate the Effect size (h)
# Calculate Simulated N for each data set 

effectSize <- rep(0, 8)
RRTenPercentNBF <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionBF[i], 
                        ControlProportionBF[i] + 0.10)
  try(
    RRTenPercentNBF[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                       power=0.8, alternative = "two.sided")$n)
}


####################### Between two proportions 15% difference

# Shannon Diversity
# Calculate the Effect size (h)
# Calculate Simulated N for each data set

effectSize <- rep(0, 8)
RRFifteenPercentNH <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionH[i], 
                        ControlProportionH[i] + 0.15)
  try(
    RRFifteenPercentNH[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                     power=0.8, alternative = "two.sided")$n)
}

# BF Ratio
# Calculate the Effect size (h)
# Calculate Simulated N for each data set 

effectSize <- rep(0, 8)
RRFifteenPercentNBF <- rep(0, 8)
for(i in 1:8){
  effectSize[i] <- ES.h(ControlProportionBF[i], 
                        ControlProportionBF[i] + 0.15)
  try(
    RRFifteenPercentNBF[i] <- pwr.2p.test(h=abs(effectSize[i]), sig.level=0.05, 
                                      power=0.8, alternative = "two.sided")$n)
}


############ Write Data Tables to file ###################################

RRSimNNeededH <- cbind(RRActSimNH, RROnePercentNH, RRFivePercentNH, 
                       RRTenPercentNH, RRFifteenPercentNH)
rownames(RRSimNNeededH) <- Study
write.csv(RRSimNNeededH, "results/tables/denovoRRSimNNeededH.csv")

RRSimNNeededBF <- cbind(RRActSimNBF, RROnePercentNBF, RRFivePercentNBF, 
                       RRTenPercentNBF, RRFifteenPercentNBF)
rownames(RRSimNNeededBF) <- Study
write.csv(RRSimNNeededBF, "results/tables/denovoRRSimNNeededBF.csv")






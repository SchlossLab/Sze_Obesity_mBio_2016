## Power Calculations Variations
## Obesesity and the bacterial microbiome
## Marc Sze
## March 18, 2016

library(statmod)
library(pwr)
source("code/UsedFunctions.R")

ShannonRRTable <- read.csv("results/tables/denovoShannonRRTable.csv")
BFRatioRRTable <- read.csv("results/tables/denovoBFRatioRRTable.csv")
PowerTable <- read.csv("results/tables/denovoPowerTable.csv")
SummaryStatsByObeseGroup <- read.csv("results/tables/denovoSSbyObeseGroup.csv")
SummaryStatsByStudy <- read.csv("results/tables/denovoSSbyStudy.csv")

Study <- as.character(PowerTable$Study)



########## ONE PERCENT #####################################################

OnePercentDiff <- SummaryStatsByStudy$averageStudyH * 0.01
H <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdH[i] 
    # This needs to change based on study
  
  H[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
                       )$power 
  }


OnePercentDiff <- SummaryStatsByStudy$averageStudyS * 0.01
S <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  S[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

OnePercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.01
J <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  J[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


OnePercentDiff <- (SummaryStatsByStudy$averageStudyB) * 0.01
bacter <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  bacter[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

OnePercentDiff <- (SummaryStatsByStudy$averageStudyF) * 0.01
firm <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  firm[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

OnePercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.01
BFratio <- rep(0, 8)

for(i in 1:8){
  cohensD <- OnePercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratio[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

ModelledPower1PerDiff <- as.data.frame(cbind(bacter, firm, BFratio, H, S, J))
rownames(ModelledPower1PerDiff) <- Study

########## FIVE PERCENT #####################################################

FivePercentDiff <- SummaryStatsByStudy$averageStudyH * 0.05
H <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  H[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


FivePercentDiff <- SummaryStatsByStudy$averageStudyS * 0.05
S <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  S[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

FivePercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.05
J <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  J[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


FivePercentDiff <- (SummaryStatsByStudy$averageStudyB) * 0.05
bacter <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  bacter[i] <- pwr.t2n.test(d=cohensD,
                            n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                            n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

FivePercentDiff <- (SummaryStatsByStudy$averageStudyF) * 0.05
firm <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  firm[i] <- pwr.t2n.test(d=cohensD,
                          n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                          n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

FivePercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.05
BFratio <- rep(0, 8)

for(i in 1:8){
  cohensD <- FivePercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratio[i] <- pwr.t2n.test(d=cohensD,
                             n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                             n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

ModelledPower5PerDiff <- as.data.frame(cbind(bacter, firm, BFratio, H, S, J))
rownames(ModelledPower5PerDiff) <- Study

########## TEN PERCENT #####################################################

TenPercentDiff <- SummaryStatsByStudy$averageStudyH * 0.1

H <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  H[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


TenPercentDiff <- SummaryStatsByStudy$averageStudyS * 0.1
S <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  S[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

TenPercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.1
J <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  J[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


TenPercentDiff <- (SummaryStatsByStudy$averageStudyB) * 0.1
bacter <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  bacter[i] <- pwr.t2n.test(d=cohensD,
                            n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                            n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

TenPercentDiff <- (SummaryStatsByStudy$averageStudyF) * 0.1
firm <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  firm[i] <- pwr.t2n.test(d=cohensD,
                          n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                          n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

TenPercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.1
BFratio <- rep(0, 8)

for(i in 1:8){
  cohensD <- TenPercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratio[i] <- pwr.t2n.test(d=cohensD,
                             n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                             n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

ModelledPower10PerDiff <- as.data.frame(cbind(bacter, firm, BFratio, H, S, J))
rownames(ModelledPower10PerDiff) <- Study

########## FIFTEEN PERCENT ##################################################

FifteenPercentDiff <- SummaryStatsByStudy$averageStudyH * 0.15

H <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdH[i] 
  # This needs to change based on study
  
  H[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


FifteenPercentDiff <- SummaryStatsByStudy$averageStudyS * 0.15
S <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdS[i] 
  # This needs to change based on study
  
  S[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

FifteenPercentDiff <- SummaryStatsByStudy$averageStudyJ * 0.15
J <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdJ[i] 
  # This needs to change based on study
  
  J[i] <- pwr.t2n.test(d=cohensD,
                       n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                       n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}


FifteenPercentDiff <- (SummaryStatsByStudy$averageStudyB) * 0.15
bacter <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdB[i] 
  # This needs to change based on study
  
  bacter[i] <- pwr.t2n.test(d=cohensD,
                            n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                            n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

FifteenPercentDiff <- (SummaryStatsByStudy$averageStudyF) * 0.15
firm <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdF[i] 
  # This needs to change based on study
  
  firm[i] <- pwr.t2n.test(d=cohensD,
                          n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                          n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

FifteenPercentDiff <- SummaryStatsByStudy$averageStudyBF * 0.15
BFratio <- rep(0, 8)

for(i in 1:8){
  cohensD <- FifteenPercentDiff[i] / SummaryStatsByStudy$sdBF[i] 
  # This needs to change based on study
  
  BFratio[i] <- pwr.t2n.test(d=cohensD,
                             n1=ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                             n2=ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i]
  )$power 
}

ModelledPower15PerDiff <- as.data.frame(cbind(bacter, firm, BFratio, H, S, J))
rownames(ModelledPower15PerDiff) <- Study

########## Actual Difference ################################################

H <- rep(0, 8)
for(i in 1:8){
  H[i] <- ((abs(SummaryStatsByObeseGroup$MeanNonObeseH[i] - 
                  SummaryStatsByObeseGroup$MeanObeseH[i])) / 
             SummaryStatsByStudy$averageStudyH[i]) * 100
}

S <- rep(0, 8)
for(i in 1:8){
  S[i] <- ((abs(SummaryStatsByObeseGroup$MeanNonObeseS[i] - 
                  SummaryStatsByObeseGroup$MeanObeseS[i])) / 
             SummaryStatsByStudy$averageStudyS[i]) * 100
}


J <- rep(0, 8)
for(i in 1:8){
  J[i] <- ((abs(SummaryStatsByObeseGroup$MeanNonObeseJ[i] - 
                  SummaryStatsByObeseGroup$MeanObeseJ[i])) / 
             SummaryStatsByStudy$averageStudyJ[i]) * 100
}


bacter <- rep(0, 8)
for(i in 1:8){
  bacter[i] <- ((abs(SummaryStatsByObeseGroup$MeanNonObeseB[i] - 
                       SummaryStatsByObeseGroup$MeanObeseB[i])) / 
                  SummaryStatsByStudy$averageStudyB[i]) * 100
}

firm <- rep(0, 8)
for(i in 1:8){
  firm[i] <- ((abs(SummaryStatsByObeseGroup$MeanNonObeseF[i] - 
                       SummaryStatsByObeseGroup$MeanObeseF[i])) / 
                SummaryStatsByStudy$averageStudyF[i]) * 100
}

BFratio <- rep(0, 8)
for(i in 1:8){
  BFratio[i] <- ((abs(SummaryStatsByObeseGroup$MeanNonObeseBF[i] - 
                             SummaryStatsByObeseGroup$MeanObeseBF[i])) / 
                        SummaryStatsByStudy$averageStudyBF[i]) * 100
}

ActualDiff <- as.data.frame(cbind(bacter, firm, BFratio, H, S, J))
rownames(ActualDiff) <- Study

#######################################################################################

                               
### Writing out Respective Tables 

write.csv(ModelledPower1PerDiff, "results/tables/denovoModelled1PerDiff.csv")
write.csv(ModelledPower5PerDiff, "results/tables/denovoModelled5PerDiff.csv")
write.csv(ModelledPower10PerDiff, "results/tables/denovoModelled10PerDiff.csv")
write.csv(ModelledPower15PerDiff, "results/tables/denovoModelled15PerDiff.csv")
write.csv(ActualDiff, "results/tables/denovoActualDiff.csv")


###############################################################################
########## Model power for Fisher

## Used the mean prevalence of obesity as a basis point for varying proportions.
## So average obesity for all groups combined is 38%.  So this value was always
## used for the control group and the treatment group would vary by:
## 1%, 5%, 10%, 15%

#One percent

OnePercentHRR <- rep(0, 8)
for(i in 1:8){
  OnePercentHRR[i] <- power.fisher.test(0.38, 0.39, 
                              ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                              ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}

OnePercentBFRR <- rep(0, 8)
for(i in 1:8){
  OnePercentBFRR[i] <- power.fisher.test(0.38, 0.39, 
                              BFRatioRRTable$tnegBF[i] + BFRatioRRTable$cnegBF[i], 
                              BFRatioRRTable$tposBF[i] + BFRatioRRTable$cposBF[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}


#Five percent

FivePercentHRR <- rep(0, 8)
for(i in 1:8){
  FivePercentHRR[i] <- power.fisher.test(0.38, 0.43, 
                              ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                              ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}

FivePercentBFRR <- rep(0, 8)
for(i in 1:8){
  FivePercentBFRR[i] <- power.fisher.test(0.38, 0.43, 
                              BFRatioRRTable$tnegBF[i] + BFRatioRRTable$cnegBF[i], 
                              BFRatioRRTable$tposBF[i] + BFRatioRRTable$cposBF[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}


#Ten percent

TenPercentHRR <- rep(0, 8)
for(i in 1:8){
  TenPercentHRR[i] <- power.fisher.test(0.38, 0.48, 
                              ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                              ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}

TenPercentBFRR <- rep(0, 8)
for(i in 1:8){
  TenPercentBFRR[i] <- power.fisher.test(0.38, 0.48, 
                              BFRatioRRTable$tnegBF[i] + BFRatioRRTable$cnegBF[i], 
                              BFRatioRRTable$tposBF[i] + BFRatioRRTable$cposBF[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}


#Fifteen percent

FifteenPercentHRR <- rep(0, 8)
for(i in 1:8){
  FifteenPercentHRR[i] <- power.fisher.test(0.38, 0.53, 
                              ShannonRRTable$tnegH[i] + ShannonRRTable$cnegH[i], 
                              ShannonRRTable$tposH[i] + ShannonRRTable$cposH[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}

FifteenPercentBFRR <- rep(0, 8)
for(i in 1:8){
  FifteenPercentBFRR[i] <- power.fisher.test(0.38, 0.53, 
                              BFRatioRRTable$tnegBF[i] + BFRatioRRTable$cnegBF[i], 
                              BFRatioRRTable$tposBF[i] + BFRatioRRTable$cposBF[i], 
                              alpha=0.05, nsim=1000, alternative="two.sided")
}


##Actual Percent differences

TreatPropH <- ShannonRRTable$tposH / (ShannonRRTable$tposH + ShannonRRTable$tnegH)
ContPropH <- ShannonRRTable$cposH / (ShannonRRTable$cposH + ShannonRRTable$cnegH)
TreatPropBF <- BFRatioRRTable$tposBF / (BFRatioRRTable$tposBF + BFRatioRRTable$tnegBF)
ContPropBF <- BFRatioRRTable$cposBF / (BFRatioRRTable$cposBF + BFRatioRRTable$cnegBF)

ActualDiffHRR <- abs(TreatPropH - ContPropH) 
ActualDiffBFRR <- abs(TreatPropBF - ContPropBF)

PowerModelHRR <- as.data.frame(cbind(
  ActualDiffHRR, OnePercentHRR, FivePercentHRR, TenPercentHRR, FifteenPercentHRR
))
rownames(PowerModelHRR) <- Study
write.csv(PowerModelHRR, "results/tables/denovoPowerModelHRR.csv")

PowerModelBFRR <- as.data.frame(cbind(
  ActualDiffBFRR, OnePercentBFRR, FivePercentBFRR, TenPercentBFRR, FifteenPercentBFRR
))
rownames(PowerModelBFRR) <- Study

write.csv(PowerModelBFRR, "results/tables/denovoPowerModelBFRR.csv")




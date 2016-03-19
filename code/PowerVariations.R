## Power Calculations Variations
## Obesesity and the bacterial microbiome
## Marc Sze
## March 18, 2016

library(statmod)
source("code/UsedFunctions.R")

ShannonRRTable <- read.csv("results/tables/denovoShannonRRTable.csv")
BFRatioRRTable <- read.csv("results/tables/denovoBFRatioRRTable.csv")
PowerTable <- read.csv("results/tables/denovoPowerTable.csv")
Study <- as.character(PowerTable$Study)

FifteenPerDiff <- mean(PowerTable$averageStudyH) * 0.15
TenPercentDiff <- mean(PowerTable$averageStudyH) * 0.10 
FivePercentDiff <- mean(PowerTable$averageStudyH) * 0.05
OnePercentDiff <- mean(PowerTable$averageStudyH) * 0.01

## Generate Estimates for Shannon Diversity

sdH <- PowerTable$sdH

#Baxter
OneDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyH[1], 
                  PowerTable$averageStudyH[1] - OnePercentDiff, 
                  ShannonRRTable$tnegH[1] + ShannonRRTable$cnegH[1], 
                  ShannonRRTable$tposH[1] + ShannonRRTable$cposH[1], 
                  sdH[1], iters=1000)

FiveDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyH[1], 
                   PowerTable$averageStudyH[1] - FivePercentDiff, 
                   ShannonRRTable$tnegH[1] + ShannonRRTable$cnegH[1], 
                   ShannonRRTable$tposH[1] + ShannonRRTable$cposH[1], 
                   sdH[1], iters=1000)

TenDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyH[1], 
                  PowerTable$averageStudyH[1] - TenPercentDiff, 
                  ShannonRRTable$tnegH[1] + ShannonRRTable$cnegH[1], 
                  ShannonRRTable$tposH[1] + ShannonRRTable$cposH[1], 
                  sdH[1], iters=1000)

FifteenDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyH[1], 
                  PowerTable$averageStudyH[1] - FifteenPerDiff, 
                  ShannonRRTable$tnegH[1] + ShannonRRTable$cnegH[1], 
                  ShannonRRTable$tposH[1] + ShannonRRTable$cposH[1], 
                  sdH[1], iters=1000)

ActualDiff <- ((abs(PowerTable$MeanNonObeseH[1] - PowerTable$MeanObeseH[1])) / 
  PowerTable$averageStudyH[1]) * 100


#Ross
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[2], 
  PowerTable$averageStudyH[2] - OnePercentDiff, 
  ShannonRRTable$tnegH[2] + ShannonRRTable$cnegH[2], 
  ShannonRRTable$tposH[2] + ShannonRRTable$cposH[2], 
  sdH[2], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[2], 
  PowerTable$averageStudyH[2] - FivePercentDiff, 
  ShannonRRTable$tnegH[2] + ShannonRRTable$cnegH[2], 
  ShannonRRTable$tposH[2] + ShannonRRTable$cposH[2], 
  sdH[2], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[2], 
  PowerTable$averageStudyH[2] - TenPercentDiff, 
  ShannonRRTable$tnegH[2] + ShannonRRTable$cnegH[2], 
  ShannonRRTable$tposH[2] + ShannonRRTable$cposH[2], 
  sdH[2], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[2], 
  PowerTable$averageStudyH[2] - FifteenPerDiff, 
  ShannonRRTable$tnegH[2] + ShannonRRTable$cnegH[2], 
  ShannonRRTable$tposH[2] + ShannonRRTable$cposH[2], 
  sdH[2], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[2] - PowerTable$MeanObeseH[2])) 
                             / PowerTable$averageStudyH[2]) * 100)

#Goodrich
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[3], 
  PowerTable$averageStudyH[3] - OnePercentDiff, 
  ShannonRRTable$tnegH[3] + ShannonRRTable$cnegH[3], 
  ShannonRRTable$tposH[3] + ShannonRRTable$cposH[3], 
  sdH[3], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[3], 
  PowerTable$averageStudyH[3] - FivePercentDiff, 
  ShannonRRTable$tnegH[3] + ShannonRRTable$cnegH[3], 
  ShannonRRTable$tposH[3] + ShannonRRTable$cposH[3], 
  sdH[3], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[3], 
  PowerTable$averageStudyH[3] - TenPercentDiff, 
  ShannonRRTable$tnegH[3] + ShannonRRTable$cnegH[3], 
  ShannonRRTable$tposH[3] + ShannonRRTable$cposH[3], 
  sdH[3], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[3], 
  PowerTable$averageStudyH[3] - FifteenPerDiff, 
  ShannonRRTable$tnegH[3] + ShannonRRTable$cnegH[3], 
  ShannonRRTable$tposH[3] + ShannonRRTable$cposH[3], 
  sdH[3], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[3] - PowerTable$MeanObeseH[3])) 
                             / PowerTable$averageStudyH[3]) * 100)

#Escobar
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[4], 
  PowerTable$averageStudyH[4] - OnePercentDiff, 
  ShannonRRTable$tnegH[4] + ShannonRRTable$cnegH[4], 
  ShannonRRTable$tposH[4] + ShannonRRTable$cposH[4], 
  sdH[4], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[4], 
  PowerTable$averageStudyH[4] - FivePercentDiff, 
  ShannonRRTable$tnegH[4] + ShannonRRTable$cnegH[4], 
  ShannonRRTable$tposH[4] + ShannonRRTable$cposH[4], 
  sdH[4], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[4], 
  PowerTable$averageStudyH[4] - TenPercentDiff, 
  ShannonRRTable$tnegH[4] + ShannonRRTable$cnegH[4], 
  ShannonRRTable$tposH[4] + ShannonRRTable$cposH[4], 
  sdH[4], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[4], 
  PowerTable$averageStudyH[4] - FifteenPerDiff, 
  ShannonRRTable$tnegH[4] + ShannonRRTable$cnegH[4], 
  ShannonRRTable$tposH[4] + ShannonRRTable$cposH[4], 
  sdH[4], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[4] - PowerTable$MeanObeseH[4])) 
                             / PowerTable$averageStudyH[4]) * 100)

#Zupancic
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[5], 
  PowerTable$averageStudyH[5] - OnePercentDiff, 
  ShannonRRTable$tnegH[5] + ShannonRRTable$cnegH[5], 
  ShannonRRTable$tposH[5] + ShannonRRTable$cposH[5], 
  sdH[5], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[5], 
  PowerTable$averageStudyH[5] - FivePercentDiff, 
  ShannonRRTable$tnegH[5] + ShannonRRTable$cnegH[5], 
  ShannonRRTable$tposH[5] + ShannonRRTable$cposH[5], 
  sdH[5], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[5], 
  PowerTable$averageStudyH[5] - TenPercentDiff, 
  ShannonRRTable$tnegH[5] + ShannonRRTable$cnegH[5], 
  ShannonRRTable$tposH[5] + ShannonRRTable$cposH[5], 
  sdH[5], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[5], 
  PowerTable$averageStudyH[5] - FifteenPerDiff, 
  ShannonRRTable$tnegH[5] + ShannonRRTable$cnegH[5], 
  ShannonRRTable$tposH[5] + ShannonRRTable$cposH[5], 
  sdH[5], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[5] - PowerTable$MeanObeseH[5])) 
                             / PowerTable$averageStudyH[5]) * 100)

#HMP
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[6], 
  PowerTable$averageStudyH[6] - OnePercentDiff, 
  ShannonRRTable$tnegH[6] + ShannonRRTable$cnegH[6], 
  ShannonRRTable$tposH[6] + ShannonRRTable$cposH[6], 
  sdH[6], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[6], 
  PowerTable$averageStudyH[6] - FivePercentDiff, 
  ShannonRRTable$tnegH[6] + ShannonRRTable$cnegH[6], 
  ShannonRRTable$tposH[6] + ShannonRRTable$cposH[6], 
  sdH[6], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[6], 
  PowerTable$averageStudyH[6] - TenPercentDiff, 
  ShannonRRTable$tnegH[6] + ShannonRRTable$cnegH[6], 
  ShannonRRTable$tposH[6] + ShannonRRTable$cposH[6], 
  sdH[6], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[6], 
  PowerTable$averageStudyH[6] - FifteenPerDiff, 
  ShannonRRTable$tnegH[6] + ShannonRRTable$cnegH[6], 
  ShannonRRTable$tposH[6] + ShannonRRTable$cposH[6], 
  sdH[6], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[6] - PowerTable$MeanObeseH[6])) 
                             / PowerTable$averageStudyH[6]) * 100)

#Wu
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[7], 
  PowerTable$averageStudyH[7] - OnePercentDiff, 
  ShannonRRTable$tnegH[7] + ShannonRRTable$cnegH[7], 
  ShannonRRTable$tposH[7] + ShannonRRTable$cposH[7], 
  sdH[7], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[7], 
  PowerTable$averageStudyH[7] - FivePercentDiff, 
  ShannonRRTable$tnegH[7] + ShannonRRTable$cnegH[7], 
  ShannonRRTable$tposH[7] + ShannonRRTable$cposH[7], 
  sdH[7], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[7], 
  PowerTable$averageStudyH[7] - TenPercentDiff, 
  ShannonRRTable$tnegH[7] + ShannonRRTable$cnegH[7], 
  ShannonRRTable$tposH[7] + ShannonRRTable$cposH[7], 
  sdH[7], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[7], 
  PowerTable$averageStudyH[7] - FifteenPerDiff, 
  ShannonRRTable$tnegH[7] + ShannonRRTable$cnegH[7], 
  ShannonRRTable$tposH[7] + ShannonRRTable$cposH[7], 
  sdH[7], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[7] - PowerTable$MeanObeseH[7])) 
                             / PowerTable$averageStudyH[7]) * 100)

#Turnbaugh
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[8], 
  PowerTable$averageStudyH[8] - OnePercentDiff, 
  ShannonRRTable$tnegH[8] + ShannonRRTable$cnegH[8], 
  ShannonRRTable$tposH[8] + ShannonRRTable$cposH[8], 
  sdH[8], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[8], 
  PowerTable$averageStudyH[8] - FivePercentDiff, 
  ShannonRRTable$tnegH[8] + ShannonRRTable$cnegH[8], 
  ShannonRRTable$tposH[8] + ShannonRRTable$cposH[8], 
  sdH[8], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[8], 
  PowerTable$averageStudyH[8] - TenPercentDiff, 
  ShannonRRTable$tnegH[8] + ShannonRRTable$cnegH[8], 
  ShannonRRTable$tposH[8] + ShannonRRTable$cposH[8], 
  sdH[8], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyH[8], 
  PowerTable$averageStudyH[8] - FifteenPerDiff, 
  ShannonRRTable$tnegH[8] + ShannonRRTable$cnegH[8], 
  ShannonRRTable$tposH[8] + ShannonRRTable$cposH[8], 
  sdH[8], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseH[8] - PowerTable$MeanObeseH[8])) 
                             / PowerTable$averageStudyH[8]) * 100)

### Putting it all together 

ModelledPowerH <- as.data.frame(
  cbind(OneDiffPower, FiveDiffPower, TenDiffPower, FifteenDiffPower, format(round(
    ActualDiff, 1), nsmall=1)))
colnames(ModelledPowerH) <- c("Pow1Diff", "Pow5Diff", "Pow10Diff", "Pow15Diff", "ActualPerDiff")
ModelledPowerH$Study <- Study


## Generate Estimates for BF Ratio

FifteenPerDiff <- mean(PowerTable$averageStudyBF) * 0.15
TenPercentDiff <- mean(PowerTable$averageStudyBF) * 0.10 
FivePercentDiff <- mean(PowerTable$averageStudyBF) * 0.05
OnePercentDiff <- mean(PowerTable$averageStudyBF) * 0.01



sdBF <- PowerTable$sdBF

#Baxter
OneDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyBF[1], 
                                   PowerTable$averageStudyBF[1] - OnePercentDiff, 
                                   BFRatioRRTable$tnegBF[1] + BFRatioRRTable$cnegBF[1], 
                                   BFRatioRRTable$tposBF[1] + BFRatioRRTable$cposBF[1], 
                                   sdH[1], iters=1000)

FiveDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyBF[1], 
                                    PowerTable$averageStudyBF[1] - FivePercentDiff, 
                                    BFRatioRRTable$tnegBF[1] + BFRatioRRTable$cnegBF[1], 
                                    BFRatioRRTable$tposBF[1] + BFRatioRRTable$cposBF[1], 
                                    sdH[1], iters=1000)

TenDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyBF[1], 
                                   PowerTable$averageStudyBF[1] - TenPercentDiff, 
                                   BFRatioRRTable$tnegBF[1] + BFRatioRRTable$cnegBF[1], 
                                   BFRatioRRTable$tposBF[1] + BFRatioRRTable$cposBF[1], 
                                   sdH[1], iters=1000)

FifteenDiffPower <- DefNonParaPowerSim(PowerTable$averageStudyBF[1], 
                                       PowerTable$averageStudyBF[1] - FifteenPerDiff, 
                                       BFRatioRRTable$tnegBF[1] + BFRatioRRTable$cnegBF[1], 
                                       BFRatioRRTable$tposBF[1] + BFRatioRRTable$cposBF[1], 
                                       sdH[1], iters=1000)

ActualDiff <- ((abs(PowerTable$MeanNonObeseBF[1] - PowerTable$MeanObeseBF[1])) / 
                 PowerTable$averageStudyBF[1]) * 100


#Ross
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[2], 
  PowerTable$averageStudyBF[2] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[2] + BFRatioRRTable$cnegBF[2], 
  BFRatioRRTable$tposBF[2] + BFRatioRRTable$cposBF[2], 
  sdH[2], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[2], 
  PowerTable$averageStudyBF[2] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[2] + BFRatioRRTable$cnegBF[2], 
  BFRatioRRTable$tposBF[2] + BFRatioRRTable$cposBF[2], 
  sdH[2], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[2], 
  PowerTable$averageStudyBF[2] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[2] + BFRatioRRTable$cnegBF[2], 
  BFRatioRRTable$tposBF[2] + BFRatioRRTable$cposBF[2], 
  sdH[2], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[2], 
  PowerTable$averageStudyBF[2] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[2] + BFRatioRRTable$cnegBF[2], 
  BFRatioRRTable$tposBF[2] + BFRatioRRTable$cposBF[2], 
  sdH[2], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[2] - PowerTable$MeanObeseBF[2])) / 
                               PowerTable$averageStudyBF[2]) * 100)

#Goodrich
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[3], 
  PowerTable$averageStudyBF[3] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[3] + BFRatioRRTable$cnegBF[3], 
  BFRatioRRTable$tposBF[3] + BFRatioRRTable$cposBF[3], 
  sdH[3], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[3], 
  PowerTable$averageStudyBF[3] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[3] + BFRatioRRTable$cnegBF[3], 
  BFRatioRRTable$tposBF[3] + BFRatioRRTable$cposBF[3], 
  sdH[3], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[3], 
  PowerTable$averageStudyBF[3] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[3] + BFRatioRRTable$cnegBF[3], 
  BFRatioRRTable$tposBF[3] + BFRatioRRTable$cposBF[3], 
  sdH[3], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[3], 
  PowerTable$averageStudyBF[3] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[3] + BFRatioRRTable$cnegBF[3], 
  BFRatioRRTable$tposBF[3] + BFRatioRRTable$cposBF[3], 
  sdH[3], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[3] - PowerTable$MeanObeseBF[3])) / 
                               PowerTable$averageStudyBF[3]) * 100)

#Escobar
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[4], 
  PowerTable$averageStudyBF[4] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[4] + BFRatioRRTable$cnegBF[4], 
  BFRatioRRTable$tposBF[4] + BFRatioRRTable$cposBF[4], 
  sdH[4], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[4], 
  PowerTable$averageStudyBF[4] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[4] + BFRatioRRTable$cnegBF[4], 
  BFRatioRRTable$tposBF[4] + BFRatioRRTable$cposBF[4], 
  sdH[4], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[4], 
  PowerTable$averageStudyBF[4] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[4] + BFRatioRRTable$cnegBF[4], 
  BFRatioRRTable$tposBF[4] + BFRatioRRTable$cposBF[4], 
  sdH[4], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[4], 
  PowerTable$averageStudyBF[4] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[4] + BFRatioRRTable$cnegBF[4], 
  BFRatioRRTable$tposBF[4] + BFRatioRRTable$cposBF[4], 
  sdH[4], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[4] - PowerTable$MeanObeseBF[4])) / 
                               PowerTable$averageStudyBF[4]) * 100)

#Zupancic
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[5], 
  PowerTable$averageStudyBF[5] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[5] + BFRatioRRTable$cnegBF[5], 
  BFRatioRRTable$tposBF[5] + BFRatioRRTable$cposBF[5], 
  sdH[5], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[5], 
  PowerTable$averageStudyBF[5] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[5] + BFRatioRRTable$cnegBF[5], 
  BFRatioRRTable$tposBF[5] + BFRatioRRTable$cposBF[5], 
  sdH[5], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[5], 
  PowerTable$averageStudyBF[5] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[5] + BFRatioRRTable$cnegBF[5], 
  BFRatioRRTable$tposBF[5] + BFRatioRRTable$cposBF[5], 
  sdH[5], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[5], 
  PowerTable$averageStudyBF[5] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[5] + BFRatioRRTable$cnegBF[5], 
  BFRatioRRTable$tposBF[5] + BFRatioRRTable$cposBF[5], 
  sdH[5], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[5] - PowerTable$MeanObeseBF[5])) / 
                               PowerTable$averageStudyBF[5]) * 100)

#HMP
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[6], 
  PowerTable$averageStudyBF[6] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[6] + BFRatioRRTable$cnegBF[6], 
  BFRatioRRTable$tposBF[6] + BFRatioRRTable$cposBF[6], 
  sdH[6], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[6], 
  PowerTable$averageStudyBF[6] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[6] + BFRatioRRTable$cnegBF[6], 
  BFRatioRRTable$tposBF[6] + BFRatioRRTable$cposBF[6], 
  sdH[6], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[6], 
  PowerTable$averageStudyBF[6] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[6] + BFRatioRRTable$cnegBF[6], 
  BFRatioRRTable$tposBF[6] + BFRatioRRTable$cposBF[6], 
  sdH[6], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[6], 
  PowerTable$averageStudyBF[6] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[6] + BFRatioRRTable$cnegBF[6], 
  BFRatioRRTable$tposBF[6] + BFRatioRRTable$cposBF[6], 
  sdH[6], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[6] - PowerTable$MeanObeseBF[6])) / 
                               PowerTable$averageStudyBF[6]) * 100)

#Wu
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[7], 
  PowerTable$averageStudyBF[7] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[7] + BFRatioRRTable$cnegBF[7], 
  BFRatioRRTable$tposBF[7] + BFRatioRRTable$cposBF[7], 
  sdH[7], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[7], 
  PowerTable$averageStudyBF[7] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[7] + BFRatioRRTable$cnegBF[7], 
  BFRatioRRTable$tposBF[7] + BFRatioRRTable$cposBF[7], 
  sdH[7], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[7], 
  PowerTable$averageStudyBF[7] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[7] + BFRatioRRTable$cnegBF[7], 
  BFRatioRRTable$tposBF[7] + BFRatioRRTable$cposBF[7], 
  sdH[7], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[7], 
  PowerTable$averageStudyBF[7] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[7] + BFRatioRRTable$cnegBF[7], 
  BFRatioRRTable$tposBF[7] + BFRatioRRTable$cposBF[7], 
  sdH[7], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[7] - PowerTable$MeanObeseBF[7])) / 
                               PowerTable$averageStudyBF[7]) * 100)

#Turnbaugh
OneDiffPower <- c(OneDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[8], 
  PowerTable$averageStudyBF[8] - OnePercentDiff, 
  BFRatioRRTable$tnegBF[8] + BFRatioRRTable$cnegBF[8], 
  BFRatioRRTable$tposBF[8] + BFRatioRRTable$cposBF[8], 
  sdH[8], iters=1000))

FiveDiffPower <- c(FiveDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[8], 
  PowerTable$averageStudyBF[8] - FivePercentDiff, 
  BFRatioRRTable$tnegBF[8] + BFRatioRRTable$cnegBF[8], 
  BFRatioRRTable$tposBF[8] + BFRatioRRTable$cposBF[8], 
  sdH[8], iters=1000))

TenDiffPower <- c(TenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[8], 
  PowerTable$averageStudyBF[8] - TenPercentDiff, 
  BFRatioRRTable$tnegBF[8] + BFRatioRRTable$cnegBF[8], 
  BFRatioRRTable$tposBF[8] + BFRatioRRTable$cposBF[8], 
  sdH[8], iters=1000))

FifteenDiffPower <- c(FifteenDiffPower, DefNonParaPowerSim(
  PowerTable$averageStudyBF[8], 
  PowerTable$averageStudyBF[8] - FifteenPerDiff, 
  BFRatioRRTable$tnegBF[8] + BFRatioRRTable$cnegBF[8], 
  BFRatioRRTable$tposBF[8] + BFRatioRRTable$cposBF[8], 
  sdH[8], iters=1000))

ActualDiff <- c(ActualDiff, ((abs(PowerTable$MeanNonObeseBF[8] - PowerTable$MeanObeseBF[8])) / 
                               PowerTable$averageStudyBF[8]) * 100)

### Putting it all together 

ModelledPowerBF <- as.data.frame(
  cbind(OneDiffPower, FiveDiffPower, TenDiffPower, FifteenDiffPower, format(round(
    ActualDiff, 1), nsmall=1)))
colnames(ModelledPowerBF) <- c("Pow1Diff", "Pow5Diff", "Pow10Diff", "Pow15Diff", "ActualPerDiff")
ModelledPowerBF$Study <- Study


#######################################################################################



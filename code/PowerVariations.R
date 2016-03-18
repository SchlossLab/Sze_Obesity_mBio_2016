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

## Generate Estimates for Shannon Diversity

sdH <- PowerTable$sdH

#Baxter
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

#Ross
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

#Goodrich
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


#Escobar
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

#Zupancic
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

#HMP
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

#Wu
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

#Turnbaugh
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

ModelledPower <- as.data.frame(
  cbind(FiveDiffPower, TenDiffPower, FifteenDiffPower))
colnames(ModelledPower) <- c("5% Mean Difference", "10% Mean Difference", 
                             "15% Mean Difference")
ModelledPower$Study <- Study

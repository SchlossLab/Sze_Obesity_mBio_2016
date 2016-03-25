## Make Figures for Manuscript
## Obesesity and the bacterial microbiome
## Marc Sze
## March 23, 2016

library(ggplot2)
library(AUCRF)
library(pROC)
source("code/StrategydenovoAnalysisMetaHitDel.R")

rm(alpha.test, AUCRFDataTable, BFRatio, BFRatioRRTable, data, combinedData, 
   demographicsTable, metadata, microbiome, overallPTable, phyla.table.rel.abund, 
   PowerTable, ShannonRRTable, SummaryStatsByStudy, SummaryStatsByObeseGroup, 
   testset, ancestry, averageStudyB, averageStudyBF, averageStudyF, 
   averageStudyH, averageStudyJ, averageStudyS, Bacter, bacter, baxAUC, 
   baxKopt, BaxterBFRR, BaxterHRR, BFratio, bmi, cnegBF, cnegH, cposBF, cposH, 
   escoAUC, EscobarBFRR, EscobarHRR, escoKopt, females, Firm, firm, goodAUC, 
   goodKopt, H, highBF, highH, HMPAUC, HMPBFRR, HMPHRR, HMPKopt, J, KoptAll, 
   lowBF, lowH, males, maxBMI, meanAge, meanBMI, MeanNonObeseB, MeanNonObeseBF, 
   MeanNonObeseF, MeanNonObeseH, MeanNonObeseJ, MeanNonObeseS, MeanObeseB, 
   MeanObeseBF, MeanObeseF, MeanObeseH, MeanObeseJ, MeanObeseS, microb.rownames, 
   minBMI, OOBAUCAll, rossAUC, RossBFRR, RossHRR, rossKopt, RRBF, RRH, S, SDAge, 
   sdB, sdBF, SDBMI, sdF, sdH, sdJ, SDNonObeseB, SDNonObeseBF, SDNonObeseF, 
   SDNonObeseH, SDNonObeseJ, SDNonObeseS, SDObeseB, SDObeseBF, SDObeseF, 
   SDObeseH, SDObeseJ, SDObeseS, sdS, StudyPowerBF, StudyPowerBFRR, 
   StudyPowerH, StudyPowerHRR, temporary, tnegBF, tnegH, totalN, tposBF, tposH, 
   TurnAUC, TurnbaughBFRR, TurnbaughHRR, TurnKopt, WuAUC, WuBFRR, WuHRR, WuKopt, 
   ZupancicBFRR, ZupancicHRR, zupAUC, zupKopt, GoodrichHRR, GoodrichBFRR)

# Make sensitivity and specificity tables for each data set

##### Baxter
test <- predict(baxterAUCfit$RFopt, type='prob')[, 2]
baxter <- roc(demographics$obese~test)
BaxterROCData <- as.data.frame(
  cbind(baxter$sensitivities, baxter$specificities))
BaxterROCData$Study <- 'Baxter'

##### Ross
RossData <- AddBMIClass(RossData, "BMI", numbers=TRUE)
test <- predict(rossAUCFit$RFopt, type='prob')[, 2]
ross <- roc(RossData$obese~test)
RossROCData <- as.data.frame(
  cbind(ross$sensitivities, ross$specificities))
RossROCData$Study <- 'Ross'

##### Goodrich
GoodrichData <- AddBMIClass(GoodrichData, "BMI", numbers=TRUE)
test <- predict(goodrichAUCFit$RFopt, type='prob')[, 2]
goodrich <- roc(GoodrichData$obese~test)
GoodrichROCData <- as.data.frame(
  cbind(goodrich$sensitivities, goodrich$specificities))
GoodrichROCData$Study <- 'Goodrich'

##### Escobar
EscobarData <- AddBMIClass(EscobarData, "BMI", numbers=TRUE)
test <- predict(escobarAUCFit$RFopt, type='prob')[, 2]
escobar <- roc(EscobarData$obese~test)
EscobarROCData <- as.data.frame(
  cbind(escobar$sensitivities, escobar$specificities))
EscobarROCData$Study <- 'Escobar'

##### Zupancic
ZupancicData <- AddBMIClass(ZupancicData, "BMI", numbers=TRUE)
test <- predict(zupancicAUCFit$RFopt, type='prob')[, 2]
zupancic <- roc(ZupancicData$obese~test)
ZupancicROCData <- as.data.frame(
  cbind(zupancic$sensitivities, zupancic$specificities))
ZupancicROCData$Study <- 'Zupancic'

##### HMP
HMPData <- AddBMIClass(HMPData, "BMI", numbers=TRUE)
test <- predict(HMPAUCFit$RFopt, type='prob')[, 2]
hmp <- roc(HMPData$obese~test)
HMPROCData <- as.data.frame(
  cbind(hmp$sensitivities, hmp$specificities))
HMPROCData$Study <- 'HMP'

##### Wu
WuData <- AddBMIClass(WuData, "BMI", numbers=TRUE)
test <- predict(WuAUCFit$RFopt, type='prob')[, 2]
wu <- roc(WuData$obese~test)
WuROCData <- as.data.frame(
  cbind(wu$sensitivities, wu$specificities))
WuROCData$Study <- 'Wu'

##### Turnbaugh
TurnbaughData <- AddBMIClass(TurnbaughData, "BMICat", numbers=FALSE)
test <- predict(TurnbaughAUCFit$RFopt, type='prob')[, 2]
turnbaugh <- roc(TurnbaughData$obese~test)
TurnbaughROCData <- as.data.frame(
  cbind(turnbaugh$sensitivities, turnbaugh$specificities))
TurnbaughROCData$Study <- 'Turnbaugh'

combinedData <- as.data.frame(rbind(
  BaxterROCData, RossROCData, GoodrichROCData, 
                      EscobarROCData, ZupancicROCData, HMPROCData, WuROCData, 
  TurnbaughROCData), stringsAsFactors = FALSE)
colnames(combinedData) <- c("Sensitivity", "Specificity", "Study")
combinedData$Study <- as.factor(combinedData$Study)

write.csv(combinedData, "results/tables/ROCcombinedData.csv")

diag <- data.frame(x=seq(1, 0), y = seq(0, 1))
diag$Study <- "Control"

rocPlot <- ggplot(data=combinedData, 
                  aes(x = Sensitivity, y = Specificity))
rocPlot + geom_line(aes(group=Study, color=Study)) + 
  scale_colour_brewer(palette="Set1") + 
  scale_x_continuous(trans = "reverse") + 
  geom_line(data=diag, aes(x=x, y=y), linetype="dashed") + 
  theme(axis.line=element_line(colour="black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(), 
       legend.key = element_rect(fill = "white"))


# Compare each curve to one another for any difference between them

# Baxter first

baxterTest <- rep(0, 8)
baxterTest[1] <- NA
baxterTest[2] <- roc.test(baxter, escobar)$p.value
baxterTest[3] <- roc.test(baxter, goodrich)$p.value
baxterTest[4] <- roc.test(baxter, hmp)$p.value
baxterTest[5] <- roc.test(baxter, ross)$p.value
baxterTest[6] <- roc.test(baxter, turnbaugh)$p.value
baxterTest[7] <- roc.test(baxter, wu)$p.value
baxterTest[8] <- roc.test(baxter, zupancic)$p.value

# Escobar Second

escobarTest <- rep(0, 8)
escobarTest[1] <- NA
escobarTest[2] <- NA
escobarTest[3] <- roc.test(escobar, goodrich)$p.value
escobarTest[4] <- roc.test(escobar, hmp)$p.value
escobarTest[5] <- roc.test(escobar, ross)$p.value
escobarTest[6] <- roc.test(escobar, turnbaugh)$p.value
escobarTest[7] <- roc.test(escobar, wu)$p.value
escobarTest[8] <- roc.test(escobar, zupancic)$p.value


# Goodrich Third

goodrichTest <- rep(0, 8)
goodrichTest[1] <- NA
goodrichTest[2] <- NA
goodrichTest[3] <- NA
goodrichTest[4] <- roc.test(goodrich, hmp)$p.value
goodrichTest[5] <- roc.test(goodrich, ross)$p.value
goodrichTest[6] <- roc.test(goodrich, turnbaugh)$p.value
goodrichTest[7] <- roc.test(goodrich, wu)$p.value
goodrichTest[8] <- roc.test(goodrich, zupancic)$p.value

# HMP Fourth

hmpTest <- rep(0, 8)
hmpTest[1] <- NA
hmpTest[2] <- NA
hmpTest[3] <- NA
hmpTest[4] <- NA
hmpTest[5] <- roc.test(hmp, ross)$p.value
hmpTest[6] <- roc.test(hmp, turnbaugh)$p.value
hmpTest[7] <- roc.test(hmp, wu)$p.value
hmpTest[8] <- roc.test(hmp, zupancic)$p.value

# Ross Fifth

rossTest <- rep(0, 8)
rossTest[1] <- NA
rossTest[2] <- NA
rossTest[3] <- NA
rossTest[4] <- NA
rossTest[5] <- NA
rossTest[6] <- roc.test(ross, turnbaugh)$p.value
rossTest[7] <- roc.test(ross, wu)$p.value
rossTest[8] <- roc.test(ross, zupancic)$p.value

# Turnbaugh Sixth
TurnbaughTest <- rep(0, 8)
TurnbaughTest[1] <- NA
TurnbaughTest[2] <- NA
TurnbaughTest[3] <- NA
TurnbaughTest[4] <- NA
TurnbaughTest[5] <- NA
TurnbaughTest[6] <- NA
TurnbaughTest[7] <- roc.test(turnbaugh, wu)$p.value
TurnbaughTest[8] <- roc.test(turnbaugh, zupancic)$p.value

# Wu Seventh
WuTest <- rep(0, 8)
WuTest[1] <- NA
WuTest[2] <- NA
WuTest[3] <- NA
WuTest[4] <- NA
WuTest[5] <- NA
WuTest[6] <- NA
WuTest[7] <- NA
WuTest[8] <- roc.test(wu, zupancic)$p.value

namesToTest <- c("escobar", "goodrich", "hmp", "ross", "turnbaugh", "wu", 
                 "zupancic")


ROCPValueTable <- as.data.frame(cbind(
  baxterTest, escobarTest, goodrichTest, hmpTest, rossTest, TurnbaughTest, 
  WuTest), stringsAsFactors = FALSE)
rownames(ROCPValueTable) <- c("baxter", "escobar", "goodrich", "hmp", "ross", 
                              "turnbaugh", "wu", "zupancic")


### Get the adjusted P-values for all of them and make another table.

rawPvalues <- c(baxterTest, escobarTest, goodrichTest, hmpTest, rossTest, 
                TurnbaughTest, WuTest)
adjustedPvalues <- p.adjust(rawPvalues, method="bonferroni", n=28)

ROCAdjustPvalueTable <- as.data.frame(cbind(
  adjustedPvalues[1:8], adjustedPvalues[9:16], adjustedPvalues[17:24], 
  adjustedPvalues[25:32], adjustedPvalues[33:40], adjustedPvalues[41:48], 
  adjustedPvalues[49:56]), stringsAsFactors = FALSE)




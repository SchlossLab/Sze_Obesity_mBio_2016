## Strategy 1
## Obesesity and the bacterial microbiome
## Marc Sze
## March 9, 2016

library(statmod)
library(vegan)
library(epiR)
library(AUCRF)
source("code/UsedFunctions.R")

############# BAXTER ######################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

# Reading in the necessary Data
demographics <- read.csv("data/process/Baxter/demographics.v2.csv")
microbiome <- read.csv("data/process/Baxter/data1.subsample.otus.csv")

# Minor modifications and Rownames adustments of data tables
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
demographics <- EditTable(demographics, 1, delRange=1)

# Create Table for phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Baxter/taxonomyKey.txt")

#Add phyla together that are not very abundant and delete them from the table
phyla.table$other <- apply(phyla.table[, c("Deferribacteres", "Deinococcus", "Fusobacteria", 
                                           "Lentisphaerae", "Spirochaetes", 
                                           "Synergistetes", "Tenericutes", "TM7")], 1, sum)

phyla.table <- phyla.table[, -c(3:4, 6:7, 9:12)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.total, phyla.table)

#Generate Alpha Diversity Measures and alpha diversity table
alpha.test <- makeAlphaTable(microbiome)

#Create BMI Classifications needed for analysis
demographics <- AddBMIClass(demographics, "BMI.classification", numbers=FALSE)

bmi <- demographics$BMI
obese <- factor(demographics$obese.num)

######################################################################################## 
########## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################


##Test BMI versus alpha diversity and phyla

baxterH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.065
baxterS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.03923
baxterJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.1254

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

baxterBacter <- wilcox.test(bacter ~ obese) #P-value=0.3175
baxterFirm <- wilcox.test(firm ~ obese) #P-value=0.6817
baxterBF <- wilcox.test(BFratio ~ obese) #P-value=0.6305

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
baxter2 <- adonis(microbiome ~ obese, permutations=1000)
baxterPERM <- baxter2$aov.tab
#Not significant PERMANOVA = 0.06394, pseudo-F = 1.4246


###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
BaxterHRR <- RunRR(alpha.test, demographics, "obese", "H")
## Risk Ratio = 1.14
## CI = 0.70, 1.85
## p-value = 0.608

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
BaxterBFRR <- RunRR(BFRatio, demographics, "obese", "BFRatio")
## Risk Ratio = 1.04
## CI = 0.64, 1.70
## p-value = 0.864

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)

#Need to add phyla and alpha diversity measures to the dataset
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(obese, testset)

#Try AUCRF with default measures provided in readme
set.seed(3)
baxterAUCfit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
baxAUC <- baxterAUCfit$`OOB-AUCopt`
baxKopt <- baxterAUCfit$Kopt


###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

BaxterZH <- scale(H)
BaxterZLogBF <- scale(log(BFratio + 1))
BaxterBMI <- bmi

###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- length(rownames(demographics))
meanAge <- format(round(mean(demographics$Age), 2), nsmall = 2)
SDAge <- format(round(sd(demographics$Age), 2), nsmall = 2)
temporary <- table(demographics$Gender)
males <- temporary[names(temporary) == "m"]
females <- temporary[names(temporary) == "f"]
temporary <- table(demographics$White)
ancestry <- format(round(temporary[names(temporary) == 1] / sum(temporary), 2), nsmall = 2)
meanBMI <- format(round(mean(demographics$BMI), 2), nsmall = 2)
SDBMI <- format(round(sd(demographics$BMI), 2), nsmall = 2)
minBMI <- format(round(min(demographics$BMI), 2), nsmall = 2)
maxBMI <- format(round(max(demographics$BMI), 2), nsmall = 2)

# Get data for P-value table
overallPTable <- as.data.frame(t(c(baxterBacter$p.value, baxterFirm$p.value, 
                                   baxterBF$p.value, baxterH$p.value, 
                                   baxterS$p.value, baxterJ$p.value, 
                                   baxterPERM[1,6])))
colnames(overallPTable) <- c("Bacteroidetes", "Firmicutes", "BFRatio", "Shannon", "OTURich", "Evenness", "BrayC")

OOBAUCAll <- format(round(baxAUC, 2), nsmall = 2)
KoptAll <- format(round(baxKopt, 2), nsmall = 2)

# Get data for the combined analysis with Zscores
combinedData <- as.data.frame(cbind(BaxterZH, BaxterZLogBF, BaxterBMI))
colnames(combinedData) <- c("ZH", "ZLogBF", "BMI")
combinedData$Study <- "Baxter"

# Get data for power simulation
StudyPowerH <- NonParaPowerSim(
  alpha.test, demographics, "H", "obese", n=1000)

set.seed(3)
StudyPowerHRR <- power.fisher.test(
  BaxterHRR$tpos/(BaxterHRR$tpos+BaxterHRR$tneg), 
  BaxterHRR$cpos/(BaxterHRR$cpos+BaxterHRR$cneg), 
  BaxterHRR$tpos+BaxterHRR$tneg, BaxterHRR$cpos+BaxterHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100

StudyPowerBF <- NonParaPowerSim(
  BFRatio, demographics, "BFRatio", "obese", n=1000)

set.seed(3)
StudyPowerBFRR <- power.fisher.test(
  BaxterBFRR$tpos/(BaxterBFRR$tpos+BaxterBFRR$tneg), 
  BaxterBFRR$cpos/(BaxterBFRR$cpos+BaxterBFRR$cneg), 
  BaxterBFRR$tpos+BaxterBFRR$tneg, BaxterBFRR$cpos+BaxterBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100

rm(baxterBacter, baxterFirm, baxterBF, baxterH, baxterS, baxterJ, baxterPERM, 
   BaxterZH, BaxterZLogBF, BaxterBMI, baxter2)



################ ROSS #####################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

microbiome <- read.table("data/process/Ross/RossGoodSub.shared", header=T)
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
metadata <- read.csv("data/process/Ross/s40168-015-0072-y-s1.csv")
sample.match <- read.csv("data/process/Ross/Hispanic_dataset.csv")

#Create a microbiome data table with sample names from metadata
test <- cbind(as.character(sample.match$Run_s), as.character(sample.match$Library_Name_s))
test <- test[order(test[, 1]), ]
rownames(test) <- test[, 1]
keep <- rownames(microbiome)
test <- test[keep, ]
test2 <- cbind(test, microbiome)
microbiome <- EditTable(test2, 2, delRange=c(1:2))
metadata <- EditTable(metadata, 1, delRange=c(1:5))

#Create a metadata file in the order of the microbiome data
order1 <- rownames(microbiome)
metadata <- metadata[order1, ]
rm (sample.match, keep, order1, test, test2)

#Get alpha diversity of the samples
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Ross/taxonomyKey.txt")

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Deinococcus", "Elusimicrobia", "Fusobacteria", 
                                           "Lentisphaerae", 
                                           "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(3:4, 6:7, 9:10)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.total, phyla.table)

#Create BMI Classifications needed for analysis
metadata <- AddBMIClass(metadata, "BMI", numbers=TRUE)

#Get paitent demographics to be tested
bmi <- metadata$BMI
obese <- factor(metadata$obese)

###########################################################################
############# First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

rossH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.2848
rossS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.2492
rossJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.3826

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

rossBacter <- wilcox.test(bacter ~ obese) #P-value=0.2036
rossFirm <- wilcox.test(firm ~ obese) #P-value=0.3799
rossBF <- wilcox.test(BFratio ~ obese) #P-value=0.2207

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
ross2 <- adonis(microbiome ~ obese, permutations=1000)
rossPERM <- ross2$aov.tab
#PERMANOVA=0.8232, pseudo-F=0.6725

###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
RossHRR <- RunRR(alpha.test, metadata, "obese", "H")
## Risk Ratio = 1.20
## CI = 0.80, 1.80
## p-value = 0.382

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
RossBFRR <- RunRR(BFRatio, metadata, "obese", "BFRatio")
## Risk Ratio = 0.87
## CI = 0.58, 1.30
## p-value = 0.503


###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
obese <- factor(metadata$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
rossAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
rossAUC <- rossAUCFit$`OOB-AUCopt`
rossKopt <- rossAUCFit$Kopt
# list of 8 Measures, AUCopt = 0.7410526

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

RossZH <- scale(H)
RossZLogBF <- scale(log(BFratio + 1))
RossBMI <- bmi

###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, 
             format(round(
               mean(metadata$age_at_visit), 2), nsmall = 2))
             
SDAge <- c(SDAge, 
           format(round(
             sd(metadata$age_at_visit), 2), nsmall = 2))

temporary <- table(metadata$sex)
males <- unname(c(males, temporary[names(temporary) == "M"]))
females <- unname(c(females, temporary[names(temporary) == "F"]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, format(round(0, 2), nsmall = 2)))
meanBMI <- c(meanBMI, format(round(mean(metadata$BMI), 2), nsmall = 2))
SDBMI <- c(SDBMI, format(round(sd(metadata$BMI), 2), nsmall = 2))
minBMI <- c(minBMI, format(round(min(metadata$BMI), 2), nsmall = 2))
maxBMI <- c(maxBMI, format(round(max(metadata$BMI), 2), nsmall = 2))

Ross <- c(rossBacter$p.value, rossFirm$p.value, rossBF$p.value, 
          rossH$p.value, rossS$p.value, rossJ$p.value, rossPERM[1,6])
overallPTable <- rbind(overallPTable, Ross)

OOBAUCAll <- c(OOBAUCAll, 
               format(round(rossAUC, 2), nsmall = 2))

KoptAll <- c(KoptAll, 
             format(round(rossKopt, 2), nsmall = 2))


RossData <- as.data.frame(cbind(RossZH, RossZLogBF, RossBMI))
RossData$Study <- "Ross"
colnames(RossData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, RossData)

# Get data for power simulation
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, metadata, "H", "obese", n=1000))
set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  RossHRR$tpos/(RossHRR$tpos+RossHRR$tneg), 
  RossHRR$cpos/(RossHRR$cpos+RossHRR$cneg), 
  RossHRR$tpos+RossHRR$tneg, RossHRR$cpos+RossHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, metadata, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  RossBFRR$tpos/(RossBFRR$tpos+RossBFRR$tneg), 
  RossBFRR$cpos/(RossBFRR$cpos+RossBFRR$cneg), 
  RossBFRR$tpos+RossBFRR$tneg, RossBFRR$cpos+RossBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(rossBacter, rossFirm, rossBF, rossH, rossS, rossJ, rossPERM,RossZH, 
   RossZLogBF, RossBMI, ross2)

############## GOODRICH ##################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

#Read in and match metadata to microbiome data
metadata <- read.csv("data/process/Goodrich/TwinsUKStudy2.csv")
metadata <- EditTable(metadata, 7, delRange=0)
microbiome <- read.table("data/process/Goodrich/GoodrichGoodSub.shared", header=T)
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
keep  <- which(metadata$body_mass_index_s != "<not provided>")
microbiome<- microbiome[keep, ]
metadata <- metadata[keep, ]
rm(keep)

#generate alpha diversity measures with vegan
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Goodrich/taxonomyKey.txt")

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Elusimicrobia", 
                                           "Fusobacteria", "Lentisphaerae", 
                                           "Spirochaetes", "SR1", "Synergistetes", 
                                           "Tenericutes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4, 6:7, 9:13)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.total, phyla.table)

#Create BMI Classifications needed for analysis
metadata$body_mass_index_s <- as.numeric(as.character(metadata$body_mass_index_s))
metadata <- AddBMIClass(metadata, "body_mass_index_s", numbers=TRUE)

#Get paitent demographics to be tested
bmi <- as.numeric(as.character(metadata$body_mass_index_s))
obese <- factor(metadata$obese)


###########################################################################
######### First Level Analysis & Alpha Diversity with BMI #################
###########################################################################

##Test BMI versus alpha diversity and phyla

goodrichH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.614
goodrichS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.3244
goodrichJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.8362

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

goodrichBacter <- wilcox.test(bacter ~ obese) #P-value=0.4353
goodrichFirm <- wilcox.test(firm ~ obese) #P-value=0.2896
goodrichBF <- wilcox.test(BFratio ~ obese) #P-value=0.2702

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
goodrich2 <- adonis(microbiome ~ obese, permutations=1000)
goodrichPERM <- goodrich2$aov.tab
#PERMANOVA=0.3257, pseudo-F=1.0833

###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
GoodrichHRR <- RunRR(alpha.test, metadata, "obese", "H")
## Risk Ratio = 0.836
## CI = 0.592, 1.18
## p-value = 0.31

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
GoodrichBFRR <- RunRR(BFRatio, metadata, "obese", "BFRatio")
## Risk Ratio = 1.02
## CI = 0.72, 1.43
## p-value = 0.93


###########################################################################
############ Classification using AUCRF####################################
###########################################################################

#Create Obese.num group
obese <- factor(metadata$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
goodrichAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
goodAUC <- goodrichAUCFit$`OOB-AUCopt`
goodKopt <- goodrichAUCFit$Kopt
# list of 41 Measures, AUCopt = 0.7662693

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

GoodrichZH <- scale(H)
GoodrichZLogBF <- scale(log(BFratio + 1))
GoodrichBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, format(round(mean(metadata$age_s), 2), nsmall = 2))
SDAge <- c(SDAge, sd(metadata$age_s))
temporary <- table(metadata$sex_s)
males <- unname(c(males, temporary[names(temporary) == 47]))
females <- unname(c(females, temporary[names(temporary) == 48]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, NA))
meanBMI <- c(meanBMI, 
             format(round(
               mean(as.numeric(as.character(metadata$body_mass_index_s))), 2), nsmall = 2))
SDBMI <- c(SDBMI, 
           format(round(
             sd(as.numeric(as.character(metadata$body_mass_index_s))), 2), nsmall = 2))
minBMI <- c(minBMI, 
            format(round(
              min(as.numeric(as.character(metadata$body_mass_index_s))), 2), nsmall = 2))
maxBMI <- c(maxBMI, 
            format(round(
              max(as.numeric(as.character(metadata$body_mass_index_s))), 2), nsmall = 2))

Goodrich <- c(goodrichBacter$p.value, goodrichFirm$p.value, 
              goodrichBF$p.value, goodrichH$p.value, goodrichS$p.value, 
              goodrichJ$p.value, goodrichPERM[1,6])

overallPTable <- rbind(overallPTable, Goodrich)

OOBAUCAll <- c(OOBAUCAll, 
               format(round(goodAUC, 2), nsmall = 2))
KoptAll <- c(KoptAll, 
             format(round(goodKopt, 2), nsmall = 2))


GoodrichData <- as.data.frame(cbind(GoodrichZH, GoodrichZLogBF, GoodrichBMI))
GoodrichData$Study <- "Goodrich"
colnames(GoodrichData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, GoodrichData)

# Get data for power simulation for obese versus non-obese
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, metadata, "H", "obese", n=1000))

set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  GoodrichHRR$tpos/(GoodrichHRR$tpos+GoodrichHRR$tneg), 
  GoodrichHRR$cpos/(GoodrichHRR$cpos+GoodrichHRR$cneg), 
  GoodrichHRR$tpos+GoodrichHRR$tneg, GoodrichHRR$cpos+GoodrichHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, metadata, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  GoodrichBFRR$tpos/(GoodrichBFRR$tpos+GoodrichBFRR$tneg), 
  GoodrichBFRR$cpos/(GoodrichBFRR$cpos+GoodrichBFRR$cneg), 
  GoodrichBFRR$tpos+GoodrichBFRR$tneg, GoodrichBFRR$cpos+GoodrichBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(goodrichBacter, goodrichFirm, goodrichBF, goodrichH, goodrichS, goodrichJ, 
   goodrichPERM, GoodrichZH, GoodrichZLogBF, GoodrichBMI, goodrich2, 
   Ross, Goodrich)


########### ESCOBAR ######################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################
            
microbiome <- read.table("data/process/Escobar/EscobarGoodSub.shared", header=T)
metadata <- read.csv("data/process/Escobar/columbian_dataset.csv")

#Organize the microbiome shared data and metadata
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
metadata <- EditTable(metadata, 5, delRange=c(1:7, 9, 13:15, 16:52))

#Sort metadata into the same order as the microbiome data and get sex information
order1 <- rownames(microbiome)
metadata <- metadata[order1, ]
metadata$sex[metadata$description_s == "Adequate weight male stool sample" | 
               metadata$description_s == "Obese male stool sample" | 
               metadata$description_s == "Overweight male stool sample"] <- "M"
metadata$sex[metadata$description_s == "Adequate weight female stool sample" | 
               metadata$description_s == "Obese female stool sample" | 
               metadata$description_s == "Overweight female stool sample"] <- "F"

rm(order1)

#Get alpha diversity of the samples
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Escobar/taxonomyKey.txt")

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Lentisphaerae", 
                                           "Spirochaetes", "Synergistetes", "TM7", 
                                           "Acidobacteria", "Deinococcus")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4, 6:7, 9:11)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.table, phyla.total)

#Create BMI Classifications needed for analysis
metadata <- AddBMIClass(metadata, "body_mass_index_s", numbers=TRUE)

#Get paitent demographics data to be tested
bmi <- metadata$body_mass_index_s
obese <- factor(metadata$obese)

###########################################################################
############# First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

escobarH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.9483
escobarS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.2307
escobarJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.6187

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

escobarBacter <- wilcox.test(bacter ~ obese) #P-value=0.05563
escobarFirm <- wilcox.test(firm ~ obese) #P-value=0.1307
escobarBF <- wilcox.test(BFratio ~ obese) #P-value=0.08221

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
escobar2 <- adonis(microbiome ~ obese, permutations=1000)
escobarPERM <- escobar2$aov.tab
#PERMANOVA=0.07393, pseudo-F=1.3756

###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
EscobarHRR <- RunRR(alpha.test, metadata, "obese", "H")
## Risk Ratio = 1
## CI = 0.37, 2.70
## p-value = 1

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
EscobarBFRR <- RunRR(BFRatio, metadata, "obese", "BFRatio")
## Risk Ratio = 0.429
## CI = 0.137, 1.23
## p-value = 0.121

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
obese <- factor(metadata$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
escobarAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
escoAUC <- escobarAUCFit$`OOB-AUCopt`
escoKopt <- escobarAUCFit$Kopt
# list of 15 Measures, AUCopt = 0.925

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

EscobarZH <- scale(H)
EscobarZLogBF <- scale(log(BFratio + 1))
EscobarBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, format(round(mean(metadata$age_s), 2), nsmall = 2))
SDAge <- c(SDAge, 
           format(round(sd(metadata$age_s), 2), nsmall = 2))
temporary <- table(metadata$sex)
males <- unname(c(males, temporary[names(temporary) == "M"]))
females <- unname(c(females, temporary[names(temporary) == "F"]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, format(round(0, 2), nsmall = 2)))
meanBMI <- c(meanBMI, 
             format(round(mean(metadata$body_mass_index_s), 2), nsmall = 2))
SDBMI <- c(SDBMI, 
           format(round(
             sd(metadata$body_mass_index_s), 2), nsmall = 2))
minBMI <- c(minBMI, 
            format(round(
              min(metadata$body_mass_index_s), 2), nsmall = 2))
maxBMI <- c(maxBMI, 
            format(round(
              max(metadata$body_mass_index_s), 2), nsmall = 2))


Escobar <- c(escobarBacter$p.value, escobarFirm$p.value, 
             escobarBF$p.value, escobarH$p.value, escobarS$p.value, 
             escobarJ$p.value, escobarPERM[1,6])

overallPTable <- rbind(overallPTable, Escobar)

OOBAUCAll <- c(OOBAUCAll, 
               format(round(escoAUC, 2), nsmall = 2))
KoptAll <- c(KoptAll, 
             format(round(escoKopt, 2), nsmall = 2))


EscobarData <- as.data.frame(cbind(EscobarZH, EscobarZLogBF, EscobarBMI))
EscobarData$Study <- "Escobar"
colnames(EscobarData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, EscobarData)

# Get data for power simulation for obese versus non-obese
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, metadata, "H", "obese", n=1000))

set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  EscobarHRR$tpos/(EscobarHRR$tpos+EscobarHRR$tneg), 
  EscobarHRR$cpos/(EscobarHRR$cpos+EscobarHRR$cneg), 
  EscobarHRR$tpos+EscobarHRR$tneg, EscobarHRR$cpos+EscobarHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, metadata, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  EscobarBFRR$tpos/(EscobarBFRR$tpos+EscobarBFRR$tneg), 
  EscobarBFRR$cpos/(EscobarBFRR$cpos+EscobarBFRR$cneg), 
  EscobarBFRR$tpos+EscobarBFRR$tneg, EscobarBFRR$cpos+EscobarBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(escobarBacter, escobarFirm, escobarBF, escobarH, escobarS, escobarJ, 
   escobarPERM, EscobarZH, EscobarZLogBF, EscobarBMI, escobar2, Escobar)


######## ZUPANCIC #########################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

#Read in and match metadata to microbiome data
metadata <- read.csv("data/process/Zupancic/amish_obesity_table2.csv")
metadata2 <- read.csv("data/process/Zupancic/amish.metadata.csv")
metadata <- metadata[!duplicated(metadata$submitted_sample_id_s), ]
rownames(metadata) <- metadata[, 4]
metadata2 <- metadata2[!duplicated(metadata2$SUBJID), ]
rownames(metadata2) <- metadata2[, 1]
microbiome <- read.table("data/process/Zupancic/ZupancicGoodSub.shared", header=T)
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
keep  <- rownames(microbiome)

metadata <- cbind(metadata[keep, ], metadata2[keep, ])
metadata <- metadata[!duplicated(metadata$submitted_subject_id_s), ]
keep  <- rownames(metadata)
microbiome <- microbiome[keep, ]

metadata <- metadata[complete.cases(metadata), ]
keep <- rownames(metadata)
microbiome <- microbiome[keep, ]
rm(keep, metadata2)

#generate alpha diversity measures with vegan
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Zupancic/taxonomyKey.txt")

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Spirochaetes", 
                                           "Elusimicrobia", "Lentisphaerae", 
                                           "Synergistetes")], 1, sum)
phyla.table <- phyla.table[, -c(3, 5:6, 8:9)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.table, phyla.total)

#Create BMI Classifications needed for analysis
metadata <- AddBMIClass(metadata, "BMI", numbers=TRUE)

#create groups to be used
bmi <- metadata$BMI
obese <- factor(metadata$obese)

###########################################################################
############# First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

zupancicH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.3108
zupancicS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.161
zupancicJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.4407

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

zupancicBacter <- wilcox.test(bacter ~ obese) #P-value=0.5674
zupancicFirm <- wilcox.test(firm ~ obese) #P-value=0.6016
zupancicBF <- wilcox.test(BFratio ~ obese) #P-value=0.5919

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
zupancic2 <- adonis(microbiome ~ obese, permutations=1000)
zupancicPERM <- zupancic2$aov.tab
#PERMANOVA=0.1788, pseudo-F=1.1332

###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
ZupancicHRR <- RunRR(alpha.test, metadata, "obese", "H")
## Risk Ratio = 1.37
## CI = 0.94, 2.01
## p-value = 0.104

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
ZupancicBFRR <- RunRR(BFRatio, metadata, "obese", "BFRatio")
## Risk Ratio = 0.972
## CI = 0.669, 1.41
## p-value = 0.883

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
obese <- factor(metadata$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
zupancicAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
zupAUC <- zupancicAUCFit$`OOB-AUCopt`
zupKopt <- zupancicAUCFit$Kopt
# list of 42 Measures, AUCopt = 0.7310296


###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

ZupancicZH <- scale(H)
ZupancicZLogBF <- scale(log(BFratio + 1))
ZupancicBMI <- bmi

###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, 
             format(round(
               mean(as.numeric(as.character(metadata$AGE))), 2), nsmall = 2))
SDAge <- c(SDAge, 
           format(round(
             sd(as.numeric(as.character(metadata$AGE))), 2), nsmall = 2))
temporary <- table(metadata$sex_s)
males <- unname(c(males, temporary[names(temporary) == "male"]))
females <- unname(c(females, temporary[names(temporary) == "female"]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, format(round(1, 2), nsmall = 2)))
meanBMI <- c(meanBMI, 
             format(round(mean(metadata$BMI), 2), nsmall = 2))
SDBMI <- c(SDBMI, 
           format(round(sd(metadata$BMI), 2), nsmall = 2))
minBMI <- c(minBMI, 
            format(round(min(metadata$BMI), 2), nsmall = 2))
maxBMI <- c(maxBMI, 
            format(round(max(metadata$BMI), 2), nsmall = 2))

Zupancic <- c(zupancicBacter$p.value, zupancicFirm$p.value, 
              zupancicBF$p.value, zupancicH$p.value, zupancicS$p.value, 
              zupancicJ$p.value, zupancicPERM[1,6])

overallPTable <- rbind(overallPTable, Zupancic)

OOBAUCAll <- c(OOBAUCAll, 
               format(round(zupAUC, 2), nsmall = 2))

KoptAll <- c(KoptAll, 
               format(round(zupKopt, 2), nsmall = 2))


ZupancicData <- as.data.frame(cbind(ZupancicZH, ZupancicZLogBF, ZupancicBMI))
ZupancicData$Study <- "Zupancic"
colnames(ZupancicData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, ZupancicData)

# Get data for power simulation for obese versus non-obese
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, metadata, "H", "obese", n=1000))

set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  ZupancicHRR$tpos/(ZupancicHRR$tpos+ZupancicHRR$tneg), 
  ZupancicHRR$cpos/(ZupancicHRR$cpos+ZupancicHRR$cneg), 
  ZupancicHRR$tpos+ZupancicHRR$tneg, ZupancicHRR$cpos+ZupancicHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, metadata, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  ZupancicBFRR$tpos/(ZupancicBFRR$tpos+ZupancicBFRR$tneg), 
  ZupancicBFRR$cpos/(ZupancicBFRR$cpos+ZupancicBFRR$cneg), 
  ZupancicBFRR$tpos+ZupancicBFRR$tneg, ZupancicBFRR$cpos+ZupancicBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(zupancicBacter, zupancicFirm, zupancicBF, zupancicH, zupancicS, zupancicJ, 
   zupancicPERM, ZupancicZH, ZupancicZLogBF, ZupancicBMI, zupancic2, 
   Zupancic)


############### HMP HMP HMP ##############################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

##Read in relevant data
microbiome <- read.table("data/process/HMP/Stool.an.0.03.subsample.shared", header=T)
meta.cat <- read.table("data/process/HMP/categorical.metadata", header=T)
meta.cont <- read.table("data/process/HMP/continuous.metadata", header=T)

#Only interested in first visit microbiome data
#need to create microbiome data set for only that
microbiome <- microbiome[grep("\\.01\\.", microbiome$Group), ] #selecting by .01.

#extracting only numeric before first "."
microb.rownames <- gsub("([0-9]+).*", "\\1", microbiome$Group) 
rownames(microbiome) <- microb.rownames

#Get rid of information not used in downstram analysis
microbiome <- microbiome[, -c(1:3)]

#Subset data to the microbiome total n
meta.cat <- meta.cat[microb.rownames, ]
meta.cont <- meta.cont[microb.rownames, ]

#generate alpha diversity measures with vegan
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/HMP/phyla.data.csv", csv=T)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Deinococcus-Thermus", 
                                           "Fusobacteria", "Lentisphaerae", 
                                           "Spirochaetes", "Synergistetes", 
                                           "Tenericutes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4, 6, 7, 9, 10:12)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.total, phyla.table)

#Create BMI Classifications needed for analysis
meta.cont <- AddBMIClass(meta.cont, "DTPBMI", numbers=TRUE)

#Get paitent demographics data to be tested
bmi <- meta.cont$DTPBMI
obese <- factor(meta.cont$obese)

###########################################################################
############# First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

HMPH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.5772
HMPS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.9654
HMPJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.4667


#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

HMPBacter <- wilcox.test(bacter ~ obese) #P-value=0.459
HMPFirm <- wilcox.test(firm ~ obese) #P-value=0.7864
HMPBF <- wilcox.test(BFratio ~ obese) #P-value=0.5858

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
HMP2 <- adonis(microbiome ~ obese, permutations=1000)
HMPPERM <- HMP2$aov.tab
#PERMANOVA=0.8112, pseudo-F=0.7024

###########################################################################
############# Relative Risk ###############################################
###########################################################################

# Run Shannon Diversity RR Test
HMPHRR <- RunRR(alpha.test, meta.cont, "obese", "H")
## Risk Ratio = 1.60
## CI = 0.77, 3.35
## p-value = 0.214

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
HMPBFRR <- RunRR(BFRatio, meta.cont, "obese", "BFRatio")
## Risk Ratio = 1.00
## CI = 0.49, 2.04
## p-value = 1

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
obese <- factor(meta.cont$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
HMPAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
HMPAUC <- HMPAUCFit$`OOB-AUCopt`
HMPKopt <- HMPAUCFit$Kopt
# list of 6 Measures, AUCopt = 0.7031773

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

HMPZH <- scale(H)
HMPZLogBF <- scale(log(BFratio + 1))
HMPBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(meta.cat)))
meanAge <- c(meanAge, 
             format(round(mean(meta.cont$AGEENR), 2), nsmall = 2))
SDAge <- c(SDAge, 
           format(round(sd(meta.cont$AGEENR), 2), nsmall = 2))
temporary <- table(meta.cat$GENDER_C)
males <- unname(c(males, temporary[names(temporary) == "Male"]))
females <- unname(c(females, temporary[names(temporary) == "Female"]))
temporary <- table(meta.cat$WHITE_C)
ancestry <- unname(c(ancestry, temporary[names(temporary) == "Yes"] / sum(temporary)))
meanBMI <- c(meanBMI, 
             format(round(mean(meta.cont$DTPBMI), 2), nsmall = 2))
SDBMI <- c(SDBMI, 
           format(round(sd(meta.cont$DTPBMI), 2), nsmall = 2))
minBMI <- c(minBMI, 
            format(round(min(meta.cont$DTPBMI), 2), nsmall = 2))
maxBMI <- c(maxBMI, 
            format(round(max(meta.cont$DTPBMI), 2), nsmall = 2))

HMP <- c(HMPBacter$p.value, HMPFirm$p.value, 
         HMPBF$p.value, HMPH$p.value, HMPS$p.value, 
         HMPJ$p.value, HMPPERM[1,6])

overallPTable <- rbind(overallPTable, HMP)


OOBAUCAll <- c(OOBAUCAll, 
               format(round(HMPAUC, 2), nsmall = 2))
KoptAll <- c(KoptAll, 
             format(round(HMPKopt, 2), nsmall = 2))


HMPData <- as.data.frame(cbind(HMPZH, HMPZLogBF, HMPBMI))
HMPData$Study <- "HMP"
colnames(HMPData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, HMPData)

# Get data for power simulation
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, meta.cont, "H", "obese", n=1000))

set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  HMPHRR$tpos/(HMPHRR$tpos+HMPHRR$tneg), 
  HMPHRR$cpos/(HMPHRR$cpos+HMPHRR$cneg), 
  HMPHRR$tpos+HMPHRR$tneg, HMPHRR$cpos+HMPHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, meta.cont, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  HMPBFRR$tpos/(HMPBFRR$tpos+HMPBFRR$tneg), 
  HMPBFRR$cpos/(HMPBFRR$cpos+HMPBFRR$cneg), 
  HMPBFRR$tpos+HMPBFRR$tneg, HMPBFRR$cpos+HMPBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(HMPBacter, HMPFirm, HMPBF, HMPH, HMPS, HMPJ, HMPPERM, HMPZH, HMPZLogBF, 
   HMPBMI, HMP2, HMP, meta.cat, meta.cont)


############## Wu Wu Wu ###################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

microbiome <- read.table("data/process/Wu/WuGoodSub.shared", header=T)
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
metadata <- read.table("data/process/Wu/bmi_info.txt", header=T)

#Match the metadata now with the microbiome data
namesToKeep <- rownames(microbiome)
metadata <- metadata[namesToKeep, ]
rm(namesToKeep)

#Get alpha diversity of the samples
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Wu/taxonomyKey.txt")

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", 
                                           "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(4, 6:7)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:6)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(keep, phyla.table, phyla.total)

#Create BMI Classifications needed for analysis
metadata <- AddBMIClass(metadata, "bmi", numbers=TRUE)

#Get paitent demographics to be tested
bmi <- metadata$bmi
obese <- factor(metadata$obese)

###########################################################################
######### First Level Analysis & Alpha Diversity with BMI #################
###########################################################################

##Test BMI versus alpha diversity and phyla

WuH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.9089
WuS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.2798
WuJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.3804

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

WuBacter <- wilcox.test(bacter ~ obese) #P-value=0.7124
WuFirm <- wilcox.test(firm ~ obese) #P-value=0.7506
WuBF <- wilcox.test(BFratio ~ obese) #P-value=0.9899

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
Wu2 <- adonis(microbiome ~ obese, permutations=1000)
WuPERM <- Wu2$aov.tab
#PERMANOVA=0.5335, pseudo-F=0.92196


###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
WuHRR <- RunRR(alpha.test, metadata, "obese", "H")
## Risk Ratio = 0.65
## CI = 0.135, 3.05
## p-value = 0.615

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
WuBFRR <- RunRR(BFRatio, metadata, "obese", "BFRatio")
## Risk Ratio = 1.45
## CI = 0.308, 6.95
## p-value = 0.668

###########################################################################
############ Classification using AUCRF####################################
###########################################################################

#Create Obese.num group
obese <- factor(metadata$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
WuAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
WuAUC <- WuAUCFit$`OOB-AUCopt`
WuKopt <- WuAUCFit$Kopt
# list of 15 Measures, AUCopt = 0.8965517

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

WuZH <- scale(H)
WuZLogBF <- scale(log(BFratio + 1))
WuBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, format(round(mean(metadata$age), 2), nsmall = 2))
SDAge <- c(SDAge, format(round(sd(metadata$age), 2), nsmall = 2))
temporary <- table(metadata$sex1m2f)
males <- unname(c(males, temporary[names(temporary) == 1]))
females <- unname(c(females, temporary[names(temporary) == 2]))
#temporary <- table(select.meta.cat$WHITE_C)
ancestry <- unname(c(ancestry, NA))
meanBMI <- c(meanBMI, format(round(mean(metadata$bmi), 2), nsmall = 2))
SDBMI <- c(SDBMI, format(round(sd(metadata$bmi), 2), nsmall = 2))
minBMI <- c(minBMI, format(round(min(metadata$bmi), 2), nsmall = 2))
maxBMI <- c(maxBMI, format(round(max(metadata$bmi), 2), nsmall = 2))

Wu <- c(WuBacter$p.value, WuFirm$p.value, 
        WuBF$p.value, WuH$p.value, WuS$p.value, 
        WuJ$p.value, WuPERM[1,6])

overallPTable <- rbind(overallPTable, Wu)

OOBAUCAll <- c(OOBAUCAll, 
               format(round(WuAUC, 2), nsmall = 2))
KoptAll <- c(KoptAll, 
             format(round(WuKopt, 2), nsmall = 2))


WuData <- as.data.frame(cbind(WuZH, WuZLogBF, WuBMI))
WuData$Study <- "Wu"
colnames(WuData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, WuData)

#Create BMI groups
combinedData$BMICat[combinedData$BMI<=24] <- "Normal"
combinedData$BMICat[combinedData$BMI>24 & 
                      combinedData$BMI<30] <- "Overweight"
combinedData$BMICat[combinedData$BMI>=30] <- "Obese"

#Create Obese Yes/No groups
combinedData$Obese[combinedData$BMICat=="Normal" | 
                     combinedData$BMICat=="Overweight"] <- "No"
combinedData$Obese[combinedData$BMICat=="Obese"] <- "Yes"

# Get data for power simulation
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, metadata, "H", "obese", n=1000))

set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  WuHRR$tpos/(WuHRR$tpos+WuHRR$tneg), 
  WuHRR$cpos/(WuHRR$cpos+WuHRR$cneg), 
  WuHRR$tpos+WuHRR$tneg, WuHRR$cpos+WuHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, metadata, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  WuBFRR$tpos/(WuBFRR$tpos+WuBFRR$tneg), 
  WuBFRR$cpos/(WuBFRR$cpos+WuBFRR$cneg), 
  WuBFRR$tpos+WuBFRR$tneg, WuBFRR$cpos+WuBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(WuBacter, WuFirm, WuBF, WuH, WuS, WuJ, WuPERM, WuZH, WuZLogBF, 
   WuBMI, Wu2, Wu)



########### Turnbaugh ####################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

microbiome <- read.table("data/process/Turnbaugh/TurnbaughSub.shared", header=T)
microbiome <- EditTable(microbiome, 2, delRange=c(1:3))
metadata <- read.csv("data/process/Turnbaugh/turnbaugh.metadata.csv")
metadata <- EditTable(metadata, 4, delRange=4)

keep1 <- which(metadata$Sample == 1) # n =146

#use the first sampling
#subset data for only the first sampling
metadata <- metadata[keep1, ]
microbiome <- microbiome[keep1, ]
rm(keep1)

#generate alpha diversity measures with vegan
alpha.test <- makeAlphaTable(microbiome)

#Get phyla information
phyla.table <- MakePhylaTable(microbiome, "data/process/Turnbaugh/phyla.csv", 
                              csv = T)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Synergistetes", "TM7", 
                                           "Lentisphaerae", "Spirochaetes")], 1, sum)
phyla.table <- phyla.table[, -c(4:5, 7:9)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.table, phyla.total)

#Create BMI Classifications needed for analysis
metadata <- AddBMIClass(metadata, "BMI.category", numbers=FALSE)

#create groups to be used
obese <- factor(metadata$obese)
bmi <- metadata$BMI.category

####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

turnbaughH <- wilcox.test(alpha.test$H ~ obese) #P-value=0.1162
turnbaughS <- wilcox.test(alpha.test$S ~ obese) #P-value=0.05479
turnbaughJ <- wilcox.test(alpha.test$J ~ obese) #P-value=0.1748

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

turnbaughBacter <- wilcox.test(bacter ~ obese) #P-value=0.8048
turnbaughFirm <- wilcox.test(firm ~ obese) #P-value=0.8587
turnbaughBF <- wilcox.test(BFratio ~ obese) #P-value=0.7886

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
turnbaugh2 <- adonis(microbiome ~ obese, permutations=1000)
turnbaughPERM <- turnbaugh2$aov.tab
#PERMANOVA=0.09491, pseudo-F=1.2114

###########################################################################
############ Relative Risk#################################################
###########################################################################

# Run Shannon Diversity RR Test
TurnbaughHRR <- RunRR(alpha.test, metadata, "obese", "H")
## Risk Ratio = 1.11
## CI = 0.883, 1.4
## p-value = 0.376

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)
TurnbaughBFRR <- RunRR(BFRatio, metadata, "obese", "BFRatio")
## Risk Ratio = 1.02
## CI = 0.812, 1.28
## p-value = 0.859


###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
obese <- factor(metadata$obese.num)
H <- alpha.test$H
S <- alpha.test$S
J <- alpha.test$J

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
TurnbaughAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
TurnAUC <- TurnbaughAUCFit$`OOB-AUCopt`
TurnKopt <- TurnbaughAUCFit$Kopt
# list of 11 Measures, AUCopt = 0.7764883

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

TurnbaughZH <- scale(H)
TurnbaughZLogBF <- scale(log(BFratio + 1))
TurnbaughObese <- as.character(metadata$obese)
TurnbaughBMICat <- as.character(bmi)


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, "21-32")
SDAge <- c(SDAge, NA)
#temporary <- table(metadata$Sex)
males <- unname(c(males, NA))
females <- unname(c(females, NA))
temporary <- table(metadata$Ancestry)
ancestry <- unname(c(ancestry, 
                     format(round(
                       temporary[names(temporary) == "EA"] / sum(temporary), 2), 
                       nsmall = 2)))
meanBMI <- c(meanBMI, NA)
SDBMI <- c(SDBMI, NA)
minBMI <- c(minBMI, NA)
maxBMI <- c(maxBMI, NA)

Turnbaugh <- c(turnbaughBacter$p.value, turnbaughFirm$p.value, 
               turnbaughBF$p.value, turnbaughH$p.value, turnbaughS$p.value, 
               turnbaughJ$p.value, turnbaughPERM[1,6])

overallPTable <- rbind(overallPTable, Turnbaugh)

OOBAUCAll <- c(OOBAUCAll, 
               format(round(TurnAUC, 2), nsmall = 2))
KoptAll <- c(KoptAll, 
             format(round(TurnKopt, 2), nsmall = 2))


TurnbaughData <- as.data.frame(cbind(TurnbaughZH, TurnbaughZLogBF))
TurnbaughData$Study <- "Turnbaugh"
TurnbaughData <- as.data.frame(cbind(TurnbaughData, TurnbaughObese, TurnbaughBMICat))
colnames(TurnbaughData) <- c("ZH", "ZLogBF", "Study", "Obese", "BMICat")

# Get data for power simulation
StudyPowerH <- c(StudyPowerH, 
                 NonParaPowerSim(
                   alpha.test, metadata, "H", "obese", n=1000))

set.seed(3)
StudyPowerHRR <- c(StudyPowerHRR, power.fisher.test(
  TurnbaughHRR$tpos/(TurnbaughHRR$tpos+TurnbaughHRR$tneg), 
  TurnbaughHRR$cpos/(TurnbaughHRR$cpos+TurnbaughHRR$cneg), 
  TurnbaughHRR$tpos+TurnbaughHRR$tneg, TurnbaughHRR$cpos+TurnbaughHRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

StudyPowerBF <- c(StudyPowerBF, 
                  NonParaPowerSim(
                    BFRatio, metadata, "BFRatio", "obese", n=1000))

set.seed(3)
StudyPowerBFRR <- c(StudyPowerBFRR, power.fisher.test(
  TurnbaughBFRR$tpos/(TurnbaughBFRR$tpos+TurnbaughBFRR$tneg), 
  TurnbaughBFRR$cpos/(TurnbaughBFRR$cpos+TurnbaughBFRR$cneg), 
  TurnbaughBFRR$tpos+TurnbaughBFRR$tneg, 
  TurnbaughBFRR$cpos+TurnbaughBFRR$cneg, 
  alpha=0.05, nsim=1000, alternative="two.sided")*100)

rm(turnbaughBacter, turnbaughFirm, turnbaughBF, turnbaughH, turnbaughS, 
   turnbaughJ, turnbaughPERM, TurnbaughZH, TurnbaughZLogBF, TurnbaughObese, 
   TurnbaughBMICat, turnbaugh2, Turnbaugh)


####################################################################################
############################ Output Tables #########################################
####################################################################################

demographicsTable <- as.data.frame(cbind(totalN, meanAge, SDAge, males, females, 
                                         ancestry, meanBMI, SDBMI, minBMI, maxBMI))
demographicsTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Turnbaugh")
write.csv(demographicsTable, "results/tables/denovodemographicsTable.csv")


write.csv(combinedData, "results/tables/denovoCombinedData.csv")

write.csv(TurnbaughData, "results/tables/denovoTurnbaughData.csv")

rownames(overallPTable) <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                             "HMP", "Wu", "Turnbaugh")
write.csv(overallPTable, "results/tables/denovoOverallPTable.csv")

tposH <- c(BaxterHRR$tpos, RossHRR$tpos, GoodrichHRR$tpos, EscobarHRR$tpos, 
           ZupancicHRR$tpos, HMPHRR$tpos, WuHRR$tpos, TurnbaughHRR$tpos)
tnegH <- c(BaxterHRR$tneg, RossHRR$tneg, GoodrichHRR$tneg, EscobarHRR$tneg, 
           ZupancicHRR$tneg, HMPHRR$tneg, WuHRR$tneg, TurnbaughHRR$tneg)
cposH <- c(BaxterHRR$cpos, RossHRR$cpos, GoodrichHRR$cpos, EscobarHRR$cpos, 
           ZupancicHRR$cpos, HMPHRR$cpos, WuHRR$cpos, TurnbaughHRR$cpos)
cnegH <- c(BaxterHRR$cneg, RossHRR$cneg, GoodrichHRR$cneg, EscobarHRR$cneg, 
           ZupancicHRR$cneg, HMPHRR$cneg, WuHRR$cneg, TurnbaughHRR$cneg)
RRH <- c(BaxterHRR$RR, RossHRR$RR, GoodrichHRR$RR, EscobarHRR$RR, 
         ZupancicHRR$RR, HMPHRR$RR, WuHRR$RR, TurnbaughHRR$RR)
lowH <- c(BaxterHRR$lowCI, RossHRR$lowCI, GoodrichHRR$lowCI, EscobarHRR$lowCI, 
          ZupancicHRR$lowCI, HMPHRR$lowCI, WuHRR$lowCI, TurnbaughHRR$lowCI)
highH <- c(BaxterHRR$highCI, RossHRR$highCI, GoodrichHRR$highCI, 
           EscobarHRR$highCI, ZupancicHRR$highCI, HMPHRR$highCI, WuHRR$highCI, 
           TurnbaughHRR$highCI)

ShannonRRTable <- as.data.frame(cbind(tposH, tnegH, cposH, cnegH, 
                                      RRH, lowH, highH))
ShannonRRTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Turnbaugh")
write.csv(ShannonRRTable, "results/tables/denovoShannonRRTable.csv")

tposBF <- c(BaxterBFRR$tpos, RossBFRR$tpos, GoodrichBFRR$tpos, 
            EscobarBFRR$tpos, ZupancicBFRR$tpos, HMPBFRR$tpos, 
            WuBFRR$tpos, TurnbaughBFRR$tpos)
tnegBF <- c(BaxterBFRR$tneg, RossBFRR$tneg, GoodrichBFRR$tneg, 
            EscobarBFRR$tneg, ZupancicBFRR$tpos, HMPBFRR$tneg, 
            WuBFRR$tneg, TurnbaughBFRR$tneg)
cposBF <- c(BaxterBFRR$cpos, RossBFRR$cpos, GoodrichBFRR$cpos, 
            EscobarBFRR$cpos, ZupancicBFRR$cpos, HMPBFRR$cpos, 
            WuBFRR$cpos, TurnbaughBFRR$cpos)
cnegBF <- c(BaxterBFRR$cneg, RossBFRR$cneg, GoodrichBFRR$cneg, 
            EscobarBFRR$cneg, ZupancicBFRR$cneg, HMPBFRR$cneg, 
            WuBFRR$cneg, TurnbaughBFRR$cneg)
RRBF <- c(BaxterBFRR$RR, RossBFRR$RR, GoodrichBFRR$RR, EscobarBFRR$RR, 
          ZupancicBFRR$RR, HMPBFRR$RR, WuBFRR$RR, TurnbaughBFRR$RR)
lowBF <- c(BaxterBFRR$lowCI, RossBFRR$lowCI, GoodrichBFRR$lowCI, 
           EscobarBFRR$lowCI, ZupancicBFRR$lowCI, HMPBFRR$lowCI, 
           WuBFRR$lowCI, TurnbaughBFRR$lowCI)
highBF <- c(BaxterBFRR$highCI, RossBFRR$highCI, GoodrichBFRR$highCI, 
            EscobarBFRR$highCI, ZupancicBFRR$highCI, HMPBFRR$highCI, 
            WuBFRR$highCI, TurnbaughBFRR$highCI)

BFRatioRRTable <- as.data.frame(cbind(tposBF, tnegBF, cposBF, cnegBF, RRBF, 
                                      lowBF, highBF))
BFRatioRRTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Turnbaugh")
write.csv(BFRatioRRTable, "results/tables/denovoBFRatioRRTable.csv")

AUCRFDataTable <- as.data.frame(cbind(OOBAUCAll, KoptAll))
AUCRFDataTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Turnbaugh")
write.csv(AUCRFDataTable, "results/tables/denovoAUCRFDataTable.csv")

PowerTable <- as.data.frame(cbind(StudyPowerH, StudyPowerHRR, 
                                  StudyPowerBFRR, StudyPowerBF))
PowerTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Turnbaugh")
write.csv(PowerTable, "results/tables/denovoPowerTable.csv")


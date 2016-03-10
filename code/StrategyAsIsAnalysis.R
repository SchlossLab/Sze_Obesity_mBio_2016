## Strategy 1
## Obesesity and the bacterial microbiome
## Marc Sze
## March 9, 2016

library(vegan)
library(epiR)
library(AUCRF)

############# BAXTER ######################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

# Reading in the necessary Data
setwd("C:/Users/marc/Desktop/obesity2/data1")
demographics <- read.csv("demographics.v2.csv")
microbiome <- read.csv("data1.subsample.otus.csv")
phylogenetic.info <- read.csv("data1.summary.taxonomy.csv")
taxonomy <- read.csv("data1.taxonomy.csv")

# Minor modifications and Rownames adustments of data tables
rownames(demographics) <- demographics[,1]
rownames(microbiome) <- microbiome[,2]
demographics <- demographics[,-1]
microbiome <- microbiome[,-2]
test.samples <- rownames(demographics)
microbiome <- microbiome[,-c(1:2)]

# Create Table for phyla information
phylogenetic.info <- phylogenetic.info[-c(1:12),]
phyla <- which(phylogenetic.info$taxlevel == 2)
phyla.table <- phylogenetic.info[phyla, ]
phyla.table <- phyla.table[,-c(1:2,4:5)]
phyla.names <- as.character(phyla.table[,1])
phyla.table <- phyla.table[,-1]
phyla.table <- as.data.frame(t(phyla.table))
colnames(phyla.table) <- phyla.names
rownames(phyla.table) <- rownames(demographics)


#Add phyla together that are not very abundant and delete them from the table
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Deferribacteres", "Deinococcus-Thermus", "Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "Tenericutes", "Cyanobacteria_Chloroplast", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4:7, 9:11, 13, 15)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]

#Generate Alpha Diversity Measures and alpha diversity table
H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Create age quartiles

demographics$age.quartile <- with(demographics, cut(Age, breaks=quantile(Age, probs=seq(0,1, by=0.25)), include.lowest=TRUE))

#Create Obese Yes/No groups
demographics$obese[demographics$BMI.classification=="Normal" | demographics$BMI.classification=="Overweight"] <- "No"
demographics$obese[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Yes"

#Create Obese.num groups
demographics$obese.num[demographics$obese=="No"] <- 0
demographics$obese.num[demographics$obese=="Yes"] <- 1
obese <- factor(demographics$obese.num)

######################################################################################## 
########## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

#Create a column with obese and extreme obese as single entity
demographics$BMIclass2[demographics$BMI.classification=="Normal"] <- "Normal"
demographics$BMIclass2[demographics$BMI.classification=="Overweight"] <- "Overweight"                     
demographics$BMIclass2[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Obese"


##Test BMI versus alpha diversity and phyla

bmi <- demographics$BMI

baxterH <- wilcox.test(H ~ obese) #P-value=0.065
baxterS <- wilcox.test(S ~ obese) #P-value=0.03923
baxterJ <- wilcox.test(J ~ obese) #P-value=0.1254

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

#Generate median values and put them into existing alpha.test dataframe

#Shannon diversity
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

##Run the RR for Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
BaxtotalN <- length(orderedHCat)
BaxHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
BaxHighShannonGroup <- as.data.frame(table(test3[c(1:BaxHHighTotal), 2]))
BaxLowShannonGroup <- as.data.frame(table(test3[c((BaxHHighTotal + 1):BaxtotalN), 2]))
#Group1 (Higher than median), obese = 22 and non-obese = 64
#Group2 (Lower than median), obese = 25 and non-obese = 61

group1 <- c(BaxHighShannonGroup[2, 2], BaxHighShannonGroup[1, 2])
group2 <- c(BaxLowShannonGroup[2, 2], BaxLowShannonGroup[1, 2])
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

baxterHEpi <- epi.2by2(r.test, method="cohort.count")
baxterHMassoc <- baxterHEpi$massoc
baxterHRR <- baxterHMassoc$RR.strata.score
baxterHRRsig <- baxterHMassoc$chisq.strata
## Risk Ratio = 1.14
## CI = 0.70, 1.85
## p-value = 0.608

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
BaxBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
BaxHighBFGroup <- as.data.frame(table(test4[c(1:BaxBFHighTotal), 2]))
BaxLowBFGroup <- as.data.frame(table(test4[c((BaxBFHighTotal + 1):BaxtotalN), 2]))
#Group1 (Higher than median), obese = 23 and non-obese = 63
#Group2 (Lower than median), obese = 24 and non-obese = 62

group1 <- c(BaxHighBFGroup[2, 2], BaxHighBFGroup[1, 2])
group2 <- c(BaxLowBFGroup[2, 2], BaxLowBFGroup[1, 2])
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

baxterBFEpi <- epi.2by2(r.test, method="cohort.count")
baxterBFMassoc <- baxterBFEpi$massoc
baxterBFRR <- baxterBFMassoc$RR.strata.score
baxterBFRRsig <- baxterBFMassoc$chisq.strata
## Risk Ratio = 1.04
## CI = 0.64, 1.70
## p-value = 0.864

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#generate test set
# get rid of values that only have 0 and something else
testset <- Filter(function(x)(length(unique(x))>2), microbiome)
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)

#Need to add phyla and alpha diversity measures to the dataset

testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(obese, testset)

#Try AUCRF with default measures provided in readme
#set.seed(3)
#baxterAUCfit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

BaxterZH <- scale(H)
BaxterZLogBF <- scale(log(BFratio))
BaxterBMI <- bmi

###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- length(rownames(demographics))
meanAge <- mean(demographics$Age)
SDAge <- sd(demographics$Age)
temporary <- table(demographics$Gender)
males <- temporary[names(temporary) == "m"]
females <- temporary[names(temporary) == "f"]
temporary <- table(demographics$White)
ancestry <- temporary[names(temporary) == 1] / sum(temporary)
meanBMI <- mean(demographics$BMI)
SDBMI <- sd(demographics$BMI)
minBMI <- min(demographics$BMI)
maxBMI <- max(demographics$BMI)

# Get data for P-value table
overallPTable <- as.data.frame(t(c(baxterBacter$p.value, baxterFirm$p.value, 
                                   baxterBF$p.value, baxterH$p.value, 
                                   baxterS$p.value, baxterJ$p.value, 
                                   baxterPERM[1,6])))
colnames(overallPTable) <- c("Bacteroidetes", "Firmicutes", "BFRatio", "Shannon", "OTURich", "Evenness", "BrayC")

# Get Data for Forest Plots
tposH <- BaxLowShannonGroup[2, 2]
tnegH <- BaxLowShannonGroup[1, 2]
cposH <- BaxHighShannonGroup[2, 2]
cnegH <- BaxHighShannonGroup[1, 2]
RRH <- baxterHRR[1,1]
lowH <- baxterHRR[1,2]
highH <- baxterHRR[1,3]

tposBF <- BaxLowBFGroup[2, 2]
tnegBF <- BaxLowBFGroup[1, 2]
cposBF <- BaxHighBFGroup[2, 2]
cnegBF <- BaxHighBFGroup[1, 2]
RRBF <- baxterBFRR[1,1]
lowBF <- baxterBFRR[1,2]
highBF <- baxterBFRR[1,3]

# Get data for the combined analysis with Zscores
combinedData <- as.data.frame(cbind(BaxterZH, BaxterZLogBF, BaxterBMI))
colnames(combinedData) <- c("ZH", "ZLogBF", "BMI")
combinedData$Study <- "Baxter"

rm(BaxLowShannonGroup, BaxHighShannonGroup, baxterHRR, baxterBacter, 
   baxterFirm, baxterBF, baxterH, baxterS, baxterJ, baxterPERM, baxterHEpi, 
   BaxterZH, BaxterZLogBF, BaxterBMI, baxter2)



################ ROSS #####################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/hispanic/")

hispanic.microb <- read.table("RossPhylotypeSub.shared", header=T)
rownames(hispanic.microb) <- hispanic.microb[, 2]
hispanic.microb <- hispanic.microb[, -c(1:3)]

metadata <- read.csv("s40168-015-0072-y-s1.csv")
sample.match <- read.csv("Hispanic_dataset.csv")

#Create a microbiome data table with sample names from metadata
test <- cbind(as.character(sample.match$Run_s), as.character(sample.match$Library_Name_s))
test <- test[order(test[, 1]), ]
rownames(test) <- test[, 1]
keep <- rownames(hispanic.microb)
test <- test[keep, ]
test2 <- cbind(test, hispanic.microb)
rownames(test2) <- test2[, 2]
his.microb.edit <- test2[, -c(1:2)]
edit.metadata <- metadata[, -c(2:5)]
rownames(edit.metadata) <- edit.metadata[, 1]
edit.metadata <- edit.metadata[, -1]

#Create a metadata file in the order of the microbiome data
order1 <- rownames(his.microb.edit)
edit.metadata2 <- edit.metadata[order1, ]

#Get alpha diversity of the samples
H <- diversity(his.microb.edit)
S <- specnumber(his.microb.edit)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.table("taxonomyKey.txt", header=T)
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(his.microb.edit)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- his.microb.edit
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("p__", "p__Cyanobacteria", "p__Elusimicrobia", "p__Fusobacteria", "p__Lentisphaerae", "p__Synergistetes", "p__Tenericutes", "p__TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4:5, 7:8, 10:12)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100

#Create BMI groups
edit.metadata2$BMI.class[edit.metadata2$BMI<=24] <- "Normal"
edit.metadata2$BMI.class[edit.metadata2$BMI>24 & edit.metadata2$BMI<30] <- "Overweight"
edit.metadata2$BMI.class[edit.metadata2$BMI>=30 & edit.metadata2$BMI<40] <- "Obese"
edit.metadata2$BMI.class[edit.metadata2$BMI>=40] <- "Extreme Obesity"


#Create Obese Yes/No groups
edit.metadata2$obese[edit.metadata2$BMI.class=="Normal" | edit.metadata2$BMI.class=="Overweight"] <- "No"
edit.metadata2$obese[edit.metadata2$BMI.class=="Obese" | edit.metadata2$BMI.class=="Extreme Obesity"] <- "Yes"

#Create BMI class 2 groups
edit.metadata2$BMI.class2[edit.metadata2$BMI<=24] <- "Normal"
edit.metadata2$BMI.class2[edit.metadata2$BMI>24 & edit.metadata2$BMI<30] <- "Overweight"
edit.metadata2$BMI.class2[edit.metadata2$BMI>=30] <- "Obese"

#Get paitent demographics to be tested

bmi <- edit.metadata2$BMI
obese <- factor(edit.metadata2$obese)

######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

rossH <- wilcox.test(H ~ obese) #P-value=0.6667
rossS <- wilcox.test(S ~ obese) #P-value=0.9559
rossJ <- wilcox.test(J ~ obese) #P-value=0.5462

#B and F tests against obesity
bacter <- phyla.table.rel.abund$p__Bacteroidetes
firm <- phyla.table.rel.abund$p__Firmicutes
BFratio <- bacter/firm

rossBacter <- wilcox.test(bacter ~ obese) #P-value=0.0867
rossFirm <- wilcox.test(firm ~ obese) #P-value=0.1576
rossBF <- wilcox.test(BFratio ~ obese) #P-value=0.096

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
ross2 <- adonis(his.microb.edit ~ obese, permutations=1000)
rossPERM <- ross2$aov.tab
#PERMANOVA=0.88259, pseudo-F=0.5025

###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
RosstotalN <- length(orderedHCat)
RossHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
RossHighShannonGroup <- as.data.frame(table(test3[c(1:RossHHighTotal), 2]))
RossLowShannonGroup <- as.data.frame(table(test3[c((RossHHighTotal + 1):RosstotalN), 2]))
group1 <- c(RossHighShannonGroup[2, 2], RossHighShannonGroup[1, 2])
group2 <- c(RossLowShannonGroup[2, 2], RossLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 18 and non-obese = 13
#Group2 (Lower than median), obese = 20 and non-obese = 12

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

rossHEpi <- epi.2by2(r.test, method="cohort.count")
rossHMassoc <- rossHEpi$massoc
rossHRR <- rossHMassoc$RR.strata.score
rossHRRsig <- rossHMassoc$chisq.strata
## Risk Ratio = 1.00
## CI = 0.65, 1.54
## p-value = 1.00


##Run the RR for B/F ratio
BFRatio <- as.data.frame(BFratio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFratio <= median(BFratio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
RossBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
RossHighBFGroup <- as.data.frame(table(test4[c(1:RossBFHighTotal), 2]))
RossLowBFGroup <- as.data.frame(table(test4[c((RossBFHighTotal + 1):RosstotalN), 2]))
group1 <- c(RossHighBFGroup[2, 2], RossHighBFGroup[1, 2])
group2 <- c(RossLowBFGroup[2, 2], RossLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 22 and non-obese = 9
#Group2 (Lower than median), obese = 16 and non-obese = 16


r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

rossBFEpi <- epi.2by2(r.test, method="cohort.count")
rossBFMassoc <- rossBFEpi$massoc
rossBFRR <- rossBFMassoc$RR.strata.score
rossBFRRsig <- rossBFMassoc$chisq.strata
## Risk Ratio = 0.70
## CI = 0.45, 1.10
## p-value = 0.11

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
edit.metadata2$obese.num[edit.metadata2$obese=="No"] <- 0
edit.metadata2$obese.num[edit.metadata2$obese=="Yes"] <- 1
obese <- factor(edit.metadata2$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), his.microb.edit)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
rossAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 5 Measures, AUCopt = 0.810049

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

RossZH <- scale(H)
RossZLogBF <- scale(log(BFratio))
RossBMI <- bmi

###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(edit.metadata2)))
meanAge <- c(meanAge, mean(edit.metadata2$age_at_visit))
SDAge <- c(SDAge, sd(edit.metadata2$age_at_visit))
temporary <- table(edit.metadata2$sex)
males <- unname(c(males, temporary[names(temporary) == "M"]))
females <- unname(c(females, temporary[names(temporary) == "F"]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, 0))
meanBMI <- c(meanBMI, mean(edit.metadata2$BMI))
SDBMI <- c(SDBMI, sd(edit.metadata2$BMI))
minBMI <- c(minBMI, min(edit.metadata2$BMI))
maxBMI <- c(maxBMI, max(edit.metadata2$BMI))

Ross <- c(rossBacter$p.value, rossFirm$p.value, rossBF$p.value, 
          rossH$p.value, rossS$p.value, rossJ$p.value, rossPERM[1,6])
overallPTable <- rbind(overallPTable, Ross)

tposH <- c(tposH, RossLowShannonGroup[2, 2])
tnegH <- c(tnegH, RossLowShannonGroup[1, 2])
cposH <- c(cposH, RossHighShannonGroup[2, 2])
cnegH <- c(cnegH, RossHighShannonGroup[1, 2])
RRH <- c(RRH, rossHRR[1,1])
lowH <- c(lowH, rossHRR[1,2])
highH <- c(highH, rossHRR[1,3])


tposBF <- c(tposBF, RossLowBFGroup[2, 2])
tnegBF <- c(tnegBF, RossLowBFGroup[1, 2])
cposBF <- c(cposBF, RossHighBFGroup[2, 2])
cnegBF <- c(cnegBF, RossHighBFGroup[1, 2])
RRBF <- c(RRBF, rossBFRR[1,1])
lowBF <- c(lowBF, rossBFRR[1,2])
highBF <- c(highBF, rossBFRR[1,3])

RossData <- as.data.frame(cbind(RossZH, RossZLogBF, RossBMI))
RossData$Study <- "Ross"
colnames(RossData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, RossData)

rm(RossLowShannonGroup, RossHighShannonGroup, rossHRR, rossBacter, 
   rossFirm, rossBF, rossH, rossS, rossJ, rossPERM, rossHEpi, 
   RossZH, RossZLogBF, RossBMI, ross2)

############## GOODRICH ##################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


setwd("C:/users/marc/Desktop/obesity2/twinsUK")

#Read in and match metadata to microbiome data
metadata <- read.csv("TwinsUKStudy2.csv")
rownames(metadata) <- metadata[, 7]
shared.data <- read.table("combined.tx.1.subsample.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
microbiome <- shared.data[, -c(1:3)]
keep  <- which(metadata$body_mass_index_s != "<not provided>")

microbiome<- microbiome[keep, ]
metadata <- metadata[keep, ]

rm(shared.data, keep)


#generate alpha diversity measures with vegan
H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity <- as.data.frame(cbind(H, S, J))
alpha.test <- alpha.diversity

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.table("phyla.txt", header=T)
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(microbiome)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- microbiome
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Chloroflexi", "Deferribacteres", "Elusimicrobia", "Fusobacteria", "Gemmatimonadetes", "Lentisphaerae", "OD1", "Planctomycetes", "Spirochaetes", "SR1", "Synergistetes", "Tenericutes")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4:6, 8:12, 14:17)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100

#Create BMI groups
metadata$BMI.class[as.numeric(as.character(metadata$body_mass_index_s))<=24] <- "Normal"
metadata$BMI.class[as.numeric(as.character(metadata$body_mass_index_s))>24 & 
                     as.numeric(as.character(metadata$body_mass_index_s))<30] <- "Overweight"
metadata$BMI.class[as.numeric(as.character(metadata$body_mass_index_s))>=30 & 
                     as.numeric(as.character(metadata$body_mass_index_s))<40] <- "Obese"
metadata$BMI.class[as.numeric(as.character(metadata$body_mass_index_s))>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
metadata$obese[metadata$BMI.class=="Normal" | metadata$BMI.class=="Overweight"] <- "No"
metadata$obese[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Yes"

#Create a column with obese and extreme obese as single entity
metadata$BMIclass2[metadata$BMI.class=="Normal"] <- "Normal"
metadata$BMIclass2[metadata$BMI.class=="Overweight"] <- "Overweight"                     
metadata$BMIclass2[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Obese"

#Get paitent demographics to be tested
bmi <- as.numeric(as.character(metadata$body_mass_index_s))
obese <- factor(metadata$obese)



###########################################################################
######### First Level Analysis & Alpha Diversity with BMI #################
###########################################################################

##Test BMI versus alpha diversity and phyla

goodrichH <- wilcox.test(H ~ obese) #P-value=0.6671
goodrichS <- wilcox.test(S ~ obese) #P-value=0.9073
goodrichJ <- wilcox.test(J ~ obese) #P-value=0.6319

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

goodrichBacter <- wilcox.test(bacter ~ obese) #P-value=0.7008
goodrichFirm <- wilcox.test(firm ~ obese) #P-value=0.997
goodrichBF <- wilcox.test(BFratio ~ obese) #P-value=0.7445

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
goodrich2 <- adonis(microbiome ~ obese, permutations=1000)
goodrichPERM <- goodrich2$aov.tab
#PERMANOVA=0.004995, pseudo-F=2.9308

###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})


##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
GoodtotalN <- length(orderedHCat)
GoodHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
GoodHighShannonGroup <- as.data.frame(table(test3[c(1:GoodHHighTotal), 2]))
GoodLowShannonGroup <- as.data.frame(table(test3[c((GoodHHighTotal + 1):GoodtotalN), 2]))
group1 <- c(GoodHighShannonGroup[2, 2], GoodHighShannonGroup[1, 2])
group2 <- c(GoodLowShannonGroup[2, 2], GoodLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 47 and non-obese = 206
#Group2 (Lower than median), obese = 56 and non-obese = 198

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

goodrichHEpi <- epi.2by2(r.test, method="cohort.count")
goodrichHMassoc <- goodrichHEpi$massoc
goodrichHRR <- goodrichHMassoc$RR.strata.score
goodrichHRRsig <- goodrichHMassoc$chisq.strata
## Risk Ratio = 1.19
## CI = 0.84, 1.68
## p-value = 0.332

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
GoodBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
GoodHighBFGroup <- as.data.frame(table(test4[c(1:GoodBFHighTotal), 2]))
GoodLowBFGroup <- as.data.frame(table(test4[c((GoodBFHighTotal + 1):GoodtotalN), 2]))
group1 <- c(GoodHighBFGroup[2, 2], GoodHighBFGroup[1, 2])
group2 <- c(GoodLowBFGroup[2, 2], GoodLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 52 and non-obese = 201
#Group2 (Lower than median), obese = 51 and non-obese = 203

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

goodrichBFEpi <- epi.2by2(r.test, method="cohort.count")
goodrichBFMassoc <- goodrichBFEpi$massoc
goodrichBFRR <- goodrichBFMassoc$RR.strata.score
goodrichBFRRsig <- goodrichBFMassoc$chisq.strata
## Risk Ratio = 0.98
## CI = 0.69, 1.38
## p-value = 0.894

###########################################################################
############ Classification using AUCRF####################################
###########################################################################

#Create Obese.num group
metadata$obese.num[metadata$obese=="No"] <- 0
metadata$obese.num[metadata$obese=="Yes"] <- 1
obese <- factor(metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
goodrichAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 8 Measures, AUCopt = 0.6767159

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

GoodrichZH <- scale(H)
GoodrichZLogBF <- scale(log(BFratio))
GoodrichBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, mean(metadata$age_s))
SDAge <- c(SDAge, sd(metadata$age_s))
temporary <- table(metadata$sex_s)
males <- unname(c(males, temporary[names(temporary) == 47]))
females <- unname(c(females, temporary[names(temporary) == 48]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, NA))
meanBMI <- c(meanBMI, mean(as.numeric(as.character(metadata$body_mass_index_s))))
SDBMI <- c(SDBMI, sd(as.numeric(as.character(metadata$body_mass_index_s))))
minBMI <- c(minBMI, min(as.numeric(as.character(metadata$body_mass_index_s))))
maxBMI <- c(maxBMI, max(as.numeric(as.character(metadata$body_mass_index_s))))

Goodrich <- c(goodrichBacter$p.value, goodrichFirm$p.value, 
              goodrichBF$p.value, goodrichH$p.value, goodrichS$p.value, 
              goodrichJ$p.value, goodrichPERM[1,6])

overallPTable <- rbind(overallPTable, Goodrich)

tposH <- c(tposH, GoodLowShannonGroup[2, 2])
tnegH <- c(tnegH, GoodLowShannonGroup[1, 2])
cposH <- c(cposH, GoodHighShannonGroup[2, 2])
cnegH <- c(cnegH, GoodHighShannonGroup[1, 2])
RRH <- c(RRH, goodrichHRR[1,1])
lowH <- c(lowH, goodrichHRR[1,2])
highH <- c(highH, goodrichHRR[1,3])

tposBF <- c(tposBF, GoodLowBFGroup[2, 2])
tnegBF <- c(tnegBF, GoodLowBFGroup[1, 2])
cposBF <- c(cposBF, GoodHighBFGroup[2, 2])
cnegBF <- c(cnegBF, GoodHighBFGroup[1, 2])
RRBF <- c(RRBF, goodrichBFRR[1,1])
lowBF <- c(lowBF, goodrichBFRR[1,2])
highBF <- c(highBF, goodrichBFRR[1,3])

GoodrichData <- as.data.frame(cbind(GoodrichZH, GoodrichZLogBF, GoodrichBMI))
GoodrichData$Study <- "Goodrich"
colnames(GoodrichData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, GoodrichData)

rm(GoodLowShannonGroup, GoodHighShannonGroup, goodrichHRR, goodrichBacter, 
   goodrichFirm, goodrichBF, goodrichH, goodrichS, goodrichJ, 
   goodrichPERM, goodrichHEpi, GoodrichZH, GoodrichZLogBF, 
   GoodrichBMI, goodrich2, Ross, Goodrich)

########### ESCOBAR ######################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/columbian/")

columbian.microb <- read.table("EscobarPhylotypeSub.shared", header=T)
metadata <- read.csv("columbian_dataset.csv")

#Organize the microbiome shared data
rownames(columbian.microb) <- columbian.microb[, 2]
columbian.microb <- columbian.microb[, -c(1:3)]

#Organize the metadata
rownames(metadata) <- metadata[, 5]

#Get only data that we are interested in
edit.metadata <- metadata[, -c(1:7, 9, 13:15, 16:52)]

#Sort metadata into the same order as the microbiome data and get sex information
order1 <- rownames(columbian.microb)
edit.metadata2 <- edit.metadata[order1, ]
edit.metadata2$sex[edit.metadata2$description_s == "Adequate weight male stool sample" | 
                     edit.metadata2$description_s == "Obese male stool sample" | 
                     edit.metadata2$description_s == "Overweight male stool sample"] <- "M"
edit.metadata2$sex[edit.metadata2$description_s == "Adequate weight female stool sample" | 
                     edit.metadata2$description_s == "Obese female stool sample" | 
                     edit.metadata2$description_s == "Overweight female stool sample"] <- "F"


#Get alpha diversity of the samples
H <- diversity(columbian.microb)
S <- specnumber(columbian.microb)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.table("taxonomyKey.txt")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(columbian.microb)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- columbian.microb
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("p__", "p__Cyanobacteria", "p__Elusimicrobia", 
                                           "p__Fusobacteria", "p__Lentisphaerae", 
                                           "p__Spirochaetes", "p__Synergistetes", 
                                           "p__Tenericutes", "p__TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4:5, 7:8, 10:13)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100

#Create Obese Yes/No groups
edit.metadata2$obese[edit.metadata2$description_s=="Adequate weight male stool sample" | edit.metadata2$description_s=="Overweight male stool sample" | edit.metadata2$description_s=="Adequate weight female stool sample" | edit.metadata2$description_s=="Overweight female stool sample"] <- "No"
edit.metadata2$obese[edit.metadata2$description_s=="Obese male stool sample" | edit.metadata2$description_s=="Obese female stool sample"] <- "Yes"

#Create BMI groups
edit.metadata2$BMI.class[edit.metadata2$body_mass_index_s<=24] <- "Normal"
edit.metadata2$BMI.class[edit.metadata2$body_mass_index_s>24 & edit.metadata2$body_mass_index_s<30] <- "Overweight"
edit.metadata2$BMI.class[edit.metadata2$body_mass_index_s>=30 & edit.metadata2$body_mass_index_s<40] <- "Obese"
edit.metadata2$BMI.class[edit.metadata2$body_mass_index_s>=40] <- "Extreme Obesity"

#Create BMI groups2
edit.metadata2$BMI.class2[edit.metadata2$body_mass_index_s<=24] <- "Normal"
edit.metadata2$BMI.class2[edit.metadata2$body_mass_index_s>24 & edit.metadata2$body_mass_index_s<30] <- "Overweight"
edit.metadata2$BMI.class2[edit.metadata2$body_mass_index_s>=30] <- "Obese"

#Get paitent demographics data to be tested
bmi <- edit.metadata2$body_mass_index_s
obese <- factor(edit.metadata2$obese)

######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

escobarH <- wilcox.test(H ~ obese) #P-value=0.5588
escobarS <- wilcox.test(S ~ obese) #P-value=0.3667
escobarJ <- wilcox.test(J ~ obese) #P-value=0.4223

#B and F tests against obesity
bacter <- phyla.table.rel.abund$p__Bacteroidetes
firm <- phyla.table.rel.abund$p__Firmicutes
BFratio <- bacter/firm

escobarBacter <- wilcox.test(bacter ~ obese) #P-value=0.02442
escobarFirm <- wilcox.test(firm ~ obese) #P-value=0.09958
escobarBF <- wilcox.test(BFratio ~ obese) #P-value=0.04392

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
escobar2 <- adonis(columbian.microb ~ obese, permutations=1000)
escobarPERM <- escobar2$aov.tab
#PERMANOVA=0.07892, pseudo-F=1.6925

###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
EscototalN <- length(orderedHCat)
EscoHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
EscoHighShannonGroup <- as.data.frame(table(test3[c(1:EscoHHighTotal), 2]))
EscoLowShannonGroup <- as.data.frame(table(test3[c((EscoHHighTotal + 1):EscototalN), 2]))
group1 <- c(EscoHighShannonGroup[2, 2], EscoHighShannonGroup[1, 2])
group2 <- c(EscoLowShannonGroup[2, 2], EscoLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 4 and non-obese = 11
#Group2 (Lower than median), obese = 6 and non-obese = 9

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

escobarHEpi <- epi.2by2(r.test, method="cohort.count")
escobarHMassoc <- escobarHEpi$massoc
escobarHRR <- escobarHMassoc$RR.strata.score
escobarHRRsig <- escobarHMassoc$chisq.strata
## Risk Ratio = 1.50
## CI = 0.53, 4.26
## p-value = 0.439

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$p__Bacteroidetes
Firm = phyla.table.rel.abund$p__Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
EscoBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
EscoHighBFGroup <- as.data.frame(table(test4[c(1:EscoBFHighTotal), 2]))
EscoLowBFGroup <- as.data.frame(table(test4[c((EscoBFHighTotal + 1):EscototalN), 2]))
group1 <- c(EscoHighBFGroup[2, 2], EscoHighBFGroup[1, 2])
group2 <- c(EscoLowBFGroup[2, 2], EscoLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 7 and non-obese = 8
#Group2 (Lower than median), obese = 3 and non-obese = 12

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

escobarBFEpi <- epi.2by2(r.test, method="cohort.count")
escobarBFMassoc <- escobarBFEpi$massoc
escobarBFRR <- escobarBFMassoc$RR.strata.score
escobarBFRRsig <- escobarBFMassoc$chisq.strata
## Risk Ratio = 0.25
## CI = 0.06, 0.99
## p-value = 0.02

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
edit.metadata2$obese.num[edit.metadata2$obese=="No"] <- 0
edit.metadata2$obese.num[edit.metadata2$obese=="Yes"] <- 1
obese <- factor(edit.metadata2$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), columbian.microb)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
escobarAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 3 Measures, AUCopt = 0.95

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

EscobarZH <- scale(H)
EscobarZLogBF <- scale(log(BFratio))
EscobarBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(edit.metadata2)))
meanAge <- c(meanAge, mean(edit.metadata2$age_s))
SDAge <- c(SDAge, sd(edit.metadata2$age_s))
temporary <- table(edit.metadata2$sex)
males <- unname(c(males, temporary[names(temporary) == "M"]))
females <- unname(c(females, temporary[names(temporary) == "F"]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, 0))
meanBMI <- c(meanBMI, mean(edit.metadata2$body_mass_index_s))
SDBMI <- c(SDBMI, sd(edit.metadata2$body_mass_index_s))
minBMI <- c(minBMI, min(edit.metadata2$body_mass_index_s))
maxBMI <- c(maxBMI, max(edit.metadata2$body_mass_index_s))


Escobar <- c(escobarBacter$p.value, escobarFirm$p.value, 
             escobarBF$p.value, escobarH$p.value, escobarS$p.value, 
             escobarJ$p.value, escobarPERM[1,6])

overallPTable <- rbind(overallPTable, Escobar)

tposH <- c(tposH, EscoLowShannonGroup[2, 2])
tnegH <- c(tnegH, EscoLowShannonGroup[1, 2])
cposH <- c(cposH, EscoHighShannonGroup[2, 2])
cnegH <- c(cnegH, EscoHighShannonGroup[1, 2])
RRH <- c(RRH, escobarHRR[1,1])
lowH <- c(lowH, escobarHRR[1,2])
highH <- c(highH, escobarHRR[1,3])

tposBF <- c(tposBF, EscoLowBFGroup[2, 2])
tnegBF <- c(tnegBF, EscoLowBFGroup[1, 2])
cposBF <- c(cposBF, EscoHighBFGroup[2, 2])
cnegBF <- c(cnegBF, EscoHighBFGroup[1, 2])
RRBF <- c(RRBF, escobarBFRR[1,1])
lowBF <- c(lowBF, escobarBFRR[1,2])
highBF <- c(highBF, escobarBFRR[1,3])

EscobarData <- as.data.frame(cbind(EscobarZH, EscobarZLogBF, EscobarBMI))
EscobarData$Study <- "Escobar"
colnames(EscobarData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, EscobarData)

rm(EscoLowShannonGroup, EscoHighShannonGroup, escobarHRR, escobarBacter, 
   escobarFirm, escobarBF, escobarH, escobarS, escobarJ, 
   escobarPERM, escobarHEpi, EscobarZH, EscobarZLogBF, 
   EscobarBMI, escobar2, Escobar)


######## ZUPANCIC #########################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


setwd("C:/users/marc/Desktop/obesity2/Amish/updatedGOOD")

#Read in and match metadata to microbiome data
metadata <- read.csv("amish_obesity_table2.csv")
metadata2 <- read.csv("amish.metadata.csv")
test <- metadata[!duplicated(metadata$submitted_sample_id_s), ]
rownames(test) <- test[, 4]
test2 <- metadata2[!duplicated(metadata2$SUBJID), ]
rownames(test2) <- test2[, 1]
shared.data <- read.table("ZupancicPhylotypeSub.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
shared.data <- shared.data[, -c(1:3)]
keep  <- rownames(shared.data)

test3 <- test[keep, ]
good.metadata <- test3
test4 <- test2[keep, ]
good.metadata2 <- test4

microbiome <- shared.data[keep, ]
metadata <- cbind(good.metadata, good.metadata2)

rm(good.metadata, shared.data, keep, test, test2, test3, test4, good.metadata2, metadata2)

test <- metadata[!duplicated(metadata$submitted_subject_id_s), ]
keep  <- rownames(test)
metadata <- test
microbiome <- microbiome[keep, ]

rm(test, keep)

metadata <- metadata[complete.cases(metadata), ]
keep <- rownames(metadata)
microbiome <- microbiome[keep, ]

rm(keep)

#generate alpha diversity measures with vegan
library(vegan)
H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))
s1.alpha.diversity <- as.data.frame(select.alpha.diversity)
alpha.test <- s1.alpha.diversity

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.table("taxonomyKey.txt", header=T)
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(microbiome)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- microbiome
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("p__Elusimicrobia", 
                                           "p__Fusobacteria", 
                                           "p__Synergistetes")], 1, sum)
phyla.table <- phyla.table[, -c(3, 5, 7)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:6)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100

#Create BMI groups
metadata$BMI.class[metadata$BMI<=24] <- "Normal"
metadata$BMI.class[metadata$BMI>24 & metadata$BMI<30] <- "Overweight"
metadata$BMI.class[metadata$BMI>=30 & metadata$BMI<40] <- "Obese"
metadata$BMI.class[metadata$BMI>=40] <- "Extreme Obesity"


#Create Obese Yes/No groups
metadata$obese[metadata$BMI.class=="Normal" | metadata$BMI.class=="Overweight"] <- "No"
metadata$obese[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Yes"

#Create BMI class 2 groups
metadata$BMI.class2[metadata$BMI<=24] <- "Normal"
metadata$BMI.class2[metadata$BMI>24 & metadata$BMI<30] <- "Overweight"
metadata$BMI.class2[metadata$BMI>=30] <- "Obese"

#create groups to be used
bmi <- metadata$BMI
obese <- factor(metadata$obese)


####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

zupancicH <- wilcox.test(H ~ obese) #P-value=0.9207
zupancicS <- wilcox.test(S ~ obese) #P-value=0.5007
zupancicJ <- wilcox.test(J ~ obese) #P-value=0.4922

#B and F tests against obesity
bacter <- phyla.table.rel.abund$p__Bacteroidetes
firm <- phyla.table.rel.abund$p__Firmicutes
BFratio <- bacter/firm

zupancicBacter <- wilcox.test(bacter ~ obese) #P-value=0.8712
zupancicFirm <- wilcox.test(firm ~ obese) #P-value=0.5536
zupancicBF <- wilcox.test(BFratio ~ obese) #P-value=0.9491

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
zupancic2 <- adonis(microbiome ~ obese, permutations=1000)
zupancicPERM <- zupancic2$aov.tab
#PERMANOVA=0.8012, pseudo-F=0.48248

###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})


##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
ZupatotalN <- length(orderedHCat)
ZupaHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
ZupaHighShannonGroup <- as.data.frame(table(test3[c(1:ZupaHHighTotal), 2]))
ZupaLowShannonGroup <- as.data.frame(table(test3[c((ZupaHHighTotal + 1):ZupatotalN), 2]))
group1 <- c(ZupaHighShannonGroup[2, 2], ZupaHighShannonGroup[1, 2])
group2 <- c(ZupaLowShannonGroup[2, 2], ZupaLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 31 and non-obese = 72
#Group2 (Lower than median), obese = 47 and non-obese = 57

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

zupancicHEpi <- epi.2by2(r.test, method="cohort.count")
zupancicHMassoc <- zupancicHEpi$massoc
zupancicHRR <- zupancicHMassoc$RR.strata.score
zupancicHRRsig <- zupancicHMassoc$chisq.strata
## Risk Ratio = 1.09
## CI = 0.75, 1.58
## p-value = 0.658


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$p__Bacteroidetes
Firm = phyla.table.rel.abund$p__Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
ZupaBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
ZupaHighBFGroup <- as.data.frame(table(test4[c(1:ZupaBFHighTotal), 2]))
ZupaLowBFGroup <- as.data.frame(table(test4[c((ZupaBFHighTotal + 1):ZupatotalN), 2]))
group1 <- c(ZupaHighBFGroup[2, 2], ZupaHighBFGroup[1, 2])
group2 <- c(ZupaLowBFGroup[2, 2], ZupaLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 41 and non-obese = 62
#Group2 (Lower than median), obese = 37 and non-obese = 67

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

zupancicBFEpi <- epi.2by2(r.test, method="cohort.count")
zupancicBFMassoc <- zupancicBFEpi$massoc
zupancicBFRR <- zupancicBFMassoc$RR.strata.score
zupancicBFRRsig <- zupancicBFMassoc$chisq.strata
## Risk Ratio = 0.97
## CI = 0.67, 1.41
## p-value = 0.883

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
metadata$obese.num[metadata$obese=="No"] <- 0
metadata$obese.num[metadata$obese=="Yes"] <- 1
obese <- factor(metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
zupancicAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
#zupAUC <- fit$`OOB-AUCopt`
# list of 3 Measures, AUCopt = 0.5746806


###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

ZupancicZH <- scale(H)
ZupancicZLogBF <- scale(log(BFratio))
ZupancicBMI <- bmi

###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, mean(as.numeric(as.character(metadata$AGE))))
SDAge <- c(SDAge, sd(as.numeric(as.character(metadata$AGE))))
temporary <- table(metadata$sex_s)
males <- unname(c(males, temporary[names(temporary) == "male"]))
females <- unname(c(females, temporary[names(temporary) == "female"]))
#temporary <- table(demographics$White)
ancestry <- unname(c(ancestry, 100))
meanBMI <- c(meanBMI, mean(metadata$BMI))
SDBMI <- c(SDBMI, sd(metadata$BMI))
minBMI <- c(minBMI, min(metadata$BMI))
maxBMI <- c(maxBMI, max(metadata$BMI))

Zupancic <- c(zupancicBacter$p.value, zupancicFirm$p.value, 
              zupancicBF$p.value, zupancicH$p.value, zupancicS$p.value, 
              zupancicJ$p.value, zupancicPERM[1,6])

overallPTable <- rbind(overallPTable, Zupancic)

tposH <- c(tposH, ZupaLowShannonGroup[2, 2])
tnegH <- c(tnegH, ZupaLowShannonGroup[1, 2])
cposH <- c(cposH, ZupaHighShannonGroup[2, 2])
cnegH <- c(cnegH, ZupaHighShannonGroup[1, 2])
RRH <- c(RRH, zupancicHRR[1,1])
lowH <- c(lowH, zupancicHRR[1,2])
highH <- c(highH, zupancicHRR[1,3])

tposBF <- c(tposBF, ZupaLowBFGroup[2, 2])
tnegBF <- c(tnegBF, ZupaLowBFGroup[1, 2])
cposBF <- c(cposBF, ZupaHighBFGroup[2, 2])
cnegBF <- c(cnegBF, ZupaHighBFGroup[1, 2])
RRBF <- c(RRBF, zupancicBFRR[1,1])
lowBF <- c(lowBF, zupancicBFRR[1,2])
highBF <- c(highBF, zupancicBFRR[1,3])

ZupancicData <- as.data.frame(cbind(ZupancicZH, ZupancicZLogBF, ZupancicBMI))
ZupancicData$Study <- "Zupancic"
colnames(ZupancicData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, ZupancicData)

rm(ZupaLowShannonGroup, ZupaHighShannonGroup, zupancicHRR, zupancicBacter, 
   zupancicFirm, zupancicBF, zupancicH, zupancicS, zupancicJ, 
   zupancicPERM, zupancicHEpi, ZupancicZH, ZupancicZLogBF, 
   ZupancicBMI, zupancic2, Zupancic)


############### HMP HMP HMP ##############################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

##Read in relevant data

setwd("C:/users/marc/Desktop/obesity2/HMP.analysis")
microbiome <- read.table("Stool.an.0.03.subsample.shared", header=T)
meta.cat <- read.table("categorical.metadata", header=T)
meta.cont <- read.table("continuous.metadata", header=T)

#Only interested in first visit microbiome data
#need to create microbiome data set for only that

test <- microbiome[grep("\\.01\\.", microbiome$Group), ] #selecting by .01.

microb.rownames <- gsub("([0-9]+).*", "\\1", test$Group) #extracting only numeric before first "."
rownames(test) <- microb.rownames

#Get rid of information not used in downstram analysis
test <- test[, -c(1:3)]

#Subset data to the microbiome total n
select.meta.cat <- meta.cat[microb.rownames, ]
select.meta.cont <- meta.cont[microb.rownames, ]

#generate alpha diversity measures with vegan
H <- diversity(test)
S <- specnumber(test)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))
alpha.test <- select.alpha.diversity

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(test)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- test
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Deinococcus-Thermus", "Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "Tenericutes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4, 6, 7, 9, 10:12)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]


#Create Obese Yes/No groups
select.meta.cat$obese[select.meta.cat$BMI_C=="normal" | select.meta.cat$BMI_C=="overweight"] <- "No"
select.meta.cat$obese[select.meta.cat$BMI_C=="obese" | select.meta.cat$BMI_C=="extreme obesity"] <- "Yes"

#Create BMI groups2
select.meta.cat$BMI.class2[select.meta.cont$DTPBMI<=24] <- "Normal"
select.meta.cat$BMI.class2[select.meta.cont$DTPBMI>24 & select.meta.cont$DTPBMI<30] <- "Overweight"
select.meta.cat$BMI.class2[select.meta.cont$DTPBMI>=30] <- "Obese"

#Get paitent demographics data to be tested
bmi <- select.meta.cont$DTPBMI
obese <- factor(select.meta.cat$obese)

######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

HMPH <- wilcox.test(H ~ obese) #P-value=0.5772
HMPS <- wilcox.test(S ~ obese) #P-value=0.9654
HMPJ <- wilcox.test(J ~ obese) #P-value=0.4667


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
HMP2 <- adonis(test ~ obese, permutations=1000)
HMPPERM <- HMP2$aov.tab
#PERMANOVA=0.8112, pseudo-F=0.7024


######################################################################################## Relative Risk ###############################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})


##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
HMPtotalN <- length(orderedHCat)
HMPHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
HMPHighShannonGroup <- as.data.frame(table(test3[c(1:HMPHHighTotal), 2]))
HMPLowShannonGroup <- as.data.frame(table(test3[c((HMPHHighTotal + 1):HMPtotalN), 2]))
group1 <- c(HMPHighShannonGroup[2, 2], HMPHighShannonGroup[1, 2])
group2 <- c(HMPLowShannonGroup[2, 2], HMPLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 10 and non-obese = 118
#Group2 (Lower than median), obese = 16 and non-obese = 112

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

HMPHEpi <- epi.2by2(r.test, method="cohort.count")
HMPHMassoc <- HMPHEpi$massoc
HMPHRR <- HMPHMassoc$RR.strata.score
HMPHRRsig <- HMPHMassoc$chisq.strata
## Risk Ratio = 1.60
## CI = 0.75, 3.39
## p-value = 0.214

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
HMPBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
HMPHighBFGroup <- as.data.frame(table(test4[c(1:HMPBFHighTotal), 2]))
HMPLowBFGroup <- as.data.frame(table(test4[c((HMPBFHighTotal + 1):HMPtotalN), 2]))
group1 <- c(HMPHighBFGroup[2, 2], HMPHighBFGroup[1, 2])
group2 <- c(HMPLowBFGroup[2, 2], HMPLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 13 and non-obese = 115
#Group2 (Lower than median), obese = 13 and non-obese = 115

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")


HMPBFEpi <- epi.2by2(r.test, method="cohort.count")
HMPBFMassoc <- HMPBFEpi$massoc
HMPBFRR <- HMPBFMassoc$RR.strata.score
HMPBFRRsig <- HMPBFMassoc$chisq.strata
## Risk Ratio = 1.00
## CI = 0.48, 2.07
## p-value = 1.00

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
select.meta.cat$obese.num[select.meta.cat$obese=="No"] <- 0
select.meta.cat$obese.num[select.meta.cat$obese=="Yes"] <- 1
obese <- factor(select.meta.cat$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), test)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
#set.seed(3)
#HMPAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 6 Measures, AUCopt = 0.7031773

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

HMPZH <- scale(H)
HMPZLogBF <- scale(log(BFratio))
HMPBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(select.meta.cat)))
meanAge <- c(meanAge, mean(select.meta.cont$AGEENR))
SDAge <- c(SDAge, sd(select.meta.cont$AGEENR))
temporary <- table(select.meta.cat$GENDER_C)
males <- unname(c(males, temporary[names(temporary) == "Male"]))
females <- unname(c(females, temporary[names(temporary) == "Female"]))
temporary <- table(select.meta.cat$WHITE_C)
ancestry <- unname(c(ancestry, temporary[names(temporary) == "Yes"] / sum(temporary)))
meanBMI <- c(meanBMI, mean(select.meta.cont$DTPBMI))
SDBMI <- c(SDBMI, sd(select.meta.cont$DTPBMI))
minBMI <- c(minBMI, min(select.meta.cont$DTPBMI))
maxBMI <- c(maxBMI, max(select.meta.cont$DTPBMI))

HMP <- c(HMPBacter$p.value, HMPFirm$p.value, 
         HMPBF$p.value, HMPH$p.value, HMPS$p.value, 
         HMPJ$p.value, HMPPERM[1,6])

overallPTable <- rbind(overallPTable, HMP)

tposH <- c(tposH, HMPLowShannonGroup[2, 2])
tnegH <- c(tnegH, HMPLowShannonGroup[1, 2])
cposH <- c(cposH, HMPHighShannonGroup[2, 2])
cnegH <- c(cnegH, HMPHighShannonGroup[1, 2])
RRH <- c(RRH, HMPHRR[1,1])
lowH <- c(lowH, HMPHRR[1,2])
highH <- c(highH, HMPHRR[1,3])

tposBF <- c(tposBF, HMPLowBFGroup[2, 2])
tnegBF <- c(tnegBF, HMPLowBFGroup[1, 2])
cposBF <- c(cposBF, HMPHighBFGroup[2, 2])
cnegBF <- c(cnegBF, HMPHighBFGroup[1, 2])
RRBF <- c(RRBF, HMPBFRR[1,1])
lowBF <- c(lowBF, HMPBFRR[1,2])
highBF <- c(highBF, HMPBFRR[1,3])


HMPData <- as.data.frame(cbind(HMPZH, HMPZLogBF, HMPBMI))
HMPData$Study <- "HMP"
colnames(HMPData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, HMPData)

rm(HMPLowShannonGroup, HMPHighShannonGroup, HMPHRR, HMPBacter, 
   HMPFirm, HMPBF, HMPH, HMPS, HMPJ, 
   HMPPERM, HMPHEpi, HMPZH, HMPZLogBF, 
   HMPBMI, HMP2, HMP)


############## Wu Wu Wu ###################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


setwd("C:/users/marc/Desktop/obesity2/COMBO/")

microbiome <- read.table("WuPhylotypeSub.shared", header=T)
rownames(microbiome) <- microbiome$Group
microbiome <- microbiome[, -c(1:3)]
seqData <- read.csv("COMBO_data_table.csv", header=T)
seqData <- seqData[-1, ]
rownames(seqData) <- seqData$Run_s
metadata <- read.table("bmi_info.txt", header=T)

#get seqData set in line with microbiome data
namesToKeep <- rownames(microbiome)
test <- seqData[namesToKeep, ]
seqData <- test
rm(test)

#change all rownames to match sample IDs
rownames(seqData) <- seqData$submitted_subject_id_s
rownames(microbiome) <- rownames(seqData)


#Match the metadata now with the microbiome data
namesToKeep <- rownames(microbiome)
test <- metadata[namesToKeep, ]
metadata <- test
rm(test)
rm(seqData, namesToKeep)

#Get alpha diversity of the samples
H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.table("taxonomyKey.txt", header=T)
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(microbiome)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- microbiome
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("p__Cyanobacteria", 
                                           "p__Fusobacteria", 
                                           "p__Lentisphaerae", 
                                           "p__Synergistetes", 
                                           "p__Tenericutes", 
                                           "p__TM7")], 1, sum)
phyla.table <- phyla.table[, -c(3, 5:6, 8:10)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:6)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(keep, phyla.names, phyla.total, phylogenetic.info, phyla.table, phyla.good)

#Create BMI groups
metadata$BMI.class[metadata$bmi<=24] <- "Normal"
metadata$BMI.class[metadata$bmi>24 & metadata$bmi<30] <- "Overweight"
metadata$BMI.class[metadata$bmi>=30 & metadata$bmi<40] <- "Obese"
metadata$BMI.class[metadata$bmi>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
metadata$obese[metadata$BMI.class=="Normal" | metadata$BMI.class=="Overweight"] <- "No"
metadata$obese[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Yes"

#Create a column with obese and extreme obese as single entity
metadata$BMIclass2[metadata$BMI.class=="Normal"] <- "Normal"
metadata$BMIclass2[metadata$BMI.class=="Overweight"] <- "Overweight"                     
metadata$BMIclass2[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Obese"

#Get paitent demographics to be tested
bmi <- metadata$bmi
obese <- factor(metadata$obese)

###########################################################################
######### First Level Analysis & Alpha Diversity with BMI #################
###########################################################################

##Test BMI versus alpha diversity and phyla

WuH <- wilcox.test(H ~ obese) #P-value=0.8462
WuS <- wilcox.test(S ~ obese) #P-value=0.9447
WuJ <- wilcox.test(J ~ obese) #P-value=0.6981

#B and F tests against obesity
bacter <- phyla.table.rel.abund$p__Bacteroidetes
firm <- phyla.table.rel.abund$p__Firmicutes
BFratio <- bacter/firm

WuBacter <- wilcox.test(bacter ~ obese) #P-value=0.7395
WuFirm <- wilcox.test(firm ~ obese) #P-value=0.8031
WuBF <- wilcox.test(BFratio ~ obese) #P-value=0.9118

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
Wu2 <- adonis(microbiome ~ obese, permutations=1000)
WuPERM <- Wu2$aov.tab
#PERMANOVA=0.0.4735, pseudo-F=0.89514


###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})


##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
WutotalN <- length(orderedHCat)
WuHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
WuHighShannonGroup <- as.data.frame(table(test3[c(1:WuHHighTotal), 2]))
WuLowShannonGroup <- as.data.frame(table(test3[c((WuHHighTotal + 1):WutotalN), 2]))
group1 <- c(WuHighShannonGroup[2, 2], WuHighShannonGroup[1, 2])
group2 <- c(WuLowShannonGroup[2, 2], WuLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 3 and non-obese = 26
#Group2 (Lower than median), obese = 2 and non-obese = 27

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

WuHEpi <- epi.2by2(r.test, method="cohort.count")
WuHMassoc <- WuHEpi$massoc
WuHRR <- WuHMassoc$RR.strata.score
WuHRRsig <- WuHMassoc$chisq.strata
## Risk Ratio = 1.50
## CI = 0.27, 8.32
## p-value = 0.64

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$p__Bacteroidetes
Firm = phyla.table.rel.abund$p__Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
WuBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
WuHighBFGroup <- as.data.frame(table(test4[c(1:WuBFHighTotal), 2]))
WuLowBFGroup <- as.data.frame(table(test4[c((WuBFHighTotal + 1):WutotalN), 2]))
group1 <- c(WuHighBFGroup[2, 2], WuHighBFGroup[1, 2])
group2 <- c(WuLowBFGroup[2, 2], WuLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 3 and non-obese = 26
#Group2 (Lower than median), obese = 2 and non-obese = 27

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

WuBFEpi <- epi.2by2(r.test, method="cohort.count")
WuBFMassoc <- WuBFEpi$massoc
WuBFRR <- WuBFMassoc$RR.strata.score
WuBFRRsig <- WuBFMassoc$chisq.strata
## Risk Ratio = 1.50
## CI = 0.27, 8.32
## p-value = 0.64

###########################################################################
############ Classification using AUCRF####################################
###########################################################################

#Create Obese.num group
metadata$obese.num[metadata$obese=="No"] <- 0
metadata$obese.num[metadata$obese=="Yes"] <- 1
obese <- factor(metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
WuAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 15 Measures, AUCopt = 0.7075472

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

WuZH <- scale(H)
WuZLogBF <- scale(log(BFratio))
WuBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, mean(metadata$age))
SDAge <- c(SDAge, sd(metadata$age))
temporary <- table(metadata$sex1m2f)
males <- unname(c(males, temporary[names(temporary) == 1]))
females <- unname(c(females, temporary[names(temporary) == 2]))
#temporary <- table(select.meta.cat$WHITE_C)
ancestry <- unname(c(ancestry, NA))
meanBMI <- c(meanBMI, mean(metadata$bmi))
SDBMI <- c(SDBMI, sd(metadata$bmi))
minBMI <- c(minBMI, min(metadata$bmi))
maxBMI <- c(maxBMI, max(metadata$bmi))

Wu <- c(WuBacter$p.value, WuFirm$p.value, 
        WuBF$p.value, WuH$p.value, WuS$p.value, 
        WuJ$p.value, WuPERM[1,6])

overallPTable <- rbind(overallPTable, Wu)

tposH <- c(tposH, WuLowShannonGroup[2, 2])
tnegH <- c(tnegH, WuLowShannonGroup[1, 2])
cposH <- c(cposH, WuHighShannonGroup[2, 2])
cnegH <- c(cnegH, WuHighShannonGroup[1, 2])
RRH <- c(RRH, WuHRR[1,1])
lowH <- c(lowH, WuHRR[1,2])
highH <- c(highH, WuHRR[1,3])

tposBF <- c(tposBF, WuLowBFGroup[2, 2])
tnegBF <- c(tnegBF, WuLowBFGroup[1, 2])
cposBF <- c(cposBF, WuHighBFGroup[2, 2])
cnegBF <- c(cnegBF, WuHighBFGroup[1, 2])
RRBF <- c(RRBF, WuBFRR[1,1])
lowBF <- c(lowBF, WuBFRR[1,2])
highBF <- c(highBF, WuBFRR[1,3])

WuData <- as.data.frame(cbind(WuZH, WuZLogBF, WuBMI))
WuData$Study <- "Wu"
colnames(WuData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, WuData)

rm(WuLowShannonGroup, WuHighShannonGroup, WuHRR, WuBacter, 
   WuFirm, WuBF, WuH, WuS, WuJ, 
   WuPERM, WuHEpi, WuZH, WuZLogBF, 
   WuBMI, Wu2, Wu)


############ Aruguman ####################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/MetaHit/")

microb <- read.csv("finalized_metaHit_shared.csv")
rownames(microb) <- microb[, 1]
microb <- microb[,-1]

metadata <- read.csv("MetaHit_metadata.csv")
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]

phyla <- read.csv("phyla.csv")
rownames(phyla) <- phyla[, 1]
phyla <- phyla[, -1]

# Generate Overall Relative Abundance - normalize to a value of 100
# Needs to be done since not all measures have 100 bacteria from metaphlan2
overall <- rowSums(microb)
microb.norm <- (microb / overall) * 100

phyla.table <- read.csv("phyla.csv", header = T)
rownames(phyla.table) <- phyla.table[, 1]
phyla.table <- phyla.table[, -1]
phyla.table <- as.data.frame(t(phyla.table))

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:10)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.table, phyla.info)

#Get alpha diversity of the samples
H <- diversity(microb.norm)
S <- specnumber(microb.norm)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)
# Seems to look okay....

#Create BMI groups
metadata$BMI.class[metadata$BMI<=24] <- "Normal"
metadata$BMI.class[metadata$BMI>24 & metadata$BMI<30] <- "Overweight"
metadata$BMI.class[metadata$BMI>=30 & metadata$BMI<40] <- "Obese"
metadata$BMI.class[metadata$BMI>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
metadata$obese[metadata$BMI.class=="Normal" | metadata$BMI.class=="Overweight"] <- "No"
metadata$obese[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Yes"


#Create Obese.num group
metadata$obese.num[metadata$obese=="No"] <- 0
metadata$obese.num[metadata$obese=="Yes"] <- 1


#Create a column with obese and extreme obese as single entity
metadata$BMIclass2[metadata$BMI<=24] <- "Normal"
metadata$BMIclass2[metadata$BMI>24 & metadata$BMI<30] <- "Overweight"      
metadata$BMIclass2[metadata$BMI>=30] <- "Obese"


#Get paitent demographics to be tested
bmi <- metadata$BMI
obese <- factor(metadata$obese)

####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

ArumugamH <- wilcox.test(H ~ obese) #P-value=0.4296
ArumugamS <- wilcox.test(S ~ obese) #P-value=0.9752
ArumugamJ <- wilcox.test(J ~ obese) #P-value=0.3516

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

ArumugamBacter <- wilcox.test(bacter ~ obese) #P-value=0.2792
ArumugamFirm <- wilcox.test(firm ~ obese) #P-value=0.2381
ArumugamBF <- wilcox.test(BFratio ~ obese) #P-value=0.2562

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
Arumugam2 <- adonis(microb.norm ~ obese, permutations=1000)
ArumugamPERM <- Arumugam2$aov.tab
#PERMANOVA=0.05495, pseudo-F=1.6047


###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity


alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
ArutotalN <- length(orderedHCat)
AruHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
AruHighShannonGroup <- as.data.frame(table(test3[c(1:AruHHighTotal), 2]))
AruLowShannonGroup <- as.data.frame(table(test3[c((AruHHighTotal + 1):ArutotalN), 2]))
group1 <- c(AruHighShannonGroup[2, 2], AruHighShannonGroup[1, 2])
group2 <- c(AruLowShannonGroup[2, 2], AruLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 17 and non-obese = 25
#Group2 (Lower than median), obese = 20 and non-obese = 23

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

ArumugamHEpi <- epi.2by2(r.test, method="cohort.count")
ArumugamHMassoc <- ArumugamHEpi$massoc
ArumugamHRR <- ArumugamHMassoc$RR.strata.score
ArumugamHRRsig <- ArumugamHMassoc$chisq.strata
## Risk Ratio = 1.15
## CI = 0.71, 1.87
## p-value = 0.575


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
AruBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
AruHighBFGroup <- as.data.frame(table(test4[c(1:AruBFHighTotal), 2]))
AruLowBFGroup <- as.data.frame(table(test4[c((AruBFHighTotal + 1):ArutotalN), 2]))
group1 <- c(AruHighBFGroup[2, 2], AruHighBFGroup[1, 2])
group2 <- c(AruLowBFGroup[2, 2], AruLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 20 and non-obese = 22
#Group2 (Lower than median), obese = 17 and non-obese = 26

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

ArumugamBFEpi <- epi.2by2(r.test, method="cohort.count")
ArumugamBFMassoc <- ArumugamBFEpi$massoc
ArumugamBFRR <- ArumugamBFMassoc$RR.strata.score
ArumugamBFRRsig <- ArumugamBFMassoc$chisq.strata
## Risk Ratio = 0.83
## CI = 0.51, 1.35
## p-value = 0.452

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
metadata$obese.num[metadata$obese=="No"] <- 0
metadata$obese.num[metadata$obese=="Yes"] <- 1
obese <- factor(metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microb.norm)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
#set.seed(3)
#ArumugamAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 19 Measures, AUCopt = 0.7649212

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

ArumugamZH <- scale(H)
ArumugamZLogBF <- scale(log(BFratio))
ArumugamBMI <- bmi


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(metadata)))
meanAge <- c(meanAge, mean(metadata$Age, na.rm = TRUE))
SDAge <- c(SDAge, sd(metadata$Age, na.rm = TRUE))
temporary <- table(metadata$Sex)
males <- unname(c(males, temporary[names(temporary) == "male"]))
females <- unname(c(females, temporary[names(temporary) == "female"]))
#temporary <- table(select.meta.cat$WHITE_C)
ancestry <- unname(c(ancestry, NA))
meanBMI <- c(meanBMI, mean(metadata$BMI))
SDBMI <- c(SDBMI, sd(metadata$BMI))
minBMI <- c(minBMI, min(metadata$BMI))
maxBMI <- c(maxBMI, max(metadata$BMI))

Arumugam <- c(ArumugamBacter$p.value, ArumugamFirm$p.value, 
        ArumugamBF$p.value, ArumugamH$p.value, ArumugamS$p.value, 
        ArumugamJ$p.value, ArumugamPERM[1,6])

overallPTable <- rbind(overallPTable, Arumugam)

tposH <- c(tposH, AruLowShannonGroup[2, 2])
tnegH <- c(tnegH, AruLowShannonGroup[1, 2])
cposH <- c(cposH, AruHighShannonGroup[2, 2])
cnegH <- c(cnegH, AruHighShannonGroup[1, 2])
RRH <- c(RRH, ArumugamHRR[1,1])
lowH <- c(lowH, ArumugamHRR[1,2])
highH <- c(highH, ArumugamHRR[1,3])

tposBF <- c(tposBF, AruLowBFGroup[2, 2])
tnegBF <- c(tnegBF, AruLowBFGroup[1, 2])
cposBF <- c(cposBF, AruHighBFGroup[2, 2])
cnegBF <- c(cnegBF, AruHighBFGroup[1, 2])
RRBF <- c(RRBF, ArumugamBFRR[1,1])
lowBF <- c(lowBF, ArumugamBFRR[1,2])
highBF <- c(highBF, ArumugamBFRR[1,3])

ArumugamData <- as.data.frame(cbind(ArumugamZH, ArumugamZLogBF, ArumugamBMI))
ArumugamData$Study <- "Arumugam"
colnames(ArumugamData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, ArumugamData)

rm(AruLowShannonGroup, AruHighShannonGroup, ArumugamHRR, ArumugamBacter, 
   ArumugamFirm, ArumugamBF, ArumugamH, ArumugamS, ArumugamJ, 
   ArumugamPERM, ArumugamHEpi, ArumugamZH, ArumugamZLogBF, 
   ArumugamBMI, Arumugam2, Arumugam)



########### Turnbaugh ####################################################


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/turnbaugh.twins")
shared.data <- read.table("test.unique.good.filter.unique.precluster.pick.pick.an.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
shared.data <- shared.data[, -c(1:3)]

subsample.data <- read.table("test.unique.good.filter.unique.precluster.pick.pick.an.0.03.subsample.shared", header=T)
rownames(subsample.data) <- subsample.data[, 2]
subsample.data <- subsample.data[, -c(1:3)]


metadata <- read.csv("turnbaugh.metadata.csv")
rownames(metadata) <- metadata[, 4]
metadata <- metadata[, -4]

keep1 <- which(metadata$Sample == 1) # n =146

#use the first sampling
#subset data for only the first sampling
s1.metadata <- metadata[keep1, ]
s1.subsample.data <- subsample.data[keep1, ]


#generate alpha diversity measures with vegan
H <- diversity(s1.subsample.data)
S <- specnumber(s1.subsample.data)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))
s1.alpha.diversity <- as.data.frame(select.alpha.diversity)
alpha.test <- s1.alpha.diversity

#Generate phyla table data

phyla.info <- read.csv("phyla.csv")
phyla.table <- shared.data
colnames(phyla.table) <- phyla.info[, 2]


#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(4:5, 7:9)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
s1.phyla.rel.abund <- phyla.table.rel.abund[rownames(shared.data), ]
s1.phyla.rel.abund <- s1.phyla.rel.abund[keep1, ]

#Create Obese Yes/No groups
s1.metadata$obese[s1.metadata$BMI.category=="Lean" | s1.metadata$BMI.category=="Overweight"] <- "No"
s1.metadata$obese[s1.metadata$BMI.category=="Obese"] <- "Yes"

#Create Obese.num group
s1.metadata$obese.num[s1.metadata$obese=="No"] <- 0
s1.metadata$obese.num[s1.metadata$obese=="Yes"] <- 1

#create groups to be used
obese <- factor(s1.metadata$obese)
bmi <- s1.metadata$BMI.category

####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

turnbaughH <- wilcox.test(H ~ obese) #P-value=0.1162
turnbaughS <- wilcox.test(S ~ obese) #P-value=0.05479
turnbaughJ <- wilcox.test(J ~ obese) #P-value=0.1748

#B and F tests against obesity
bacter <- s1.phyla.rel.abund$Bacteroidetes
firm <- s1.phyla.rel.abund$Firmicutes
BFratio <- bacter/firm

turnbaughBacter <- wilcox.test(bacter ~ obese) #P-value=0.3998
turnbaughFirm <- wilcox.test(firm ~ obese) #P-value=0.2245
turnbaughBF <- wilcox.test(BFratio ~ obese) #P-value=0.3354

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
turnbaugh2 <- adonis(s1.subsample.data ~ obese, permutations=1000)
turnbaughPERM <- turnbaugh2$aov.tab
#PERMANOVA=0.09491, pseudo-F=1.2114

###########################################################################
############ Relative Risk#################################################
###########################################################################

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
orderedHCat <- test3[, 1]
TurntotalN <- length(orderedHCat)
TurnHHighTotal <- length(orderedHCat[orderedHCat=="higher"])
TurnHighShannonGroup <- as.data.frame(table(test3[c(1:TurnHHighTotal), 2]))
TurnLowShannonGroup <- as.data.frame(table(test3[c((TurnHHighTotal + 1):TurntotalN), 2]))
group1 <- c(TurnHighShannonGroup[2, 2], TurnHighShannonGroup[1, 2])
group2 <- c(TurnLowShannonGroup[2, 2], TurnLowShannonGroup[1, 2])
#Group1 (Higher than median), obese = 47 and non-obese = 26
#Group2 (Lower than median), obese = 52 and non-obese = 21

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

turnbaughHEpi <- epi.2by2(r.test, method="cohort.count")
turnbaughHMassoc <- turnbaughHEpi$massoc
turnbaughHRR <- turnbaughHMassoc$RR.strata.score
turnbaughHRRsig <- turnbaughHMassoc$chisq.strata
## Risk Ratio = 1.11
## CI = 0.88, 1.38
## p-value = 0.376


##Run the RR for B/F ratio
Bacter = s1.phyla.rel.abund$Bacteroidetes
Firm = s1.phyla.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
orderedBFCat <- test4[, 1]
TurnBFHighTotal <- length(orderedBFCat[orderedBFCat=="higher"])
TurnHighBFGroup <- as.data.frame(table(test4[c(1:TurnBFHighTotal), 2]))
TurnLowBFGroup <- as.data.frame(table(test4[c((TurnBFHighTotal + 1):TurntotalN), 2]))
group1 <- c(TurnHighBFGroup[2, 2], TurnHighBFGroup[1, 2])
group2 <- c(TurnLowBFGroup[2, 2], TurnLowBFGroup[1, 2])
#Group1 (Higher than median), obese = 54 and non-obese = 19
#Group2 (Lower than median), obese = 45 and non-obese = 28

r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

turnbaughBFEpi <- epi.2by2(r.test, method="cohort.count")
turnbaughBFMassoc <- turnbaughBFEpi$massoc
turnbaughBFRR <- turnbaughBFMassoc$RR.strata.score
turnbaughBFRRsig <- turnbaughBFMassoc$chisq.strata
## Risk Ratio = 0.83
## CI = 0.66, 1.05
## p-value = 0.111

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

#Create Obese.num group
#s1.metadata$obese.num[s1.metadata$obese=="No"] <- 0
#s1.metadata$obese.num[s1.metadata$obese=="Yes"] <- 1
#obese <- factor(s1.metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
#testset <- Filter(function(x)(length(unique(x))>5), s1.subsample.data)
#testset <- cbind(obese, testset)
#colnames(testset)[1] <- "obese"
#testset <- cbind(testset, H, S, J, s1.phyla.rel.abund)

#Try AUCRF with default measures provided in readme
#set.seed(3)
#TurnbaughAUCFit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
# list of 11 Measures, AUCopt = 0.7764883

###########################################################################
############ Z-score Data Preparation ###################################
###########################################################################

TurnbaughZH <- scale(H)
TurnbaughZLogBF <- scale(log(BFratio))
TurnbaughObese <- as.character(obese)
TurnbaughBMICat <- as.character(bmi)


###########################################################################
############ Combining Data Together ######################################
###########################################################################

#Get data for demographics table
totalN <- c(totalN, length(rownames(s1.metadata)))
meanAge <- c(meanAge, "21-32")
SDAge <- c(SDAge, NA)
#temporary <- table(metadata$Sex)
males <- unname(c(males, NA))
females <- unname(c(females, NA))
temporary <- table(s1.metadata$Ancestry)
ancestry <- unname(c(ancestry, temporary[names(temporary) == "EA"] / sum(temporary)))
meanBMI <- c(meanBMI, NA)
SDBMI <- c(SDBMI, NA)
minBMI <- c(minBMI, NA)
maxBMI <- c(maxBMI, NA)

Turnbaugh <- c(turnbaughBacter$p.value, turnbaughFirm$p.value, 
               turnbaughBF$p.value, turnbaughH$p.value, turnbaughS$p.value, 
               turnbaughJ$p.value, turnbaughPERM[1,6])

overallPTable <- rbind(overallPTable, Turnbaugh)

tposH <- c(tposH, TurnLowShannonGroup[2, 2])
tnegH <- c(tnegH, TurnLowShannonGroup[1, 2])
cposH <- c(cposH, TurnHighShannonGroup[2, 2])
cnegH <- c(cnegH, TurnHighShannonGroup[1, 2])
RRH <- c(RRH, turnbaughHRR[1,1])
lowH <- c(lowH, turnbaughHRR[1,2])
highH <- c(highH, turnbaughHRR[1,3])

tposBF <- c(tposBF, TurnLowBFGroup[2, 2])
tnegBF <- c(tnegBF, TurnLowBFGroup[1, 2])
cposBF <- c(cposBF, TurnHighBFGroup[2, 2])
cnegBF <- c(cnegBF, TurnHighBFGroup[1, 2])
RRBF <- c(RRBF, turnbaughBFRR[1,1])
lowBF <- c(lowBF, turnbaughBFRR[1,2])
highBF <- c(highBF, turnbaughBFRR[1,3])

TurnbaughData <- as.data.frame(cbind(TurnbaughZH, TurnbaughZLogBF, TurnbaughObese, TurnbaughBMICat))
TurnbaughData$Study <- "Turnbaugh"
colnames(TurnbaughData) <- c("ZH", "ZLogBF", "Obese","BMICat", "Study")


rm(TurnLowShannonGroup, TurnHighShannonGroup, turnbaughHRR, turnbaughBacter, 
   turnbaughFirm, turnbaughBF, turnbaughH, turnbaughS, turnbaughJ, 
   turnbaughPERM, turnbaughHEpi, TurnbaughZH, TurnbaughZLogBF, TurnbaughObese, 
   TurnbaughBMICat, turnbaugh2, Turnbaugh)


####################################################################################
############################ Output Tables #########################################
####################################################################################

setwd("C:/users/marc/Desktop")

demographicsTable <- as.data.frame(cbind(totalN, meanAge, SDAge, males, females, 
                                         ancestry, meanBMI, SDBMI, minBMI, maxBMI))
demographicsTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Arumugam", "Turnbaugh")
write.csv(demographicsTable, "demographicsTable.csv")

write.csv(combinedData, "combinedData.csv")
write.csv(TurnbaughData, "TurnbaughData.csv")

rownames(overallPTable) <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                             "HMP", "Wu", "Arumugam", "Turnbaugh")
write.csv(overallPTable, "overallPTable.csv")

ShannonRRTable <- as.data.frame(cbind(tposH, tnegH, cposH, cnegH, RRH, lowH, highH))
ShannonRRTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Arumugam", "Turnbaugh")
write.csv(ShannonRRTable, "ShannonRRTable.csv")

BFRatioRRTable <- as.data.frame(cbind(tposBF, tnegBF, cposBF, cnegBF, RRBF, lowBF, highBF))
BFRatioRRTable$Study <- c("Baxter", "Ross", "Goodrich", "Escobar", "Zupancic", 
                          "HMP", "Wu", "Arumugam", "Turnbaugh")
write.csv(BFRatioRRTable, "BFRatioRRTable.csv")





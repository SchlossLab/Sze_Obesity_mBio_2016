## Strategy 1
## Obesesity and the bacterial microbiome
## Marc Sze
## March 9, 2016

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
First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

#Create a column with obese and extreme obese as single entity
demographics$BMIclass2[demographics$BMI.classification=="Normal"] <- "Normal"
demographics$BMIclass2[demographics$BMI.classification=="Overweight"] <- "Overweight"                     
demographics$BMIclass2[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Obese"


##Test BMI versus alpha diversity and phyla

bmi <- demographics$BMI

baxterH <- anova(lm(H ~ obese)) #P-value=0.01643
baxterS <- anova(lm(S ~ obese)) #P-value=0.01571
baxterJ <- anova(lm(J ~ obese)) #P-value=0.03241

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

baxterBacter <- anova(lm(bacter ~ obese)) #P-value=0.3175
baxterFirm <- anova(lm(firm ~ obese)) #P-value=0.6817
baxterBF <- anova(lm(BFratio ~ obese)) #P-value=0.6305

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


overallPTable <- as.data.frame(t(c(baxterBacter[1,5], baxterFirm[1,5], baxterBF[1,5], baxterH[1,5], baxterS[1,5], baxterJ[1,5], baxterPERM[1,6])))
colnames(overallPTable) <- c("Bacteroidetes", "Firmicutes", "BFRatio", "Shannon", "OTURich", "Evenness", "BrayC")

tpos <- BaxLowShannonGroup[2, 2]
tneg <- BaxLowShannonGroup[1, 2]

cpos <- BaxHighShannonGroup[2, 2]
cneg <- BaxHighShannonGroup[1, 2]

RR <- baxterHRR[1,1]

low <- baxterHRR[1,2]

high <- baxterHRR[1,3]


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

rossH <- anova(lm(H ~ obese)) #P-value=0.3154
rossS <- anova(lm(S ~ obese)) #P-value=0.8179
rossJ <- anova(lm(J ~ obese)) #P-value=0.2269

#B and F tests against obesity
bacter <- phyla.table.rel.abund$p__Bacteroidetes
firm <- phyla.table.rel.abund$p__Firmicutes
BFratio <- bacter/firm

rossBacter <- anova(lm(bacter ~ obese)) #P-value=0.06935
rossFirm <- anova(lm(firm ~ obese)) #P-value=0.09236
rossBF <- anova(lm(BFratio ~ obese)) #P-value=0.09534

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

Ross <- c(rossBacter[1,5], rossFirm[1,5], rossBF[1,5], rossH[1,5], rossS[1,5], rossJ[1,5], rossPERM[1,6])
overallPTable <- rbind(overallPTable, Ross)

tpos <- c(tpos, RossLowShannonGroup[2, 2])
tneg <- c(tneg, RossLowShannonGroup[1, 2])

cpos <- c(cpos, RossHighShannonGroup[2, 2])
cneg <- c(cneg, RossHighShannonGroup[1, 2])

RR <- c(RR, rossHRR[1,1])

low <- c(low, rossHRR[1,2])

high <- c(high, rossHRR[1,3])

RossData <- as.data.frame(cbind(RossZH, RossZLogBF, RossBMI))
RossData$Study <- "Ross"
colnames(RossData) <- c("ZH", "ZLogBF", "BMI", "Study")
combinedData <- rbind(combinedData, RossData)

rm(RossLowShannonGroup, RossHighShannonGroup, rossHRR, rossBacter, 
   rossFirm, rossBF, rossH, rossS, rossJ, rossPERM, rossHEpi, 
   RossZH, RossZLogBF, RossBMI, ross2)





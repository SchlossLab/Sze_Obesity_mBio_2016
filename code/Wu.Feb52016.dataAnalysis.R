#Wu et al.
# Single study analysis
#Feb 5th, 2016

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


setwd("C:/users/marc/Desktop/obesity2/COMBO/")

microbiome <- read.table("combined.good.subsample.shared", header=T)
rownames(microbiome) <- microbiome$Group
microbiome <- microbiome[, -c(1:3)]
seqData <- read.csv("COMBO_data_table.csv", header=T)
seqData <- seqData[-1, ]
rownames(seqData) <- seqData$Run_s
metadata <- read.table("BMI_and_Counts.txt", header=T)
metadata <- metadata[, c(1, 2)]

#get seqData set in line with microbiome data
namesToKeep <- rownames(microbiome)
test <- seqData[namesToKeep, ]
seqData <- test
rm(test)

#change all rownames to match sample IDs
rownames(seqData) <- seqData$submitted_subject_id_s
rownames(microbiome) <- rownames(seqData)
rownames(metadata) <- metadata[, 1]

#Match the metadata now with the microbiome data
namesToKeep <- rownames(microbiome)
test <- metadata[namesToKeep, ]
metadata <- test
rm(test)
rm(seqData, keep, namesToKeep)

#Get alpha diversity of the samples
library(vegan)

H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.table("combined.good.phyla.txt", header=T)
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
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(4, 6:7)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:6)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(keep, phyla.names, phyla.total, phylogenetic.info, phyla.table, phyla.good)

#Create BMI groups
metadata$BMI.class[metadata$BMI<=24] <- "Normal"
metadata$BMI.class[metadata$BMI>24 & metadata$BMI<30] <- "Overweight"
metadata$BMI.class[metadata$BMI>=30 & metadata$BMI<40] <- "Obese"
metadata$BMI.class[metadata$BMI>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
metadata$obese[metadata$BMI.class=="Normal" | metadata$BMI.class=="Overweight"] <- "No"
metadata$obese[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Yes"

#Create a column with obese and extreme obese as single entity
metadata$BMIclass2[metadata$BMI.class=="Normal"] <- "Normal"
metadata$BMIclass2[metadata$BMI.class=="Overweight"] <- "Overweight"                     
metadata$BMIclass2[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Obese"

#Get paitent demographics to be tested
bmi <- metadata$BMI
obese <- factor(metadata$obese)
bmi.class <- factor(metadata$BMI.class)
bmi.class2 <- factor(metadata$BMIclass2)


###########################################################################
######### First Level Analysis & Alpha Diversity with BMI #################
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.7493
anova(lm(H ~ bmi.class)) #P-value=0.2959
anova(lm(H ~ bmi.class2)) #P-value=0.2089

anova(lm(S ~ obese)) #P-value=0.7511
anova(lm(S ~ bmi.class)) #P-value=0.5128
anova(lm(S ~ bmi.class2)) #P-value=0.3891

anova(lm(J ~ obese)) #P-value=0.5557
anova(lm(J ~ bmi.class)) #P-value=0.2647
anova(lm(J ~ bmi.class2)) #P-value=0.1818

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.8276, R2=0.0008541
summary(lm(S ~ bmi)) #P-value=0.8195, R2=0.0009377
summary(lm(J ~ bmi)) #P-value=0.8521, R2=0.0006262

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.9173
anova(lm(bacter ~ bmi.class)) #P-value=0.5554
anova(lm(bacter ~ bmi.class2)) #P-value=0.6846

anova(lm(firm ~ obese)) #P-value=0.5765
anova(lm(firm ~ bmi.class)) #P-value=0.2157
anova(lm(firm ~ bmi.class2)) #P-value=0.2817

anova(lm(BFratio ~ obese)) #P-value=0.6735
anova(lm(BFratio ~ bmi.class)) #P-value=0.6208
anova(lm(BFratio ~ bmi.class2)) #P-value=0.4681

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(microbiome ~ obese, permutations=1000) 
#PERMANOVA=0.8921, pseudo-F=0.66065

set.seed(3)
adonis(microbiome ~ bmi.class, permutations=1000) 
#PERMANOVA=0.9321, pseudo-F=0.73921

set.seed(3)
adonis(microbiome ~ bmi.class2, permutations=1000) 
#PERMANOVA=0.75922, pseudo-F=0.8791

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 3.000655
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 98.5
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6432727
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:29), 4])
table(test3[c(30:58), 4])
#Group1 (Higher than median), obese = 3 and non-obese = 26
#Group2 (Lower than median), obese = 2 and non-obese = 27

group1 <- c(3, 26)
group2 <- c(2, 27)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.67
## CI = 0.12, 3.70
## p-value = 0.64

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 2.438335
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:29), 2]) #higher group
table(test4[c(30:58), 2]) #lower group
#Group1 (Higher than median), obese = 3 and non-obese = 26
#Group2 (Lower than median), obese = 2 and non-obese = 27

group1 <- c(3, 26)
group2 <- c(2, 27)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.67
## CI = 0.12, 3.70
## p-value = 0.64

###########################################################################
############ Classification using AUCRF####################################
###########################################################################

library(AUCRF)

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
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 15 Measures, AUCopt = 0.7075472
plot(fit)
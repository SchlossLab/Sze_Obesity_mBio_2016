#Zupancic et al.
# Single study analysis
#Feb 24th, 2016

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
shared.data <- read.table("combined.0.03.subsample.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
shared.data <- shared.data[, -c(1:3)]
BactFamily <- read.csv("BacterialFamilyData.csv")
rownames(BactFamily) <- BactFamily[, 1]
BactFamily <- BactFamily[, -1]
keep  <- rownames(shared.data)

test3 <- test[keep, ]
good.metadata <- test3
test4 <- test2[keep, ]
good.metadata2 <- test4

family <- BactFamily[keep, ]
microbiome <- shared.data[keep, ]
metadata <- cbind(good.metadata, good.metadata2)

rm(good.metadata, shared.data, keep, test, test2, test3, test4, BactFamily, good.metadata2, metadata2)

test <- metadata[!duplicated(metadata$submitted_subject_id_s), ]
keep  <- rownames(test)
metadata <- test
family <- family[keep, ]
microbiome <- microbiome[keep, ]

rm(test, keep)

metadata <- metadata[complete.cases(metadata), ]
keep <- rownames(metadata)
family <- family[keep, ]
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
phylogenetic.info <- read.csv("phyla.data.csv")
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
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Spirochaetes")], 1, sum)
phyla.table <- phyla.table[, -c(4, 6)]

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

anova(lm(H ~ obese)) #P-value=0.2671
anova(lm(S ~ obese)) #P-value=0.05423
anova(lm(J ~ obese)) #P-value=0.6229

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.8592
anova(lm(firm ~ obese)) #P-value=0.2564
anova(lm(BFratio ~ obese)) #P-value=0.989

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(microbiome ~ obese, permutations=1000) 
#PERMANOVA=0.6863, pseudo-F=0.92888

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

median(H) # 3.058573
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 92
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.677342
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:99), 4])
table(test3[c(100:199), 4])
#Group1 (Higher than median), obese = 28 and non-obese = 71
#Group2 (Lower than median), obese = 42 and non-obese = 58

group1 <- c(28, 71)
group2 <- c(42, 58)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.49
## CI = 1.01, 2.19
## p-value = 0.043


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.5414258
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:99), 2]) #higher group
table(test4[c(100:199), 2]) #lower group
#Group1 (Higher than median), obese = 36 and non-obese = 63
#Group2 (Lower than median), obese = 34 and non-obese = 66

group1 <- c(36, 63)
group2 <- c(34, 66)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.94
## CI = 0.64, 1.36
## p-value = 0.727

# Ruminococcaceae RR analysis
rumi <- as.data.frame(family$Ruminococcaceae)
colnames(rumi) <- "rumi"
median(rumi$rumi) # 9.2
rumi <- within(rumi, {rumiTest = ifelse(rumi <= median(rumi), "less", "higher")})

bmi.cat <- as.character(obese)
test3 <- cbind(rumi, bmi.cat)
test3 <- test3[order(test3$rumiTest), ]
table(test3[c(1:99), 3])
table(test3[c(100:199), 3])
#Group1 (Higher than median), obese = 32 and non-obese = 67
#Group2 (Lower than median), obese = 38 and non-obese = 62

group1 <- c(32, 67)
group2 <- c(38, 62)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.18
## CI = 0.80, 1.72
## p-value = 0.402


###########################################################################
############ Classification using AUCRF ###################################
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
summary(fit) # list of 40 Measures, AUCopt = 0.6365534
zupAUC <- fit$`OOB-AUCopt`
plot(fit)

#Try AUCRF with only phyla and alpha diversity measures in the testset
testset <- cbind(obese, H, S, J, phyla.table.rel.abund)
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) 
# list of 3 Measures, AUCopt = 0.600477
# S, Actinobacteria, H
plot(fit)

# Generate data for the ZScore analysis

H.corr <- H - mean(H)
sd(H) #0.5932651

BF.ratio.corr <- log(BFratio) - log(mean(BFratio))
ZBF <- scale(BFratio)
Zbacter <- scale(bacter)
Zfirm <- scale(firm)
sd(BF.ratio.corr) # 1.463952

HNormData <- cbind(H.corr, bmi)
BFNormData <- cbind(ZBF, Zbacter, Zfirm, BF.ratio.corr, bmi)

write.csv(HNormData, "Amish.H.norm.csv")
write.csv(BFNormData, "Amish.BF.norm.csv")

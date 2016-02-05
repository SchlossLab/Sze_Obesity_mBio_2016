#Arumugam et al.
# Single study analysis
#Feb 5th, 2016

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
library(vegan)

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
bmi.class <- factor(metadata$BMI.class)
bmi.class2 <- factor(metadata$BMIclass2)

####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.2174
anova(lm(H ~ bmi.class)) #P-value=0.5339
anova(lm(H ~ bmi.class2)) #P-value=0.3628

anova(lm(S ~ obese)) #P-value=0.9878
anova(lm(S ~ bmi.class)) #P-value=0.7083
anova(lm(S ~ bmi.class2)) #P-value=0.4985

anova(lm(J ~ obese)) #P-value=0.1874
anova(lm(J ~ bmi.class)) #P-value=0.3772
anova(lm(J ~ bmi.class2)) #P-value=0.2326

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.1022, R2=0.03185
summary(lm(S ~ bmi)) #P-value=0.7642, R2=0.001091
summary(lm(J ~ bmi)) #P-value=0.06682, R2=0.0399

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.1911
anova(lm(bacter ~ bmi.class)) #P-value=0.2975
anova(lm(bacter ~ bmi.class2)) #P-value=0.2263

anova(lm(firm ~ obese)) #P-value=0.2925
anova(lm(firm ~ bmi.class)) #P-value=0.3654
anova(lm(firm ~ bmi.class2)) #P-value=0.2998

anova(lm(BFratio ~ obese)) #P-value=0.05254
anova(lm(BFratio ~ bmi.class)) #P-value=0.1808
anova(lm(BFratio ~ bmi.class2)) #P-value=0.1284

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(microb.norm ~ obese, permutations=1000) 
#PERMANOVA=0.05495, pseudo-F=1.6047

set.seed(3)
adonis(microb.norm ~ bmi.class, permutations=1000) 
#PERMANOVA=0.1798, pseudo-F=1.176

set.seed(3)
adonis(microb.norm ~ bmi.class2, permutations=1000) 
#PERMANOVA=0.2138, pseudo-F=1.1686

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

median(H) # 2.886851
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 85
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.664365
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:42), 4])
table(test3[c(43:85), 4])
#Group1 (Higher than median), obese = 17 and non-obese = 25
#Group2 (Lower than median), obese = 20 and non-obese = 23

group1 <- c(17, 25)
group2 <- c(20, 23)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.15
## CI = 0.71, 1.87
## p-value = 0.575


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.9224021
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:42), 2]) #higher group
table(test4[c(43:85), 2]) #lower group
#Group1 (Higher than median), obese = 20 and non-obese = 22
#Group2 (Lower than median), obese = 17 and non-obese = 26

group1 <- c(20, 22)
group2 <- c(17, 26)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.83
## CI = 0.51, 1.35
## p-value = 0.452

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
testset <- Filter(function(x)(length(unique(x))>5), microb.norm)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 19 Measures, AUCopt = 0.7649212
plot(fit)
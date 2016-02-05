#Ross et al.
# Single study analysis
#Feb 5th, 2016

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/hispanic/")

hispanic.microb <- read.table("hispanic.subsample.shared", header=T)

metadata <- read.csv("s40168-015-0072-y-s1.csv")
sample.match <- read.csv("Hispanic_dataset.csv")

#Create a microbiome data table with sample names from metadata
test <- cbind(as.character(sample.match$Run_s), as.character(sample.match$Library_Name_s))
test <- test[order(test[, 1]), ]

test2 <- cbind(test, hispanic.microb)
rownames(test2) <- test2[, 2]
his.microb.edit <- test2[, -c(1:5)]
edit.metadata <- metadata[, -c(2:5)]
rownames(edit.metadata) <- edit.metadata[, 1]
edit.metadata <- edit.metadata[, -1]

#Create a metadata file in the order of the microbiome data
order1 <- rownames(his.microb.edit)
edit.metadata2 <- edit.metadata[order1, ]

#Get alpha diversity of the samples
library(vegan)

H <- diversity(his.microb.edit)
S <- specnumber(his.microb.edit)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
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
phyla.table$other <- apply(phyla.table[, c("Deinococcus-Thermus", "Elusimicrobia", "Fusobacteria", "Lentisphaerae", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(3:4, 6:7, 9:10)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- as.data.frame(phyla.table.rel.abund[, -8])


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
bmi.class <- factor(edit.metadata2$BMI.class)
bmi.class2 <- factor(edit.metadata2$BMI.class2)


######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla


anova(lm(H ~ obese)) #P-value=0.1657
anova(lm(H ~ bmi.class)) #P-value=0.5448
anova(lm(H ~ bmi.class2)) #P-value=0.3615

anova(lm(S ~ obese)) #P-value=0.1841
anova(lm(S ~ bmi.class)) #P-value=0.61
anova(lm(S ~ bmi.class2)) #P-value=0.4

anova(lm(J ~ obese)) #P-value=0.2538
anova(lm(J ~ bmi.class)) #P-value=0.6431
anova(lm(J ~ bmi.class2)) #P-value=0.492

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.343, R2=0.01476
summary(lm(S ~ bmi)) #P-value=0.345, R2=0.01466
summary(lm(J ~ bmi)) #P-value=0.472, R2=0.008502

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.2092
anova(lm(bacter ~ bmi.class)) #P-value=0.5436
anova(lm(bacter ~ bmi.class2)) #P-value=0.4374

anova(lm(firm ~ obese)) #P-value=0.2569
anova(lm(firm ~ bmi.class)) #P-value=0.1741
anova(lm(firm ~ bmi.class2)) #P-value=0.321

anova(lm(BFratio ~ obese)) #P-value=0.2508
anova(lm(BFratio ~ bmi.class)) #P-value=0.4966
anova(lm(BFratio ~ bmi.class2)) #P-value=0.506


###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(his.microb.edit ~ obese, permutations=1000) 
#PERMANOVA=0.7153, pseudo-F=0.73249

set.seed(3)
adonis(his.microb.edit ~ bmi.class, permutations=1000) 
#PERMANOVA=0.3357, pseudo-F=1.0772

set.seed(3)
adonis(his.microb.edit ~ bmi.class2, permutations=1000) 
#PERMANOVA=0.3017, pseudo-F=1.1061


###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 2.747097
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 56
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6819666
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:31), 4])
table(test3[c(32:63), 4])
#Group1 (Higher than median), obese = 18 and non-obese = 13
#Group2 (Lower than median), obese = 20 and non-obese = 12

group1 <- c(18, 13)
group2 <- c(20, 12)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.08
## CI = 0.72, 1.61
## p-value = 0.719


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 1.543046
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:31), 2]) #higher group
table(test4[c(32:63), 2]) #lower group
#Group1 (Higher than median), obese = 22 and non-obese = 9
#Group2 (Lower than median), obese = 16 and non-obese = 16

group1 <- c(22, 9)
group2 <- c(16, 16)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.70
## CI = 0.47, 1.07
## p-value = 0.089

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

library(AUCRF)

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
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 3 Measures, AUCopt = 0.7526
plot(fit)


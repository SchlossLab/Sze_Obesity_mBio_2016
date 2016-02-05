#Escobar et al.
# Single study analysis
#Feb 5th, 2016

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/columbian/")

columbian.microb <- read.table("columbian.0.03.subsample.shared", header=T)
metadata <- read.csv("columbian_dataset.csv")

#Organize the microbiome shared data
rownames(columbian.microb) <- columbian.microb[, 2]
columbian.microb <- columbian.microb[, -c(1:3)]

#Organize the metadata
rownames(metadata) <- metadata[, 5]

#Get only data that we are interested in
edit.metadata <- metadata[, -c(1:7, 9, 13:15, 16:52)]

#Sort metadata into the same order as the microbiome data
order1 <- rownames(columbian.microb)
edit.metadata2 <- edit.metadata[order1, ]


#Get alpha diversity of the samples
library(vegan)

H <- diversity(columbian.microb)
S <- specnumber(columbian.microb)
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
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(4:5, 7:9)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]


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
bmi.class <- factor(edit.metadata2$BMI.class)
bmi.class2 <- factor(edit.metadata2$BMI.class2)

######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.8636
anova(lm(H ~ bmi.class)) #P-value=0.9584
anova(lm(H ~ bmi.class2)) #P-value=0.9584

anova(lm(S ~ obese)) #P-value=0.1947
anova(lm(S ~ bmi.class)) #P-value=0.1432
anova(lm(S ~ bmi.class2)) #P-value=0.1432

anova(lm(J ~ obese)) #P-value=0.4762
anova(lm(J ~ bmi.class)) #P-value=0.764
anova(lm(J ~ bmi.class2)) #P-value=0.764

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.7895, R2=0.002587
summary(lm(S ~ bmi)) #P-value=0.4474, R2=0.02076
summary(lm(J ~ bmi)) #P-value=0.9794, R2=2.423x10-5

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.09253
anova(lm(bacter ~ bmi.class)) #P-value=0.1742
anova(lm(bacter ~ bmi.class2)) #P-value=0.1742

anova(lm(firm ~ obese)) #P-value=0.3602
anova(lm(firm ~ bmi.class)) #P-value=0.5409
anova(lm(firm ~ bmi.class2)) #P-value=0.5409

anova(lm(BFratio ~ obese)) #P-value=0.4663
anova(lm(BFratio ~ bmi.class)) #P-value=0.7545
anova(lm(BFratio ~ bmi.class2)) #P-value=0.7545

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################


set.seed(3)
adonis(columbian.microb ~ obese, permutations=1000) 
#PERMANOVA=0.08891, pseudo-F=1.3797

set.seed(3)
adonis(columbian.microb ~ bmi.class, permutations=1000) 
#PERMANOVA=0.2358, pseudo-F=1.1142

set.seed(3)
adonis(columbian.microb ~ bmi.class2, permutations=1000) 
#PERMANOVA=0.2358, pseudo-F=1.1142

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 4.080262
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 362
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6896595
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:15), 4])
table(test3[c(16:30), 4])
#Group1 (Higher than median), obese = 4 and non-obese = 11
#Group2 (Lower than median), obese = 6 and non-obese = 9

group1 <- c(4, 11)
group2 <- c(6, 9)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.50
## CI = 0.53, 4.26
## p-value = 0.439

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.2516729
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:15), 2]) #higher group
table(test4[c(16:30), 2]) #lower group
#Group1 (Higher than median), obese = 7 and non-obese = 8
#Group2 (Lower than median), obese = 3 and non-obese = 12

group1 <- c(7, 8)
group2 <- c(3, 12)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.43
## CI = 0.14, 1.35
## p-value = 0.121

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
testset <- Filter(function(x)(length(unique(x))>5), columbian.microb)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 3 Measures, AUCopt = 0.95
plot(fit)
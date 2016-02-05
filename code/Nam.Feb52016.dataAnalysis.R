#Nam et al.
# Single study analysis
#Feb 5th, 2016


# It should be noted that for the grouping part of this analysis it was done 
# in relationship to being overweight instead of obese. But for the combined 
# portion of the analysis or third arm (all data sets taken together) it was 
# used.  Conclusions on the other arms are not dependent on this specific 
# data set.

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/korean/")

microb <- read.table("DRR000776.0.03.subsample.shared", header=T)
rownames(microb) <- microb[, 2]
microb <- microb[,-c(1:3)]

metadata <- read.csv("korean.metadata.csv")

#Need to get rid of the repeated samples and the children
#First the metadata
metadata2 <- metadata[-c(2:3, 5:6, 8:9, 11:12, 14:15, 17:18, 19:24), ]
rownames(metadata2) <- metadata2[, 2]
metadata2 <- metadata2[, -c(1:2)]

#Second the microbiome data set
microb2 <- microb[-c(2:3, 5:6, 8:9, 11:12, 14:15, 17:18, 19:24), ]
rownames(microb2) <- rownames(metadata2)

#Get alpha diversity of the samples
library(vegan)

H <- diversity(microb2)
S <- specnumber(microb2)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)
# Seems to look okay....

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(microb2)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- microb2
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Lentisphaerae", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(4:5, 7:8)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]


#Create BMI groups
metadata2$BMI.class[metadata2$BMI<=24] <- "Normal"
metadata2$BMI.class[metadata2$BMI>24 & metadata2$BMI<30] <- "Overweight"
metadata2$BMI.class[metadata2$BMI>=30 & metadata2$BMI<40] <- "Obese"
metadata2$BMI.class[metadata2$BMI>=40] <- "Extreme Obesity"

#Create Overweight Yes/No groups
metadata2$overweight[metadata2$BMI.class=="Normal"] <- "No"
metadata2$overweight[metadata2$BMI.class=="Obese" | metadata2$BMI.class=="Extreme Obesity" | metadata2$BMI.class=="Overweight"] <- "Yes"

#Get paitent demographics data to be tested 
bmi <- metadata2$BMI
overweight <- factor(metadata2$overweight)


######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ overweight)) #P-value=0.961
anova(lm(S ~ overweight)) #P-value=0.5165
anova(lm(J ~ overweight)) #P-value=0.6547

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.8838, R2=0.001376
summary(lm(S ~ bmi)) #P-value=0.3153, R2=0.06293
summary(lm(J ~ bmi)) #P-value=0.7022, R2=0.009386

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ overweight)) #P-value=0.4382
anova(lm(firm ~ overweight)) #P-value=0.2609
anova(lm(BFratio ~ overweight)) #P-value=0.4089

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(microb2 ~ overweight, permutations=1000) 
#PERMANOVA=0.4236, pseudo-F=0.97425

######################################################################################## Relative Risk ###############################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 3.888484
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 156.5
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.7715805
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(overweight)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:9), 4])
table(test3[c(10:18), 4])
#Group1 (Higher than median), overweight = 4 and normal = 5
#Group2 (Lower than median), overweight = 2 and normal = 7

group1 <- c(4, 5)
group2 <- c(2, 7)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Overweight", "normal")
rownames(r.test) <- c("group1", "group2")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.50
## CI = 0.12, 2.08
## p-value = 0.317

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.4243408
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, overweight)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:9), 2]) #higher group
table(test4[c(10:18), 2]) #lower group
#Group1 (Higher than median), overweight = 2 and normal = 7
#Group2 (Lower than median), overweight = 4 and normal = 5

group1 <- c(2, 7)
group2 <- c(4, 5)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 2.00
## CI = 0.48, 8.31
## p-value = 0.317


###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

library(AUCRF)

#Create Obese.num group
metadata2$over.num[metadata2$overweight=="No"] <- 0
metadata2$over.num[metadata2$overweight=="Yes"] <- 1
overweight <- factor(metadata2$over.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microb2)
testset <- cbind(overweight, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=500, nodesize=1)
#Since this data set was so small had to modify the standardized        
#parameters.
summary(fit) # list of 2 Measures, AUCopt = 1.00
plot(fit)
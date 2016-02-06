#Goodrich et al.
# Single study analysis
#Feb 6th, 2016

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


setwd("C:/users/marc/Desktop/obesity2/twinsUK")

#Read in and match metadata to microbiome data
metadata <- read.csv("TwinsUKStudy_metadata.csv")
rownames(metadata) <- metadata[, 1]
shared.data <- read.table("combined.tx.1.subsample.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
microbiome <- shared.data[, -c(1:3)]
keep  <- rownames(metadata)

microbiome<- microbiome[keep, ]

rm(good.metadata, shared.data, keep)


#generate alpha diversity measures with vegan
library(vegan)
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
metadata$BMI.class[metadata$body_mass_index_s<=24] <- "Normal"
metadata$BMI.class[metadata$body_mass_index_s>24 & metadata$body_mass_index_s<30] <- "Overweight"
metadata$BMI.class[metadata$body_mass_index_s>=30 & metadata$body_mass_index_s<40] <- "Obese"
metadata$BMI.class[metadata$body_mass_index_s>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
metadata$obese[metadata$BMI.class=="Normal" | metadata$BMI.class=="Overweight"] <- "No"
metadata$obese[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Yes"

#Create a column with obese and extreme obese as single entity
metadata$BMIclass2[metadata$BMI.class=="Normal"] <- "Normal"
metadata$BMIclass2[metadata$BMI.class=="Overweight"] <- "Overweight"                     
metadata$BMIclass2[metadata$BMI.class=="Obese" | metadata$BMI.class=="Extreme Obesity"] <- "Obese"

#Get paitent demographics to be tested
bmi <- metadata$body_mass_index_s
obese <- factor(metadata$obese)
bmi.class <- factor(metadata$BMI.class)
bmi.class2 <- factor(metadata$BMIclass2)


###########################################################################
######### First Level Analysis & Alpha Diversity with BMI #################
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.6671
anova(lm(H ~ bmi.class)) #P-value=0.9732
anova(lm(H ~ bmi.class2)) #P-value=0.8987

anova(lm(S ~ obese)) #P-value=0.9073
anova(lm(S ~ bmi.class)) #P-value=0.9124
anova(lm(S ~ bmi.class2)) #P-value=0.8202

anova(lm(J ~ obese)) #P-value=0.6319
anova(lm(J ~ bmi.class)) #P-value=0.943
anova(lm(J ~ bmi.class2)) #P-value=0.8394

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.8662, R2=5.63x10-5
summary(lm(S ~ bmi)) #P-value=0.9142, R2=2.302x10-5
summary(lm(J ~ bmi)) #P-value=0.8582, R2=6.332x10-5

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.6661
anova(lm(bacter ~ bmi.class)) #P-value=0.3184
anova(lm(bacter ~ bmi.class2)) #P-value=0.3184

anova(lm(firm ~ obese)) #P-value=0.9688
anova(lm(firm ~ bmi.class)) #P-value=0.399
anova(lm(firm ~ bmi.class2)) #P-value=0.2455

anova(lm(BFratio ~ obese)) #P-value=0.8709
anova(lm(BFratio ~ bmi.class)) #P-value=0.7915
anova(lm(BFratio ~ bmi.class2)) #P-value=0.5976

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(microbiome ~ obese, permutations=1000) 
#PERMANOVA=0.004995, pseudo-F=2.9308

set.seed(3)
adonis(microbiome ~ bmi.class, permutations=1000) 
#PERMANOVA=0.005994, pseudo-F=1.8853

set.seed(3)
adonis(microbiome ~ bmi.class2, permutations=1000) 
#PERMANOVA=0.007992, pseudo-F=2.1453

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 2.725646
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 73
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6311328
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:253), 4])
table(test3[c(254:507), 4])
#Group1 (Higher than median), obese = 47 and non-obese = 206
#Group2 (Lower than median), obese = 56 and non-obese = 198

group1 <- c(47, 206)
group2 <- c(56, 198)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.19
## CI = 0.84, 1.68
## p-value = 0.332

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.5304444
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:253), 2]) #higher group
table(test4[c(254:507), 2]) #lower group
#Group1 (Higher than median), obese = 52 and non-obese = 201
#Group2 (Lower than median), obese = 51 and non-obese = 203

group1 <- c(52, 201)
group2 <- c(51, 203)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.98
## CI = 0.69, 1.38
## p-value = 0.894

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
summary(fit) # list of 8 Measures, AUCopt = 0.6767159
plot(fit)


#Try AUCRF with only phyla and alpha diversity measures in the testset
testset <- cbind(obese, H, S, J, phyla.table.rel.abund)
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) 
# list of 7 Measures, AUCopt = 0.538138
# Proteobacteria, unclassified, Verrucomicrobia
plot(fit)


## Generate Data for combined overall analysis
H.corr <- H - mean(H)
sd(H) #0.3489582

H.corr.data <- cbind(H.corr, metadata$body_mass_index_s)
write.csv(H.corr.data, "H.corr.data.csv")

ZBF <- scale(BFratio)
Zbacter <- scale(Bacter)
Zfirm <- scale(Firm)

BF.ratio.corr <- log(BFratio) - log(mean(BFratio))
sd(BFratio) #0.9802656

BF.corr.data <- cbind(ZBF, Zbacter, Zfirm, BF.ratio.corr, metadata$body_mass_index_s)
write.csv(BF.corr.data, "BF.corr.data.csv")











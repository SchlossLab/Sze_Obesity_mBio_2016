#Human Microbiome Project
# Single study analysis
#Feb 5th, 2016

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
library(vegan)
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
bmi.class <- factor(select.meta.cat$BMI_C)
bmi.class2 <- factor(select.meta.cat$BMI.class2)

######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.585
anova(lm(H ~ bmi.class)) #P-value=0.322
anova(lm(H ~ bmi.class2)) #P-value=0.322

anova(lm(S ~ obese)) #P-value=0.845
anova(lm(S ~ bmi.class)) #P-value=0.2051
anova(lm(S ~ bmi.class2)) #P-value=0.2051

anova(lm(J ~ obese)) #P-value=0.3559
anova(lm(J ~ bmi.class)) #P-value=0.3501
anova(lm(J ~ bmi.class2)) #P-value=0.3501

##linear correlations
summary(lm(H ~ bmi)) #P-value=0.7041, R2=0.0005689
summary(lm(S ~ bmi)) #P-value=0.2304, R2=0.005657
summary(lm(J ~ bmi)) #P-value=0.962, R2=8.932x10-6

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.3518
anova(lm(bacter ~ bmi.class)) #P-value=0.1343
anova(lm(bacter ~ bmi.class2)) #P-value=0.1343

anova(lm(firm ~ obese)) #P-value=0.636
anova(lm(firm ~ bmi.class)) #P-value=0.07813
anova(lm(firm ~ bmi.class2)) #P-value=0.07813

anova(lm(BFratio ~ obese)) #P-value=0.5688
anova(lm(BFratio ~ bmi.class)) #P-value=0.6807
anova(lm(BFratio ~ bmi.class2)) #P-value=0.6807

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(test ~ obese, permutations=1000) 
#PERMANOVA=0.8112, pseudo-F=0.7024

set.seed(3)
adonis(test ~ bmi.class, permutations=1000) 
#PERMANOVA=0.8022, pseudo-F=0.77896

set.seed(3)
adonis(test ~ bmi.class2, permutations=1000) 
#PERMANOVA=0.8022, pseudo-F=77896

######################################################################################## Relative Risk ###############################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 2.541824
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 57
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6305742
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:128), 4])
table(test3[c(129:256), 4])
#Group1 (Higher than median), obese = 10 and non-obese = 118
#Group2 (Lower than median), obese = 16 and non-obese = 112

group1 <- c(10, 118)
group2 <- c(16, 112)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.60
## CI = 0.75, 3.39
## p-value = 0.214

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 2.529168
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:128), 2]) #higher group
table(test4[c(129:256), 2]) #lower group
#Group1 (Higher than median), obese = 13 and non-obese = 115
#Group2 (Lower than median), obese = 13 and non-obese = 115

group1 <- c(13, 115)
group2 <- c(13, 115)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.00
## CI = 0.48, 2.07
## p-value = 1.00

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

library(AUCRF)

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
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 6 Measures, AUCopt = 0.7031773
plot(fit)
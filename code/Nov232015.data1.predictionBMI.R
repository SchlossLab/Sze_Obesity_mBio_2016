##Data set 1 prediction of BMI
## Nov 23, 2015


library(vegan)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(scales)
library(Boruta)
library(pheatmap)
library(randomForest)
library(AUCRF)


setwd("C:/Users/marc/Desktop/obesity2/data1")
demographics <- read.csv("demographics.v2.csv")
microbiome <- read.csv("data1.subsample.otus.csv")
phylogenetic.info <- read.csv("data1.summary.taxonomy.csv")
taxonomy <- read.csv("data1.taxonomy.csv")

rownames(demographics) <- demographics[,1]
rownames(microbiome) <- microbiome[,2]
demographics <- demographics[,-1]
microbiome <- microbiome[,-2]

test.samples <- rownames(demographics)
microbiome <- microbiome[,-c(1:2)]

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

#Generate Alpha Diversity Measures

H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)


#Create age quartiles

demographics$age.quartile <- with(demographics, cut(Age, breaks=quantile(Age, probs=seq(0,1, by=0.25)), include.lowest=TRUE))

#Create Obese Yes/No groups
demographics$obese[demographics$BMI.classification=="Normal" | demographics$BMI.classification=="Overweight"] <- "No"
demographics$obese[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Yes"

#Create Obese.num groups
demographics$obese.num[demographics$obese=="No"] <- 0
demographics$obese.num[demographics$obese=="Yes"] <- 1
obese <- factor(demographics$obese.num)

#generate test set
  # get rid of values that only have 0 and something else
testset <- Filter(function(x)(length(unique(x))>2), microbiome)
  # get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microbiome)

#Need to add phyla and alpha diversity measures to the dataset

testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(obese, testset)
         
#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 14 Measures, AUCopt = 0.7553
plot(fit)

#Try AUCRF with changing ntree 
#decrease first
set.seed(3)
fit2 <- AUCRF(obese ~ ., data=testset, ntree=800, nodesize=20)
summary(fit2) # list of 23 Measures, AUCopt = 0.7602
plot(fit2)

set.seed(3)
fit3 <- AUCRF(obese ~ ., data=testset, ntree=600, nodesize=20)
summary(fit3) # list of 14 Measures, AUCopt = 0.7489
plot(fit3)

#increase second
set.seed(3)
fit4 <- AUCRF(obese ~ ., data=testset, ntree=1200, nodesize=20)
summary(fit4) # list of 14 Measures, AUCopt = 0.7593
plot(fit4)
rm(testset)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), microbiome)
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(demographics$BMI, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset, ntree=1000, doTrace=2)
getSelectedAttributes(bmi.bortuaRegress)
  #OTU 00050, 00075, 00119, 00145, 00222, 00320, 00452
test <- randomForest(testset[,getSelectedAttributes(bmi.bortuaRegress)], testset$bmi)
  # No. variables at each split = 2, No. trees = 500, % Var.explained = 20.39, mean of squared residuals = 22.48

#See if using only caucasians improves any of the signal

caucasian <- which(demographics$White == 1)
cauc.microb <- microbiome[caucasian, ]
cauc.demo <- demographics[caucasian, ]
cauc.phyla <- phyla.table.rel.abund[caucasian, ]
alpha.diversity <- as.data.frame(alpha.diversity.shannon)
cauc.alpha <- alpha.diversity[caucasian, ]
rm(testset)

cauc.obese <- factor(cauc.demo$obese.num)

#generate caucasian only test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), cauc.microb)

#Need to add phyla and alpha diversity measures to the dataset

testset <- cbind(testset, cauc.alpha, cauc.phyla)
testset <- cbind(cauc.obese, testset)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(cauc.obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 3 Measures, AUCopt = 0.7847
plot(fit)

#Try AUCRF with changing ntree 
#decrease first
set.seed(3)
fit2 <- AUCRF(cauc.obese ~ ., data=testset, ntree=800, nodesize=20)
summary(fit2) # list of 13 Measures, AUCopt = 0.7900
plot(fit2)

set.seed(3)
fit3 <- AUCRF(cauc.obese ~ ., data=testset, ntree=600, nodesize=20)
summary(fit3) # list of 7 Measures, AUCopt = 0.7869
plot(fit3)
rm(testset)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), cauc.microb)
testset <- cbind(testset, cauc.alpha, cauc.phyla)
testset <- cbind(cauc.demo$BMI, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset, ntree=1000, doTrace=2)
getSelectedAttributes(bmi.bortuaRegress)
#OTU 00010, 00016, 00025, 00036, 00051, 00119, 00145, 00175, 00180, 00212, 
#cont. H, S, Unclassified
test <- randomForest(testset[,getSelectedAttributes(bmi.bortuaRegress)], testset$bmi)
# No. variables at each split = 4, No. trees = 500, % Var.explained = 29.08, mean of squared residuals = 19.41



##Korean Obesity (PLOS One 2011) data analysis
##Nov 23, 2015
##There are no obese Koreans in this data set
##So have to look at normal vs overweight
##Look at prediction

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

#Create Obese Yes/No groups
metadata2$obese[metadata2$BMI.class=="Normal" | metadata2$BMI.class=="Overweight"] <- "No"
metadata2$obese[metadata2$BMI.class=="Obese" | metadata2$BMI.class=="Extreme Obesity"] <- "Yes"

#Create Overweight Yes/No groups
metadata2$overweight[metadata2$BMI.class=="Normal"] <- "No"
metadata2$overweight[metadata2$BMI.class=="Obese" | metadata2$BMI.class=="Extreme Obesity" | metadata2$BMI.class=="Overweight"] <- "Yes"

#Create overweight.num group
metadata2$over.num[metadata2$overweight=="No"] <- 0
metadata2$over.num[metadata2$overweight=="Yes"] <- 1
overweight <- factor(metadata2$over.num)


#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microb2)
testset <- cbind(overweight, testset)
colnames(testset)[1] <- "overweight"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with small ntree, froze program when default was tried
set.seed(3)
fit <- AUCRF(overweight ~ ., data=testset)
summary(fit) # list of 2 Measures, AUCopt = 1.00
plot(fit)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), microb2)
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(metadata2$BMI, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset)
getSelectedAttributes(bmi.bortuaRegress)
#Variables Otu00057, Actinobacteria
test <- randomForest(testset[,getSelectedAttributes(bmi.bortuaRegress)], testset$bmi)
# No. variables at each split = 1, No. trees = 500, % Var.explained = 41.02, mean of squared residuals = 4.2288


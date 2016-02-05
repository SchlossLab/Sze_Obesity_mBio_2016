##Data Analysis
##HMP and Obesity
##Nov 23, 2015
##Looking at prediction



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
alpha.diversity <- as.data.frame(select.alpha.diversity)

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

#generate categorical variables of interest
bmi <- select.meta.cont$DTPBMI
bmi.cat <- select.meta.cat$BMI_C

#Create Obese Yes/No groups
select.meta.cat$obese[select.meta.cat$BMI_C== "normal" | select.meta.cat$BMI_C== "overweight"] <- "No"
select.meta.cat$obese[select.meta.cat$BMI_C== "obese"] <- "Yes"

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
summary(fit) # list of 6 Measures, AUCopt = 0.7032
plot(fit)

#Try AUCRF with changing ntree 
#decrease first
set.seed(3)
fit2 <- AUCRF(obese ~ ., data=testset, ntree=800, nodesize=20)
summary(fit2) # list of 6 Measures, AUCopt = 0.7054
plot(fit2)

set.seed(3)
fit3 <- AUCRF(obese ~ ., data=testset, ntree=600, nodesize=20)
summary(fit3) # list of 6 Measures, AUCopt = 0.7094
plot(fit3)

set.seed(3)
fit4 <- AUCRF(obese ~ ., data=testset, ntree=400, nodesize=20)
summary(fit4) # list of 8 Measures, AUCopt = 0.6747
plot(fit4)
rm(testset)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), test)
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(select.meta.cont$DTPBMI, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset, ntree=1000, doTrace=2)
getSelectedAttributes(bmi.bortuaRegress)
# 0 Vairables were important
#Try random forest without selection
randomForest(bmi ~ ., data=testset, ntree=1000, nodesize=20)
# get -7.59 percent of variation explained

##Try to see if using only Caucasians improves any of the readings.

caucasians <- which(meta.cat$WHITE_C == "Yes" & meta.cat$ASIAN_C == "No" & meta.cat$BLACK_C == "No")
cauc.alpha <- alpha.diversity[caucasians, ]
cauc.alpha <- cauc.alpha[-c(206:237),]
cauc.microb <- test[caucasians, ]
cauc.microb <- cauc.microb[-c(206:237),]
cauc.phyla <- phyla.table[caucasians, ]
cauc.phyla <- cauc.phyla[-c(206:237),]
cauc.meta.cat <- select.meta.cat[caucasians, ]
cauc.meta.cont <- select.meta.cont[caucasians, ]
cauc.meta.cat <- cauc.meta.cat[-c(206:237), ]
cauc.meta.cont <- cauc.meta.cont[-c(206:237), ]

cauc.obese <- factor(cauc.meta.cat$obese.num)
cauc.bmi <- cauc.meta.cont$DTPBMI

#generate caucasian only test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), cauc.microb)

#Need to add phyla and alpha diversity measures to the dataset

testset <- cbind(testset, cauc.alpha, cauc.phyla)
testset <- cbind(cauc.obese, testset)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(cauc.obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 8 Measures, AUCopt = 0.5984
plot(fit)

#Try AUCRF with changing ntree 
#decrease first
set.seed(3)
fit2 <- AUCRF(cauc.obese ~ ., data=testset, ntree=800, nodesize=20)
summary(fit2) # list of 6 Measures, AUCopt = 0.5849
plot(fit2)

set.seed(3)
fit3 <- AUCRF(cauc.obese ~ ., data=testset, ntree=600, nodesize=20)
summary(fit3) # list of 11 Measures, AUCopt = 0.5876
plot(fit3)
rm(testset)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), cauc.microb)
testset <- cbind(testset, cauc.alpha, cauc.phyla)
testset <- cbind(cauc.bmi, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset, ntree=1000, doTrace=2)
getSelectedAttributes(bmi.bortuaRegress)
#0 Attributes


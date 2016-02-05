##Hispanic Obesity data analysis
##Nov 11, 2015

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
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]

#Get paitent demographics to be tested

bmi <- edit.metadata2$BMI

#Create BMI groups
edit.metadata2$BMI.class[edit.metadata2$BMI<=24] <- "Normal"
edit.metadata2$BMI.class[edit.metadata2$BMI>24 & edit.metadata2$BMI<30] <- "Overweight"
edit.metadata2$BMI.class[edit.metadata2$BMI>=30 & edit.metadata2$BMI<40] <- "Obese"
edit.metadata2$BMI.class[edit.metadata2$BMI>=40] <- "Extreme Obesity"


#Create Obese Yes/No groups
edit.metadata2$obese[edit.metadata2$BMI.class=="Normal" | edit.metadata2$BMI.class=="Overweight"] <- "No"
edit.metadata2$obese[edit.metadata2$BMI.class=="Obese" | edit.metadata2$BMI.class=="Extreme Obesity"] <- "Yes"

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

#Try AUCRF with changing ntree 
#decrease first
set.seed(3)
fit2 <- AUCRF(obese ~ ., data=testset, ntree=800, nodesize=20)
summary(fit2) # list of 3 Measures, AUCopt = 0.7705
plot(fit2)

set.seed(3)
fit3 <- AUCRF(obese ~ ., data=testset, ntree=600, nodesize=20)
summary(fit3) # list of 3 Measures, AUCopt = 0.7574
plot(fit3)

set.seed(3)
fit4 <- AUCRF(obese ~ ., data=testset, ntree=400, nodesize=20)
summary(fit4) # list of 6 Measures, AUCopt = 0.7479
plot(fit4)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), his.microb.edit)
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)
testset <- cbind(edit.metadata2$BMI, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset, ntree=1000, doTrace=2)
getSelectedAttributes(bmi.bortuaRegress)
#OTU 0012, 0037
test <- randomForest(testset[,getSelectedAttributes(bmi.bortuaRegress)], testset$bmi)
# No. variables at each split = 1, No. trees = 500, % Var.explained = 12.25, mean of squared residuals = 23.95










##MetaHit (Danish) Obesity data analysis
##Nov 23, 2015
##Invesitgating prediction

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

# Generate Overall Relative Abundance for phyla - normalize to a value of 100
# Needs to be done since not all measures have 100 bacteria from metaphlan2
t.phyla <- as.data.frame(t(phyla))
t.phyla$other <- apply(t.phyla[, c("Acidobacteria", "Candidatus_Saccharibacteria", "Chlorobi", "Fusobacteria", "Spirochaetes")], 1, sum)
t.phyla <- t.phyla[, -c(1, 4:5, 7, 9)]
overall2 <- rowSums(t.phyla)
phyla.norm <- (t.phyla / overall) * 100

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
obese <- factor(metadata$obese.num)

#generate test set
  # get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), microb.norm)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.norm)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 15 Measures, AUCopt = 0.7635
plot(fit)

#Try AUCRF with changing ntree 
#decrease first
set.seed(3)
fit2 <- AUCRF(obese ~ ., data=testset, ntree=800, nodesize=20)
summary(fit2) # list of 24 Measures, AUCopt = 0.7652
plot(fit2)

set.seed(3)
fit3 <- AUCRF(obese ~ ., data=testset, ntree=600, nodesize=20)
summary(fit3) # list of 24 Measures, AUCopt = 0.7669
plot(fit3)

set.seed(3)
fit4 <- AUCRF(obese ~ ., data=testset, ntree=400, nodesize=20)
summary(fit4) # list of 15 Measures, AUCopt = 0.7903
plot(fit4)

set.seed(3)
fit5 <- AUCRF(obese ~ ., data=testset, ntree=200, nodesize=20)
summary(fit5) # list of 11 Measures, AUCopt = 0.7658
plot(fit5)

#increase second
set.seed(3)
fit6 <- AUCRF(obese ~ ., data=testset, ntree=1200, nodesize=20)
summary(fit6) # list of 15 Measures, AUCopt = 0.7764
plot(fit6)

set.seed(3)
fit7 <- AUCRF(obese ~ ., data=testset, ntree=1400, nodesize=20)
summary(fit7) # list of 31 Measures, AUCopt = 0.7720
plot(fit7)
rm(testset)

#Use random forest to look at regression
testset <- Filter(function(x)(length(unique(x))>5), microb.norm)
testset <- cbind(testset, H, S, J, phyla.norm)
testset <- cbind(metadata$BMI, testset)
colnames(testset)[1] <- "bmi"

set.seed(3)
bmi.bortuaRegress <- Boruta(bmi ~ ., data=testset, ntree=1000, doTrace=2)
getSelectedAttributes(bmi.bortuaRegress)
#OTU36, OTU88, OTU246, OTU349, Verrucomicrobia
test <- randomForest(testset[,getSelectedAttributes(bmi.bortuaRegress)], testset$bmi)
# No. variables at each split = 1, No. trees = 500, % Var.explained = 23.64, mean of squared residuals = 25.70




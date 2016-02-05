##Turnbaugh Twinset data
# Feb 3rd, 2016

#Put into desired format for StatisticsProcessing.Rmd

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

setwd("C:/users/marc/Desktop/obesity2/turnbaugh.twins")
shared.data <- read.table("test.unique.good.filter.unique.precluster.pick.pick.an.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
shared.data <- shared.data[, -c(1:3)]

subsample.data <- read.table("test.unique.good.filter.unique.precluster.pick.pick.an.0.03.subsample.shared", header=T)
rownames(subsample.data) <- subsample.data[, 2]
subsample.data <- subsample.data[, -c(1:3)]


metadata <- read.csv("turnbaugh.metadata.csv")
rownames(metadata) <- metadata[, 4]
metadata <- metadata[, -4]

keep1 <- which(metadata$Sample == 1) # n =146

#use the first sampling
#subset data for only the first sampling
s1.metadata <- metadata[keep1, ]
s1.subsample.data <- subsample.data[keep1, ]


#generate alpha diversity measures with vegan
library(vegan)
H <- diversity(s1.subsample.data)
S <- specnumber(s1.subsample.data)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))
s1.alpha.diversity <- as.data.frame(select.alpha.diversity)
alpha.test <- s1.alpha.diversity

#Generate phyla table data

phyla.info <- read.csv("phyla.csv")
phyla.table <- shared.data
colnames(phyla.table) <- phyla.info[, 2]


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
s1.phyla.rel.abund <- phyla.table.rel.abund[rownames(shared.data), ]
s1.phyla.rel.abund <- s1.phyla.rel.abund[keep1, ]

#Create Obese Yes/No groups
s1.metadata$obese[s1.metadata$BMI.category=="Lean" | s1.metadata$BMI.category=="Overweight"] <- "No"
s1.metadata$obese[s1.metadata$BMI.category=="Obese"] <- "Yes"

#Create Obese.num group
s1.metadata$obese.num[s1.metadata$obese=="No"] <- 0
s1.metadata$obese.num[s1.metadata$obese=="Yes"] <- 1

#create groups to be used
obese <- factor(s1.metadata$obese)
bmi <- s1.metadata$BMI.category

####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.06976
anova(lm(S ~ obese)) #P-value=0.02435
anova(lm(J ~ obese)) #P-value=0.1253


##linear correlations
summary(lm(H ~ bmi)) #P-value=0.1926, R2=0.02277
summary(lm(S ~ bmi)) #P-value=0.03434, R2=0.04606
summary(lm(J ~ bmi)) #P-value=0.2988, R2=0.01675

#B and F tests against obesity
bacter <- s1.phyla.rel.abund$Bacteroidetes
firm <- s1.phyla.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.5771
anova(lm(firm ~ obese)) #P-value=0.4073
anova(lm(BFratio ~ obese)) #P-value=0.2216

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(s1.subsample.data ~ obese, permutations=1000) 
#PERMANOVA=0.09491, pseudo-F=1.2114

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

median(H) # 4.327308
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 290.5
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.7611711
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:73), 4])
table(test3[c(74:146), 4])
#Group1 (Higher than median), obese = 47 and non-obese = 26
#Group2 (Lower than median), obese = 52 and non-obese = 21

group1 <- c(47, 26)
group2 <- c(52, 21)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.11
## CI = 0.88, 1.38
## p-value = 0.376


##Run the RR for B/F ratio
Bacter = s1.phyla.rel.abund$Bacteroidetes
Firm = s1.phyla.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.4764057
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:73), 2]) #higher group
table(test4[c(74:146), 2]) #lower group
#Group1 (Higher than median), obese = 54 and non-obese = 19
#Group2 (Lower than median), obese = 45 and non-obese = 28

group1 <- c(54, 19)
group2 <- c(45, 28)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.83
## CI = 0.66, 1.05
## p-value = 0.111

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

library(AUCRF)

#Create Obese.num group
s1.metadata$obese.num[s1.metadata$obese=="No"] <- 0
s1.metadata$obese.num[s1.metadata$obese=="Yes"] <- 1
obese <- factor(s1.metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), s1.subsample.data)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, s1.phyla.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 11 Measures, AUCopt = 0.7764883
plot(fit)





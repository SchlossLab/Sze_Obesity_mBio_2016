#Zupancic et al.
# Single study analysis
#Feb 5th, 2016

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


setwd("C:/users/marc/Desktop/obesity2/Amish/updatedGOOD")

#Read in and match metadata to microbiome data
metadata <- read.csv("amish_obesity_table2.csv")
test <- metadata[!duplicated(metadata$submitted_sample_id_s), ]
rownames(test) <- test[, 4]
shared.data <- read.table("combined.0.03.subsample.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
shared.data <- shared.data[, -c(1:3)]
keep  <- rownames(shared.data)

test2 <- test[keep, ]
good.metadata <- test2
rm(metadata, test, test2)

metadata <- good.metadata[!grepl("_s2", good.metadata$submitted_sample_id_s),]
keep <- rownames(metadata)
microbiome <- shared.data[keep, ]
rm(good.metadata, shared.data, keep)


#generate alpha diversity measures with vegan
library(vegan)
H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))
s1.alpha.diversity <- as.data.frame(select.alpha.diversity)
alpha.test <- s1.alpha.diversity

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
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
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Spirochaetes")], 1, sum)
phyla.table <- phyla.table[, -c(4, 6)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:6)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100

#Create Obese Yes/No groups
# One problem is that we are assuming those with "not provided" are defaulting to the not affected group
metadata$obese[metadata$subject_is_affected_s=="<not provided>" | metadata$subject_is_affected_s=="No"] <- "No"
metadata$obese[metadata$subject_is_affected_s=="Yes"] <- "Yes"


#Create Obese.num group
metadata$obese.num[metadata$obese=="No"] <- 0
metadata$obese.num[metadata$obese=="Yes"] <- 1

#create groups to be used
obese <- factor(metadata$obese)


####################################################################################### First Level Analysis & Alpha Diversity with BMI #############
###########################################################################
###########################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.2145
anova(lm(S ~ obese)) #P-value=0.05925
anova(lm(J ~ obese)) #P-value=0.5093

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.9485
anova(lm(firm ~ obese)) #P-value=0.2379
anova(lm(BFratio ~ obese)) #P-value=0.9277

####################################################################################### NMDS and PERMANOVA Analysis###################################
###########################################################################

set.seed(3)
adonis(microbiome ~ obese, permutations=1000) 
#PERMANOVA=0.5504, pseudo-F=0.97273

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity

median(H) # 3.073974
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 93
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6790902
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:100), 4])
table(test3[c(101:201), 4])
#Group1 (Higher than median), obese = 29 and non-obese = 71
#Group2 (Lower than median), obese = 45 and non-obese = 56

group1 <- c(29, 71)
group2 <- c(45, 56)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.54
## CI = 1.05, 2.24
## p-value = 0.022


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.537415
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:100), 2]) #higher group
table(test4[c(101:201), 2]) #lower group
#Group1 (Higher than median), obese = 38 and non-obese = 62
#Group2 (Lower than median), obese = 36 and non-obese = 65

group1 <- c(38, 62)
group2 <- c(36, 65)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.94
## CI = 0.65, 1.35
## p-value = 0.729

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

library(AUCRF)

#Create Obese.num group
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
summary(fit) # list of 32 Measures, AUCopt = 0.6679613
plot(fit)
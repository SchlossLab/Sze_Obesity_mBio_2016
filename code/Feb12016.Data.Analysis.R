##Amish Obesity data analysis
## Using a subsample of 500 instead (Go from 68 to 80 total samples)
##Feb 1, 2016

setwd("C:/users/marc/Desktop/obesity2/Amish")

#Read in and match metadata to microbiome data
metadata <- read.csv("amish_obesity_table2.csv")
test <- metadata[!duplicated(metadata$submitted_sample_id_s), ]
rownames(test) <- test[, 5]
shared.data <- read.table("Amish.subsample500.shared", header=T)
rownames(shared.data) <- shared.data[, 2]
shared.data <- shared.data[, -c(1:3)]
keep  <- rownames(shared.data)

test2 <- test[keep, ]
good.metadata <- test2
rm(metadata, test, test2)

#Create alpha diversity table
library(vegan)
H <- diversity(shared.data)
S <- specnumber(shared.data)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))
alpha.diversity <- as.data.frame(select.alpha.diversity)
alpha.test <- as.data.frame(alpha.diversity)
rm(select.alpha.diversity)


#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.info.csv")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(shared.data)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- shared.data
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Deferribacteres", "Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(3, 5:6, 8:10)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.table, phyla.info)


#Create Obese Yes/No groups
# One problem is that we are assuming those with "not provided" are defaulting to the not affected group
good.metadata$obese[good.metadata$subject_is_affected_s=="<not provided>" | good.metadata$subject_is_affected_s=="No"] <- "No"
good.metadata$obese[good.metadata$subject_is_affected_s=="Yes"] <- "Yes"

# First test of interest is BMI categories and alpha diversity measures
obese <- good.metadata$obese


######################################################################################## 
######### First Level Analysis & Alpha Diversity with BMI ##############################
########################################################################################

##Test BMI versus alpha diversity and phyla

anova(lm(H ~ obese)) #P-value=0.7507
anova(lm(S ~ obese)) #P-value=0.7528
anova(lm(J ~ obese)) #P-value=0.780


#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm
BFratio[67] <- NA # Made value NA since Firmicute value was 0

anova(lm(bacter ~ obese)) #P-value=0.5248
anova(lm(firm ~ obese)) #P-value=0.6641
anova(lm(BFratio ~ obese)) #P-value=0.7315

############################################################################
############ NMDS and PERMANOVA Analysis####################################
############################################################################

set.seed(3)
adonis(shared.data ~ obese, permutations=1000) 
#PERMANOVA=0.5734, pseudo-F=0.92921


############################################################################
############ Relative Risk##################################################
############################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe
#Shannon diversity
median(H) # 4.509466
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 210
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.8470437
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:40), 4])
table(test3[c(41:80), 4])
#Group1 (Higher than median), obese = 19 and non-obese = 21
#Group2 (Lower than median), obese = 16 and non-obese = 24

group1 <- c(19, 21)
group2 <- c(16, 24)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.84
## CI = 0.51, 1.39
## p-value = 0.499


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.3443504
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:40), 2]) #higher group
table(test4[c(41:80), 2]) #lower group
#Group1 (Higher than median), obese = 16 and non-obese = 24
#Group2 (Lower than median), obese = 19 and non-obese = 21

group1 <- c(16, 24)
group2 <- c(19, 21)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.19
## CI = 0.72, 1.96
## p-value = 0.499

############################################################################
############ Classification using AUCRF ####################################
############################################################################

library(AUCRF)

#Create Obese.num group
good.metadata$obese.num[good.metadata$obese=="No"] <- 0
good.metadata$obese.num[good.metadata$obese=="Yes"] <- 1
obese <- factor(good.metadata$obese.num)

#generate test set
# get rid of those with 0 and only 4 other values
testset <- Filter(function(x)(length(unique(x))>5), shared.data)
testset <- cbind(obese, testset)
colnames(testset)[1] <- "obese"
testset <- cbind(testset, H, S, J, phyla.table.rel.abund)

#Try AUCRF with default measures provided in readme
set.seed(3)
fit <- AUCRF(obese ~ ., data=testset, ntree=1000, nodesize=20)
summary(fit) # list of 4 Measures, AUCopt = 0.6247619
plot(fit)

## Generate Data for combined overall analysis
H.corr <- H - mean(H)
sd(H) #0.661348

H.corr.data <- cbind(H.corr, good.metadata$BMI)


ZBF <- scale(BFratio)
Zbacter <- scale(Bacter)
Zfirm <- scale(Firm)

test <- log(BFratio + 1)
#cannot have a value for log 0 so I used a +1 transformation
test2 <- BFratio + 1
meanTest <- mean(test2[-67])

BF.corr <- test - log(meanTest)


sd(test2[-67]) #9.865909

BF.corr.data <- cbind(ZBF, Zbacter, Zfirm, BF.corr, good.metadata$BMI)


dataset <- cbind(good.metadata, H.corr, BF.ratio.corr, Zbacter, Zfirm, ZBF)
write.csv(dataset, "Amish.normalized500.csv")






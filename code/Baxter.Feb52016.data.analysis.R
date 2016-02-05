#Baxter, et al
# Individual R file
# Feb 5, 2016


###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################

# Reading in the necessary Data
setwd("C:/Users/marc/Desktop/obesity2/data1")
demographics <- read.csv("demographics.v2.csv")
microbiome <- read.csv("data1.subsample.otus.csv")
phylogenetic.info <- read.csv("data1.summary.taxonomy.csv")
taxonomy <- read.csv("data1.taxonomy.csv")

# Minor modifications and Rownames adustments of data tables
rownames(demographics) <- demographics[,1]
rownames(microbiome) <- microbiome[,2]
demographics <- demographics[,-1]
microbiome <- microbiome[,-2]
test.samples <- rownames(demographics)
microbiome <- microbiome[,-c(1:2)]

# Create Table for phyla information
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

#Generate Alpha Diversity Measures and alpha diversity table
library(vegan)
H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)

#Create age quartiles

demographics$age.quartile <- with(demographics, cut(Age, breaks=quantile(Age, probs=seq(0,1, by=0.25)), include.lowest=TRUE))

#Create Obese Yes/No groups
demographics$obese[demographics$BMI.classification=="Normal" | demographics$BMI.classification=="Overweight"] <- "No"
demographics$obese[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Yes"

#Create Obese.num groups
demographics$obese.num[demographics$obese=="No"] <- 0
demographics$obese.num[demographics$obese=="Yes"] <- 1
obese <- factor(demographics$obese.num)

######################################################################################## First Level Analysis & Alpha Diversity with BMI #############
###########################################################################

#Create a column with obese and extreme obese as single entity
demographics$BMIclass2[demographics$BMI.classification=="Normal"] <- "Normal"
demographics$BMIclass2[demographics$BMI.classification=="Overweight"] <- "Overweight"                     
demographics$BMIclass2[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Obese"


##Test BMI versus alpha diversity and phyla

bmi <- demographics$BMI
bmi.class <- factor(demographics$BMI.classification)
bmi.class2 <- factor(demographics$BMIclass2)

anova(lm(H ~ obese)) #P-value=0.01643
anova(lm(H ~ bmi.class)) #P-value=0.06
anova(lm(H ~ bmi.class2)) #P-value=0.04471
a2 <- aov(H ~ bmi.class2)
TukeyHSD(x=a2, 'bmi.class2', conf.level=0.95)
# Obese-Norm (p-value=0.201), Over-Norm (p-value=0.769), Over-Obese (p-value=0.036)

anova(lm(S ~ obese)) #P-value=0.01571
anova(lm(S ~ bmi.class)) #P-value=0.06
anova(lm(S ~ bmi.class2)) #P-value=0.03
a2 <- aov(S ~ bmi.class2)
TukeyHSD(x=a2, 'bmi.class2', conf.level=0.95)
# Obese-Norm (p-value=0.025), Over-Norm (p-value=0.546), Over-Obese (p-value=0.178)

anova(lm(J ~ obese)) #P-value=0.03241
anova(lm(J ~ bmi.class)) #P-value=0.0272
a2 <- aov(J ~ bmi.class)
TukeyHSD(x=a2, 'bmi.class', conf.level=0.95)
# Obese-Norm (p-value=0.02), Over-Norm (p-value=0.047), Over-Obese (p-value=0.06)
anova(lm(J ~ bmi.class2)) #P-value=0.02
a2 <- aov(J ~ bmi.class2)
TukeyHSD(x=a2, 'bmi.class2', conf.level=0.95)
# Obese-Norm (p-value=0.630), Over-Norm (p-value=0.159), Over-Obese (p-value=0.019)

##linear correlations
H.bmi.test <- lm(H ~ bmi)
summary(H.bmi.test)  
#R2 = 0.01246, P-value=0.0774

S.bmi.test <- lm(S ~ bmi)
summary(S.bmi.test)
#R2 = 0.03222, P-value=0.0105

J.bmi.test <- lm(J ~ bmi)
summary(J.bmi.test)
#R2 = -0.0003052, P-value=0.332

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

anova(lm(bacter ~ obese)) #P-value=0.3175
anova(lm(bacter ~ bmi.class)) #P-value=0.6457
anova(lm(bacter ~ bmi.class2)) #P-value=0.5993

anova(lm(firm ~ obese)) #P-value=0.6817
anova(lm(firm ~ bmi.class)) #P-value=0.9293
anova(lm(firm ~ bmi.class2)) #P-value=0.8936

anova(lm(BFratio ~ obese)) #P-value=0.6305
anova(lm(BFratio ~ bmi.class)) #P-value=0.8334
anova(lm(BFratio ~ bmi.class2)) #P-value=0.881

###########################################################################
############ NMDS and PERMANOVA Analysis###################################
###########################################################################


table(demographics[,23])
###Normal=54, Overweight=71, Obese=45, Extreme Obesity=2
###Need to omit extreme obesity
#Removing rows with extreme obesity
test <- which(demographics$BMI.classification == "Extreme Obesity")
BMI.class.demo <- demographics[-test,]
BMI.microbiome <- microbiome[-test,]
bmi.class.group <- factor(BMI.class.demo$BMI.classification)

set.seed(3)
adonis(BMI.microbiome ~ bmi.class.group, permutations=1000)
#Not significant PERMANOVA = 0.2278, pseudo-F = 1.1323

set.seed(3)
adonis(microbiome ~ bmi.class, permutations=1000)
#Not significant PERMANOVA = 0.1848, pseudo-F = 1.1368

set.seed(3)
adonis(microbiome ~ bmi.class2, permutations=1000)
#Not significant PERMANOVA = 0.1359, pseudo-F = 1.1937

###########################################################################
############ Relative Risk#################################################
###########################################################################

library(epiR)

#Generate median values and put them into existing alpha.test dataframe

#Shannon diversity
median(H) # 3.67469
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 198
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.692766
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Run the RR for Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:86), 4])
table(test3[c(87:172), 4])
#Group1 (Higher than median), obese = 22 and non-obese = 64
#Group2 (Lower than median), obese = 25 and non-obese = 61

group1 <- c(22, 64)
group2 <- c(25, 61)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.14
## CI = 0.70, 1.85
## p-value = 0.608

##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.3621609
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:86), 2]) #higher group
table(test4[c(87:172), 2]) #lower group
#Group1 (Higher than median), obese = 23 and non-obese = 63
#Group2 (Lower than median), obese = 24 and non-obese = 62

group1 <- c(23, 63)
group2 <- c(24, 62)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.04
## CI = 0.64, 1.70
## p-value = 0.864

###########################################################################
############ Classification using AUCRF ###################################
###########################################################################

library(AUCRF)

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




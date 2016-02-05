##MetaHit (Danish) Obesity data analysis
##Nov 17, 2015

setwd("C:/users/marc/Desktop/obesity2/MetaHit/")

microb <- read.csv("finalized_metaHit_shared.csv")
rownames(microb) <- microb[, 1]
microb <- microb[,-1]

metadata <- read.csv("MetaHit_metadata.csv")
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]

# Generate Overall Relative Abundance - normalize to a value of 100
# Needs to be done since not all measures have 100 bacteria from metaphlan2
overall <- rowSums(microb)
microb.norm <- (microb / overall) * 100

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

#Plot NMDS and test with PERMANOVA
#Non-Obese versus Obese

obese <- factor(metadata$obese)
cols.obese <- c('blue', 'red')[obese]
obese.nmds <- capscale(microb.norm ~ 1, distance = "bray")
plot(obese.nmds, type="n")
points(obese.nmds, display="sites", pch=16, col=cols.obese)
legend("bottomright", c("Yes", "No"), pch=16, col=c('red', 'blue'))
set.seed(3)
adonis(microb.norm ~ obese, permutations=1000) 
#PERMANOVA=0.05495, pseudo-F=1.6047

## Test out if there is any correlations with BMI and alpha diversity measures

# First test of interest is BMI categories and alpha diversity measures

# Shannon Diversity
bmi <- metadata$BMI

anova(lm(H ~ obese)) #P-value=0.2174

summary(lm(H ~ bmi)) #P-value=0.1022, R2=0.02018
plot(H ~ bmi)
abline(lm(H ~ bmi), col="red")

# OTU Richness
anova(lm(S ~ obese)) #P-value=0.9878

summary(lm(S ~ bmi)) #P-value=0.7642, R2=-0.01094
plot(S ~ bmi)
abline(lm(S ~ bmi), col="red")

# Evenness
anova(lm(J ~ obese)) #P-value=0.1874

summary(lm(J ~ bmi)) #P-value=0.06682, R2=0.02833
plot(J ~ bmi)
abline(lm(J ~ bmi), col="red")


##Lets look at median to bring it inline with the other stuff

#Generate median values and put them into existing alpha.test dataframe

#Shannon diversity
median(H) # 2.886851
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 85
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.664365
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})


## Actually test the RR now

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:42), 4])
table(test3[c(43:85), 4])
#Group1 (Higher than median), obese = 17 and non-obese = 25
#Group2 (Lower than median), obese = 20 and non-obese = 23

group1 <- c(17, 25)
group2 <- c(20, 23)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.15
## CI = 0.71, 1.87
## p-value = 0.575


##Run the RR for B/F ratio
phyla.table <- read.csv("phyla.csv", header = T)
rownames(phyla.table) <- phyla.table[, 1]
phyla.table <- phyla.table[, -1]
phyla.table <- as.data.frame(t(phyla.table))

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:10)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.table, phyla.info)

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

ZscoreBacter <- scale(bacter)
ZScoreFirm <- scale(firm)
ZBF <- scale(BFratio)
metaHit.BandF.scales <- cbind(rownames(phyla.table.rel.abund), ZscoreBacter, ZScoreFirm, ZBF)
colnames(metaHit.BandF.scales) <- c("Sample", "Zbacter", "Zfirm", "ZBF")
write.csv(metaHit.BandF.scales, "metaHit.BF.phyla.stand.csv")


BFNorm = (log(BFratio) - log(1.143849))/log(0.9707577)

anova(lm(bacter ~ obese)) #P-value=0.1911
anova(lm(firm ~ obese)) #P-value=0.2925
anova(lm(BFratio ~ obese)) #P-value=0.05254


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

#Plot NMDS and test with PERMANOVA
#Non-Obese versus Obese

obese <- factor(edit.metadata2$obese)
cols.obese <- c('blue', 'red')[obese]
obese.nmds <- capscale(his.microb.edit ~ 1, distance = "bray")
plot(obese.nmds, type="n")
points(obese.nmds, display="sites", pch=16, col=cols.obese)
legend("bottomright", c("Yes", "No"), pch=16, col=c('red', 'blue'))
set.seed(3)
adonis(his.microb.edit ~ obese, permutations=1000) 
  #PERMANOVA=0.7153, pseudo-F=0.73249


## Test out if there is any correlations with BMI and alpha diversity measures

# First test of interest is BMI categories and alpha diversity measures

# Shannon Diversity
anova(lm(H ~ obese)) #P-value=0.1657

summary(lm(H ~ bmi)) #P-value=0.343, R2=-0.001388
plot(H ~ bmi)
abline(lm(H ~ bmi), col="red")

# OTU Richness
anova(lm(S ~ obese)) #P-value=0.1841

summary(lm(S ~ bmi)) #P-value=0.345, R2=-0.001497
plot(S ~ bmi)
abline(lm(S ~ bmi), col="red")

# Evenness
anova(lm(J ~ obese)) #P-value=0.2538

summary(lm(J ~ bmi)) #P-value=0.472, R2=-0.007752
plot(J ~ bmi)
abline(lm(J ~ bmi), col="red")

##Lets look at median to bring it inline with the other stuff

#Generate median values and put them into existing alpha.test dataframe

#Shannon diversity
median(H) # 2.747097
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 56
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6819666
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})


## Actually test the RR now

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:31), 4])
table(test3[c(32:63), 4])
#Group1 (Higher than median), obese = 21 and non-obese = 10
#Group2 (Lower than median), obese = 17 and non-obese = 15

group1 <- c(21, 10)
group2 <- c(17, 15)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.78
## CI = 0.52, 1.18
## p-value = 0.236





















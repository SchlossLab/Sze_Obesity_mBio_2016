##Korean Obesity (PLOS One 2011) data analysis
##Nov 17, 2015
##There are no obese Koreans in this data set
##So have to look at normal vs overweight

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

#Plot NMDS and test with PERMANOVA
#Non-Obese versus Obese

overweight <- factor(metadata2$overweight)
cols.overweight <- c('blue', 'red')[overweight]
overweight.nmds <- capscale(microb2 ~ 1, distance = "bray")
plot(overweight.nmds, type="n")
points(overweight.nmds, display="sites", pch=16, col=cols.overweight)
legend("bottomright", c("Yes", "No"), pch=16, col=c('red', 'blue'))
set.seed(3)
adonis(microb2 ~ overweight, permutations=1000) 
#PERMANOVA=0.4236, pseudo-F=0.97425

## Test out if there is any correlations with BMI and alpha diversity measures

# First test of interest is BMI categories and alpha diversity measures

# Shannon Diversity
bmi <- metadata2$BMI

anova(lm(H ~ overweight)) #P-value=0.961

summary(lm(H ~ bmi)) #P-value=0.8838, R2=-0.06104
plot(H ~ bmi)
abline(lm(H ~ bmi), col="red")

# OTU Richness
anova(lm(S ~ overweight)) #P-value=0.5165

summary(lm(S ~ bmi)) #P-value=0.3153, R2=0.004366
plot(S ~ bmi)
abline(lm(S ~ bmi), col="red")

# Evenness
anova(lm(J ~ overweight)) #P-value=0.6547

summary(lm(J ~ bmi)) #P-value=0.7022, R2=-0.05253
plot(J ~ bmi)
abline(lm(J ~ bmi), col="red")


##Lets look at median to bring it inline with the other stuff

#Generate median values and put them into existing alpha.test dataframe

#Shannon diversity
median(H) # 3.888484
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 156.5
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.7715805
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})


## Actually test the RR now

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
over.cat <- as.character(overweight)
test3 <- cbind(H.cat, S.cat, J.cat, over.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:9), 4])
table(test3[c(10:18), 4])
#Group1 (Higher than median), overweight = 4 and normal = 5
#Group2 (Lower than median), overweight = 2 and normal = 7

group1 <- c(4, 5)
group2 <- c(2, 7)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Overweight", "Normal")
rownames(r.test) <- c("group2", "group1")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.50
## CI = 0.12, 2.08
## p-value = 0.317
##Columbian Obesity data analysis
##Nov 10, 2015

setwd("C:/users/marc/Desktop/obesity2/columbian/")

columbian.microb <- read.table("columbian.0.03.subsample.shared", header=T)
metadata <- read.csv("columbian_dataset.csv")

#Organize the microbiome shared data
rownames(columbian.microb) <- columbian.microb[, 2]
columbian.microb <- columbian.microb[, -c(1:3)]

#Organize the metadata
rownames(metadata) <- metadata[, 5]

#Get only data that we are interested in
edit.metadata <- metadata[, -c(1:7, 9, 13:15, 16:52)]

#Sort metadata into the same order as the microbiome data
order1 <- rownames(columbian.microb)
edit.metadata2 <- edit.metadata[order1, ]


#Get alpha diversity of the samples
library(vegan)

H <- diversity(columbian.microb)
S <- specnumber(columbian.microb)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)
alpha.test <- as.data.frame(alpha.diversity.shannon)


#Get paitent demographics data to be tested 

bmi <- edit.metadata2$body_mass_index_s

#Create Obese Yes/No groups
edit.metadata2$obese[edit.metadata2$description_s=="Adequate weight male stool sample" | edit.metadata2$description_s=="Overweight male stool sample" | edit.metadata2$description_s=="Adequate weight female stool sample" | edit.metadata2$description_s=="Overweight female stool sample"] <- "No"
edit.metadata2$obese[edit.metadata2$description_s=="Obese male stool sample" | edit.metadata2$description_s=="Obese female stool sample"] <- "Yes"


#Plot NMDS and test with PERMANOVA
#Non-Obese versus Obese

obese <- factor(edit.metadata2$obese)
cols.obese <- c('blue', 'red')[obese]
obese.nmds <- capscale(columbian.microb ~ 1, distance = "bray")
plot(obese.nmds, type="n")
points(obese.nmds, display="sites", pch=16, col=cols.obese)
legend("bottomright", c("Yes", "No"), pch=16, col=c('red', 'blue'))
set.seed(3)
adonis(columbian.microb ~ obese, permutations=1000) 
#PERMANOVA=1.3797, pseudo-F=0.08891


## Test out if there is any correlations with BMI and alpha diversity measures

# First test of interest is BMI categories and alpha diversity measures

# Shannon Diversity
anova(lm(H ~ obese)) #P-value=0.8636

summary(lm(H ~ bmi)) #P-value=0.7895, R2=-0.03304
plot(H ~ bmi)
abline(lm(H ~ bmi), col="red")

# OTU Richness
anova(lm(S ~ obese)) #P-value=0.1947

summary(lm(S ~ bmi)) #P-value=0.4474, R2=-0.01421
plot(S ~ bmi)
abline(lm(S ~ bmi), col="red")

# Evenness
anova(lm(J ~ obese)) #P-value=0.4762

summary(lm(J ~ bmi)) #P-value=0.9794, R2=-0.03569
plot(J ~ bmi)
abline(lm(J ~ bmi), col="red")



##Lets look at median to bring it inline with the other stuff

#Generate median values and put them into existing alpha.test dataframe

#Shannon diversity
median(H) # 4.080262
alpha.test <- within(alpha.test, {shannon.cat = ifelse(H <= median(H), "less", "higher")})

#OTU Richness
median(S) # 362
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= median(S), "less", "higher")})

#Evenness
median(J) # 0.6896595
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})


## Actually test the RR now

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(obese)
test3 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test3 <- test3[order(H.cat), ]
table(test3[c(1:15), 4])
table(test3[c(16:30), 4])
#Group1 (Higher than median), obese = 4 and non-obese = 11
#Group2 (Lower than median), obese = 6 and non-obese = 9

group1 <- c(4, 11)
group2 <- c(6, 9)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.50
## CI = 0.53, 4.26
## p-value = 0.439





##Data Analysis
##HMP and Obesity
##Nov 3, 2015
##Investigating Relative Risk as a means of looking at obesity and microbiome



##Read in relevant data

setwd("C:/users/marc/Desktop/HMP.analysis")
microbiome <- read.table("Stool.an.0.03.subsample.shared", header=T)
meta.cat <- read.table("categorical.metadata", header=T)
meta.cont <- read.table("continuous.metadata", header=T)

#Only interested in first visit microbiome data
#need to create microbiome data set for only that

test <- microbiome[grep("\\.01\\.", microbiome$Group), ] #selecting by .01. 

microb.rownames <- gsub("([0-9]+).*", "\\1", test$Group) #extracting only numeric before first "."
rownames(test) <- microb.rownames

#Get rid of information not used in downstram analysis
test <- test[, -c(1:3)]

#Subset data to the microbiome total n
select.meta.cat <- meta.cat[microb.rownames, ]
select.meta.cont <- meta.cont[microb.rownames, ]

#generate alpha diversity measures with vegan
library(vegan)
H <- diversity(test)
S <- specnumber(test)
J <- H/log(S)
select.alpha.diversity <- as.data.frame(cbind(H, S, J))

#generate categorical variables of interest
bmi <- select.meta.cont$DTPBMI
bmi.cat <- select.meta.cat$BMI_C
caucasian <- select.meta.cat$WHITE_C
asian <-select.meta.cat$ASIAN_C
black <- select.meta.cat$BLACK_C

#Add bacterial microbiome categorical variables based on median values

#Shannon Diversity
median(H) #2.541824
alpha.test <- within(select.alpha.diversity, {shannon.cat = ifelse(H <= 2.541824, "less", "higher")})

#OTU Richness
median(S) #57
alpha.test <- within(alpha.test, {S.cat = ifelse(S <= 57, "less", "higher")})

#Evenness
median(J)
alpha.test <- within(alpha.test, {J.cat = ifelse(J <= median(J), "less", "higher")})

##Shannon Diversity
H.cat <- alpha.test$shannon.cat
S.cat <- alpha.test$S.cat
J.cat <- alpha.test$J.cat
bmi.cat <- as.character(select.meta.cat$BMI_C)
test2 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test2 <- test2[order(H.cat), ]
table(test2[c(1:128), 4])
table(test2[c(129:256), 4])
#Group1 (Higher than median), obese = 10 and non-obese = 118
#Group2 (Lower than median), obese = 16 and non-obese = 112

group1 <- c(10, 118)
group2 <- c(16, 112)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

library(epiR)
epi.2by2(r.test, method="cohort.count")
  ## Risk Ratio = 1.6
  ## CI = 0.75, 3.39
  ## p-value = 0.214


##OTU Richness

test2 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test2 <- test2[order(S.cat), ]
table(test2[c(1:128), 4])
table(test2[c(129:256), 4])
#Group1 (Higher than median), obese = 12 and non-obese = 116
#Group2 (Lower than median), obese = 14 and non-obese = 114

group1 <- c(12, 116)
group2 <- c(14, 114)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

epi.2by2(r.test, method="cohort.count")
  ## Risk Ratio = 1.17
  ## CI = 0.53, 2.68
  ## p-value = 0.679


##Evenness

test2 <- cbind(H.cat, S.cat, J.cat, bmi.cat)
test2 <- test2[order(J.cat), ]
table(test2[c(1:128), 4])
table(test2[c(129:256), 4])
#Group1 (Higher than median), obese = 12 and non-obese = 116
#Group2 (Lower than median), obese = 14 and non-obese = 114

group1 <- c(12, 116)
group2 <- c(14, 114)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group1", "group2")

epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.17
## CI = 0.53, 2.68
## p-value = 0.679

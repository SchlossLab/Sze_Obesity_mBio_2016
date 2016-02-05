##Data Analysis
##HMP and Obesity
##Nov 23, 2015
##Looking at prediction



##Read in relevant data

setwd("C:/users/marc/Desktop/obesity2/HMP.analysis")
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
alpha.diversity <- as.data.frame(select.alpha.diversity)

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(test)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- test
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Deinococcus-Thermus", "Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "Tenericutes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4, 6, 7, 9, 10:12)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]

#generate categorical variables of interest
bmi <- select.meta.cont$DTPBMI
bmi.cat <- select.meta.cat$BMI_C

#Create Obese Yes/No groups
select.meta.cat$obese[select.meta.cat$BMI_C== "normal" | select.meta.cat$BMI_C== "overweight"] <- "No"
select.meta.cat$obese[select.meta.cat$BMI_C== "obese"] <- "Yes"

#Create Obese.num group
select.meta.cat$obese.num[select.meta.cat$obese=="No"] <- 0
select.meta.cat$obese.num[select.meta.cat$obese=="Yes"] <- 1
obese <- factor(select.meta.cat$obese.num)

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

BFNorm <- (log(BFratio) - log(4.946849))/log(8.301864)

ZscoreBacter <- scale(bacter)
ZScoreFirm <- scale(firm)
ZBF <- scale(BFratio)
HMP.BandF.scales <- cbind(rownames(phyla.table.rel.abund), ZscoreBacter, ZScoreFirm, ZBF)
colnames(HMP.BandF.scales) <- c("Sample", "Zbacter", "Zfirm", "ZBF")
write.csv(HMP.BandF.scales, "HMP.BF.phyla.stand.csv")

anova(lm(bacter ~ obese)) #P-value=0.3518
anova(lm(firm ~ obese)) #P-value=0.636
anova(lm(BFratio ~ obese)) #P-value=0.5688


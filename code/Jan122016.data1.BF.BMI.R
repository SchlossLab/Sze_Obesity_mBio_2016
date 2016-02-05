##Data set 1 BF Relative Risk BMI
## Jan 12, 2016


setwd("C:/Users/marc/Desktop/obesity2/data1")
demographics <- read.csv("demographics.v2.csv")
microbiome <- read.csv("data1.subsample.otus.csv")
phylogenetic.info <- read.csv("data1.summary.taxonomy.csv")
taxonomy <- read.csv("data1.taxonomy.csv")

rownames(demographics) <- demographics[,1]
rownames(microbiome) <- microbiome[,2]
demographics <- demographics[,-1]
microbiome <- microbiome[,-2]

test.samples <- rownames(demographics)
microbiome <- microbiome[,-c(1:2)]

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

#Generate Alpha Diversity Measures

H <- diversity(microbiome)
S <- specnumber(microbiome)
J <- H/log(S)
alpha.diversity.shannon <- cbind(H,S,J)


#Create age quartiles

demographics$age.quartile <- with(demographics, cut(Age, breaks=quantile(Age, probs=seq(0,1, by=0.25)), include.lowest=TRUE))

#Create Obese Yes/No groups
demographics$obese[demographics$BMI.classification=="Normal" | demographics$BMI.classification=="Overweight"] <- "No"
demographics$obese[demographics$BMI.classification=="Obese" | demographics$BMI.classification=="Extreme Obesity"] <- "Yes"

#Create Obese.num groups
demographics$obese.num[demographics$obese=="No"] <- 0
demographics$obese.num[demographics$obese=="Yes"] <- 1
obese <- factor(demographics$obese.num)


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

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 1.04
## CI = 0.64, 1.70
## p-value = 0.864

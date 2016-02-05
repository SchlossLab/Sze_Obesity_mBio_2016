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

#B and F tests against obesity
bacter <- phyla.table.rel.abund$Bacteroidetes
firm <- phyla.table.rel.abund$Firmicutes
BFratio <- bacter/firm

BFNorm <- (log(BFratio) - log(0.4769099))/log(0.4451074)

ZscoreBacter <- scale(bacter)
ZScoreFirm <- scale(firm)
ZBF <- scale(BFratio)
data1.BandF.scales <- cbind(rownames(phyla.table.rel.abund), ZscoreBacter, ZScoreFirm, ZBF)
colnames(data1.BandF.scales) <- c("Sample", "Zbacter", "Zfirm", "ZBF")
write.csv(data1.BandF.scales, "data1.BF.phyla.stand.csv")

anova(lm(bacter ~ obese)) #P-value=0.3175
anova(lm(firm ~ obese)) #P-value=0.6817
anova(lm(BFratio ~ obese)) #P-value=0.6305

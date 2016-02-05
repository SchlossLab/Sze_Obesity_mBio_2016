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

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(columbian.microb)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- columbian.microb
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing)
#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(4:5, 7:9)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
phyla.table.rel.abund <- phyla.table.rel.abund[, -8]


#Get paitent demographics data to be tested 

bmi <- edit.metadata2$body_mass_index_s

#Create Obese Yes/No groups
edit.metadata2$obese[edit.metadata2$description_s=="Adequate weight male stool sample" | edit.metadata2$description_s=="Overweight male stool sample" | edit.metadata2$description_s=="Adequate weight female stool sample" | edit.metadata2$description_s=="Overweight female stool sample"] <- "No"
edit.metadata2$obese[edit.metadata2$description_s=="Obese male stool sample" | edit.metadata2$description_s=="Obese female stool sample"] <- "Yes"

#Create Obese.num group
edit.metadata2$obese.num[edit.metadata2$obese=="No"] <- 0
edit.metadata2$obese.num[edit.metadata2$obese=="Yes"] <- 1
obese <- factor(edit.metadata2$obese.num)


##Run the RR for B/F ratio
Bacter = phyla.table.rel.abund$Bacteroidetes
Firm = phyla.table.rel.abund$Firmicutes
BFRatio = Bacter/Firm
BFRatio <- as.data.frame(BFRatio)

median(BFRatio$BFRatio) # 0.2516729
BFRatio <- within(BFRatio, {BFRatio.cat = ifelse(BFRatio <= median(BFRatio), "less", "higher")})

BFRatio.cat <- BFRatio$BFRatio.cat
test4 <- cbind(BFRatio.cat, obese)
test4 <- test4[order(BFRatio.cat), ]
table(test4[c(1:15), 2]) #higher group
table(test4[c(16:30), 2]) #lower group
#Group1 (Higher than median), obese = 7 and non-obese = 8
#Group2 (Lower than median), obese = 3 and non-obese = 12

group1 <- c(7, 8)
group2 <- c(3, 12)
r.test <- rbind(group2, group1)
colnames(r.test) <- c("Obese", "Not.Obese")
rownames(r.test) <- c("group2", "group1")

library(epiR)
epi.2by2(r.test, method="cohort.count")
## Risk Ratio = 0.43
## CI = 0.14, 1.35
## p-value = 0.121

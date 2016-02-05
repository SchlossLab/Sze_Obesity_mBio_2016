# Testing the accuracy of metaphlan2 and 16S analysis
# linear correlation, bland-altman plot, and bar plots
# January 8th, 2016

#Phyla first then will explore other measures if data is promising

setwd("C:/Users/marc/Desktop/obesity2/HMP.analysis")

metaphlanData <- read.table("HMP_metaphlan2_merged_abundance_table.txt", header=T)
IDconversion <- read.csv("HMP.stool.IDs.csv")
microbiome <- read.table("Stool.an.0.03.subsample.shared", header=T)

#Only interested in first visit microbiome data
#need to create microbiome data set for only that

Data16S <- microbiome[grep("\\.01\\.", microbiome$Group), ] #selecting by .01. 

microb.rownames <- gsub("([0-9]+).*", "\\1", Data16S$Group) #extracting only numeric before first "."
rownames(Data16S) <- microb.rownames

#Get rid of information not used in downstram analysis
Data16S <- Data16S[, -c(1:3)]

#Get phyla information
#Edited out non phyla information first with sed in linux
#combined new labels with previous taxonomy file with excel
phylogenetic.info <- read.csv("phyla.data.csv")
rownames(phylogenetic.info) <- phylogenetic.info[,1]
phylogenetic.info <- phylogenetic.info[,-c(1)]
phyla.names <- as.character(phylogenetic.info$Taxonomy)
keep <- colnames(Data16S)
phyla.good <- phylogenetic.info[keep, ]
phyla.names <- as.character(phyla.good[,2])
phyla.table <- Data16S
colnames(phyla.table) <- phyla.names
#add all the same columns up and then return the sum
testing <- t(rowsum(t(phyla.table), group = rownames(t(phyla.table))))
phyla.table <- as.data.frame(testing)
rm(testing, keep, microb.rownames, phyla.names)

#combine phyla that are not that abundant
phyla.table$other <- apply(phyla.table[, c("Acidobacteria", "Deinococcus-Thermus", "Fusobacteria", "Lentisphaerae", "Spirochaetes", "Synergistetes", "Tenericutes", "TM7")], 1, sum)
phyla.table <- phyla.table[, -c(1, 4, 6, 7, 9, 10:12)]

#Create a relative abundance table for phyla
phyla.total <- apply(phyla.table[, c(1:7)], 1, sum)
phyla.table.rel.abund <- (phyla.table/phyla.total)*100
rm(phyla.total)

# Create comparable data set with 
metaphlanPhlyaData <- read.csv("HMP.metaphlan.phyla.data.csv", header=T)
rownames(metaphlanPhlyaData) <- metaphlanPhlyaData$ID
metaphlanPhlyaData <- metaphlanPhlyaData[,-1]
metaphlanPhlyaData <- as.data.frame(t(metaphlanPhlyaData))

#combine phyla that are not that abundant
metaphlanPhlyaData$other <- apply(metaphlanPhlyaData[, c("Acidobacteria", "Candidatus_Saccharibacteria", "Deinococcus_Thermus", "Fusobacteria", "Spirochaetes", "Synergistetes", "Tenericutes")], 1, sum)
metaphlanPhlyaData <- metaphlanPhlyaData[, -c(1, 4:5, 7, 9:11)]

#Create a relative abundance table for phyla
metaphlanPhlyaData.total <- apply(metaphlanPhlyaData[, c(1:6)], 1, sum)
metaphlanPhlyaData.rel.abund <- (metaphlanPhlyaData/metaphlanPhlyaData.total)*100

#Need to get sampleIDs to match 16S data
metaphlanNames <- rownames(metaphlanPhlyaData)
rownames(IDconversion) <- IDconversion[,1]
matchedData <- IDconversion[metaphlanNames, c(1,3)]
IndIDmatched <- matchedData[,2]

#Need to find and remove duplicated values (No idea which sample goes with which)
test1 <- as.data.frame(table(matchedData$Individual.ID))
good <- which(test1$Freq == 1)
unique <- test1[good, ]
uniqueNames <- as.character(unique$Var1)

test <- phyla.table.rel.abund[rownames(phyla.table.rel.abund) %in% uniqueNames, ]
finalUniqueNames <- rownames(test)
UpdatedMatchedData <- matchedData[matchedData$Individual.ID %in% finalUniqueNames, ]

goodMetaphlanPhylaData.rel.abund <- metaphlanPhlyaData.rel.abund[rownames(metaphlanPhlyaData.rel.abund) %in% rownames(UpdatedMatchedData), ]
rownames(goodMetaphlanPhylaData.rel.abund) <- UpdatedMatchedData[, 2]
testing <- goodMetaphlanPhylaData.rel.abund[order(rownames(goodMetaphlanPhylaData.rel.abund)), ]

phyla16SData <- test
phylaMetaphlanData <- testing
rm(finalUniqueNames, good, IndIDmatched, metaphlanNames, metaphlanPhlyaData.total, uniqueNames, goodMetaphlanPhylaData.rel.abund, matchedData, metaphlanPhlyaData, metaphlanPhlyaData.rel.abund, phyla.good, phyla.table, phyla.table.rel.abund, phylogenetic.info, test, test1, testing, unique, UpdatedMatchedData)

# Can start off with a general bar plot of the two side by side 
library(qdap)
library(plyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra) #using grid.arrange() can add different plots together
library(plotrix)

#phyla16SData <- phyla16SData[,-5] only if want to ignore unclassified
phylaMetaphlanData$unclassified <- 0
phylaMetaphlanData <- subset(phylaMetaphlanData, select=c(Actinobacteria:Proteobacteria, unclassified, Verrucomicrobia, other))
phyla16SData$dataset <- 1 
phylaMetaphlanData$dataset <- 2

testset <- rbind(phyla16SData, phylaMetaphlanData)
testset$dataset <- factor(testset$dataset)
melt.testset <- melt(testset, id.vars=c("dataset"))
colnames(melt.testset)[2] <- "phyla"
barplotAnalysisData <- ddply(melt.testset, c("dataset", "phyla"), summarise, mean = mean(value), sd=sd(value), sem=sd(value/sqrt(length(value))))
  
bpg1 <- ggplot(barplotAnalysisData, aes(x=dataset, y=mean, fill=phyla)) + geom_bar(position = "fill", stat="identity")
bpg2 <- bpg1 + scale_fill_brewer(palette="Set1") + scale_y_continuous(label=percent)


# Now can look at direct linear correlation between same sample
# Between metaphlan2 and 16S amplicon sequencing

#Actinobacteria
summary(lm(phyla16SData$Actinobacteria ~ phylaMetaphlanData$Actinobacteria)) 
#P-value 0.5543, R2= 0.006899
plot(phyla16SData$Actinobacteria ~ phylaMetaphlanData$Actinobacteria)
abline(lm(phyla16SData$Actinobacteria ~ phylaMetaphlanData$Actinobacteria), col="red")

#Bacteroidetes
summary(lm(phyla16SData$Bacteroidetes ~ phylaMetaphlanData$Bacteroidetes)) 
#P-value 0.0002277, R2= 0.2358
plot(phyla16SData$Bacteroidetes ~ phylaMetaphlanData$Bacteroidetes)
abline(lm(phyla16SData$Bacteroidetes ~ phylaMetaphlanData$Bacteroidetes), col="red")

#Firmicutes
summary(lm(phyla16SData$Firmicutes ~ phylaMetaphlanData$Firmicutes)) 
#P-value 0.1223, R2= 0.04617
plot(phyla16SData$Firmicutes ~ phylaMetaphlanData$Firmicutes)
abline(lm(phyla16SData$Firmicutes ~ phylaMetaphlanData$Firmicutes), col="red")

#Proteobacteria
summary(lm(phyla16SData$Proteobacteria ~ phylaMetaphlanData$Proteobacteria)) 
#P-value 0.004252, R2= 0.1494
plot(phyla16SData$Proteobacteria ~ phylaMetaphlanData$Proteobacteria)
abline(lm(phyla16SData$Proteobacteria ~ phylaMetaphlanData$Proteobacteria), col="red")

#Verrucomicrobia
summary(lm(phyla16SData$Verrucomicrobia ~ phylaMetaphlanData$Verrucomicrobia)) 
#P-value 0.1761, R2= 0.03559
plot(phyla16SData$Verrucomicrobia ~ phylaMetaphlanData$Verrucomicrobia)
abline(lm(phyla16SData$Verrucomicrobia ~ phylaMetaphlanData$Verrucomicrobia), col="red")

#other
summary(lm(phyla16SData$other ~ phylaMetaphlanData$other)) 
#P-value 0.8054, R2= 0.001201
plot(phyla16SData$other ~ phylaMetaphlanData$other)
abline(lm(phyla16SData$other ~ phylaMetaphlanData$other), col="red")

# Create a bland Altman of the best correlation
  # bias (difference of metaphlan) based on 16S analysis  
library(BlandAltmanLeh)
bland.altman.plot(phyla16SData$Bacteroidetes, phylaMetaphlanData$Bacteroidetes, main="Bland Altman Plot of Bacteroidetes", xlab="Means", ylab="Differences")




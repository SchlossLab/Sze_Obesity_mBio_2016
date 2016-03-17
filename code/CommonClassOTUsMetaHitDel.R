## Strategy 2: Investigating Common Classification OTUs
## Obesesity and the bacterial microbiome
## Marc Sze
## March 15, 2016


library(vegan)
library(epiR)
library(AUCRF)
source("code/UsedFunctions.R")

###########################################################################
############ Preparing Data Tables for Analysis ###########################
###########################################################################


####### Baxter #########

# Reading in the necessary Data
BaxterDemo <- read.csv("data/process/Baxter/demographics.v2.csv")
BaxterMicrobiome <- read.csv("data/process/Baxter/data1.subsample.otus.csv")

# Minor modifications and Rownames adustments of data tables
rownames(BaxterDemo) <- BaxterDemo[,1]
rownames(BaxterMicrobiome) <- BaxterMicrobiome[,2]
BaxterDemo <- BaxterDemo[,-1]
BaxterMicrobiome <- BaxterMicrobiome[,-2]
test.samples <- rownames(BaxterDemo)
microbiome <- BaxterMicrobiome[,-c(1:2)]
#Create Obese Yes/No groups
BaxterDemo$obese[BaxterDemo$BMI.classification=="Normal" | 
                   BaxterDemo$BMI.classification=="Overweight"] <- "No"
BaxterDemo$obese[BaxterDemo$BMI.classification=="Obese" | 
                   BaxterDemo$BMI.classification=="Extreme Obesity"] <- "Yes"

#Create Obese.num groups
BaxterDemo$obese.num[BaxterDemo$obese=="No"] <- 0
BaxterDemo$obese.num[BaxterDemo$obese=="Yes"] <- 1

#Create a column with obese and extreme obese as single entity
BaxterDemo$BMIclass2[BaxterDemo$BMI.classification=="Normal"] <- "Normal"
BaxterDemo$BMIclass2[BaxterDemo$BMI.classification=="Overweight"] <- "Overweight"                     
BaxterDemo$BMIclass2[BaxterDemo$BMI.classification=="Obese" | 
                       BaxterDemo$BMI.classification=="Extreme Obesity"] <- "Obese"

rm(microbiome, test.samples)

####### Ross #########

RossMicrobiome <- read.table("data/process/Ross/RossGoodSub.shared", header=T)
rownames(RossMicrobiome) <- RossMicrobiome[, 2]
RossMicrobiome <- RossMicrobiome[, -c(1:3)]

metadata <- read.csv("data/process/Ross/s40168-015-0072-y-s1.csv")
sample.match <- read.csv("data/process/Ross/Hispanic_dataset.csv")

#Create a microbiome data table with sample names from metadata
test <- cbind(as.character(sample.match$Run_s), as.character(sample.match$Library_Name_s))
test <- test[order(test[, 1]), ]
rownames(test) <- test[, 1]
keep <- rownames(RossMicrobiome)
test <- test[keep, ]
test2 <- cbind(test, RossMicrobiome)
rownames(test2) <- test2[, 2]
RossMicrobiome <- test2[, -c(1:2)]
edit.metadata <- metadata[, -c(2:5)]
rownames(edit.metadata) <- edit.metadata[, 1]
edit.metadata <- edit.metadata[, -1]

#Create a metadata file in the order of the microbiome data
order1 <- rownames(RossMicrobiome)
RossDemo <- edit.metadata[order1, ]

#Create BMI groups
RossDemo$BMI.class[RossDemo$BMI<=24] <- "Normal"
RossDemo$BMI.class[RossDemo$BMI>24 & RossDemo$BMI<30] <- "Overweight"
RossDemo$BMI.class[RossDemo$BMI>=30 & RossDemo$BMI<40] <- "Obese"
RossDemo$BMI.class[RossDemo$BMI>=40] <- "Extreme Obesity"


#Create Obese Yes/No groups
RossDemo$obese[RossDemo$BMI.class=="Normal" | RossDemo$BMI.class=="Overweight"] <- "No"
RossDemo$obese[RossDemo$BMI.class=="Obese" | RossDemo$BMI.class=="Extreme Obesity"] <- "Yes"

#Create BMI class 2 groups
RossDemo$BMI.class2[RossDemo$BMI<=24] <- "Normal"
RossDemo$BMI.class2[RossDemo$BMI>24 & RossDemo$BMI<30] <- "Overweight"
RossDemo$BMI.class2[RossDemo$BMI>=30] <- "Obese"

rm (test, test2, metadata, keep, order1)


####### Goodrich #########

#Read in and match metadata to microbiome data
GoodrichDemo <- read.csv("data/process/Goodrich/TwinsUKStudy2.csv")
rownames(GoodrichDemo) <- GoodrichDemo[, 7]
GoodrichMicrobiome <- read.table("data/process/Goodrich/GoodrichGoodSub.shared", header=T)
rownames(GoodrichMicrobiome) <- GoodrichMicrobiome[, 2]
GoodrichMicrobiome <- GoodrichMicrobiome[, -c(1:3)]
keep  <- which(GoodrichDemo$body_mass_index_s != "<not provided>")

GoodrichMicrobiome<- GoodrichMicrobiome[keep, ]
GoodrichDemo <- GoodrichDemo[keep, ]

#Create BMI groups
GoodrichDemo$BMI.class[as.numeric(as.character(GoodrichDemo$body_mass_index_s))<=24] <- "Normal"
GoodrichDemo$BMI.class[as.numeric(as.character(GoodrichDemo$body_mass_index_s))>24 & 
                     as.numeric(as.character(GoodrichDemo$body_mass_index_s))<30] <- "Overweight"
GoodrichDemo$BMI.class[as.numeric(as.character(GoodrichDemo$body_mass_index_s))>=30 & 
                     as.numeric(as.character(GoodrichDemo$body_mass_index_s))<40] <- "Obese"
GoodrichDemo$BMI.class[as.numeric(as.character(GoodrichDemo$body_mass_index_s))>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
GoodrichDemo$obese[GoodrichDemo$BMI.class=="Normal" | 
                     GoodrichDemo$BMI.class=="Overweight"] <- "No"
GoodrichDemo$obese[GoodrichDemo$BMI.class=="Obese" | 
                     GoodrichDemo$BMI.class=="Extreme Obesity"] <- "Yes"

#Create a column with obese and extreme obese as single entity
GoodrichDemo$BMIclass2[GoodrichDemo$BMI.class=="Normal"] <- "Normal"
GoodrichDemo$BMIclass2[GoodrichDemo$BMI.class=="Overweight"] <- "Overweight"                     
GoodrichDemo$BMIclass2[GoodrichDemo$BMI.class=="Obese" | 
                         GoodrichDemo$BMI.class=="Extreme Obesity"] <- "Obese"
colnames(GoodrichDemo)[12] <- "BMI"
GoodrichDemo$BMI <- as.numeric(as.character(GoodrichDemo$BMI))
rm(keep)


####### Escobar #########

EscobarMicrobiome <- read.table("data/process/Escobar/EscobarGoodSub.shared", header=T)
EscobarDemo <- read.csv("data/process/Escobar/columbian_dataset.csv")

#Organize the microbiome shared data
rownames(EscobarMicrobiome) <- EscobarMicrobiome[, 2]
EscobarMicrobiome <- EscobarMicrobiome[, -c(1:3)]

#Organize the metadata
rownames(EscobarDemo) <- EscobarDemo[, 5]

#Get only data that we are interested in
EscobarDemo <- EscobarDemo[, -c(1:7, 9, 13:15, 16:52)]

#Sort metadata into the same order as the microbiome data and get sex information
order1 <- rownames(EscobarMicrobiome)
EscobarDemo <- EscobarDemo[order1, ]
EscobarDemo$sex[EscobarDemo$description_s == "Adequate weight male stool sample" | 
                  EscobarDemo$description_s == "Obese male stool sample" | 
                  EscobarDemo$description_s == "Overweight male stool sample"] <- "M"
EscobarDemo$sex[EscobarDemo$description_s == "Adequate weight female stool sample" | 
                  EscobarDemo$description_s == "Obese female stool sample" | 
                  EscobarDemo$description_s == "Overweight female stool sample"] <- "F"

rm(order1)

#Create Obese Yes/No groups
EscobarDemo$obese[EscobarDemo$description_s=="Adequate weight male stool sample" | 
                    EscobarDemo$description_s=="Overweight male stool sample" | 
                    EscobarDemo$description_s=="Adequate weight female stool sample" | 
                    EscobarDemo$description_s=="Overweight female stool sample"] <- "No"
EscobarDemo$obese[EscobarDemo$description_s=="Obese male stool sample" | 
                    EscobarDemo$description_s=="Obese female stool sample"] <- "Yes"

#Create BMI groups
EscobarDemo$BMI.class[EscobarDemo$body_mass_index_s<=24] <- "Normal"
EscobarDemo$BMI.class[EscobarDemo$body_mass_index_s>24 & 
                        EscobarDemo$body_mass_index_s<30] <- "Overweight"
EscobarDemo$BMI.class[EscobarDemo$body_mass_index_s>=30 & 
                        EscobarDemo$body_mass_index_s<40] <- "Obese"
EscobarDemo$BMI.class[EscobarDemo$body_mass_index_s>=40] <- "Extreme Obesity"

#Create BMI groups2
EscobarDemo$BMI.class2[EscobarDemo$body_mass_index_s<=24] <- "Normal"
EscobarDemo$BMI.class2[EscobarDemo$body_mass_index_s>24 & EscobarDemo$body_mass_index_s<30] <- "Overweight"
EscobarDemo$BMI.class2[EscobarDemo$body_mass_index_s>=30] <- "Obese"

####### Zupancic #########

#Read in and match metadata to microbiome data
metadata <- read.csv("data/process/Zupancic/amish_obesity_table2.csv")
metadata2 <- read.csv("data/process/Zupancic/amish.metadata.csv")
test <- metadata[!duplicated(metadata$submitted_sample_id_s), ]
rownames(test) <- test[, 4]
test2 <- metadata2[!duplicated(metadata2$SUBJID), ]
rownames(test2) <- test2[, 1]
ZupancicMicrobiome <- read.table("data/process/Zupancic/ZupancicGoodSub.shared", header=T)
rownames(ZupancicMicrobiome) <- ZupancicMicrobiome[, 2]
ZupancicMicrobiome <- ZupancicMicrobiome[, -c(1:3)]
keep  <- rownames(ZupancicMicrobiome)

test3 <- test[keep, ]
good.metadata <- test3
test4 <- test2[keep, ]
good.metadata2 <- test4

ZupancicMicrobiome <- ZupancicMicrobiome[keep, ]
metadata <- cbind(good.metadata, good.metadata2)
test <- metadata[!duplicated(metadata$submitted_subject_id_s), ]
keep  <- rownames(test)
ZupancicDemo <- test
ZupancicMicrobiome <- ZupancicMicrobiome[keep, ]
ZupancicDemo <- ZupancicDemo[complete.cases(ZupancicDemo), ]
keep <- rownames(ZupancicDemo)
ZupancicMicrobiome <- ZupancicMicrobiome[keep, ]

#Create BMI groups
ZupancicDemo$BMI.class[ZupancicDemo$BMI<=24] <- "Normal"
ZupancicDemo$BMI.class[ZupancicDemo$BMI>24 & ZupancicDemo$BMI<30] <- "Overweight"
ZupancicDemo$BMI.class[ZupancicDemo$BMI>=30 & ZupancicDemo$BMI<40] <- "Obese"
ZupancicDemo$BMI.class[ZupancicDemo$BMI>=40] <- "Extreme Obesity"


#Create Obese Yes/No groups
ZupancicDemo$obese[ZupancicDemo$BMI.class=="Normal" | 
                     ZupancicDemo$BMI.class=="Overweight"] <- "No"
ZupancicDemo$obese[ZupancicDemo$BMI.class=="Obese" | 
                     ZupancicDemo$BMI.class=="Extreme Obesity"] <- "Yes"

#Create BMI class 2 groups
ZupancicDemo$BMI.class2[ZupancicDemo$BMI<=24] <- "Normal"
ZupancicDemo$BMI.class2[ZupancicDemo$BMI>24 & ZupancicDemo$BMI<30] <- "Overweight"
ZupancicDemo$BMI.class2[ZupancicDemo$BMI>=30] <- "Obese"

rm(good.metadata, good.metadata2, metadata, metadata2, test, test2, test3, test4, keep)


####### HMP #########

##Read in relevant data
HMPmicrobiome <- read.table("data/process/HMP/Stool.an.0.03.subsample.shared", header=T)
meta.cat <- read.table("data/process/HMP/categorical.metadata", header=T)
meta.cont <- read.table("data/process/HMP/continuous.metadata", header=T)

#Only interested in first visit microbiome data
#need to create microbiome data set for only that

test <- HMPmicrobiome[grep("\\.01\\.", HMPmicrobiome$Group), ] #selecting by .01.

microb.rownames <- gsub("([0-9]+).*", "\\1", test$Group) #extracting only numeric before first "."
rownames(test) <- microb.rownames

#Get rid of information not used in downstram analysis
HMPmicrobiome <- test[, -c(1:3)]

#Subset data to the microbiome total n
select.meta.cat <- meta.cat[microb.rownames, ]
select.meta.cont <- meta.cont[microb.rownames, ]
HMPDemo <- cbind(select.meta.cat, select.meta.cont)

#Create Obese Yes/No groups
HMPDemo$obese[HMPDemo$BMI_C=="normal" | HMPDemo$BMI_C=="overweight"] <- "No"
HMPDemo$obese[HMPDemo$BMI_C=="obese" | HMPDemo$BMI_C=="extreme obesity"] <- "Yes"

#Create BMI groups2
HMPDemo$BMI.class2[HMPDemo$DTPBMI<=24] <- "Normal"
HMPDemo$BMI.class2[HMPDemo$DTPBMI>24 & HMPDemo$DTPBMI<30] <- "Overweight"
HMPDemo$BMI.class2[HMPDemo$DTPBMI>=30] <- "Obese"

rm(meta.cat, meta.cont, select.meta.cat, select.meta.cont, test, microb.rownames)


####### Wu #########

WuMicrobiome <- read.table("data/process/Wu/WuGoodSub.shared", header=T)
rownames(WuMicrobiome) <- WuMicrobiome$Group
WuMicrobiome <- WuMicrobiome[, -c(1:3)]
WuDemo <- read.table("data/process/Wu/bmi_info.txt", header=T)

#Match the metadata now with the microbiome data
namesToKeep <- rownames(WuMicrobiome)
WuDemo <- WuDemo[namesToKeep, ]

#Create BMI groups
WuDemo$BMI.class[WuDemo$bmi<=24] <- "Normal"
WuDemo$BMI.class[WuDemo$bmi>24 & WuDemo$bmi<30] <- "Overweight"
WuDemo$BMI.class[WuDemo$bmi>=30 & WuDemo$bmi<40] <- "Obese"
WuDemo$BMI.class[WuDemo$bmi>=40] <- "Extreme Obesity"

#Create Obese Yes/No groups
WuDemo$obese[WuDemo$BMI.class=="Normal" | WuDemo$BMI.class=="Overweight"] <- "No"
WuDemo$obese[WuDemo$BMI.class=="Obese" | WuDemo$BMI.class=="Extreme Obesity"] <- "Yes"

#Create a column with obese and extreme obese as single entity
WuDemo$BMIclass2[WuDemo$BMI.class=="Normal"] <- "Normal"
WuDemo$BMIclass2[WuDemo$BMI.class=="Overweight"] <- "Overweight"                     
WuDemo$BMIclass2[WuDemo$BMI.class=="Obese" | WuDemo$BMI.class=="Extreme Obesity"] <- "Obese"

rm(namesToKeep)


####### Turnbaugh #########

TurnbaughMicrobiome <- read.table("data/process/Turnbaugh/TurnbaughSub.shared", header=T)
rownames(TurnbaughMicrobiome) <- TurnbaughMicrobiome[, 2]
TurnbaughMicrobiome <- TurnbaughMicrobiome[, -c(1:3)]


TurnbuaghDemo <- read.csv("data/process/Turnbaugh/turnbaugh.metadata.csv")
rownames(TurnbuaghDemo) <- TurnbuaghDemo[, 4]
TurnbuaghDemo <- TurnbuaghDemo[, -4]

keep1 <- which(TurnbuaghDemo$Sample == 1) # n =146

#use the first sampling
#subset data for only the first sampling
TurnbuaghDemo <- TurnbuaghDemo[keep1, ]
TurnbaughMicrobiome <- TurnbaughMicrobiome[keep1, ]

#Create Obese Yes/No groups
TurnbuaghDemo$obese[TurnbuaghDemo$BMI.category=="Lean" | 
                      TurnbuaghDemo$BMI.category=="Overweight"] <- "No"
TurnbuaghDemo$obese[TurnbuaghDemo$BMI.category=="Obese"] <- "Yes"

rm(keep1)


###########################################################################
########## Visualizing How the Different OTUs look ########################
################### Clostridiales OTUs ####################################
###########################################################################

library(gridExtra)

a1 <- ggBaseBox(BaxterMicrobiome, BaxterDemo, "Otu00479", "Baxter", log=TRUE)
b1 <- ggBaseBox(RossMicrobiome, RossDemo, "Otu00102", "Ross", log=TRUE)
c1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000287", "Goodrich", log=TRUE)
d1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000662", "Goodrich", log=TRUE)
e1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000145", "Goodrich", log=TRUE)
f1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000244", "Goodrich", log=TRUE)
g1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000442", "Goodrich", log=TRUE)
h1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000322", "Zupancic", log=TRUE)
i1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000062", "Zupancic", log=TRUE)
h1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000322", "Zupancic", log=TRUE)
j1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000256", "Zupancic", log=TRUE)
k1 <- ggBaseBox(HMPmicrobiome, HMPDemo, "Otu000293", "HMP", log=TRUE)
l1 <- ggBaseBox(WuMicrobiome, WuDemo, "Otu00191", "Wu", log=TRUE)

png("results/figures/suppFigClostridiales.png")
grid.arrange(a1, b1, c1, d1, e1, f1, g1, h1, i1, h1, j1, k1, l1)
dev.off()


###########################################################################
########## Visualizing How the Different OTUs look ########################
################### Lachnospiraceae OTUs ##################################
###########################################################################

library(gridExtra)

a1 <- ggBaseBox(BaxterMicrobiome, BaxterDemo, "Otu00036", "Baxter", log=TRUE)
b1 <- ggBaseBox(RossMicrobiome, RossDemo, "Otu00177", "Ross", log=TRUE)
c1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000180", "Goodrich", log=TRUE)
d1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000075", "Goodrich", log=TRUE)
e1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000122", "Goodrich", log=TRUE)
f1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000092", "Goodrich", log=TRUE)
g1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000673", "Goodrich", log=TRUE)
h1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000523", "Goodrich", log=TRUE)
i1 <- ggBaseBox(EscobarMicrobiome, EscobarDemo, 
                "Otu00067", "Escobar", log=TRUE)
h1 <- ggBaseBox(EscobarMicrobiome, EscobarDemo, 
                "Otu00074", "Escobar", log=TRUE)
j1 <- ggBaseBox(EscobarMicrobiome, EscobarDemo, 
                "Otu00080", "Escobar", log=TRUE)
k1 <- ggBaseBox(WuMicrobiome, WuDemo, "Otu00075", "Wu", log=TRUE)
l1 <- ggBaseBox(WuMicrobiome, WuDemo, "Otu00229", "Wu", log=TRUE)
m1 <- ggBaseBox(WuMicrobiome, WuDemo, "Otu00005", "Wu", log=TRUE)
n1 <- ggBaseBox(WuMicrobiome, WuDemo, "Otu00067", "Wu", log=TRUE)

png("results/figures/suppFigLachno.png")
grid.arrange(a1, b1, c1, d1, e1, f1, g1, h1, i1, h1, j1, k1, l1, m1, n1)
dev.off()



###########################################################################
########## Visualizing How the Different OTUs look ########################
################### Ruminocoaccaceae OTUs ##################################
###########################################################################

library(gridExtra)

a1 <- ggBaseBox(BaxterMicrobiome, BaxterDemo, "Otu00320", "Baxter", log=TRUE)
b1 <- ggBaseBox(RossMicrobiome, RossDemo, "Otu00022", "Ross", log=TRUE)
c1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000068", "Goodrich", log=TRUE)
d1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000255", "Goodrich", log=TRUE)
e1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu001145", "Goodrich", log=TRUE)
f1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000558", "Goodrich", log=TRUE)
g1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000297", "Goodrich", log=TRUE)
h1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000006", "Goodrich", log=TRUE)
i1 <- ggBaseBox(GoodrichMicrobiome, GoodrichDemo, 
                "Otu000132", "Goodrich", log=TRUE)
h1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000398", "Zupancic", log=TRUE)
j1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000850", "Zupancic", log=TRUE)
k1 <- ggBaseBox(ZupancicMicrobiome, ZupancicDemo, 
                "Otu000817", "Zupancic", log=TRUE)
l1 <- ggBaseBox(HMPmicrobiome, HMPDemo, "Otu000252", "HMP", log=TRUE)
m1 <- ggBaseBox(WuMicrobiome, WuDemo, "Otu00123", "Wu", log=TRUE)

png("results/figures/suppFigRumi.png")
grid.arrange(a1, b1, c1, d1, e1, f1, g1, h1, i1, h1, j1, k1, l1, m1)
dev.off()








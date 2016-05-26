# Table 1 Creation 
# Make all the measurements for the table here so it just needs to be read into
# the Rmd.


# Read in necessary tables
meta.bax <- read.table("../data/baxter/baxter.metadata", header = T)
meta.esc <- read.table("../data/escobar/escobar.metadata", header = T)
meta.goo <- read.table("../data/goodrich/goodrich.metadata", header = T)
meta.hmp <- read.table("../data/hmp/hmp.metadata", header = T)
meta.ros <- read.table("../data/ross/ross.metadata", header = T)
meta.sch <- read.table("../data/schubert/schubert.metadata", header = T)
meta.tur <- read.table("../data/turnbaugh/turnbaugh.metadata", header = T)
meta.wu <- read.table("../data/wu/wu.metadata", header = T)
meta.zee <- read.table("../data/zeevi/zeevi.metadata", header = T)
meta.zup <- read.table("../data/zupancic/zupancic.metadata", header = T)
beta_tests.summary <- read.table("../data/process/beta_tests.summary", 
                                 header=T, row.names = 1)

# Create data table entries
Study <- as.character(rownames(beta_tests.summary))

Age <- round(c(mean(meta.bax$age), mean(meta.esc$age), mean(meta.hmp$age), 
               mean(meta.ros$age), mean(meta.sch$age), mean(meta.tur$age), 
               mean(meta.wu$age), mean(as.numeric(as.character(meta.zup$age)), 
                                       na.rm=TRUE), mean(meta.goo$age), 
               mean(meta.zee$age, na.rm=TRUE)), 2)

AgeSD <- round(c(sd(meta.bax$age), sd(meta.esc$age), sd(meta.hmp$age), 
                 sd(meta.ros$age), sd(meta.sch$age), sd(meta.tur$age), 
                 sd(meta.wu$age), sd(as.numeric(as.character(meta.zup$age)), 
                                     na.rm=TRUE), sd(meta.goo$age), 
               sd(meta.zee$age, na.rm = TRUE)), 2)

AgeTable <- paste(as.character(Age)," (",as.character(AgeSD),")", sep="")

# No female/male information from Turnbaugh study so have to use NA 
females <- c(
  length(meta.bax$sex[meta.bax$sex == 'f']), 
  length(meta.esc$sex[meta.esc$sex == 'f']), 
  length(meta.hmp$sex[meta.hmp$sex == 'f']), 
  length(meta.ros$sex[meta.ros$sex == 'F']), 
  length(meta.sch$sex[meta.sch$sex == 'f']), NA, 
  length(meta.wu$sex[meta.wu$sex == 'f']), 
  length(meta.zup$sex[meta.zup$sex == 'f']), 
  length(meta.goo$sex[meta.goo$sex == 'f']), 
  length(meta.zee$sex[meta.zee$sex == 'f']))

# No female/male information from Turnbaugh study so have to use NA 
males <- c(
  length(meta.bax$sex[meta.bax$sex == 'm']), 
  length(meta.esc$sex[meta.esc$sex == 'm']), 
  length(meta.hmp$sex[meta.hmp$sex == 'm']), 
  length(meta.ros$sex[meta.ros$sex == 'M']), 
  length(meta.sch$sex[meta.sch$sex == 'm']), NA, 
  length(meta.wu$sex[meta.wu$sex == 'm']), 
  length(meta.zup$sex[meta.zup$sex == 'm']), 
  length(meta.goo$sex[meta.goo$sex == 'm']), 
  length(meta.zee$sex[meta.zee$sex == 'm']))

Sex <- paste(as.character(females), "|", as.character(males))

white <- c(length(meta.bax$white[meta.bax$white == TRUE]), 
           length(meta.esc$white[meta.esc$white == TRUE]), 
           length(meta.hmp$white[meta.hmp$white == TRUE]), 
           length(meta.ros$white[meta.ros$white == TRUE]), 
           length(meta.sch$white[meta.sch$white == TRUE]), 
           length(meta.tur$white[meta.tur$white == TRUE]), 
           length(meta.wu$white[meta.wu$white == TRUE]), 
           length(meta.zup$white[meta.zup$white == TRUE]), 
           length(meta.goo$white[meta.goo$white == TRUE]), NA)

nonwhite <- c(length(meta.bax$white[meta.bax$white == FALSE]), 
              length(meta.esc$white[meta.esc$white == FALSE]), 
              length(meta.hmp$white[meta.hmp$white == FALSE]), 
              length(meta.ros$white[meta.ros$white == FALSE]), 
              length(meta.sch$white[meta.sch$white == FALSE]), 
              length(meta.tur$white[meta.tur$white == FALSE]), 
              length(meta.wu$white[meta.wu$white == FALSE]), 
              length(meta.zup$white[meta.zup$white == FALSE]), 
              length(meta.goo$white[meta.zee$white == FALSE]), NA)

total <- c(length(rownames(meta.bax)), length(rownames(meta.esc)), 
           length(rownames(meta.hmp)), length(rownames(meta.ros)), 
           length(rownames(meta.sch)), length(rownames(meta.tur)), 
           length(rownames(meta.wu)), length(rownames(meta.zup)), 
           length(rownames(meta.zee)))

# white percentage of study population (white | nonwhite)
Ancestry <- paste(as.character(round((white/total)*100), 2), " | ", 
                  as.character(round((nonwhite/total)*100), 2), sep="")

ObeseTotal <- c(length(meta.bax$obese[meta.bax$obese == TRUE]), 
                length(meta.esc$obese[meta.esc$obese == TRUE]), 
                length(meta.hmp$obese[meta.hmp$obese == TRUE]), 
                length(meta.ros$obese[meta.ros$obese == TRUE]), 
                length(meta.sch$obese[meta.sch$obese == TRUE]), 
                length(meta.tur$obese[meta.tur$obese == TRUE]), 
                length(meta.wu$obese[meta.wu$obese == TRUE]), 
                length(meta.zup$obese[meta.zup$obese == TRUE]), 
                length(meta.goo$obese[meta.goo$obese == TRUE]), 
                length(meta.zee$obese[meta.zee$obese == TRUE]))

NonObeseTotal <- c(length(meta.bax$obese[meta.bax$obese == FALSE]), 
                   length(meta.esc$obese[meta.esc$obese == FALSE]), 
                   length(meta.hmp$obese[meta.hmp$obese == FALSE]), 
                   length(meta.ros$obese[meta.ros$obese == FALSE]), 
                   length(meta.sch$obese[meta.sch$obese == FALSE]), 
                   length(meta.tur$obese[meta.tur$obese == FALSE]), 
                   length(meta.wu$obese[meta.wu$obese == FALSE]), 
                   length(meta.zup$obese[meta.zup$obese == FALSE]), 
                   length(meta.goo$obese[meta.goo$obese == FALSE]), 
                   length(meta.zee$obese[meta.zee$obese == FALSE]))

ObBreakdown <- paste(as.character(round((ObeseTotal/total)*100), 2), " | ", 
                     as.character(round((NonObeseTotal/total)*100), 2))

# For BMI stuff Turnbaugh did not have data available so have to use NA
meanBMI <- as.character(round(c(mean(meta.bax$bmi), mean(meta.esc$bmi), 
                                mean(meta.hmp$bmi), mean(meta.ros$bmi), 
                                mean(meta.sch$bmi), NA, mean(meta.wu$bmi), 
                                mean(as.numeric(as.character(meta.zup$bmi)), 
                                     na.rm=TRUE), mean(meta.goo$bmi), 
                                mean(meta.zee$bmi, na.rm=TRUE)), 2))

Min <- as.character(round(c(min(meta.bax$bmi), min(meta.esc$bmi), 
                            min(meta.hmp$bmi), min(meta.ros$bmi), 
                            min(meta.sch$bmi), NA, min(meta.wu$bmi), 
                            min(as.numeric(as.character(meta.zup$bmi)), 
                                na.rm=TRUE), min(meta.goo$bmi), 
                            min(meta.zee$bmi, na.rm=TRUE)), 2))

Max <- as.character(round(c(max(meta.bax$bmi), max(meta.esc$bmi), 
                            max(meta.hmp$bmi), max(meta.ros$bmi), 
                            max(meta.sch$bmi), NA, max(meta.wu$bmi), 
                            max(as.numeric(as.character(meta.zup$bmi)), 
                                na.rm=TRUE), max(meta.goo$bmi), 
                            max(meta.zee$bmi, na.rm=TRUE)), 2))

table1 <- cbind(Study, AgeTable, Sex, Ancestry, ObBreakdown, meanBMI, Min, Max)
table1[table1=="NA | NA" | table1=="NA (NA)"] <- "Not Available"
table1[which(table1[, 'Study'] == "turnbaugh"), 'meanBMI'] <- "Not Available"
table1[which(table1[, 'Study'] == "turnbaugh"), 'Min'] <- "Not Available"
table1[which(table1[, 'Study'] == "turnbaugh"), 'Max'] <- "Not Available"

write.table(table1, file="../results/tables/table1.demographics", sep='\t', 
            row.names=F)
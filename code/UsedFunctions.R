## Functions For Meta Analysis Paper on
## Obesesity and the bacterial microbiome
## Marc Sze
## March 16, 2016

## Function to make table for RR Analysis and run RR Test.
## Needed Libraries include epiR.  
## An Example of how to use:
  ## RunRR(alpha.test, demographics, "obese", "H")
RunRR <- function(testData, demo, testVar1, testVar2){
  
  # Check for necessary Libraries and load in not already attached
  deps = c("epiR");
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
  
  # Generate categorical Data
  Var1 <- demo[testVar1]
  
  # Generate Median Categorical Data
  VarOfInt <- testData[testVar2]
  VarOfInt <- as.data.frame(within(VarOfInt, 
                     {Var.cat = ifelse(
                       VarOfInt[,1] <= median(VarOfInt[,1]), 
                       "less", "higher")}))
  VarOfInt <- cbind(VarOfInt, Var1)
  VarOfInt <- VarOfInt[order(VarOfInt[, 2]), ]
  
  # Generate Test table
  totalN <- length(VarOfInt[, 2])
  HighTotal <- length(which(VarOfInt[, 2]=="higher"))
  HighMedianGroup <- as.data.frame(table(
    VarOfInt[c(1:HighTotal), 3]))
  LowMedianGroup <- as.data.frame(table(
    VarOfInt[c((HighTotal + 1):totalN), 3]))
  ConGroup <- c(HighMedianGroup[2, 2], HighMedianGroup[1, 2])
  TestGroup <- c(LowMedianGroup[2, 2], LowMedianGroup[1, 2])
  r.test <- rbind(TestGroup, ConGroup)
  colnames(r.test) <- c("Yes", "No")
  rownames(r.test) <- c("TestGroup", "ConGroup")
  
  # Run test and store important variables
  EpiTest <- epi.2by2(r.test, method="cohort.count")
  MassocTest <- EpiTest$massoc
  RRValues <- MassocTest$RR.strata.score
  RRsignif <- MassocTest$chisq.strata
  
  # Create Other variables to store for future use if needed
  tpos <- LowMedianGroup[2, 2]
  tneg <- LowMedianGroup[1, 2]
  cpos <- HighMedianGroup[2, 2]
  cneg <- HighMedianGroup[1, 2]
  RR <- RRValues[1,1]
  lowCI <- RRValues[1,2]
  highCI <- RRValues[1,3]
  
  output <- list(totalN, tpos, tneg, cpos, cneg, RR, lowCI, highCI, RRsignif[1, 3])
  names(output) <- c("totalN", "tpos", "tneg", "cpos", "cneg", 
                     "RR", "lowCI", "highCI", "pValue")
  
  return(output)
  
}

## Function to get the Simulated Power of Study.
## Need to make sure that the TestData and demo tables rownames match
## Example on how to use:
## NonParaPowerSim(alpha.test, s1.metadata, "H", "obese", n=1000)
NonParaPowerSim <- function(TestData, demo, Var1, Var2, n=1000){
  keep <- which(demo[Var2] == "Yes")
  TestData <- TestData[Var1]
  VarOfInt <- TestData[, 1]
  set.seed(3)
  PvalSim <- replicate(n, wilcox.test(
    rnorm(length(VarOfInt[keep]), mean(VarOfInt[keep]), sd(VarOfInt[keep])), 
    rnorm(length(VarOfInt[-keep]), mean(VarOfInt[-keep]), 
          sd(VarOfInt[-keep])))$p.value)
  return((sum(PvalSim < 0.05)) / 10)
}

## Function to graph OTU versus BMI linear regression with log scale option
## Default for log is FALSE. Necessary library is ggplot2. 
## Example on how to use:
## q1 <- ggBaseLinear(RossMicrobiome, RossDemo, "Otu00102", log=TRUE)
ggBaseLinear <- function(microbiome, demo, otuNum, log=FALSE){
  
  # Check for needed library and load if necessary
  deps = c("ggplot2");
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
  
  # Run Log or normal linear regression
  if(log==FALSE){
    OTUData <- as.data.frame(microbiome[otuNum])
  } else{
    OTUData <- as.data.frame(log(microbiome[otuNum] + 1))
  }
  bmi <- demo$BMI
  combined <- cbind(OTUData, bmi)
  ggplot(combined, aes(x=combined[, 1], y=combined[, 2])) + 
    geom_point(shape=1) + 
    geom_smooth(method=lm) + xlab(otuNum) + ylab(colnames(combined)[2])
}


## Function to graph OTU versus BMI  boxplot with log scale option
## Default for log is FALSE. Necessary library is ggplot2 
## Example on how to use:
## q1 <- ggBaseLinear(RossMicrobiome, RossDemo, "Otu00102", "Ross", log=TRUE)
ggBaseBox <- function(microbiome, demo, otuNum, title, log=FALSE){
  # Check for needed libraries and load if necessary
  deps = c("ggplot2", "RColorBrewer");
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
  
  # Run log or normal box plot graph
  if(log==FALSE){
    OTUData <- as.data.frame(microbiome[otuNum])
  } else{
    OTUData <- as.data.frame(log(microbiome[otuNum] + 1))
  }
  obese <- as.factor(demo$obese)
  combined <- cbind(OTUData, obese)
  ggplot(combined, aes(obese, combined[, 1])) + 
    geom_boxplot(aes(fill = obese)) + 
    scale_fill_brewer(palette="Set1") + 
    xlab(colnames(combined)[2]) + ylab(colnames(combined)[1]) + 
    ggtitle(title)
  
}


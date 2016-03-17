## Functions For Meta Analysis Paper on
## Obesesity and the bacterial microbiome
## Marc Sze
## March 16, 2016


## Function to graph OTU versus BMI linear regression with log scale option
## Default for log is FALSE. 
## Example:
## q1 <- ggBaseLinear(RossMicrobiome, RossDemo, "Otu00102", log=TRUE)
ggBaseLinear <- function(microbiome, demo, otuNum, log=FALSE){
  deps = c("ggplot2");
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
  
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
## Default for log is FALSE. 
## Example:
## q1 <- ggBaseLinear(RossMicrobiome, RossDemo, "Otu00102", "Ross", log=TRUE)
ggBaseBox <- function(microbiome, demo, otuNum, title, log=FALSE){
  deps = c("ggplot2", "RColorBrewer");
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
  
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


## Power Calculations Variations: Estimating needed n
## Obesesity and the bacterial microbiome
## Marc Sze
## March 19, 2016

library(statmod)
source("code/UsedFunctions.R")

ShannonRRTable <- read.csv("results/tables/denovoShannonRRTable.csv")
BFRatioRRTable <- read.csv("results/tables/denovoBFRatioRRTable.csv")
PowerTable <- read.csv("results/tables/denovoPowerTable.csv")
Study <- as.character(PowerTable$Study)


### Generating predicted n to reach 80% power for Shannon Diversity.

# Optimization for n Baxter
# Use the means and sds for the actual non-obese and obese groups
## Setting variables to fixed amounts execept n1 and n2 for optimization
## Need to manualy adjust parameters for each data set.

#ForOpNonParaPowerSim <- function(nValues){

#  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[1], PowerTable$MeanObeseH[1], 
#                      nValues[1], nValues[2], 
#                      PowerTable$SDNonObeseH[1], PowerTable$SDObeseH[1])

#}

# Varying both n, not really sure I like this approach for this question
#optim(par = c(15,15), ForOpNonParaPowerSim, method="BFGS", 
#      control=list(maxit=10000, fnscale=-1, trace=TRUE))

# Hold control n constant and vary obese n to give an idea how many obese needed 
## In optimize function maximum is the minimum n needed to get the max power (objective)
## If n needed went above 2000 increased control n by 100 until a reasonable number was 
## found

# For current control n
ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[1], PowerTable$MeanObeseH[1], 
                      n1, ShannonRRTable$tnegH[1] + ShannonRRTable$cnegH[1], 
                      PowerTable$SDNonObeseH[1], PowerTable$SDObeseH[1])
  
}

optimize(ForOpNonParaPowerSim, lower=10, upper=100, maximum=TRUE)
# Max = 89.5, $Objective = 80.3, | control n= 125, obese = 47


# Optimization for n Ross
# Not possible to get a p=0.8 with a reasonable n with study control n (greater than 2000)
# Thus set control n at 100 and looked for obese n to get to a p=0.8

# For current control n
# obese n > 2000
ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[2], PowerTable$MeanObeseH[2], 
                      n1, ShannonRRTable$tnegH[2] + ShannonRRTable$cnegH[2], 
                      PowerTable$SDNonObeseH[2], PowerTable$SDObeseH[2])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)

# For control n = 100
ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[2], PowerTable$MeanObeseH[2], 
                      n1, 100, 
                      PowerTable$SDNonObeseH[2], PowerTable$SDObeseH[2])
  
}

optimize(ForOpNonParaPowerSim, lower=200, upper=260, maximum=TRUE)
# Max = 222, $Objective = 80.4, | control n= set to 100, obese = 38

# Optimization for n Goodrich
## In optimize function maximum is the minimum n needed to get the max power (objective)
## If n needed went above 2000 increased control n by 100 upto 1000.  The total
## control n was not increased once power=0.8 and obese n < 2000 or if 1000 control n
## was reached first.

# For current control n
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[3], PowerTable$MeanObeseH[3], 
                      n1, ShannonRRTable$tnegH[3] + ShannonRRTable$cnegH[3], 
                      PowerTable$SDNonObeseH[3], PowerTable$SDObeseH[3])
  
}

optimize(ForOpNonParaPowerSim, lower=1500, upper=2000, maximum=TRUE)
# Max = 1788, $Objective = 36.6, | control n= 404, obese = 103
# To get to 80 with the given numbers you would need more than 2000 obese individuals

# For control n = 500
# obese n > 2000

# For control n = 600
# obese n > 2000

# For control n = 700
# obese n > 2000

# For control n = 800
# obese n > 2000

# For control n = 900
# obese n > 2000

# For control n = 1000
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[3], PowerTable$MeanObeseH[3], 
                      n1, 1000, 
                      PowerTable$SDNonObeseH[3], PowerTable$SDObeseH[3])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)



# Optimization for n Escobar
# Not possible to get a p=0.8 with a reasonable n with study control n
# Thus set control n at 100 at looked for obese n to get to a p=0.8

# For current control n
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[4], PowerTable$MeanObeseH[4], 
                      n1, ShannonRRTable$tnegH[4] + ShannonRRTable$cnegH[4], 
                      PowerTable$SDNonObeseH[4], PowerTable$SDObeseH[4])
  
}

# For control n = 1000
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[4], PowerTable$MeanObeseH[4], 
                      n1, 1000, 
                      PowerTable$SDNonObeseH[4], PowerTable$SDObeseH[4])
  
}


optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# Optimization for n Zupancic
# For current control n
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[5], PowerTable$MeanObeseH[5], 
                      n1, ShannonRRTable$tnegH[5] + ShannonRRTable$cnegH[5], 
                      PowerTable$SDNonObeseH[5], PowerTable$SDObeseH[5])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# For control n = 1000
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[4], PowerTable$MeanObeseH[4], 
                      n1, 1000, 
                      PowerTable$SDNonObeseH[4], PowerTable$SDObeseH[4])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# Optimization for n HMP
# For current control n
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[6], PowerTable$MeanObeseH[6], 
                      n1, ShannonRRTable$tnegH[6] + ShannonRRTable$cnegH[6], 
                      PowerTable$SDNonObeseH[6], PowerTable$SDObeseH[6])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)

# For n = 1000
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[6], PowerTable$MeanObeseH[6], 
                      n1, 1000, 
                      PowerTable$SDNonObeseH[6], PowerTable$SDObeseH[6])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# Optimization for n Wu
# For current control n
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[7], PowerTable$MeanObeseH[7], 
                      n1, ShannonRRTable$tnegH[7] + ShannonRRTable$cnegH[7], 
                      PowerTable$SDNonObeseH[7], PowerTable$SDObeseH[7])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# For n = 1000
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[7], PowerTable$MeanObeseH[7], 
                      n1, 1000, 
                      PowerTable$SDNonObeseH[7], PowerTable$SDObeseH[7])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# Optimization for n Turnbaugh
# For current control n
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[8], PowerTable$MeanObeseH[8], 
                      n1, ShannonRRTable$tnegH[8] + ShannonRRTable$cnegH[8], 
                      PowerTable$SDNonObeseH[8], PowerTable$SDObeseH[8])
  
}

optimize(ForOpNonParaPowerSim, lower=2100, upper=2200, maximum=TRUE)


# For n = 1000
# obese n > 2000

ForOpNonParaPowerSim <- function(n1){
  
  DefNonParaPowerSim2(PowerTable$MeanNonObeseH[7], PowerTable$MeanObeseH[7], 
                      n1, 1000, 
                      PowerTable$SDNonObeseH[7], PowerTable$SDObeseH[7])
  
}

optimize(ForOpNonParaPowerSim, lower=1900, upper=2000, maximum=TRUE)
# Using a quick for loop to generate the necessary effect size information
# This would be to generate actual effect sizes observed.

# Read in required data

alpha_power.obs <- read.table("../data/process/alpha_power.observed", 
                              header = T, row.names = 1)
beta_tests.summary <- read.table("../data/process/beta_tests.summary", 
                                 header=T, row.names = 1)
alpha.data.ALL <- read.table("../data/process/alpha.data", header = T)

# Create variables needed to run the for loop
variableTested <- seq(1:6)
variables <- colnames(alpha_power.obs)[1:6]
alpha_effectSize <- matrix(nrow=length(rownames(beta_tests.summary)), ncol=length(variableTested))
studies <- rownames(beta_tests.summary)

# For loops to create an effect size table
for(i in variableTested){
  tempdata <- c()
  for(j in studies){
    tempdata <- rbind(tempdata, abs((mean(alpha.data.ALL[which(alpha.data.ALL$dataset == j & alpha.data.ALL$obese == FALSE), variables[i]]) - mean(alpha.data.ALL[which(alpha.data.ALL$dataset == j & alpha.data.ALL$obese == TRUE), variables[i]]))/mean(alpha.data.ALL[, variables[i]]) * 100))
  }
  alpha_effectSize[, i] <- tempdata
}

colnames(alpha_effectSize) <- variables
rownames(alpha_effectSize) <- studies


write.table(alpha_effectSize, file="../data/process/alpha_EffectSize.summary", 
            quote = F, sep='\t', row.names=T)



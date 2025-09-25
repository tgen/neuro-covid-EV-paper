

library(tidyverse)
library(RnaSeqSampleSize)


countTable <- read_tsv("data/countTable.txt")
sampleData <- read_tsv("data/sampleData.txt")

control <- sampleData$Sample[sampleData$GroupCode == "out_pt_neg"]
non_icu <- sampleData$Sample[sampleData$GroupCode == "non_icu"]
outpt <- sampleData$Sample[sampleData$GroupCode == "out_pt_pos"]
icuA <- sampleData$Sample[sampleData$GroupCode == "icu" & sampleData$Visit == "A"]
icuB <- sampleData$Sample[sampleData$GroupCode == "icu" & sampleData$Visit == "B"]

conSub <- grepl(paste(control, collapse = "|"), names(countTable))
treatSub <- grepl(paste(icuB, collapse = "|"), names(countTable))

dataMatrix <- cbind(countTable[,conSub], countTable[,treatSub])
dataMatrix <- dataMatrix[rowSums(dataMatrix) > 0,]

write_tsv(dataMatrix, "data/dataMatrix.txt")

# dataMatrix is the raw count table with 10 samples in columns. Genes at rows. First five samples are control and last five samples are treatment.

dataMatrixDistribution <- est_count_dispersion(dataMatrix, group=c(rep(0,5), rep(1,5)))



# n is proposed sample size, f is FDR, rho is minimal fold change. use repNumber=5 for test purpose. Using default for other parameters, total gene m=10,000,  proposed differential gene m1=100.

est_power_distribution(n=114,f=0.05,rho=2,distributionObject=dataMatrixDistribution,repNumber=5, minAveCount = 10, m = 200)



# repNumber=1000 to get estimation results. This is a little slow.

set.seed(123)

est_power_distribution(n=114,f=0.05,rho=2,distributionObject=dataMatrixDistribution,repNumber=1000, minAveCount = 10, m = 200)









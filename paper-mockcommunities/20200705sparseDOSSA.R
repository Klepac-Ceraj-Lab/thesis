library(BiocManager)
library("sparseDOSSA")

n.microbes <- 150
n.samples <- 50

sparseDOSSA::sparseDOSSA( number_features = n.microbes, 
                          number_samples = n.samples)



setwd("/Users/danielle/Documents/thesis/paper-mockcommunities")

mock <- read.csv("SyntheticMicrobiome.csv", header=TRUE)

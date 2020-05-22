library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)

setwd("/Users/danielle/Documents/thesis/analysis")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)


# 8 samples that are replicates of same fecal sample, just different storage method
sampleid <- unique(df$sampleid)

samplename <- gsub("-.*", "", sampleid)
length(samplename)
length(unique(samplename))

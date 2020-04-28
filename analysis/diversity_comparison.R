library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(sna)
library(reshape2)

setwd("/Users/danielle/Documents/thesis/analysis")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

df[is.na(df)] <- 0
df["shannon"] <- diversity(df[,6:307], "shannon")
df$dev_stage[df$AgeMonths <= 15] <- "less than 15 months"
df$dev_stage[df$AgeMonths > 15 & df$AgeMonths <= 30] <- "15 to 30 months"
df$dev_stage[df$AgeMonths > 30] <- "older than 30 months"

df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"),ordered = TRUE)

p1 <- ggplot(df, aes(x = dev_stage, y = shannon)) + geom_boxplot(aes(colour = method))
p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("alpha diversity")

p1

# shannon diversity for different ages based on profiling method
mean(df[df$dev_stage == "less than 15 months" && df$method == "amp",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

# t-test

df_stage1 <- df[df$dev_stage == "less than 15 months",]
t.test(as.numeric(df_stage1$shannon) ~ df_stage1$method)

df_stage2 <- df[df$dev_stage == "15 to 30 months",]
t.test(as.numeric(df_stage2$shannon) ~ df_stage2$method)

df_stage3 <- df[df$dev_stage == "older than 30 months",]
t.test(as.numeric(df_stage3$shannon) ~ df_stage3$method)

anova_dev_stage <- aov(as.numeric(df$shannon) ~ df$dev_stage)
summary(anova_dev_stage)
posthoc <- TukeyHSD(anova_dev_stage, 'df$dev_stage', conf.level=0.95)
posthoc

# bray curtis dissimilarity for samples

abund_table <- df[,6:307]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,6:307])!=0)[,2:5]

bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))

# upper left side of matrix is bc differences for mgx samples: mgx samples
mgx_bc <- bc_matrix[1:99,1:99]
mgx_bc[lower.tri(mgx_bc)] <- NA
mgx_bc[mgx_bc == 0.0] <- NA
mgx_bc_vec <- na.omit(as.vector(mgx_bc))
mgx_bc_vec

# lower right side of matrix = BC for 16S samples: 16S samples

amp_bc <- bc_matrix[100:198,100:198]
amp_bc[lower.tri(mgx_bc)] <- NA
amp_bc[amp_bc == 0.0] <- NA
amp_bc_vec <- na.omit(as.vector(amp_bc))
amp_bc_vec

# diagonal of the upper right side of matrix = BC for same samples (16S vs mgx)
paired_bc_vec <- as.vector(diag(bc_matrix[1:99, 100:198]))
length(paired_bc_vec) <- length(amp_bc_vec)


# upper right side of matrix excluding diagonal = BC for different samples (16S vs mgx)
unpaired_bc <- diag.remove(bc_matrix[1:99, 100:198], remove.val=NA)
unpaired_bc[lower.tri(unpaired_bc)] <- NA
unpaired_bc_vec <- na.omit(as.vector(unpaired_bc))

bc_df <- data.frame(mgx_bc_vec, amp_bc_vec, paired_bc_vec, unpaired_bc_vec)

p2 <- ggplot(melt(bc_df), aes(variable, value)) + geom_boxplot(aes(fill=variable))
p2

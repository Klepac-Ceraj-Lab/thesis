library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(reshape2)

setwd(Sys.getenv("REPOPATH"))

# download data as described in `diabimmune/data/README.md`

# creates "metadata"
load("diabimmune/data/DIABIMMUNE_Karelia_metadata.RData")
# creates data_16s
load("diabimmune/data/DIABIMMUNE_Karelia_16S_data.RData")
data_mgx <- read.csv("diabimmune/data/DIABIMMUNE_karelia_metaphlan_table.txt", sep='\t')

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

df$Escherichia_combined = sum(df$Escherichia, df$Escherichia.Shigella)
df <- subset(df, select = -c(Escherichia,Escherichia.Shigella ))

df$dev_stage[df$AgeMonths <= 15] <- "less than 15 months"
df$dev_stage[df$AgeMonths > 15 & df$AgeMonths <= 30] <- "15 to 30 months"
df$dev_stage[df$AgeMonths > 30] <- "older than 30 months"

df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"),ordered = TRUE)

df_stage12 <- df[df$dev_stage == "less than 15 months" | df$dev_stage ==  "15 to 30 months",]

df_mgx <- df_stage12[df_stage12$method == "mgx",6:306]
df_mgx_na <- Filter(function(x)!all(!is.na(x)), df_mgx)

df_stage1_mgx_not1 <- df_mgx[colMeans(df_mgx, na.rm = TRUE)>0,]

df_stage1_amp <- df[df$dev_stage == "less than 15 months" & df$method == "amp",6:306]
df_stage1_amp_not0 <- df_stage1_amp[df_stage1_amp.columns[df_stage1_amp.mean(axis=0) > 0.2]]
df_stage1_amp_na <- Filter(function(x)!all(!is.na(x)), df_stage1_amp)

bugs_16S <- merge(df_stage1_mgx_na, )
bugs_mgx <- unique(df[(!is.na(abund$mgx_abund)) & is.na(df$amplicon_abund),]) 

# looking for kids with greatest bray-curtis dissimilarities

meta_table <- df_stage12[1:38,c(3,5)]

abund_table <- df_stage12[,6:306]
abund_table[is.na(abund_table)] <- 0


bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))

# pairwise comparisons of 16s-amp for same samples

paired_bc_vec <- as.vector(diag(bc_matrix[1:38,39:76]))
bc_dist <- cbind(meta_table, paired_bc_vec)
bc_dist <- bc_dist[order(bc_dist$paired_bc_vec, decreasing = TRUE),]

# samples with top 10 greatest differences
largest_diff <- as.vector(bc_dist[1:10,]$sampleid)

# looking at relative abundance for these samples
largest_diff_df <- df[is.element(df$sampleid, largest_diff),]
largest_diff_df <- subset(largest_diff_df, select = -c(uid, X, AgeMonths, dev_stage))
largest_diff_df[is.na(largest_diff_df)] <-  0

write.csv(largest_diff_df, "youngchildren_bug_diff.csv")

diff_df_melt <- melt(largest_diff_df)

plt <- ggplot(data = diff_df_melt, aes(x = variable, y = value))
plt + geom_boxplot() + theme_minimal() + labs(x = "Title", y = "x")

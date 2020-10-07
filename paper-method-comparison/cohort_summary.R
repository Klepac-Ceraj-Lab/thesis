library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(grid)
library(reshape2)

setwd("/Users/danielle/Documents/thesis/paper-abundance-tables")

df <- read.csv("mgx_abund_df.csv", header=TRUE)

# calculating shannon diversity

df["shannon"] <- diversity(df[,4:143], "shannon")

df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"), ordered = TRUE)
# age statistics

df_stage1 <- df[df$dev_stage == "less than 15 months",]
df_stage2 <- df[df$dev_stage == "15 to 30 months",]
df_stage3 <- df[df$dev_stage == "older than 30 months",]

# number of kids in each age group

nrow(df_stage1) # 85
nrow(df_stage2) # 15
nrow(df_stage3) # 60

nrow(df)

# age statistics

mean(df_stage1$AgeMonths)
sqrt(var(df_stage1$AgeMonths))

mean(df_stage2$AgeMonths)
sqrt(var(df_stage2$AgeMonths))

mean(df_stage3$AgeMonths)
sqrt(var(df_stage3$AgeMonths))

# only counting 1 profiling method (to get number of kids, not profiles)

ageplot <- ggplot(df, aes(x = AgeMonths, fill = dev_stage)) + 
  geom_histogram(data = df, 
                 position = "identity", 
                 alpha = 0.7, 
                 binwidth = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + 
  ylab("# of samples") + 
  xlab("Age (months)") +
  labs(title = "", tag = "") + 
  theme(legend.position = c(0.8, 0.8))

ageplot


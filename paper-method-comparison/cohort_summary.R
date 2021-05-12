library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(grid)
library(reshape2)

setwd("/Users/danielle/Documents/thesis/paper-method-comparison")
df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

total_columns <- ncol(df)
total_samples <- nrow(df)


# adding developmental stage

df$dev_stage[df$AgeMonths<15 ] <- "less than 15 months"
df$dev_stage[(df$AgeMonths>=15 & df$AgeMonths <30) ] <- "15 to 30 months"
df$dev_stage[df$AgeMonths>30 ] <- "older than 30 months"


df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"), ordered = TRUE)
# number of kids in each age group

nrow(df[df$dev_stage == "less than 15 months",])/2 # 204 kids <15 months
nrow(df[df$dev_stage == "15 to 30 months",])/2 # 41 kids >15 months, <30 months
nrow(df[df$dev_stage == "older than 30 months",])/2 # 193 >30 months

# age statistics, calculating mean and standard deviation of ages in each
# developmental stage

mean(df[df$dev_stage == "less than 15 months",]$AgeMonths)
sqrt(var(df[df$dev_stage == "less than 15 months",]$AgeMonths))

mean(df[df$dev_stage == "15 to 30 months",]$AgeMonths)
sqrt(var(df[df$dev_stage == "15 to 30 months",]$AgeMonths))

mean(df[df$dev_stage == "older than 30 months",]$AgeMonths)
sqrt(var(df[df$dev_stage == "older than 30 months",]$AgeMonths))

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


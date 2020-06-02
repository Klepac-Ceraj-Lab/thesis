library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(grid)
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


# 8 samples that are replicates of same fecal sample, just different storage method
sampleid <- unique(df$sampleid)

samplename <- gsub("-.*", "", sampleid)
length(samplename)
length(unique(samplename))

# identifying samples that occured more than once
n_occur <- data.frame(table(samplename))
duplicates <- n_occur[n_occur$Freq > 1,]

startsWith(duplicates$samplename, sampleid)

sampleid.exists(startsWith(duplicates$samplename))

# age statistics

df_stage1 <- df[df$dev_stage == "less than 15 months",]
mean(df_stage1[df_stage1$method == "mgx", ]$shannon)
mean(df_stage1[df_stage1$method == "amp", ]$shannon)
t.test(as.numeric(df_stage1$shannon)~ df_stage1$method)

df_stage2 <- df[df$dev_stage == "15 to 30 months",]
kruskal.test(as.numeric(df_stage2$shannon) ~ df_stage2$method)

df_stage3 <- df[df$dev_stage == "older than 30 months",]
kruskal.test(as.numeric(df_stage3$shannon) ~ df_stage3$method)

anova_dev_stage <- aov(as.numeric(df$shannon) ~ df$dev_stage)
summary(anova_dev_stage)
posthoc <- TukeyHSD(anova_dev_stage, 'df$dev_stage', conf.level=0.95)
posthoc

# number of kids in each age group

nrow(df_stage1)
nrow(df_stage2)
nrow(df_stage3)

nrow(df)

# age statistics

mean(df_stage1$AgeMonths)
sqrt(var(df_stage1$AgeMonths))

mean(df_stage2$AgeMonths)
sqrt(var(df_stage2$AgeMonths))

mean(df_stage3$AgeMonths)
sqrt(var(df_stage3$AgeMonths))

# only counting 1 profiling method (to get number of kids, not profiles)
df_stage1_mgx <- df_stage1[df_stage1$method == "mgx",]
df_stage2_mgx <- df_stage2[df_stage2$method == "mgx",]
df_stage3_mgx <- df_stage3[df_stage3$method == "mgx",]

df_ages <- df[,c("sampleid","method","dev_stage", "AgeMonths")]
df_ages <- df_ages[df_ages$method == "mgx",]


ageplot<- ggplot(df, aes(x = AgeMonths)) + 
  geom_histogram(data=df_stage1_mgx, fill = "#481567FF", alpha=0.7, binwidth = 10)+
  geom_histogram(data=df_stage2_mgx, fill = "#238A8DFF", alpha=0.7, binwidth = 10)+
  geom_histogram(data=df_stage3_mgx, fill = "#FDE725FF", alpha=0.7, binwidth = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="bottom") + ylab("# of samples") + xlab("Age (months)") +
  labs(title = "", tag = "B")+ theme(legend.position = c(0.5, 0.5))+
  scale_color_manual(values=c("less than 15 months"="#481567FF", "15 to 30 months"="#238A8DFF",
                              "older than 30 months"="#FDE725FF"))

ageplot <- ggplot(df_ages, aes(x = AgeMonths, fill = dev_stage)) + 
  geom_histogram(data = df_ages, position = "identity", alpha = 0.7, binwidth = 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + ylab("# of samples") + xlab("Age (months)") +
  labs(title = "", tag = "B")+ theme(legend.position = c(0.8, 0.8))

ageplot


library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)


setwd("/Users/danielle/Documents/thesis/subsampled_analysis")

df <- read.csv("subsampled_df.csv", header=TRUE)
df["shannon"] <- diversity(df[,8:95], "shannon")
df["evenness"] <- diversity(df[,8:95], "simpson")
df["richness"] <- apply(df[,8:95]>0,1,sum)

plot1 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(shape=dev_stage, color = sampling_cat))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="richness", y="evenness", 
       color="read depth", shape="developmental stage")
plot1

plot2 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(shape=sampling_cat, color = sampleid)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="richness", y="evenness", color="sample", shape = "sampling depth")

plot2

# pcoa plots

abund <- df[,8:95] # only contains abundances, no metadata
abund <- abund[apply(abund[,-1], 1, function(x) !all(x==0)),]
ord <-  metaMDS(abund, distance = "bray")

plot(ord,type = "n")
points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow",
       col = df[sampleid])

varespec %>%
  +     metaMDS(trace = F) %>%
  +     ordiplot(type = "none") %>%
  +     text("sites")


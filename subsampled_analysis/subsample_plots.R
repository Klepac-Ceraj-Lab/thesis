library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)


setwd("/Users/danielle/Documents/thesis/subsampled_analysis")

df <- read.csv("subsampled_df.csv", header=TRUE)
df["shannon"] <- diversity(df[,8:95], "shannon")
df["evenness"] <- diversity(df[,8:95], "simpson")
df["richness"] <- apply(df[,8:95]>0,1,sum)

# subsample children statistics

mean(df[df$dev_stage == "less than 15 months",]$AgeMonths)
nrow(df[df$dev_stage == "less than 15 months",])
nrow(df[df$dev_stage == "15 to 30 months",])
nrow(df[df$dev_stage == "older than 30 months",])


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

# pcoa plots with vegan

meta <-df[,2:7]
abund <- df[,8:95] # only contains abundances, no metadata
abund <- abund[apply(abund[,-1], 1, function(x) !all(x==0)),]
ord <-  metaMDS(abund, distance = "bray")

mod <- rda(abund, scale = TRUE)
plot(mod)
with(meta, levels(sampling_cat))
colvec <- c("red2", "green4", "mediumblue", "yellow")
with(meta, points(mod, display = "sites", col = colvec[sampling_cat],
                      scaling = scl, pch = 21, bg = colvec[sampling_cat]))
with(meta, legend("topright", legend = levels(sampling_cat), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec))

# trying again to make pcoa plots

abund_table <- df[,8:95]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,8:95])!=0)[,2:7]

sol<-cca(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, meta_table$sampling_cat, meta_table$dev_stage)
colnames(df_sites)<-c("x,","y","sampleid", "sampling_cat", "dev_stage")

# coloring by sampleid
p1<-ggplot()
p1<-p1+geom_point(data=df_sites,aes(x,y,colour=sampleid))+theme_bw()
p1

# coloring by sampling depth
p2 <- ggplot()
p2 <- p2+geom_point(data=df_sites,aes(x,y,colour=sampleid, shape = dev_stage, 
      size = 0.1, stroke = 0)) +
  theme_bw()+labs(x="axis 1, 14.22%", y="axis 2, 11.79%", 
 color="sample id", shape = "developmental stage")
p2

# coloring by dev stage and sampling depth
p3 <- ggplot()
p3 <- p3+geom_point(data=df_sites,aes(x,y,colour=sampling_cat, shape = dev_stage, 
                                      size = 0.1, stroke = 0)) +
  theme_bw()+labs(x="axis 1, 14.22%", y="axis 2, 11.79%", 
                  color="sampling depth", shape = "developmental stage")
p3

# calculating statistics

mean((df[df$sampling_cat == "10000.0" ,]$richness), na.rm = TRUE)
IQR((df[df$sampling_cat == "10000.0" ,]$richness), na.rm = TRUE)
mean((df[df$sampling_cat == "original depth",]$richness), na.rm = TRUE)
IQR((df[df$sampling_cat == "original depth",]$richness), na.rm = TRUE)

anova <- aov(df$richness~df$sampling_cat)
summary(anova)
posthoc <- TukeyHSD(anova,conf.level=0.95)
posthoc
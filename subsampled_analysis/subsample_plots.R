library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(ggpubr)


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

plot0 <- ggplot(df, aes(sampling_cat, shannon)) + 
  geom_point(aes(shape=dev_stage, color = sampling_cat))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="sampling depth", y="Shannon Index", 
       color="read depth", shape="developmental stage")+
  labs(title = "", tag = "A")
plot0


plot1 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(shape=dev_stage, color = sampling_cat))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="richness", y="evenness", 
       color="read depth", shape="developmental stage")+
  labs(title = "", tag = "A")
plot1

plot2 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(shape=sampling_cat, color = sampleid)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="richness", y="evenness", color="sample", shape = "sampling depth")

plot2


gl <- list(plot1, p3, plot2)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,2, 2),
                        c(3, 2,2))
)

# pcoa plots
# http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

abund_table <- df[,8:97]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,8:97])!=0)[,2:7]

sol<-cca(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, meta_table$sampling_cat, meta_table$dev_stage)
colnames(df_sites)<-c("x,","y","sampleid", "sampling_cat", "dev_stage")

# coloring by dev stage and sampling depth

df_sites$dev_stage <- factor(df_sites$dev_stage,
                            levels = c("less than 15 months", 
                                       "15 to 30 months", 
                                       "older than 30 months"),ordered = TRUE)


p1<-ggplot(df_sites, aes(df_sites$x, df_sites$y, colour=sampling_cat, shape = dev_stage)) +
  geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Axis 1, 14.22%", y="Axis 2, 11.79%", 
                  color="sampling depth", shape = "developmental stage")+ 
  theme(legend.position = "none")+
  labs(title = "", tag = "B")
p1

# calculating statistics
my_comparisons <- list( c("10000.0", "100000.0"), c("100000.0", "1000000.0"), c("1000000.0", "original depth") )


p4 <- ggplot(df, aes(x = sampling_cat, y = richness)) + geom_boxplot()
p4 <- p4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("richness")+ xlab("sampling depth") +
  labs(title = "", tag = "C")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 50)
p4

gl <- list(plot1, p4, p1)

grid.arrange(
  grobs = gl,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2),
                        c(3,2))
)  

mean((df[df$sampling_cat == "10000.0" ,]$richness), na.rm = TRUE)
IQR((df[df$sampling_cat == "10000.0" ,]$richness), na.rm = TRUE)

mean((df[df$sampling_cat == "100000.0" ,]$richness), na.rm = TRUE)
IQR((df[df$sampling_cat == "100000.0" ,]$richness), na.rm = TRUE)

mean((df[df$sampling_cat == "1000000.0" ,]$richness), na.rm = TRUE)
IQR((df[df$sampling_cat == "1000000.0" ,]$richness), na.rm = TRUE)

mean((df[df$sampling_cat == "original depth",]$richness), na.rm = TRUE)
IQR((df[df$sampling_cat == "original depth",]$richness), na.rm = TRUE)

anova <- aov(df$richness~df$sampling_cat)
summary(anova)
posthoc <- TukeyHSD(anova,conf.level=0.95)
posthoc
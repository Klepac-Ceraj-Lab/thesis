library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(ggpubr)

# index of where bug abundances columns start and end
start_bugs <- 7
end_bugs <- 104

setwd("/Users/danielle/Documents/thesis/subsampled_analysis")

df <- read.csv("subsampled_df.csv", header=TRUE)

df[is.na(df)] <- 0
df["shannon"] <- diversity(df[,start_bugs:end_bugs], "shannon")
df["evenness"] <- diversity(df[,start_bugs:end_bugs], "simpson")
df["richness"] <- apply(df[,start_bugs:end_bugs]>0,1,sum)

# subsample children statistics

# average age of children less than 15 months old
mean(df[df$dev_stage == "less than 15 months",]$AgeMonths)

# number of kids in each age group
nrow(df[df$dev_stage == "less than 15 months",])
nrow(df[df$dev_stage == "15 to 30 months",])
nrow(df[df$dev_stage == "older than 30 months",])

# plots evenness and richness by developmental stage and read depth
plot0 <- ggplot(df, aes(sampling_cat, shannon)) + 
  geom_point(aes(shape=dev_stage, color = sampling_cat)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x="sampling depth", y="Shannon Index", 
       color="read depth", shape="developmental stage") +
  labs(title = "", tag = "A")
plot0

df$dev_stage <- factor(df$dev_stage,
                            levels = c("less than 15 months", 
                                       "15 to 30 months", 
                                       "older than 30 months"),ordered = TRUE)

plot1 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(shape=dev_stage, color = sampling_cat))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
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


gl <- list(plot0, plot1, plot2)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,2, 2),
                        c(3, 2,2))
)

# pcoa plots
# http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

abund_table <- df[,start_bugs:end_bugs]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,start_bugs:end_bugs])!=0)[,1:start_bugs-1]

sol<-rda(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, meta_table$sampling_cat, meta_table$dev_stage)
colnames(df_sites)<-c("x,","y","sampleid", "sampling_cat", "dev_stage")

# finding main axes
axis.expl <- function(mod, axes = 1:2) {
  
  if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
  
}

axis.expl(sol)

# coloring by dev stage and sampling depth



df_sites$dev_stage <- factor(df_sites$dev_stage,
                            levels = c("less than 15 months", 
                                       "15 to 30 months", 
                                       "older than 30 months"),ordered = TRUE)


p1<-ggplot(df_sites, aes(df_sites$x, df_sites$y, colour=sampling_cat, shape = dev_stage)) +
  geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="RDA 1, 37.01%", y="RDA 2, 21.67%", 
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


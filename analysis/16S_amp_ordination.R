library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)

setwd("/Users/danielle/Documents/thesis/paper-abundance-tables")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

abund_table <- df[,5:209]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)

meta_table <- subset(df,rowSums(df[,5:209])!=0)[,1:4]

sol<-rda(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, meta_table$AgeMonths, meta_table$dev_stage, meta_table$method)
colnames(df_sites)<-c("x","y","sampleid", "AgeMonths", "dev_stage", "method")

df_sites$dev_stage <- factor(df_sites$dev_stage,
                               levels = c("less than 15 months", 
                                          "15 to 30 months", 
                                          "older than 30 months"),ordered = TRUE)

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


# coloring by sampling depth
p1 <- ggplot(data = df_sites, aes(x,y,colour=dev_stage, shape = method, group = sampleid))
p1 <- p1+geom_point(aes(colour=dev_stage, shape = method, group = sampleid), size = 3, alpha = 0.7) + geom_line() +
theme_bw()+labs(x="RDA 1, 35.21% ", y="RDA 2, 21.84%", color="developmental stage", shape = "profiling method")+
theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
axis.text.x = element_blank(), axis.text.y = element_blank()) 

p1


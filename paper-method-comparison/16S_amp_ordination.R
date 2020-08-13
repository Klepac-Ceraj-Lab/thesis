library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)

setwd("/Users/danielle/Documents/thesis/paper-abundance-tables")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

abund_table <- df[,5:206]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)

meta_table <- subset(df,rowSums(df[,5:206])!=0)[,1:4]

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
pcoa1 <- ggplot(data = df_sites, aes(x,y,colour=dev_stage, shape = method, group = sampleid))
pcoa1 <- pcoa1+geom_point(aes(colour=dev_stage, shape = method, group = sampleid), size = 3, alpha = 0.7) + 
geom_line(size = 0.5) +
theme_bw()+labs(x="RDA 1, 35.21% ", y="RDA 2, 21.84%", color="developmental stage", shape = "profiling method")+
theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(legend.position = c(0.2, 0.8)) + theme(legend.title=element_blank()) + 
  labs(title = "", tag = "D")

pcoa1

# only amp profiles (to layer for powerpoint presentation)

amp_df_sites <- df_sites[df_sites$method == "amp",]

pcoa2 <- ggplot(data = amp_df_sites, aes(x,y,colour=dev_stage, group = sampleid))
pcoa2 <- pcoa2+geom_point(aes(colour=dev_stage, shape = method, group = sampleid), size = 3, alpha = 0.7) + geom_line() +
  theme_bw()+labs(x="RDA 1, 35.21% ", y="RDA 2, 21.84%", color="developmental stage", shape = "profiling method")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  xlim(-0.5, 0.5)+ ylim(-0.5, 1)

pcoa1

grid.arrange(pcoa1, pcoa2, ncol = 2)


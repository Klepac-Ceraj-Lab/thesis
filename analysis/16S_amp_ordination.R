library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)

setwd("/Users/danielle/Documents/thesis/analysis")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

# combining Escherichia and E. Shigella into 1

df$Escherichia_combined = sum(df$Escherichia, df$Escherichia.Shigella)
df <- subset(df, select = -c(Escherichia,Escherichia.Shigella ))

df[is.na(df)] <- 0

abund_table <- df[,6:306]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)

meta_table <- subset(df,rowSums(df[,6:306])!=0)[,2:5]
meta_table$dev_stage[meta_table$AgeMonths <= 15] <- "less than 15 months"
meta_table$dev_stage[meta_table$AgeMonths > 15 & meta_table$AgeMonths <= 30] <- "15 to 30 months"
meta_table$dev_stage[meta_table$AgeMonths > 30] <- "older than 30 months"


sol<-cca(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, meta_table$method, meta_table$dev_stage)
colnames(df_sites)<-c("x","y","sampleid", "method", "dev_stage")

df_sites$dev_stage <- factor(df_sites$dev_stage,
                               levels = c("less than 15 months", 
                                          "15 to 30 months", 
                                          "older than 30 months"),ordered = TRUE)

# finding main axes



# coloring by sampling depth
p1 <- ggplot(data = df_sites, aes(x,y,colour=dev_stage, shape = method, group = sampleid))
#p1 <- p1+geom_point(aes(colour=dev_stage, shape = method, group = sampleid), size = 3, alpha = 0.7) +
#  theme_bw()+labs(x="Axis 1, 9.775% ", y="Axis 2, 8.815%", color="developmental stage", shape = "profiling method")+
#  xlim(-1.5, 0.5) + ylim(-1.5, 3)+ theme_bw() + 
#  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#axis.text.x = element_blank(), axis.text.y = element_blank())


p1 <- p1+geom_point(size = 3, alpha = 0.7) + geom_line() +
  labs(x="Axis 1, 9.775% ", y="Axis 2, 8.815%", color="developmental stage", shape = "profiling method")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1


axis.text.x = element_blank(), axis.text.y = element_blank()
xlim(-1.5, 1.5) + ylim(-1.5, 0.5) +

# triangle without line, don't know why...all sampleid's are matching pairs
df_mgx <- subset(df, df$method=="mgx")
df_amp <- subset(df, df$method == "amp")

mgx_unique <- unique(df_mgx$sampleid)
amp_unique <- unique(df_amp$sampleid)

mgx_unique == amp_unique
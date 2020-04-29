library(ggplot2)
library(vegan)

# looking at shannnon data from cohort

setwd("/Users/danielle/Documents/thesis/theoretical")
babies <- read.csv("theoretical_babies_df.csv")

babies <- na.omit(babies)

babies["shannon"] <- diversity(babies[,6:134], index="shannon")

babies$dev_stage <- factor(babies$dev_stage,
                             levels = c("less than 15 months", 
                                        "15 to 30 months", 
                                        "older than 30 months"),ordered = TRUE)

# plotting age vs Shannon diversity


plot<- ggplot(babies, aes(log(AgeMonths), shannon))

plot1 <- plot + geom_point(aes(color = dev_stage),alpha = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") + ylab("Shannon diversity")
plot1

# boxplots of Shannon diversity given developmental stage

plot2 <- ggplot(babies, aes(x = dev_stage, y = shannon)) + geom_boxplot(aes(colour=dev_stage))
plot2 <- plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = "none") +
  ylab("Shannon diversity") + xlab("developmental stage") 

plot2

mean(babies[babies$dev_stage == "less than 15 months",]$shannon)
mean(babies[babies$dev_stage == "15 to 30 months",]$shannon)
mean(babies[babies$dev_stage == "older than 30 months",]$shannon)

# bray curtis distance matrix

babies2 <- na.omit(babies2)

abund_table <- babies2[,6:134]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)

meta_table <- subset(babies2,rowSums(babies2[,6:134])!=0)[,1:5]

sol<-cca(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, meta_table$dev_stage)
colnames(df_sites)<-c("x","y","sampleid", "dev_stage")

df_sites$dev_stage <- factor(df_sites$dev_stage,
                           levels = c("less than 15 months", 
                                      "15 to 30 months", 
                                      "older than 30 months"),ordered = TRUE)

# CCA plot of all samples

plot3 <- ggplot()
plot3 <- plot3+geom_point(data=df_sites, aes(x,y,colour=dev_stage), alpha = 0.5) +
  theme_bw()+labs(x="Axis 1, 13.52% ", y="Axis 2, 11.32", color="developmental stage")+
  xlim(-1.5, 1.0) + ylim(-1.25, 1.5)+theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot3

gl <- list(plot1, plot3)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,1,NA,NA),
                        c(3, 3,3,NA))
)  



library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(sna)
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

p1 <- ggplot(df, aes(x = dev_stage, y = shannon)) + geom_boxplot(aes(colour = method))
p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Shannon diversity")+labs(title = "", tag = "A") + xlab("developmental stage") +
  theme(legend.position = c(0.8, 0.2)) + theme(legend.title=element_blank())
p1

# shannon diversity for different ages based on profiling method
mean(df[df$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)


# bray curtis dissimilarity for samples

abund_table <- df[,6:307]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,6:307])!=0)[,2:5]

bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))

# upper left side of matrix is bc differences for mgx samples: mgx samples
mgx_bc <- bc_matrix[1:99,1:99]
mgx_bc[lower.tri(mgx_bc)] <- NA
mgx_bc[mgx_bc == 0.0] <- NA
mgx_bc_vec <- na.omit(as.vector(mgx_bc))
mgx_bc_vec

# lower right side of matrix = BC for 16S samples: 16S samples

amp_bc <- bc_matrix[100:198,100:198]
amp_bc[lower.tri(mgx_bc)] <- NA
amp_bc[amp_bc == 0.0] <- NA
amp_bc_vec <- na.omit(as.vector(amp_bc))
amp_bc_vec

# diagonal of the upper right side of matrix = BC for same samples (16S vs mgx)
paired_bc_vec <- as.vector(diag(bc_matrix[1:99, 100:198]))
length(paired_bc_vec) <- length(amp_bc_vec)


# upper right side of matrix excluding diagonal = BC for different samples (16S vs mgx)
unpaired_bc <- diag.remove(bc_matrix[1:99, 100:198], remove.val=NA)
unpaired_bc[lower.tri(unpaired_bc)] <- NA
unpaired_bc_vec <- na.omit(as.vector(unpaired_bc))

bc_df <- data.frame(paired_bc_vec, unpaired_bc_vec)

# for paper, add color to first two boxplots in figure
p2 <- ggplot(melt(bc_df), aes(variable, value)) + geom_boxplot()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("between sample beta diversity")+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("16S & mgx, paired", "16S & mgx, unpaired"))+labs(title = "", tag = "C")
  
p2

# statistics

ks.test(bc_df$mgx_bc_vec, bc_df$amp_bc_vec)

# no significant difference between paired and unpaired samples
ks.test(bc_df$paired_bc_vec, bc_df$unpaired_bc_vec) 

# BC distance for kids of different ages

# less than 15, 16S
amp_under15 <- df[df$dev_stage == "less than 15 months"
                  & df$method == "amp",]
amp.15.table <- amp_under15[,6:307]
amp.15.table<-subset(amp.15.table,rowSums(amp.15.table)!=0)
amp.15.matrix <- as.matrix(vegdist(amp.15.table, "bray", diag=FALSE,upper=FALSE))
amp.15.matrix[amp.15.matrix == 0.0] <- NA
amp.15.vec <- na.omit(as.vector(amp.15.matrix))
amp.15.vec
amp.15.df <- data.frame(amp.15.vec)
colnames(amp.15.df)<-c("BC_dist")
amp.15.df["method"] <- "amp"
amp.15.df["dev_stage"] <- "less than 15 months"
  
mgx_under15 <- df[df$dev_stage == "less than 15 months"
                                & df$method == "mgx",]
mgx.15.table <- mgx_under15[,6:307]
mgx.15.table<-subset(mgx.15.table,rowSums(mgx.15.table)!=0)
mgx.15.matrix <- as.matrix(vegdist(mgx.15.table, "bray", diag=FALSE,upper=FALSE))
mgx.15.matrix[mgx.15.matrix == 0.0] <- NA
mgx.15.vec <- na.omit(as.vector(mgx.15.matrix))
mgx.15.vec
mgx.15.df <- data.frame(mgx.15.vec)
colnames(mgx.15.df)<-c("BC_dist")
mgx.15.df["method"] <- "mgx"
mgx.15.df["dev_stage"] <- "less than 15 months"


amp_15to30 <- df[df$dev_stage == "15 to 30 months"
                & df$method == "amp",]
amp.1530.table <- amp_15to30[,6:307]
amp.1530.table <-subset(amp.1530.table,rowSums(amp.1530.table)!=0)
amp.1530.matrix <- as.matrix(vegdist(amp.1530.table, "bray", diag=FALSE,upper=FALSE))
amp.1530.matrix[amp.1530.matrix == 0.0] <- NA
amp.1530.vec <- na.omit(as.vector(amp.1530.matrix))
amp.1530.vec
amp.1530.df <- data.frame(amp.1530.vec)
colnames(amp.1530.df)<-c("BC_dist")
amp.1530.df["method"] <- "amp"
amp.1530.df["dev_stage"] <- "15 to 30 months"

mgx_15to30 <- df[df$dev_stage == "15 to 30 months"
                 & df$method == "mgx",]
mgx.1530.table <- mgx_15to30[,6:307]
mgx.1530.table <-subset(mgx.1530.table,rowSums(mgx.1530.table)!=0)
mgx.1530.matrix <- as.matrix(vegdist(mgx.1530.table, "bray", diag=FALSE,upper=FALSE))
mgx.1530.matrix[mgx.1530.matrix == 0.0] <- NA
mgx.1530.vec <- na.omit(as.vector(mgx.1530.matrix))
mgx.1530.vec
mgx.1530.df <- data.frame(mgx.1530.vec)
colnames(mgx.1530.df)<-c("BC_dist")
mgx.1530.df["method"] <- "mgx"
mgx.1530.df["dev_stage"] <- "15 to 30 months"
  
amp_over30 <-df[df$dev_stage == "older than 30 months"
                & df$method == "amp",]
amp.over30.table <- amp_over30[,6:307]
amp.over30.table <-subset(amp.over30.table,rowSums(amp.over30.table)!=0)
amp.over30.matrix <- as.matrix(vegdist(amp.over30.table, "bray", diag=FALSE,upper=FALSE))
amp.over30.matrix[amp.over30.matrix == 0.0] <- NA
amp.over30.vec <- na.omit(as.vector(amp.over30.matrix))
amp.over30.vec
amp.over30.df <- data.frame(amp.over30.vec)
colnames(amp.over30.df)<-c("BC_dist")
amp.over30.df["method"] <- "amp"
amp.over30.df["dev_stage"] <- "older than 30 months"

  
mgx_over30 <- df[df$dev_stage == "older than 30 months"
                 & df$method == "mgx",]
mgx.over30.table <- mgx_over30[,6:307]
mgx.over30.table <-subset(mgx.over30.table,rowSums(mgx.over30.table)!=0)
mgx.over30.matrix <- as.matrix(vegdist(mgx.over30.table, "bray", diag=FALSE,upper=FALSE))
mgx.over30.matrix[mgx.over30.matrix == 0.0] <- NA
mgx.over30.vec <- na.omit(as.vector(mgx.over30.matrix))
mgx.over30.vec
mgx.over30.df <- data.frame(mgx.over30.vec)
colnames(mgx.over30.df)<-c("BC_dist")
mgx.over30.df["method"] <- "mgx"
mgx.over30.df["dev_stage"] <- "older than 30 months"
  

bc_method_df <- rbind(amp.15.df, mgx.15.df, amp.1530.df, mgx.1530.df, amp.over30.df, mgx.over30.df)

bc_method_df$dev_stage <- factor(bc_method_df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"),ordered = TRUE)

p3 <- ggplot(bc_method_df, aes(x = dev_stage, y = BC_dist)) + geom_boxplot(aes(colour = method))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("developmental stage")+
  labs(title = "", tag = "B")+
  theme(legend.position = "none")
p3

gl <- list(p1, p3, p2)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,2, 2),
                        c(3, 3,3))
)  
# figure for thesis defense
grid.arrange(
  grobs =gl, 
  layout_matrix = rbind(c(1,2,3)))

# statistics

under_30 <- as.vector(as.matrix(bc_method_df[1:4]))
over_30 <- as.vector(as.matrix(bc_method_df[5:6]))

t.test(under_30, over_30)

mean(na.omit(under_30))
mean(na.omit(over_30))

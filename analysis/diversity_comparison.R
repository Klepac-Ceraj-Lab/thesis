library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(sna)
library(reshape2)
library(gtools)

setwd("/Users/danielle/Documents/thesis/paper-abundance-tables")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

df["shannon"] <- diversity(df[,5:206], "shannon")
df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"),ordered = TRUE)

my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )

p1 <- ggplot(df, aes(x = dev_stage, y = shannon)) + geom_boxplot(aes(colour = method))
p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Shannon diversity")+labs(title = "", tag = "A") + xlab("developmental stage") +
  theme(legend.position = c(0.8, 0.2)) + theme(legend.title=element_blank())+
  stat_compare_means(aes(group = dev_stage, method), label = "p.format")

p1

# shannon diversity for different ages based on profiling method
mean(df[df$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

# statistics

# paired t-test

t.test((df[df$dev_stage == "less than 15 months" & df$method == "amp",]$shannon),
       df[df$dev_stage == "less than 15 months" & df$method == "mgx",]$shannon,paired=TRUE)
t.test((df[df$dev_stage == "15 to 30 months" & df$method == "amp",]$shannon),
       df[df$dev_stage == "15 to 30 months" & df$method == "mgx",]$shannon,paired=TRUE)
t.test((df[df$dev_stage == "older than 30 months" & df$method == "amp",]$shannon),
       df[df$dev_stage == "older than 30 months" & df$method == "mgx",]$shannon,paired=TRUE)


# bray curtis dissimilarity for samples

abund_table <- df[,5:206]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,5:206])!=0)[,1:4]

bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))

# upper left side of matrix is bc differences for mgx samples: mgx samples
mgx_bc <- bc_matrix[1:130,1:130]
mgx_bc[lower.tri(mgx_bc)] <- NA
mgx_bc[mgx_bc == 0.0] <- NA
mgx_bc_vec <- na.omit(as.vector(mgx_bc))
mgx_bc_vec

# lower right side of matrix = BC for 16S samples: 16S samples

amp_bc <- bc_matrix[131:260,131:260]
amp_bc[lower.tri(mgx_bc)] <- NA
amp_bc[amp_bc == 0.0] <- NA
amp_bc_vec <- na.omit(as.vector(amp_bc))
amp_bc_vec

# diagonal of the upper right side of matrix = BC for same samples (16S vs mgx)
paired_bc_vec <- as.vector(diag(bc_matrix[1:130, 131:260]))
length(paired_bc_vec) <- length(amp_bc_vec)

# samples with largest differences
paired_distances <- data.frame(distances = (na.omit(paired_bc_vec)), 
                               sampleid = as.vector(unique(df$sampleid)))
paired_distances$distances <- as.numeric(paired_distances$distances)
                        
# upper right side of matrix excluding diagonal = BC for different samples (16S vs mgx)
unpaired_bc <- diag.remove(bc_matrix[1:130, 131:260], remove.val=NA)
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

  stat_compare_means(method = "t.test", label = "p.signif", label.x = 2, label.y = 1)
  
p2

# statistics

t.test(bc_df$paired_bc_vec, bc_df$unpaired_bc_vec)
mean(bc_df$unpaired_bc_vec) - mean(bc_df$paired_bc_vec, na.rm = TRUE)


# BC distance for kids of different ages

# less than 15, 16S
amp_under15 <- df[df$dev_stage == "less than 15 months"
                  & df$method == "amp",]
amp.15.table <- amp_under15[,5:206]
amp.15.table <- amp.15.table[,colSums(amp.15.table) > 0]
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
mgx.15.table <- mgx_under15[,5:206]
mgx.15.table<- mgx.15.table[,colSums(mgx.15.table) > 0]
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
amp.1530.table <- amp_15to30[,5:206]
amp.1530.table <- amp.1530.table[,colSums(amp.1530.table) > 0]
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
mgx.1530.table <- mgx_15to30[,5:206]
mgx.1530.table <- mgx.1530.table[,colSums(mgx.1530.table) > 0]
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
amp.over30.table <- amp_over30[,5:206]
amp.over30.table <- amp.over30.table[,colSums(amp.over30.table) > 0]
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
mgx.over30.table <- mgx_over30[,5:206]
mgx.over30.table <- mgx.over30.table[,colSums(mgx.over30.table) > 0]
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

t.test(amp.15.vec, mgx.15.vec)
t.test(amp.1530.vec, mgx.1530.vec)
t.test(amp.over30.vec, mgx.over30.vec)


# paired vs unpaired distances for kids over different ages

# 15 months
combined.15 <- smartbind(amp.15.table, mgx.15.table)
combined.15[is.na(combined.15)] <- 0

combined.15.matrix <- as.matrix(vegdist(combined.15, "bray", diag=FALSE,upper=FALSE))
combined.15.matrix[combined.15.matrix == 0.0] <- NA

paired.15 <- as.vector(diag(combined.15.matrix[1:70, 71:140]))

unpaired.15 <- diag.remove(combined.15.matrix[1:70, 71:140], remove.val=NA)
unpaired.15[lower.tri(unpaired.15)] <- NA
unpaired.15 <- na.omit(as.vector(unpaired.15))
              
paired.15.df <- data.frame(paired.15)
colnames(paired.15.df)<-c("BC")
paired.15.df["method"] <- "paired"
paired.15.df["dev_stage"] <- "less than 15 months"

unpaired.15.df <- data.frame(unpaired.15)
colnames(unpaired.15.df)<-c("BC")
unpaired.15.df["method"] <- "unpaired"
unpaired.15.df["dev_stage"] <- "less than 15 months"

length(paired.15) <- length(unpaired.15)

# 15 to 30 months
combined.1530 <- smartbind(amp.1530.table, mgx.1530.table)
combined.1530[is.na(combined.1530)] <- 0

combined.1530.matrix <- as.matrix(vegdist(combined.1530, "bray", diag=FALSE,upper=FALSE))
combined.1530.matrix[combined.1530.matrix == 0.0] <- NA

paired.1530 <- as.vector(diag(combined.1530.matrix[1:15, 16:30]))

unpaired.1530 <- diag.remove(combined.1530.matrix[1:15, 16:30], remove.val=NA)
unpaired.1530[lower.tri(unpaired.1530)] <- NA
unpaired.1530 <- na.omit(as.vector(unpaired.1530))
length(unpaired.1530) <- length(unpaired.15)

paired.1530 <- as.vector(diag(combined.1530.matrix[1:15, 16:30]))
length(paired.1530) <- length(unpaired.15)

paired.1530.df <- data.frame(paired.1530)
colnames(paired.1530.df)<-c("BC")
paired.1530.df["method"] <- "paired"
paired.1530.df["dev_stage"] <- "15 to 30 months"

unpaired.1530.df <- data.frame(unpaired.1530)
colnames(unpaired.1530.df)<-c("BC")
unpaired.1530.df["method"] <- "unpaired"
unpaired.1530.df["dev_stage"] <- "15 to 30 months"

# over 30 months
combined.over30 <- smartbind(amp.over30.table, mgx.over30.table)
combined.over30[is.na(combined.over30)] <- 0

combined.over30.matrix <- as.matrix(vegdist(combined.over30, "bray", diag=FALSE,upper=FALSE))
combined.over30.matrix[combined.over30.matrix == 0.0] <- NA

unpaired.over30 <- diag.remove(combined.over30.matrix[1:45, 46:90], remove.val=NA)
unpaired.over30[lower.tri(unpaired.over30)] <- NA
unpaired.over30 <- na.omit(as.vector(unpaired.over30))
length(unpaired.over30) <- length(unpaired.15)

paired.over30 <- as.vector(diag(combined.over30.matrix[1:45, 46:90]))
length(paired.over30) <- length(unpaired.15)

paired.over30.df <- data.frame(paired.over30)
colnames(paired.over30.df)<-c("BC")
paired.over30.df["method"] <- "paired"
paired.over30.df["dev_stage"] <- "older than 30 months"

unpaired.over30.df <- data.frame(unpaired.over30)
colnames(unpaired.over30.df)<-c("BC")
unpaired.over30.df["method"] <- "unpaired"
unpaired.over30.df["dev_stage"] <- "older than 30 months"


age_bc_df <- rbind(paired.15.df, unpaired.15.df, 
                   paired.1530.df, unpaired.1530.df, 
                   paired.over30.df, unpaired.over30.df)

age_bc_df$dev_stage <- factor(age_bc_df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"),ordered = TRUE)

# for paper, add color to first two boxplots in figure
p4 <- ggplot(age_bc_df, aes(dev_stage, BC)) + geom_boxplot(aes(fill=method))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("between sample beta diversity") +
  scale_fill_manual(values=c("#69b3a2", "grey"))

p4


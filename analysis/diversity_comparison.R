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
  ylab("Shannon diversity")

p1

# shannon diversity for different ages based on profiling method
mean(df[df$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

# statistics

df_stage1 <- df[df$dev_stage == "less than 15 months",]
kruskal.test(as.numeric(df_stage1$shannon)~ df_stage1$method)

df_stage2 <- df[df$dev_stage == "15 to 30 months",]
kruskal.test(as.numeric(df_stage2$shannon) ~ df_stage2$method)

df_stage3 <- df[df$dev_stage == "older than 30 months",]
kruskal.test(as.numeric(df_stage3$shannon) ~ df_stage3$method)

anova_dev_stage <- aov(as.numeric(df$shannon) ~ df$dev_stage)
summary(anova_dev_stage)
posthoc <- TukeyHSD(anova_dev_stage, 'df$dev_stage', conf.level=0.95)
posthoc




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

bc_df <- data.frame(mgx_bc_vec, amp_bc_vec, paired_bc_vec, unpaired_bc_vec)

p2 <- ggplot(melt(bc_df), aes(variable, value)) + geom_boxplot(aes(fill=variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("relative Bray Curtis distance") + xlab("between sample beta diversity")+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c('mgx, all','16S, all', 
    "16S & mgx, paired", "16S & mgx, unpaired"))
  
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
  
mgx_under15 <- df[df$dev_stage == "less than 15 months"
                                & df$method == "mgx",]
mgx.15.table <- mgx_under15[,6:307]
mgx.15.table<-subset(mgx.15.table,rowSums(mgx.15.table)!=0)
mgx.15.matrix <- as.matrix(vegdist(mgx.15.table, "bray", diag=FALSE,upper=FALSE))
mgx.15.matrix[mgx.15.matrix == 0.0] <- NA
mgx.15.vec <- na.omit(as.vector(mgx.15.matrix))
mgx.15.vec
  
amp_15to30 <- df[df$dev_stage == "15 to 30 months"
                & df$method == "amp",]
amp.1530.table <- amp_15to30[,6:307]
amp.1530.table <-subset(amp.1530.table,rowSums(amp.1530.table)!=0)
amp.1530.matrix <- as.matrix(vegdist(amp.1530.table, "bray", diag=FALSE,upper=FALSE))
amp.1530.matrix[amp.1530.matrix == 0.0] <- NA
amp.1530.vec <- na.omit(as.vector(amp.1530.matrix))
amp.1530.vec

mgx_15to30 <- df[df$dev_stage == "15 to 30 months"
                 & df$method == "mgx",]
mgx.1530.table <- mgx_15to30[,6:307]
mgx.1530.table <-subset(mgx.1530.table,rowSums(mgx.1530.table)!=0)
mgx.1530.matrix <- as.matrix(vegdist(mgx.1530.table, "bray", diag=FALSE,upper=FALSE))
mgx.1530.matrix[mgx.1530.matrix == 0.0] <- NA
mgx.1530.vec <- na.omit(as.vector(mgx.1530.matrix))
mgx.1530.vec
  
amp_over30 <-df[df$dev_stage == "older than 30 months"
                & df$method == "amp",]
amp.over30.table <- amp_over30[,6:307]
amp.over30.table <-subset(amp.over30.table,rowSums(amp.over30.table)!=0)
amp.over30.matrix <- as.matrix(vegdist(amp.over30.table, "bray", diag=FALSE,upper=FALSE))
amp.over30.matrix[amp.over30.matrix == 0.0] <- NA
amp.over30.vec <- na.omit(as.vector(amp.over30.matrix))
amp.over30.vec
  
mgx_over30 <- df[df$dev_stage == "older than 30 months"
                 & df$method == "mgx",]
mgx.over30.table <- mgx_over30[,6:307]
mgx.over30.table <-subset(mgx.over30.table,rowSums(mgx.over30.table)!=0)
mgx.over30.matrix <- as.matrix(vegdist(mgx.over30.table, "bray", diag=FALSE,upper=FALSE))
mgx.over30.matrix[mgx.over30.matrix == 0.0] <- NA
mgx.over30.vec <- na.omit(as.vector(mgx.over30.matrix))
mgx.over30.vec

vec_list <- list(amp.15.vec, mgx.15.vec, amp.1530.vec, mgx.1530.vec,amp.over30.vec, mgx.over30.vec )
  
  
bc_method_df <- data.frame(lapply(vec_list, "length<-", max(lengths(vec_list))))
colnames(bc_method_df)<-c("amp, under 15", "mgx, under 15",
                          "amp, 15-30", "mgx, 15-30",
                          "amp, over 30", "mgx, over 30")

p3 <- ggplot(melt(bc_method_df), aes(variable, value)) + geom_boxplot(aes(fill=variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("relative Bray Curtis distance") + xlab("beta diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3

gl <- list(p1, p3, p2)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,2, 2),
                        c(3, 3,3))
)  


# statistics

under_30 <- as.vector(as.matrix(bc_method_df[1:4]))
over_30 <- as.vector(as.matrix(bc_method_df[5:6]))

ks.test(under_30, over_30)

mean(na.omit(under_30))
mean(na.omit(over_30))

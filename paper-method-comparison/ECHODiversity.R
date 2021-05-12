library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(sna)
library(reshape2)
library(gtools)
library(ape)

# This R code file investigates and visualizes the differences in diversity 
# between child gut microbiome samples taxonomically profiled with 16S rRNA 
# amplicon gene and shotgun metagenomics sequencing. Starting with a .csv file 
# containing relative abundances of 130 samples profiled with the two methods, 
# it compares ecological diversity metrics between children of different ages.

# The code investigates alpha (within sample) diversity between the two 
# profiling methods and among the three developmental stages. Next, it analyzes 
# beta (between sample) diversity between the two methods and age groups.

# Afterwards, the code explores the beta diversity between profiles from the 
# same fecal sample (paired) and from random, different fecal samples (unpaired) 
# for each age group. If 16S and shotgun metagenomics profiled microbial 
# communities in the same way, we would expect the paired dissimilarity to be 0. 
# Lastly, the code visualizes the beta diversity between 16S and shotgun 
# metagenomics profiles with a PCoA plot.  

# change working directory if necessary
setwd("/Users/danielle/Documents/thesis/paper-abundance-tables")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

total_columns <- ncol(df)
abund_table <- as.matrix(df[,5:total_columns])


df["shannon"] <- diversity(abund_table, "shannon")
df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"), ordered = TRUE)

p1 <- ggplot(df, aes(x = dev_stage, y = shannon)) + 
  geom_boxplot(aes(colour = method))
p1 <- p1 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("Shannon diversity")+
  labs(title = "", tag = "A") + 
  xlab("developmental stage") +
  theme(legend.position = "none")

p1

# shannon diversity as children age
mean(df[df$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df[df$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

# paired t-test comparing shannon index as children age & between the two methods
t.test((df[df$dev_stage == "less than 15 months" & 
             df$method == "amp",]$shannon),
       df[df$dev_stage == "less than 15 months" & 
            df$method == "mgx",]$shannon,paired=TRUE)
t.test((df[df$dev_stage == "15 to 30 months" & 
             df$method == "amp",]$shannon),
       df[df$dev_stage == "15 to 30 months" & 
            df$method == "mgx",]$shannon,paired=TRUE)
t.test((df[df$dev_stage == "older than 30 months" & 
             df$method == "amp",]$shannon),
       df[df$dev_stage == "older than 30 months" & 
            df$method == "mgx",]$shannon,paired=TRUE) ## this number is in the paper


# bray curtis dissimilarity for samples

# BC distance for kids of different ages
### bug_start = first column of data frame that includes relative abundance data
### bug_end = last column of data frame that incldues relative abundance data
### dev_stage_arg = developmental stage (<15, 15-30, >30 months)
### method_arg = method to profile microbial community (16S or metagenomics)
### returns a data frame of bray-curtis dissimilarity values for each age 
### group/profiling method

calc_bc_dist <- function(bug_start, bug_end, dev_stage_arg, method_arg) {
  df_slice <- df[df$dev_stage == dev_stage_arg
                 & df$method == method_arg,]
  
  abund_table <- df_slice[,bug_start:bug_end]
  abund_table <- abund_table[,colSums(abund_table) > 0]
  bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))
  bc_matrix[bc_matrix == 0.0] <- NA
  bc_vec <- na.omit(as.vector(bc_matrix))
  bc_df <- data.frame(bc_vec)
  colnames(bc_df)<-c("BC_dist")
  bc_df["method"] <- method_arg
  bc_df["dev_stage"] <- dev_stage_arg
  
  return(bc_df)
}

# calculating bray curtis distance amongst samples from same age and profiled 
# with same method
amp_under15 <- calc_bc_dist(5, 206, "less than 15 months", "amp")
mgx_under15 <- calc_bc_dist(5, 206, "less than 15 months", "mgx")
amp_15to30 <- calc_bc_dist(5, 206, "15 to 30 months", "amp")
mgx_15to30 <- calc_bc_dist(5, 206, "15 to 30 months", "mgx")
amp_over30 <- calc_bc_dist(5, 206, "older than 30 months", "amp")
mgx_over30 <- calc_bc_dist(5, 206, "older than 30 months", "mgx")

bc_method_df <- rbind(amp_under15, mgx_under15, amp_15to30, mgx_15to30, 
                      amp_over30, mgx_over30)

bc_method_df$dev_stage <- factor(bc_method_df$dev_stage,
                                 levels = c("less than 15 months", 
                                            "15 to 30 months", 
                                            "older than 30 months"),
                                 ordered = TRUE)

# visualizing bray curtis distances
p2 <- ggplot(bc_method_df, aes(x = dev_stage, y = BC_dist)) + 
  geom_boxplot(aes(colour = method)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("developmental stage") +
  labs(title = "", tag = "B") +
  theme(legend.position = c(0.9, 0.95)) + theme(legend.title=element_blank())
p2


# calculating paired & unpaired Bray Curtis distances for kids of different ages
### see above function for arguments

calc_paired_BC <- function(bug_start, bug_end, dev_stage_arg){
  
  df_amp <- df[df$dev_stage == dev_stage_arg
               & df$method == "amp",]
  df_mgx <- df[df$dev_stage == dev_stage_arg
               & df$method == "mgx",]
  amp_abund <- df_amp[,bug_start:bug_end]
  mgx_abund <- df_mgx[,bug_start:bug_end]
  
  n_samples <- nrow(df_amp)
  
  combined_abund <- smartbind(amp_abund, mgx_abund)
  combined_abund[is.na(combined_abund)] <- 0
  combined_abund_matrix <- as.matrix(vegdist(combined_abund, "bray", 
                                             diag=FALSE,upper=FALSE))
  combined_abund_matrix[combined_abund_matrix == 0.0] <- NA
  
  paired <- as.vector(diag(combined_abund_matrix[1:n_samples, 
                                                 (n_samples+1):(2*n_samples)]))
  unpaired <- diag.remove(combined_abund_matrix[1:n_samples, 
                                                (n_samples+1):(2*n_samples)], 
                          remove.val=NA)
  unpaired[lower.tri(unpaired)] <- NA
  unpaired <- na.omit(as.vector(unpaired))
  
  paired_df <- data.frame(paired)
  colnames(paired_df)<-c("BC_dist")
  paired_df["method"] <- "paired"
  paired_df["dev_stage"] <- dev_stage_arg
  
  unpaired_df <- data.frame(unpaired)
  colnames(unpaired_df)<-c("BC_dist")
  unpaired_df["method"] <- "unpaired"
  unpaired_df["dev_stage"] <- dev_stage_arg
  return(smartbind(paired_df, unpaired_df))
}

paired_under15 <- calc_paired_BC(5, 206, "less than 15 months")
paired_15to30 <- calc_paired_BC(5, 206, "15 to 30 months")
paired_over30 <- calc_paired_BC(5, 206, "older than 30 months")

age_bc_df <- rbind(paired_under15, paired_15to30, paired_over30)
age_bc_df$dev_stage <- factor(age_bc_df$dev_stage,
                              levels = c("less than 15 months", 
                                         "15 to 30 months", 
                                         "older than 30 months"),ordered = TRUE)

# plotting paired bray curtis distances, broken down by age and method
p3 <- ggplot(age_bc_df, aes(dev_stage, BC_dist)) + geom_boxplot(aes(fill=method))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("between sample beta diversity") +
  scale_fill_manual(values=c("#69b3a2", "grey")) + 
  theme(legend.title=element_blank()) + theme(legend.position = c(0.9, 0.95)) + 
  labs(title = "", tag = "C")
p3

# making a PCoA plot 

abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,5:206])!=0)[,1:4]

### Pcoa plot
dist <- vegdist(abund_table,  method = "bray")
PCOA <- pcoa(dist)

pcoa_vectors <- as.data.frame(PCOA$vectors)

pcoa_2_vectors <- cbind.data.frame(pcoa_vectors$Axis.1, pcoa_vectors$Axis.2)
pcoa_plot <- cbind(pcoa_2_vectors, meta_table$sampleid, meta_table$AgeMonths, 
                   meta_table$dev_stage, meta_table$method)

colnames(pcoa_plot)<-c("x","y","sampleid", "AgeMonths", "dev_stage", "method")
pcoa_plot$dev_stage <- factor(df_sites$dev_stage,
                              levels = c("less than 15 months", 
                                         "15 to 30 months", 
                                         "older than 30 months"),ordered = TRUE)

axis1 <- paste("PCoA 1, ", 
               round(PCOA$values$Relative_eig[1],4)*100, "%", sep = "")
axis2 <- paste("PCoA 2, ", 
               round(PCOA$values$Relative_eig[2],4)*100, "%", sep = "")

pcoa1 <- ggplot(data = pcoa_plot, aes(x,y,colour=dev_stage, 
                                      shape = method, 
                                      group = sampleid))

pcoa1 <- pcoa1+geom_point(aes(colour=dev_stage, shape = method, 
                              group = sampleid), size = 3, alpha = 0.7) + 
  geom_line(size = 0.5) + theme_bw()+labs(x=axis1, 
                                          y=axis2, 
                                          color="developmental stage", 
                                          shape = "profiling method") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  theme(legend.position = c(0.9, 0.95),
        legend.text=element_text(size=10),
        legend.title=element_blank()) + 
  labs(title = "", tag = "D") +
  theme(legend.key.size = unit(0.35, "cm"))
pcoa1


# variance explained by first 10 principal components
top_10_components <- (cumsum(PCOA$values$Relative_eig[1:10]))
names(top_10_components) <- seq(1, 10, by=1)

barplot(top_10_components,
        ylab="cummulative percent explained (%)",
        xlab="principal component #")


# aggregating all the plots together

gl <- list(p1, p2, p3, pcoa1)
grid.arrange(grobs = gl,widths = c(2, 1, 1), 
             layout_matrix = rbind(c(1,2, 2), c(3, 4, 4)))

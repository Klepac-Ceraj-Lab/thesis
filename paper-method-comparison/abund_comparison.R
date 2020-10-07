library(ggplot2)
library(gridExtra)
library(formattable)
library(phyloseq)
library(plyr)
library(ggpubr)

setwd("/Users/danielle/Documents/thesis/paper-abundance-tables")

abund <- read.csv("paper_abund_df.csv", na.strings=c("","NA"))
abund$taxa <- as.character(abund$taxa)

all_bugs_16S_df <- abund[,c("taxa","amplicon_abund")]
all_bugs_16S_df["method"] <- "amp"
all_bugs_16S_df <- rename(all_bugs_16S_df, c("amplicon_abund"="abund"))

all_bugs_mgx_df <- abund[,c("taxa","mgx_abund")] 
all_bugs_mgx_df["method"] <- "mgx"
all_bugs_mgx_df <- rename(all_bugs_mgx_df, c("mgx_abund"="abund"))

all_bugs_abund <- rbind(all_bugs_16S_df, all_bugs_mgx_df)

# remove bugs with abundance of 0
all_bugs_abund[all_bugs_abund==0] <- NA 
all_bugs_abund <- all_bugs_abund[!is.na(all_bugs_abund$abund), ]


# plot relative abundance of each method of all bugs
p0 <- ggplot(all_bugs_abund, 
             aes(x = method, 
                 y = log(abund))) + 
  geom_boxplot()

p0 <- p0 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("log(relative abundance)") + 
  xlab("profiling method") +
  theme(text = element_text(size=20)) +
  stat_compare_means(method = "t.test", 
                     label = "p.signif", 
                     label.x = 2, 
                     label.y = 10)
p0

mean(log(all_bugs_abund$abund[all_bugs_abund$method=="amp"]))
mean(log(all_bugs_abund$abund[all_bugs_abund$method=="mgx"]))

all_bugs_abund$log_abund <- log(all_bugs_abund$abund)
t.test(abund~method, data = all_bugs_abund)

# read in files with bugs unique to each method
bugs_16S <- scan("unique_amplicon.txt", 
                 what="character", sep=",", 
                 strip.white = TRUE)

bugs_mgx <- scan("unique_mgx.txt", 
                 what="character", 
                 sep=",", 
                 strip.white = TRUE)

# bugs found in both methods (intersection)
bugs_intersect <- scan("all_taxa.txt", 
                       what = "character", 
                       sep=",", 
                       strip.white = TRUE)

# filter out bugs not unique to either method
filtered_abund <- abund[abund$taxa %in% bugs_16S | abund$taxa %in% bugs_mgx, ]

# plotting average relative abundance of bugs unique to each sample

bugs_16S_df <- filtered_abund[,c("taxa","amplicon_abund")]
bugs_16S_df["method"] <- "amp"
bugs_16S_df <- rename(bugs_16S_df, c("amplicon_abund"="abund"))

bugs_mgx_df <- filtered_abund[,c("taxa","mgx_abund")] 
bugs_mgx_df["method"] <- "mgx"
bugs_mgx_df <- rename(bugs_mgx_df, c("mgx_abund"="abund"))

unique_bugs_abund <- rbind(bugs_16S_df, bugs_mgx_df)

# remove bugs with abundance of 0
unique_bugs_abund[unique_bugs_abund==0] <- NA 
unique_bugs_abund <- unique_bugs_abund[!is.na(unique_bugs_abund$abund), ]

p1 <- ggplot(unique_bugs_abund, aes(x = method, y = log(abund))) + 
  geom_boxplot()

p1 <- p1 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("log(relative abundance)") + xlab("profiling method") +
  theme(text = element_text(size=20)) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 2, 
                     label.y = 10)

p1

mean((unique_bugs_abund$abund[unique_bugs_abund$method=="amp"]))
mean((unique_bugs_abund$abund[unique_bugs_abund$method=="mgx"]))

unique_bugs_abund$log_abund <- log(unique_bugs_abund$abund)
t.test(abund~method, data = unique_bugs_abund)

grid.arrange(p0, p1, ncol = 2)

# histograms of the distribution of relative abundances of bugs unique to each 
# sample
# important to make sure that we understand distribution before applying any 
# statistical tests :)

par(mfrow=c(2,1))
hist(all_bugs_16S_df$abund,
     xlab = "relative abundances of bugs uniquely found by amplicon sequencing",
     main = NULL)
hist(all_bugs_mgx_df$abund,
     xlab = "relative abundances of bugs uniquely found by metagenomic sequencing",
     main = NULL)

# plotting absolute difference in relative abundance by taxa

# calculating mean absolute difference
all_taxa <- unique(abund$taxa)
avg_absdiff <- aggregate(abund$abs_diff, by=list(abund$taxa), 
                         FUN = mean, na.rm=TRUE)

# find largest relative abundance for each sample
abund %>% group_by(taxa) %>% summarise(B = sum(B))

taxa_order <- avg_absdiff[order(avg_absdiff$x, decreasing = TRUE),]
top_taxa <- taxa_order[1:5,1] #taxa with largest differences between 16S and mgx
largest_diff_df <- abund[abund$taxa %in% top_taxa,]

p2 <- ggplot(largest_diff_df, aes(x = taxa, y = abs_diff)) + geom_boxplot()
p2 <- p2 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black")) +
  ylab("absolute difference in relative abundance") +
  labs(title = "", tag = "A")

p2


# plotting total difference in relative abundance by taxa

avg_totdiff <- aggregate(abund$tot_diff, 
                         by=list(abund$taxa), 
                         FUN = mean, na.rm=TRUE)

taxa_order_totaldiff <- avg_totdiff[order(abs(avg_totdiff$x), decreasing = TRUE),]
top_taxa_totaldiff <- taxa_order_totaldiff[1:5,1] #taxa with largest differences between 16S and mgx
largest_totaldiff_df <- abund[abund$taxa %in% top_taxa_totaldiff,]

p3 <- ggplot(largest_totaldiff_df, aes(x = taxa, y = tot_diff)) + 
  geom_boxplot()
p3 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("total difference in relative abundance") +
  labs(title = "", tag = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
p3
grid.arrange(p2, p3, ncol = 2)

# calculate average total difference for Bacteroides
mean(largest_totaldiff_df$tot_diff[largest_totaldiff_df$taxa == "Bacteroides"])

intersect_df <- abund[abund$taxa %in% bugs_intersect,]

# 16S abundance was greater
intersect_df$amplicon_greater[intersect_df$tot_diff > 0 ] <- 1
# mgx abundance was greater
intersect_df$mgx_greater[intersect_df$tot_diff < 0 ] <- 1 

# Mann-whitney test

# remove bugs that are not present in both methods
for (bug in unique(intersect_df$taxa)) {
  bug_df <- subset(intersect_df, intersect_df$taxa == bug)
  if (mean(bug_df$amplicon_abund) == 0 | mean(bug_df$mgx_abund) == 0) {
    intersect_df <- subset(intersect_df, intersect_df$taxa != bug)}}

wilcox_df <- data.frame(matrix(ncol = 4, 
                               nrow = length(unique(intersect_df$taxa))))
colnames(wilcox_df) <- c("taxa","wilcox", "mean_16S", "mean_mgx")
 wilcox_df$taxa <- unique(intersect_df$taxa)

for (bug in wilcox_df$taxa) {
  # add wilcox test results to dataframe
  wilcox_df[wilcox_df$taxa == bug,]$wilcox <-(wilcox.test(
    intersect_df$amplicon_abund[intersect_df$taxa==bug], 
    intersect_df$mgx_abund[intersect_df$taxa==bug])$p.value)
  
  # add mean relative abundance values to dataframe
  wilcox_df[wilcox_df$taxa == bug,]$mean_16S <- 
    mean(intersect_df$amplicon_abund[intersect_df$taxa==bug])
  
  wilcox_df[wilcox_df$taxa == bug,]$mean_mgx <- 
    mean(intersect_df$mgx_abund[intersect_df$taxa==bug])
}

 wilcox16S <- wilcox_df[wilcox_df$mean_16S > wilcox_df$mean_mgx, ]
 wilcox16S <- (wilcox16S[order(wilcox16S$wilcox,
                             decreasing = FALSE),][1:20,c(1:4)])
 row.names(wilcox16S) <- NULL
 formattable(wilcox16S)
 
 wilcoxmgx <- wilcox_df[wilcox_df$mean_mgx > wilcox_df$mean_16S, ]
 wilcoxmgx <- (wilcoxmgx[order(wilcoxmgx$wilcox,
                               decreasing = FALSE),][1:20,c(1:4)])
 row.names(wilcoxmgx) <- NULL
 formattable(wilcoxmgx)
 

 
 
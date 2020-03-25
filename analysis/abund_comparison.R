library(ggplot2)
library(gridExtra)
library(formattable)
library(phyloseq)

setwd("/Users/danielle/Documents/thesis/analysis")

abund <- read.csv("taxa_abundance_comparison.csv")

# find bugs unique to each method
bugs_16S <- unique(abund[(!is.na(abund$amplicon_abund)) & is.na(abund$mgx_abund),]) #159
bugs_mgx <- unique(abund[(!is.na(abund$mgx_abund)) & is.na(abund$amplicon_abund),]) #59

# find taxa unique to each method 
unique_16S <- bugs_16S$taxa
unique_mgx <- bugs_mgx$taxa

# filter out bugs not unique to either method

filtered_abund <- abund[abund$taxa %in% bugs_16S | abund$taxa %in% bugs_mgx, ]

# plotting average relative abundance of bugs unique to each sample
boxplot(log(bugs_16S$amplicon_abund[bugs_16S$amplicon_abund!=0]), 
        log(bugs_mgx$mgx_abund[bugs_mgx$mgx_abund!=0]),
        names = c("16S rRNA", "shotgun metagenomics"),
        xlab="sequencing method",
        ylab="log(relative abundance)")

mean((bugs_16S$amplicon_abund[bugs_16S$amplicon_abund!=0]))
mean((bugs_mgx$mgx_abund[bugs_mgx$mgx_abund!=0]))
t.test(log(bugs_16S$amplicon_abund[bugs_16S$amplicon_abund!=0]), log(bugs_mgx$mgx_abund[bugs_mgx$mgx_abund!=0]))
ks.test(bugs_16S$amplicon_abund[bugs_16S$amplicon_abund!=0], bugs_mgx$mgx_abund[bugs_mgx$mgx_abund!=0])

# histograms of the distribution of relative abundances of bugs unique to each sample
# important to make sure that we understand distribution before applying any statistical tests :)

par(mfrow=c(2,1))
hist(bugs_16S$amplicon_abund[bugs_16S$amplicon_abund!=0],
     xlab = "relative abundances of bugs uniquely found by amplicon sequencing",
     main = NULL)
hist(bugs_mgx$mgx_abund[bugs_mgx$mgx_abund!=0],
     xlab = "relative abundances of bugs uniquely found by metagenomic sequencing",
     main = NULL)

# plotting absolute difference in relative abundance by taxa

# calculating mean absolute difference
all_taxa <- unique(abund$taxa)
avg_absdiff <- aggregate(abund$abs_diff, by=list(abund$taxa), FUN = mean, na.rm=TRUE)

taxa_order <- avg_diff[order(avg_absdiff$x, decreasing = TRUE),]
top_taxa <- taxa_order[1:5,1] #taxa with largest differences between 16S and mgx
largest_diff_df <- abund[abund$taxa %in% top_taxa,]

p2 <- ggplot(largest_diff_df, aes(x = taxa, y = abs_diff)) + geom_boxplot()
p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("absolute difference in relative abundance")

# plotting total difference in relative abundance by taxa

avg_totdiff <- aggregate(abund$tot_diff, by=list(abund$taxa), FUN = mean, na.rm=TRUE)
taxa_order_totaldiff <- avg_totdiff[order(abs(avg_totdiff$x), decreasing = TRUE),]
top_taxa_totaldiff <- taxa_order_totaldiff[1:5,1] #taxa with largest differences between 16S and mgx
largest_totaldiff_df <- abund[abund$taxa %in% top_taxa_totaldiff,]

p3 <- ggplot(largest_totaldiff_df, aes(x = taxa, y = tot_diff)) + geom_boxplot()
p3 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("total difference in relative abundance")

grid.arrange(p2, p3, ncol = 2)

# calculating times difference was positive or negative

abund$amplicon_greater[abund$tot_diff > 0 ] <- 1 # 16S abundance was greater
abund$mgx_greater[abund$tot_diff < 0 ] <- 1 # mgx abundance was greater

total_amplicon_greater <- aggregate(abund$amplicon_greater, by=list(abund$taxa), FUN = sum, na.rm=TRUE)
total_mgx_greater <- aggregate(abund$mgx_greater, by=list(abund$taxa), FUN = sum, na.rm=TRUE)

ratio_df <- total_amplicon_greater
colnames(ratio_df)[1] <- "taxa"
colnames(ratio_df)[2] <- "amplicon"
ratio_df <- cbind(ratio_df, total_mgx_greater$x)
colnames(ratio_df)[3] <- "mgx"

ratio_df$"ratio: 16S/mgx" <- ratio_df$amplicon/ratio_df$mgx
ratio_df$"ratio: mgx/16S" <- ratio_df$mgx/ratio_df$amplicon

hist(ratio_df$"ratio: 16S/mgx")
hist(ratio_df$"ratio: mgx/16S")

# making tables of the bugs with highest ratio of 16S or mgx

formattable(ratio_df[order(ratio_df$"ratio: 16S/mgx", decreasing = TRUE),][1:15,c(1:4)])
formattable(ratio_df[order(ratio_df$"ratio: mgx/16S", decreasing = TRUE),][1:15,c(1,2,3,5)])


# calculating percentage of bugs present uniquely in mgx, 16S, or both by sample

taxa_diff <- read.csv("taxa_difference.csv")


taxa_diff$amplicon[taxa_diff$amp_avg_abund > 0 ] <- 1
taxa_diff$amplicon[taxa_diff$amp_avg_abund == 0 | is.na(taxa_diff$amp_avg_abund)] <- 0
taxa_diff$mgx[taxa_diff$mgx_avg_abund > 0 ] <- 1
taxa_diff$mgx[taxa_diff$mgx_avg_abund == 0 | is.na(taxa_diff$mgx_avg_abund)] <- 0

taxa_diff <- taxa_diff[order("amplicon"),]

write.csv(taxa_diff,'Documents/thesis/analysis/taxa_difference_bools.csv')

# histogram of abundance difference: 16S-mgx
hist(abund$diff_weight,
     main=NULL,
     xlab= "differences")

# identify bugs in tails

# 



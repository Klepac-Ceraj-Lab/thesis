library(stringr)
library(ggplot2)
library(gridExtra)
library(formattable)
library(phyloseq)
library(plyr)
library(ggpubr)

setwd("/Users/danielle/Documents/thesis/paper-primers")
df <- read.csv("arb-silva.de_testprime_taxlist_929741.csv", header=TRUE, sep=';')

# only keep genera
df$genus <- sapply(strsplit(as.character(df$taxonomy), "\\;"), `[`, 7)
df_genus <- df[!is.na(df$genus),]

# remove Candidates, unknowns

df_genus <- df_genus[!grepl("Candidatus", df_genus$genus),]
df_genus <- df_genus[!grepl("uncultured", df_genus$genus),]

# sort by coverage perc

df_genus <- df_genus[order(df_genus$coverage),]

# find average coverage for genus

avg_coverage <- aggregate(df_genus$coverage,by=list(df_genus$genus),data=df_genus,FUN=mean)

# genera only found by 16S

bugs_mgx <- scan("/Users/danielle/Documents/thesis/paper-abundance-tables/unique_mgx.txt", 
                 what="character", sep=",", 
                 strip.white = TRUE)
bugs_mgx <- gsub("\\[|\\]", "", bugs_mgx)
bugs_mgx_df <- avg_coverage[avg_coverage$Group.1 %in% bugs_mgx,]
not_aligned_bugs <- data.frame("Group.1" = 
                                 c("Agathobaculum", "Ruthenibacterium",
                                    "Abisella", "Lawsonibacter",
                                    "Massilomicrobiota", "Gemmiger",
                                    "Turicimonas", "Metakosakonia",
                                    "Anaeromassilibacillus", "Anaerotignum"),
                               "x" = rep(0,10))
bugs_mgx_df <- rbind(bugs_mgx_df, not_aligned_bugs)
bugs_mgx_df$method <- "mgx"

bugs_amp <- bugs_16S <- scan("/Users/danielle/Documents/thesis/paper-abundance-tables/unique_amplicon.txt", 
                             what="character", sep=",", 
                             strip.white = TRUE)
bugs_amp <- gsub("\\[|\\]", "", bugs_amp)
bugs_amp_df <- avg_coverage[avg_coverage$Group.1 %in% bugs_amp,]
bugs_amp_df$method <- "amp"

bugs_intersection <- scan("/Users/danielle/Documents/thesis/paper-abundance-tables/intersection.txt", 
                          what="character", sep=",", 
                          strip.white = TRUE)
intersection_df <- avg_coverage[avg_coverage$Group.1 %in% bugs_intersection,]
intersection_df$method <- "both"

bugs_method_melt <- rbind(bugs_mgx_df, bugs_amp_df, intersection_df)
# compare coverage of bugs unique to each method

## Kruskal Wallace test

kruskal.test(x ~ method, data = bugs_method_melt)
pairwise.wilcox.test(bugs_method_melt$x, bugs_method_melt$method,
                     p.adjust.method = "BH")

## Violin plot of differences

p1 <- ggplot(bugs_method_melt, aes(x = method, y = x)) + 
  geom_violin()

p1 <- p1 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("primer coverage (%)") + xlab("profiling method") +
  theme(text = element_text(size=20)) +
  stat_compare_means()
p1

# add quartiles and minimum



# average coverage found by mgx
mean(bugs_mgx_df$x)

# look in tree for neighbors that are hit by 16S 
# (mgx is labelling some stuff as one genus, but may be found as something else in silva)
# some bugs that have fine coverage
# others have very low (Pediococcus)
# find taxa not in table bc not in data base (include in table as 0)

# bee swarm plot or violin plot



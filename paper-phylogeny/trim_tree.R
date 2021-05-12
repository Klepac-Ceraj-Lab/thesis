library(ape)

setwd("/Users/danielle/Documents/thesis/paper-phylogeny")

# read in tree with all the species
phylo <- read.tree("genus_union.nwk")
# remove _ in between genus and species
phylo$tip.label <- gsub("_", " ", phylo$tip.label, fixed=TRUE)

# read in list of species we want to remove
remove_tips <- scan("remove_species.txt", what = 'character', sep = ",")
remove_tips <- gsub("\\[|\\]", "", remove_tips)
remove_tips <- trimws(remove_tips)

additional_remove_tips <- c("", 
                            "Streptomycesviridis(Lombardo-Pellegrino1903)Waksman1953")


# drop tips we wish to remove
phylo_drop <- drop.tip(phylo, tip = remove_tips)


# replace species with genera
phylo_drop$tip.label <- unlist((lapply(phylo_drop$tip.label, 
                                      function(x) sub(" .*", "", x))))

phylo_drop <- drop.tip(phylo_drop, tip = additional_remove_tips)

write.tree(phylo_drop, "trimmed_tree_new.nwk")

intersection <- scan("intersection.txt", what = 'character', sep = ",")
unique_mgx <- scan("mgx_only.txt", what = 'character', sep = ",")
unique_amp <- scan("amplicon_only.txt", what = 'character', sep = ",")

# taxa present in metaphlan3 database
present_metaphlan3 <- scan("present_metaphlan3.txt", what = 'character', sep = ",")

tip_labels <- data.frame("taxa" = phylo_drop$tip.label)
tip_labels['method'] <- NA
tip_labels['in_db'] <- NA

tip_labels[tip_labels$taxa %in% intersection,]$method <- "both"
tip_labels[tip_labels$taxa %in% intersection,]$in_db <- 0

tip_labels[tip_labels$taxa %in% unique_mgx,]$method <- "amp"

tip_labels[tip_labels$taxa %in% unique_amp,]$method <- "mgx"
tip_labels[(tip_labels$taxa %in% present_metaphlan3 & tip_labels$taxa %in% unique_amp),]$in_db <- 0
# tip_labels[(!tip_labels$taxa %in% present_metaphlan3 & tip_labels$taxa %in% unique_amp),]$in_db <- 1
# no taxa that were uniquely found by 16S that weren't present in the metaphlan3 database

unique_amp %in% present_metaphlan3

write.csv(tip_labels, file = "tip_labels.csv")



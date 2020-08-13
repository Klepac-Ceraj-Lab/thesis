library(ape)

setwd("/Users/danielle/Documents/thesis/paper-phylogeny")

# replace species with genera

phylo <- read.tree("species_timetree.nwk")
phylo$tip.label <- unlist(lapply(phylo$tip.label, function(x) gsub("_.*","",x)))
write.tree(phylo, "trimmed_tree.nwk")


library("ggtree")
library("ggplot2")
library("phyloseq")
library("treeio")
library("microbiome")

# Set this to your current path
thesis.path <- "/Users/ksb/repos/lab/danielle_thesis"

phylo.path <- file.path(thesis.path, "phylogeny")
taxa.path <- file.path(thesis.path, "analysis/taxa_difference.csv")

setwd(phylo.path)


# creating dataframe for phylogenetic tree metadata
taxa_diff <- read.csv("/Users/danielle/Documents/thesis/analysis/taxa_difference.csv")
tree <- taxa_diff
tree$sequencing[!is.na(taxa_diff$mgx_avg_abund)] <- "mgx"
tree$sequencing[!is.na(taxa_diff$amp_avg_abund)] <- "amp"
tree$sequencing[!is.na(taxa_diff$mgx_avg_abund)&!is.na(taxa_diff$amp_avg_abund) ] <- "both"
tree <- tree[,c(2,7)]
write.csv(tree,"phylogeny.csv")

otu <- as.matrix(read.csv("phyloseq_otu.csv"))
otu[,2:ncol(otu)] <- as.numeric(otu[,2:ncol(otu)])
OTU = otu_table(otu, taxa_are_rows = TRUE)


phy <- read.tree("phyliptree.phy")
metadata <- read.csv("phyloseq_metadata.csv")

physeq1 = merge_phyloseq(otu, metadata)
pseq1 <- read_phyloseq(otu, metadata, type = "simple")

read_csv2phyloseq(
  otu.file = otu,
  metadata.file = metadata,
  sep = ","
)



plot_tree(phy)
plot_tree(physeq1, label.tips="taxa_names", color="sequencing.type") + coord_polar(theta="y")  + 
  geom_point(color="sequencing")


phy <- read.newick("phyliptree.phy")
full_join(phy, metadata, by="label")
ggtree(beast_tree, aes(color=rate))




#test

otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))

OTU = otu_table(otumat, taxa_are_rows = TRUE)

sampledata = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  row.names=rownames(otumat),
  stringsAsFactors=FALSE
))

random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))

physeq <- merge_phyloseq(OTU, sampledata)

plot_tree(physeq, color="sequencing.type", label.tips="taxa_names", ladderize="left", plot.margin=0.3)

# playing with Global Patterns

physeq = prune_taxa(taxa_names(GlobalPatterns)[1:50], GlobalPatterns)
plot_tree(physeq)
physeq <- merge_phyloseq(GlobalPatterns@otu_table, sampledata, random_tree)
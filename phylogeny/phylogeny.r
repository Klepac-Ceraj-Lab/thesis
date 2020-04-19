library("ggtree")
library("ggplot2")
library("phyloseq")
library("treeio")
library("microbiome")

# Set this to your current path
thesis.path <- "/Users/danielle/Documents/thesis/"

phylo.path <- file.path(thesis.path, "phylogeny")
taxa.path <- file.path(thesis.path, "analysis/taxa_difference.csv")

setwd(phylo.path)


# creating dataframe for phylogenetic tree metadata
taxa_diff <- read.csv(taxa.path)
tree <- taxa_diff
tree$sequencing[!is.na(taxa_diff$mgx_avg_abund)] <- "mgx"
tree$sequencing[!is.na(taxa_diff$amp_avg_abund)] <- "amp"
tree$sequencing[!is.na(taxa_diff$mgx_avg_abund)&!is.na(taxa_diff$amp_avg_abund) ] <- "both"
tree <- tree[,c(2,7)]
write.csv(tree,"phylogeny.csv")

otu <- as.data.frame(read.csv("phyloseq_otu.csv"))
# convert NA to 0
otu[is.na(otu)] <- 0
# find smallest non-zero number
min.non0 <- min(otu[,-1][otu[,-1] > 0])

# convert everything to ints
otu[,-1] <- sapply(otu[,-1], function(x) {x * (1/min.non0)})
otu.matrix = as.matrix(otu[-1])

sd <- sample_data(otu)

OTU <- otu_table(otu.matrix, taxa_are_rows=TRUE)

phy <- read.tree("phyliptree.phy")
metadata <- read.csv("phyloseq_metadata.csv")

physeq1 = merge_phyloseq(OTU, metadata, phy)

plot_tree(phy)
plot_tree(physeq1, label.tips="taxa_names", color="sequencing.type") + coord_polar(theta="y")  + 
  geom_point(color="sequencing")





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

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq = phyloseq(OTU, TAX)

physeq1 = merge_phyloseq(OTU, sampledata, random_tree)

plot_tree(physeq, color="sequencing.type", label.tips="taxa_names", ladderize="left", plot.margin=0.3)


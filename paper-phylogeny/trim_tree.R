library(ape)

setwd("/Users/danielle/Documents/thesis/paper-phylogeny")

# read in original phylogenetic tree
phylo <- read.tree("all_taxa.nwk")

# read in tip names (species names)
tip <- c(scan("remove_species.txt", what="character", sep=",", 
                            strip.white = TRUE))

# replace "_" with empty spaces in species names (like how they are in phylo tree)
tip <- lapply(tip, function(x) gsub(" ", "_", x))

plot(drop.tip(phylo, tip, trim.internal = FALSE))
write.tree((drop.tip(phylo, tip, trim.internal = FALSE)), "trimmed_tree.nwk")

# example
data(bird.families)
plot(bird.families)
example <- c(
  "Eopsaltriidae", "Acanthisittidae", "Pittidae", "Eurylaimidae",
  "Philepittidae", "Tyrannidae", "Thamnophilidae", "Furnariidae",
  "Formicariidae", "Conopophagidae", "Rhinocryptidae", "Climacteridae",
  "Menuridae", "Ptilonorhynchidae", "Maluridae", "Meliphagidae",
  "Pardalotidae", "Petroicidae", "Irenidae", "Orthonychidae",
  "Pomatostomidae", "Laniidae", "Vireonidae", "Corvidae",
  "Callaeatidae", "Picathartidae", "Bombycillidae", "Cinclidae",
  "Muscicapidae", "Sturnidae", "Sittidae", "Certhiidae",
  "Paridae", "Aegithalidae", "Hirundinidae", "Regulidae",
  "Pycnonotidae", "Hypocoliidae", "Cisticolidae", "Zosteropidae",
  "Sylviidae", "Alaudidae", "Nectariniidae", "Melanocharitidae",
  "Paramythiidae","Passeridae", "Fringillidae")
plot(drop.tip(bird.families, tip))
plot(drop.tip(bird.families, tip, trim.internal = FALSE))
data(bird.orders)
plot(drop.tip(bird.orders, 6:23, subtree = TRUE))
plot(drop.tip(bird.orders, c(1:5, 20:23), subtree = TRUE))
plot(drop.tip(bird.orders, c(1:20, 23), subtree = TRUE))
plot(drop.tip(bird.orders, c(1:20, 23), subtree = TRUE, rooted = FALSE))
### Examples of the use of `root.edge'
tr <- read.tree(text = "(A:1,(B:1,(C:1,(D:1,E:1):1):1):1):1;")
drop.tip(tr, c("A", "B"), root.edge = 0) # = (C:1,(D:1,E:1):1);
drop.tip(tr, c("A", "B"), root.edge = 1) # = (C:1,(D:1,E:1):1):1;
drop.tip(tr, c("A", "B"), root.edge = 2) # = (C:1,(D:1,E:1):1):2;
drop.tip(tr, c("A", "B"), root.edge = 3) # = (C:1,(D:1,E:1):1):3;
# }

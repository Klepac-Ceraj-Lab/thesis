library("ggtree")
library("ggplot2")
library("phyloseq")
library("treeio")
library("microbiome")
library("vegan")
library("ggdendro")


# Set this to your current path
thesis.path <- "/Users/danielle/Documents/thesis/"

phylo.path <- file.path(thesis.path, "phylogeny")
taxa.path <- file.path(thesis.path, "analysis/taxa_difference.csv")

setwd(phylo.path)


# code that may work for making phylogenetic tree!

otu <- read.csv("/Users/danielle/Documents/thesis/analysis/transposed_mgxamp_df.csv")
otu[is.na(otu)] <- 0

meta <- otu[1:4]

betad<-vegdist(otu[,5:306] ,method="bray")
hc <- hclust(betad)
hc_d <- dendro_data(as.dendrogram(hc))
hc_d$labels$Type<-meta[as.character(hc_d$labels$label),4]
gg_color_hue<-function(n){
  hues=seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}

cols=gg_color_hue(length(unique(hc_d$labels$Type)))
hc_d$labels$color=cols[hc_d$labels$Type]

## Plot clusters
p1 <- ggplot(data = segment(hc_d)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_x_discrete(labels=label(hc_d)$label) +
  ylab("Distance (beta diversity = bray)") + theme_bw()+
  theme(axis.text.y = element_text(color = hc_d$labels$color),
        axis.title.y = element_blank())
p1 <- p1 + geom_point(data=hc_d$label, aes(x = x, y = y, color = Type), inherit.aes =F, alpha = 0)
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=1, alpha = 10)))+
  scale_color_manual(values = cols)
print(p1)



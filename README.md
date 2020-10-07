# Danielle Peterson Thesis

## Dependencies

```R
install.packages("ggplot2") # 3.3.2
install.packages("ape") # 5.4
install.packages("gridExtra") # 2.3
install.packages("formattable") # 0.2.0.1
install.packages("phyloseq") # 1.30.0
install.packages("plyr") # 1.8.6
install.packages("dplyr") # 1.0.0
install.packages("ggpubr") # 0.4.0
install.packages("vegan") # 2.5.6
install.packages("sna") # 2.5
install.packages("reshape2") # 1.4.4
install.packages("gtools") # 3.8.2
install.packages("sparseDOSSA") # 1.10.0
# ... others
```

### Bioconductor

Several packages require bioconductor # 3.10

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

```R
library(BiocManager)
BiocManager::install("ggtree") # 2.0.4
BiocManager::install("treeio") # 1.10.0
BiocManager::install("microbiome") # 1.8.0
BiocManager::install("phyloseq") # 1.30.0
# ... others
```


## FASTQ Subsampling code

- see https://github.com/kescobo/dpthesis-subsample
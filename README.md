# Danielle Peterson Thesis

## Contents

- include order for going through notebooks
  - Alternatively (or in addition) you can number directories (eg `1_16S`)

## Dependencies

```R
install.packages("ggplot2") # add version as comment
# ... others
```

### Bioconductor

Several packages require bioconductor

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

```R
library(BiocManager)
BiocManager::install("ggtree") # add version as comment
BiocManager::install("treeio") # add version as comment
BiocManager::install("microbiome") # add version as comment
BiocManager::install("phyloseq") # add version as comment
# ... others
```


## FASTQ Subsampling code

- see https://github.com/kescobo/dpthesis-subsample
  - commit hash used
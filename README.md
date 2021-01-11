# Danielle Peterson Thesis

## Dependencies

The R environment is maintained by [`renv`](https://blog.rstudio.com/2019/11/06/renv-project-environments-for-r/).

```R
install.packages("renv")
```

```R
library(renv)
renv::init()
```

Note - if this is the first time using `renv`, this might take a while.

Project dependencies and versions can be found in `renv.lock`
and are automatically loaded / installed when you run `renv::init()`

They include:

- ggplot2
- ape
- gridExtra
- formattable
- phyloseq
- plyr
- dplyr
- ggpubr
- vegan
- sna
- reshape2
- gtools
```

### Bioconductor

Several packages require bioconductor # 3.10

This is how they were installed initially,
but renv should handle this in the future.

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
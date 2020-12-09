# GSEABenchmarkeR

The [GSEABenchmarkeR](https://bioconductor.org/packages/GSEABenchmarkeR)
package implements an extendable framework for 
reproducible evaluation of set- and network-based methods for enrichment 
analysis of gene expression data. This includes support for the efficient 
execution of these methods on comprehensive real data compendia (microarray and 
RNA-seq) using parallel computation on standard workstations and 
institutional computer grids. Methods can then be assessed with respect to 
runtime, statistical significance, and relevance of the results for the 
phenotypes investigated.
    
## Installation

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GSEABenchmarkeR")
```


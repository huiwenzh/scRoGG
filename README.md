# scRoGG: identifying RObust Gene-Gene correlation for single-cell RNA-sequencing data (scRoGG)

## Installation                                            <img src="https://user-images.githubusercontent.com/46465953/211182218-b577f94b-6b44-4c11-83aa-7d2ba20406e2.png" width=10% height=12%> 

For installation please use the following codes in R
```
devtools::install_github("huiwenzh/scRoGG")
library(scRoGG)
# requires the installastion of foreach, igraph, fgsea, msigdb, ggplot2 and propr
```
We welcome any suggestions on the package! Please submit a ticket to issues if you have encountered any techinical problems. For more details for the method, please refer to our manuscript. 

## 1. Background

Single-cell RNA-sequencing data is highly sparse and heterogenous, which largely impact the measurement of coexpression. Traditional methods such as Pearson and Spearman correlation fail to capture the gene pair correlation and often return low correlation coefficient even with significantly correlated gene pairs. As suggested in the review paper, proportional correlation showed best performance in measuring the assoiation in single-cell RNA-seq data [1]. Further extensive discussion suggested the improvement of proportional correlation can be achieved by providing a set of reference [2]. From our previous work, scRoGG implements essential genes that hve been tested in single-cell atlas (scEssentials) that shown high expression, detection rates and tight coexpression level.

## 2. Main functions
The input gene count matrix can be either raw or normalised but not log data, which can be switched on by parameter **normalised**. As the assumption of the asymptotic test we implemented, we suggest having at least 30 scEssentials identified for each test. Consider increasing **ES_number** to obtain a sufficient number for the downstream analysis. 

scRoGG has three main functions, which are:
1. Identify robust coexpressed genes in a single condition. Here we use **naive T cell** from healthy PMBC data that contains 711 cells and ~20k genes [3].  

```
library(scRoGG)
library(Seurat)
# extract raw count for naive CD4+ T cells
naiveT <- subset(pbmc, subset=name =='naive CD T' )
naiveT <- as.matrix(naiveT@assays$RNA@counts)
# run scRoGG with default settings
naiveT_cor <- scRoGG(naiveT) 
naiveT_sig_cor <- robustness(naiveT_cor) # 1884 signficantly robustly coexpressed genes
# with extremely large dataset, one can turn stats = F so that top 0.1% quantile of abs(RS) gene pairs will be returned.
naiveT_sig_cor1 <- robustness(naiveT_cor,stats = T)
```
 







## Reference
[1] Skinnider, M.A., Squair, J.W. & Foster, L.J. Evaluating measures of association for single-cell transcriptomics. Nat Methods 16, 381–386 (2019). https://doi.org/10.1038/s41592-019-0372-4.
[2] Erb, I., Notredame, C. How should we measure proportionality on relative gene expression data?. Theory Biosci. 135, 21–36 (2016). https://doi.org/10.1007/s12064-015-0220-8.
[3] Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017). https://doi.org/10.1038/ncomms14049.



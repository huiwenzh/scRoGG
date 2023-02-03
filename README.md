# scRoGG- identifying RObust Gene-Gene correlation for single-cell RNA-sequencing data (scRoGG)

## Installation                                            <img src="https://user-images.githubusercontent.com/46465953/211182218-b577f94b-6b44-4c11-83aa-7d2ba20406e2.png" width=15% height=12%> 

```
devtools::install_github("huiwenzh/scRoGG")
library(scRoGG)
# requires the installastion of foreach, igraph, fgsea and propr
```
We welcome any suggestions on the package! Please submit an issue  
## 1. Background

Single-cell RNA-sequencing data is highly sparse and heterogenous, which largely impact the measurement of coexpression. Traditional methods such as Pearson and Spearman correlation fail to capture the gene pair correlation and often return low correlation coefficient even with significantly correlated gene pairs. Inspired by [review paper](https://www.nature.com/articles/s41592-019-0372-4), proportional correlation showed best performance in measuring the assoiation in single-cell RNA-seq data. Further extensive discussion suggested the improvement of proportional correlation can be achieved by providing a set of reference [1]. From our previous work, scRoGG implements essential genes that hve been tested in single-cell atlas (scEssentials) that shown high expression, detection rates and tight coexpression level.









### Reference
[1] Erb, I., Notredame, C. How should we measure proportionality on relative gene expression data?. Theory Biosci. 135, 21–36 (2016). https://doi.org/10.1007/s12064-015-0220-8.


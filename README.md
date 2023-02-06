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

scRoGG has **four** main functions, which are:
1. Identify robust coexpressed genes in a single condition. Here I use **naive T cell** from healthy PMBC data that contains 711 cells and ~20k genes [3].

```
library(scRoGG)
library(Seurat)
# extract raw count for naive CD4+ T cells
naiveT <- readRDS('~ Example_data/naiveT.rds')
naiveT <- as.matrix(naiveT@assays$RNA@counts)
# run scRoGG with default settings
naiveT_cor <- scRoGG(naiveT) 
naiveT_sig_cor <- robustness(naiveT_cor) # 1884 signficantly robustly coexpressed genes
# with extremely large dataset, one can turn stats = F so that top 0.1% quantile of abs(RS) gene pairs will be returned.
naiveT_sig_cor1 <- robustness(naiveT_cor,stats = T)

# To visualise a specific gene pair

plotCoExp(naiveT_cor$transformed_data,'CCL7','NKG7')
```
2. Identify dofferentailly coexpressed gene pairs between multiple conditions. Here we use three endothelial cells from **mouse brain arteriovenous region**, including arterial endothelial (aEC), capillary endothelial (cEC) and venous endothelial (vEC). Each cell type contains 397, 405, and 298 cells [4].
```
library(scRoGG)
aEC <- brain_3EC[,brain_3EC$annoated.cell.types=='aEC']
#Smartseq2 data is less sparse than 10X genomics data, so filter has increased to 0.2 to reduce the computational burden. For mouse data, org sets to 'mmu'. 
aEC_cor  <- scRoGG(dat=assay(aEC),filter = 0.2,ES_number = 200,org = 'mmu') 

capEC <- brain_3EC[,brain_3EC$annoated.cell.types=='capilEC']
capEC_cor <- scRoGG(dat=assay(capEC),filter = 0.2,ES_number = 200,org = 'mmu') 

vEC <- brain_3EC[,brain_3EC$annoated.cell.types=='vEC']
vEC_cor <- scRoGG(dat=assay(vEC),filter = 0.2,ES_number = 200,org = 'mmu')

# To test DC genes, robustness2 is applied
EC_list <- list(aEC_cor,capEC_cor,vEC_cor)
ECs_cordiff <- robustness2(EC_list,p.adj = 0.1) #221 significantly DC genes identified
```
3. Constructing coexpression network with community detection by Ledien clustering, returing *igraph* object for visualisation transformation. This method can be used for any types of edgelist object (3 columns: gene1; gene2; association). Here, I use the result from naiveT cells as an example. For more details, please refere to our vignette.

```
library(scRoGG)
naiveT_net <- coExp_network(naiveT_sig_cor[,1:3],n_networks = 1) # To plot the largest sub-communities
# naiveT_net is a list with 'all' for the full network and sub_1 as the largest sub-communities by size.
```

4. Performing network-based GSEA (net-GSEA). In scRoGG, we implement a novel method to perform gene-set enrichment analysis for network object (igraph). To obtain gene set information, scRoGG implements the syntax from msigd. Please refer to our manuscript for how the analysis is performed.

```
naiveT_gsea <- net_gsea(network = naiveT_net[['all']],species = 'Homo sapiens', category = 'C5',subcategory = 'GO:BP',minSize = 20,sign = F) # against GO BP pathways
naiveT_gsea1 <- net_gsea(network = naiveT_net[['all']],species = 'Homo sapiens', category = 'C2',subcategory = 'CP:KEGG',minSize = 10) # against KEGG pathways

# you can visualise the significantly entiched pathways in barplot.
library(ggplot2)
plotPways(naiveT_gsea)+ggtitle('Naive CD4+ T enriched pathways')
```


## Reference
[1] Skinnider, M.A., Squair, J.W. & Foster, L.J. Evaluating measures of association for single-cell transcriptomics. Nat Methods 16, 381–386 (2019). https://doi.org/10.1038/s41592-019-0372-4.

[2] Erb, I., Notredame, C. How should we measure proportionality on relative gene expression data?. Theory Biosci. 135, 21–36 (2016). https://doi.org/10.1007/s12064-015-0220-8.

[3] Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017). https://doi.org/10.1038/ncomms14049.

[4] Vanlandewijck, M., He, L., Mäe, M. A., Andrae, J., Ando, K., Del Gaudio, F., Nahar, K., Lebouvier, T., Laviña, B., Gouveia, L., Sun, Y., Raschperger, E., Räsänen, M., Zarb, Y., Mochizuki, N., Keller, A., Lendahl, U., & Betsholtz, C. (2018). A molecular atlas of cell types and zonation in the brain vasculature. Nature, 554(7693), 475–480. https://doi.org/10.1038/nature25739.




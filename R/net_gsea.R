#' Perform network-based geneset enrichment analysis
#' @author Huiwen Zheng
#' @param species Species information to msigdbr package.
#' @param category Category information to msigdbr package.
#' @param subcategory Subcategory information to msigdbr package.
#' @param minSize Minimum size in each pathway.
#' @param maxSize Maximum size in each pathway.
#' @param BPPARAM Parellalisation
#' @usage net_gsea(G, species='Homo sapiens',category='H', minSize=5)
#' @export

net_gsea <- function(network,species = "Homo sapiens", category , subcategory = NULL, minSize=5,sign=F,maxSize=500,BPPARAM=NULL){
  msigdbr_df <- msigdbr::msigdbr(species = species, category = category,subcategory = subcategory )
  pathways <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
  # calculate average weight
  if(sign==F){
    igraph::E(network)$weight <- abs(igraph::E(network)$weight)
    satistic <- cbind(igraph::as_edgelist(network),igraph::E(network)$weight)
    RCC_ave <- c()
    for (ver in union(satistic[,1],satistic[,2])){
      val  <-  mean(as.numeric(satistic[satistic[,1]==ver|satistic[,2]==ver,3]))
      names(val) <- ver
      RCC_ave <- c(RCC_ave,val)
    }
    stat <- RCC_ave + igraph::hub_score(network)$vector
  }
  if(sign==T){
    satistic <- cbind(igraph::as_edgelist(network),igraph::E(network)$weight)
    RCC_ave <- c()
    for (ver in union(satistic[,1],satistic[,2])){
      val  <-  mean(as.numeric(satistic[satistic[,1]==ver|satistic[,2]==ver,3]))
      names(val) <- ver
      RCC_ave <- c(RCC_ave,val)
    }
    stat <- RCC_ave + igraph::hub_score(network)$vector
  }

result.GSEA <- fgsea::fgsea(pathways, stat,  minSize = minSize,scoreType="pos", maxSize = maxSize, BPPARAM = BPPARAM)
result.GSEA
}

# consider scoreType="pos" to 'std' to group-wise comparison

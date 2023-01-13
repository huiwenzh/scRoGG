#' Visualise pathway enrichment results
#' @author Huiwen Zheng
#' @param pways_df Output from pathway analysis.
#' @param n Number of pathways plotted in the barplot.
#' @param p.adj Adjusted p value to select significantly enriched genes
#' @export
plotPways <- function(pways_df, n=10, p.adj = 0.1){
  pways_df <- pways_df[order(pways_df$padj),]
  if (nrow(pways_df)>=n){n=n}
  if (nrow(pways_df)<n){n=nrow(pways_df)}
  pways_sub <- pways_df[1:n,]
  pways_sub <- pways_sub[pways_sub$padj<p.adj,]
  #pways_sub$size <- as.numeric(pways_sub$size)
  pways_sub <- pways_sub[order(pways_sub$size, decreasing = F),]
  ggplot2::ggplot(pways_sub,aes(reorder(pathway, size), y= size, fill=padj))+
    geom_col()+
    scale_x_discrete(labels = function(x)
      stringr::str_wrap(x, width = 16))+
    scale_fill_gradient2(low = "yellow",mid = 'yellow', high = "red")+coord_flip()+xlab('')+theme_classic()
}

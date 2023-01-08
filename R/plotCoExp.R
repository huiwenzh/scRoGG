#' Visualise correlated gene pairs with transformed data
#' @author Huiwen Zheng
#' @param dat Transformed data, encourage to use the transformed data from prop_distribution function but can be applied to any data.
#' @param gene1 A gene to plot
#' @param gene2 A gene to plot
#' @export

plotCoExp <- function(dat, gene1, gene2){
df <- as.data.frame(cbind(dat[,gene1],dat[,gene2]))
colnames(df) <- c(gene1, gene2)
p <- ggpubr::ggscatter(df, x = gene1, y = gene2,
                       size = 1.3, alpha = 0.6)
ggExtra::ggMarginal(p, type = "density",xparams = list(fill = "darkgreen"),yparams = list(fill = "orange"))
}

#'

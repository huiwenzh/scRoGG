#' Visualisation of selected coexpressed gene pairs' correlation distribution
#' @author Huiwen Zheng
#' @param prop_list A list of prop_distribution dataframes.
#' @param gene_name1 Gene1 to plot
#' @param gene_name2 Gene2 to plot
#' @return A ggplot object with violin plot
#' @export
#' 
plotCoDis <- function(dat_list,gene_name1, gene_name2){
if (length(dat_list)==1){
  dat <- dat_list[[1]]
  dat2 <- dat[dat$gene1==gene_name1&dat$gene2==gene_name2,-1:-3]
  null_dist <- colMeans(dat[,-1:-3], na.rm = T)
  p_data <- data.frame(value=c(as.numeric(na.omit(dat2)),null_dist), group = rep(c("Sample",'Null'),c(length(na.omit(dat2)),length(null_dist))))
  ggplot(p_data, aes(x=group, y=value, fill=group)) + 
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, fill="white")+
    labs(title=paste0("Correlation distribution (",gene_name1,' , ', gene_name2,")"),x="", y = "Correlation under scEssential")+ scale_fill_brewer(palette="Dark2") + theme(legend.position = "none")+ theme_bw()+ theme(legend.position = "none")      
}
else{
  dat <- sapply(dat_list,function(l)as.numeric(na.omit(l[l$gene1==gene_name1&l$gene2==gene_name2,-1:-3]))*l[l$gene1==gene_name1&l$gene2==gene_name2,3])
  p_data <- data.frame(value=unlist(dat), group = rep(paste0('sample',1:length(dat_list)),unlist(sapply(dat_list,function(s)ncol(s)-3))))
  
  ggplot(p_data, aes(x=group, y=value, fill=group)) + 
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.2, fill="white")+
    labs(title=paste0("Correlation distribution (",gene_name1,' , ', gene_name2,")"),x="", y = "Correlation under scEssential")+ scale_fill_brewer(palette="Dark2") + theme(legend.position = "none")+ theme_bw()+ theme(legend.position = "none")      
}
}


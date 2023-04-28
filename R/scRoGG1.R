#'Generate proportional distribution based on essential genes for identifying robust gene-gene correlation in scRNA-seq data (with updated scEssentials)
#' @author Huiwen Zheng
#' @param dat A cell*gene matrix that can be raw or normliased without logarithm.
#' @param normalised By default, performs normalisation via scattle package.     If changes to FALSE, the input data will be used in the following steps.
#' @param filter By default, genes that expressed in less than 0.1 of total cells will be removed. You may want to increase the filtering threshold with more in-depth sequencing platforms.
#' @param ES_number Number of the ES genes used in estimating proportional correlation distribution.
#' @param org Selects certain ES genes under the organism. Currently only hsa and mmu are available.
#' @param bino Bivariate distribution probability.
#' @param details By default, transformed data will be returned and can be used to visualise  correlation.
#' @return A list with three components:      1. bi_zero: joint zero for each gene pair.     2. prop_distribution: correlation distribution for each gene pair.    3.  transformed_data: transformed data using the mean of essential genes.
#' @export
#'
scRoGG1 <- function(dat,normalised = TRUE, filter = 0.1, ES_number = 1000, org = 'hsa', bino = 0.05,details = T){
  if (normalised==T) {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=dat))
    sce <- scuttle::logNormCounts(sce[,colSums(SingleCellExperiment::counts(sce)) > 0], log=F)
    dat_keep <- SingleCellExperiment::normcounts(sce)
  }
  index <-  apply(dat_keep,1,function(y)sum(y>0))
  # at least present in more than 10%
  dat_keep <- dat_keep[index>(filter*ncol(dat_keep)),]
  # remove specific gene names - mito
  MT <- grep(pattern = "^MT-|^mt-", x = rownames(dat_keep), value = TRUE)
  dat_keep <- dat_keep[setdiff(rownames(dat_keep),MT),]
  message(paste('Removing', length(MT) , 'mitochondrial genes'))
  dat.coprob = dismay::binomial(t(dat_keep))
  # less than 0.05
  dat.coprob[lower.tri(dat.coprob,diag=TRUE)] <- 0
  cor_df <- reshape2::melt(dat.coprob)
  cor_df_sig <- cor_df[cor_df$value>(-log10(bino)),]

  dat_p <- dat_keep[rownames(dat_keep)%in%union(cor_df_sig$Var1,cor_df_sig$Var2),]

  message(paste0(nrow(cor_df_sig),' gene pairs pass the bivariate distribution probability test'))

  # calculate proportionality against per essential genes
  ct <- t(dat_p)

  # get each essential genes
  if (org=='hsa'){
    ES <- scEssential_hsa
  }
  else if(org=='mmu'){
    ES <- scEssential_mmu
  }
  else {warning('No data exist! Please try keyword: hsa or mmu')}
  ES <- ES[order(ES$score, decreasing = T),]
  with_ES <- ES$gene[1:ES_number][ES$gene[1:ES_number] %in% colnames(ct)]
  with_ES <- colnames(ct[,colnames(ct) %in% with_ES])[colSums(ct[,colnames(ct) %in% with_ES]>0)/nrow(ct)>0.5]
  message(paste0('Using ',length(with_ES),' essential genes to estimate proportional distribution'))
  if (length(with_ES)<30){
    warning('Warining! Number of essential genes used to estimate proportional distribution is less than 30, consider to increase ES_number!')
  }

  #replace zeros with next smallest value
  if(any(as.matrix(ct) == 0)){
    zeros <- ct == 0
    ct[zeros] <- min(ct[!zeros])
  }
  x <- list()
  if (details==T){
    use <- propr::ivar2index(ct, ivar = with_ES)
    logX <- log(ct)
    logSet <- logX[, use, drop = FALSE]
    ref <- rowMeans(logSet)
    lr <- sweep(logX, 1, ref, "-")
    x[['transformed_data']] <- lr
  }
  # remove ribosomal genes
  cor_df_sig <- cor_df_sig[!is.infinite(cor_df_sig$value),]
  cor_phs <- lapply(with_ES, function(a){
    use <- propr::ivar2index(ct, ivar = a)
    logX <- log(ct)
    logSet <- logX[, use, drop = FALSE]
    ref <- rowMeans(logSet)
    lr <- sweep(logX, 1, ref, "-")
    phs_cor <- apply(cor_df_sig,1,function(s){
      if(s[1]==a|s[2]==a) NA
      else
        phs <-  var(lr[,s[1]]-lr[,s[2]])/ var(lr[,s[1]]+lr[,s[2]])
    })
    (max(phs_cor, na.rm = T)-phs_cor)/max(phs_cor, na.rm = T)

  })
  pearson_sign <- sign(apply(cor_df_sig[,1:2],1,function(b){cor(dat_p[b[1],],dat_p[b[2],])}))

  cors <- cbind(cor_df_sig[,1:2],pearson_sign,cor_phs)
  colnames(cors)[1:3] <- c('gene1','gene2','sign')
  x[["prop_distribution"]] <- cors
  bi_nonzero <- apply(cor_df_sig[,1:2],1,function(b){
    datt <- sum((dat_p[rownames(dat_p)==b[1],])*(dat_p[rownames(dat_p)==b[2],])==0)
    datt_pct <- (ncol(dat_p)-datt)/ncol(dat_p)
  })

  x[['bi_nonzero']] <- bi_nonzero
  x
}

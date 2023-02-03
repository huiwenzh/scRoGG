# figure plotting 
# Fig 1 #####
# find one example 
dat_cv <- data.frame(value=c(as.numeric(na.omit(a)),null_dist), group = rep(c("sample1",'sample2'),c(length(na.omit(a)),length(null_dist))))
library(ggthemr)
dust_theme <- ggthemr(palette = "dust", set_theme = FALSE,spacing = 1, type = 'outer')

ggplot(dat_cv, aes(x=group, y=value, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="                          Correlation CV (𝑔𝑖, 𝑔𝑗)",x="", y = "Correlation")+ scale_fill_brewer(palette="Dark2") +
  dust_theme$theme  + theme(legend.position = "none")

ggthemr_reset()

# Fig 2 ######
# number of sig. correlated gene pairs and cell number 
pmbc_df <- data.frame(value = c(nrow(naiveT_sig_cor),nrow(CD14_Mono_sig_cor),nrow(MemoryT_sig_cor), nrow(B_sig_cor), nrow(CD8T_sig_cor), nrow(FCGR3A_sig_cor),nrow(NK_sig_cor)),
                      type = c('Naive CD4+ T','CD14+ Mono','Memory CD4+ T','B','CD8+ T','FCRG3A+ Mono','NK'),group = c('Adaptive','Innate','Adaptive','Adaptive','Adaptive','Innate','Innate'), 
                      number= c(711,480,472,344,279,162,144))

pmbc_df$type <- factor(pmbc_df$type,  
                  levels = pmbc_df$type[order(pmbc_df$value, decreasing = TRUE)])

fig2a <- ggplot(pmbc_df, aes(x=type, y=value, fill=group))+
  geom_bar(stat="identity", color="black")+
  scale_fill_brewer(palette="Dark2")+
  geom_text(aes(label=value), vjust=-1, size=3.5, col= 'black') +
  theme_classic()+labs(x="", y = "Number of significantly correlated gene pairs")
wilcox.test(c(1884, 10288,154,14608),c(69003,57566,16158))
 
library(ggpubr)
ggboxplot(pmbc_df, x = "group", y = "value",
          fill = "group", palette = "Dark2", add = "jitter")+ stat_compare_means(label.y = 75000,)+
  labs(x="", y = "Number of significantly correlated gene pairs")+theme_bw()+ theme(legend.position = "none")

cor_list <- list("Naive CD4+ T" = paste0(naiveT_sig_cor$gene1,'-',naiveT_sig_cor$gene2),"CD14+ Mono" = paste0(CD14_Mono_sig_cor$gene1,'-',CD14_Mono_sig_cor$gene2),"Memory CD4+ T" = paste0(MemoryT_sig_cor$gene1,'-',MemoryT_sig_cor$gene2),'B' = paste0(B_sig_cor$gene1,'-',B_sig_cor$gene2),'CD8+ T' = paste0(CD8T_sig_cor$gene1,'-',CD8T_sig_cor$gene2),'FCRG3A+ Mono'= paste0(FCGR3A_sig_cor$gene1,'-',FCGR3A_sig_cor$gene2),'NK'= paste0(NK_sig_cor$gene1,'-',NK_sig_cor$gene2))
Reduce(intersection, cor_list)
library(UpSetR)
cor_share <- fromList(cor_list)
UpSetR::upset(cor_share,nsets = 7)
library(ggVennDiagram)
ggVennDiagram(
  cor_list,label_alpha = 0, label = 'count')+ scale_fill_gradient(low="white",high = "red")
library(dendextend)
dend <- t(cor_share) %>% 
  dist(method = 'binary') %>% # calculate a distance matrix, 
  hclust(method = "centroid") %>% # on it compute hierarchical clustering using the "average" method, 
  as.dendrogram # and lastly, turn that object into a dendrogram.

dend %>% plot(ylim=c(0,1))

# top 10 cor from each cell type 
# for supplementary figure
library("dplyr") 
topcor <- function(x, n=10){
  data_new <- x %>%              
    arrange(desc(abs(RS))) %>% 
    slice(1:n)
  as.data.frame(data_new)
}
naiveT_10 <- topcor(naiveT_sig_cor)
CD14_Mono_10 <- topcor(CD14_Mono_sig_cor)
MemoryT_10 <- topcor(MemoryT_sig_cor)
B_10 <- topcor(B_sig_cor)
CD8T_10 <- topcor(CD8T_sig_cor)
FCRG3A_Mono_10 <- topcor(FCGR3A_sig_cor)
NK_10 <- topcor(NK_sig_cor)

library(scRoGG)
p1 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[1]),as.character(FCRG3A_Mono_10$gene2[1]))
p2 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[2]),as.character(FCRG3A_Mono_10$gene2[2]))
p3 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[3]),as.character(FCRG3A_Mono_10$gene2[3]))
p4 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[4]),as.character(FCRG3A_Mono_10$gene2[4]))
p5 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[5]),as.character(FCRG3A_Mono_10$gene2[5]))
p6 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[6]),as.character(FCRG3A_Mono_10$gene2[6]))
p7 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[7]),as.character(FCRG3A_Mono_10$gene2[7]))
p8 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[8]),as.character(FCRG3A_Mono_10$gene2[8]))
p9 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[9]),as.character(FCRG3A_Mono_10$gene2[9]))
p10 <- plotCoExp(FCGR3A_Mono$transformed_data, as.character(FCRG3A_Mono_10$gene1[10]),as.character(FCRG3A_Mono_10$gene2[10]))
P <- ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol = 5,nrow = 2,labels = 'auto') 
annotate_figure(P, top = text_grob("FCGR3A+ Mono", face = "bold", size = 14))

# selected gene pair for plotting Fig2b
# specific gene pair for CD4+ T and memeory CD4+ T
fig2b1 <- plotCoExp(naiveT$transformed_data,'NKG7','CCL5')+ggtitle('Naive CD4+ T')
fig2b2 <- plotCoExp(MemoryT$transformed_data,'GZMA','GZMK')
# shared gene pair COX6A1 and COX6A1P2
fig2c1 <- plotCoExp(naiveT$transformed_data,'COX6A1','COX6A1P2')
fig2c2 <- plotCoExp(CD8T$transformed_data,'COX6A1','COX6A1P2')
fig2bc <- ggpubr::ggarrange(fig2b1,fig2b2,fig2c1,fig2c1, ncol = 2, nrow = 2, labels =c("b","c","d","e"))
ggpubr::ggarrange(fig2a,fig2bc, ncol=2,nrow = 1,labels ='g')

# supp. to show the biggest subcommunity
# restricted to top 1000, if the number of correlated gene pairs is greater than 20,000.
naiveT_net <- coExp_network(naiveT_sig_cor[,1:3],n_networks = 1)
CD14_Mono_net <- coExp_network(CD14_Mono_sig_cor[1:1000,1:3],n_networks = 1)
MemoryT_net <- coExp_network(MemoryT_sig_cor[,1:3],n_networks = 1)
B_net <- coExp_network(B_sig_cor[,1:3],n_networks = 1)
CD8T_net <- coExp_network(CD8T_sig_cor[1:1000,1:3],n_networks = 1)
FCGR3A_net <- coExp_network(FCGR3A_sig_cor[1:1000,1:3],n_networks = 1)
NK_net <- coExp_network(NK_sig_cor[1:1000,1:3],n_networks = 1)
gd <- edge_density(MemoryT_net[['sub_1']])

g.random <- erdos.renyi.game(n = gorder(MemoryT_net[['sub_1']]), p.or.m = gd, type = "gnp")
plot(g.random)
plot(MemoryT_net[['sub_1']])
# net-based GSEA
MemoryT_gsea <- net_gsea(network = MemoryT_net[['all']],species = 'Homo sapiens', category = 'C5',subcategory = 'GO:BP',minSize = 20)
plotPways(MemoryT_gsea,p.adj = 0.05)+ggtitle('Memory CD4+ T cells')

# Fig 3 healthy beta cell #########
# results for healthy pancreas
# filter = 0.7
beta_cor  <- robustness(filter_0.7_ES_1000)
beta_network  <- scRoGG::coExp_network(beta_cor[,1:3])
beta_pathway <- scRoGG::net_gsea(beta_network[['all']],species='Homo sapiens',category='H')
pway1 <- plotPways(beta_pathway)+ggtitle('scRoGG')
beta_pathway1 <- beta_pathway[beta_pathway$padj<0.2,]
beta_pathway1$leadingEdge <- sapply(beta_pathway1$leadingEdge , function(a)paste0(a,collapse='/'))
#write.csv(beta_pathway1,'scRoGG_beta_pathways.csv')
# keep as same as what scRoGG used 
index <-  apply(beta_healthy,1,function(y)sum(y>0))
beta_healthy_sub <- beta_healthy[index>(0.7*ncol(beta_healthy)),]
# remove specific gene names - mito
MT <- grep(pattern = "^MT-|^mt-", x = rownames(beta_healthy_sub), value = TRUE)
beta_healthy_sub <- beta_healthy_sub[setdiff(rownames(beta_healthy_sub),MT),]
# pearson + raw
library(Hmisc)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  x = data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
  x$p.adj <- p.adjust(x$p)
  x
}
pearson_beta <- rcorr(t(beta_healthy_sub)) 
pearson_beta <- flattenCorrMatrix(pearson_beta$r, pearson_beta$P)
pearson_beta_final <- pearson_beta[pearson_beta$p.adj<0.05,]
library('scRoGG')
pearson_beta_network <- coExp_network(pearson_beta_final[,1:3],n_networks = 1)
pearson_beta_pways <- net_gsea(pearson_beta_network[['all']],species='Homo sapiens',category='H')
pway2 <- plotPways(pearson_beta_pways)+ggtitle('Pearson_Raw')
beta_pathway2<- pearson_beta_pways[pearson_beta_pways$padj<0.2,]
beta_pathway2$leadingEdge <- sapply(beta_pathway2$leadingEdge , function(a)paste0(a,collapse='/'))
write.csv(beta_pathway2,'Pearson_Raw_pathways.csv')
# pearson + scTransformed
pan3 <- CreateSeuratObject(t(beta_raw),meta.data =beta_meta)

pan3_healthy <- subset(pan3, subset = diabetes =='ND')
pan3_healthy <- SCTransform(pan3_healthy, vars.to.regress = 'dataset')

SCT_beta <- as.matrix(pan3_healthy@assays$SCT@data)
index <-  apply(SCT_beta,1,function(y)sum(y>0))
SCT_beta <- SCT_beta[index>(0.7*ncol(SCT_beta)),]
# remove specific gene names - mito
MT <- grep(pattern = "^MT-|^mt-", x = rownames(SCT_beta), value = TRUE)
SCT_beta <- SCT_beta[setdiff(rownames(SCT_beta),MT),]

pearsonSCT_beta <- rcorr(t(SCT_beta)) 
pearsonSCT_beta <- flattenCorrMatrix(pearsonSCT_beta$r, pearsonSCT_beta$P)
pearsonSCT_beta_final <- pearsonSCT_beta[pearsonSCT_beta$p.adj<0.05,]
  
pearsonSCT_beta_network <- coExp_network(pearsonSCT_beta_final[,1:3],n_networks = 1)
pearsonSCT_beta_pways <- net_gsea(pearsonSCT_beta_network[['all']],species='Homo sapiens',category='H')

# compare three results
# compare pathways
# write.csv(beta_cor2,'scRoGG_beta_correlated_genepair.csv')
# write.csv(pearson_beta_final,'PearsonRaw_beta_correlated_genepair.csv')
# write.csv(pearsonSCT_beta_final,'PearsonBC_beta_correlated_genepair.csv')

ggpubr::ggarrange(pway1,pway2, common.legend = T, legend = 'right')

# extract top 500 from Pearson data
pearson_500 <- pearson_beta_final[order(abs(pearson_beta_final$cor),decreasing=T),][1:500,]
pearson_500_net <- coExp_network(pearson_500[,1:3])
pearsonSCT_500 <- pearsonSCT_beta_final[order(abs(pearsonSCT_beta_final$cor),decreasing=T),][1:500,]
pearsonSCT_500_net <- coExp_network(pearsonSCT_500[,1:3])

# overlap of gene pairs under different number of ES genes
scRoGG_500 <- apply(beta_cor[,1:2],1,function(s)paste0(s[1],'_',s[2]))
pearson_500n <- apply(pearson_500[,1:2],1,function(s)paste0(s[1],'_',s[2]))
pearsonSCT_500n <- apply(pearsonSCT_500[,1:2],1,function(s)paste0(s[1],'_',s[2]))

library(UpSetR)
pair_list <- list(scRoGG = scRoGG_500,
                     Pearson_Raw = pearson_500n,
                     Pearson_BC = pearsonSCT_500n)
upset(fromList(pair_list), order.by = 'freq', text.scale	= 2,main.bar.color=c('gray23','gray23','gray23','gray23','brown','darkgreen')) 
# overlap between scRoGG and Pearson_BC
shared <- intersect(scRoGG_500,pearsonSCT_500n)
shared <- t(sapply(shared,function(d)stringr::str_split(d,'_')[[1]])
)
# compare networks scRoGG and Pearson_BC
scRoGG_net <- beta_network[['all']] # actually use beta_network2
PearsonRaw_net <- pearson_500_net[['all']]
PearsonBC_net <- pearsonSCT_500_net[['all']]

library(igraph)
assortativity.degree(scRoGG_net, directed = F)
assortativity.degree(PearsonRaw_net, directed = F)
assortativity.degree(PearsonBC_net, directed = F)

edge_density(scRoGG_net,loops = F)
edge_density(PearsonRaw_net,loops = F)
edge_density(PearsonBC_net,loops = F)

diameter(scRoGG_net, weights = NA)
diameter(PearsonRaw_net, weights = NA)
diameter(PearsonBC_net, weights = NA)

boxplot(degree(scRoGG_net,normalized = T),degree(PearsonRaw_net,normalized = T),degree(PearsonBC_net,normalized = T))
 
# compare different number of ES
corpair_list <- list(ES_50 = corpair_50,
                     ES_150 = corpair_150,
                     ES_350 = corpair_350)
upset(fromList(corpair_list), order.by = 'freq',main.bar.color=c('gray23','gray23', '#8E1600','gray23','gray23'),text.scale	=2)

shared_corpair <- Reduce(intersect,corpair_list)
shared_pair <- sapply(shared_corpair,function(d)stringr::str_split(d,'_')[[1]])
shared_corpair_df <- apply(shared_pair,2,function(z)beta_cor2[beta_cor2$gene1==z[1]&beta_cor2$gene2==z[2],])
shared_corpair_df <- do.call('rbind',shared_corpair_df)

share_network <- scRoGG::coExp_network(shared_corpair_df[,1:3])
share_pathway <- scRoGG::net_gsea(share_network[['all']],species='Homo sapiens',category='H')
plotPways(share_pathway,p.adj = 0.2)+scale_fill_gradient2(low = "black",mid = '#8E1600', high = "red")
write.csv(beta_cor,'scRoGG_ES50_correlated_genepair.csv')
write.csv(beta_cor1,'scRoGG_ES150_correlated_genepair.csv')
write.csv(beta_cor2,'scRoGG_ES350_correlated_genepair.csv')
share_pathway$leadingEdge <- sapply(share_pathway$leadingEdge , function(a)paste0(a,collapse='/'))
write.csv(share_pathway,'scRoGG_ESshare_pathways.csv')

# Fig 3 healthy vs T2D beta cell #####
names(filter_0.7_ES_1000)[3] <- 'bi_nonzero'
names(betaT2D_cor)[3] <- 'bi_nonzero'
beta_T2D_list <- list(filter_0.7_ES_1000,betaT2D_cor)
beta_T2D <- scRoGG::robustness2(beta_T2D_list, p.adj=0.1) 
write.csv(beta_T2D,'scRoGG_differentially_correlated_T2D.csv')
 
beta_T2D_net <- scRoGG::coExp_network(beta_T2D[,c(1,2,4)],n_networks = 2)
p3c  <-plotCoExp(dat = filter_0.7_ES_1000$transformed_data,'CNTN1','APLP1') 
p3c1 <- plotCoExp(dat = betaT2D_cor$transformed_data,'CNTN1','APLP1')
ggpubr::ggarrange(p3c,p3c1,labels = c('Healthy','T2D'))

p3d  <-plotCoExp(dat = filter_0.7_ES_1000$transformed_data,'TULP4','RAB21') 
p3d1 <- plotCoExp(dat = betaT2D_cor$transformed_data,'TULP4','RAB21')
ggpubr::ggarrange(p3d,p3d1,labels = c('He ','T2D'))
# to construct network 
betandT2D_pathway <- net_gsea(network = beta_T2D_net[['all']],species = 'Homo sapiens', category = 'C5',subcategory = 'GO:BP',minSize = 10)
plotPways(betandT2D_pathway)
betandT2D_pathway.s<- betandT2D_pathway[betandT2D_pathway$padj<0.2,]
betandT2D_pathway.s$leadingEdge <- sapply(betandT2D_pathway.s$leadingEdge , function(a)paste0(a,collapse='/'))
#write.csv(betandT2D_pathway.s,'scRoGG_healthyT2D_pathways.csv')

# run EdgeR 
# Create single cell experiment object
index <-  apply(betaT2D,1,function(y)sum(y>0))
beta_T2D_sub <- betaT2D[index>(0.7*ncol(betaT2D)),]
# remove specific gene names - mito
MT <- grep(pattern = "^MT-|^mt-", x = rownames(betaT2D), value = TRUE)
beta_T2D_sub <- beta_T2D_sub[setdiff(rownames(beta_T2D_sub),MT),]
s_genes <- intersect(rownames(beta_healthy_sub),rownames(beta_T2D_sub))
beta_all<- cbind(beta_healthy_sub[match(s_genes,rownames(beta_healthy_sub)),],beta_T2D_sub[match(s_genes,rownames(beta_T2D_sub)),])

beta_meta <- beta_meta[match(colnames(beta_all),beta_meta$X),]
rownames(beta_meta) <- beta_meta$X
sce_beta <- SingleCellExperiment(assays = list(counts = beta_all),
                            colData = beta_meta)

# # Identify groups for aggregation of counts
# groups <- colData(sce_beta)[, c( "diabetes",'dataset')]
# summed <- aggregateAcrossCells(sce_beta, 
#                                ids=groups)
# summed
library(edgeR)
y.all <- DGEList(counts(sce_beta), samples = colData(sce_beta),group = colData(sce_beta)$diabetes)
y.all <- calcNormFactors(y.all)
design <- model.matrix(~colData(sce_beta)$diabetes+colData(sce_beta)$ dataset)
y.all<- estimateDisp(y.all,design)
fit <- glmFit(y.all, design)
lrt <- glmLRT(fit,coef=2)
EdgeR <- topTags(lrt,sort.by="logFC", n=Inf)
EdgeR_DE_genes <- EdgeR$table[EdgeR$table$FDR<0.1,]
#write.csv(EdgeR_DE_genes,'EdgeR_DE_genes.csv')
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(sce_beta),
                              colData = colData(sce_beta),
                              design= ~ diabetes+dataset)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res_DESeq2 <- as.data.frame(results(dds, name="diabetes_T2D_vs_ND"))

DESeq2_DE_genes <- res_DESeq2[res_DESeq2$padj<0.1,]
#write.csv(DESeq2_DE_genes,'DESeq2_DE_genes.csv')

Gene_list <- list(scRoGG = union(beta_T2D$gene1,beta_T2D$gene2),
                  edgeR = rownames(EdgeR_DE_genes),
                  DEseq2 = rownames(DESeq2_DE_genes))
#write.csv(union(beta_T2D$gene1,beta_T2D$gene2),'scRoGG_DC_genes.csv')
library(UpSetR)
upset(fromList(Gene_list), order.by = 'freq', text.scale = 2)
Reduce(intersect,Gene_list)

CNTN1_edgeR <- cpm(y.all, log=T, prior.count=1)[rownames(y.all)=='CNTN1',]
APLP1_edgeR <- cpm(y.all, log=T, prior.count=1)[rownames(y.all)=='APLP1',]
CNTN1_DEseq2 <-log2(counts(dds, normalized=TRUE)[rownames(y.all)=='CNTN1',]+1)
APLP1_DEseq2 <-log2(counts(dds, normalized=TRUE)[rownames(y.all)=='APLP1',]+1)
# supp. data to showcase the difference between DE and DC
DE.df <- data.frame(value = c(CNTN1_edgeR,APLP1_edgeR,CNTN1_DEseq2,APLP1_DEseq2),group = beta_meta$diabetes, method = rep(c('edegR','DEseq2'), each = 801*2),gene = rep(c('CNTN1','APLP1','CNTN1','APLP1'), each = 801))
library(ggplot2)
ggplot(DE.df, aes(x = group, y = value, fill=group))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  facet_grid(gene~method, scales = "free")+theme_bw()+xlab('')+ylab('Log-normalised gene expression')+ scale_fill_brewer(palette="Dark2")

# Fig 4 mouse endothelial arteriovenous zonation ####
# visualise in a tSNE
brain_3EC = readRDS('brain_3EC.rds')
library(scran)
brain_3EC <- computeSumFactors(brain_3EC)
brain_3EC <- logNormCounts(brain_3EC)
dec <- modelGeneVar(brain_3EC) 
sced <- denoisePCA(brain_3EC, dec, subset.row=getTopHVGs(dec, prop=0.1)) # 741 HVGs
ncol(reducedDim(sced, "PCA")) #30 
brain_3EC <- fixedPCA(brain_3EC, subset.row=getTopHVGs(dec, prop=0.1))
reducedDim(brain_3EC, "PCAsub") <- reducedDim(brain_3EC, "PCA")[,1:30]
library(scater)
brain_3EC <- runTSNE(brain_3EC, dimred="PCAsub")
plotTSNE(brain_3EC, colour_by="annoated.cell.types", text_by="annoated.cell.types")+ scale_color_manual(values=c("#cd9797",  "#b09db9","#9fa4c6"))                                      

#identify patterns
library(scRoGG)
aEC <- brain_3EC[,brain_3EC$annoated.cell.types=='aEC']
aEC_cor <- scRoGG(dat=assay(aEC),filter = 0.2,ES_number = 200,org = 'mmu')
capEC <- brain_3EC[,brain_3EC$annoated.cell.types=='capilEC']
capEC_cor <- scRoGG(dat=assay(capEC ),filter = 0.2,ES_number = 200,org = 'mmu')
vEC <- brain_3EC[,brain_3EC$annoated.cell.types=='vEC']
vEC_cor <- scRoGG(dat=assay(vEC),filter = 0.2,ES_number = 200,org = 'mmu')

EC_list <- list(aEC_cor1,capEC_cor1,vEC_cor1)
ECs_cordiff <- robustness2(EC_list,p.adj = 0.1)

# different patterns
# sign difference 
pattern5 <- ECs_cordiff[ECs_cordiff$type=='sign', ]
# DC change 
ECs_cordiff_DC <- ECs_cordiff[ECs_cordiff$type=='DC',]
# extra column: 3/2
ECs_cordiff_DC3 <- ECs_cordiff_DC$diff2/ECs_cordiff_DC$diff1
#group1<group2<group3
pattern1 <- ECs_cordiff_DC[ECs_cordiff_DC$diff1>1&ECs_cordiff_DC3>1, ]
coExp_network(pattern1[,c(1,2,5)],n_networks = 1)
#group1<group2>group3
pattern2 <- ECs_cordiff_DC[ECs_cordiff_DC$diff1>1&ECs_cordiff_DC3<1, ]
#group1>group2<group3
pattern3 <- ECs_cordiff_DC[ECs_cordiff_DC$diff1<1&ECs_cordiff_DC3>1, ]
#group1>group2>group3
pattern4 <- ECs_cordiff_DC[ECs_cordiff_DC$diff1<1&ECs_cordiff_DC3<1, ]

# Over-representation analysis with chi-square test 
library(GO.db)
library(org.Mm.eg.db)
keys = keys(org.Mm.eg.db)
columns(org.Mm.eg.db)
GO_info = AnnotationDbi::select(org.Mm.eg.db, keys=keys, columns = c("SYMBOL", "GO"))
keep = GO_info$SYMBOL %in% rownames(brain_3EC)
table(keep)
GO_info_filt = GO_info[keep,]
# at least 10 genes in the term and less than 500
allTerms = names(which(table(GO_info_filt$GO) >= 5 & table(GO_info_filt$GO) <= 500))

GO_info_terms = AnnotationDbi::select(GO.db, columns = columns(GO.db), keys = allTerms)
rownames(GO_info_terms) <- GO_info_terms$GOID

allTermNames = GO_info_terms[allTerms, "TERM"]
names(allTermNames) <- allTerms

GO_list = sapply(allTerms, function(term) {
  genes = GO_info_filt[GO_info_filt$GO == term, "SYMBOL"]
  return(sort(unique(genes[!is.na(genes)])))
}, simplify = FALSE)
names(GO_list) <- allTermNames
GO_list_for_patterns = GO_list[unlist(lapply(GO_list, function(x) any(x %in% union(ECs_cordiff$gene1,ECs_cordiff$gene2))))]
genesetGOtest = function(set, universe, termsList, n_plot=10) {
  # set is the character vector of genes in geneset
  # universe is character vector of genes to include in universe
  # termsList is a list of the pathways
  
  termsListfiltered = lapply(termsList, function(x) return(intersect(universe, x)))
  keepForTesting = unlist(lapply(termsListfiltered, length)) >= 5 & unlist(lapply(termsListfiltered, length)) <= 500
  
  pval = sapply(1:length(termsList), function(i){
    
    if (!keepForTesting[i]) return(1)
    
    termGenes = termsListfiltered[[i]]
    return(fisher.test(table(universe %in% set, universe %in% termGenes), alt = "g")$p.value)
  })
  names(pval) <- names(termsListfiltered)
  name <- sapply(termsListfiltered,function(y){ 
    intersect(y,set) 
  })
  df = data.frame(pathway = factor(names(pval), levels = c(names(pval), "")),  size = sapply(name, length),name = sapply(name, function(s)paste0(s,collapse = '/')),
                  pval = pval,
                  padj = p.adjust(pval, method = "BH"))
  
  df_sorted = reshape::sort_df(df, "pval")[1:n_plot,]
  df_sorted$pathway = factor(df_sorted$pathway, levels =  rev(df_sorted$pathway))
  rownames(df_sorted) <- NULL
  return(df_sorted)
}
# identify enriched pathways
pattern1_pway <- genesetGOtest(union(pattern1$gene1,pattern1$gene2), rownames(brain_3EC), GO_list_for_patterns)
pattern2_pway <- genesetGOtest(union(pattern2$gene1,pattern2$gene2), rownames(brain_3EC), GO_list_for_patterns)
pattern3_pway <- genesetGOtest(union(pattern3$gene1,pattern3$gene2), rownames(brain_3EC), GO_list_for_patterns)
pattern4_pway <- genesetGOtest(union(pattern4$gene1,pattern4$gene2), rownames(brain_3EC), GO_list_for_patterns)
pattern5_pway <- genesetGOtest(union(pattern5$gene1,pattern5$gene2), rownames(brain_3EC), GO_list_for_patterns)
library(ggplot2)
library('scRoGG')
plotPways(pattern1_pway,p.adj = 0.1)+ ggtitle('Pattern I')
plotPways(pattern2_pway,p.adj = 0.1) # no significance
pattern3.plot <- plotPways(pattern23_pway,p.adj = 0.1)+ ggtitle('Pattern III')  
pattern4.plot <-plotPways(pattern4_pway,p.adj = 0.1)+ ggtitle('Pattern IV') 
pattern5.plot <-plotPways(pattern5_pway,p.adj = 0.1)+ ggtitle('Pattern V') 
ggpubr::ggarrange(pattern3.plot,pattern4.plot,pattern5.plot, nrow = 1, ncol = 3, common.legend = T, legend = 'right')

# visualise the overall coexpression network
ECs_net <- scRoGG::coExp_network(ECs_cordiff[,c(1,2,4)],n_networks = 1,min_size = 3)
library(igraph)
graph = ECs_net[['all']]
edgelist= as_long_data_frame(graph)
library(RColorBrewer)
brewer.pal(n = 5, name = "Dark2")
type <- apply(edgelist,1,function(s){
  ss <- paste0(s[4],'_',s[5])
  if(ss %in% apply(pattern1,1,function(o)paste0(o[1],'_',o[2]))){
    'P1'
  }else if (ss %in% apply(pattern2,1,function(o)paste0(o[1],'_',o[2]))){
    'P2'
  }else if (ss %in% apply(pattern3,1,function(o)paste0(o[1],'_',o[2]))){
    'P3'
  }else if (ss %in% apply(pattern4,1,function(o)paste0(o[1],'_',o[2]))){
    'P4'} else {'P5'}
})
E(graph)$type = brewer.pal(n = 5, name = "Dark2")[as.factor(type)]
plot(graph,layout=layout.fruchterman.reingold, edge.width=4, edge.color= E(graph)$type,vertex.size=10,vertex.label.color="black", vertex.color = "#cc6600b3",          # Node color
     vertex.frame.color = "white",
     vertex.label.cex=0.7)

# supplementary figures
pattern1_net <- coExp_network(pattern1[c(1,2,4)])[['all']]
plot(pattern1_net,layout=layout.fruchterman.reingold, edge.width=4, edge.color= '#D4B1A4',vertex.size=10,vertex.label.color="black", vertex.color = "#C7D0EE",          # Node color
     vertex.frame.color = "white",
     vertex.label.cex=0.7)

pattern2_net <- coExp_network(pattern2[c(1,2,4)])[['all']]
plot(pattern2_net,layout=layout.fruchterman.reingold, edge.width=4, edge.color= '#D4B1A4',vertex.size=10,vertex.label.color="black", vertex.color = "#C7D0EE",          # Node color
     vertex.frame.color = "white",
     vertex.label.cex=0.7)

pattern3_net <- coExp_network(pattern3[c(1,2,4)])[['all']]
plot(pattern3_net,layout=layout.fruchterman.reingold, edge.width=4, edge.color= '#D4B1A4',vertex.size=10,vertex.label.color="black", vertex.color = "#C7D0EE",          # Node color
     vertex.frame.color = "white",
     vertex.label.cex=0.7)

pattern4_net <- coExp_network(pattern4[c(1,2,4)])[['all']]
plot(pattern4_net,layout=layout.fruchterman.reingold, edge.width=4, edge.color= '#D4B1A4',vertex.size=10,vertex.label.color="black", vertex.color = "#C7D0EE",          # Node color
     vertex.frame.color = "white",
     vertex.label.cex=0.7)

pattern5_net <- coExp_network(pattern5[c(1,2,4)])[['all']]
plot(pattern5_net,layout=layout.fruchterman.reingold, edge.width=4, edge.color= '#D4B1A4',vertex.size=10,vertex.label.color="black", vertex.color = "#C7D0EE",          # Node color
     vertex.frame.color = "white",
     vertex.label.cex=0.7)

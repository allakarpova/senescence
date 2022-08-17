# ssgsea scores for ccrcc
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)
library(GSVA)

signatures <- fread('~/lab_Ding/work/single_cell/senescence/Papers/Unmasking_transcriptional_heterogeneity/processed_sen_signatures_Segura_paper.txt', data.table = F)
signatures.pos <- subset(signatures, log2FC > 0)
signatures.neg <- subset(signatures, log2FC < 0)
wd <- '~/lab_Ding/work/single_cell/senescence/ccrcc/normal_samples/ssgsea'
dir.create(wd, showWarnings = F, recursive = T)
setwd(wd)

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res.7.Normal epithelial cells_ccrcc.2.sample.RDS')

DefaultAssay(panc) <- 'RNA'
FeaturePlot(object = panc, features = 'CCND1')
DimPlot(object = panc)

panc <- FindClusters(object = panc, resolution = 2)
cell.type = 'Normal_epithelial'


# make it readable for gsva function
gene.sets.list <- split (signatures.pos$Gene_names, signatures.pos$ont)
markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), 
                           Cell_cycle = c('CCND1', 'CDKN1A','CDKN2A', 'TP53'), 
                           DDR = c('H2AFX','H2AFY', 'TP53BP1'),
                           Proliferation = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1'),
                           SASP = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'))
gene.sets.list <- c(gene.sets.list, markers.from.email)

expr.mat <- FetchData(panc, vars = rownames(panc), slot = 'data') %>% t() %>% scale()

###### 
ss.results.prot <- gsva(expr = as.matrix(expr.mat), 
                        gset.idx.list = gene.sets.list,
                        method=c( "ssgsea"),
                        kcdf=c("Gaussian"),
                        abs.ranking=FALSE,
                        min.sz=1,
                        max.sz=Inf,
                        parallel.sz=0,
                        parallel.type="SOCK",
                        mx.diff=TRUE,
                        tau=0.25,
                        ssgsea.norm=T,
                        verbose=TRUE)
fwrite(data.frame(ss.results.prot), 'Ccrcc_normal_epith_senescence_signatures_upregulated_genes_scaled.txt', sep = '\t', row.names = T)

# make it readable for gsva function
gene.sets.list <- split (signatures.neg$Gene_names, signatures.neg$ont)

###### 
ss.results.prot <- gsva(expr = as.matrix(expr.mat), 
                        gset.idx.list = gene.sets.list,
                        method=c( "ssgsea"),
                        kcdf=c("Gaussian"),
                        abs.ranking=FALSE,
                        min.sz=1,
                        max.sz=Inf,
                        parallel.sz=0,
                        parallel.type="SOCK",
                        mx.diff=TRUE,
                        tau=0.25,
                        ssgsea.norm=T,
                        verbose=TRUE)
fwrite(data.frame(ss.results.prot), 'Ccrcc_normal_epith_senescence_signatures_downregulated_genes_scaled.txt', sep = '\t', row.names = T)


pos.scores <- fread('Ccrcc_normal_epith_senescence_signatures_upregulated_genes.txt', data.table = F) %>% data.frame(row.names = 1); neg.scores <- fread('Ccrcc_normal_epith_senescence_signatures_downregulated_genes.txt' )%>% data.frame(row.names = 1)
pos.scores <- cbind(Embeddings(panc, reduction = 'umap'),seurat_clusters = panc$seurat_clusters, data.frame(t(pos.scores)), dir = 'pos')
neg.scores <- cbind(Embeddings(panc, reduction = 'umap'),seurat_clusters = panc$seurat_clusters, data.frame(t(neg.scores)), dir = 'neg')
toplot.melt <- rbind(melt(pos.scores, id.vars = c('seurat_clusters', 'UMAP_1', 'UMAP_2', 'dir')), 
                     melt(neg.scores, id.vars = c('seurat_clusters', 'UMAP_1', 'UMAP_2', 'dir')))


#toplot <- cbind(Embeddings(panc, reduction = 'umap'),seurat_clusters = panc$seurat_clusters, data.frame(t(ss.results.prot)))
#toplot.melt <- melt(toplot, id.vars = c('seurat_clusters', 'UMAP_1', 'UMAP_2'))
colnames(toplot)

ggplot(data = toplot.melt, aes(x = value, fill = dir)) +
  geom_density(alpha = 0.5) +
  facet_grid(seurat_clusters~variable, scales = 'free') +
  theme_classic()# +
#scale_x_continuous( breaks = c(-0.3, 0, 0.3, 0.6), )
ggsave('Density_ccrcc_plot_signatures_positive_negative_genes.pdf', useDingbats = F, width = 15, height = 10)


toplot <- pos.scores
toplot.melt <- melt(toplot, id.vars = c('seurat_clusters', 'UMAP_1', 'UMAP_2'))
for (s in names(gene.sets.list)) {
  ggplot(data = toplot,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = s)) +
    geom_point( size = 1) +
    theme_cowplot() +
    #ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = paste(s, "\nssggsea score"), low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c')
  ggsave(paste0('Featureplot_ssgsea_pos_', s, '.pdf'), useDingbats = F, width = 7, height = 4.5)
  
}

getwd()
toplot <- neg.scores
toplot.melt <- melt(toplot, id.vars = c('seurat_clusters', 'UMAP_1', 'UMAP_2'))
for (s in names(gene.sets.list)) {
  ggplot(data = toplot,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = s)) +
    geom_point( size = 1) +
    theme_cowplot() +
    #ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = paste(s, "\nssggsea score"), low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c')
  ggsave(paste0('Featureplot_ssgsea_neg_', s, '.pdf'), useDingbats = F, width = 7, height = 4.5)
  
}

colnames(neg.scores)[4:(ncol(neg.scores)-1)] <- paste0(colnames(neg.scores)[4:(ncol(neg.scores)-1)], '.neg')
panc <- AddMetaData(panc, neg.scores[4:(ncol(neg.scores)-1)])

colnames(pos.scores)[4:(ncol(pos.scores)-1)] <- paste0(colnames(pos.scores)[4:(ncol(pos.scores)-1)], '.pos')
panc <- AddMetaData(panc, pos.scores[4:(ncol(pos.scores)-1)])

FeaturePlot(panc, features = c('Universal_sen.neg', 'Universal_sen.pos'), blend = T, cols = c('#1f78b4', '#e41a1c'), order = T)




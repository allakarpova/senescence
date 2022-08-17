library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)
library(ggplot2)

filter_counts <- function( counts, cell.pct = 0.05){
  cells.id <- counts@Dimnames[[2]]
  nCells <- length(cells.id)
  
  Cellmin <- cell.pct * nCells
  nonzero <- counts > 0
  # Sums all TRUE values and returns TRUE if more than Cellmin TRUE values per gene
  keep_genes <- Matrix::rowSums(counts) >= Cellmin
  keep_genes <- names(keep_genes)[keep_genes]
  return(keep_genes)
}


#panc - object
# s - signature name to correlate with
# ct - cell type that this object contains
do.cor.good.cover <- function(panc, s, ct, dot.size = 1, cell.pct.cutoff = 0.05, cor.adjust.method = 'fdr') {
  print(s)
  #s = 'Fibroblast_replicative_sen'
  genes.in.signature <- subset(signatures, ont == s)$Gene_names
  #get genes with good coverage
  counts <- GetAssayData(object = panc, slot = "counts")
  keep_genes <- filter_counts(counts, cell.pct = cell.pct.cutoff)
  genes.in.signature.good.coverage <- intersect(genes.in.signature, keep_genes)
  cat(paste0('Genes from a signature total ', length(genes.in.signature), '\n'))
  cat(paste0('Genes from a signature with good coverage ', length(genes.in.signature.good.coverage), '\n'))
  #get expression data and scale within each cell
  expr.mat <- FetchData(panc, vars = genes.in.signature.good.coverage, slot = 'data') %>% t()# %>% scale()
  #dim(expr.mat)
  fc.sign <- subset(signatures, ont == s)
  rownames(fc.sign) <- fc.sign$Gene_names
  fc.sign <- fc.sign[genes.in.signature.good.coverage,]
  #dim(fc.sign)
  if(dim(expr.mat)[1] != dim(fc.sign)[1]) {
    fc.sign <- subset(fc.sign, Gene_names %in% rownames(expr.mat))
    stopifnot(fc.sign$Gene_names== rownames(expr.mat))
  }
  
  cor.result <- psych::corr.test(expr.mat, fc.sign$log2FC, method = 'spearman', adjust = cor.adjust.method)
  df <- cbind('rho' = cor.result$r , 'p.val' = cor.result$p) %>% data.frame(check.rows = F, check.names = F)
  colnames(df) <- c('rho', 'p.val')
  
  toplot <- cbind(Embeddings(panc, reduction = 'umap'), df)
  max.cor <- max(abs(toplot$rho))
  
  ggplot(data = toplot,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho')) +
    geom_point(
      size = dot.size) +
    theme_cowplot() +
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0('Featureplot_', ct, '_cor_', s, '.pdf'), useDingbats = F, width = 5, height = 4.5)
  
  #color dots only if p.value is < 0.05
  ggplot() +
    geom_point(data = subset(toplot, p.val > 0.05), 
               aes_string(x = 'UMAP_1', y = 'UMAP_2'), 
               size = 0.5, color = 'lightgrey') +
    geom_point(data = subset(toplot, p.val < 0.05),  
               aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho'),
               size = dot.size) +
    theme_cowplot() +
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0('Featureplot_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 5, height = 4.5)
  
}

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/'
dir.create(wd, showWarnings = F)
setwd(wd)

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res.7.mCAF_htan.brca.sc.v10.RDS')
panc <- FindClusters(panc, resolution = 2)
cell.type = 'mCAF'
Idents(panc) <- panc$cell_type_specific
DefaultAssay(panc) <- 'RNA'


pdf(paste0('Dimplots_',cell.type,'.pdf'), width = 8.5, height =6, useDingbats = F)
print(DimPlot(panc, label = T, group.by = 'seurat_clusters') & coord_fixed())
print(DimPlot(panc, group.by = 'cell_type_specific', label = T)& coord_fixed())
print(DimPlot(panc, group.by = 'sample', label = T)& coord_fixed())
dev.off()


signatures <- fread('~/lab_Ding/work/single_cell/senescence/Papers/Unmasking_transcriptional_heterogeneity/processed_sen_signatures_Segura_paper.txt', data.table = F)
signatures$ont %>% unique
signatures %>% group_by(ont) %>% tally()
intersect(subset(signatures, ont== 'Fibroblasts_sen')$Gene_names, 
          subset(signatures, ont== 'Fibroblast_meta_sen')$Gene_names)




for(s in signatures$ont %>% unique) {
  do.cor (panc, s, ct = cell.type, dot.size = 0.5)
}

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/with_fdr_good_coverage'
dir.create(wd, showWarnings = F)
setwd(wd)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc, s, ct = cell.type, dot.size = 0.5, cell.pct.cutoff = 0.05)
}

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/no_fdr_good_coverage1pct'
dir.create(wd, showWarnings = F)
setwd(wd)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc, s, ct = cell.type, dot.size = 0.5, cell.pct.cutoff = 0.01, cor.adjust.method = 'none')
}

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/no_fdr_good_coverage1pct_scaled_within_cell'
dir.create(wd, showWarnings = F)
setwd(wd)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc, s, ct = cell.type, dot.size = 0.5, cell.pct.cutoff = 0.01, cor.adjust.method = 'none')
}

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/no_fdr_good_coverage1pct_not_scaled'
dir.create(wd, showWarnings = F)
setwd(wd)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc, s, ct = cell.type, dot.size = 0.5, cell.pct.cutoff = 0.01, cor.adjust.method = 'none')
}

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/no_fdr_not_scaled'
dir.create(wd, showWarnings = F)
setwd(wd)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc, s, ct = cell.type, dot.size = 0.5, cell.pct.cutoff = 0.0, cor.adjust.method = 'none')
}
pdf(paste0('Featurplots_',cell.type,'.pdf'), width = 12, height =6, useDingbats = F)
print(FeaturePlot(panc, features = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), order = T, ncol = 4))
print(FeaturePlot(panc, features = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'), order = T, ncol = 4))
print(FeaturePlot(panc, features = c('LMNB1'), order = T)&coord_fixed())
dev.off()

panc$SCT_snn_res.2
getwd()
pdf(paste0('Featurplots_',cell.type,'.pdf'), width = 12, height =6, useDingbats = F)
DimPlot(panc, group.by = 'SCT_snn_res.2', label = T)&coord_fixed()
VlnPlot(object = panc, features = c('CDKN1A','CDKN2A', 'TP53', 'CCND1'), pt.size = 0.1, group.by = 'seurat_clusters', ncol = 3)
VlnPlot(object = panc, features = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1'), pt.size = 0.9, group.by = 'seurat_clusters', ncol = 3)
VlnPlot(object = panc, features = c('H2AFX','H2AFY', 'TP53BP1'), pt.size = 0.9, group.by = 'seurat_clusters', ncol = 3)
VlnPlot(object = panc, features = c('SERPINE1','LMNB1'), pt.size = 0.9, group.by = 'seurat_clusters', ncol = 3)
VlnPlot(object = panc, features = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'), pt.size = 0.9, group.by = 'seurat_clusters', ncol = 3)
dev.off()

getwd()
pdf(paste0('VLNplot_','toplot_mCAF','.pdf'), width = 12, height =6, useDingbats = F)
print(VlnPlot(object = panc, features = c('CDKN2A','H2AFX', 'SERPINE1', 'IL6','CCL2','IGFBP7'), pt.size = 0.1, group.by = 'seurat_clusters', ncol = 3))
dev.off()


cpdb.split.df <- fread('~/R_working_dir/data/ConsensusPathDB/CPDB_pathways_genes.2column.txt', data.table = F)
cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), 
                           Cell_cycle = c('CCND1', 'CDKN1A','CDKN2A', 'TP53'), 
                           DDR = c('H2AFX','H2AFY', 'TP53BP1'),
                           Proliferation = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1'),
                           SASP = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'))
db.sen <- split(cpdb.split.df.sen$gene, cpdb.split.df.sen$ont)
total.sen <- c(db.sen, markers.from.email)
names(total.sen) <- make.names(names(total.sen))
names(total.sen) <- gsub('_', '-', names(total.sen))

cat('calculate modules\n')
panc <- AddModuleScore(
  object = panc,
  assay = 'RNA',
  features = total.sen,
  ctrl = 5,
  name = 'senescence'
)
cat('done\n')

n.list <- length(total.sen)
colnames(panc@meta.data)[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)] <- names(total.sen)

FeaturePlot(panc, features = c('Generic', 'SASP'), order = T, min.cutoff = 0)
ggsave(paste0('Featureplot_',cell.type,'_Zhang_markers_module.pdf'), useDingbats = F, width = 8, height = 4)

FeaturePlot(panc, features = names(total.sen)[grepl('DAMAG', names(total.sen))], order = T, min.cutoff = 0, ncol = 5)
ggsave(paste0('Featureplot_',cell.type,'_all_damage_modules.pdf'), useDingbats = F, width = 30, height = 40, limitsize = F)

FeaturePlot(panc, features = names(total.sen)[grepl('SENESCE', names(total.sen))], order = T, min.cutoff = 0)
ggsave(paste0('Featureplot_',cell.type,'_all_senescence_modules.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)

FeaturePlot(panc, features = names(total.sen)[grepl('INFLAMMAT', names(total.sen))], order = T, min.cutoff = 0)
ggsave(paste0('Featureplot_',cell.type,'_all_inflammatory_modules.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)



### dropout rates
panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res.7.htan.brca.sc.v10.RDS')
wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/dropout'
dir.create(wd, showWarnings = F)
setwd(wd)
counts <- GetAssayData(object = panc, slot = "counts")
str(counts)
counts@Dim
zero.mat <- counts == 0
keep_genes <- data.frame(zeros = Matrix::rowSums(zero.mat), total = length(counts@Dimnames[[2]])) %>% mutate(pct = zeros/total, gene = rownames(.))
fwrite(keep_genes, 'Dropout_rate_htan_breast_norm.txt', sep = '\t')
hist(keep_genes$pct, xlim = c(0,0.5), ylim = c(0,100), breaks = 50, xlab = 'Dropout', ylab = 'Number of genes')


ggplot() +
  geom_histogram(data = keep_genes, aes (pct), color = 'white') +
  theme_cowplot() +
  labs(x = 'dropout', y = 'Number of genes', title = 'Dropout distribution in breast normal cells') +
  geom_point(data = subset(keep_genes, gene %in% markers.from.email$Generic | gene %in% markers.from.email$SASP), 
             aes(x = pct, y = -100), color = 'red', size = 0.1)
ggsave('Histogram_dropout_all_genes_breast.pdf', height = 5, width = 6, useDingbats = F)

fwrite(subset(keep_genes, gene %in% markers.from.email$Generic | gene %in% markers.from.email$SASP),
       'Dropout_rate_htan_breast_norm_specific_markers.txt', sep = '\t')





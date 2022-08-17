library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)

filter_counts <- function( counts, cell.pct = 0.05){
  cells.ct <- counts@Dimnames[[2]]
  nCells <- length(cells.ct)
  Cellmin <- cell.pct * nCells
  counts.filtered <- counts[,which(colnames(counts) %in% cells.ct)]
  nonzero <- counts.filtered > 0
  # Sums all TRUE values and returns TRUE if more than Cellmin TRUE values per gene
  keep_genes <- Matrix::rowSums(counts.filtered) >= Cellmin
  keep_genes <- names(keep_genes)[keep_genes]
  return(keep_genes)
}

#panc - object
# s - signature name to correlate with
# ct - cell type that this object contains
do.cor <- function(panc, s, ct, dot.size = 1) {
  print(s)
  genes.in.signature <- subset(signatures, ont == s)$Gene_names
  counts <- GetAssayData(object = panc, slot = "counts")
  str(counts)
  counts@Dimnames[[2]]
  
  
  
  
  
  
  expr.mat <- FetchData(panc, vars = genes.in.signature, slot = 'scale.data') %>% t()
  #dim(expr.mat)
  fc.sign <- subset(signatures, ont == s)
  #dim(fc.sign)
  if(dim(expr.mat)[1] != dim(fc.sign)[1]) {
    fc.sign <- subset(fc.sign, Gene_names %in% rownames(expr.mat))
  }
  cor.result <- psych::corr.test(expr.mat, fc.sign$log2FC, method = 'spearman', adjust = 'fdr')
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

wd <- '~/lab_Ding/work/single_cell/senescence/gbm/normal_samples/with_fdr'
dir.create(wd, showWarnings = F, recursive = T)
setwd(wd)

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res.7.htan.brca.sc.v10.RDS')
cell.type = 'Neurons'

pdf(paste0('Dimplots_',cell.type,'.pdf'), width = 8.5, height =6, useDingbats = F)
print(DimPlot(panc, label = T) & coord_fixed())
print(DimPlot(panc, group.by = 'cell_type_v20210128', label = T)& coord_fixed())
print(DimPlot(panc, group.by = 'case_id', label = T)& coord_fixed())
dev.off()


signatures <- fread('~/lab_Ding/work/single_cell/senescence/Papers/Unmasking_transcriptional_heterogeneity/processed_sen_signatures_Segura_paper.txt', data.table = F)
signatures$ont %>% unique

Idents(panc) <- panc$cell_type_v20210128
DefaultAssay(panc) <- 'RNA'

for(s in signatures$ont %>% unique) {
  do.cor (panc, s, ct = cell.type, dot.size = 0.5)
}


pdf(paste0('Featurplots_',cell.type,'.pdf'), width = 12, height =6, useDingbats = F)
print(FeaturePlot(panc, features = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), order = T, ncol = 4))
print(FeaturePlot(panc, features = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'), order = T, ncol = 4))
print(FeaturePlot(panc, features = c('LMNB1'), order = T)&coord_fixed())
dev.off()


# do gene set scoring as we did before
cpdb.split.df <- fread('~/R_working_dir/data/ConsensusPathDB/CPDB_pathways_genes.2column.txt', data.table = F)
cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), 
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



panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res1.gbm.v2.0.RDS')
wd <- '~/lab_Ding/work/single_cell/senescence/gbm/normal_samples/dropout'
dir.create(wd, showWarnings = F)
setwd(wd)
counts <- GetAssayData(object = panc, slot = "counts")
str(counts)
counts@Dimnames[[1]]
zero.mat <- counts == 0
keep_genes <- data.frame(zeros = Matrix::rowSums(zero.mat), total = length(counts@Dimnames[[2]])) %>% mutate(pct = zeros/total, gene = rownames(.))
fwrite(keep_genes, 'Dropout_rate_brain_normal.samples.txt', sep = '\t')
hist(keep_genes$pct, xlim = c(0,0.5), ylim = c(0,100), breaks = 50, xlab = 'Dropout', ylab = 'Number of genes')


ggplot() +
  geom_histogram(data = keep_genes, aes (pct), color = 'white') +
  theme_cowplot() +
  labs(x = 'dropout', y = 'Number of genes', title = 'Dropout distribution in brain normal cells') +
  geom_point(data = subset(keep_genes, gene %in% markers.from.email$Generic | gene %in% markers.from.email$SASP), 
             aes(x = pct, y = -100), color = 'red', size = 0.1)
ggsave('Histogram_dropout_brain_all_genes.pdf', height = 5, width = 6, useDingbats = F)

fwrite(subset(keep_genes, gene %in% markers.from.email$Generic | gene %in% markers.from.email$SASP),
       'Dropout_rate_htan_brain_norm_specific_markers.txt', sep = '\t')



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
do.cor.good.cover <- function(panc, s, ct, dot.size = 1, cell.pct.cutoff = 0.05, cor.adjust.method = 'fdr', folder ) {
  print(s)
  #s = 'Fibroblast_replicative_sen'
  genes.in.signature <- subset(signatures, ont == s)$Gene_names %>% tolower() %>% Caps()
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
  rownames(fc.sign) <- fc.sign$Gene_names %>% tolower() %>% Caps()
  fc.sign <- fc.sign[genes.in.signature.good.coverage,]
  #dim(fc.sign)
  if(dim(expr.mat)[1] != dim(fc.sign)[1]) {
    fc.sign <- subset(fc.sign, Gene_names %in% rownames(expr.mat))
    stopifnot(fc.sign$Gene_names== rownames(expr.mat))
  }
  
  cor.result <- psych::corr.test(expr.mat, fc.sign$log2FC, method = 'spearman', adjust = cor.adjust.method)
  df <- cbind('rho' = cor.result$r , 'p.val' = cor.result$p) %>% data.frame(check.rows = F, check.names = F)
  colnames(df) <- c('rho', 'p.val')
  panc <- AddMetaData(panc, df)
  
  toplot <- cbind(Embeddings(panc, reduction = 'umap'), orig.ident = panc$orig.ident, df)
  max.cor <- max(abs(toplot$rho))
  
  ggplot(data = toplot,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho')) +
    geom_point(
      size = dot.size) +
    theme_cowplot() +
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0(folder, 'Featureplot_', ct, '_cor_', s, '.pdf'), useDingbats = F, width = 5.5, height = 4.5)
  
  ggplot(data = toplot,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho')) +
    geom_point(
      size = dot.size) +
    theme_cowplot() +
    facet_wrap(~orig.ident, ncol = 3 )+
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0(folder, 'Featureplot_', ct, '_cor_', s, '_splitted.pdf'), useDingbats = F, width = 13, height = 5.5)
  
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
  ggsave(paste0(folder, 'Featureplot_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 5, height = 4.5)
  
  ggplot() +
    geom_point(data = subset(toplot, p.val > 0.05), 
               aes_string(x = 'UMAP_1', y = 'UMAP_2'), 
               size = 0.5, color = 'lightgrey') +
    geom_point(data = subset(toplot, p.val < 0.05),  
               aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho'),
               size = dot.size) +
    theme_cowplot() +
    facet_wrap(~orig.ident, ncol = 3 )+
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0(folder, 'Featureplot_', ct, '_cor_', s, '_significant_only_splitetd.pdf'), useDingbats = F, width = 15, height = 5.5)
  
  
  VlnPlot(object = panc, features = 'rho', pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
          assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
  ggsave(paste0(folder, 'Vlnplot_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 18, height = 8)
  
  VlnPlot(object = panc, features = 'rho', pt.size = 0.9, group.by = 'seurat_clusters', ncol = 1,split.by = 'orig.ident',
          assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
  ggsave(paste0(folder, 'Vlnplot_cliusters_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 18, height = 8)
  
}
Caps <- function(x) {
  s <- x
  toreturn = sapply(s, function(x) paste0(toupper(substring(x, 1,1)), substring(x, 2)) )
  return(as.character(toreturn))
}


signatures <- fread('~/lab_Ding/work/single_cell/senescence/Papers/Unmasking_transcriptional_heterogeneity/processed_sen_signatures_Segura_paper.txt', data.table = F)
panc <- readRDS('~/lab_Ding/work/single_cell/senescence/snRNA_combo/objects/mouse_liver_by_cell_type/Kuppfer.object.res.1.mouse_liver_old_young.rds')
cell.type = 'Kuppfer'

wd <- '~/lab_Ding/work/single_cell/senescence/snRNA_combo/senescence_correlation/Kuppfer/'
dir.create(wd, recursive = T)
setwd(wd)
panc$orig.ident
DimPlot(panc,label = T)
ggsave('Dimplot_cluster.pdf', width = 5.5, height = 5, useDingbats = F)

DimPlot(panc, group.by ='Cell_type',  label = T)
ggsave('Dimplot_cell_type.pdf', width = 6, height = 5, useDingbats = F)

DimPlot(panc, group.by ='Age',  label = T, cols = 'Paired')
ggsave('Dimplot_age.pdf', width = 6, height = 5, useDingbats = F)


fd <- 'no_fdr_good_coverage/'
dir.create(fd, showWarnings = F, recursive = T)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc = panc, s = s, ct = cell.type, dot.size = 1, cell.pct.cutoff = 0.05,cor.adjust.method = 'none', folder = fd)
}

fd <- 'no_fdr_coverage1pct/'
dir.create(fd, showWarnings = F, recursive = T)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc = panc, s = s, ct = cell.type, dot.size = 1, cell.pct.cutoff = 0.01,cor.adjust.method = 'none', folder = fd)
}

fd <- 'with_fdr_good_coverage/'
dir.create(fd, showWarnings = F, recursive = T)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc = panc, s = s, ct = cell.type, dot.size = 0.5, cell.pct.cutoff = 0.05,cor.adjust.method = 'fdr', folder = fd)
}

fd <- 'with_fdr_coverage1pct/'
dir.create(fd, showWarnings = F, recursive = T)
for(s in signatures$ont %>% unique) {
  do.cor.good.cover (panc = panc, s = s, ct = cell.type, dot.size = 1, cell.pct.cutoff = 0.01,cor.adjust.method = 'fdr', folder = fd)
}






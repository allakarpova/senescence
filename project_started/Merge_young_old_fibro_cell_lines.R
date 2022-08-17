#senescence for GBM non tumor cells 
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)
library(paletteer)
library(viridis)
library(tidyverse)


markers.from.email <- list(Generic = c('CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), 
                           SASP = c('IL1A', 'IL6', 'CXCL8', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'))
senmayo <- fread('~/lab_Ding/work/Sennet/Papers/SenMayo/SenMayo_geneset.txt')

obj.y <- '~/lab_Ding/work/Sennet/rds_objects/fixed/20211214-WS1_PD45N1Z1_1Bn1_1/20211214-WS1_PD45N1Z1_1Bn1_1_processed.rds'
obj.o <- '~/lab_Ding/work/Sennet/rds_objects/fixed/20211110-WS1_PD102N1Z1_1Bn1_1/20211110-WS1_PD102N1Z1_1Bn1_1_processed.rds'


out_path <- '~/lab_Ding/work/Sennet/analysis/test_cell_line_snRNA_old_vs_young_fixed/'
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)
add_filename <- 'snRNA_old_vs_young'


cat('opening object...\n')
obj.y <- readRDS(obj.y)
obj.y$Age <- 'Young'
obj.o <- readRDS(obj.o)
obj.o$Age <- 'Old'
cat('done\n')

DefaultAssay(obj.y) <- 'RNA'
DefaultAssay(obj.o) <- 'RNA'

obj.y <- DietSeurat(obj.y, assays = 'RNA', counts = T, data = T, scale.data = F)
obj.o <- DietSeurat(obj.o, assays = 'RNA', counts = T, data = T, scale.data = F)

merged <- merge(obj.y, obj.o, add.cell.ids = c('Young', 'Old'))
rm(obj.y, obj.o)
gc()

merged <- merged %>% 
  NormalizeData(assay = 'RNA') %>%
  SCTransform(assay = "RNA", verbose = TRUE, conserve.memory = T) %>%
  RunPCA(assay = "SCT", verbose = TRUE) %>%
  #RunHarmony("Slice", assay.use = "SCT", plot_convergence = TRUE, max.iter.harmony = 15) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(verbose = TRUE, resolution = 0.8) %>%
  RunUMAP(reduction = "pca",
          dims = 1:30)


DimPlot(merged, group.by = 'Age', label = T, cols = 'Paired') + 
  DimPlot(merged, label = T) 
ggsave(paste0('DimPlot_clusters', add_filename,'.pdf'), width = 10, height = 5)
ggsave(paste0('DimPlot_clusters', add_filename,'.png'), width = 10, height = 5)


p <- DotPlot(merged, features =markers.from.email$Generic, group.by = 'Age', assay = 'RNA', scale = FALSE ) 
p <- p  + RotatedAxis()
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + paletteer::scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0(out_path,'Dotplot_generic_markers_', add_filename,'.pdf'), useDingbats = F, width = 5, height = 3, limitsize = F)
ggsave(paste0(out_path,'Dotplot_generic_markers_', add_filename,'.png'), width = 5, height = 3, limitsize = F)


p <- DotPlot(merged, features =markers.from.email$SASP, group.by = 'Age', assay = 'RNA', scale = FALSE ) 
p <- p  + RotatedAxis()
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + paletteer::scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0(out_path,'Dotplot_SASP_markers_', add_filename,'.pdf'), useDingbats = F, width = 5, height = 3, limitsize = F)
ggsave(paste0(out_path,'Dotplot_SASP_markers_', add_filename,'.png'),  width = 5, height = 3, limitsize = F)


p <- DotPlot(merged, features =c('IL1A', 'IL6', 'CXCL8', 'CCL2'), group.by = 'Age', assay = 'RNA', scale = FALSE ) 
p <- p  + RotatedAxis()
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + paletteer::scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0(out_path,'Dotplot_SASP_chemokines_markers_', add_filename,'.pdf'), useDingbats = F, width = 5, height = 3, limitsize = F)
ggsave(paste0(out_path,'Dotplot_SASP_chemokines_markers_', add_filename,'.png'), width = 5, height = 3, limitsize = F)


FeaturePlot(merged, features = 'MKI67', order = T)

VlnPlot(merged, markers.from.email$Generic, group.by = 'Age', ncol = 4, assay = 'RNA')
ggsave(paste0(out_path,'VlnPlot_Generic_markers_', add_filename,'.pdf'), useDingbats = F, width = 10, height = 4, limitsize = F)

VlnPlot(merged, markers.from.email$SASP, group.by = 'Age', ncol = 8, assay = 'RNA')
ggsave(paste0(out_path,'VlnPlot_SASP_markers_', add_filename,'.pdf'), useDingbats = F, width = 19, height = 4, limitsize = F)


markers.from.email$Generic %>% walk(function(gene) {
  FeaturePlot(merged, features = gene, order = T, ncol = 1, pt.size = 0.5) +
    scale_color_viridis_c(option = "magma", direction = 1)
ggsave(paste0(out_path,'Featureplot_generic_markers_',gene,'_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4)
ggsave(paste0(out_path,'Featureplot_generic_markers_',gene,'_', add_filename,'.png'), width = 4.5, height = 4)})


markers.from.email$SASP %>% walk(function(gene) {
  FeaturePlot(merged, features = gene, order = T, ncol = 1, pt.size = 0.5) +
    scale_color_viridis_c(option = "magma", direction = 1)
ggsave(paste0(out_path,'Featureplot_SASP_markers_genes_',gene,'_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4) 
ggsave(paste0(out_path,'Featureplot_SASP_markers_genes_',gene,'_', add_filename,'.png'), width = 4.5, height = 4) })


markers <- FindMarkers(object = merged, group.by = 'Age', assay = 'SCT', test.use = 'LR', ident.1 = 'Old')
fwrite(markers, paste0(out_path,'Markers_old_vs_young_LR_', add_filename,'.txt') , row.names = T, sep = '\t')


toplot <- arrange(markers, avg_log2FC) %>% filter(abs(avg_log2FC) >= 0.75)

toplot <- arrange(markers, avg_log2FC) %>% filter(abs(avg_log2FC) >= 0.75)

p <- DotPlot(merged, features =rownames(toplot), group.by = 'Age', assay = 'RNA', scale = FALSE ) 
p <- p  + RotatedAxis()
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + paletteer::scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0(out_path,'Dotplot_SASP_DEGs_old_vs_young_', add_filename,'.pdf'), useDingbats = F, width = 10, height = 3, limitsize = F)
ggsave(paste0(out_path,'Dotplot_SASP_DEGs_old_vs_young', add_filename,'.png'), width = 10, height = 3, limitsize = F)


FeaturePlot(merged, features = 'MIF', order = T, ncol = 1, pt.size = 0.5) +
  scale_color_viridis_c(option = "magma", direction = 1)

cat('calculate senMayo modules\n')
merged <- AddModuleScore(
  object = merged,
  assay = 'RNA',
  features = list(SenMayo = senmayo$Gene),
  ctrl = 15,
  name = 'SenMayo'
)
cat('done\n')

merged@meta.data
FeaturePlot(merged, features = 'SenMayo1', order = T, min.cutoff = 0, pt.size = 0.5) +
  scale_color_viridis_c(option = "magma", direction = 1)
ggsave(paste0(out_path,'Featureplot_SenMayo_geneset_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4)

VlnPlot(merged, 'SenMayo1', group.by = 'Age', cols = brewer.pal(n = 9, 'Dark2')) +
  geom_boxplot(outlier.shape = NA, inherit.aes = FALSE,
               mapping = aes_string(x = 'ident', y = 'SenMayo1', fill = 'ident'),
               alpha = 0.5,
               color = "black") +
  ggpubr::stat_compare_means(label = 'p.signif')
ggsave(paste0(out_path,'VlnPlot_SenMayo_scores_', add_filename,'.pdf'), useDingbats = F, width = 3, height = 6, limitsize = F)

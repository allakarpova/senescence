# find signatures of senescence per cluster
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)
library(viridis)
library(RColorBrewer)


wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/cluster7_signature'
dir.create(wd, showWarnings = F, recursive = T)
setwd(wd)

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res.7.mCAF_htan.brca.sc.v10.RDS')
panc <- FindClusters(object = panc, resolution = 2)
pdf('../Dimplot_clusters_res2.pdf', width = 5, height = 5, useDingbats =F )
print(DimPlot(panc, label = T)&coord_fixed())
dev.off()
cell.type = 'mCAF'
DefaultAssay(panc) <- 'RNA'

# find DEGs
panc$for_test <- case_when (panc$seurat_clusters == '7' ~ 'cluster 7',
                            panc$seurat_clusters == '3' ~ 'cluster 3',
                            TRUE ~ 'the rest')

panc$for_test <- factor (panc$for_test , levels = c('cluster 7', 'cluster 3','the rest' ))
Idents(panc) <- 'for_test'
lr.deg <- FindMarkers(panc, test.use = 'LR',  ident.1 = 'cluster 7', ident.2 = 'the rest')
fwrite(lr.deg, 'DEGs_cluster7_vs_other_no3.txt', row.names = T, sep = '\t')

lr.deg3 <- FindMarkers(panc, test.use = 'LR',  ident.1 = 'cluster 3', ident.2 = 'the rest')
fwrite(lr.deg3, 'DEGs_cluster3_vs_other_no7.txt', row.names = T, sep = '\t')

cluster3.deg <- rownames(subset(lr.deg3, p_val_adj < 0.05 & avg_log2FC < 0))

cluster7.unique.deg <- subset(lr.deg, p_val_adj < 0.05 & avg_log2FC < 0) %>%
  dplyr::filter (!(rownames(.) %in% cluster3.deg))
fwrite(cluster7.unique.deg, 'DEGs_cluster7_unique_down_vs_other_no3.txt', row.names = T, sep = '\t')


cluster7.unique.deg <- fread('DEGs_cluster7_unique_up_vs_other_no3.txt', data.table = F)
DoHeatmap(panc, features = cluster7.unique.deg$V1, raster = T) +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(10, "RdBu")))(1024))
ggsave('Heatmap_cluster7_unique_up_degs_all_blue.pdf', width = 15, height = 8, useDingbats = F)


cluster7.unique.deg <- fread('DEGs_cluster7_unique_down_vs_other_no3.txt', data.table = F)
DoHeatmap(panc, features = cluster7.unique.deg$V1, raster = T) +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(10, "RdBu")))(1024))
ggsave('Heatmap_cluster7_unique_down_degs_all_blue.pdf', width = 15, height = 6, useDingbats = F)


cluster7.top.deg <- subset(lr.deg, p_val_adj < 0.05 & avg_log2FC < 0)

DoHeatmap(panc, features = rownames(cluster7.top.deg)[1:50], raster = F, ) +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(10, "RdBu")))(1024))
ggsave('Heatmap_cluster7_top_down_degs.pdf', width = 15, height = 8, useDingbats = F)


subset(lr.deg, p_val_adj < 0.05 & avg_log2FC > 0) %>% subset (rownames(.) %in% toupper(sasp.mouse))
cluster7.unique.deg %>% subset (V1 %in% toupper(sasp.mouse))

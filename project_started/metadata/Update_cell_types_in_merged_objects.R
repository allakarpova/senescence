suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readxl))

setwd('/Users/allakarpova/lab_Ding/work/Sennet/analysis/merged/merge_liver_combo_1')
panc <- readRDS('2_Merged_normalized_liver_combo.rds')
meta <- fread('2_Merged_normalized_liver_combo.metadata.tsv')

dir.create('upd.cell.type')
setwd('upd.cell.type')

DefaultAssay(panc) = 'SCT'
colnames(panc@meta.data)

tb <- table(panc$seurat_clusters, panc$cell_type)
cluster.match.celltype <- apply(tb, 1, function(x) {
  colnames(tb)[which.max (x)]
})
panc$cell_type <- cluster.match.celltype[as.character(panc$seurat_clusters)]
panc$cell_type[panc$seurat_clusters==30] <- 'Pericytes'
panc$cell_type[panc$seurat_clusters==27] <- 'Inflammatory macs'

DimPlot(panc, group.by = 'cell_type', reduction = 'wnn.umap', label = TRUE)
ggsave('Dimplot_cell_type_updated.pdf')
FeaturePlot(panc, 'RGS5' ,order = TRUE, reduction = 'wnn.umap')

panc$Age.group <- ifelse(panc$Age<50, 'Young', 'Old')
DimPlot(panc, group.by = 'Age.group', reduction = 'wnn.umap', label = TRUE)



####### human liver single cell paper
degs <- read_excel('~/lab_Ding/work/single_cell/senescence/snATAC_combo/papers/41467_2018_6318_MOESM5_ESM_formatted.xlsx', col_names = T)

genes2plot <- degs$Gene %>% unique()
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA', cluster.idents = TRUE)
p$data <- merge(p$data, degs[c(1,3)], by.x = 'features.plot', by.y = 'Gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell_type , scales = "free",  drop = T, ncol = 3)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_viridis_c("viridis", direction = 1)

ggsave(paste0( "Dotplot_20human_liver_paper_markers_gene_expression_RNA.pdf"),plot = p, height=50,width=20,useDingbats=FALSE,limitsize = FALSE)


# plot dotplot markers MINE
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v04192022.txt', data.table = F, header = T)
genes2plot <- myeloid.genes$Gene %>% unique()

p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA', cluster.idents = TRUE)

p$data <- merge(p$data, myeloid.genes[1:3], by.x = 'features.plot', by.y = 'Gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~Gene_set_group + Gene_set , scales = "free",  drop = T, ncol = 9)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_viridis_c("viridis", direction = 1)

ggsave(paste0( "Dotplot_marker_gene_expression_RNA.pdf"),height=80,width=50,useDingbats=FALSE,limitsize = FALSE)

fwrite(panc@meta.data, paste0('2_Merged_normalized_liver_combo.metadata_05232022.tsv'), row.names = T, sep = '\t')









library(Seurat)
library(dplyr)
library(umap)
library(ggplot2)
library(cowplot)
library(readxl)
library(data.table)
library(paletteer)



panc1 = readRDS('~/lab_Ding/work/single_cell/senescence/objects/Second_DOXO/Second_DOXO_processed.rds')
panc2 = readRDS('~/lab_Ding/work/single_cell/senescence/objects/Second_Old/Second_Old_processed.rds')
panc3 = readRDS('~/lab_Ding/work/single_cell/senescence/objects/Second_Young/Second_Young_processed.rds')


int_st <- merge(panc1, y = c(panc2, panc3) , add.cell.ids = c('Second_DOXO', 'Second_Old', 'Second_Young'), project = "Second_round_mouse")
rm(panc1, panc2, panc3)

int_st <- int_st %>%
  SCTransform(
    assay = 'RNA',
    vars.to.regress = c("nCount_RNA", "percent.mito"),
    #variable.features.n = opt$top_variable_features,
    return.only.var.genes = FALSE
  ) %>%
  RunPCA(assay = 'SCT', do.print = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = "pca") %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = 1)

int_st <- NormalizeData(int_st, assay = 'RNA')
int_st <- ScaleData(int_st, assay = 'RNA')

saveRDS(int_st, '~/lab_Ding/work/single_cell/senescence/objects/Secound_mouse_integrated.rds')

pdf('~/lab_Ding/work/single_cell/senescence/mouse_BM/cell_type_annotation/Dimplots_second.pdf', useDingbats = F, width = 6, height = 5)
print(DimPlot(int_st, group.by = 'orig.ident')&coord_fixed())
print(DimPlot(int_st, label = T)&coord_fixed())
dev.off()

# process marker genes from https://www.sciencedirect.com/science/article/pii/S0092867419304593?via%3Dihub#mmc1
BM.stroma.markers <- read_excel('~/lab_Ding/work/single_cell/senescence/Papers/Mouse_BM_stroma/1-s2.0-S0092867419304593-mmc1.xlsx', skip = 1)
BM.stroma.markers <- BM.stroma.markers %>% mutate (Cell_type = case_when(BM.stroma.markers$cluster==0 ~ 'EC-sinusoidal',
                                                                         BM.stroma.markers$cluster==1 ~ 'Lepr-MSC',
                                                                         BM.stroma.markers$cluster==2 ~ 'Chondro-hyper',
                                                                         BM.stroma.markers$cluster==3 ~ 'Fibro-4',
                                                                         BM.stroma.markers$cluster==4 ~ 'Chondro-progen',
                                                                         BM.stroma.markers$cluster==5 ~ 'Fibro-5',
                                                                         BM.stroma.markers$cluster==6 ~ 'EC-arteriolar',
                                                                         BM.stroma.markers$cluster==7 ~ 'OLC-1',
                                                                         BM.stroma.markers$cluster==8 ~ 'OLC-2',
                                                                         BM.stroma.markers$cluster==9 ~ 'Fibro-1',
                                                                         BM.stroma.markers$cluster==10 ~ 'Chondro-prol.rest',
                                                                         BM.stroma.markers$cluster==11 ~ 'EC-arterial',
                                                                         BM.stroma.markers$cluster==12 ~ 'Pericytes',
                                                                         BM.stroma.markers$cluster==13 ~ 'Chondro',
                                                                         BM.stroma.markers$cluster==15 ~ 'Fibro-2',
                                                                         BM.stroma.markers$cluster==16 ~ 'Fibro-3',
                                                                         BM.stroma.markers$cluster==17 ~ 'Chondro-prehyper-2'))

top10.stroma.markers <- BM.stroma.markers %>% filter(p_val_adj < 0.05) %>% group_by(Cell_type) %>% top_n (n = 20, wt = avg_logFC)


#process markers from https://www.nature.com/articles/s41556-019-0439-6#Sec27
BM.cells.markers <-  read_excel('~/lab_Ding/work/single_cell/senescence/Papers/Mouse_BM_single_cell/41556_2019_439_MOESM3_ESM.xlsx', sheet = 2)[-1,]
BM.cells.markers <- BM.cells.markers[1:20,]
top20.cell.markers <- reshape2::melt(as.matrix(BM.cells.markers))[2:3] %>% data.frame() 
colnames(top20.cell.markers) <- c('Cell_type', 'gene')


wd <- '~/lab_Ding/work/single_cell/senescence/mouse_BM/cell_type_annotation/integrated'
dir.create(wd, recursive = T)
setwd(wd)
add_filename <- 'integrated_second'
#####
genes2plot <- top10.stroma.markers$gene %>% unique()
p <- DotPlot(object = int_st, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'SCT')

p$data <- merge(p$data, top10.stroma.markers[7:8], by.x = 'features.plot', by.y = 'gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell_type , scales = "free",  drop = T, ncol = 6)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0( "Dotplot_20stroma_gene_expression_", add_filename, "_SCT.pdf"),height=20,width=40,useDingbats=FALSE,limitsize = FALSE)



#####
genes2plot <- top20.cell.markers$gene %>% unique()
p <- DotPlot(object = int_st, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'SCT')

p$data <- merge(p$data, top20.cell.markers, by.x = 'features.plot', by.y = 'gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell_type , scales = "free",  drop = T, ncol = 6)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0( "Dotplot_20cell.type_gene_expression_", add_filename, "_SCT.pdf"),height=35,width=40,useDingbats=FALSE,limitsize = FALSE)

markers <- FindAllMarkers(int_st, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25, assay = 'SCT')
fwrite(markers, 'Markers.cluster.integrated_second.txt', sep ='\t')

# update cell type annotation
int_st <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Secound_mouse_integrated.rds')
cell.types <- fread('Cell_clusters.txt', data.table = F)
metadata <- int_st@meta.data %>% mutate(Row.names = rownames(.))
metadata.merged <- merge(metadata, cell.types, by.x = 'seurat_clusters', by.y = 'Cluster', all.x = T )
rownames(metadata.merged) <- metadata.merged$Row.names
int_st <- AddMetaData(int_st, metadata.merged[-10])
getwd()
DimPlot(int_st, group.by = 'Cell_type', label = T)
ggsave('Dimplot_cell_type.pdf', width = 10, height = 7, useDingbats = F)
saveRDS(int_st, '~/lab_Ding/work/single_cell/senescence/objects/Secound_mouse_integrated.rds')



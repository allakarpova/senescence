#senescence for fetal fibroblast cells 
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

#obj.path <- '~/lab_Ding/work/Sennet/rds_objects/20211214-WS1_PD45G1Z1_1Bn1_1/20211214-WS1_PD45G1Z1_1Bn1_1_processed.rds'
obj.path <- '~/lab_Ding/work/Sennet/rds_objects/SN001H1-Md1A3Y1K2N1/SN001H1-Md1A3Y1K2N1_processed_multiomic.rds'

add_filename <- 'SN001H1-Md1A3Y1K2N1'
out_path <- paste0('~/lab_Ding/work/Sennet/analysis/', add_filename, '/')
dir.create(out_path, showWarnings = F, recursive = T)


cat('opening object...\n')
panc <- readRDS(obj.path)
cat('done\n')
# 
meta <- fread(paste0('~/lab_Ding/work/Sennet/cell_typing/manual/',add_filename,'/',add_filename,'_processed_multiomic_cellTyped.meta.data')) %>% 
   data.frame(row.names = 1) 

panc <- AddMetaData(panc, meta)

dim(panc)

sum(panc@assays$RNA@counts['CDKN2A', ] > 0)

DefaultAssay(panc) <- 'SCT'

red.touse <- ifelse("wnn.umap" %in% Reductions(panc), "wnn.umap","umap")

DimPlot(object = panc, group.by = 'predicted_doublet',cols = c('black', 'yellow'), 
        reduction = red.touse,label=TRUE,label.size=6)

DimPlot(panc, label = T)


markers.from.email$Generic %>% walk(function(gene) {
  FeaturePlot(panc, features = gene, order = T, ncol = 1, pt.size = 0.5, reduction = red.touse) +
  scale_color_viridis_c(option = "magma", direction = 1)
  ggsave(paste0(out_path,'Featureplot_generic_markers_',gene,'_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4)
  ggsave(paste0(out_path,'Featureplot_generic_markers_',gene,'_', add_filename,'.png'), width = 4.5, height = 4)})

markers.from.email$SASP %>% walk(function(gene) {
  FeaturePlot(panc, features = gene, order = T, ncol = 1, pt.size = 0.5, reduction = red.touse) +
  scale_color_viridis_c(option = "magma", direction = 1)
  ggsave(paste0(out_path,'Featureplot_SASP_markers_genes_',gene,'_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4) 
  ggsave(paste0(out_path,'Featureplot_SASP_markers_genes_',gene,'_', add_filename,'.png'),  width = 4.5, height = 4) })


cat('prepare gene lists\n')
cpdb.split.df <- fread('~/R_working_dir/data/ConsensusPathDB/CPDB_pathways_genes.2column.txt', data.table = F)
cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

db.sen <- sapply(split(cpdb.split.df.sen, cpdb.split.df.sen$ont), function(x) {
  y <- as.character(x$gene)
  return(y)
})

total.sen <- db.sen
names(total.sen) <- make.names(names(total.sen))
names(total.sen) <- gsub('_', '-', names(total.sen))

cat('calculate modules\n')
panc <- AddModuleScore(
  object = panc,
  assay = 'RNA',
  features = total.sen,
  ctrl = 15,
  name = 'senescence'
)
cat('done\n')

n.list <- length(total.sen)
colnames(panc@meta.data)[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)] <- names(total.sen)


#Add senescence module as an assay
cat('add modules in a separate assay\n')
panc[['senescence.module']] <- CreateAssayObject(data = t(x = FetchData(object = panc, vars = names(total.sen)) ))

FeaturePlot(panc, features = names(total.sen)[grepl('SENESCE', names(total.sen))], order = T, min.cutoff = 0, pt.size = 0.5, reduction = "wnn.umap")
ggsave(paste0(out_path,'Featureplot_all_senescence_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)

FeaturePlot(panc, features = names(total.sen)[grepl('DAMAGE', names(total.sen))], order = T, min.cutoff = 0, pt.size = 0.5, reduction = "wnn.umap")
ggsave(paste0(out_path,'Featureplot_all_DAMAGE_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 30, limitsize = F)

FeaturePlot(panc, features = names(total.sen)[grepl('INFLAMMAT', names(total.sen))], order = T, min.cutoff = 0, pt.size = 0.5, reduction = "wnn.umap")
ggsave(paste0(out_path,'Featureplot_all_INFLAMMAT_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)

dim(panc)

cat('calculate senMayo modules\n')
panc <- AddModuleScore(
  object = panc,
  assay = 'RNA',
  features = list(SenMayo = senmayo$Gene),
  ctrl = 15,
  name = 'SenMayo'
)
cat('done\n')

panc@meta.data
FeaturePlot(panc, features = 'SenMayo1', order = T, min.cutoff = 0, pt.size = 0.5, reduction = "wnn.umap") +
  scale_color_viridis_c(option = "magma", direction = 1)
ggsave(paste0(out_path,'Featureplot_SenMayo_geneset_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4)

panc$Sample <- 'All cells'
VlnPlot(panc, markers.from.email$Generic, group.by = 'cell_type', ncol = 4, assay = 'RNA')
ggsave(paste0(out_path,'VlnPlot_Generic_markers_', add_filename,'.pdf'), useDingbats = F, width = 15, height = 4, limitsize = F)

VlnPlot(panc, markers.from.email$SASP, group.by = 'cell_type', ncol = 7, assay = 'RNA')
ggsave(paste0(out_path,'VlnPlot_SASP_markers_', add_filename,'.pdf'), useDingbats = F, width = 25, height = 4, limitsize = F)


p <- DotPlot(panc, features =markers.from.email$Generic, group.by = 'cell_type', assay = 'SCT', scale = FALSE, scale.max = 100, scale.min = 0 , cluster.idents = TRUE ) 
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
ggsave(paste0(out_path,'Dotplot_generic_markers_', add_filename,'.pdf'), useDingbats = F, width = 6, height = 5, limitsize = F)
ggsave(paste0(out_path,'Dotplot_generic_markers_', add_filename,'.png'),  width = 6, height = 5, limitsize = F)


p <- DotPlot(panc, features =markers.from.email$SASP, group.by = 'cell_type', assay = 'SCT', scale = FALSE, 
             scale.max = 100, scale.min = 0, cluster.idents = TRUE) 
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
ggsave(paste0(out_path,'Dotplot_SASP_markers_', add_filename,'.pdf'), useDingbats = F, width = 6, height = 5, limitsize = F)
ggsave(paste0(out_path,'Dotplot_SASP_markers_', add_filename,'.png'), width = 6, height = 5, limitsize = F)


p <- DotPlot(panc, features = 'CDKN2A', group.by = 'cell_type', assay = 'SCT', scale = FALSE,
              cluster.idents = FALSE ) 
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
ggsave(paste0(out_path,'Dotplot_CDKN2A_', add_filename,'.pdf'), useDingbats = F, width = 5, height = 5, limitsize = F)
ggsave(paste0(out_path,'Dotplot_CDKN2A_', add_filename,'.png'),  width = 5, height = 5, limitsize = F)
*

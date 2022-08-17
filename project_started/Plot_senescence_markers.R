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
obj.path <- '~/lab_Ding/work/Sennet/rds_objects/fixed/20211110-WS1_PD102N1Z1_1Bn1_1/20211110-WS1_PD102N1Z1_1Bn1_1_processed.rds'

out_path <- '~/lab_Ding/work/Sennet/analysis/Test_cell_line_old_snRNA_fixed/'
dir.create(out_path, showWarnings = F, recursive = T)
add_filename <- 'Test_cell_line_old_snRNA'


cat('opening object...\n')
panc <- readRDS(obj.path)
cat('done\n')


dim(panc)

sum(panc@assays$RNA@counts['CDKN2A', ] > 0)

DefaultAssay(panc) <- 'SCT'

DimPlot(object = panc, group.by = 'predicted_doublet',cols = c('black', 'yellow'), 
        reduction = "umap",label=TRUE,label.size=6)

DimPlot(panc, label = T)

markers.from.email$Generic %>% walk(function(gene) {
  FeaturePlot(panc, features = gene, order = T, ncol = 1, pt.size = 0.5) +
  scale_color_viridis_c(option = "magma", direction = 1)
  ggsave(paste0(out_path,'Featureplot_generic_markers_',gene,'_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4)
  ggsave(paste0(out_path,'Featureplot_generic_markers_',gene,'_', add_filename,'.png'), width = 4.5, height = 4)})

markers.from.email$SASP %>% walk(function(gene) {
  FeaturePlot(panc, features = gene, order = T, ncol = 1, pt.size = 0.5) +
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

FeaturePlot(panc, features = names(total.sen)[grepl('SENESCE', names(total.sen))], order = T, min.cutoff = 0, pt.size = 0.5)
ggsave(paste0(out_path,'Featureplot_all_senescence_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)

FeaturePlot(panc, features = names(total.sen)[grepl('DAMAGE', names(total.sen))], order = T, min.cutoff = 0, pt.size = 0.5)
ggsave(paste0(out_path,'Featureplot_all_DAMAGE_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 30, limitsize = F)

FeaturePlot(panc, features = names(total.sen)[grepl('INFLAMMAT', names(total.sen))], order = T, min.cutoff = 0, pt.size = 0.5)
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
FeaturePlot(panc, features = 'SenMayo1', order = T, min.cutoff = 0, pt.size = 0.5) +
  scale_color_viridis_c(option = "magma", direction = 1)
ggsave(paste0(out_path,'Featureplot_SenMayo_geneset_', add_filename,'.pdf'), useDingbats = F, width = 4.5, height = 4)

panc$Sample <- 'All cells'
VlnPlot(panc, markers.from.email$Generic, group.by = 'Sample', ncol = 5, assay = 'RNA')
ggsave(paste0(out_path,'VlnPlot_Generic_markers_', add_filename,'.pdf'), useDingbats = F, width = 10, height = 6, limitsize = F)

VlnPlot(panc, markers.from.email$SASP, group.by = 'Sample', ncol = 8, assay = 'RNA')
ggsave(paste0(out_path,'VlnPlot_SASP_markers_', add_filename,'.pdf'), useDingbats = F, width = 16, height = 6, limitsize = F)


p <- DotPlot(panc, features =markers.from.email$Generic, group.by = 'seurat_clusters', assay = 'SCT', scale = FALSE, scale.max = 100, scale.min = 0  ) 
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
ggsave(paste0(out_path,'Dotplot_generic_markers_', add_filename,'.pdf'), useDingbats = F, width = 5, height = 5, limitsize = F)
ggsave(paste0(out_path,'Dotplot_generic_markers_', add_filename,'.png'),  width = 5, height = 5, limitsize = F)


p <- DotPlot(panc, features =markers.from.email$SASP, group.by = 'seurat_clusters', assay = 'SCT', scale = FALSE, scale.max = 100, scale.min = 0 ) 
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
ggsave(paste0(out_path,'Dotplot_SASP_markers_', add_filename,'.pdf'), useDingbats = F, width = 5, height = 5, limitsize = F)
ggsave(paste0(out_path,'Dotplot_SASP_markers_', add_filename,'.png'), width = 5, height = 5, limitsize = F)


FeaturePlot(panc, features = c('CDKN1A','CDKN2A'), blend = T, order = T, min.cutoff = 0, pt.size = 0.5, 
            #cols = c('grey', '#ffff33', '#377eb8'), 
            blend.threshold = 0.1)
ggsave(paste0(out_path,'Featureplot_CDKN1A_CDKN2A_blend_', add_filename,'.pdf'), useDingbats = F, width = 12, height = 4, limitsize = F)

FeaturePlot(panc, features = c("nCount_RNA","percent.mt"), order = T, min.cutoff = 0, pt.size = 0.5)

DefaultAssay(panc) <- 'RNA'
FeatureScatter(panc, feature1 = 'CDKN1A', feature2 = 'TP53', group.by = 'Sample', slot = 'data') +
  geom_density2d()

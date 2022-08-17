library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)
library(paletteer)
library(viridis)

doSenescence <- function(input.path, out_path, add_filename) {
  cat('opening object...\n')
  panc <- readRDS(input.path)
  cat('done\n')
  
  fwrite(data.frame(Embeddings(panc, reduction = 'umap')), paste0(out_path,'Integrated_object_embeddings_', add_filename,'.txt'), sep = '\t', row.names = T)
  
  
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
  #Add senescence module as an assay
  cat('add modules in a separate assay\n')
  panc[['senescence.module']] <- CreateAssayObject(data = t(x = FetchData(object = panc, vars = names(total.sen)) ))
  
  
  
  fwrite(panc@meta.data[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)], 
         paste0(out_path,'Senescense_scores_', add_filename,'.txt'), sep = '\t', row.names = T)
  
  fwrite(panc@meta.data, 
         paste0(out_path,'Metadata_', add_filename,'.txt'), sep = '\t', row.names = T)
  
  cat('make plots\n')
  DimPlot(object = panc, group.by = 'cell_type_specific', label = T)
  ggsave(paste0(out_path,'Dimplot_cell_types_', add_filename,'.pdf'), useDingbats = F, width = 12, height = 7.5)
  
  DimPlot(object = panc, label = T)
  ggsave(paste0(out_path,'Dimplot_clusters_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 7.5)
  
  DimPlot(object = panc, group.by = 'tissue_type', label = T)
  ggsave(paste0(out_path,'Dimplot_tissue_type_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 7.5)
  
  FeaturePlot(panc, features = c('Generic', 'SASP'), order = T, min.cutoff = 0)
  ggsave(paste0(out_path,'Featureplot_Zhang_markers_module_', add_filename,'.pdf'), useDingbats = F, width = 8, height = 4)
  
  FeaturePlot(panc, features = names(total.sen)[grepl('DAMAG', names(total.sen))], order = T, min.cutoff = 0)
  ggsave(paste0(out_path,'Featureplot_all_damage_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 40, limitsize = F)
  
  FeaturePlot(panc, features = names(total.sen)[grepl('SENESCE', names(total.sen))], order = T, min.cutoff = 0)
  ggsave(paste0(out_path,'Featureplot_all_senescence_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)
  
  FeaturePlot(panc, features = names(total.sen)[grepl('INFLAMMAT', names(total.sen))], order = T, min.cutoff = 0)
  ggsave(paste0(out_path,'Featureplot_all_inflammatory_modules_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 20, limitsize = F)
  
  FeaturePlot(panc, features = markers.from.email$Generic, order = T, ncol = 3)
  ggsave(paste0(out_path,'Featureplot_Zhang_generic_markers_genes_', add_filename,'.pdf'), useDingbats = F, width = 12, height = 8)
  
  FeaturePlot(panc, features = markers.from.email$SASP, order = T, ncol = 3)
  ggsave(paste0(out_path,'Featureplot_Zhang_SASP_markers_genes_', add_filename,'.pdf'), useDingbats = F, width = 12, height = 12)
  
  FeaturePlot(panc, features = 'LMNB1', order = T)
  ggsave(paste0(out_path,'Featureplot_LMNB1_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 6)
  FeaturePlot(panc, features = 'CGAS', order = T)
  ggsave(paste0(out_path,'Featureplot_CGAS_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 6)
  #FeaturePlot(panc, features = 'STING1', order = T)
  #ggsave(paste0(out_path,'Featureplot_STING1_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 6)
  
  
  DoHeatmap(object = panc, features = names(total.sen), assay = 'senescence.module', slot = 'data', group.by = 'cell_type_specific', label=FALSE)
  ggsave(paste0(out_path,'Heatmap_by_cell_type_scored_senesence_DNA-damage.paths_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 8)
  
  DotPlot(object = panc, features = markers.from.email$SASP, assay = 'RNA', group.by = 'cell_type_specific')
  ggsave(paste0(out_path,'Dotplot_Zhang_SASP_markers_genes_cell_type_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 15)
  
  DotPlot(object = panc, features = markers.from.email$SASP)
  ggsave(paste0(out_path,'Dotplot_Zhang_SASP_markers_genes_cluster_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 15)
  
  DotPlot(object = panc, features = markers.from.email$Generic, assay = 'RNA', group.by = 'cell_type_specific')
  ggsave(paste0(out_path,'Dotplot_Zhang_Generic_markers_genes_cell_type_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 15)
  
  DotPlot(object = panc, features = markers.from.email$Generic)
  ggsave(paste0(out_path,'Dotplot_Zhang_Generic_markers_genes_cluster_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 15)
  
}


cat('prepare gene lists\n')
cpdb.split.df <- fread('/diskmnt/Projects/Users/allakarpova/Projects/senescence/CPDB_pathways_genes.2column.txt', data.table = F)
cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), SASP = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'))
db.sen <- split(cpdb.split.df.sen$gene, cpdb.split.df.sen$ont)
total.sen <- c(db.sen, markers.from.email)
names(total.sen) <- make.names(names(total.sen))
names(total.sen) <- gsub('_', '-', names(total.sen))

input.path <- '/diskmnt/Projects/Users/allakarpova/Projects/senescence/brca/non.tumor.cells/Non.tumor.object.res1.htan.brca.sn.v14.RDS'
out_path <- '/diskmnt/Projects/Users/allakarpova/Projects/senescence/brca/non.tumor.cells/Plots/'
dir.create(out_path, showWarnings = F)
add_filename <- 'non_tumor_htan_brca_sn_v14'
doSenescence(input.path, out_path, add_filename)


input.path <- '/diskmnt/Projects/Users/allakarpova/Projects/senescence/brca/non.tumor.cells/Non.tumor.object.res1.htan.brca.sc.v10.RDS'
out_path <- '/diskmnt/Projects/Users/allakarpova/Projects/senescence/brca/non.tumor.cells/Plots/'
dir.create(out_path, showWarnings = F)
add_filename <- 'non_tumor_htan_brca_sc_v10'

doSenescence(input.path, out_path, add_filename)






#senescence for GBM non tumor cells 
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)
library(paletteer)
library(viridis)


cat('prepare gene lists\n')
cpdb.split.df <- fread('/diskmnt/Projects/Users/allakarpova/Projects/senescence/CPDB_pathways_genes.2column.txt', data.table = F)
cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1'), 
                           SASP = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7'))
db.sen <- sapply(split(cpdb.split.df.sen, cpdb.split.df.sen$ont), function(x) {
  y <- as.character(x$gene)
  return(y)
})

total.sen <- c(db.sen, markers.from.email)
names(total.sen) <- make.names(names(total.sen))
names(total.sen) <- gsub('_', '-', names(total.sen))

input.path <- '/diskmnt/Projects/Users/allakarpova/Projects/senescence/gbm/non.tumor.cells/Non.tumor.object.res1.gbm.v2.0.RDS'
out_path <- '/diskmnt/Projects/Users/allakarpova/Projects/senescence/gbm/non.tumor.cells/Plots/'
dir.create(out_path, showWarnings = F, recursive = T)
add_filename <- 'non.tumor_gbm_28_sample'

# input.path <- '~/lab_Ding/work/single_cell/kidney/ccRCC/objects/integration.202002012.v3.immune.lymphoid_reclustered20200707.RDS'
# out_path <- '/diskmnt/Projects/cptac_scratch_4/CPTAC3_GBM/snRNA_seq/snRNA_analysis/'
# add_filename <- 'all_genes_no_doublets_18sample'

cat('opening object...\n')
panc <- readRDS(input.path)
cat('done\n')

DefaultAssay(panc) <- 'RNA'

my.metadata <- fread('/diskmnt/Projects/cptac_scratch_4/CPTAC3_GBM/single_cell_data_freeze/v2/snRNA/snRNA_metadata.v2.1.tsv.gz', data.table = F) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)
panc <- AddMetaData(panc, metadata = my.metadata)

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

fwrite(data.frame(Embeddings(panc, reduction = 'umap')), paste0(out_path,'Integrated_object_embeddings_', add_filename,'.txt'), sep = '\t', row.names = T)
fwrite(panc@meta.data[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)], 
       paste0(out_path,'Senescense_scores_', add_filename,'.txt'), sep = '\t', row.names = T)

cat('make plots\n')
DimPlot(object = panc, group.by = 'cell_type_v20210128')
ggsave(paste0(out_path,'Dimplot_cell_types_', add_filename,'.pdf'), useDingbats = F, width = 12, height = 7.5)

DimPlot(object = panc, group.by = 'case_id')
ggsave(paste0(out_path,'Dimplot_case_id_', add_filename,'.pdf'), useDingbats = F, width = 12, height = 7.5)

DimPlot(object = panc, label = T)
ggsave(paste0(out_path,'Dimplot_clusters_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 7.5)



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


DoHeatmap(object = panc, features = names(total.sen), assay = 'senescence.module', slot = 'data', group.by = 'cell_type_v20210128', label=FALSE)
ggsave(paste0(out_path,'Heatmap_by_cell_type_scored_senesence_DNA-damage.paths_', add_filename,'.pdf'), useDingbats = F, width = 30, height = 8)

DotPlot(object = panc, features = markers.from.email$SASP, assay = 'RNA', group.by = 'cell_type_v20210128')
ggsave(paste0(out_path,'Dotplot_Zhang_SASP_markers_genes_cell_type_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 8)

DotPlot(object = panc, features = markers.from.email$SASP)
ggsave(paste0(out_path,'Dotplot_Zhang_SASP_markers_genes_cluster_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 12)

DotPlot(object = panc, features = markers.from.email$Generic, assay = 'RNA', group.by = 'cell_type_v20210128')
ggsave(paste0(out_path,'Dotplot_Zhang_Generic_markers_genes_cell_type_', add_filename,'.pdf'), useDingbats = F, width = 9, height = 8)

DotPlot(object = panc, features = markers.from.email$Generic)
ggsave(paste0(out_path,'Dotplot_Zhang_Generic_markers_genes_cluster_', add_filename,'.pdf'), useDingbats = F, width = 7, height = 12)




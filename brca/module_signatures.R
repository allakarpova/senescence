library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)
library(ggplot2)


panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Normal.object.res.7.mCAF_htan.brca.sc.v10.RDS')
panc <- FindClusters(panc, resolution = 2)
DimPlot(panc, group.by = 'Cell_type')
ggsave('~/lab_Ding/work/single_cell/senescence/mouse_BM/senescence_correlation/Osteo/Dimplot_celltype.pdf', width = 5, height = 4.5, useDingbats = F)
DimPlot(panc, group.by = 'orig.ident', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('~/lab_Ding/work/single_cell/senescence/mouse_BM/senescence_correlation/Osteo/Dimplot_sample.pdf', width = 5, height = 4.5, useDingbats = F)
DimPlot(panc, label = T)
ggsave('~/lab_Ding/work/single_cell/senescence/mouse_BM/senescence_correlation/Osteo/Dimplot_cluster.pdf', width = 5, height = 4.5, useDingbats = F)

cell.type <- 'mCAF'
source('~/R_working_dir/scripts/senescence/markers.R')

# cpdb.split.df <- fread('~/R_working_dir/data/ConsensusPathDB/CPDB_pathways_genes.2column.txt', data.table = F)
# cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
# cpdb.split.df.sen$gene <- cpdb.split.df.sen$gene %>% tolower() %>% Caps()
# cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1') %>% tolower() %>% Caps(), 
                           Cell_cycle = c('CCND1', 'CDKN1A','CDKN2A', 'TP53') %>% tolower() %>% Caps(), 
                           DDR = c('H2AFX','H2AFY', 'TP53BP1') %>% tolower() %>% Caps(),
                           Proliferation = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1') %>% tolower() %>% Caps(),
                           SASP = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7') %>% tolower() %>% Caps())
markers.from.paper <- list (Sen.core = sen.core.mouse %>% toupper(), 
                            Sen.effector = sen.effector.mouse  %>% toupper(), 
                            sasp = sasp.mouse  %>% toupper(),
                            DDR = c('H2AFX','H2AFY', 'TP53BP1'), 
                            Proliferation = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1'))

# db.sen <- split(cpdb.split.df.sen$gene, cpdb.split.df.sen$ont)
# total.sen <- c(db.sen, markers.from.email)
# names(total.sen) <- make.names(names(total.sen))
# names(total.sen) <- gsub('_', '-', names(total.sen))
total.sen <- markers.from.paper

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
panc@meta.data

wd <- '~/lab_Ding/work/single_cell/senescence/brca/normal_samples/module_scores/'
dir.create(wd, recursive = T)
setwd(wd)
for(s in names(total.sen)) {
  #s <- 'Sen.core'
  
  VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'seurat_clusters', ncol = 1,
          assay = 'RNA')
  ggsave(paste0( 'Vlnplot_', cell.type, '_clusters_scores_', s, '.pdf'), useDingbats = F, width = 10, height = 6)
  
  FeaturePlot(panc, features = s, order = T, min.cutoff = 0)
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '.pdf'), useDingbats = F, width = 6, height = 5)
  
  FeaturePlot(panc, features = total.sen[[s]], order = T, min.cutoff = 0, ncol = 3)
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '_gene_level.pdf'), useDingbats = F, width = 10, height = length(total.sen[[s]]), limitsize = F)
  
  ggplot (panc@meta.data, aes_string (x = s, fill = 'SCT_snn_res.2')) +
    geom_density(alpha = 0.5) +
    theme_cowplot() +
    geom_vline(xintercept = 0) +
    facet_grid(SCT_snn_res.2 ~., scales = 'free') +
    ggtitle(s) 
    #scale_fill_manual(values =c( '#fc8d62','#8da0cb', '#66c2a5'))
  ggsave(paste0('Densityplot_',cell.type, '_scores_', s, '_filled.pdf'), useDingbats = F, width =5, height = 25)
}





VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
        assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5')) + ggpubr::stat_compare_means(method = 'wilcox.test',  hide.ns = T, label = "p.signif")

VlnPlot(object = panc, features = c('Generic', 'SASP'), pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
        assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave(paste0(folder, 'Vlnplot_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 18, height = 8)

FeaturePlot(panc, features = c('Generic', 'SASP'), order = T, min.cutoff = 0)
ggsave(paste0('Featureplot_',cell.type,'_Zhang_markers_module.pdf'), useDingbats = F, width = 8, height = 4)


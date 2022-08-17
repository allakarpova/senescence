library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)
library(ggplot2)

source('~/R_working_dir/scripts/senescence/markers.R')
Caps <- function(x) {
  s <- x
  toreturn = sapply(s, function(x) paste0(toupper(substring(x, 1,1)), substring(x, 2)) )
  return(as.character(toreturn))
}

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/snATAC_RNA/objects/Endothelial.rds')
panc <- NormalizeData(panc, assay = 'RNA')
cell.type <- 'Endothelial'
#panc@meta.data <- panc@meta.data[-(ncol(panc@meta.data)-2):-ncol(panc@meta.data)]
panc@meta.data$orig.ident <- panc@meta.data$dataset
panc@meta.data$Cell_type <- panc@meta.data$Cell_type_macs2

wd <- paste0('~/lab_Ding/work/single_cell/senescence/snATAC_RNA/module_signatures/', cell.type)
dir.create(wd, recursive = T)
setwd(wd)

getwd()

DimPlot(panc,label = T)
ggsave('Dimplot_cluster.pdf', width = 5.5, height = 5, useDingbats = F)

DimPlot(panc, group.by ='cell_type',  label = T)
ggsave('Dimplot_cell_type.pdf', width = 6, height = 5, useDingbats = F)

DimPlot(panc, group.by ='age',  label = T, cols = 'Paired')
ggsave('Dimplot_age.pdf', width = 6, height = 5, useDingbats = F)



# cpdb.split.df <- fread('~/R_working_dir/data/ConsensusPathDB/CPDB_pathways_genes.2column.txt', data.table = F)
# cpdb.split.df.sen <- filter(cpdb.split.df, grepl('SENESCE|DAMAGE|INFLAMMAT', ont))
# cpdb.split.df.sen$gene <- cpdb.split.df.sen$gene %>% tolower() %>% Caps()
# cpdb.split.df.sen$ont <- as.character(cpdb.split.df.sen$ont)

markers.from.email <- list(Generic = c('GLB1', 'CDKN1A','CDKN2A', 'TRP53', 'SERPINE1') %>% tolower() %>% Caps(), 
                           Cell_cycle = c('CCND1', 'CDKN1A','CDKN2A', 'TP53') %>% tolower() %>% Caps(), 
                           DDR = c('H2AFX','H2AFY', 'TP53BP1') %>% tolower() %>% Caps(),
                           Proliferation = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1') %>% tolower() %>% Caps(),
                           SASP.short = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7') %>% tolower() %>% Caps())
markers.from.paper <- list (Sen.core = sen.core.mouse, Sen.effector = sen.effector.mouse, sasp = sasp.mouse)

# db.sen <- split(cpdb.split.df.sen$gene, cpdb.split.df.sen$ont)
# total.sen <- c(db.sen, markers.from.email)
# names(total.sen) <- make.names(names(total.sen))
# names(total.sen) <- gsub('_', '-', names(total.sen))
total.sen <- c(markers.from.email, markers.from.paper)

cat('calculate modules\n')
DefaultAssay(panc) <- 'RNA'
panc <- AddModuleScore(
  object = panc,
  assay = 'RNA',
  features = total.sen,
  ctrl = 10,
  name = 'senescence'
)
cat('done\n')

n.list <- length(total.sen)
colnames(panc@meta.data)[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)] <- names(total.sen)

DefaultAssay(panc) <- 'RNA'

for(s in names(total.sen)) {
  #s <- 'sasp'
  VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'cell_type', ncol = 1,split.by = 'age',
          assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
  ggsave(paste0( 'Vlnplot_', cell.type, '_scores_', s, '.pdf'), useDingbats = F, width = 18, height = 8)
  
  VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'seurat_clusters', ncol = 1, split.by = 'age',
          assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
  ggsave(paste0( 'Vlnplot_', cell.type, '_clusters_scores_', s, '.pdf'), useDingbats = F, width = 18, height = 8)
  
  FeaturePlot(panc, features = s, order = T, min.cutoff = 0)
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '.pdf'), useDingbats = F, width = 6, height = 5)
  
  FeaturePlot(panc, features = s, order = T, min.cutoff = 0,split.by = 'age')
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '_splitted.pdf'), useDingbats = F, width = 15, height = 5.8)
  
  FeaturePlot(panc, features = total.sen[[s]], order = T, min.cutoff = 0, split.by = 'age')
  ggsave(paste0('Featureplot_',cell.type, '_scores_', s, '_gene_level_splitted.pdf'), useDingbats = F, width = 10, height = 3*length(total.sen[[s]]), limitsize = F)
  
  ggplot (panc@meta.data, aes_string (x = s, fill = 'age')) +
    geom_density(alpha = 0.5) +
    theme_cowplot() +
    facet_wrap(~cell_type, scales = 'free') +
    ggtitle(s) +
    scale_fill_manual(values =c( '#fc8d62','#8da0cb', '#66c2a5'))
  ggsave(paste0('Densityplot_',cell.type, '_scores_', s, '_filled.pdf'), useDingbats = F, width = 8, height = 5.8)
}





VlnPlot(object = panc, features = s, pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
        assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5')) + ggpubr::stat_compare_means(method = 'wilcox.test',  hide.ns = T, label = "p.signif")

VlnPlot(object = panc, features = c('Generic', 'SASP'), pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
        assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave(paste0(folder, 'Vlnplot_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 18, height = 8)

FeaturePlot(panc, features = c('Generic', 'SASP'), order = T, min.cutoff = 0)
ggsave(paste0('Featureplot_',cell.type,'_Zhang_markers_module.pdf'), useDingbats = F, width = 8, height = 4)


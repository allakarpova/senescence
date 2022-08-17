

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readxl)
library(data.table)
library(paletteer)

# nothing good
Caps <- function(x) {
  s <- x
  toreturn = sapply(s, function(x) paste0(toupper(substring(x, 1,1)), substring(x, 2)) )
  return(as.character(toreturn))
}


source('~/R_working_dir/scripts/senescence/markers.R')

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/snRNA_combo/objects/mouse_liver_by_cell_type/Hepathocytes.object.res.1.mouse_liver_old_young.rds')
cell.type <- 'Hepathocytes'
wd <- paste0('~/lab_Ding/work/single_cell/senescence/snRNA_combo/marker_expr/', cell.type)
dir.create(wd, recursive = T)
setwd(wd)



DimPlot(panc)
ggsave('Dimplot_sample.pdf',width = 9, height = 7, useDingbats = F)

VlnPlot(object = panc, features = c('CDKN1A','CDKN2A', 'TRP53', 'CCND1')%>% tolower() %>% Caps(), pt.size = 0.1, 
        group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident', assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Vlnplot_cell_cycle.pdf', width = 10, height = 15, useDingbats = F)

VlnPlot(object = panc, features = c('CDKN1A','CDKN2A', 'TRP53', 'CCND1')%>% tolower() %>% Caps(), pt.size = 0.1, 
        ncol = 1, assay = 'RNA')
ggsave('Vlnplot_cell_cycle_by_cluster.pdf', width = 6, height = 15, useDingbats = F)

VlnPlot(object = panc, features = c('Lepr', 'Cxcl12')%>% tolower() %>% Caps(), pt.size = 0.1, group.by = 'Cell_type',
        ncol = 1, assay = 'RNA')
ggsave('Vlnplot_Lepr_Cxcl12.pdf', width = 7, height = 11, useDingbats = F)



VlnPlot(object = panc, features = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1') %>% tolower() %>% Caps(), pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
        assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Vlnplot_proliferation.pdf', width = 10, height = 25, useDingbats = F)
VlnPlot(object = panc, features = c('MKI67', 'TOP2A', 'MYBL2', 'BUB1', 'PLK1')%>% tolower() %>% Caps(), pt.size = 0.1, 
        ncol = 1, assay = 'RNA')
ggsave('Vlnplot_proliferation_by_cluster.pdf', width = 7, height = 20, useDingbats = F)


VlnPlot(object = panc, features = c('H2AFX','H2AFY', 'TP53BP1') %>% tolower() %>% Caps(), pt.size = 0.9, group.by = 'Cell_type', ncol = 1,split.by = 'orig.ident',
        assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Vlnplot_DNA_damage.pdf', width = 10, height = 10, useDingbats = F)
VlnPlot(object = panc, features = c('H2AFX','H2AFY', 'TP53BP1')%>% tolower() %>% Caps(), pt.size = 0.1, 
        ncol = 1, assay = 'RNA')
ggsave('Vlnplot_DNA_damage_by_cluster.pdf', width = 7, height = 10, useDingbats = F)



VlnPlot(object = panc, features = c('SERPINE1','LMNB1'), pt.size = 0.9, group.by = 'Cell_type', ncol = 2,split.by = 'orig.ident')

VlnPlot(object = panc, features = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7')%>% tolower() %>% Caps(), pt.size = 0.9, group.by = 'Cell_type',
        ncol = 1,split.by = 'orig.ident',assay = 'RNA', cols =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Vlnplot_SASP.pdf', width = 13, height = 33, useDingbats = F)
VlnPlot(object = panc, features = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7') %>% tolower() %>% Caps(), pt.size = 0.1, 
        ncol = 1, assay = 'RNA')
ggsave('Vlnplot_SASP_by_cluster.pdf', width = 7, height = 25, useDingbats = F)





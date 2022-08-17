library(spatstat)
library(Signac)
library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(paletteer)
library(readxl)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/snATAC_RNA/objects/Kuppfer.rds')
panc <- NormalizeData(panc, assay = 'RNA')
cell.type <- 'kupffer'

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
# add the gene information to the object
DefaultAssay(panc) <- 'peaksinters'
Annotation(panc) <- annotations	


wd <- paste0('~/lab_Ding/work/single_cell/senescence/snATAC_RNA/markers/', cell.type)
dir.create(wd, recursive = T)
setwd(wd)

DefaultAssay(panc) <- 'RNA'
FeaturePlot(panc, features = c('Cdkn1a', 'Cdkn2a'), order = T)
ggsave('Featureplot_CDKN2A_nosplit.pdf', width =8, height = 4)

FeaturePlot(panc, features = c('Cdkn1a', 'Cdkn2a'), order = T, split.by = 'age')
ggsave('Featureplot_CDKN2A_splitted.pdf', width =8, height = 8)

FeaturePlot(panc, features = sen.core.mouse, order = T, split.by = 'sen.anno')
ggsave('Featureplot_sen.core_splitted.pdf', width =8, height = 30)

FeaturePlot(panc, features = 'Lmnb1', order = T, split.by = 'sen.anno')
VlnPlot(panc, sasp.mouse, group.by = 'sen.anno')


panc$p16.status <- ifelse(FetchData(panc, 'Cdkn2a')==0, 'p16-', 'p16+')
panc$p16.status.age <- paste(panc$p16.status, panc$age, sep = '_')
panc$p16.status.age %>% unique


FeatureScatter(panc, feature1 = 'Cdkn1a', feature2 = 'Cdkn2a')

DefaultAssay(panc) <- 'peaksinters'
panc <- RegionStats(panc, genome = BSgenome.Mmusculus.UCSC.mm10)
panc <- LinkPeaks(panc, peak.assay = 'peaksinters', expression.assay = 'RNA', genes.use = c('Cdkn2a' ,'Cdkn1a'))

gene = 'Cdkn2a'
CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'p16.status', region =gene,  extend.upstream = 2000, extend.downstream = 1000)
ggsave(paste0('CoveragePlot_kupffer_', gene, '_RNA_', cell.type, '_p16.status_linked.pdf'), useDingbats = F, width = 10, height = 5)

CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'p16.status.age', region =gene, extend.upstream = 2000, extend.downstream = 1000)
ggsave(paste0('CoveragePlot_kupffer_', gene, '_RNA_', cell.type, '_p16.status.age_linked.pdf'), useDingbats = F, width = 10, height = 8)

CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'age', region =gene, features = gene, expression.assay = "RNA", extend.upstream = 2000, extend.downstream = 1000)
ggsave(paste0('CoveragePlot_kupffer_', gene, '_RNA_', cell.type, '_age_linked.pdf'), useDingbats = F, width = 10, height = 5)

DefaultAssay(panc) <- 'peaksinters'
peak.df <- data.frame (peakAnno)
peak.df$m_coords=paste(peak.df$seqnames,"-",peak.df$start,"-",peak.df$end,sep='')


peak.df %>% subset(SYMBOL == 'Ccl2' & annotation == 'Promoter') %>% dplyr::select (m_coords) %>% unlist()

FeaturePlot(panc, features = 'chr4-89281833-89282555', order = T, split.by = 'age')
ggsave('Featureplot_p16_promoter_peak.pdf', width = 8, height = 4)

FeaturePlot(panc, features = 'chr4-89294185-89295752', order = T, split.by = 'age')
ggsave('Featureplot_p14_promoter_peak.pdf', width = 8, height = 4)

FeaturePlot(panc, features = 'chr17-29093476-29095781', order = T, split.by = 'age')
ggsave('Featureplot_p21_promoter_peak.pdf', width = 8, height = 4)

FeaturePlot(panc, features = 'chr11-82035073-82035747', order = T, split.by = 'age')
ggsave('Featureplot_ccl2_promoter_peak.pdf', width = 8, height = 4)


p16.promoter.peak <- FetchData(panc, vars = c('chr4-89281833-89282555', 'age')) 
head(p16.promoter.peak)
ggplot(data = p16.promoter.peak, aes ((`chr4-89281833-89282555`), fill = age)) +
  geom_histogram(alpha = 0.5) +
  theme_cowplot() +
  scale_fill_brewer(palette = 'Set2') +
  scale_x_continuous(breaks = seq(0,5,0.25))
ggsave('DensityPlot_chr4-89281833-89282555_p16_promoter.pdf', useDingbats = F, width = 6, height = 5)

panc$p16.promoter_open <- ifelse (p16.promoter.peak[,1] > 1, 'open', 'close')
panc$p16.promoter_open2 <- ifelse (p16.promoter.peak[,1] > 0.3, 'open', 'close')

CoveragePlot(object = panc, assay = 'peaksinters', group.by = 'p16.promoter_open', region =gene, features = gene, expression.assay = "RNA", extend.upstream = 2000, extend.downstream = 1000)
ggsave(paste0('CoveragePlot_kupffer_', gene, '_RNA_', cell.type, '_p16.promoter_open2_linked.pdf'), useDingbats = F, width = 10, height = 5)

DefaultAssay(panc) <- 'RNA'
FeaturePlot(panc, features = c('Cdkn1a', 'Cdkn2a'), order = T, split.by = c('age'))
ggsave('Featureplot_CDKN2A_by_p16_promoter_open.pdf', width =8, height = 8)

FeaturePlot(panc, features = sen.core.mouse, order = T, split.by = c('p16.promoter_open'))
ggsave('Featureplot_sen.core_by_p16_promoter.pdf', width =8, height = 30)

FeaturePlot(panc, features = sen.core.mouse, order = T, split.by = c('p16.status'))
ggsave('Featureplot_sen.core_by_p16.status.pdf', width =8, height = 30)

panc@meta.data %>% group_by(p16.status,p16.promoter_open, age) %>% tally()
DimPlot(panc, group.by = 'age', cols =c('#fc8d62', '#66c2a5'))
ggsave('Dimplot_age.pdf', useDingbats = F, width = 5, height = 4)

DefaultAssay(panc) <- 'peaksinters'
p14.promoter.peak <- FetchData(panc, vars = c('chr4-89294185-89295752', 'age'))
ggplot(data = p14.promoter.peak, aes ((`chr4-89294185-89295752`), fill = age)) +
  geom_density(alpha = 0.5) +
  theme_cowplot() +
  scale_fill_brewer(palette = 'Set2') +
  scale_x_continuous(breaks = seq(0,5,0.25))
ggsave('DensityPlot_chr4-89294185-8929575_p16_promoter.pdf', useDingbats = F, width = 6, height = 5)

panc$p14.promoter_open <- ifelse (p14.promoter.peak[,1] > 0.5, 'open', 'close')

panc@meta.data %>% group_by(age,p16.promoter_open,p16.status) %>% tally()

sen.cells.data <- FetchData(panc, vars = c('p16.promoter_open','p16.status', 'age')) 
sen.cells.old <- sen.cells.data %>% subset(p16.promoter_open == 'open' | p16.status == 'p16+') %>% subset(age == 'Old')
sen.cells.old <- rownames(sen.cells.old)
norm.cells.old <- sen.cells.data %>% subset(p16.promoter_open != 'open' & p16.status != 'p16+') %>% subset(age == 'Old')
norm.cells.old <- rownames(norm.cells.old)


panc$sen.anno.age <- case_when (rownames(panc@meta.data) %in% sen.cells.old ~ 'sen.cells.old',
                            rownames(panc@meta.data) %in% norm.cells.old ~ 'norm.cells.old',
                            TRUE ~ 'young')

sen.cells <- sen.cells.data %>% subset(p16.promoter_open == 'open' | p16.status == 'p16+')
sen.cells <- rownames(sen.cells)
norm.cells <- sen.cells.data %>% subset(p16.promoter_open != 'open' & p16.status != 'p16+')
norm.cells <- rownames(norm.cells)

panc$sen.anno <- case_when (rownames(panc@meta.data) %in% sen.cells ~ 'sen.cells',
                            rownames(panc@meta.data) %in% norm.cells ~ 'norm.cells')

DefaultAssay(panc) <- 'RNA'
Idents(panc) <- 'sen.anno.age'
lr.deg <- FindMarkers(panc, test.use = 'LR', assay = 'RNA', ident.1 = 'sen.cells.old', ident.2 = 'norm.cells.old')
wilc.deg <- FindMarkers(panc, assay = 'RNA', ident.1 = 'sen.cells.old', ident.2 = 'norm.cells.old')

Idents(panc) <- 'sen.anno'
wilc.deg <- FindMarkers(panc, assay = 'RNA', ident.1 = 'sen.cells', ident.2 = 'norm.cells')


DefaultAssay(panc) <- 'peaksinters'
Idents(panc) <- 'sen.anno.age'
lr.deg <- FindMarkers(panc, test.use = 'LR',  ident.1 = 'sen.cells.old', ident.2 = 'norm.cells.old', latent.vars = 'nCount_peaks')
head(lr.deg)
lr.deg <- merge (lr.deg, peak.df[c('annotation', 'SYMBOL', 'm_coords')], by.x = 0, by.y = 'm_coords', all.x = T)

lr.deg %>% dplyr::filter (SYMBOL %in% (SASP.human %>% tolower () %>% Caps())) %>% arrange (p_val)
df.text <- lr.deg %>% dplyr::filter ()
ggplot(data = lr.deg, aes(x = avg_log2FC, y = -log(p_val_adj), color = annotation)) +
  geom_point() +
  geom_hline(yintercept = -log(0.05))

Idents(panc) <- 'age'
lr.deg <- FindMarkers(panc, test.use = 'LR',  ident.1 = 'Old', ident.2 = 'Young', latent.vars = 'nCount_peaks')

####### plot UMAPs p16 expression and promoter opening
p16.promoter.peak <- cbind (Embeddings (panc, reduction = 'umap'), FetchData(panc, vars = c( 'p16.promoter_open', 'p16.status', 'age')) )
p16.promoter.peak$promoter_gene <- paste(p16.promoter.peak$p16.promoter_open, p16.promoter.peak$p16.status, sep = '_')
ggplot() +
  geom_point(data = subset(p16.promoter.peak, promoter_gene == 'close_p16-'), 
             aes_string(x = 'UMAP_1', y = 'UMAP_2'), 
             size = 1, color = 'lightgrey', shape = 16) +
  geom_point(data = subset(p16.promoter.peak, promoter_gene != 'close_p16-'),  
             aes_string(x = 'UMAP_1', y = 'UMAP_2', fill = 'promoter_gene' ),
             size = 3, shape = 21, stroke = 0.15, color = 'white') +
  theme_cowplot() +
  facet_wrap(~age, ncol = 2, scales = 'free')+
  scale_fill_manual(values = c('#e7298a', '#7570b3', '#66a61e'))
ggsave(paste0( 'Featureplot_Cdkn2a_promoter_expression.pdf'), useDingbats = F, width = 8.5, height = 4)


getwd()

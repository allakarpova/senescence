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
library(Signac)
library(Seurat)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

sen.core.mouse <- c('Cdkn2a', 'Bmi1', 'Trp53', 'Hmga1', 'Chek1', 'Chek2', 'Prodh', 'Tnfrsf10b', 'Cdkn1a', 'Dao')
sen.effector.mouse <- c('Ppp1ca', 'Ahcy', 'Brf1', 'Map2k3', 'Map2k6', 'Smurf2', 'Tgfb1i1', 'Srsf1', 'Angptl2')
sasp.mouse <- c('Ccl2', 'Ccl24', 'Ccl3', 'Ccl5', 'Ctnnb1', 'Cxcl1', 'Cxcl10', 'Cxcl12', 'Cxcl2', 'Cxcl16', 'Hgf', 'Hmgb1', 'Icam1', 'Igfbp2', 'Igfbp3', 'Igfbp4', 'Igfbp5', 'Igfbp6', 'Igfbp7',
                'Il15', 'Il18', 'Il1a', 'Il1b', 'Il2', 'Il6', 'Mif', 'Mmp12', 'Mmp13', 'Mmp14', 'Pgf', 'Plat', 'Timp2', 'Serpine1', 'Ccl3', 'Ccl4', 'Ang', 'Csf2', 'Kitl', 'Serpine2', 'Tnfrsf1a', 'Hgf', 'Nrg1', 'Ereg', 'Areg')
Generic.human = c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1')
SASP.human = c('IL1A', 'IL6', 'CXCL8', 'TNF', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7')


setwd('/diskmnt/Projects/SenNet_analysis/plots/mouse_liver/')
panc <- readRDS ('/diskmnt/Projects/SenNet_analysis/merged_obj/ATAC.only/mouse_liver_by_cell_type/Hepathocytes.snATAC.object.res.1.mouse_liver_old_young.rds')
DefaultAssay(panc) <- 'peaksMACS2'
cell.type <- 'Hepathocytes'
Idents(panc) <- 'Age'

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(panc) <- annotations	


gene <- 'Fgf21'
for (gene in sen.core.mouse) {
  CoveragePlot(
    object = panc,
    region = gene,
    extend.upstream = 2000,
    extend.downstream = 1000
  )
  ggsave(paste0('CoveragePlot.', cell.type, '_', gene, '.pdf'), width  = 15, height = 5, useDingbats = F)
}


for (gene in sasp.mouse) {
  CoveragePlot(
    object = panc,
    region = gene,
    extend.upstream = 2000,
    extend.downstream = 1000
  )
  ggsave(paste0('CoveragePlot.', cell.type, '_', gene, '.pdf'), width  = 15, height = 5, useDingbats = F)
}
gene = 'chr8-70831521-70831814'
CoveragePlot(
  object = panc,
  region = gene,
  extend.upstream = 2000,
  extend.downstream = 1000
)
ggsave(paste0('CoveragePlot.', cell.type, '_', gene, '.pdf'), width  = 15, height = 5, useDingbats = F)


gene <- 'Fgf21'
CoveragePlot(
  object = panc,
  region = gene,
  extend.upstream = 2000,
  extend.downstream = 1000
)
ggsave(paste0('CoveragePlot.', cell.type, '_', gene, '.pdf'), width  = 15, height = 5, useDingbats = F)


### joint RNA + ATAC object
obj <- readRDS ('/diskmnt/Projects/SenNet_analysis/analysis/dan/combo/SM001H1-Md/merge/SM001_merged_both_combo_v3.rds')
ibj.hep <- NormalizeData(ibj.hep)
DefaultAssay(ibj.hep) <- 'peaksinters'
cell.type <- 'Hepatocytes'
setwd('/diskmnt/Projects/SenNet_analysis/plots/mouse_liver_joint/')

# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Mmusculus.UCSC.mm10)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(ibj.hep) <- annotations	

gene = 'Hgf'
# link peaks to genes
ibj.hep <- LinkPeaks(
  object = ibj.hep,
  score_cutoff = 0,
  pvalue_cutoff = 0.07,
  peak.assay = "peaksinters",
  expression.assay = "RNA",
  genes.use = gene
)


#ibj.hep$anno = paste(ibj.hep$cell_type,obj$age,sep="_")

Idents(obj) = 'age'
CoveragePlot(
  object = obj,
  region = gene,
  features = gene,
  expression.assay = "SCT",
  extend.upstream = 2000,
  extend.downstream = 1000
)
ggsave(paste0('CoveragePlot.joint.', cell.type, '_', gene, '.pdf'), width  = 15, height = 5, useDingbats = F)

GetGroups(
  object = obj,
  group.by = 'age',
  idents = 'Hepathocytes'
)



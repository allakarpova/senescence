### subset and recluster cells from normal samples CCRCC
#!/usr/bin/env Rscript --vanilla
# Alla Karpova
# 102/15/2021
# load required libraries

library(optparse)
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(uwot)

option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to the integrated objected",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- '~/lab_Ding/work/single_cell/senescence/objects/Secound_mouse_integrated.rds'
out_path <- '~/lab_Ding/work/single_cell/senescence/objects/'
add_filename <- 'Mouse_BM_second'
cell.types <- c('Fibro', 'MSC', 'Osteo|OLC', 'Chondro')

#read in input
cat('opening object...\n')
panc.my <- readRDS(input.path)
cat('done\n')

for (ct in cell.types) {
  cat('Subset not tumor cells cells\n')
  panc.my.ct <- subset(x = panc.my, cells = rownames(subset(panc.my@meta.data, grepl(ct, Cell_type))))
  cat('SCTransform\n')
  panc.my.ct <- SCTransform(panc.my.ct, vars.to.regress = c("percent.mito","nCount_RNA"), return.only.var.genes = F)
  # These are now standard steps in the Seurat workflow for visualization and clustering
  panc.my.ct <- RunPCA(panc.my.ct,verbose = FALSE)
  panc.my.ct <- RunUMAP(panc.my.ct, dims = 1:30, reduction="pca")
  panc.my.ct <- FindNeighbors(panc.my.ct, dims = 1:30, reduction="pca",force.recalc=TRUE)
  # resolution 1 seems to work better than resolution 0.5
  panc.my.ct <- FindClusters(panc.my.ct, resolution = 1)
  panc.my.ct <- NormalizeData(panc.my.ct, assay = 'RNA', normalization.method = "LogNormalize", scale.factor = 10000)
  panc.my.ct <- ScaleData(panc.my.ct, assay = 'RNA')
  saveRDS(panc.my.ct, paste0(out_path, "Cell_type.object.res.1.",make.names(ct),'_', add_filename, ".RDS"))
  #  find myeloid markers
  cat('Find markers\n')
  myeloid.markers <- FindAllMarkers(panc.my.ct, only.pos = F, min.pct = 0.2, logfc.threshold = 0.2, assay = 'RNA')
  fwrite(myeloid.markers, paste0(out_path, 'Markers.normal.',make.names(ct),'_', add_filename,'.txt'), sep = '\t')
  
}






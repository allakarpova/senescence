### subset and recluster cells 
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
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
cell.types <- c('Hepathocytes', 'Kuppfer', 'Endothelial', 'stellate', 'Cholangiocytes')

dir.create (out_path, showWarnings = F, recursive = T)
setwd(out_path)
#read in input
cat('opening object...\n')
panc.my <- readRDS(input.path)
cat('done\n')

for (ct in cell.types) {
  cat(paste0('Subset ',ct,' cells\n'))
  cells.subset <- rownames(subset(panc.my@meta.data, grepl(ct, Cell_type)))
  head(cells.subset)
  panc.my.ct <- subset(x = panc.my, cells = cells.subset)
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
  saveRDS(panc.my.ct, paste0(make.names(ct), ".object.res.1.", add_filename, ".rds"))
  #  find myeloid markers
  cat('Find markers\n')
  myeloid.markers <- FindAllMarkers(panc.my.ct, only.pos = F, min.pct = 0.2, logfc.threshold = 0.2, assay = 'RNA')
  fwrite(myeloid.markers, paste0('Markers.normal.',make.names(ct),'_', add_filename,'.txt'), sep = '\t')
  
}






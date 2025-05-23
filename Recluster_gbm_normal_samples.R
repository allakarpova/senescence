### subset and recluster not tumor cells for snRNA-seq object GBM
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
  # make_option(c("-m", "--metadata"),
  #             type="character",
  #             default=NULL, 
  #             help="output folder path",
  #             metavar="character"),
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

#read in input
cat('opening object...\n')
panc <- readRDS(input.path)
cat('done\n')

my.metadata <- fread('/diskmnt/Projects/cptac_scratch_4/CPTAC3_GBM/single_cell_data_freeze/v2/snRNA/snRNA_metadata.v2.1.tsv.gz', data.table = F) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)
panc <- AddMetaData(panc, metadata = my.metadata)

cat('Subset normal cells\n')
panc.my <- subset(x = panc, cells = rownames(subset(panc@meta.data, grepl('PT-', case_id))))

cat('SCTransform\n')
panc.my <- SCTransform(panc.my, vars.to.regress = c("percent.mito","nCount_RNA"), return.only.var.genes = F)
# These are now standard steps in the Seurat workflow for visualization and clustering
panc.my <- RunPCA(panc.my,verbose = FALSE)
panc.my <- RunUMAP(panc.my, dims = 1:30, reduction="pca")
panc.my <- FindNeighbors(panc.my, dims = 1:30, reduction="pca",force.recalc=TRUE)
# resolution 1 seems to work better than resolution 0.5
panc.my <- FindClusters(panc.my, resolution = 1)
panc.my <- NormalizeData(panc.my, assay = 'RNA', normalization.method = "LogNormalize", scale.factor = 10000)
panc.my <- ScaleData(panc.my, assay = 'RNA')

saveRDS(panc.my, paste0(out_path, "Normal.object.res1.", add_filename, ".RDS"))

#  find myeloid markers
cat('Find markers\n')
myeloid.markers <- FindAllMarkers(panc.my, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25, assay = 'RNA')
fwrite(myeloid.markers, paste0(out_path, 'Markers.normal.', add_filename,'.txt'), sep = '\t')


for (ct in c('OPC', 'Neuron', 'Oligodendrocyte', "Pericyte/vSMC")) {
  cat(paste('Subset' ,ct, 'cells\n'))
  panc.my.ct <- subset(x = panc.my, cells = rownames(subset(panc.my@meta.data, cell_type_v20210128 == ct)))
  cat('SCTransform\n')
  panc.my.ct <- SCTransform(panc.my.ct, vars.to.regress = c("percent.mito","nCount_RNA"), return.only.var.genes = F)
  # These are now standard steps in the Seurat workflow for visualization and clustering
  panc.my.ct <- RunPCA(panc.my.ct,verbose = FALSE)
  panc.my.ct <- RunUMAP(panc.my.ct, dims = 1:30, reduction="pca")
  panc.my.ct <- FindNeighbors(panc.my.ct, dims = 1:30, reduction="pca",force.recalc=TRUE)
  # resolution 1 seems to work better than resolution 0.5
  panc.my.ct <- FindClusters(panc.my.ct, resolution = 0.7)
  panc.my.ct <- NormalizeData(panc.my.ct, assay = 'RNA', normalization.method = "LogNormalize", scale.factor = 10000)
  panc.my.ct <- ScaleData(panc.my.ct, assay = 'RNA')
  saveRDS(panc.my.ct, paste0(out_path, "Normal.object.res.7.",ct,'_', add_filename, ".RDS"))
  #  find myeloid markers
  cat('Find markers\n')
  myeloid.markers <- FindAllMarkers(panc.my.ct, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25, assay = 'RNA')
  fwrite(myeloid.markers, paste0(out_path, 'Markers.normal.',ct,'_', add_filename,'.txt'), sep = '\t')
  
}






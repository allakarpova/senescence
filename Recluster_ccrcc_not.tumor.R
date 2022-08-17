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

my.metadata <- fread('/diskmnt/Projects/Users/allakarpova/Projects/ccRCC/ligand_receptor/ccRCC_metadata_AK_v20210125_scrublet_fixed_integr_barcocde.txt', data.table = F) %>% 
  data.frame(row.names = 2, check.rows = F, check.names = F)
panc <- AddMetaData(panc, metadata = my.metadata)

cat('Subset non tumor cells\n')
panc.my <- subset(x = panc, cells = rownames(subset(panc@meta.data, !grepl('umor|oublet', Cell_type.shorter))))

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

saveRDS(panc.my, paste0(out_path, "Non.tumor.object.res1.", add_filename, ".RDS"))

#  find myeloid markers
cat('Find markers\n')
myeloid.markers <- FindAllMarkers(panc.my, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25, assay = 'RNA')
fwrite(myeloid.markers, paste0(out_path, 'Markers.not.tumor.', add_filename,'.txt'), sep = '\t')



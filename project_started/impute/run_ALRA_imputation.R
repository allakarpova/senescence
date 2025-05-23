#run ALRA imputation
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readxl))
suppressMessages(library(optparse))
library(SeuratWrappers)

option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object to score",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output

dir.create(out_path, showWarnings = F, recursive = TRUE)
setwd(out_path)

panc <- readRDS(input.path)

DefaultAssay(panc) <- 'SCT'
panc <- RunALRA(panc)

panc <- panc %>% 
  FindVariableFeatures(assay = 'alra') %>%
  ScaleData() %>%
  RunPCA(assay = 'alra', reduction.name = 'pca.alra') %>%
  RunUMAP( dims = 1:30, assay = 'alra', reduction = 'pca.alra', 
           reduction.name = "umap.alra", 
           reduction.key = "alraUMAP_", verbose = F) 

base.name <- tools::file_path_sans_ext(basename(input.path))
saveRDS(panc, paste0(base.name, '.imputed.rds'))



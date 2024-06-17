#run ALRA imputation
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readxl))
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
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="foo",
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-m","--metadata"),
              type="character",
              default='/diskmnt/Projects/HTAN_analysis_2/BRCA/Analyses/Alla/DCIS_project/cell_typing/st_obj/',
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-r","--reduction"),
              type="character",
              default='wnn.umap',
              help = "wnn.umap or umap",
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

base.name <- tools::file_path_sans_ext(basename(input.path))
saveRDS(panc, paste0(base.name, '.imputed.rds'))



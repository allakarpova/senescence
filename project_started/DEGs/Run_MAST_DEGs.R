# Fubd DEGs by cell types senescent phenotypes
###libraries
##################

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
library(MAST)

################################

###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object",
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
              metavar="character"),
  make_option(c("--metadata"),
              type="character",
              default=NULL, 
              help="path to metadata"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
panc <- readRDS(input.path)

library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024 ^ 3)

if(!is.null(meta.path)) {
  meta <- fread(meta.path, header = TRUE) %>% 
    column_to_rownames('V1') %>%
    dplyr::select(all_of(cell_column))
  panc <- AddMetaData(panc, meta)
  
}

print(head(panc@meta.data))
Idents(panc) <- cell_column
DefaultAssay(panc) <- 'RNA'
unique(Idents(panc))

deg <- FindAllMarkers(panc, assay = "RNA", 
                      logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                      test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F)

fwrite(deg, paste0('DEG_findAllMarkers_by_',cell_column,'_RNA_MAST_',add_filename,'.txt'), sep = '\t', row.names = T)











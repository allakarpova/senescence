# Fubd daps by cell types senescent phenotypes
###libraries
##################
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(EnsDb.Hsapiens.v100))

################################

###options###
######################
option_list = list(
  
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
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
input.path <- '/diskmnt/Projects/SenNet_analysis/Main.analysis/data_freeze/v1.1/combo_RNA_ATAC/mCRC_SenNet/each_cell_type/34_liver_14_mCRC_combo_fibroblasts.stellate.chromvar.rds'
out_path <- opt$output
add_filename <- 'Stellate.fibroblasts'
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
DefaultAssay(panc) <- 'ATAC_merged'
unique(Idents(panc))


dap1 <- FindMarkers(panc, ident.1 = 'PS fibroblasts CDKN1A CCL2',
                    ident.2 = c('Portal fibroblasts'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR',only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS fibroblasts CDKN1A CCL2',
         ident.2 = 'Portal fibroblasts')

dap2 <- FindMarkers(panc,  ident.1 = 'PS fibroblasts CDKN1A CCL2',
                    ident.2 = c('Activated HSCs', "Activating HSCs", 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR',  only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS fibroblasts CDKN1A CCL2',
         ident.2 = 'HSCs')


fwrite(rbind(dap1, dap2), paste0('DAP_findMarkers_PS_fibro_ATAC_merged_LR_no_latent.txt'), sep = '\t', row.names = F)


dap1 <- FindMarkers(panc, ident.1 = c('Portal fibroblasts'),
                    ident.2 = c('Activated HSCs', "Activating HSCs", 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR',  only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Portal fibroblasts',
         ident.2 =  "HSCs")

dap2 <- FindMarkers(panc,  ident.1 = c('Activated HSCs'),
                    ident.2 = c("Activating HSCs", 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Activated HSCs',
         ident.2 =  "Other HSCs")

dap3 <- FindMarkers(panc,ident.1 = c('Activating HSCs'),
                    ident.2 = c( 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Activating HSCs',
         ident.2 =  "Quiescent HSCs")

dap4 <- FindMarkers(panc,  ident.1 = c( 'Quiescent HSCs'),
                    ident.2 = c( 'Activated HSCs', "Activating HSCs"),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR',  only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Quiescent HSCs',
         ident.2 =  "Other HSCs")


fwrite(rbind(dap1, dap2,dap3, dap4), paste0('DAP_findMarkers_normal_fibro_HSC_ATAC_merged_LR_no_latent.txt'), sep = '\t', row.names = F)



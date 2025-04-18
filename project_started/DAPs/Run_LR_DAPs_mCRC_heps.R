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
input.path <- '/diskmnt/Projects/SenNet_analysis/Main.analysis/data_freeze/v1.1/combo_RNA_ATAC/mCRC_SenNet/each_cell_type/34_liver_14_mCRC_combo_Hepatocytes.chromvar.rds'
out_path <- opt$output
add_filename <- 'Hepatocytes'
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


dap1 <- FindMarkers(panc, ident.1 = 'PS Hepatocytes CDKN1A SERPINE1',
                    ident.2 = c('Hepatocytes Zone 1', "Hepatocytes Zone 2/3"),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS Hepatocytes CDKN1A SERPINE1',
         ident.2 = 'Hepatocytes')

dap2 <- FindMarkers(panc, ident.1 = 'PS Hepatocytes CDKN1A',
                    ident.2 = c('Hepatocytes Zone 1', "Hepatocytes Zone 2/3"),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS Hepatocytes CDKN1A',
         ident.2 = 'Hepatocytes')

dap3 <- FindMarkers(panc,  ident.1 = 'Hepatocytes CRP',
                    ident.2 = c('Hepatocytes Zone 1', "Hepatocytes Zone 2/3"),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Hepatocytes CRP',
         ident.2 = 'Hepatocytes')

fwrite(rbind(dap1, dap2, dap3), paste0('DAP_findMarkers_PS_Heps_ATAC_merged_LR.txt'), sep = '\t', row.names = F)


dap1 <- FindMarkers(panc,  ident.1 = 'Hepatocytes Zone 1',
                    ident.2 = "Hepatocytes Zone 2/3",
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'LR', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Hepatocytes Zone 1',
         ident.2 =  "Hepatocytes Zone 2/3")


fwrite(dap1, paste0('DAP_findMarkers_Heps_zone1_vs_3_ATAC_merged_LR.txt'), sep = '\t', row.names = F)



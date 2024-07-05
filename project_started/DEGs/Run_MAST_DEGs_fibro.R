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
input.path <- '/diskmnt/Projects/SenNet_analysis/Main.analysis/data_freeze/v1.1/RNA_only/each_cell_type_harmony/37_Merged_normalized_liver_snRNA_no_doublets_PF_HSC_NoVSMC.rds'
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
DefaultAssay(panc) <- 'RNA'
unique(Idents(panc))


deg1 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS fibroblasts CDKN1A CCL2',
                    ident.2 = c('Portal fibroblasts'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS fibroblasts CDKN1A CCL2',
         ident.2 = 'Portal fibroblasts')

deg2 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS fibroblasts CDKN1A CCL2',
                    ident.2 = c('Activated HSCs', "Activating HSCs", 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS fibroblasts CDKN1A CCL2',
         ident.2 = 'HSCs')


fwrite(rbind(deg1, deg2), paste0('DEG_findMarkers_PS_fibro_RNA_MAST.txt'), sep = '\t', row.names = F)


deg1 <- FindMarkers(panc, assay = "RNA", ident.1 = c('Portal fibroblasts'),
                    ident.2 = c('Activated HSCs', "Activating HSCs", 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Portal fibroblasts',
         ident.2 =  "HSCs")

deg2 <- FindMarkers(panc, assay = "RNA", ident.1 = c('Activated HSCs'),
                    ident.2 = c("Activating HSCs", 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Activated HSCs',
         ident.2 =  "Other HSCs")

deg3 <- FindMarkers(panc, assay = "RNA", ident.1 = c('Activating HSCs'),
                    ident.2 = c( 'Quiescent HSCs'),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Activating HSCs',
         ident.2 =  "Quiescent HSCs")

deg4 <- FindMarkers(panc, assay = "RNA", ident.1 = c( 'Quiescent HSCs'),
                    ident.2 = c( 'Activated HSCs', "Activating HSCs"),
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Quiescent HSCs',
         ident.2 =  "Other HSCs")


fwrite(rbind(deg1, deg2,deg3, deg4), paste0('DEG_findMarkers_normal_fibro_HSC_RNA_MAST.txt'), sep = '\t', row.names = F)



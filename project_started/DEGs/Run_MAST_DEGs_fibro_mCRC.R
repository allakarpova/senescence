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
input.path <- '/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/merge_SenNet_mCRC_snRNA_2/37_liver_18_mCRC_snRNA_Fibroblasts_harmony.rds'
out_path <- opt$output
add_filename <- 'Fibroblasts'
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
                    ident.2 = 'Portal fibroblasts',
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

deg3 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS fibroblasts CDKN1A CCL2',
                    ident.2 = 'PS tumor fibroblasts CDKN1A SERPINE1',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = "PS fibroblasts CDKN1A CCL2",
         ident.2 = "PS tumor fibroblasts CDKN1A SERPINE1")

deg7 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS fibroblasts CDKN1A CCL2',
                    ident.2 = 'PS tumor fibroblasts CDKN2A/2B',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = "PS fibroblasts CDKN1A CCL2",
         ident.2 = "PS tumor fibroblasts CDKN2A/2B")


deg4 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS tumor fibroblasts CDKN2A/2B',
                    ident.2 = 'PS tumor fibroblasts CDKN1A SERPINE1',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = "PS tumor fibroblasts CDKN2A/2B",
         ident.2 = "PS tumor fibroblasts CDKN1A SERPINE1")

deg5 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS tumor fibroblasts CDKN2A/2B',
                    ident.2 = 'Tumor fibroblasts',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = "PS tumor fibroblasts CDKN2A/2B",
         ident.2 = "Tumor fibroblasts")

deg6 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS tumor fibroblasts CDKN1A SERPINE1',
                    ident.2 = 'Tumor fibroblasts',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = "PS tumor fibroblasts CDKN1A SERPINE1",
         ident.2 = "Tumor fibroblasts")



fwrite(rbind(deg1, deg2, deg3,deg7, deg4, deg5,deg6), paste0('DEG_findMarkers_PS_fibros_RNA_MAST.txt'), sep = '\t', row.names = F)


deg1 <- FindMarkers(panc, assay = "RNA", subset.ident = 'PS fibroblasts CDKN1A CCL2',
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = 'PS fibroblasts CDKN1A CCL2',
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')

deg2 <- FindMarkers(panc, assay = "RNA", subset.ident = 'Activated HSCs',
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = 'Activated HSCs',
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')

deg3 <- FindMarkers(panc, assay = "RNA", subset.ident = 'Activating HSCs',
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = 'Activating HSCs',
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')

deg4 <- FindMarkers(panc, assay = "RNA", subset.ident = "Quiescent HSCs",
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = "Quiescent HSCs",
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')
fwrite(rbind(deg1, deg2, deg3, deg4), paste0('DEG_findMarkers_mCRC_vs_normal_liver_all_liver_fibro_RNA_MAST.txt'), sep = '\t', row.names = F)

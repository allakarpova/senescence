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
input.path <- '/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/merge_SenNet_mCRC_snRNA_2/37_liver_18_mCRC_snRNA_Hepatocytes_harmony.rds'
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


deg1 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS Hepatocytes CDKN1A SERPINE1',
                    ident.2 = 'Hepatocytes',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS Hepatocytes CDKN1A SERPINE1',
         ident.2 = 'Hepatocytes')

deg2 <- FindMarkers(panc, assay = "RNA", ident.1 = 'PS Hepatocytes CDKN1A',
                    ident.2 = 'Hepatocytes',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'PS Hepatocytes CDKN1A',
         ident.2 = 'Hepatocytes')

deg3 <- FindMarkers(panc, assay = "RNA", ident.1 = 'Hepatocytes CRP',
                    ident.2 = 'Hepatocytes',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = 'Hepatocytes CRP',
         ident.2 = 'Hepatocytes')


deg4 <- FindMarkers(panc, assay = "RNA", ident.1 = "PS Hepatocytes CDKN2A/2B",
                    ident.2 = 'Hepatocytes',
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.02,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         ident.1 = "PS Hepatocytes CDKN2A/2B",
         ident.2 = 'Hepatocytes')

fwrite(rbind(deg1, deg2, deg3, deg4), paste0('DEG_findMarkers_PS_Heps_RNA_MAST.txt'), sep = '\t', row.names = F)


deg1 <- FindMarkers(panc, assay = "RNA", subset.ident = 'PS Hepatocytes CDKN1A SERPINE1',
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.0,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
        celltype = 'PS Hepatocytes CDKN1A SERPINE1',
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')

deg2 <- FindMarkers(panc, assay = "RNA", subset.ident = 'PS Hepatocytes CDKN1A',
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.0,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = 'PS Hepatocytes CDKN1A',
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')

deg3 <- FindMarkers(panc, assay = "RNA", subset.ident = 'Hepatocytes CRP',
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.0,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = 'Hepatocytes CRP',
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')

deg4 <- FindMarkers(panc, assay = "RNA", subset.ident = "PS Hepatocytes CDKN2A/2B",
                    group.by = 'Cohort',  ident.1 = 'mCRC liver', ident.2 = 'Normal liver', 
                    logfc.threshold = 0.1, min.pct = 0, min.diff.pct = 0.0,
                    test.use = 'MAST', latent.vars = 'Patient_ID', only.pos = F) %>%
  mutate(gene=rownames(.),
         celltype = "PS Hepatocytes CDKN2A/2B",
         ident.1 = 'mCRC liver',
         ident.2 = 'Normal liver')
fwrite(rbind(deg1, deg2, deg3, deg4), paste0('DEG_findMarkers_mCRC_vs_normal_liver_PS_Heps_RNA_MAST.txt'), sep = '\t', row.names = F)

# subset immune cells from the ATAC object and recall peaks with MACS2
###libraries
##################

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(org.Hs.eg.db))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

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
  make_option(c("-m","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file

dir.create(out_path, showWarnings = F)
setwd(out_path)

cat('opening object \n')
panc <- readRDS(input.path)

my.metadata <- fread(meta.path, data.table = F) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F) %>%
  dplyr::select(-seurat_clusters)

panc <- AddMetaData(panc, my.metadata)

filter <- dplyr::filter
select <- dplyr::select

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

Idents(panc) <- 'cell_type'
DefaultAssay(panc) <- 'SCT'
print(head(panc@meta.data))

print(unique(Idents(panc)))

degs <- unique(as.character(Idents(panc))) %>% map (function(ct) {
  
  print(ct)
  old.cells <- rownames(panc@meta.data %>% filter(cell_type==ct) %>% filter(Age.group == 'Old'))
  young.cells <- rownames(panc@meta.data %>% filter(cell_type==ct) %>% filter(Age.group == 'Young'))
  print(length(old.cells))
  print(length(young.cells))
  degs <- NA
  if ((ct != 'Doublet') & (length(old.cells)>50) & (length(young.cells)>50)) {
    degs <- FindMarkers(
      object = panc,
      subset.ident = ct,group.by = 'Age.group', ident.1 = 'Old',
      #cells.1 = old.cells, cells.2 = young.cells,
      only.pos = F,
      test.use = 'LR',
      latent.vars = 'nCount_SCT', 
      pseudocount.use = 0.1
    )
    degs$cell_type <- ct
    degs$ident.1 <- 'Old'
    degs$ident.2 <- 'Young'
    degs$gene <- rownames(degs)
    
  } 
  return(degs)
  
})
print(str(degs))
degs <- degs[!is.na(degs)]
degs <- degs  %>% rbindlist()
fwrite(degs, paste0('DEG_findAllMarkers_old_vs_young_cell_types_',add_filename,'.txt'), sep = '\t', row.names = F)


if ("ATAC_merged" %in% Assays(panc)) {
  
  DefaultAssay(panc) <- 'ATAC_merged'
  unique(Idents(panc))
  daps <- unique(as.character(Idents(panc))) %>% map (function(ct) {
    
    print(ct)
    old.cells <- rownames(panc@meta.data %>% filter(cell_type==ct) %>% filter(Age.group == 'Old'))
    young.cells <- rownames(panc@meta.data %>% filter(cell_type==ct) %>% filter(Age.group == 'Young'))
    
    if ((ct != 'Doublet') & (length(old.cells)>50) & (length(young.cells)>50)) {
      degs <- FindMarkers(
        object = panc,
        subset.ident = ct,group.by = 'Age.group', ident.1 = 'Old',
        #cells.1 = old.cells, cells.2 = young.cells,
        only.pos = F,
        test.use = 'LR',
        latent.vars = 'nCount_ATAC_merged', 
        pseudocount.use = 0.1
      )
      degs$cell_type <- ct
      degs$ident.1 <- 'Old'
      degs$ident.2 <- 'Young'
      degs$peak <- rownames(degs)
    }
    
    return(degs)
    
  })
  
  daps <- daps[!is.na(daps)]
  daps <- daps  %>% rbindlist()
  fwrite(daps, paste0('DAP_findAllMarkers_old_vs_young_cell_types_',add_filename,'.txt'), sep = '\t', row.names = F)

}









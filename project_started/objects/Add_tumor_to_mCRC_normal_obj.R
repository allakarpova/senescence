# add some tumor cells to the mCRC object for Xenium label transfer
###libraries
##################
library(future)

#plan("multicore", workers = 20)
#options(future.globals.maxSize = 100 * 1024 ^ 3)

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
#suppressMessages(library(doParallel))



################################

#####################################
####### FUNCTIONS ##################
####################################



############################################

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
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

NormalizeRNA <- function(obj){
  ######## Normalize RNA
  DefaultAssay(obj) <- 'RNA'
  cat('normalizing RNA\n')
  obj <- obj %>%
    NormalizeData(assay = 'RNA') %>%
    CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F) %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", 
                          "percent.mt", 
                          "S.Score", 
                          "G2M.Score"
      ),
      return.only.var.genes = TRUE, verbose = T) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
    RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = T) %>%
    FindNeighbors( dims = 1:50) %>%
    FindClusters(resolution = 2, verbose = FALSE)
  
  return(obj)
}

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
  column_to_rownames('V1') 

panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)

# Recluster cell types 
###libraries
##################
library(future)

#plan("multicore", workers = 20)
#options(future.globals.maxSize = 100 * 1024 ^ 3)

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
#suppressMessages(library(doParallel))



################################

#####################################
####### FUNCTIONS ##################
####################################



############################################

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
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

NormalizeRNA <- function(obj){
  ######## Normalize RNA
  DefaultAssay(obj) <- 'RNA'
  cat('normalizing RNA\n')
  obj <- obj %>%
    NormalizeData(assay = 'RNA') %>%
    CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F) %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", 
                          "percent.mt", 
                          "S.Score", 
                          "G2M.Score"
      ),
      return.only.var.genes = TRUE, verbose = T) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
    RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = T) %>%
    FindNeighbors( dims = 1:50) %>%
    FindClusters(resolution = 2, verbose = FALSE)
  
  return(obj)
}

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

if(!file.exists(paste0(add_filename,".rds"))) {
  my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
    column_to_rownames('V1') 
  
  panc.my <- readRDS(input.path)
  panc.my <- AddMetaData(panc.my, my.metadata)
  
  relevant.patients <- panc.my@meta.data %>% pull(Patient_ID) %>% unique()
  
  tumor.panc <- readRDS('/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/merge_mCRC_snRNA_2/30_Merged_normalized_mCRC_snRNA.rds')
  
  set.seed(666)
  
  cells.toadd <- tumor.panc@meta.data %>% 
    filter(Patient_ID %in% relevant.patients) %>%
    filter(merged_cell_type == 'Tumor') %>%
    dplyr::sample_n(size=5000) %>%
    rownames(.)
  
  tumor.panc.5k <- subset(tumor.panc, cells = cells.toadd)
  tumor.panc.5k@meta.data[[cell_column]] <- 'Tumor'
  
  common.columns <- intersect(colnames(tumor.panc.5k@meta.data), colnames(panc.my@meta.data))
  tumor.panc.5k@meta.data <- tumor.panc.5k@meta.data[,common.columns]
  
  
  int.sub <- merge(panc.my, tumor.panc.5k)
  
  print(dim(int.sub))

  int.sub <- NormalizeRNA(int.sub)
  
  saveRDS(int.sub,  paste0(add_filename,".rds"))
  
  
} else {
  int.sub <- readRDS(paste0(add_filename,".rds"))
  
}
  

int.sub@meta.data <- int.sub@meta.datamutate(cell_type_broad = case_when(grepl('Hepatocytes', cell_type_sen_mCRC) ~ 'Hepatocytes',
                                   grepl('fibroblas', cell_type_sen_mCRC) ~ 'Portal fibroblasts',
                                   grepl('HSC', cell_type_sen_mCRC) ~ 'Hepatic stellate cells',
                                   grepl('PS LSEC', cell_type_sen_mCRC) ~ 'Mid lobular LSECs',
                                   grepl('Cholangiocytes', cell_type_sen_mCRC) ~ 'Cholangiocytes',
                                   TRUE ~ cell_type_sen_mCRC))

fwrite(int.sub@meta.data,  paste0(add_filename,".metadata.tsv"), sep='\t', row.names = TRUE)
DimPlot(int.sub, group.by = cell_column, label = TRUE)
ggsave(paste0('Dimplot_', cell_column,'.pdf'), width = 15, height = 5, useDingbats = F)












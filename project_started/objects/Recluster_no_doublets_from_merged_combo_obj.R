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
suppressMessages(library(doParallel))



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

NormalizeAll <- function(obj){
  ######## Normalize RNA
  DefaultAssay(obj) <- 'RNA'
  cat('normalizing RNA\n')
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-") 
  obj <- obj %>%
    NormalizeData(assay = 'RNA') %>%
    CellCycleScoring(s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
  
  obj <- obj %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      return.only.var.genes = FALSE, verbose = F) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = F) %>%
    RunUMAP(dims = 1:30,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = F) 
  
  ######## Normalize ATAC
  DefaultAssay(obj) <- "ATAC_merged"
  
  Annotation(obj) <- annotations
  
  cat('normalizing ATAC\n')
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 'q10') %>%
    RunSVD(verbose = F) %>%
    RunUMAP(reduction = 'lsi', 
            dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = F)
  
  # do wnn analysis
  cat('doing WNN\n')
  obj <- FindMultiModalNeighbors(obj, 
                                     reduction.list = list("pca", "lsi"), 
                                     dims.list = list(1:30, 2:30)) %>%
    RunUMAP( nn.name = "weighted.nn", 
             reduction.name = "wnn.umap", 
             reduction.key = "wnnUMAP_")  %>%
    FindClusters(graph.name = "wsnn", algorithm = 3,resolution = 2, verbose = T)
  
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

my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)

panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)

panc.my <- subset(x = panc.my, cells = rownames(dplyr::filter(panc.my@meta.data, !grepl('Doubl', .data[[cell_column]]))))

panc.my <- NormalizeAll(panc.my)

saveRDS(panc.my,  paste0(add_filename,".rds"))









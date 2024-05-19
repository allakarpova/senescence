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

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

panc.my <- readRDS(input.path)

panc@meta.data$Sample[is.na(panc@meta.data$Sample)] <- 'NA'
panc <- subset(panc, Sample=='NA' , invert = TRUE)

samples.id <- panc@meta.data %>% pull(Sample) %>% unique()

panc.my <- NormalizeRNA(panc.my)

saveRDS(panc.my,  paste0(length(samples.id),"_Merged_normalized_", add_filename,".rds"))









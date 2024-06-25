# Merge same cell types from mCRC and liver datasets
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
suppressMessages(library(harmony))
suppressMessages(library(googlesheets4))
gs4_deauth()

################################

#####################################
####### FUNCTIONS ##################
####################################



############################################

###options###
######################
option_list = list(
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

runHarmonyNormalization <- function(obj, dims=30, column = 'Patient_ID') {
  
  obj <- obj %>%
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    FindNeighbors(reduction = "harmony", dims = 1:dims) %>%
    FindClusters(verbose = FALSE, resolution = 1, 
                 algorithm = 1) %>%
    RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:dims)
  
  return(obj)
  
}

# read in initial arguments

out_path <- opt$output
add_filename <- opt$extra

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

obj.paths <- read_sheet("https://docs.google.com/spreadsheets/d/1VeWme__vvVHAhHaQB3wCvAGuq-w3WrhZ5cT-7Mh7Sr0/edit#gid=0", 
                         sheet = "liver_mCRC_cell_type_obj")

pwalk(list(obj.paths$Cell_type, obj.paths$liver, obj.paths$mCRC), function(ct, l.path, m.path) {
  obj1 <- readRDS(l.path)
  obj2 <- readRDS(m.path)
  
  obj1@meta.data$Cohort <- 'Normal liver'
  obj2@meta.data$Cohort <- 'mCRC liver'
  
  int.sub <- merge(obj1, obj2)

  print(dim(int.sub))
  int.sub <- NormalizeRNA(int.sub)
  int.sub <- runHarmonyNormalization(int.sub)
  
  saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), "_harmony.rds"))
  
  ct <- make.names(ct)
  
  DimPlot(int.sub, reduction='umap.harmony', group.by = 'Patient_ID')
  ggsave(glue::glue("Dimplot_{ct}_Patient.id.pdf"), width = 8, height = 5)
  
  DimPlot(int.sub, reduction='umap.harmony', group.by = 'Cohort')
  ggsave(glue::glue("Dimplot_{ct}_Cohort.pdf"), width = 8, height = 5)
  
  DimPlot(int.sub, reduction='umap.harmony', group.by = 'seurat_clusters')
  ggsave(glue::glue("Dimplot_{ct}_seurat_clusters.pdf"), width = 6.5, height = 5)
    

  
})









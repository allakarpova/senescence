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

runAllNormalization <- function(int.sub, assay = "ATAC_merged") {
  int.sub <- int.sub %>%
    NormalizeData(assay = 'RNA') %>%
    CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  
  int.sub <- int.sub %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      return.only.var.genes = TRUE, verbose = F) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = F) %>%
    RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = F) 
  
  ######## Normalize ATAC
  DefaultAssay(int.sub) <- assay
  
  Annotation(int.sub) <- annotations
  
  cat('normalizing ATAC\n')
  int.sub <- int.sub %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 'q10') %>%
    RunSVD(verbose = F) %>%
    RunUMAP(reduction = 'lsi', 
            dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = F)
  
  # do wnn analysis
  cat('doing WNN\n')
  int.sub <- FindMultiModalNeighbors(int.sub, 
                                     reduction.list = list("pca", "lsi"), 
                                     dims.list = list(1:50, 2:50)) %>%
    RunUMAP( nn.name = "weighted.nn", 
             reduction.name = "wnn.umap", 
             reduction.key = "wnnUMAP_")  %>%
    FindClusters(graph.name = "wsnn", algorithm = 1, resolution = 1, verbose = T)
  return(int.sub)
}

runAllChromvar <- function(obj, assay = 'ATAC_merged') {
  DefaultAssay(obj) <- assay
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
  )
  
  # Scan the DNA sequence of each peak for the presence of each motif
  motif.matrix <- CreateMotifMatrix(
    features = granges(obj),
    pwm = pfm,
    genome = 'BSgenome.Hsapiens.UCSC.hg38',
    use.counts = FALSE
  )
  
  # Create a new Motif object to store the results
  motif <- CreateMotifObject(
    data = motif.matrix,
    pwm = pfm
  )
  
  # Add the Motif object to the assay
  obj <- SetAssayData(
    object = obj,
    assay = assay,
    slot = 'motifs',
    new.data = motif
  )
  
  cat('doing chromvar\n')
  obj <- RegionStats(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  obj <- RunChromVAR(
    object = obj,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  DefaultAssay(obj) <- 'chromvar'
  obj@assays$chromvar@scale.data <- obj@assays$chromvar@data
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
                        sheet = "liver_mCRC_cell_type_obj") %>%
    filter(datatype=='combo')

pwalk(list(obj.paths$Cell_type, obj.paths$liver, obj.paths$mCRC), function(ct, l.path, m.path) {
  obj1 <- readRDS(l.path)
  obj2 <- readRDS(m.path)
  
  obj1@meta.data$Cohort <- 'Normal liver'
  obj2@meta.data$Cohort <- 'mCRC liver'
  
  DefaultAssay(obj1) <- 'RNA'
  obj1 <- DietSeurat(obj1, assays = c('RNA', 'ATAC_merged'))
  DefaultAssay(obj2) <- 'RNA'
  obj2 <- DietSeurat(obj2, assays = c('RNA', 'ATAC_merged'))
  int.sub <- merge(obj1, obj2)
  
  print(dim(int.sub))
  int.sub <- runAllNormalization(int.sub, assay = assay.towork)
  
  
  saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".rds"))
  
  ct <- make.names(ct)
  
  DimPlot(int.sub, reduction='wnn.umap', group.by = 'Patient_ID')
  ggsave(glue::glue("Dimplot_{ct}_Patient.id.pdf"), width = 8, height = 5)
  
  DimPlot(int.sub, reduction='wnn.umap', group.by = 'Cohort')
  ggsave(glue::glue("Dimplot_{ct}_Cohort.pdf"), width = 8, height = 5)
  
  DimPlot(int.sub, reduction='wnn.umap', group.by = 'seurat_clusters')
  ggsave(glue::glue("Dimplot_{ct}_seurat_clusters.pdf"), width = 6.5, height = 5)
  
  int.sub <- runAllChromvar(int.sub)
  saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".chromvar.rds"))
  
  int.sub@meta.data %>% fwrite(paste0(add_filename,"_",make.names(ct), ".metadata.tsv"), sep='\t', row.names = TRUE)
})









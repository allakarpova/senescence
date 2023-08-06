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
DefaultAssay(panc.my) <- 'RNA'
panc.my <- DietSeurat(panc.my, assays = 'RNA')

int.sub <- subset(x = panc.my, cells = rownames(dplyr::filter(panc.my@meta.data, !grepl('Doubl', .data[[cell_column]]))))

ct<- unique(int.sub@meta.data[[cell_column]])
print(ct)
print(dim(int.sub))


all.rna.list <- SplitObject(int.sub, split.by = 'Patient_ID')

all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  
  x <- x %>% SCTransform(
    assay = 'RNA',
    vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
    return.only.var.genes = FALSE, verbose = F)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = all.rna.list, nfeatures = 3000)
all.rna.list <- PrepSCTIntegration(object.list = all.rna.list, anchor.features = features)
all.rna.list <- lapply(X = all.rna.list, FUN = RunPCA, features = features)

rna.anchors <- FindIntegrationAnchors(object.list = all.rna.list, normalization.method = "SCT",
                                      anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
int <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:30)
int <- RunPCA(int, verbose = FALSE)
int <- RunUMAP(int, reduction = "pca", dims = 1:30)
int <- FindNeighbors(int, reduction = "pca", dims = 1:30)
int <- FindClusters(int, resolution = 1, algorithm = 4, method = "igraph")
int <- PrepSCTFindMarkers(int)

cat('saving the object...\n')
saveRDS(int,  saveRDS(int,  paste0("Integrated_RNA_", add_filename,"_",ct, ".rds")))





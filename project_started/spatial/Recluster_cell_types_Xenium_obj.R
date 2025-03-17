# Remove doublets in mouse samples and run harmony on c('Mouse_ID', 'Cohort')
###libraries
##################
#plan("multicore", workers = 20)

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(future))
options(future.globals.maxSize = 100 * 1024 ^ 3)

suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(harmony))
options(Seurat.object.assay.version = "v3")


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

normalizeXenium <- function(xenium.sub, num_pcs = 30) {
  xenium.sub <- SCTransform(xenium.sub, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2")
  xenium.sub <- RunPCA(xenium.sub, npcs = num_pcs, features = rownames(xenium.sub))
  xenium.sub <- RunUMAP(xenium.sub, dims = 1:num_pcs, reduction.name = paste("umap.",num_pcs,"PC", sep = ""), reduction.key = paste("UMAP",num_pcs,"PC_",sep=""))
  xenium.sub <- FindNeighbors(xenium.sub, reduction = "pca", dims = 1:num_pcs)
  xenium.sub <- FindClusters(xenium.sub, resolution = 0.3)
  return(xenium.sub)
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
  data.frame(row.names = 1, check.rows = F, check.names = F) 

panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)
panc.my@images[[1]] <- NULL

cell.types.oi <- c( 'Hepatocytes', 'Cholangiocytes', 'HSCs|fibro',  'LSECs')

cell.types.in.object <- unique(as.character(unlist(panc.my[[cell_column]])))

DefaultAssay(panc.my) <- 'Xenium'
panc.my <- DietSeurat(panc.my, assays = 'Xenium')

cell.types.oi %>% walk (function(ct) {

    print(ct)
    int.sub <- subset(x = panc.my, 
                      cells = rownames(dplyr::filter(panc.my@meta.data, 
                                                     grepl(ct, .data[[cell_column]])
                      )
                      )
    )
    print(dim(int.sub))
    
    int.sub <- normalizeXenium(int.sub)
    
    cat('saving the object...\n')
    saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".rds"))
    
    DimPlot(int.sub, group.by = 'seurat_clusters', raster = F, label = TRUE)
    ggsave(glue::glue("Dimplot_seurat_clusters_{ct}.pdf"),
           height=5,width=8,useDingbats=FALSE,limitsize = FALSE)
    
    
  }
)








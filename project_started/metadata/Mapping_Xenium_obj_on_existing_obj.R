#Map Xenium object on snRNA object and transfer labels 

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))

###options###
######################
option_list = list(
  make_option(c("-r", "--input.reference"),
              type="character",
              default=NULL, 
              help="path to reference object",
              metavar="character"),
  make_option(c("-q", "--input.samples"),
              type="character",
              default=NULL, 
              help="path to list of samples and rds paths",
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
  make_option(c("-t", "--meta"),
              type="character",
              default="./", 
              help="metadata path",
              metavar="character"),
  make_option(c("-c", "--cell.column"),
              type="character",
              default="cell_type_merged", 
              help="column with cell type to transfer",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
samples.tb.path <- opt$input.samples
ref.path <- opt$input.reference
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$meta
cell_column <- opt$cell.column

ref.obj <- readRDS(ref.path)

meta <- fread(meta.path) %>%
  column_to_rownames(var = 'V1') 
print (head(meta))

meta <- meta %>%    
  select(all_of(cell_column))
ref.obj <- AddMetaData(ref.obj, meta)


DefaultAssay(ref.obj) <- 'SCT'
ref.obj <- RunUMAP(ref.obj, dims=1:30,  
                   nn.name = "weighted.nn", 
                   reduction.name = "wnn.umap", 
                   reduction.key = "wnnUMAP_", return.model = TRUE)


sample.table <- fread(samples.tb.path, header = F)

walk2(sample.table$V1, sample.table$V2, function(sample, query.path){
    
    
    dir.create(paste(out_path, sample, sep='/'), showWarnings = F)
    setwd(paste(out_path, sample, sep='/'))
    
    set.seed(666)
    query.obj <- readRDS(query.path)
    DefaultAssay(query.obj) <- 'SCT'
    
    anchors <- FindTransferAnchors(
      reference = ref.obj,
      query = query.obj,
      normalization.method = "SCT",
      reduction = 'rpca',
      reference.reduction = "pca",
      dims = 1:50
    )
    
    query.obj <- MapQuery(
      anchorset = anchors,
      query = query.obj,
      reference = ref.obj,
      refdata = list(
        celltype.l1 = cell_column),
      reference.reduction = "pca", 
      reduction.model = "wnn.umap"
    )
    
    p1 = DimPlot(query.obj, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
    p2 = DimPlot(ref.obj, reduction = "wnn.umap", group.by = cell_column, label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
    p1+p2
    ggsave(glue::glue('Dimplot_predicted.celltype.l1_ref_{sample}.pdf'), width = 12, height = 5)
    
    DimPlot(query.obj, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
    ggsave(glue::glue('Dimplot_predicted.celltype.l1_{sample}.pdf'), width = 7, height = 5)
    
    #saveRDS(query.obj, glue::glue('Mapped_object_{sample}.rds'))
    query.obj@meta.data %>% 
      dplyr::rename(cell_type = predicted.celltype.l1) %>% 
      fwrite( glue::glue('{sample}_processed_Xenium_cellTyped.meta.data'), sep = '\t', row.names = T)
  })

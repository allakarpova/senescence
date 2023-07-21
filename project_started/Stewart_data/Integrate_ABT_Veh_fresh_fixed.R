
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(glue))


select <- dplyr::select
out_path <- '/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/Stewart_ABT_fresh_fixed/'
input_path <- '/diskmnt/Projects/SenNet_primary/rds_objects/StewartLab/'
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)

samples.to.merge <- c("ABT_CD45_neg-1",  "ABT_CD45_neg-2" , "Veh_CD45_neg-1", "Veh_CD45_neg-2")

obj <- map(samples.to.merge, function(s) {
  o <- readRDS(glue('{input_path}/{s}/{s}_processed.rds'))
  o@meta.data$sample <- s
  return(o)
})

combined <- merge(x = obj[[1]], y = obj[-1], add.cell.ids=samples.to.merge)  

combined@meta.data$Preservation <- case_when(grepl("neg-2", combined$sample) ~ 'Fixed',
                                           TRUE ~ 'Fresh')
combined@meta.data$Treatment<- case_when(grepl("ABT", combined$sample) ~ 'ABT737',
                                             TRUE ~ 'Vehicle')  
  
cat('normalizing RNA\n')
DefaultAssay(combined) <- 'RNA'

cat ('Integrate regular RNA and combo RNA by batches \n')
all.rna.list <- SplitObject(combined, split.by = 'Preservation')

all.rna.list <- lapply(X = all.rna.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")

  x <- x %>% SCTransform(
    assay = 'RNA',
    vars.to.regress =  c("nCount_RNA", "percent.mt"),
    conserve.memory = T,
    return.only.var.genes = T
  )
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
int <- FindClusters(int, resolution = 2, algorithm = 4, method = "igraph")

int <- PrepSCTFindMarkers(int)

cat('saving the object...\n')
saveRDS(int,  paste0(length(samples.to.merge),"_Stewart_ABT_Veh_Fresh_Fixed_Integrated_normalized.rds"))


DimPlot(int, group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_clusters_Stewart_ABT_Veh_Fresh_Fixed_Integrated.pdf"),height=10,width=11)
DimPlot(int, group.by = "sample", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_Sample_Stewart_ABT_Veh_Fresh_Fixed_Integrated.pdf"),height=10,width=12)

DimPlot(int, group.by = "Preservation",label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_Preservation_Stewart_ABT_Veh_Fresh_Fixed_Integrated.pdf"),height=10,width=11)
DimPlot(int,  group.by = "Treatment",  label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_Treatment_Stewart_ABT_Veh_Fresh_Fixed_Integrated.pdf"),height=10,width=11)


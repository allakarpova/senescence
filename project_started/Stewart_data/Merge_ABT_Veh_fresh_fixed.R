
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
out_path <- '/diskmnt/Projects/HTAN_analysis_2/BRCA/Analyses/Alla/DCIS_project/merged/snRNA_05182023/recluster2/HT297B1_HT323B1_HT480B1/'
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

combined <- combined %>%
  NormalizeData(assay = 'RNA')

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-") 

#print(head(combined$G2M.Score))

combined <- combined %>%
  SCTransform(
    assay = 'RNA',
    vars.to.regress = c("nCount_RNA", 
                        "percent.mt" 
                        
    ),
    return.only.var.genes = FALSE, verbose = T) %>%
  RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
  RunUMAP(dims = 1:30, verbose = T) %>%
  FindNeighbors( dims = 1:30) %>%
  FindClusters(resolution = 2, verbose = FALSE,  algorithm = 4, method = "igraph")#method igraph is needed for leiden to work on large datasets



cat('saving the object...\n')
saveRDS(combined,  paste0(length(samples.to.merge),"_Stewart_ABT_Veh_Fresh_Fixed_Merged_normalized.rds"))

DimPlot(combined, group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_clusters_Stewart_ABT_Veh_Fresh_Fixed_Merged.pdf"),height=10,width=11)
DimPlot(combined, group.by = "sample", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_Sample_Stewart_ABT_Veh_Fresh_Fixed_Merged.pdf"),height=10,width=12)

DimPlot(combined, group.by = "Preservation",label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_Preservation_Stewart_ABT_Veh_Fresh_Fixed_Merged.pdf"),height=10,width=11)
DimPlot(combined,  group.by = "Treatment",  label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("Dimplot_Treatment_Stewart_ABT_Veh_Fresh_Fixed_Merged.pdf"),height=10,width=11)


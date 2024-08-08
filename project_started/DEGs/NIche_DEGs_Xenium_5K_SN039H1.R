suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(harmony)
library(readxl)
library(ggpubr)
library(viridis)

setwd('/diskmnt/Projects/SenNet_analysis/Main.analysis/old_vs_young/xenium_validation/BANKSY/transfer_to_5K')
niche5k <- fread('/diskmnt/Projects/SenNet_analysis/Main.analysis/kgallant/xenium_analysis/banksy_analysis_5k/SN039H1-Md1Fp2U2_processed_BANKSY_xenium_20240722_post_pca/SN039H1-Md1Fp2U2_processed_BANKSY_xenium_20240722_post_pca_BANKSY_xenium_20240722_final_metadata.tsv') %>%
  select(V1, BANKSY_snn_res.0.5) %>%
  rename(niche = BANKSY_snn_res.0.5) %>%
  column_to_rownames('V1')

head(niche5k)

xenium.obj <- readRDS('/diskmnt/Projects/SenNet_analysis/Main.analysis/data_freeze/v1.1/Xenium/objects/SN039H1-Md1Fp2U2_processed.rds')


xenium.obj <- AddMetaData(xenium.obj, niche5k)

Idents(xenium.obj) <- 'niche'

deg <- FindAllMarkers(xenium.obj, logfc.threshold = 0.1, min.pct = 0.001, only.pos = F)

fwrite(deg, 'Niche_DEGs_Banksy_res0.5_SN039H1_5K.pdf', sep='\t')
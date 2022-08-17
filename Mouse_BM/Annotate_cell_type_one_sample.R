library(Seurat)
library(dplyr)
library(umap)
library(ggplot2)
library(cowplot)
library(stringr)
library(data.table)
library(paletteer)
library(biomaRt)
library(readxl)
require("biomaRt")

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  #return(humanx)
  return(genesV2)
}


Caps <- function(x) {
  s <- x
  toreturn = sapply(s, function(x) paste(toupper(substring(x, 1,1)), substring(x, 2),
                                       sep="") )
  return(as.character(toreturn))
}
# process marker genes from https://www.sciencedirect.com/science/article/pii/S0092867419304593?via%3Dihub#mmc1
BM.stroma.markers <- read_excel('~/lab_Ding/work/single_cell/senescence/Papers/Mouse_BM_stroma/1-s2.0-S0092867419304593-mmc1.xlsx', skip = 1)
BM.stroma.markers <- BM.stroma.markers %>% mutate (Cell_type = case_when(BM.stroma.markers$cluster==0 ~ 'EC-sinusoidal',
                                                                         BM.stroma.markers$cluster==1 ~ 'Lepr-MSC',
                                                                         BM.stroma.markers$cluster==2 ~ 'Chondro-hyper',
                                                                         BM.stroma.markers$cluster==3 ~ 'Fibro-4',
                                                                         BM.stroma.markers$cluster==4 ~ 'Chondro-progen',
                                                                         BM.stroma.markers$cluster==5 ~ 'Fibro-5',
                                                                         BM.stroma.markers$cluster==6 ~ 'EC-arteriolar',
                                                                         BM.stroma.markers$cluster==7 ~ 'OLC-1',
                                                                         BM.stroma.markers$cluster==8 ~ 'OLC-2',
                                                                         BM.stroma.markers$cluster==9 ~ 'Fibro-1',
                                                                         BM.stroma.markers$cluster==10 ~ 'Chondro-prol.rest',
                                                                         BM.stroma.markers$cluster==11 ~ 'EC-arterial',
                                                                        BM.stroma.markers$cluster==12 ~ 'Pericytes',
                                                   BM.stroma.markers$cluster==13 ~ 'Chondro',
                                                   BM.stroma.markers$cluster==15 ~ 'Fibro-2',
                                                   BM.stroma.markers$cluster==16 ~ 'Fibro-3',
                                                   BM.stroma.markers$cluster==17 ~ 'Chondro-prehyper-2'))

top10.stroma.markers <- BM.stroma.markers %>% filter(p_val_adj < 0.05) %>% group_by(Cell_type) %>% top_n (n = 20, wt = avg_logFC)
fwrite(top10.stroma.markers, '~/lab_Ding/work/single_cell/senescence/Papers/Mouse_BM_stroma/Top20_up_genes_stroma_BM_cells.txt', sep = '\t')

#process markers from https://www.nature.com/articles/s41556-019-0439-6#Sec27
BM.cells.markers <-  read_excel('~/lab_Ding/work/single_cell/senescence/Papers/Mouse_BM_single_cell/41556_2019_439_MOESM3_ESM.xlsx', sheet = 2)[-1,]
BM.cells.markers <- BM.cells.markers[1:20,]
top20.cell.markers <- reshape2::melt(as.matrix(BM.cells.markers))[2:3] %>% data.frame() 
colnames(top20.cell.markers) <- c('Cell_type', 'gene')
fwrite(top20.cell.markers, '~/lab_Ding/work/single_cell/senescence/Papers/Mouse_BM_single_cell/Top20_up_genes_all_cells_BM.txt', sep = '\t')


panc = readRDS('~/lab_Ding/work/single_cell/senescence/objects/Second_DOXO/Second_DOXO_processed.rds')
add_filename <- 'Second_DOXO'
wd <- '~/lab_Ding/work/single_cell/senescence/mouse_BM/cell_type_annotation/Second_DOXO'
dir.create(wd, recursive = T)
setwd(wd)

DimPlot(panc, label = T)

# plot dotplot markers
myeloid.genes <- fread('~/lab_Ding/work/single_cell/Cell_state_markers_v12182020.txt', data.table = F, header = T)
myeloid.genes$Gene <- myeloid.genes$Gene %>% tolower() %>% Caps()

genes2plot <- myeloid.genes$Gene %>% unique()
genes2plot <- 'Fn1'
DotPlot(object = panc, features = 'Fn1', )
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'RNA')
p
p$data <- merge(p$data, myeloid.genes[1:2], by.x = 'features.plot', by.y = 'Gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Gene_set , scales = "free",  drop = T, ncol = 7)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0( "Dotplot_marker_gene_expression_", add_filename, "_SCT.pdf"),height=50,width=50,useDingbats=FALSE,limitsize = FALSE)

#####
genes2plot <- top10.stroma.markers$gene %>% unique()
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'SCT')

p$data <- merge(p$data, top10.stroma.markers[7:8], by.x = 'features.plot', by.y = 'gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell_type , scales = "free",  drop = T, ncol = 6)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0( "Dotplot_20stroma_gene_expression_", add_filename, "_SCT.pdf"),height=15,width=40,useDingbats=FALSE,limitsize = FALSE)



#####
genes2plot <- top20.cell.markers$gene %>% unique()
p <- DotPlot(object = panc, group.by = 'seurat_clusters', features = genes2plot, col.min = 0, dot.min = 0.01, assay = 'SCT')

p$data <- merge(p$data, top20.cell.markers, by.x = 'features.plot', by.y = 'gene', all.x=T)
p <- p  + RotatedAxis()
p <- p + facet_wrap(~ Cell_type , scales = "free",  drop = T, ncol = 6)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5, size = 14),
               axis.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
               strip.placement = "outside")
p <- p + scale_color_paletteer_c("viridis::viridis", direction = 1)
p
ggsave(paste0( "Dotplot_20cell.type_gene_expression_", add_filename, "_SCT.pdf"),height=30,width=40,useDingbats=FALSE,limitsize = FALSE)


markers <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25, assay = 'SCT')
fwrite(markers, 'Markers.cluster.txt', sep ='\t')

panc$Cell_type <- case_when( panc$seurat_clusters %in% c('0') ~ 'Fibro-4.5',
                             panc$seurat_clusters %in% c('1') ~ 'Chondrocytes-progenitor',
                             panc$seurat_clusters %in% c('2') ~ 'Lepr−MSC.Adipo−CAR',
                             panc$seurat_clusters %in% c('3') ~ 'Chondrocytes-pre-/hypertrophic ',
                              panc$seurat_clusters %in% c('4') ~ 'Neutrophil/Monocyte progenitors',
                             panc$seurat_clusters %in% c( '5') ~ 'Erythroid/Megakaryocyte progenitors',
                             
                             panc$seurat_clusters %in% c('7') ~ 'Fibro-1',
                             
                             
                             panc$seurat_clusters %in% c('10') ~ 'Endothelial cells',
                              panc$seurat_clusters %in% c('11') ~ 'Fibro-2',
                            panc$seurat_clusters %in% c('12') ~ 'Osteoblasts',
                            panc$seurat_clusters %in% c('13') ~ 'Fibro-3',
                            panc$seurat_clusters %in% c( '14') ~ 'Eosinophil/Basophil progenitors')
                            
                                
                                 
                                 
                                 
                                 
                                 
                                  
                                  
                                 

DotPlot(object = panc, features = 'Mgp', col.min = 0, assay = 'SCT')

FeaturePlot(panc, 'Mgp', assay = 'SCT')
#######GENE LISTS###MYELOMA SPECIFIC####
plasma_only_candidates<-c("CD79A","CD79B","FCRL5","FCRL2","FCER2","CD38","GPR160","LAMP3","TNFRSF13B","RASGRP3","TNFRSF17","SLAMF7","SDC1")
erythrocyte_features<-c("HBB","HBA1","CA1","AHSP","HBA2","BLVRB","TFRC","GYPA","SLC25A37","HBM","HBD")
tcell_features <- c("MALAT1","LTB","IL32","IL7R","CD2","B2M","ACAP1","CD27","STK17A","CTSW","LDHB","CD3G","CD3E","CD3E","NOSIP")
nk_features <- c("GNLY","TYROBP","FCGR3A","KLRD1","KLRC1","HOXP","B2M","ACAP1","STK17A","CTSW")
DC_features<-c("FCER1A","CD1C","ENHO","CLEC10A","PKIB")
pDC_features<-c("CLEC4C","LILRA4") #https://www.sciencedirect.com/science/article/pii/S0022175918303387?via%3Dihub
CD34CYTL1_features<-c("CD34","CYTL1")
Plasma_features <- c("SDC1", "IGHG1", "IGHG3", "IGHG4","MZB1","TNFRSF17","TNFRSF13B","FCRL5","IGKC","SSR4","DERL3","SEC11C","IGLC2","IGHA1","IGHA2")
B_features <- c("CD79A", "CD79B","CD19", "MS4A1","IGHM","RALGPS2","HLA-DRA","HLA-DPB1","HLA-DQA1","HLA-DRB1","HLA-DPA1","CD74")
#HLA genes are also shared with other tissue types but help specific plasma and  B cell  groups
Monocyte_macrophage_features <-c("LYZ", "CD14", "S100A8","S100A9","FCN1","CST3","VCAN","AIF1","FCGR3A","TYROBP","LGALS1","MS4A7", "IFITM3","CD68","CSF1R")
#newbies
Treg<-c("CD4","IL2RA","FOXP3","CTLA4","TIGIT","TNFRSF4","LAG3","PDCD1","HAVCR2")
Tnaive<-c("CD4","CD8A","CD8B","CCR7","SELL","LEF1")
CD8T<-c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B","EVL","CD69","CCL5","ARHGDIB","CD7","IL2RB","RAC2","PTPRC")
CD4T<-c("CD3D", "CD3E", "CD3G", "CD6","IL7R", "LDHB", "NOSIP", "CD4","TBX21","PTPRC")
MacroM2<-c("IL4", "CD163", "SRSF10","MRC1","CD200R1","TGM2","TNFRSF6B","IL1R2")
MacroM1<-c("HLA-DRA","ICAM1", "LBP","CD80","CD86","CD68","IL1R1","TLR2","TLR4","NOS2","SOCS2")
Neutrophils<-c("AZU1","MPO","ELANE","CTSG","LYZ")
Monocyte<-c("LYZ", "CD14", "S100A8")
Macrophage<-c("FCGR3A", "MS4A7", "IFITM3","CD68","CSF1R")
Megakaryocyte<-c("PPBP","PF4","NRGN","RGS18","GNG11","HIST1H2AC","CAVIN2","TAGLN2","PLEK")#https://www.biorxiv.org/content/10.1101/127217v1.full.pdf
Pre_B<-c("IGLL1","VPREB1","CD24","STNM1","HIST1H1C","SOX4","MKI67","CD79B","IGHM","IRF4","CD38")
HPC<-c("SPINK2","AREG","KIAA0125","ANKRD28","PRSS57","TSC22D1","SOX4","STMN1","IGLL1","KIAA0101")
Immature_neutrophil<-c("MPO","PRTN3","ELANE","AZU1","LYZ","RNASE2","CTSG","PRSS57","RETN")
mDC<-c("HLA-DPB1","FCER1A","HLA-DPA1","HLA-DRA","CST3","HLA-DQA1","HLA-DRB1","HLA-DQB1","CD74","HLA-DRB5")
NK_strong<-c("NCAM1","FCGR3A","SPON2","KLRF1","CCL4","CCL5")
NK_weak<-c("NCAM1","GZMK","CCR7")
T_markers<-c("CD4","CD8A","CD8B", #base 
             "CD28","CCR7","SELL",#TCM
             "PRF1","GZMK","GZMB","TNF","FGFBP2","CX3CR1","FCGR3A", #CD8 effector + INFG 
             "MKI67","IFNG","CCL4","IL2RA","TBX21", #effector memory Th1
             #"IL4","IL5","IL13","GATA3", #effector memory Th2 - These data are indicative of a prevalent but slightly skewed TH2 response, consistent with the previous report that IL13 rather than IL4 is the dominant TH2 cytokine observed in CD19-BB-3z CAR-T cell activation [13]. Thus, for further analysis, IL13 rather than IL4 was chosen as the TH2 dominant cytokine gene.
             #"FOXP3","IL10", #Treg
             "PDCD1","LAG3","HAVCR2","EOMES",#exhausted
             "CTLA4","ICOS","CSF2", #upregulated in other studies
             "CD27","IL7R","KLRB1","CCR4","SLC4A10"
)


# int_st <- int_st %>%
#   SCTransform(
#     assay = 'Spatial', 
#     vars.to.regress = c("nCount_Spatial", "percent.mito"),
#     #variable.features.n = opt$top_variable_features,
#     return.only.var.genes = FALSE
#   ) %>%
#   RunPCA(assay = 'SCT', do.print = FALSE) %>%
#   RunUMAP(dims = 1:30, reduction = "pca") %>%
#   FindNeighbors(dims = 1:30, reduction = "pca") %>%
#   FindClusters(resolution = 0.8)
# 
# saveRDS(int_st, 'obj/NMK_3W_52W_MF_integrated.rds')

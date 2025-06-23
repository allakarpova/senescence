# written by Austin Southard-Smith
# 2023-10-11
# based on https://satijalab.org/seurat/articles/spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
# addresses error mentioned https://github.com/satijalab/seurat/issues/7491
library(Seurat)
library(tidyverse)
library(future)
library(optparse)

set.seed(1234)
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="path to data folder (e.g. xenium output folder: /diskmnt/primary/Xenium/data/20230803__195918__20230803_firstsamples/output-XETG00122__0010660__HT270P1-S1H1A1Fp1Us1_1__20230803__200008/)",
                metavar="character"),
    make_option(c('-s', "--sampleID"),
                type="character",
                default=NULL,
                help="SampleID that output file names will be based on",
                metavar="character"),
    make_option(c("-o", "--output"),
                type="character",
                default="./",
                help="output folder path",
                metavar="character"),
    make_option(c("--mols.qv.threshold"),
                type="integer",
                default=20,
                help="threshold for the qv score cutoff threshold. Transcripts above this will be included.",
                metavar="integer"),
    make_option(c("--pc_num"),
                type="integer",
                default=30,
                help="number of principal components to use",
                metavar="integer"),
    make_option(c("--counts_cutoff"),
                type="integer",
                default=0,
                help="retain only cells with nCounts greater than the value specified by the counts_cutoff",
                metavar="integer"),
    make_option(c("--plot_features"),
                type='character',
                default=NA,
                help = "comma separated list of features to plot when the object is generated (for pancreas this looks like 'KRT7,AMY2A,TFF2,CFTR,EPCAM,MKI67,TOP2A,PPARG,PECAM1,VWF,MS4A1,CD3D,CD3E,CD4,MYC,INS,GCG,C7,PDGFRA,SRPX,MAMDC2,VCAN,MEDAG,MFAP5,GPRC5A,CAV1,PTPRC,IL7R,TMC5,MET,EHF,CAPN8,CD163,MS4A6A,PDGFRB,MYH11,KCNMA1'",
                metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
    print_help(opt_parser)
    stop("Path to data is required (--input).n", call.=FALSE)
}
if (is.null(opt$sampleID)){
    print_help(opt_parser)
    stop("sampleID is required (--sampleID).n", call.=FALSE)
}

input = opt$input
sampleID = opt$sampleID
out_path = paste0(opt$output,"/",sampleID,"/")
dir.create(out_path)
num_pcs = opt$pc_num

plan("multicore", workers = 30)
# The below errors out because it is reading in the Cell IDs incorrectly per https://github.com/satijalab/seurat/issues/7491
# xenium.obj <- LoadXenium(data.dir=input, fov = "fov", assay = "Xenium") #this might be fixed in Seurat5.
xenium.list <- ReadXenium(data.dir=input, outs = c("matrix", "microns"), type = c("centroids", "segmentations"), mols.qv.threshold = opt$mols.qv.threshold)
plan("sequential")
assay <- "Xenium"
segmentations.data <- list(
    "centroids" = CreateCentroids(xenium.list$centroids),
    "segmentation" = CreateSegmentation(xenium.list$segmentations)
)
coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = xenium.list$microns,
    assay = assay
)
xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay)
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Probe"]])
fov <- "fov" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
xenium.obj[[fov]] <- coords
#check to make sure that plotting with transcript/cell coordinates works
pdf(paste(out_path,"VlnPlot_QC_",sampleID,".pdf", sep=''),useDingbat=FALSE, width=8, height=6)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()
# pdf(paste("ImageDimPlot_QC_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #This shows all of the dots on one plot
# ImageDimPlot(xenium.obj, fov = "fov", molecules = c("TFF2", "KRT7", "VCAN", "CFTR"), nmols = 20000)
# dev.off()
# pdf(paste("ImageDimPlot_QC_features_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# ImageFeaturePlot(xenium.obj, features = c("TFF2", "KRT7", "VCAN", "CFTR"), max.cutoff = c(20, 20, 15, 5), size = 0.25, cols = c("lightgrey", "blue", "yellow"), border.size= NA)
# dev.off()
# pdf(paste("ImageDimPlot_QC_features_2color_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# ImageFeaturePlot(xenium.obj, features = c("TFF2", "KRT7", "VCAN", "CFTR"), max.cutoff = c(20, 20, 15, 5), size = 0.25, cols = c("blue", "yellow"), border.size= NA)
# dev.off()
# save the raw
saveRDS(xenium.obj, paste(out_path,sampleID,"_raw.rds",sep=''))
# some cells will be empty so you need to remove them
counts <- xenium.obj@assays$Xenium@counts
empty_cells <- colnames(counts)[(colSums(counts) == 0)]
write.table(empty_cells,paste(out_path,"Empty_cells_",sampleID,'.tsv',sep=''),sep="\t",quote=FALSE)
cells_to_keep <- colnames(counts)[(colSums(counts) > opt$counts_cutoff)]
xenium.sub <- subset(xenium.obj, cells = cells_to_keep)
# normalize,umap,cluster, SCTv2 handles low number of genes very well compared to LogNormalize (see Ilya's results)
xenium.sub <- SCTransform(xenium.sub, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2")
xenium.sub <- RunPCA(xenium.sub, npcs = num_pcs, features = rownames(xenium.sub))
xenium.sub <- RunUMAP(xenium.sub, dims = 1:num_pcs, reduction.name = paste("umap.",num_pcs,"PC", sep = ""), reduction.key = paste("UMAP",num_pcs,"PC_",sep=""))
xenium.sub <- FindNeighbors(xenium.sub, reduction = "pca", dims = 1:num_pcs)
xenium.sub <- FindClusters(xenium.sub, resolution = 0.3)
DefaultFOV(xenium.sub, assay='SCT') <- 'fov'
# the above line results in the following warning and may cause problems down the line. to do this properly you would copy the current 'fov' slot and save it to the FOV slot for the SCT assay I think. Not 100% sure how to do this yet.
# Warning: FOV 'fov' currently associated with assay 'Xenium', changing to 'SCT'

# save the processed object.
Idents(xenium.sub) <- "seurat_clusters"
# saveRDS(xenium.sub, file='HT270P1-S1H1A1Fp1Us1_1_processed.rds')
# adding in cell labels (here they are based on the xenium browser clusters)
# can't use the traditional method for seurat cluster adding and I think it has to do with the xenium browser clusters not being a factor.
xenium_browser_clusters <- read.table(paste(input,"/analysis/clustering/gene_expression_graphclust/clusters.csv",sep=''), sep = ",", header=T)
xenium_browser_clusters <- xenium_browser_clusters[(xenium_browser_clusters$Barcode %in% cells_to_keep),]
rownames(xenium_browser_clusters) <- xenium_browser_clusters$Barcode
xenium_browser_clusters$Barcode <- NULL
xenium.sub <- AddMetaData(object = xenium.sub, metadata = xenium_browser_clusters, col.name = "xenium_browser_graph_based_clusters")
# adding in cell labels (here they are based on the xenium browser clusters)
# can't use the traditional method for seurat cluster adding and I think it has to do with the xenium browser clusters not being a factor.
# could also be due to the xenium object not containing the right cell type labels.
#cell_types <- c("your","clusters","go","here")
#xenium.sub$cell_type_individual <- NA
#xenium_clusters = sort(unique(xenium.sub$xenium_browser_graph_based_clusters))
#for (i in 1:length(xenium_clusters)) {
#    value = xenium_clusters[i]
#    cell_type = cell_types[i]
#    xenium.sub$cell_type_individual[xenium.sub$xenium_browser_graph_based_clusters == value] <- cell_type
#}
# Setting the default values to be included in plots
DefaultAssay(xenium.sub) <- "SCT"
Idents(xenium.sub) <- "seurat_clusters"
# save the processed object.
saveRDS(xenium.sub, paste(out_path,sampleID,"_processed.rds",sep=''))

print("generating_QC_plots")
DefaultAssay(xenium.sub) <- "Xenium"
DefaultFOV(xenium.sub, assay='Xenium') <- 'fov'
pdf(paste(out_path,"DimPlots_QC_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6)
DimPlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), group.by = "seurat_clusters", label=TRUE, label.size=4, raster=TRUE)
DimPlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), group.by = "xenium_browser_graph_based_clusters", label=TRUE, label.size=4, raster=FALSE)
FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = "nCount_Xenium", label=TRUE, label.size=4, raster=TRUE)
FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = "nFeature_Xenium", label=TRUE, label.size=4, raster=TRUE)
ImageFeaturePlot(xenium.sub, features = "nCount_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("blue", "yellow"), border.size= NA) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
ImageFeaturePlot(xenium.sub, features = "nFeature_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("blue", "yellow"), border.size= NA) #max.cutoff is 90th percentile
ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
ImageDimPlot(xenium.sub, fov = "fov", group.by="xenium_browser_graph_based_clusters", size = 0.15, border.size=NA)
dev.off()
DefaultAssay(xenium.sub) <- "SCT"
DefaultFOV(xenium.sub, assay='SCT') <- 'fov' 
if (is.na(opt$plot_features)){
    print("skipping plotting of features as no features were provided")
} else {
    genes_filtered <-c()
    plot_features=str_split(opt$plot_features, pattern=',')[[1]]
    plot_features=unique(plot_features)
    sample_genes <- rownames(x = xenium.sub@assays$SCT@counts)
    
    for (gene in plot_features) {
        if (gene %in% sample_genes) {
            genes_filtered <- c(genes_filtered,gene)
        }
    }
    print("the following genes are unique and are present in the SCT assay for inclusion in plots")
    print(genes_filtered)
    pdf(paste(out_path,"FeaturePlots_markers_of_interest_SCT_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6)
    for (marker in genes_filtered) {
        print(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = marker, label=TRUE, label.size=4, raster=TRUE))
    }
    dev.off()
    pdf(paste(out_path,"ImageDimPlot_QC_features_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
    for (marker in genes_filtered){
        print(ImageFeaturePlot(xenium.sub, features = marker, max.cutoff = 'q95', size = 0.15, cols = c("lightgrey", "blue", "yellow"), border.size= NA)) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
        print(ImageFeaturePlot(xenium.sub, features = marker, max.cutoff = 'q90', size = 0.15, cols = c("blue", "yellow"), border.size= NA)) #max.cutoff is 90th percentile
    }
    dev.off()
}


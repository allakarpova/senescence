#run several rounds of clustering with different resolutions
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(data.table))

suppressMessages(library(optparse))


###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to object",
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
  make_option(c("--multiome"),
              type="logical",
              default=FALSE,
              help = "is it a multiome object or not?",
              metavar="logical")
  
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

select <- dplyr::select
filter <- dplyr::filter

# read in initial arguments
input.path<- opt$input.object
out_path <- opt$output
add_filename <- opt$extra

dir.create(out_path, showWarnings = F)
setwd(out_path)

panc<- readRDS(input.path)

resolutions <- c(seq(0.1, 0.9, 0.1), seq(1, 1.8, 0.2), 2)
algs <- c(4,1,3)

markers.from.email <- c('GLB1', 'CDKN1A','CDKN2A', 'TP53', 'SERPINE1','FOS', 'JUN', 'IL1A', 'IL6', 'CXCL8', 'CCL2', 'MMP1', 'MMP3', 'IGFBP7')
more.markers <- c('CDKN2D', 'CENPB', 'JUN', 'LMNB1',  'TNF', 'TGFB1', 'CXCL1', 'NFKB1', 'RELA', 'HMGB1', 'SERPINB2',
                  'INHBA', 'GDF15', 'CCL11', 'MIF', 'CXCL2', 'IGFBP3', 'CCL24', 'MMP12', 'CXCL10', 'TP53BP1', 'BCL2', 'MKI67', 'TOP2A')
dotplot.color <- colorRampPalette(c('#eae2b7','#f3d180', '#fcbf49','#f77f00','#d62828','#6b2c39', '#003049'))(10) 

for (resol in resolutions) {
  if(opt$multiome) {
    panc <- FindClusters(panc, resolution = resol,graph.name = "wsnn", algorithm = 4, method = "igraph")
  } else {
    panc <- FindClusters(panc, resolution = resol, algorithm = 4, method = "igraph")
  }
  
  DotPlot(panc, features = c(markers.from.email, more.markers), cluster.idents = TRUE,
          group.by = 'seurat_clusters') + scale_color_gradientn(colors = dotplot.color) +
    theme(axis.text.x = element_text(angle =90, hjust=1, vjust=0.5))+
    scale_size_area(limits=c(0,20), oob=scales::squish)
  ggsave(paste0('Dotplot_', 'known_and_extra_markers','_by_seurat_clusters_',resol,'_alg4.pdf'), 
         width = length(c(markers.from.email, more.markers))*0.2 + 2, height = 5)
  
  
  print(head(panc@meta.data))
}

cluster.tb <- panc@meta.data %>% select(dplyr::contains('res.'))
fwrite(cluster.tb, paste0('Clusters_res0.1_to_2_alg4_', add_filename, '.txt'), sep='\t', row.names = TRUE)


for (resol in resolutions) {
  if(opt$multiome) {
    panc <- FindClusters(panc, resolution = resol,graph.name = "wsnn", algorithm = 1)
  } else {
    panc <- FindClusters(panc, resolution = resol, algorithm = 1)
  }
  DotPlot(panc, features = c(markers.from.email, more.markers), cluster.idents = TRUE,
          group.by = 'seurat_clusters') + scale_color_gradientn(colors = dotplot.color) +
    theme(axis.text.x = element_text(angle =90, hjust=1, vjust=0.5))+
    scale_size_area(limits=c(0,20), oob=scales::squish)
  ggsave(paste0('Dotplot_', 'known_and_extra_markers','_by_seurat_clusters_',resol,'_alg1.pdf'), 
         width = length(c(markers.from.email, more.markers))*0.2 + 2, height = 5)
  print(head(panc@meta.data))
}

cluster.tb <- panc@meta.data %>% select(dplyr::contains('res.'))
fwrite(cluster.tb, paste0('Clusters_res0.1_to_2_alg1_', add_filename, '.txt'), sep='\t', row.names = TRUE)

for (resol in resolutions) {
  if(opt$multiome) {
    panc <- FindClusters(panc, resolution = resol,graph.name = "wsnn", algorithm = 3)
  } else {
    panc <- FindClusters(panc, resolution = resol, algorithm = 3)
  }
  
  DotPlot(panc, features = c(markers.from.email, more.markers), cluster.idents = TRUE,
          group.by = 'seurat_clusters') + scale_color_gradientn(colors = dotplot.color) +
    theme(axis.text.x = element_text(angle =90, hjust=1, vjust=0.5))+
    scale_size_area(limits=c(0,20), oob=scales::squish)
  ggsave(paste0('Dotplot_', 'known_and_extra_markers','_by_seurat_clusters_',resol,'_alg3.pdf'), 
         width = length(c(markers.from.email, more.markers))*0.2 + 2, height = 5)
  print(head(panc@meta.data))
}
cluster.tb <- panc@meta.data %>% select(dplyr::contains('res.'))
fwrite(cluster.tb, paste0('Clusters_res0.1_to_2_alg3_', add_filename, '.txt'), sep='\t', row.names = TRUE)



library(clustree)
clustree.dir <- out_path

reduct <- ifelse(opt$multiome, 'wnn.umap', 'umap.harmony')
algs %>% walk( function(a) {
  rna.clust <- fread(paste0(clustree.dir,'/Clusters_res0.1_to_2_alg',a,'_',add_filename,'.txt')) %>%
    data.frame(row.names = 1)
  p2 <- clustree(rna.clust, prefix = "SCT_snn_res.",node_text_size = 5, #node_colour = "sc3_stability", 
                 layout = "sugiyama")
  
  pdf(paste0("All_clusters_plot_alg",a,".pdf"), width = 17, height = 10)
  print(p2)
  dev.off()
  panc <- AddMetaData(panc, rna.clust)
  
  resolutions %>% walk(function( res) {
    DimPlot(panc, group.by = paste0('SCT_snn_res.',res), reduction = reduct, label =TRUE)
    ggsave(paste0("Dimplot_T-cell_",res,"_alg",a,".pdf"),height=5,width=6,useDingbats=FALSE)
  })
})




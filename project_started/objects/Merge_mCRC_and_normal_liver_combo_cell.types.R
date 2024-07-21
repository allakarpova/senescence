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

suppressMessages(library(EnsDb.Hsapiens.v100))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(harmony))
suppressMessages(library(googlesheets4))
gs4_deauth()

library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(chromVAR)

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

iterative_removal <- function(all_peaks.f) {
  #print(cancer.type)
  #all_peaks.f <- all_peaks.f[order(score.norm, decreasing = T), ]
  #just load existing peaks if any
  if (file.exists(paste0('peaks/',length(samples.id),'_recentered_final.',add_filename,'.tsv'))) {
    recentered_final.f <- fread(paste0('peaks/',length(samples.id),'_recentered_final.',add_filename,'.tsv'))
  } else {
    
    recentered_p=StringToGRanges(regions = all_peaks.f$new_peak, sep = c("-", "-"))
    
    cat(paste0('finding overlapping peaks \n'))
    overlapping=as.data.table(x = findOverlaps(query = recentered_p, 
                                               subject = recentered_p)) # find which peaks overlap
    overlapping=overlapping[queryHits!=subjectHits,]
    overlapping.peak.number <- unique(x = overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks.f table
    recentered_non_overlapping=all_peaks.f[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
    fwrite(recentered_non_overlapping,paste0('peaks/',length(samples.id),'_recentered_nonOverlapping.',add_filename,'.tsv'),
           sep='\t',row.names=FALSE)
    if (length(overlapping.peak.number)>0) {
      tmp <- data.table(chr = all_peaks.f$seqnames[overlapping.peak.number], 
                        num = overlapping.peak.number)
      overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) #split peaks by chromosome 
      registerDoParallel(cores=25)
      #this is where iterative removal of peaks is done
      best_in_overlapping_num <- foreach(peak.numbers=overlapping.peak.number.split) %dopar% {
        cat('removing overlapping peaks in each chromosome\n')
        iterative_removal_core (peak.numbers = peak.numbers, overlapping.f = overlapping)
      }
      stopImplicitCluster()
      best_in_overlapping_num <- do.call('c', best_in_overlapping_num) #combine best peak numbers from all chromosomes
      best_in_overlapping_cancer <- all_peaks.f[best_in_overlapping_num,] #extract peaks themselves
      fwrite(best_in_overlapping_cancer,paste0('peaks/',length(samples.id),'_recentered_Overlapping.',add_filename,'.tsv'),
             sep='\t',row.names=FALSE)
      recentered_final.f=rbindlist(list(recentered_non_overlapping,best_in_overlapping_cancer))
    } else {
      recentered_final.f=recentered_non_overlapping
    }
    final.overlaps <-  recentered_final.f$new_peak %>% 
      unique %>% 
      StringToGRanges %>% 
      countOverlaps
    if (sum(final.overlaps>1)>0) {
      stop("Execution stopped. Overlapping peaks remained")
    }
    
  }
  return(recentered_final.f)
}

# this works like a charm
iterative_removal_core <- function(peak.numbers, overlapping.f) {
  chr = peak.numbers$chr[1]
  running.vector <- peak.numbers$num
  peaks.to.trash <- NULL
  peaks.to.keep <- NULL
  while (length(running.vector) != 0) {
    n <- running.vector[1] # this is the first and the best peak since peaks are sorted by scores
    neighbor.peaks.num.discard <- overlapping.f[queryHits==n, subjectHits] #find positions of other peaks overlapping with the first one 
    running.vector <- setdiff(running.vector, neighbor.peaks.num.discard) # remove them from the list of peaks
    running.vector <- setdiff(running.vector, n)
    peaks.to.keep <- c(peaks.to.keep, n) # add this peak to the keeping list
    peaks.to.trash <- unique(c(peaks.to.trash, neighbor.peaks.num.discard)) # add neighbors to the list of peaks to discard
  }
  cat('done\n')
  return(peaks.to.keep)
}

getFeatureMatrix <- function (obj, peaks, pro_n, assay.towork.f) {
  frag <- Fragments(obj[[assay.towork.f]])
  cat('Making a large count matrix...\n')
  matrix.counts <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    process_n = pro_n,
    sep = c("-","-"),
    cells = colnames(obj)
  )
  return(matrix.counts)
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
    filter(datatype=='combo', Cell_type=='All')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v100,  standard.chromosomes = TRUE)
ct='All'

if (!file.exists(paste0(add_filename,"_",make.names(ct), ".rds"))) {
  l.path=obj.paths$liver
  m.path=obj.paths$mCRC
  obj1 <- readRDS(l.path)
  obj2 <- readRDS(m.path)
  
  obj1.peaks <- fread('/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/merge_liver_combo_9/peaks/34_recentered_final.liver_combo.tsv', header = TRUE)
  obj2.peaks <- fread('/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/merge_mCRC_combo_2/peaks/26_recentered_final.liver_combo.tsv', header = TRUE)
  
  obj1@meta.data$Cohort <- 'Normal liver'
  obj2@meta.data$Cohort <- 'mCRC liver'
  DefaultAssay(obj1) <- 'RNA'
  obj1 <- DietSeurat(obj1, assays = c('RNA', 'ATAC_merged'))
  DefaultAssay(obj2) <- 'RNA'
  obj2 <- DietSeurat(obj2, assays = c('RNA', 'ATAC_merged'))
  
  obj <- list(obj1, obj2)
  
  DefaultAssay(obj2) <- 'ATAC_merged'
  peaks.tokeep <- AccessiblePeaks(obj2, min.cells = 10)
  dim(obj2.peaks)
  obj2.peaks.filt <- filter(obj2.peaks, new_peak %in% peaks.tokeep) # trying to remove cancer specific peaks
  dim(obj2.peaks.filt)
  
  dir.create('peaks')
  
  samples.id <- c(unique(obj1$Sample), unique(obj2$Sample))
  print(samples.id)
  
  all_peaks <- rbindlist(list(obj1.peaks, obj2.peaks.filt))
  all_peaks <- all_peaks[order(score.norm, decreasing = T), ] # order peaks by normalized scores, this is essential for filtering out overlapping peaks
  fwrite(all_peaks, paste0('peaks/',length(samples.id),"_sample_MACS2_peaks_",add_filename,".tsv"),
         sep='\t',row.names=FALSE)
  
  recentered_final <- iterative_removal(all_peaks)
  
  fwrite(recentered_final, paste0('peaks/',length(samples.id),'_recentered_final.',add_filename,'.tsv'),sep='\t',
         row.names=FALSE)
  recentered_p=StringToGRanges(unique(recentered_final$new_peak), sep = c("-", "-"))
  
  plan("multicore", workers = 20)
  options(future.globals.maxSize = 250 * 1024^3) # for 250 Gb RAM
  
  peak.number <- length(unique(recentered_final$new_peak))
  n.peaks <- round(peak.number/20)
  
  matrix.counts <- map (obj, ~getFeatureMatrix(.x, recentered_p, n.peaks, 'ATAC_merged'))
  
  registerDoParallel(cores=10)
  cat ('creating assays on common peaks\n')
  obj <- foreach (obj = obj, co = matrix.counts, .combine=c) %dopar% {
    DefaultAssay(obj) <- 'RNA'
    frag = Fragments(obj[['ATAC_merged']])
    obj[['ATAC_merged']] <- NULL
    obj[['ATAC_merged']] <- CreateChromatinAssay(counts = co,
                                                 fragments=frag, 
                                                 min.cells = -1, 
                                                 min.features = -1)
    
    return(obj)
  }
  stopImplicitCluster()
  
  
  int.sub <- merge(obj[[1]], obj[[2]])
  
  # add metadata
  ifelse('passed_filters' %in% colnames(int.sub@meta.data),
         total_fragments_cell <- int.sub$passed_filters, 
         total_fragments_cell <- int.sub$atac_fragments)
  
  peak.counts <- colSums(x = GetAssayData(int.sub, slot = 'counts'))
  frip <- peak.counts *100 / total_fragments_cell
  int.sub <- AddMetaData(object = int.sub, metadata = frip, col.name = 'pct_read_in_peaks_ATAC_merged')
  int.sub <- AddMetaData(object = int.sub, metadata = peak.counts, col.name = 'peak_region_fragments_ATAC_merged')
  
  
  print(dim(int.sub))
  int.sub <- runAllNormalization(int.sub, assay =  'ATAC_merged')
  
  
  saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".rds"))
  
} else {
  int.sub <- readRDS( paste0(add_filename,"_",make.names(ct), ".rds"))
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
#})

}







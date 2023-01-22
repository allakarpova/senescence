suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
set.seed(1234)
suppressMessages(library(future))
plan("multicore", workers =10)
options(future.globals.maxSize = 50 * 1024^3)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(future))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(AUCell))
library(msigdbr)

Run_AUCell <- function(obj, geneset.list, name) {
  
    DefaultAssay(obj) <- 'SCT'
    exprMatrix <- GetAssayData(obj, slot = 'data')
    
    cells_rankings <- AUCell_buildRankings(exprMatrix)
    cells_AUC <- AUCell_calcAUC(geneset.list, cells_rankings,  nCores=20)
    saveRDS(cells_AUC, file=glue::glue("{name}_cells_AUC.rds"))
    
    pdf(glue::glue('{name}_AUCell_histograms.pdf'))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
    print(cells_assignment)
    dev.off()
    
    saveRDS(cells_assignment, file=glue::glue("{name}_cells_assignment.rds"))
    
    obj[["AUCell"]] <- CreateAssayObject(getAUC(cells_AUC))
    
    pdf(glue::glue('{name}_AUCell_featureplots.pdf'), onefile = TRUE)
    print(DimPlot(obj, group.by = 'cell_type_merged', cols = 'Dark2', label = TRUE))
    names(geneSets) %>% map(function(f) {
      tryCatch({
        print(suppressWarnings(FeaturePlot(obj, features = f)))
      }, error = function(e) {
        message(e)
        NULL
      })
    })
    dev.off()
    dev.off()
    
}
  


Run_correlation <- function(obj, s, ct, dot.size = 1, cor.adjust.method = 'fdr') {
  print(s)
  #s = 'Fibroblast_replicative_sen'
  genes.in.signature <- subset(all.signatures, ont == s)$Gene_names
  # #get genes with good coverage
  # counts <- GetAssayData(object = panc, slot = "counts")
  # keep_genes <- filter_counts(counts, cell.pct = cell.pct.cutoff)
  # genes.in.signature.good.coverage <- intersect(genes.in.signature, keep_genes)
  cat(paste0('Genes from a signature total ', length(genes.in.signature), '\n'))
  #cat(paste0('Genes from a signature with good coverage ', length(genes.in.signature.good.coverage), '\n'))
  #get expression data and scale within each cell
  expr.mat <- FetchData(obj, vars = genes.in.signature, slot = 'data') %>% t() %>% scale()
  #dim(expr.mat)
  fc.sign <- subset(all.signatures, ont == s)
  rownames(fc.sign) <- fc.sign$Gene_names
  fc.sign <- fc.sign[genes.in.signature,]
  #dim(fc.sign)
  if(dim(expr.mat)[1] != dim(fc.sign)[1]) {
    fc.sign <- subset(fc.sign, Gene_names %in% rownames(expr.mat))
    stopifnot(fc.sign$Gene_names== rownames(expr.mat))
  }
  
  cor.result <- psych::corr.test(expr.mat, fc.sign$log2FC, method = 'spearman', adjust = cor.adjust.method)
  df <- cbind('rho' = cor.result$r , 'p.val' = cor.result$p) %>% data.frame(check.rows = F, check.names = F)
  colnames(df) <- c('rho', 'p.val')
  
  toplot <- cbind(Embeddings(obj, reduction = 'umap'), df)
  fwrite(toplot, glue::glue("{ct}_{s}_correlation_Spearman.tsv"), row.names = T, sep='\t')
  
  max.cor <- max(abs(toplot$rho))
  
  ggplot(data = toplot,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho')) +
    geom_point(
      size = dot.size) +
    theme_cowplot() +
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0('Featureplot_', ct, '_cor_', s, '.pdf'), useDingbats = F, width = 5, height = 4.5)
  
  #color dots only if p.value is < 0.05
  ggplot() +
    geom_point(data = subset(toplot, p.val > 0.05), 
               aes_string(x = 'UMAP_1', y = 'UMAP_2'), 
               size = 0.5, color = 'lightgrey') +
    geom_point(data = subset(toplot, p.val < 0.05),  
               aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'rho'),
               size = dot.size) +
    theme_cowplot() +
    ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = "Spearman's rho", low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c', limits = c(-max.cor, max.cor))
  ggsave(paste0( 'Featureplot_', ct, '_cor_', s, '_significant_only.pdf'), useDingbats = F, width = 5, height = 4.5)
  
}

###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object to score",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="foo",
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-m","--metadata"),
              type="character",
              default='/diskmnt/Projects/HTAN_analysis_2/BRCA/Analyses/Alla/DCIS_project/cell_typing/st_obj/',
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
meta.path <- opt$metadata
add_filename <- opt$extra

dir.create(out_path, showWarnings = F, recursive = TRUE)
setwd(out_path)

obj <- readRDS(input.path)
meta <- fread(meta.path, header = TRUE) %>%
  data.frame(row.names = 1)
obj <- AddMetaData(obj, meta)

cat('doing one way scoring')
message('do CP genes')
c2_gene_sets = msigdbr(species = "human", category="C2")

msig_df <- c2_gene_sets %>% dplyr::filter(grepl('SENESCE|senesce|Senesce|SASP', gs_name)) %>% 
  dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()
geneSets.msig <- split(msig_df$gene_symbol, msig_df$gs_name)

c5_gene_sets = msigdbr(species = "human", category="C5")

msig_df <- c5_gene_sets %>% dplyr::filter(grepl('SENESCE|senesce|Senesce|SASP', gs_name)) %>% 
  dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()
geneSets.msig2 <- split(msig_df$gene_symbol, msig_df$gs_name)

#add senmayo scores
senmayo <- fread('/diskmnt/Projects/SenNet_analysis/Main.analysis/gene_sets/SenMayo_geneset.txt')
senmayo  <- senmayo %>% mutate(Senmayo.sub = paste0('Senmayo_', Classification))

geneSets.senmayo <- list('SenMayo' = senmayo$Gene)
geneSets.senmayo <- c(geneSets.senmayo, split(senmayo$Gene, senmayo$Senmayo.sub))

#add sensig scores
sensig <- fread('/diskmnt/Projects/SenNet_analysis/Papers/SenSig/SkmFibGenes_SenSig_res0.05.ovl_with_cherryetal.csv')
sensig <- sensig %>% mutate(Direction = case_when(sign(logFC)==1 ~ 'SenSig_up',
                                                  TRUE ~ 'SenSig_down')) %>%
  select(HGNC.symbol,Direction ) %>%
  filter(!is.na(HGNC.symbol)) %>%
  filter((HGNC.symbol!=''))

geneSets.sensig <- split(sensig$HGNC.symbol, sensig$Direction)
geneSets.sensig

geneSets <- c(geneSets.msig, geneSets.msig2, geneSets.senmayo, geneSets.sensig)
names(geneSets) <- make.names(names(geneSets))
names(geneSets) <- str_replace_all(names(geneSets), '_', '-')

print(names(geneSets))

#Run_AUCell(obj, geneSets, name = add_filename)



cat('doing two way scoring')

signatures <- fread('/diskmnt/Projects/SenNet_analysis/Papers/Unmasking_transcriptional_heterogeneity/processed_sen_signatures_Segura_paper.txt', data.table = F)
signatures

sensig <- fread('/diskmnt/Projects/SenNet_analysis/Papers/SenSig/SkmFibGenes_SenSig_res0.05.ovl_with_cherryetal.csv')
sensig <- sensig %>% 
  select(HGNC.symbol,logFC ) %>%
  filter(!is.na(HGNC.symbol)) %>%
  filter((HGNC.symbol!='')) %>%
  mutate(ont='SenSig') %>%
  rename(Gene_names=HGNC.symbol,
         log2FC=logFC)
all.signatures <- rbind(signatures, sensig)


all.signatures$ont %>% unique %>% walk(function(s) {
  Run_correlation(obj, s, add_filename, dot.size = 0.1, cor.adjust.method = 'fdr')})
  
  
  



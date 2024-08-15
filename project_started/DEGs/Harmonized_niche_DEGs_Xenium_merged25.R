
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
set.seed(1234)
suppressMessages(library(data.table))
library(cowplot)


setwd('/diskmnt/Projects/SenNet_analysis/Main.analysis/old_vs_young/xenium_validation/cell_typing/label_transfer/anno_06182024')

merged.xen <- readRDS('/diskmnt/Projects/SenNet_analysis/Main.analysis/merged/merge_liver_Xenium_all_2/25_Merged_normalized_liver_Xenium_all.rds')

niche.labels <- fread('/diskmnt/Projects/SenNet_analysis/Main.analysis/kgallant/xenium_analysis/banksy_072024/banksy_20_harmonized_annot_louvain.tsv', header = FALSE, col.names = c('unknown', 'cell_id', 'niche')) %>%
  select(-unknown) %>%
  mutate(niche = paste('Niche', niche)) 

niche.labels

merged.xen <- AddMetaData(merged.xen, column_to_rownames(niche.labels, 'cell_id'))

merged.xen$niche[is.na(merged.xen$niche)]  <- 'irrelevant'


all.niches <- merged.xen$niche %>% unique()
all.niches <- all.niches[all.niches!='irrelevant']

deg.all <- all.niches %>% map_dfr(function(ni) {
  Idents(merged.xen) <- 'niche'
  print(ni)
  deg.mine <- FindMarkers(object = merged.xen, 
                          min.pct = 0, logfc.threshold = 0.1,latent.vars = 'Patient_ID',
                          `ident.1` = ni, `ident.2` = setdiff(all.niches, ni), 
                          test.use = 'MAST', assay = 'SCT')
  deg.mine <- deg.mine %>% mutate(gene= rownames(deg.mine),
                                  Niche = ni)
  fwrite(deg.mine, glue::glue('/diskmnt/Projects/SenNet_analysis/Main.analysis/old_vs_young/xenium_validation/BANKSY/harmonized_niche_DEGs/{ni}_MAST_DEGs_20cases_laten_var_Patient_ID.tsv'), sep='\t')
})

deg.all %>% fwrite(glue::glue('/diskmnt/Projects/SenNet_analysis/Main.analysis/old_vs_young/xenium_validation/BANKSY/harmonized_niche_DEGs/All_niches_MAST_DEGs_20cases_laten_var_Patient_ID.tsv'), sep='\t')

top3 <- deg.all %>% filter(avg_log2FC > 0.1) %>% arrange(Niche) %>% group_by(Niche) %>% top_n(n=3, wt=avg_log2FC)
top3

dotplot.color <- colorRampPalette(c('#eae2b7','#f3d180', '#fcbf49','#f77f00','#d62828','#6b2c39', '#003049'))(10) 

p <- DotPlot(merged.xen,
             features = unique(top3$gene),
             dot.min = 0,
             group.by = 'niche') + 
  scale_color_gradientn(colors = dotplot.color) +
  theme(axis.text.x = element_text(angle =90, hjust=1, vjust=0.5))+
  scale_size_area(limits=c(0,50), oob=scales::squish) +
  coord_fixed()
p
ggsave(glue::glue('/diskmnt/Projects/SenNet_analysis/Main.analysis/old_vs_young/xenium_validation/BANKSY/harmonized_niche_DEGs/Dotplot_harm.niche.DEGs_mine_merged_object_by_niches.pdf'), plot=p,
       width= 13, height = 6, useDingbats = FALSE)



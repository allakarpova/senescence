library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)
library(paletteer)
library(viridis)
library(tidyverse)
library(cowplot)

metas <- c('~/lab_Ding/work/Sennet/cell_typing/manual/SN001H1-Md1A3Y1K2N1/SN001H1-Md1A3Y1K2N1_processed_multiomic_cellTyped.meta.data',
           '~/lab_Ding/work/Sennet/cell_typing/manual/SN004H1-Md1A4K1G1/SN004H1-Md1A4K1G1_processed_multiomic_cellTyped.meta.data')

data.type <- c('multiome', "5' sc")

toplot <- map2(metas, data.type, function(x, y) {
  tb <- fread(x)
  tb$data.type <- y
  return(tb[, c('orig.ident', 'V1', 'cell_type', 'data.type')])
}) %>% rbindlist() %>%
  mutate(cell_type_broad = case_when(grepl('DC', cell_type) ~ 'DC',
                                     grepl('LSECs', cell_type) ~ 'LSECs',
                                     TRUE ~ cell_type)) %>%
  group_by(orig.ident, cell_type_broad, data.type) %>% tally(name='N')


ggplot(toplot, aes(x = data.type, y = N, fill = cell_type_broad)) +
  geom_bar(stat = 'identity') +
  theme_cowplot() +
  scale_fill_paletteer_d("ggsci::category20c_d3")

ggplot(toplot, aes(x = data.type, y = N, fill = cell_type_broad)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_cowplot() +
  scale_fill_paletteer_d("ggsci::category20c_d3")



#################
meta <- fread('~/lab_Ding/work/Sennet/analysis/merged/merge_liver_combo_2/upd.cell.type/5_Merged_normalized_liver_combo.metadata_06052022.tsv')

toplot <- meta %>% filter(cell_type_merged != 'Doublet') %>% mutate(Age.group = ifelse(Age >= 50, 'Old', 'Young'))

ggplot(toplot, aes(x = cell_type_merged,  fill = cell_type_merged)) +
  geom_bar(stat = 'count', position = 'dodge') +
  theme_cowplot() +
  facet_wrap(~ Age.group + Sample , scales = 'free', ncol = 5) +
  scale_fill_paletteer_d("ggsci::category20c_d3") +
  theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust=0.5))

ggsave('~/lab_Ding/work/Sennet/analysis/merged/merge_liver_combo_2/upd.cell.type/Cell_type_proportion.summary.pdf',
       width = 19, height = 5)



paths <- c('~/lab_Ding/work/Sennet/analysis/merged/recluster_cell_types_from_merge_liver_combo_2/Pathways_scores_differential_activity_Hepatocytes.tsv'
          )
cell.types <- c('Hepatocytes')
ucell <- fread('~/lab_Ding/work/Sennet/analysis/merged/recluster_cell_types_from_merge_liver_combo_2/Ucell_CPDB_pathways_Hepatocytes.tsv') %>% 
  group_by(Sample, cell_type_merged) %>% summarise_all( mean)

path.both <- map2(paths, cell.types, function(x, y) {
  tb <- fread(x)
  tb$cell.type <- y
  return(tb)
}) %>% rbindlist()

common.deg <- 
  path.both %>% 
  mutate(FC_sign = sign(avg_diff)) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(V1, FC_sign, cell.type) %>%
  tally %>% arrange(-n) %>% 
  filter(n >= 2)


toplot <- path.both %>% 
  mutate(FC_sign = sign(avg_diff)) %>%
  filter(p_val_adj < 0.05 & abs(avg_diff) > 0.1) %>%
  group_by( FC_sign, cell.type) %>% top_n(10, wt = -p_val_adj) %>%
  mutate(Path = str_split_fixed(V1, '-', 2)[,1],
         Path = str_split_fixed(Path, '[.]', 2)[,2],
         Path = gsub(pattern = '[.]', replacement = ' ', x = Path),
         Path = tolower(Path)) %>% 
  arrange(avg_diff) %>%
  mutate(Path = factor(Path, levels = Path))



ggplot(toplot, aes(x = Path, y= avg_diff,  fill = as.character(FC_sign))) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_cowplot() +
  #facet_wrap(~ cell.type , scales = 'free', ncol = 2) +
  #scale_fill_manual(c('blue', 'red')) +
  theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust=0.5))

ggsave('~/lab_Ding/work/Sennet/analysis/merged/recluster_cell_types_from_merge_liver_combo_2/Paths_avg_diff_Hepatocytes.pdf',
       width = 10, height = 10)



ucell %>% filter(cell_type_merged == 'Hepatocytes') %>% 
  select(Sample, str_replace(toplot$V1, '-', '_')) %>%
  melt() %>%
  mutate(Age.group = case_when(grepl('004|007', Sample) ~ 'Young', TRUE ~ 'Old')) %>%
  ggplot( aes(x = Sample, y= value,  fill = Age.group)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_cowplot() +
  facet_wrap(~ variable , scales = 'free', ncol = 10) +
  #scale_fill_manual(c('blue', 'red')) +
  theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust=0.5))



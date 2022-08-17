library(data.table)
library(dplyr)
library(readxl)

setwd('~/lab_Ding/work/single_cell/senescence/Papers/Unmasking_transcriptional_heterogeneity/')
fibro.onco <- read_excel('/1-s2.0-S0960982217308928-mmc2.xlsx', sheet = 2, skip = 2)
fibro.replic <- read_excel('1-s2.0-S0960982217308928-mmc2.xlsx', sheet = 3, skip = 2)
fibro.radia <- read_excel('1-s2.0-S0960982217308928-mmc2.xlsx', sheet = 4, skip =2)
fibro.meta <- read_excel('1-s2.0-S0960982217308928-mmc2.xlsx', sheet = 5, skip =2)

colnames(fibro.onco)
fibro.onco <- fibro.onco %>% 
  select (Gene_Names, log2FoldChange) %>% 
  rename(log2FoldChange = 'log2FC', Gene_Names = 'Gene_names') %>% 
  mutate(ont = 'Fibroblast_oncogene_ind_sen')

colnames(fibro.replic)
fibro.replic <- fibro.replic %>% 
  dplyr::select (Gene_names, `NB-GLM_FC`) %>% 
  rename(`NB-GLM_FC` = 'log2FC') %>% 
  mutate(ont = 'Fibroblast_replicative_sen')

colnames(fibro.radia)
fibro.radia <- fibro.radia %>% 
  select (Gene_names, log2FoldChange) %>% 
  rename(log2FoldChange = 'log2FC') %>%
  mutate(ont = 'Fibroblast_radiation_ind_sen')

colnames(fibro.meta)
fibro.meta <- fibro.meta %>% 
  select (Gene_names, `Meta_NB-GLM_FC`) %>% 
  rename(`Meta_NB-GLM_FC` = 'log2FC') %>%
  mutate(ont = 'Fibroblast_meta_sen')


fibro.sign <- do.call('rbind', list(fibro.onco, fibro.replic, fibro.radia, fibro.meta))


keratin <- read_excel('1-s2.0-S0960982217308928-mmc3.xlsx', sheet = 2, skip = 2)
melano <- read_excel('1-s2.0-S0960982217308928-mmc3.xlsx', sheet = 3, skip = 2)
astro <- read_excel('1-s2.0-S0960982217308928-mmc3.xlsx', sheet = 4, skip =2)
fibro <- read_excel('1-s2.0-S0960982217308928-mmc3.xlsx', sheet = 5, skip =2)
univers <- read_excel('1-s2.0-S0960982217308928-mmc3.xlsx', sheet = 6, skip =2)


colnames(keratin)
keratin <- keratin %>% 
  select (Gene_Names, log2FoldChange) %>% 
  rename(log2FoldChange = 'log2FC', Gene_Names = 'Gene_names') %>% 
  mutate(ont = 'Keratinocytes_sen')

colnames(melano)
melano <- melano %>% 
  select (Gene_Names, log2FoldChange) %>% 
  rename(log2FoldChange = 'log2FC', Gene_Names = 'Gene_names') %>% 
  mutate(ont = 'Melanocytes_sen')

colnames(astro)
astro <- astro %>% 
  select (Gene_Names, log2FoldChange) %>% 
  rename(log2FoldChange = 'log2FC', Gene_Names = 'Gene_names') %>% 
  mutate(ont = 'Actrocytes_sen')

colnames(fibro)
fibro <- fibro %>% 
  select (Gene_Names, log2FC_NB_GLM) %>% 
  rename(log2FC_NB_GLM = 'log2FC', Gene_Names = 'Gene_names') %>%
  mutate(ont = 'Fibroblasts_sen')

colnames(univers)
univers <- univers %>% 
  mutate(log2FC = rowMeans(select(., ends_with("FC")), na.rm = TRUE) ) %>%
  select (Gene_names, log2FC) %>% 
  mutate(ont = 'Universal_sen')

cell.type.sign <- do.call('rbind', list(keratin, melano, astro, fibro, univers))

fwrite(rbind(fibro.sign, cell.type.sign), 'processed_sen_signatures_Segura_paper.txt', sep = '\t')  
  
  
  
  
  
  

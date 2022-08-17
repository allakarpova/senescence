library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)
library(paletteer)
library(viridis)
library(cowplot)


########################
emb <- fread('~/lab_Ding/work/single_cell/senescence/gbm/Integrated_object_embeddings_gbm_28_sample.txt', data.table = F)
scores <- fread('~/lab_Ding/work/single_cell/senescence/gbm/Senescense_scores_gbm_28_sample.txt', data.table = F)
meta <- fread('~/lab_Ding/work/single_cell/gbm/snRNA_only_study/metadata/snRNA_metadata.v2.0.tsv.gz', data.table = F)

master <- merge(emb, meta, by.x = 'V1', by.y = 'cell_id')
master <- merge(master, scores, by.x = 'V1', by.y = 'V1')
master$sample_type <- ifelse(grepl('PT-', master$case_id), 'normal','tumor')
head(master)

toplot <- master
toplot <- toplot %>% arrange(SASP)
toplot$SASP[toplot$SASP < 0] <- 0
ggplot(data = toplot, aes(x = UMAP_1, y = UMAP_2, color = SASP)) +
  geom_point(size = 0.02) +
  theme_cowplot() +
  ggtitle ('SASP') + 
  facet_wrap(~sample_type, ncol = 2) +
  scale_color_gradient(low = 'grey', high = '#e41a1c')

ggsave('~/lab_Ding/work/single_cell/senescence/gbm/Featureplot_all_SASP_modules_gbm_28_sample.pdf', width = 10, height = 5, useDingbats = F)

toplot <- master
toplot <- toplot %>% arrange(Generic)
toplot$Generic[toplot$Generic < 0] <- 0
ggplot(data = toplot, aes(x = UMAP_1, y = UMAP_2, color = Generic)) +
  geom_point(size = 0.02) +
  theme_cowplot() +
  ggtitle ('Generic') + 
  facet_wrap(~sample_type, ncol = 2) +
  scale_color_gradient(low = 'grey', high = '#e41a1c')

ggsave('~/lab_Ding/work/single_cell/senescence/gbm/Featureplot_all_Generic_modules_gbm_28_sample.pdf', width = 10, height = 5, useDingbats = F)

master.normal <- master %>% filter (sample_type == 'normal')
which(master.normal[60:100]==max(master.normal[60:100]))

sapply(master.normal[60:100], function (x) which(x==max(master.normal[60:100])))
colnames(master.normal)
toplot <- master
toplot <- toplot %>% arrange(`Wikipathways-DNA.DAMAGE.RESPONSE..ONLY.ATM.DEPENDENT.`)
toplot$`Wikipathways-DNA.DAMAGE.RESPONSE..ONLY.ATM.DEPENDENT.`[toplot$`Wikipathways-DNA.DAMAGE.RESPONSE..ONLY.ATM.DEPENDENT.` < 0] <- 0

ggplot(data = toplot, aes(x = UMAP_1, y = UMAP_2, color = `Wikipathways-DNA.DAMAGE.RESPONSE..ONLY.ATM.DEPENDENT.`)) +
  geom_point(size = 0.02) +
  theme_cowplot() +
  ggtitle ('Wikipathways-DNA.DAMAGE.RESPONSE..ONLY.ATM.DEPENDENT.') + 
  facet_wrap(~mgmt_methyl_stp27_prediction, ncol = 3) +
  scale_color_gradient(low = 'grey', high = '#e41a1c', name = 'score')

ggsave('~/lab_Ding/work/single_cell/senescence/gbm/Featureplot_all_KReactome-Wikipathways-DNA.DAMAGE.RESPONSE..ONLY.ATM.DEPENDENT._by_MGMT_status_gbm_28_sample.pdf', width = 15, height = 5, useDingbats = F)


colnames(master.normal)
toplot <- master
toplot <- toplot %>% arrange(`Reactome-ONCOGENE.INDUCED.SENESCENCE`)
toplot$`Reactome-ONCOGENE.INDUCED.SENESCENCE`[toplot$`Reactome-ONCOGENE.INDUCED.SENESCENCE` < 0] <- 0

ggplot(data = toplot, aes(x = UMAP_1, y = UMAP_2, color = `Reactome-ONCOGENE.INDUCED.SENESCENCE`)) +
  geom_point(size = 0.02) +
  theme_cowplot() +
  ggtitle ('Reactome-ONCOGENE.INDUCED.SENESCENCE') + 
  facet_wrap(~sample_type, ncol = 2) +
  scale_color_gradient(low = 'grey', high = '#e41a1c', name = 'score')

ggsave('~/lab_Ding/work/single_cell/senescence/gbm/Featureplot_all_Reactome-ONCOGENE.INDUCED.SENESCENCE_modules_gbm_28_sample.pdf', width = 10, height = 5, useDingbats = F)

######
ggplot(data = toplot, aes(x = UMAP_1, y = UMAP_2, color = case_id)) +
  geom_point(size = 0.1) +
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=3)))
  
ggsave('~/lab_Ding/work/single_cell/senescence/gbm/Dimplot_case_id_gbm_28_sample.pdf', width = 8 , height = 5, useDingbats = F)


ggplot(data = master, aes(x = UMAP_1, y = UMAP_2, color = cell_type_v20210126)) +
  geom_point(size = 0.1) +
  theme_cowplot() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  #ggtitle ('Generic') + 
  scale_color_brewer(palette = 'Paired')

ggsave('~/lab_Ding/work/single_cell/senescence/gbm/Dimplot_celltype_gbm_28_sample.pdf', width = 7.5, height = 5, useDingbats = F)

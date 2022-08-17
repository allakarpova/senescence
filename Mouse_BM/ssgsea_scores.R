library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(stringr)
library(cowplot)
library(ggplot2)
library(GSVA)

panc <- readRDS('~/lab_Ding/work/single_cell/senescence/objects/Cell_type.object.res.7.Fibro_Mouse_BM_second.RDS')
cell.type <- 'Fibro'
source('~/R_working_dir/scripts/senescence/markers.R')
DefaultAssay(panc) <- 'RNA'

gene.sets.list <- list (Sen.core = sen.core.mouse, Sen.effector = sen.effector.mouse, sasp = sasp.mouse)
expr.mat <- FetchData(panc, vars = rownames(panc), slot = 'data') %>% t()# %>% scale()
###### 
ss.results.prot <- gsva(expr = as.matrix(expr.mat), 
                        gset.idx.list = total.sen,
                        method=c( "ssgsea"),
                        kcdf=c("Gaussian"),
                        abs.ranking=FALSE,
                        min.sz=1,
                        max.sz=Inf,
                        parallel.sz=0,
                        parallel.type="SOCK",
                        mx.diff=TRUE,
                        tau=0.25,
                        ssgsea.norm=T,
                        verbose=TRUE)

scores <- cbind(Embeddings(panc, reduction = 'umap'),seurat_clusters = panc$seurat_clusters,sample = panc$orig.ident, data.frame(t(ss.results.prot)))
toplot.melt <- reshape2::melt(scores, id.vars = c('seurat_clusters','sample', 'UMAP_1', 'UMAP_2'))

wd <- '~/lab_Ding/work/single_cell/senescence/mouse_BM/senescence_ssgsea_scores/'
dir.create(wd, recursive = T)
setwd(wd)

ggplot(data = toplot.melt, aes(x = value, fill = sample)) +
  geom_density(alpha = 0.5) +
  facet_grid(~variable, scales = 'free') +
  theme_classic() +
  scale_fill_manual(values =c( '#fc8d62','#8da0cb', '#66c2a5'))
#scale_x_continuous( breaks = c(-0.3, 0, 0.3, 0.6), )
ggsave(paste0('Density_plot_',cell.type,'_senes.sets.pdf'), useDingbats = F, width = 10, height =4)

for (s in names(gene.sets.list)) {
  scores <- scores[order(scores[s]),]
  ggplot(data = scores,  
         aes_string(x = 'UMAP_1', y = 'UMAP_2', color = s)) +
    geom_point( size = 1) +
    facet_grid(~sample, scales = 'free') +
    theme_cowplot() +
    #ggtitle (paste0('Correlation with\n', s)) +
    scale_color_gradient2(name = paste(s, "\nssggsea score"), low = '#1f78b4', mid = 'lightgrey', high = '#e41a1c')
    #scale_color_gradient(name = paste(s, "\nssggsea score"), low = 'lightgrey', high = '#e41a1c')
  ggsave(paste0('Featureplot_ssgsea_', s, '.pdf'), useDingbats = F, width = 10, height = 4.5)
  
}







################# testing
panc <- AddModuleScore(
  object = panc,
  assay = 'RNA',
  features = total.sen,
  ctrl = 10,
  name = 'senescence'
)
cat('done\n')

n.list <- length(total.sen)
colnames(panc@meta.data)[(ncol(panc@meta.data) - n.list + 1) : ncol(panc@meta.data)] <- names(total.sen)
test <- data.frame(cbind(panc@meta.data[names(total.sen)], data.frame(t(ss.results.prot))), check.rows = F, check.names = F)
test$sample = panc$orig.ident
str(test)
colnames(test)[1:3] <- paste0(colnames(test)[1:3], '.module')
head(test)

ggplot(data = test, aes(x = Sen.core.module, y = Sen.core, color = sample)) +
  geom_point(size = 0.2) +
  #facet_grid(~variable, scales = 'free') +
  theme_classic() +
  ylab ('ssgsea score') + xlab ('seurat module scores') + 
  scale_color_manual(values =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Scatterplot_Sencore.ssgse_vs_module.pdf', width = 7, height = 5.5 , useDingbats = F)

ggplot(data = test, aes(x = Sen.effector.module, y = Sen.effector, color = sample)) +
  geom_point(size = 0.2) +
  #facet_grid(~variable, scales = 'free') +
  theme_classic() +
  ylab ('ssgsea score') + xlab ('seurat module scores') + 
  scale_color_manual(values =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Scatterplot_Seneffector.ssgsea_vs_module.pdf', width = 7, height = 5.5 , useDingbats = F)

ggplot(data = test, aes(x = sasp.module, y = sasp, color = sample)) +
  geom_point(size = 0.2) +
  #facet_grid(~variable, scales = 'free') +
  theme_classic() +
  ylab ('ssgsea score') + xlab ('seurat module scores') + 
  scale_color_manual(values =c( '#fc8d62','#8da0cb', '#66c2a5'))
ggsave('Scatterplot_SASP.ssgsea_vs_module.pdf', width = 7, height = 5.5 , useDingbats = F)




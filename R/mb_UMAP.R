# UMAP


#--- FUNCTION AS IS NOW IS NOT FLEXIBLE, IT USES COLUMN NAMES FROM METADATA, PLEASE CONSIDER THIS AND EXCLUDE, OR USE color=group like in the other dimensionality reduction plotters!
library(embed)
library(tidymodels)

pan2umap <- function(pan, metadata, color="group", title_suffix="", color_scale=NULL,discrete=FALSE){
  meta.sub = metadata %>% select(genome,Genome_size,group)
  pan = pan[metadata$genome,]
  pan = pan[,which(colSums(pan)!=0)]
  pan = pan[,which(colSums(pan)!=1)]
  pan <-pan[,which(apply(pan,2,sd)>0)]

  df.pan = cbind(meta.sub,pan[meta.sub$genome,])

  umap_rec <- recipe(~.,data=df.pan) %>%
    update_role(genome, Genome_size, group, new_role = 'id' ) %>% #descide what it is...
    step_normalize(all_predictors()) %>% #  center and scale
    step_umap(all_predictors())
  #umap_rec

  #prepare recipe
  umap_prep <- prep(umap_rec)
  umap_prep

  if(discrete==FALSE){
    p <- juice(umap_prep) %>%
      ggplot(aes(umap_1,umap_2,fill=group))
  }else{
    p <- juice(umap_prep) %>%
      ggplot(aes(umap_1,umap_2,fill=as.factor(group)))
  }

  #(why only 5 dimensions?)
  p <- p +
    geom_point(shape=21,size=2)+
    theme_classic() +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    ggtitle(label=paste0(title_suffix))+
    geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
    geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
    fdb_style() + theme(legend.position = "none")

  if(!is.null(color_scale)){
    p <- p + scale_fill_manual(values=color_scale, na.value = 'white')
  }

  return(p)
}

#=============


pan2pca <- function(pan, metadata, color="group", title_suffix="", color_scale=NULL,discrete=FALSE){
  meta.sub = metadata %>% select(genome,Genome_size,group)
  pan = pan[metadata$genome,]
  pan = pan[,which(colSums(pan)!=0)]
  pan = pan[,which(colSums(pan)!=1)]
  pan <-pan[,which(apply(pan,2,sd)>0)]

  df.pan = cbind(meta.sub,pan[meta.sub$genome,])

  pca_rec <- recipe(~.,data=df.pan) %>%
    update_role(genome, Genome_size, group, new_role = 'id' ) %>% #descide what it is...
    step_normalize(all_predictors()) %>% #  center and scale
    step_pca(all_predictors())
  #umap_rec

  #prepare recipe
  pca_prep <- prep(pca_rec)
  pca_prep

  if(discrete==FALSE){
    p <- juice(pca_prep) %>%
      ggplot(aes(PC1,PC2,fill=group))
  }else{
    p <- juice(pca_prep) %>%
      ggplot(aes(PC1,PC2,fill=as.factor(group)))
  }

  #(why only 5 dimensions?)
  p <- p +
    geom_point(shape=21,size=2)+
    theme_classic() +
    xlab('PC1') +
    ylab('PC2') +
    ggtitle(label=paste0(title_suffix))+
    geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
    geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
    fdb_style() + theme(legend.position = "none")

  if(!is.null(color_scale)){
    p <- p + scale_fill_manual(values=color_scale, na.value = 'white')
  }

  return(p)
}


#=================================================================================#
# Run tSNE on all pan tables

selcol <- 'sc0.85'
colcol <- 'cc0.85'
metdad <- genome.tbl %>%
  filter(genome%in% sorted.cliques$genome) %>%
  left_join(sorted.cliques,by='genome') %>%
  filter(grepl('Marinobacter|Tamil',genome)) %>%
  left_join(out.colors, by='genome') %>%
  filter(!is.na(!!as.name(selcol))) %>%
  filter(quality_class != 'low') %>%
  data.frame()
metdad$group <- paste0('cl',metdad[,colcol])

rownames(metdad) <- metdad$genome
metdad[,selcol] <- as.factor(metdad[,selcol])
metdad[,colcol] <- as.factor(metdad[,colcol])


colors <- metdad %>% select(paste0('col.',colcol),colcol, group) %>% unique() %>% data.frame()
clique.colors <- as.character(colors[,1])
names(clique.colors) <- as.character(colors[,'group'])


plotslist.t <- list()
for(pan.name in names(pan.list)){
  message(paste0('analysing ',pan.name))
  pan <- pan.list[[pan.name]]
  pan.p.umap <- pan2umap(pan[rownames(metdad),],
                         metadata=metdad,
                         color="group",
                         title_suffix=pan.name,
                         color_scale= clique.colors,
                         discrete=TRUE)
  plotslist.t[[pan.name]] <- pan.p.umap
}

p.all.pan.umap <-grid.arrange(grobs = plotslist.t, ncol = 6)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_UMAP.pdf',plot=p.all.pan.umap, width = 12, height = 4,unit='in')



plotslist.t <- list()
for(pan.name in names(pan.list)){
  message(paste0('analysing ',pan.name))
  pan <- pan.list[[pan.name]]
  pan.p.pca <- pan2pca(pan[rownames(metdad),],
                         metadata=metdad,
                         color="group",
                         title_suffix=pan.name,
                         color_scale= clique.colors,
                         discrete=TRUE)
  plotslist.t[[pan.name]] <- pan.p.pca
}

p.all.pan.pca <-grid.arrange(grobs = plotslist.t, ncol = 6)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_PCA.pdf',plot=p.all.pan.pca, width = 12, height = 4,unit='in')


plotlist <- list('PCoA'=p.all.pan,
                 'PCA'=p.all.pan.pca,
                 'tSNE'=p.all.pan.tsne,
                 'UMAP'=p.all.pan.umap)
p.all.pan.dim_red <-grid.arrange(grobs = plotlist, nrow = 4)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_dim_red.pdf',plot=p.all.pan.dim_red, width = 12, height = 12,unit='in')

#----


p.umap.ANIb = pan2umap(pan = ANIb[rownames(metadata),rownames(metadata)],metadata=metadata,color="group",
                       title_suffix=pan.name,
                       color_scale= clique.colors,
                       discrete=TRUE)
p.umap.ANIb

p.umap.AAI = pan2umap(pan = AAI[rownames(metadata),rownames(metadata)],metadata=metadata,title_suffix="AAI")
p.umap.AAI


p.umap.aa_usage = pan2umap(pan = aa_usage[rownames(metadata),],metadata=metadata,title_suffix="aa_usage")
p.umap.aa_usage
p.umap.codon_usage = pan2umap(pan = codon_usage[rownames(metadata),],metadata=metadata,title_suffix="codon_usage")
p.umap.codon_usage
p.umap.kmer_usage = pan2umap(pan = kmer_usage[rownames(metadata),],metadata=metadata,title_suffix="kmer_usage")
p.umap.kmer_usage


p.umap.OFpan = pan2umap(pan = OG.pan,metadata=metadata,title_suffix="OF")
p.umap.OFpan

p.umap.COGpan = pan2umap(pan = COG.pan,metadata=metadata,title_suffix="COG")
p.umap.COGpan

p.umap.kfam.pan = pan2umap(pan = kfm.pan ,metadata=metadata,title_suffix="KOfam")
p.umap.kfam.pan

p.umap.cazy.pan = pan2umap(pan = cazy.pan ,metadata=metadata,title_suffix="CAZY")
p.umap.cazy.pan

p.umap.tcdb.pan = pan2umap(pan = tcdb.pan ,metadata=metadata,title_suffix="TCDB")
p.umap.tcdb.pan



p.umap.combined <- gridExtra::grid.arrange(p.umap.OFpan, p.umap.COGpan, p.umap.kfam.pan, p.umap.cazy.pan,p.umap.tcdb.pan ,nrow=1)
p.umap.combined


p.combined.dimred <- gridExtra::grid.arrange(p.pcoa.combined, p.tsne.combined,p.umap.combined, ncol=1)
p.combined.dimred



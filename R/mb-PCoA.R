# PAN GENOME BASED DIMENSIONALITY REDUCTION
#--------------------------------------------------------------#

#### FUNCTIONS ####

# spiderized PCA plot
#---------------------------------------------------------------#
#' plots PCA as spider using specified group
#'
#' @param pan pan/abundance table
#' @param metadata metadata table with rownames
#' @param color column name in metadata table used to color the points
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' plotSpiderPCA(m.pan, metadata ,color='phylogroup') +
#' scale_fill_manual(values=phylogroup.colors) +
#' scale_color_manual(values=phylogroup.colors)

plotSpiderPCA <- function(pan, metadata, color='group'){

  compPCA.scrs = c()
  compPCA.cent = c()
  compPCA.segs = c()

  pca_data=prcomp(pan[rownames(metadata),], scale = F)

  pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
  scrs <- vegan::scores(pca_data, display = 'sites')
  scrs <- cbind(as.data.frame(scrs), group = metadata[rownames(scrs),color])
  cent <- aggregate(cbind(PC1, PC2) ~ group, data = scrs, FUN = mean)
  segs <- merge(scrs, setNames(cent, c('group','oPC1','oPC2')),
                by = 'group', sort = FALSE)

  compPCA.scrs = bind_rows(compPCA.scrs, scrs)
  compPCA.cent = bind_rows(compPCA.cent, cent)
  compPCA.segs = bind_rows(compPCA.segs, segs)

  compPCA.scrs =data.frame(compPCA.scrs)
  compPCA.cent =data.frame(compPCA.cent)
  compPCA.scrs =data.frame(compPCA.scrs)
  pcaPercExpl =data.frame(pcaPercExpl)

  pca.plots = ggplot2::ggplot(compPCA.scrs, ggplot2::aes(x = PC1, y = PC2, fill=group)) +
    ggplot2::geom_segment(data = compPCA.segs,mapping = aes(xend = oPC1, yend = oPC2, colour = group),alpha=.6) + # spiders
    ggplot2::geom_point(ggplot2::aes(color=group),shape=21,alpha=0.6,size=2) +                                              # sample scores
    ggplot2::geom_point(data = compPCA.cent, size = 3.5,shape=21) +                         # centroids
    ggplot2::labs(x=paste0("PC1"), y=paste0("PC2")) +
    ggplot2::theme_classic()+
    fdb_style()

  return(pca.plots)
}

#================================================================================#

plotSpiderPCoA <- function(pan.dist, metadata, color='group'){

  compPCA.scrs = c()
  compPCA.cent = c()
  compPCA.segs = c()

  pca_data=cmdscale(pan.dist, eig = TRUE)

  #pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
  scrs <- vegan::scores(pca_data, display = 'sites')
  scrs <- cbind(as.data.frame(scrs), group = metadata[rownames(scrs),color])
  cent <- aggregate(cbind(Dim1, Dim2) ~ group, data = scrs, FUN = mean)
  segs <- merge(scrs, setNames(cent, c('group','oPC1','oPC2')),
                by = 'group', sort = FALSE)

  compPCA.scrs = bind_rows(compPCA.scrs, scrs)
  compPCA.cent = bind_rows(compPCA.cent, cent)
  compPCA.segs = bind_rows(compPCA.segs, segs)

  compPCA.scrs =data.frame(compPCA.scrs)
  compPCA.cent =data.frame(compPCA.cent)
  compPCA.scrs =data.frame(compPCA.scrs)
  pcaPercExpl =data.frame(pcaPercExpl)

  pca.plots = ggplot2::ggplot(compPCA.scrs, ggplot2::aes(x = PC1, y = PC2, fill=group)) +
    ggplot2::geom_segment(data = compPCA.segs,mapping = aes(xend = oPC1, yend = oPC2, colour = group),alpha=.6) + # spiders
    ggplot2::geom_point(ggplot2::aes(color=group),shape=21,alpha=0.6,size=2) +                                              # sample scores
    ggplot2::geom_point(data = compPCA.cent, size = 3.5,shape=21) +                         # centroids
    ggplot2::labs(x=paste0("PCoA1"), y=paste0("PCoA2")) +
    ggplot2::theme_classic()+
    fdb_style()

  return(pca.plots)
}




plotSpiderPCA(distances.tibble$distance[[1]], distances.tibble$metadata[[1]] ,color='phylogroup') +
  scale_color_manual(values=phylogroup.colors)+scale_fill_manual(values=phylogroup.colors)


i=4

betadisper.obj <- betadisper(as.dist(distances.tibble$distance[[i]]),
           group=distances.tibble$metadata[[i]]$group)

anova(betadisper.obj)

#anova(betadisper.obj)

TukeyHSD(betadisper.obj, which = "group", ordered = FALSE,
         conf.level = 0.95)



#betadisper.obj <- betadisper(as.dist(1-ANIb[rownames(distances.tibble$metadata[[i]]),rownames(distances.tibble$metadata[[i]])]),
#                             group=distances.tibble$metadata[[i]]$group)

plotSpiderPCA(1-ANIb[rownames(distances.tibble$metadata[[1]]),rownames(distances.tibble$metadata[[1]])],
              distances.tibble$metadata[[1]] ,color='phylogroup') +
  scale_color_manual(values=phylogroup.colors)+scale_fill_manual(values=phylogroup.colors)




# pan2dist
#=================================================================================#
#' Title
#'
#' @param pan
#' @param dist.method
#'
#' @return
#' @export
#'
#' @examples
pan2dist <- function(pan, dist.method="Manhattan"){
  if(dist.method=="Manhattan"){
    print("Calculating pairwise Manhattan distances")
    distPan <- micropan::distManhattan(pan)
  }
  if(dist.method=="Jaccard"){
    print("Calculating pairwise Jaccard distances")
    distPan <- micropan::distJaccard(pan)
  }
  distPan= as.matrix(distPan)
  distPan = distPan[ order(row.names(distPan)), ]
  distPan = distPan[ , order(colnames(distPan))]
  return(distPan)
}


#EXAMPLE USAGE
# metadata <- genome.tbl %>%
#   filter(genome %in% base.tree$tip.label) %>%
#   left_join(sorted.ANI.cliques,by='genome') %>%
#   data.frame()
# rownames(metadata) <- metadata$genome
#
# #pfam
# pfam.dist <- pan2dist(pfam.pan[sel.genome,], dist.method="Manhattan")
#
# pan2pcoa(dist.pan = pfam.dist,
#          metadata=metadata,
#          color="phylogroup",
#          title_suffix="pfam") + scale_fill_manual(values=phylogroup.colors)
#
# #OG
# og.dist <- pan2dist(OG.pan[sel.genome,], dist.method="Manhattan")
#
# pan2pcoa(dist.pan = og.dist,
#          metadata=metadata,
#          color="phylogroup",
#          title_suffix="OG") + scale_fill_manual(values=phylogroup.colors)
#




# bootsrtap dense clusters

#' Title
#'
#' @param metadata
#' @param sample.col
#' @param n
#'
#' @return
#' @export
#'
#' @examples
bootstrap_group <- function(metadata, sample.col='sc0.95',n=25){
  grps <- unique(metadata[,sample.col])
  grps <- grps[!is.na(grps)]

  bootstrapped <- list()
  for(it in 1:n){
    bootstrap.tips <- NULL
    for(i in grps){
      stps <- metadata[metadata[,sample.col]==i,'genome']
      stps <- stps[!is.na(stps)]
      stps <- stps %>% sample(1)
      bootstrap.tips <- c(bootstrap.tips, stps)
    }
    bootstrap.tips <- c(bootstrap.tips,
                        metadata %>% filter(is.na(sc0.95))%>%select(genome)%>%pull())
    bootstrapped[[it]] <- bootstrap.tips
  }
  return(bootstrapped)
}

bootstrap.tips <- bootstrap_group(sorted.cliques, sample.col='sc0.95',n=25)
#
# cazy.dist <- pan2dist(OG.pan[bootstrap.tips,], dist.method="Manhattan")
#
# cazy.dist <- pan2dist(OG.pan[bootstrap.tips,], dist.method="Manhattan")
# pan2pcoa(dist.pan = cazy.dist,
#          metadata=metadata[bootstrap.tips,],
#          color="phylogroup",
#          title_suffix="HGM_FDR<0.05",color_scale = phylogroup.colors)
#
#
# pan2tsne(OG.pan[bootstrap.tips,cOutSig],
#          perplexity=10, theta=0.5,
#          metadata=metadata[bootstrap.tips,],
#          color="phylogroup",
#          title_suffix="TCDB",
#          color_scale= phylogroup.colors,
#          discrete=TRUE)



# pan2pcoa
#=================================================================================#
#' Function that plots a PCoA/MDS plot based on a pan table
#'
#' @param pan input pan table (genome x features)
#' @param dist.method specifies the distance method, Manhattan or Jaccard
#' @param metadata metadata
#' @param color specify the column name in metadata used for coloring the points, can be discrete or continuous
#' @param title_suffix specifies what to put in the title of the plot
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' cazy.p.pcoa <- pan2pcoa(pan = cazy.pan,
#'     dist.method="Manhattan",
#'     metadata=metadata,
#'     color="group",
#'     title_suffix="CAZY")
#' cazy.p.pcoa
#'
pan2pcoa <- function(dist.pan, metadata=metadata, color="group",title_suffix="", color_scale=NULL,discrete=FALSE,ellipse=FALSE){

  if(is.factor(metadata[,color])){
    colorGroups = as.character(metadata[,color])
  }else{
    colorGroups = metadata[,color]
  }
  names(colorGroups) = rownames(metadata)

  print("using cmdscale to calculate MDS/PCoA")
  MDS = cmdscale(dist.pan, eig = TRUE)
  MDSdata = data.frame(MDS$points)
  MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
  MDSdata$Name = rownames(MDSdata)
  MDSdata$Group = colorGroups[MDSdata$Name]
  colnames(MDSdata) = c("x","y","Name","Group")

  if(discrete==FALSE){
    p = ggplot2::ggplot(MDSdata,aes(x=x,y=y,fill=Group,label=Name))
  }else{
    p = ggplot2::ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(Group),label=Name))
  }
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::xlab(paste("PCoA1 (", MDS.var.perc[1], "%)",sep="")) +
    ggplot2::ylab(paste("PCoA2 (", MDS.var.perc[2], "%)",sep="")) +
    ggplot2::ggtitle(label=title_suffix)+
    ggplot2::geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_point(shape = 21, size = 2) +
    fdb_style() +
    ggplot2::theme(legend.position = "none")

  if(ellipse==TRUE){
    p <- p + ggplot2::stat_ellipse(aes(fill=as.factor(Group)))
  }
  if(!is.null(color_scale)){
    p <- p + ggplot2::scale_fill_manual(values=color_scale, na.value = 'white')
  }

  return(p)
}


#some example usage of pan2dist and pan2pcoa
# sel.genome <- sorted.ANI.cliques$genome
#
# pfam.dist <- pan2dist(pfam.pan[sel.genome,],
#                       dist.method="Manhattan")
# pfam.p.pcoa <- pan2pcoa(dist.pan= pfam.dist,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="pfam",
#                         color_scale=phylogroup.colors)


#---- RUN OVER DIFFERENT LEVELS OF ANI CLIQUES
#plot.range_pcoa
#=================================================================================#
#' Title
#'
#' @param dist.matrix
#' @param pan
#' @param dist.method
#' @param sorted.cliques
#' @param columns
#' @param shrink
#'
#' @return
#' @export
#'
#' @examples
#' p.list <- plot.range_pcoa(pan = OG.pan,
#'     dist.method='Manhattan',
#'     sorted.cliques=sorted.ANI.cliques,
#'     columns=c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
#'     shrink=FALSE)

plot.range_pcoa <- function(dist.matrix, pan, dist.method = "Manhattan", sorted.cliques, columns,shrink=TRUE, ellipse=FALSE){
  plotslist.t <- list()
  sel.genome <- rownames(pan)
  if(shrink==TRUE){
    for(colcol in columns){
      selcol <- colcol

      metdad <-genome.tbl %>%
        dplyr::left_join(sorted.cliques,by=c('genome'='genome')) %>%
        dplyr::left_join(out.colors, by=c('genome'='genome')) %>%
        dplyr::filter(genome %in% sel.genome) %>%
        dplyr::filter(!is.na(!!as.name(selcol)))%>%
        data.frame()
      rownames(metdad) <- metdad$genome
      metdad$group <- paste0('cl',metdad[,colcol])
      #metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))-1)),'clNA'))

      metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))))))

      colors <- metdad%>%
        dplyr::select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>%
        #filter(group!='clNA') %>%
        unique() %>%
        data.frame()

      grid.colors <- as.character(colors[,1])
      names(grid.colors) <- as.character(colors[,'group'])
      grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
      grid.colors <- c(grid.colors,'clNA'='white')

      dist.matrix <- pan2dist(pan=pan[metdad$genome,],dist.method=dist.method)

      if(ellipse==FALSE){
        plotslist.t[[colcol]] <- pan2pcoa(dist.pan = dist.matrix,
                                          metadata=metdad[metdad$genome,],
                                          color='group',
                                          title_suffix=colcol,
                                          color_scale=grid.colors,
                                          discrete=TRUE)
      }else{
        plotslist.t[[colcol]] <- pan2pcoa(dist.pan = dist.matrix,
                                          metadata=metdad[metdad$genome,],
                                          color='group',
                                          title_suffix=colcol,
                                          color_scale=grid.colors,
                                          discrete=TRUE,
                                          ellipse=TRUE)
        }
      }
  }else{
    metdad <-genome.tbl %>%
      dplyr::left_join(sorted.cliques,by=c('genome'='genome')) %>%
      dplyr::left_join(out.colors, by=c('genome'='genome')) %>%
      dplyr::filter(genome %in% sel.genome) %>%
      data.frame()

    dist.matrix <- pan2dist(pan=pan[metdad$genome,],dist.method=dist.method)

    for(colcol in columns){
      metdad <-genome.tbl %>%
        dplyr::left_join(sorted.cliques,by=c('genome'='genome')) %>%
        dplyr::left_join(out.colors, by=c('genome'='genome')) %>%
        dplyr::filter(genome %in% sel.genome) %>%
        data.frame()

      metdad$group <- paste0('cl',metdad[,colcol])

      rownames(metdad) <- metdad$genome
      metdad$group <- paste0('cl',metdad[,colcol])
      #metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))-1)),'clNA'))

      metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))))))


      colors <- metdad%>%
        dplyr::select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>%
        dplyr::filter(group!='clNA') %>%
        unique() %>%
        data.frame()

      grid.colors <- as.character(colors[,1])
      names(grid.colors) <- as.character(colors[,'group'])
      grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
      grid.colors <- c(grid.colors,'clNA'='white')

      metdad <- metdad %>% arrange(group)
      rownames(metdad) <- metdad$genome

      if(ellipse==FALSE){
        plotslist.t[[colcol]] <- pan2pcoa(dist.pan = dist.matrix,
                                          metadata=metdad,
                                          color='group',
                                          title_suffix=colcol,
                                          color_scale=grid.colors,
                                          discrete=TRUE)
      }else{
        plotslist.t[[colcol]] <- pan2pcoa(dist.pan = dist.matrix,
                                          metadata=metdad,
                                          color='group',
                                          title_suffix=colcol,
                                          color_scale=grid.colors,
                                          discrete=TRUE,
                                          ellipse=TRUE)
      }

    }
  }
  return(plotslist.t)
}




#=================================================================================#



#=================================================================================#
# tSNE - graph based, non-linear, stochastic method, only local distances preserved, distances bewteen groups not meaningful

# library(Rtsne) # Load package
# set.seed(42) # Sets seed for reproducibility


# pan2tsne
#=================================================================================#

#' Function that plots a t-SNE plot based on a pan table
#' uses function Rtsne from library(Rtsne)
#' @param pan
#' @param autotune boolean to set auto tuning, uses optimal nr of PCAs and optimal perplexity
#' @param perplexity
#' @param theta
#' @param metadata
#' @param color
#' @param title_suffix
#' @param color_scale
#' @param discrete
#'
#' @usage pan2tsne(pan, perplexity=15, theta=0.5, metadata=metadata, color="group",title_suffix="")
#'
#' @importFrom Rtsne Rtnse
#' @return
#' @export
#'
#' @examples
#' tcdb.p.tsne <- pan2tsne(tcdb.pan,perplexity=7, theta=0.5,metadata=metadata,color="group",title_suffix="TCDB")

# auto tuning is according to some simple principles
# PCA used to assess optimal dimensions
# # https://towardsdatascience.com/how-to-tune-hyperparameters-of-tsne-7c0596a18868

pan2tsne <- function(pan,autotune=FALSE, perplexity=15, theta=0.5, metadata=metadata, color="group",title_suffix="", color_scale=NULL,discrete=FALSE){
  set.seed(123456)
  if(is.factor(metadata[,color])){
    colorGroups = as.character(metadata[,color])
  }else{
    colorGroups = metadata[,color]
  }
  names(colorGroups) = rownames(metadata)

  if(autotune==FALSE){
    tsne_model_1 = Rtsne::Rtsne(as.matrix(pan), check_duplicates=FALSE, pca=TRUE, perplexity=perplexity, theta=theta, dims=2)
  }else{
    message('assessing optimal dimensions')
    PC <- prcomp(t(pan), center=TRUE, scale=FALSE)
    expl_var <- PC$sdev^2/sum(PC$sdev^2)
    N_perm <- 25
    expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = N_perm)
    for(k in 1:N_perm)
    {
      pan_perm <- apply(pan,2,sample)
      PC_perm <- prcomp(t(pan_perm), center=TRUE, scale=FALSE)
      expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
    }
    pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
    optPC<-head(which(pval>=0.05),1)-1

    N_cells<-dim(pan)[1]
    optimal_perplexity = round(sqrt(N_cells),0)
    if(optPC>2){
      message(paste0("Optimal number of principal components = ",optPC))
      tsne_model_1<-Rtsne::Rtsne(pan,
                                 perplexity=optimal_perplexity,
                                 initial_dims=optPC,
                                 max_iter=10000)
    }else{
      if(optPC %in% c(1,2)){
        message(paste0("Optimal number of principal components = ",optPC))
        message('this seems too low, so we continue with default')
        tsne_model_1<-Rtsne::Rtsne(pan,
                                   perplexity=optimal_perplexity,
                                   max_iter=10000)
      }else{
        message('failed to detect optimal number of PCs, continue with default')
        tsne_model_1<-Rtsne::Rtsne(pan,
                                   perplexity=optimal_perplexity,
                                   max_iter=10000)
      }

    }

  }


  d_tsne_1 = as.data.frame(tsne_model_1$Y)
  d_tsne_1$Name = rownames(pan)
  d_tsne_1$Group = colorGroups[d_tsne_1$Name]
  colnames(d_tsne_1) = c("x","y","Name","Group")

  if(discrete==FALSE){
    p = ggplot2::ggplot(d_tsne_1,aes(x=x,y=y,fill=Group,label=Name))
  }else{
    p = ggplot2::ggplot(d_tsne_1,aes(x=x,y=y,fill=as.factor(Group),label=Name))
  }

  p <- p +
    ggplot2::theme_classic() +
    ggplot2::xlab(paste("tSNE1",sep="")) +
    ggplot2::ylab(paste("tSNE2",sep="")) +
    ggplot2::ggtitle(label=paste0(title_suffix))+
    ggplot2::geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_point(shape = 21, size = 2) +
    fdb_style() +
    ggplot2::theme(legend.position = "none")

  if(!is.null(color_scale)){
    p <- p +
      ggplot2::scale_fill_manual(values=color_scale, na.value = 'white')
  }
  return(p)
}

#usage



#=================================================================================#
#' filter all tables in a panlist obtject to contain only specified genomes, and relevant orthogroups
#'
#' @param panlist
#' @param genomes
#'
#' @return list of feature abundance tables/pan-tables with only the selected genomes
#' @export
#'
#' @examples
#'
panlist_filter <- function(panlist, genomes, clean = TRUE){
  outpanlist<-list()
  for(pan.name in names(pan.list)){
    if(clean==TRUE){
      outpanlist[[pan.name]] <-  pan_clean(pan=pan.list[[pan.name]][genomes,],remove_singletons=FALSE)

    }else{
      outpanlist[[pan.name]] <-  pan.list[[pan.name]][genomes,]

    }
  }
  return(outpanlist)
}




#-----------------------------------------------------------#
# UMAP
# non-linear graph-based dimension reduction method (like tSNE)
# Based on topological structure in multidimensional space
# unlike tSNE you can compute the structure at once (no randomisation),
# perserved global structrure better than tSNE
#----------------------------------------------------#
#library(umap)

pan2umap <- function(pan, metadata=metadata, color="phylogroup",title_suffix="", color_scale=NULL){
  if(is.factor(metadata[,color])){
    colorGroups = as.character(metadata[,color])
  }else{
    colorGroups = metadata[,color]
  }
  names(colorGroups) = rownames(metadata)

  #umap
  pan.umap = umap::umap(pan)

  #wrangle data
  df.umap <- pan.umap$layout %>% data.frame() %>% rownames_to_column(var="genome")
  colnames(df.umap) <- c('genome','UMAP1','UMAP2')
  df.umap <- df.umap %>% mutate(Group=colorGroups[df.umap$genome]) #left_join(metadata, by='genome')

  #visualise
  p <- df.umap %>%
    ggplot2::ggplot(aes(x=UMAP1,y=UMAP2,fill=as.factor(Group),label=genome)) +
    ggplot2::geom_point(shape=21)

  #finish up
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::xlab(paste("UMAP1",sep="")) +
    ggplot2::ylab(paste("UMAP2",sep="")) +
    ggplot2::ggtitle(label=paste0(title_suffix))+
    ggplot2::geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_point(shape = 21, size = 2) +
    fdb_style() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values=color_scale)
  return(p)
}



## END OF FUNCTIONS ####
#### START OF ANALYSIS ####
genomes.selected <- metadata %>% filter(!sourceClass=='') %>% select(genome) %>% pull

d <- og.dist[genomes.selected,genomes.selected]

distance <- d  %>%
  data.frame() %>%
  tibble::rownames_to_column(var='genome')

meta_distance <- distance %>% left_join(metadata,by='genome')

all_dist <- meta_distance %>%
  dplyr::select(dplyr::all_of(.[['genome']])) %>%
  as.dist() %>% as.matrix()

mod1 = vegan::adonis(all_dist ~ sourceClass, data = meta_distance, permutation=10000)


#=================================================================================#
### RUN ####

#### 1. PCoA on various grouping levels ####

og.dist <- pan2dist(OG.pan[sel.genome,], dist.method="Manhattan")

og.p <- pan2pcoa(dist.pan = og.dist,
         metadata=metadata,
         color="phylogroup",
         title_suffix="OG") + scale_fill_manual(values=phylogroup.colors)

#with shrink = FALSE
p.list <- plot.range_pcoa(pan = OG.pan,
                          dist.method='Manhattan',
                          sorted.cliques=sorted.ANI.cliques,
                          columns=c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
                          shrink=FALSE)

p.ani.pcoa <-gridExtra::grid.arrange(grobs = p.list, ncol = 5)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANIclique_range_PCOA_OG.pdf',plot=p.ani.pcoa, width = 12, height = 4,unit='in')


#with shrink = TRUE
p.list <- plot.range_pcoa(pan = OG.pan,
                          dist.method='Manhattan',
                          sorted.cliques=sorted.ANI.cliques,
                          columns=c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
                          shrink=TRUE)

p.ani.pcoa.shrunk <-gridExtra::grid.arrange(grobs = p.list, ncol = 5)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANIclique_shrunk_PCOA_OG.pdf',plot=p.ani.pcoa, width = 12, height = 4,unit='in')

#combine
p.ani.pcoa.combined <-grid.arrange(p.ani.pcoa, p.ani.pcoa.shrunk, ncol = 1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANIclique_combined_PCOA_OG.pdf',plot=p.ani.pcoa.combined, width = 12, height = 6,unit='in')


#add phylogroup level plot
p.list[['PG']]=og.p
p.ani.pcoa <-gridExtra::grid.arrange(grobs =  p.list[c('PG','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')], ncol = 6)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANIclique_range_PCOA_OG_pg.pdf',plot=p.ani.pcoa, width = 12, height = 4,unit='in')



### 2. in specific lineage ####
#select some genomes
# for example those in broad lineage P5, and cicking out the ones that done group at 0.8 ANI
selected.g <- metadata %>%
  dplyr::filter(phylogroup=='P5') %>%
  dplyr::filter(!is.na(sc0.8)) %>%
  dplyr::select(genome) %>% pull

#filter the pan.list object
pan.list.filtered <- panlist_filter(panlist= pan.list, genomes=selected.g)


plotslist.t<-list()
for(pan.name in names(pan.list)){
  pan <- pan.list[[pan.name]]

  og.dist <- pan2dist(pan[sel.genome,], dist.method="Manhattan")

  og.p <- pan2pcoa(dist.pan = og.dist,
                   metadata=metadata,
                   color="phylogroup",
                   title_suffix=pan.name) + scale_fill_manual(values=phylogroup.colors)


  p.list <- plot.range_pcoa(pan = pan,
                            dist.method='Manhattan',
                            sorted.cliques=sorted.ANI.cliques,
                            columns=c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
                            shrink=TRUE)

  p.list[['PG']]=og.p
  p.all.pan <-gridExtra::grid.arrange(grobs = p.list[c('PG','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')], ncol = 6)
  plotslist.t[[pan.name]] <- p.all.pan
}

p.all.pan.all.clique <-gridExtra::grid.arrange(grobs = plotslist.t, nrow = length(names(pan.list)))
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_PCOA_all_clique_PG.pdf',plot=p.all.pan.all.clique, width = 12, height = 12,unit='in')


# Same thing, could be function, but here for the filtered one
# I had to remove the >0.98 similarity, as this is not present for many
plotslist.t<-list()
for(pan.name in names(pan.list.filtered)){
  pan <- pan.list.filtered[[pan.name]]
  p.list <- plot.range_pcoa(pan = pan,
                            dist.method='Manhattan',
                            sorted.cliques=sorted.ANI.cliques,
                            columns=c('sc0.8','sc0.85','sc0.9','sc0.95'),
                            shrink=TRUE)
  p.all.pan <-grid.arrange(grobs = p.list, ncol = 4)
  plotslist.t[[pan.name]] <- p.all.pan
}

p.all.pan.all.clique <-grid.arrange(grobs = plotslist.t, nrow = length(names(pan.list)))
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_PCOA_all_clique_filtered.pdf',plot=p.all.pan.all.clique, width = 12, height = 12,unit='in')



#=================================================================================#
### 3. per pan on specific grouping level ####
## integrating tibble based storage of all intermediate data
# idea is to store things in lists, and store the lists in analysis specific tibble
# this then can be used at later stages to be analysed, say distance matrices and metadatasets, can be used for permanova analysis


selcol <- 'sc0.8'
colcol <- 'cc0.8'
metdad <- genome.tbl %>%
  filter(genome%in% sorted.cliques$genome) %>%
  left_join(sorted.cliques,by='genome') %>%
  filter(grepl('Marinobacter|Tamil',genome)) %>%
  left_join(out.colors, by='genome') %>%
  filter(quality_class!='low')%>%
  filter(!is.na(!!as.name(selcol))) %>%
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
  pan <- pan.list[[pan.name]]
  pan.dist <- pan2dist(pan[rownames(metdad),],
                        dist.method="Manhattan")
  pan.p.pcoa <- pan2pcoa(dist.pan= pan.dist,
                          metadata=metdad,
                          color="group",
                          title_suffix=pan.name,
                          color_scale=clique.colors)
  plotslist.t[[pan.name]] <- pan.p.pcoa
  }

p.all.pan <-grid.arrange(grobs = plotslist.t, ncol = 6)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_PCOA.pdf',plot=p.all.pan, width = 12, height = 4,unit='in')



#=================================================================================#
#### Run tSNE ####
#### 1. on all pan tables ####

pan.list <- list('OG'=OG.pan[sel.genome,],
                 'KO'=kfm.pan[sel.genome,],
                 'COG'=COG.pan[sel.genome,],
                 'pfam'=pfam.pan[sel.genome,],
                 'CAZy'=cazy.pan[sel.genome,],
                 'TCDB'=tcdb.pan[sel.genome,])

selcol <- 'sc0.85'
colcol <- 'cc0.85'
metdad <- genome.tbl %>%
  filter(genome%in% sorted.cliques$genome) %>%
  left_join(sorted.cliques,by='genome') %>%
  filter(grepl('Marinobacter|Tamil',genome)) %>%
  left_join(out.colors, by='genome') %>%
  filter(quality_class!='low') %>%
  filter(!is.na(!!as.name(selcol))) %>%
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
  pan <- pan.list[[pan.name]]

  pan.p.tsne <- pan2tsne(unique(pan[sel.genome,]),
                         autotune=TRUE,
                         metadata=metadata,
                         color="phylogroup",
                         title_suffix=pan.name,
                         color_scale= phylogroup.colors,
                         discrete=TRUE)
  plotslist.t[[pan.name]] <- pan.p.tsne
}



p.all.pan.tsne <-gridExtra::grid.arrange(grobs = plotslist.t, ncol = 6)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_tSNE.pdf',plot=p.all.pan.tsne, width = 12, height = 4,unit='in')

#=================================================================================#


#---- for each pangenome using PcOA
plotslist.t <- list()
for(pan.name in names(pan.list)){
  pan <- pan.list[[pan.name]]

  og.dist <- pan2dist(pan[sel.genome,], dist.method="Manhattan")

  og.p <- pan2pcoa(dist.pan = og.dist,
                   metadata=metadata,
                   color="phylogroup",
                   title_suffix=pan.name) + scale_fill_manual(values=phylogroup.colors)

  plotslist.t[[pan.name]] <- og.p
}

p.all.pan.pcoa <-gridExtra::grid.arrange(grobs = plotslist.t, ncol = 6)


#combine plots from PCoA
p.tsne.combined <- gridExtra::grid.arrange(OGpan.p.tsne, COGpan.p.tsne, kfam.p.tsne, pfam.p.tsne,tcdb.p.tsne ,nrow=1)
p.tsne.combined


p.combined.dimred <- gridExtra::grid.arrange(p.pcoa.combined, p.tsne.combined, ncol=1)
p.combined.dimred



#### UMAP ####
#### 1. on all pan ####

pan2umap(pan=OG.pan[sel.genome,],  metadata=metadata,title_suffix = "UMAP OG", color_scale = phylogroup.colors)


#sel.genome <- sco.94.concat.nuc.tree$tip.label

plotslist.t <- list()
for(pan.name in names(pan.list)){
  pan <- pan.list[[pan.name]]

  pan.p.umap <- pan2umap(pan=pan[sel.genome,],
                         metadata=metadata,
                         title_suffix = pan.name,
                         color_scale = phylogroup.colors)
  plotslist.t[[pan.name]] <- pan.p.umap
}

p.all.pan.umap <-gridExtra::grid.arrange(grobs = plotslist.t, ncol = 6)
p.all.pan.umap

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_umap.pdf',plot=p.all.pan.umap, width = 12, height = 4,unit='in')



p.combined.dimred <- gridExtra::grid.arrange(p.all.pan.pcoa,p.all.pan.umap, p.all.pan.tsne, ncol=1)
p.combined.dimred

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_dimred.pdf',plot=p.combined.dimred, width = 12, height = 8,unit='in')






#--------------------#


selcol <- 'sc0.85'
colcol <- 'sc0.85'

metdad <-genome.tbl %>%
  dplyr::left_join(sorted.cliques,by=c('genome'='genome')) %>%
  dplyr::left_join(out.colors, by=c('genome'='genome')) %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(!is.na(!!as.name(selcol)))%>%
  data.frame()
rownames(metdad) <- metdad$genome
metdad$group <- paste0('cl',metdad[,colcol])
#metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))-1)),'clNA'))

metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))))))

colors <- metdad%>%
  dplyr::select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>%
  #filter(group!='clNA') %>%
  unique() %>%
  data.frame()

grid.colors <- as.character(colors[,1])
names(grid.colors) <- as.character(colors[,'group'])
grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
grid.colors <- c(grid.colors,'clNA'='white')


shrunk.genome <- metdad$genome

#the filtered sets
p.cutoff <- 0.05
selcol <- 'sc0.85'

sig_enriched_binary <- ANI.hmg.grp.tbl %>% filter(fdr.p.value_enriched_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_enriched_raw <- ANI.hmg.grp.tbl %>% filter(fdr.p.value_enriched_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_depleted_binary <- ANI.hmg.grp.tbl %>% filter(fdr.p.value_depletion_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_depleted_raw <- ANI.hmg.grp.tbl %>% filter(fdr.value_depletion_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()

extracte.ogs <- c(sig_enriched_raw, sig_depleted_binary , sig_enriched_binary, sig_depleted_raw) %>% unique()

OGpan.p.umap.sig <- pan2umap(pan=OG.pan[sel.genome,extracte.ogs],
                             metadata=metadata,
                             title_suffix = "UMAP",
                             color_scale = phylogroup.colors)

OGpan.p.tsne.sig <- pan2tsne(OG.pan[sel.genome,extracte.ogs],
                         autotune=TRUE,
                         metadata=metadata,
                         color="phylogroup",
                         title_suffix="tSNE",
                         color_scale= phylogroup.colors,
                         discrete=TRUE)

OGpan.dist <- pan2dist(OG.pan[sel.genome,extracte.ogs],
                      dist.method="Manhattan")
OGpan.p.pcoa.sig <- pan2pcoa(dist.pan= OGpan.dist,
                        metadata=metadata,
                        color="phylogroup",
                        title_suffix="PCoA",
                        color_scale=phylogroup.colors)
p.OG.pan.dimred <-gridExtra::grid.arrange(OGpan.p.pcoa.sig,OGpan.p.umap.sig,OGpan.p.tsne.sig, ncol = 3)


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_sig_only_dim_red.pdf',plot=p.OG.pan.dimred, width = 8, height = 4,unit='in')

#--------------------


OGpan.p.umap.sig <- pan2umap(pan=OG.pan[shrunk.genome,extracte.ogs],
                             metadata=metdad,
                             color='group',
                             title_suffix = "UMAP",
                             color_scale = grid.colors)

OGpan.p.tsne.sig <- pan2tsne(OG.pan[shrunk.genome,extracte.ogs],
                             autotune=TRUE,
                             metadata=metdad,
                             color="group",
                             title_suffix="tSNE",
                             color_scale= grid.colors,
                             discrete=TRUE)

OGpan.dist <- pan2dist(OG.pan[shrunk.genome,extracte.ogs],
                       dist.method="Manhattan")

OGpan.p.pcoa.sig <- pan2pcoa(dist.pan= OGpan.dist,
                             metadata=metdad,
                             color="group",
                             title_suffix="PCoA",
                             color_scale=grid.colors)

p.OG.pan.dimred_sig <-grid.arrange(OGpan.p.pcoa.sig,OGpan.p.umap.sig,OGpan.p.tsne.sig, ncol = 3)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_sig_only_dim_red_significant.pdf',plot=p.OG.pan.dimred_sig, width = 8, height = 4,unit='in')



OGpan.p.umap.sig <- pan2umap(pan=OG.pan[shrunk.genome,],
                             metadata=metdad,
                             color='group',
                             title_suffix = "UMAP",
                             color_scale = grid.colors)

OGpan.p.tsne.sig <- pan2tsne(OG.pan[shrunk.genome,],
                             autotune=TRUE,
                             metadata=metdad,
                             color="group",
                             title_suffix="tSNE",
                             color_scale= grid.colors,
                             discrete=TRUE)

OGpan.dist <- pan2dist(OG.pan[shrunk.genome,],
                       dist.method="Manhattan")

OGpan.p.pcoa.sig <- pan2pcoa(dist.pan= OGpan.dist,
                             metadata=metdad,
                             color="group",
                             title_suffix="PCoA",
                             color_scale=grid.colors)

p.OG.pan.dimred_all <-grid.arrange(OGpan.p.pcoa.sig,OGpan.p.umap.sig,OGpan.p.tsne.sig, ncol = 3)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_sig_only_dim_red_shrunk_allgenes.pdf',plot=p.OG.pan.dimred_sig, width = 8, height = 4,unit='in')




p.both <- grid.arrange(p.OG.pan.dimred_all,p.OG.pan.dimred_sig, ncol = 1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_sig0.85_only_dim_red_shrunk_allgenes_both.pdf',plot=p.both, width = 8, height = 8,unit='in')



p.OG.pan.dimred_sig









grid.arrange(OGpan.p.pcoa.sig,OGpan.p.tsne.sig, ncol = 3)

















OGpan.p.umap <- pan2umap(pan=OG.pan[sel.genome,],
                             metadata=metadata,
                             title_suffix = "UMAP",
                             color_scale = phylogroup.colors)

OGpan.p.tsne <- pan2tsne(OG.pan[sel.genome,],
                             autotune=TRUE,
                             metadata=metadata,
                             color="phylogroup",
                             title_suffix="tSNE",
                             color_scale= phylogroup.colors,
                             discrete=TRUE)

OGpan.dist <- pan2dist(OG.pan[sel.genome,],
                       dist.method="Manhattan")
OGpan.p.pcoa <- pan2pcoa(dist.pan= OGpan.dist,
                             metadata=metadata,
                             color="phylogroup",
                             title_suffix="PCoA",
                             color_scale=phylogroup.colors)
p.OG.pan.dimred2 <-grid.arrange(OGpan.p.pcoa,OGpan.p.umap,OGpan.p.tsne, ncol = 3)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_all_dim_red.pdf',plot=p.OG.pan.dimred, width = 7, height = 4,unit='in')




OGpan.p.umap <- pan2umap(pan=pan_clean(OG.pan[sel.genome,],remove_singletons = TRUE),
                         metadata=metadata,
                         title_suffix = "UMAP",
                         color_scale = phylogroup.colors)

OGpan.p.tsne <- pan2tsne(pan_clean(OG.pan[sel.genome,],remove_singletons = TRUE),
                         autotune=TRUE,
                         metadata=metadata,
                         color="phylogroup",
                         title_suffix="tSNE",
                         color_scale= phylogroup.colors,
                         discrete=TRUE)

OGpan.dist <- pan2dist(pan_clean(OG.pan[sel.genome,],remove_singletons = TRUE),
                       dist.method="Manhattan")
OGpan.p.pcoa <- pan2pcoa(dist.pan= OGpan.dist,
                         metadata=metadata,
                         color="phylogroup",
                         title_suffix="PCoA",
                         color_scale=phylogroup.colors)
p.OG.pan.dimred3 <-gridExtra::grid.arrange(OGpan.p.pcoa,OGpan.p.umap,OGpan.p.tsne, ncol = 3)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_all_dim_red_pan_clean.pdf',plot=p.OG.pan.dimred, width = 7, height = 4,unit='in')




p.OG.pan.dimred_all <-gridExtra::grid.arrange(p.OG.pan.dimred,p.OG.pan.dimred2,p.OG.pan.dimred3, ncol = 1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ogpan_all__sigonly_all_and_cleanall.pdf',plot=p.OG.pan.dimred_all, width = 7, height = 8,unit='in')

p.OG.pan.dimred_all <-grid.arrange(p.all.pan,p.all.pan.umap,p.all.pan.tsne, ncol = 1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/all_pan_dimred.pdf',plot=p.OG.pan.dimred_all, width = 7, height = 8,unit='in')



#-----------------------------------------------------------------------#
#### Effect of enrichment analysis on PCoA ####






df.all.pcoa <- c()
#entire pangenome
og.dist <- pan2dist(out.dist.tibble$pan[[3]], dist.method="Manhattan")
p.o <- pan2pcoa(dist.pan = og.dist,
                metadata=out.dist.tibble$metadata[[3]],
                color="group",
                title_suffix="OG") + scale_fill_manual(values=col.list[['sc0.85']])

for(s.grp in unique(rf.all$group)){
#for(s.grp in c('cl1')){
  tmp.df <- p.o$data %>% mutate(selgroup=s.grp) %>% mutate(method='all')
  df.all.pcoa <- rbind(df.all.pcoa, tmp.df )

  message(s.grp)
  all250topRF <- rf.all %>% filter(group==s.grp) %>%
    group_by(group) %>%
    arrange(desc(MeanDecreaseAccuracy)) %>%
    slice(1:250) %>%
    select(feature.id) %>%
    pull()
  message(length(all250topRF))

  sigcheck <- HGM_RF_COMB %>%
    dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts,fdr.p.value_depletion_binary,fdr.value_depletion_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
    dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Depl_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Depl_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::select(feature.id.x, group.x, Enrich_bin,Depl_bin,Depl_raw, Enrich_raw) %>%
    dplyr::mutate(signEnriched = ifelse(Enrich_bin == TRUE | Enrich_raw==TRUE, TRUE, FALSE)) %>%
    dplyr::mutate(signDepleted = ifelse(Depl_bin == TRUE | Depl_raw==TRUE, TRUE, FALSE))


  lssig <- list()
  lssig[['RF']] <- all250topRF
  lssig[['sig_en_bin']] <- sigcheck %>% filter(group.x==s.grp) %>%
    dplyr::filter(Enrich_bin==TRUE) %>% select(feature.id.x) %>% pull() %>% unique() %>% as.character()
  lssig[['sig_en_raw']] <- sigcheck %>% filter(group.x==s.grp) %>%
    dplyr::filter(Enrich_raw==TRUE) %>% select(feature.id.x) %>% pull() %>% unique() %>% as.character()
  lssig[['sig_de_bin']] <- sigcheck %>% filter(group.x==s.grp) %>%
    dplyr::filter(Depl_bin==TRUE) %>% select(feature.id.x) %>% pull() %>% unique() %>% as.character()
  lssig[['sig_de_raw']] <- sigcheck %>% filter(group.x==s.grp) %>%
    dplyr::filter(Depl_raw==TRUE) %>% select(feature.id.x) %>% pull() %>% unique() %>% as.character()
  lssig[['allsig']] <- unique(unlist(lssig))




  for(sg in names(lssig)){
    if(length(lssig[[sg]]) > 10){
      og.dist <- pan2dist(out.dist.tibble$pan[[3]][,lssig[[sg]]], dist.method="Manhattan")
      p.or <- pan2pcoa(dist.pan = og.dist,
                      metadata=out.dist.tibble$metadata[[3]],
                      color="group",
                      title_suffix="OG") + scale_fill_manual(values=col.list[['sc0.85']])

      tmp.df <- p.or$data %>% mutate(selgroup=s.grp) %>% mutate(method=sg)
      df.all.pcoa <- rbind(df.all.pcoa, tmp.df )

    }
  }
}


df.all.pcoa <- df.all.pcoa %>%
  mutate(method=factor(method, levels=c('all','RF','allsig','sig_en_bin','sig_en_raw','sig_de_bin','sig_de_raw')))



plotlist <- list()
for(s.grp in unique(rf.all$group)){
  plotlist[[s.grp]] <- df.all.pcoa %>%
    filter(selgroup==s.grp) %>%
    mutate(ingroup=ifelse(Group==s.grp,s.grp,NA)) %>%
    ggplot2::ggplot(aes(x=x,y=y)) +
    ggplot2::theme_classic() +
    #ggplot2::xlab(paste("PCoA1 (", MDS.var.perc[1], "%)",sep="")) +
    #ggplot2::ylab(paste("PCoA2 (", MDS.var.perc[2], "%)",sep="")) +
    ggplot2::ggtitle(label=s.grp)+
    ggplot2::geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
    ggplot2::geom_point(shape = 21, size = 2,aes(fill=ingroup)) +
    ggforce::geom_mark_ellipse(aes(color=ingroup,fill=ingroup),alpha=.5) + #ggplot2::stat_ellipse(aes(color=ingroup)) +
    fdb_style() +
    ggplot2::theme(legend.position = "none") +
    facet_wrap(.~method,scales='free',nrow=1) +
    ggplot2::scale_fill_manual(values=col.list[['sc0.85']], na.value = 'grey') +
    ggplot2::scale_color_manual(values=col.list[['sc0.85']], na.value = 'grey')

}

#add ellipse?



p.all <-gridExtra::grid.arrange(grobs = plotlist, ncol =1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_all_associated2.pdf',plot=p.all, width = 16, height = 16,unit='in')


#==============================================================================#



intraPG <- interaDistances(df.dist = out.dist.tibble$distance[[1]],
                           var_name = "ANIb",
                           metadata = out.dist.tibble$metadata[[1]],
                           grp = "phylogroup")

p.phl <- intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(interaction(grouping,group.x), value)) +
  geom_jitter(aes(fill=group.x),shape=21,size=1,width=0.2,alpha=0.7,color='white') +
  geom_boxplot(aes(fill=group.x),outlier.shape=NA,alpha=0.7) +
  #facet_wrap(~grouping)+
  scale_color_manual(values=phylogroup.colors) +
  scale_fill_manual(values=phylogroup.colors) +
  fdb_style(aspect.ratio=0.75)

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/distances_test.pdf',plot=p.phl, width = 6, height = 6,unit='in')


#==============================================================================#




out.dist.tibble$metadata[[3]]


vegan::adonis2(as.dist(out.dist.tibble$distance[[3]]) ~ as.formula("group + sourceClass + Completeness + Contamination + Genome_size"), data = out.dist.tibble$metadata[[3]],permutation=prm)







allSig_combined <- unique(c(all250topRF,allSig_hgm))






#------------------------------#
#### ADDITIONAL CODE ####
#------------------------------#

# # Dimensionality reduction
# #---------------------------------------------------------------------#
# pan2metaMDS <- function(dist,trymax=trymax){
#   out.nmds <-  vegan::metaMDS(as.matrix(dist),trymax = trymax)
#   data.scores <- as.data.frame(vegan::scores(out.nmds))
#   data.scores$site <- rownames(data.scores)
#   return(data.scores)
# }
#
#
# MDSdata <- pan2metaMDS(ANIb[sel.genome,sel.genome], trymax=500)
#
# p = ggplot2::ggplot(MDSdata,aes(x=NMDS1,y=NMDS2))+
#   ggplot2::theme_classic() +
#   geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
#   geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
#   geom_point(shape = 21, size = 2) +
#   fdb_style() + theme(legend.position = "none")
#
# p
#
# vegan::scores(out.nmds) %>%
#   cbind(ANIb[sel.genome,sel.genome]) %>%
#   ggplot(aes(x = NMDS1, y = NMDS2)) +
#   geom_point(size=2,shape=21)+
#   ggtitle(paste0("stress: ", format(out.nmds$stress, digits = 4)))+
#   fdb_style()
#
# vegan::scores(out.nmds) %>%
#   cbind(ANIb[sel.genome,sel.genome]) %>%
#   ggplot(aes(x = NMDS1, y = NMDS2)) +
#   geom_point(aes(size = Richness, color = Group)) +
#   stat_ellipse(geom = "polygon", aes(group = Group, color = Group, fill = Group), alpha = 0.3) +
#   annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
#   theme_bw()
#
#
#
# dstMan.nmds <- vegan::metaMDS(distMan[rownames(ANIb),rownames(dstANI)],trymax = 500)
# coph.nmds <- vegan::metaMDS(COPH[rownames(ANIb),rownames(dstANI)],trymax = 500)
#
# ani.coph.procrust <- vegan::protest(ani.nmds,dstMan.nmds, symmetric=T)
# ani.coph.procrust <- vegan::protest(ani.nmds,coph.nmds, symmetric=T)
# ani.coph.procrust <- vegan::protest(coph.nmds, ani.nmds, symmetric=T)
#
# summary(ani.coph.procrust)
# #plot(ani.coph.procrust)
#
# ## Re-plotting the procrustes
# prodat <- as.data.frame(ani.coph.procrust$X)
# prodat <- cbind(prodat,ani.coph.procrust$Yrot,metadata[rownames(prodat),])
#
# colnames(prodat)[colnames(prodat)=="1"] <- "Xst"
# colnames(prodat)[colnames(prodat)=="2"] <- "Yst"
#
# #Shape 19 or 21 as first nmds (COPH) #I THINK TGHEY GOT IT WRONG
# #shape 17 or 24 as second nmds (ANI)
#
# ggplot() +
#   geom_point(data=prodat, mapping=aes(x=NMDS1, y=NMDS2,fill=group), size=2, shape=19,alpha=0.7) +#21
#   geom_point(data=prodat, mapping=aes(x=Xst, y=Yst,fill=group), size=2, shape=17,alpha=0.7) +#24
#   geom_segment(data=prodat, mapping=aes(x=Xst, y=Yst, xend=NMDS1, yend=NMDS2,color=group), size=0.5, arrow = arrow(length =  unit(0.03, "npc"))) +
#   labs(x="Procrustes axis 1", y="Procrustes axis 2") +
#   theme_classic() + fdb_style()
#
# ##### Correlating composition and function ######
# #### Procrustes analysis comparing functional and compositional NMDS plots for D63
# ani.coph.procrust <- protest(expf63n.nmds,ecofxn63.nmds, symmetric=T)
# ani.coph.procrust ## Correlation in a symmetric Procrustes rotation: 0.8991, Significance:  0.001
#
# #### Procrustes analysis comparing functional and compositional NMDS plots for D63
# ani.coph.procrust <- protest(expf63n.nmds,ecofxn63.nmds, symmetric=T)
# ani.coph.procrust ## Correlation in a symmetric Procrustes rotation: 0.8991, Significance:  0.001
# summary(ani.coph.procrust)
# plot(ani.coph.procrust)

#------------#
#not automatic
# pfam.p.tsne <- pan2tsne(pfam.pan[sel.genome,],
#                         perplexity=7, theta=0.5,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="pfam",
#                         color_scale= phylogroup.colors,
#                         discrete=TRUE)
#
# #with hyperparameter tuning
# pfam.p.tsne <- pan2tsne(pan_clean(pfam.pan[sel.genome,],remove_singletons = TRUE),
#                         autotune=TRUE,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="pfam",
#                         color_scale= phylogroup.colors,
#                         discrete=TRUE)
#
# grid.arrange(pfam.p.tsne,pfam.p.tsne.auto,ncol=2)
#
# tcdb.p.tsne <- pan2tsne(pan_clean(tcdb.pan[sel.genome,],remove_singletons = TRUE),
#                         autotune=TRUE,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="TCDB",
#                         color_scale= phylogroup.colors,
#                         discrete=TRUE)
#
# kfam.p.tsne <- pan2tsne(pan_clean(kfm.pan[sel.genome,],remove_singletons = TRUE),
#                         autotune=TRUE,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="kfam",
#                         color_scale= phylogroup.colors,
#                         discrete=TRUE)
#
# cazy.p.tsne <- pan2tsne(pan_clean(cazy.pan[sel.genome,],remove_singletons = TRUE) %>% unique(),
#                         autotune=TRUE,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="CAZY",
#                         color_scale= phylogroup.colors,
#                         discrete=TRUE)
#
# COGpan.p.tsne <- pan2tsne(COG.pan[sel.genome,],
#                         autotune=TRUE,
#                         metadata=metadata,
#                         color="phylogroup",
#                         title_suffix="COG",
#                         color_scale= phylogroup.colors,
#                         discrete=TRUE)
#
# OGpan.p.tsne <- pan2tsne(OG.pan[sel.genome,],
#                           autotune=TRUE,
#                           metadata=metadata,
#                           color="phylogroup",
#                           title_suffix="OG",
#                           color_scale= phylogroup.colors,
#                           discrete=TRUE)
#
# OGpan.p.tsne

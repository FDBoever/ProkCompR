#' Title
#'
#' @param adonis.model
#'
#' @return
#' @export
#'
#' @examples
#' AICc_PERMANOVA(adonis_phylogeny)

AICc_PERMANOVA <- function(adonis.model) {
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
  # Calculating AICc using residual sum of squares (RSS)

  # AIC : 2*k + n*ln(RSS)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],

  RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
  MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
  k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
  nn <- nrow(adonis.model$model.matrix)


  AIC <- 2*k + nn*log(RSS)
  AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
  AIC.MSE <- 2*k + nn * log(MSE)
  AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
  AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
  AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
  AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)

  out <- list("AIC" = AIC,
              "AIC.g" = AIC.g,
              "AICc" = AICc,
              "AIC.MSE" = AIC.MSE,
              "AICc.MSE" = AICc.MSE,
              "AIC.pi" = AIC.pi,
              "AICc.pi" = AICc.pi,
              "k" = k)
  return(out)
}

#usage
#AICc_PERMANOVA(adonis_phylogeny)




### removed color subtraction - need to make a smart function to do this, and then annotate it later, if needed for plotting


# do the sorted.cliques attachment here??
#' Title
#'
#' @param pan
#' @param metadata
#' @param dist.method
#' @param sorted.cliques
#' @param columns
#' @param shrink
#'
#' @return
#' @export
#'
#' @examples
create_range_distances <- function(pan, metadata, dist.method = "Manhattan", sorted.cliques, columns, shrink=TRUE){
  l.gemome.names <- list()
  l.pan <- list()
  l.metadata <- list()
  l.distace <- list()

  sel.genome <- rownames(pan)

  if(shrink==TRUE){ #remove false has been removed here
    for(colcol in columns){
      selcol <- colcol
      metdad <-genome.tbl %>%
        dplyr::left_join(sorted.cliques,by=c('genome'='genome')) %>%
        dplyr::filter(genome %in% sel.genome) %>%
        dplyr::filter(!is.na(!!as.name(selcol)))%>%
        data.frame()

      rownames(metdad) <- metdad$genome

      if(selcol=='phylogroup'){
        metdad$group <- metdad$phylogroup
        }else{
        metdad$group <- paste0('cl',metdad[,colcol])
        metdad$group <- factor(metdad$group,levels = c(paste0('cl',1:(length(unique(metdad$group))))))
        }

      #colors <- metdad%>%
      #  dplyr::select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>%
      #  unique() %>%
      #  data.frame()

      #grid.colors <- as.character(colors[,1])
      #names(grid.colors) <- as.character(colors[,'group'])
      #grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
      #grid.colors <- c(grid.colors,'clNA'='white')

      pan.cleaned <- pan_clean(pan[metdad$genome,])
      dist.matrix <- pan2dist(pan=pan.cleaned,dist.method=dist.method)

      l.gemome.names[[colcol]] <- metdad$genome
      l.pan[[colcol]] <- pan.cleaned
      l.metadata[[colcol]] <- metdad
      l.distace[[colcol]] <- dist.matrix
    }
  }

  #add the lists into a tibble object
  out.tibble <- tibble::tibble(columns=columns,
                               genomes.set = l.gemome.names,
                               pan = l.pan,
                               metadata = l.metadata,
                               distance = l.distace)
  return(out.tibble)
}


#usage
distances.tibble <- create_range_distances(
                pan = OG.pan,
                dist.method='Manhattan',
                sorted.cliques=sorted.ANI.cliques,
                columns=c('phylogroup','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
                shrink=TRUE)


#adds column with the output of vegan::adonis2 from distance tibble
adonis2_from_tibble <- function(distances.tibble){
  l.adonis <- list()
  l.adonis.res <- list()

  for(i in 1:nrow(distances.tibble)){
    message('analysing', i)
    l.adonis[[i]]<-adonis_model(distances.tibble$distance[[i]], distances.tibble$metadata[[i]])
    l.adonis.res[[i]] = cbind('terms'=rownames(data.frame(l.adonis[[i]]))[1:nrow(data.frame(l.adonis[[i]]))-1],data.frame(l.adonis[[i]])[1:nrow(data.frame(l.adonis[[i]]))-1,]) %>% mutate(model='mod1')
  }

  out.tibble <- distances.tibble %>%
    dplyr::mutate(adonis = l.adonis,
                  res.adonis = l.adonis.res)

  return(out.tibble)
}


#run over all bootstraps and calculate the adonis models
#usage
distances.adonis.tibble <- adonis2_from_tibble(distances.tibble)



adonis_model <- function(pan,metadata){
  vegan::adonis2(pan ~ as.formula("group + sourceClass + Completeness + Contamination + Genome_size"), data = metadata,permutation=prm)
}




#### RUN ####

#======================================================#
# run over different pantables and stpre all the data
# this may take a while to run...



out.dist.tibble <- c()
for(i in names(panList)){
  message(i)
  pan = panList[[i]]
  distances.tibble <- create_range_distances(
    pan = pan,
    dist.method='Manhattan',
    sorted.cliques=sorted.ANI.cliques,
    columns=c('phylogroup','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
    shrink=TRUE)
  distances.adonis.tibble <- adonis2_from_tibble(distances.tibble)
  distances.adonis.tibble <- distances.adonis.tibble %>% mutate(pan_name=i)
  out.dist.tibble <- rbind(out.dist.tibble, distances.adonis.tibble)
  }





out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))) %>%
  ggplot2::ggplot(ggplot2::aes(x = columns, y = abs(R2),fill=terms))+
  ggplot2::geom_bar(stat="identity",position="fill")+
  ggplot2::scale_y_continuous(expand = c(0, 0))+
  ggplot2::theme_bw()+
  ggplot2::theme( panel.background = ggplot2::element_rect(colour = "black", size=0.5))+
  ggplot2::geom_text(aes(label=paste0(ifelse(Pr..F.<.5,'*',''),round(R2*100,2))),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  ggplot2::scale_fill_brewer(palette = "BuPu")+
  ggplot2::facet_wrap(~ pan_name,scales='free',nrow=1)+
  ggplot2::xlab(label='model')+ggplot2::ylab(label='variance explained (R2)')+
  ggplot2::coord_flip()+
  ggplot2::ggtitle(label='PERMANOVA (adonis) model comparison')+
  ggplot2::scale_x_discrete(limits=rev)


p1 <- out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))) %>%
  ggplot(aes(x = columns, y = abs(R2),fill=terms))+
  geom_bar(stat="identity",position="fill")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=paste0(ifelse(Pr..F.<.5,'*',''),round(R2*100,2))),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  scale_fill_brewer(palette = "BuPu")+
  facet_grid(.~pan_name)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison') + theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev)#+fdb_style()

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/PERMANOVA_pans_fullModel.pdf',plot=p1, width = 12, height = 4,unit='in')





p1 <- out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))) %>%
  ggplot(aes(x = columns, y =terms))+
  geom_point(aes(size= abs(R2)),color='grey')+
  geom_text(data=.%>% dplyr::filter(Pr..F.<.5), aes(label=abs(round(R2*100,2))),color='black',size=2) + #,stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  scale_fill_brewer(palette = "BuPu")+
  facet_grid(.~pan_name)+
  xlab(label=' ')+ylab(label=' ')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison') +
  fdb_style() +
  ggplot2::theme(aspect.ratio = 1,
                 axis.text.x=element_text(angle=45, hjust=1))+
  theme( panel.background = element_rect(colour = "black", size=1))+
  scale_x_discrete(limits=rev)

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/PERMANOVA_pans_fullModel_geom_point2.pdf',plot=p1, width = 12, height = 4,unit='in')






p2 <- out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::filter(columns != 'sc0.98') %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))) %>%
  dplyr::mutate(pan_name = factor(pan_name, levels=c('OGpan','COGpan','KOpan','PFAMpan','CAZYpan','TCDBpan'))) %>%
  ggplot(aes(x = columns, y = abs(R2),fill=terms))+
  geom_bar(stat="identity",position="fill")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=paste0(ifelse(Pr..F.<.5,'*',''),round(R2*100,2))),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  scale_fill_brewer(palette = "Greys")+
  facet_grid(.~pan_name)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison') + theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev)+fdb_style() +
  ggplot2::theme(aspect.ratio = 1,
                 axis.text.x=element_text(angle=45, hjust=1))+
  theme( panel.background = element_rect(colour = "black", size=1))+
  scale_x_discrete(limits=rev)


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/PERMANOVA_pans_fullModel_cleaner.pdf',plot=p2, width = 12, height = 4,unit='in')






p2 <- out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::filter(columns != 'sc0.98') %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))) %>%
  dplyr::mutate(pan_name = factor(pan_name, levels=c('OGpan','COGpan','KOpan','PFAMpan','CAZYpan','TCDBpan'))) %>%
  dplyr::filter(pan_name %in% c('OGpan','KOpan','PFAMpan')) %>%
  ggplot(aes(x = columns, y = abs(R2),fill=terms))+
  geom_bar(stat="identity",position="fill")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=paste0(ifelse(Pr..F.<.5,'*',''),round(R2*100,2))),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  scale_fill_brewer(palette = "Greys")+
  facet_grid(.~pan_name)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison') + theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev)+fdb_style() +
  ggplot2::theme(aspect.ratio = 1,
                 axis.text.x=element_text(angle=45, hjust=1))+
  theme( panel.background = element_rect(colour = "black", size=1))+
  scale_x_discrete(limits=rev)


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/PERMANOVA_pans_fullModel_cleaner_blues3.pdf',plot=p2, width = 8, height = 4,unit='in')








out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::filter(columns != 'sc0.98') %>%
  dplyr::filter(columns == 'phylogroup') %>%
  dplyr::mutate(terms = factor(terms, levels=rev(c('Residual','Contamination','Completeness','Genome_size','sourceClass','group')))) %>%
  dplyr::mutate(pan_name = factor(pan_name, levels=c('OGpan','COGpan','KOpan','PFAMpan','CAZYpan','TCDBpan'))) %>%
  dplyr::filter(pan_name %in% c('OGpan','KOpan','PFAMpan')) %>%
  ggplot(aes(x = terms, y = abs(R2),fill=terms))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=paste0(ifelse(Pr..F.<.5,'*',''),round(R2*100,2))),stat='identity',size=2.5,color='white')+
  scale_fill_brewer(palette = "Greys",direction = -1)+
  facet_grid(.~pan_name)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison') + theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev)+fdb_style() +
  ggplot2::theme(aspect.ratio = 1,
                 axis.text.x=element_text(angle=45, hjust=1))+
  theme( panel.background = element_rect(colour = "black", size=1))+
  scale_x_discrete(limits=rev)


p.perm  <- out.dist.tibble %>% tidyr::unnest(res.adonis) %>%
  dplyr::filter(columns != 'sc0.98') %>%
  dplyr::filter(columns == 'phylogroup') %>%
  dplyr::mutate(terms = factor(terms, levels=rev(c('Residual','Contamination','Completeness','Genome_size','sourceClass','group')))) %>%
  dplyr::mutate(pan_name = factor(pan_name, levels=c('OGpan','COGpan','KOpan','PFAMpan','CAZYpan','TCDBpan'))) %>%
  dplyr::filter(pan_name %in% c('OGpan','KOpan','PFAMpan')) %>%
  ggplot(aes(x = terms, y = abs(R2),fill=terms))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=paste0(ifelse(Pr..F.<.5,'*',''),round(R2*100,2))),stat='identity',size=2.5,color='white')+
  scale_fill_brewer(palette = "Greys",direction = -1)+
  facet_grid(pan_name~.)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison') + theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev)+fdb_style() +
  ggplot2::theme(aspect.ratio = 1,
                 axis.text.x=element_text(angle=45, hjust=1))+
  theme( panel.background = element_rect(colour = "black", size=1))+
  scale_x_discrete(limits=rev)


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/PERMANOVA_pans_sift_cleaner_blues3.pdf',plot=p.perm, width = 8, height = 8,unit='in')




#-------------------------------------------------------------------------------#
# run dbRDA
dbRDA <- vegan::capscale(distances.tibble$pan[[3]] ~ group + Completeness + Contamination + sourceClass + Genome_size, data = distances.tibble$metadata[[3]],distance='manhattan')
plot(dbRDA)



dbRDA <- vegan::capscale(distances.tibble$pan[[3]] ~ group + Completeness + Contamination + sourceClass + Genome_size + GC, data = distances.tibble$metadata[[3]],distance='manhattan')
plot(dbRDA)


#dbRDA <- vegan::capscale(distances.tibble$pan[[3]] ~ group , data = distances.tibble$metadata[[3]],distance='manhattan')
#plot(dbRDA)

#dbRDA <- vegan::capscale(distances.tibble$pan[[3]] ~ sourceClass , data = distances.tibble$metadata[[3]],distance='manhattan')
#plot(dbRDA)


#significance of overall model
anova(dbRDA)

#test signifciance of axes
anova(dbRDA, by ='axis',perm.max=500)

# tes of terms
anova(dbRDA, by ='terms',permu=500)

#-------------------------------------------------------------------------------#

# calculate betadisp
ls.betasdisp
for(i in 1:nrow(out.dist.tibble)){
  message(i)
}

betadisper.obj <- vegan::betadisper(as.dist(out.dist.tibble$distance[[3]]),
                             group=out.dist.tibble$metadata[[3]]$group)




vegan::adonis2(pan ~ as.formula("group + sourceClass + Completeness + Contamination + Genome_size"), data = metadata,permutation=prm)

vegan::adonis2(as.dist(out.dist.tibble$distance[[3]]) ~ as.formula("group + sourceClass + Completeness + Contamination + Genome_size"), data = out.dist.tibble$metadata[[3]],permutation=prm)


betadisper(as.dist(out.dist.tibble$distance[[i]]),group=out.dist.tibble$metadata[[i]]$group)




metadata = cbind(metadata, mergedEnv[metadata$genome,c('SS','sourceClass','lon','lat')])


#----------------------------------------------------------------------------------------------------------------#


metadata_prep <- function(genome.tbl, sorted.cliques, selcol,remove_na=TRUE){
  mtd <- genome.tbl %>% left_join(sorted.cliques, by='genome')  %>% data.frame()
  mtd$grp <- mtd[,selcol]
  mtd <- mtd %>% mutate(group=paste0('cl',grp))
  mtd[mtd$group=='clNA','group']<-NA
  if(remove_na==TRUE){
    mtd <-  mtd %>% dplyr::filter(!is.na(group))
    }
  return(mtd)
  }



# Here we want to filter out those that do not have a SourceClass
metadata <- metadata_prep(genome.tbl, sorted.cliques=sorted.ANI.cliques, selcol='sc0.85')



metadata <- metadata %>% filter(!is.na(sourceClass))
metadata <- metadata %>% mutate(SS=as.character(SS)) %>% mutate(SS=ifelse(SS=='','unknown',SS))

#metadata <- metadata_prep(genome.tbl, sorted.cliques=sorted.ANI.cliques, selcol='sc0.85', remove_na=FALSE)



#### FFS THIS REALLY ANOYINGLY FUCKS UP LOOK BETTER AS THIS USED TO WORK
#
# panList = list('OGpan'=OG.pan[rownames(OG.pan),rownames(OG.pan)],
#                'COGpan'=COG.pan[rownames(OG.pan),rownames(OG.pan)],
#                'KOpan'=kfm.pan[rownames(OG.pan),rownames(OG.pan)],
#                'PFAMpan'=pfam.pan[rownames(OG.pan),rownames(OG.pan)],
#                'CAZYpan'=cazy.pan[rownames(OG.pan),rownames(OG.pan)],
#                'TCDBpan' = tcdb.pan[rownames(OG.pan),rownames(OG.pan)])

metadata <- metadata %>% filter(!is.na(sourceClass)) %>% filter(!is.na(phylogroup))

adonisdf=c()
adonisdfAIC =c()

#sto re distances
panDistances <- list()

for(i in names(panList)){
  message(i)
  pan = panList[[i]]
  metadf = metadata %>% data.frame()
  rownames(metadf) <- metadata$genome
  PRM_data = data.frame(pan[metadata$genome,])
  dPRM_data = micropan::distManhattan(PRM_data)
  panDistances[[i]] = dPRM_data
  #AICc_PERMANOVA(adonis_phylogeny)
  mod1 = vegan::adonis(dPRM_data ~ sourceClass, data = metadf[rownames(PRM_data),],permutation=10000)
  #mod2 = vegan::adonis(dPRM_data ~ group + sourceClass, data = metadata[rownames(PRM_data),],permutation=10000)
  adonisdf = rbind(adonisdf,
        cbind('terms'=rownames(data.frame(mod1$ aov.tab))[1:nrow(data.frame(mod1$ aov.tab))-1],data.frame(mod1$ aov.tab)[1:nrow(data.frame(mod1$ 	aov.tab))-1,]) %>% mutate(pan_name = i, model='mod1')
    )
  #adonisdfAIC = rbind(adonisdfAIC,
  #      unlist(c('model'='mod1',AICc_PERMANOVA(mod1),'pan_name'=i)))
  }


#adonisdfAIC = data.frame(adonisdfAIC)
#adonisdfAIC$AIC = as.numeric(as.character(adonisdfAIC$AIC))






# USING ADONIS2
#========================================#

adonisdf=c()
adonisdfAIC =c()

#sto re distances
panDistances <- list()

for(i in names(panList)){
  message(i)
  pan = panList[[i]]
  metadf = metadata %>% data.frame()
  rownames(metadf) <- metadata$genome
  PRM_data = data.frame(pan[metadata$genome,])
  dPRM_data = micropan::distManhattan(PRM_data)
  panDistances[[i]] = dPRM_data
  mod1 = vegan::adonis2(dPRM_data ~ group + Genome_size + Contamination + sourceClass , data = metadf[rownames(PRM_data),],permutation=10000)
  adonisdf = rbind(adonisdf,
                   cbind('terms'=rownames(data.frame(mod1)[1:(nrow(mod1)-1),]),
                         data.frame(mod1)[1:(nrow(mod1)-1),]) %>%
                     mutate(pan_name = i, model='mod1')
  )
}




p1 <- ggplot(adonisdf, aes(x = pan_name, y = R2,fill=terms))+
  geom_bar(stat="identity",position="fill")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=round(R2*100,2)),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~ model,scales='free',nrow=1)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+ggtitle(label='PERMANVA (adonis) model comparison')

# p2 <- ggplot(adonisdfAIC, aes(x =  pan_name, y = AIC,fill= k))+
#   geom_point()+
#   facet_wrap(~ model,scales='free',nrow=1)+
#   geom_bar(stat="identity")+
#   scale_y_continuous(expand = c(0, 0))+
#   theme_minimal()+
#   theme( panel.background = element_rect(colour = "black", size=0.5))+
#   xlab(label='model')+ylab(label='AIC')+
#   coord_flip()+
#   ggtitle(label='PERMANVA (adonis) model comparison')
#
# p3 <- ggplot(adonisdfAIC, aes(x =  k, y = AIC))+
#   geom_point(shape=1,size=0.5) +
#   geom_smooth(color='red',size=0.3,alpha=0.3,method='lm')+
#   facet_wrap(~ model,scales='free',nrow=1)+
#   theme_minimal()+theme( panel.background = element_rect(colour = "black", size=0.5))+ggtitle(label='AIC (adonis) model comparison - k vs AIC')




#=====
#experiment with tidyr
genome.phylo.nested <- metadata %>% group_by(phylogroup) %>% tidyr::nest()
genome.phylo.nested[[2]][[1]]
genome.phylo.nested$data[[1]]
df <-genome.phylo.nested$data[[1]]
#make a model

gss_model <- function(df){
  lm(Genome_size ~ predicted_genes, data=df)
}


#add models to the tibble
models <- bootstraps %>%
  mutate(model=purrr::map(metadata, gss_model))

models <- models %>% mutate(tidy=map(model, broom::tidy),
                            rsq = glance %>% purrr::map_dbl('r.squared'),
                            glance=map(model,broom::glance),
                            augment=map(model,broom::augment))

models  %>% ggplot(aes(rsq,reorder(bootstrap,rsq))) +
  geom_point()


models %>% unnest(tidy)
models %>% unnest(glance)
models %>% unnest(augment)

models %>% unnest(augment) %>% ggplot(aes(predicted_genes, .resid)) +
  geom_line(aes(group=bootstrap),alpha=1/3)+geom_hline(yintercept=0)+
  facet_wrap(~bootstrap)

models %>% unnest(tidy) %>% select(phylogroup, term, estimate, rsq) %>%
  spread(term,estimate) %>%
  ggplot2::ggplot(aes(`(Intercept)`,predicted_genes)) +
  geom_point(aes(size=rsq, color=phylogroup)) +
  geom_smooth(se=FALSE)+
  xlab('genome size') +
  ylab('predicted_genes')+scale_color_manual(values=phylogroup.colors)


#=== phylogroup stuff
genome.phylo.nested <- metadata %>% group_by(phylogroup) %>% tidyr::nest()
genome.phylo.nested[[2]][[1]]
genome.phylo.nested$data[[1]]
df <-genome.phylo.nested$data[[1]]

models <- genome.phylo.nested %>%
  mutate(model=purrr::map(metadata, gss_model))

models %>% unnest(tidy) %>% select(phylogroup, term, estimate, rsq) %>%
  spread(term,estimate) %>%
  ggplot2::ggplot(aes(`(Intercept)`,predicted_genes)) +
  geom_point(aes(size=rsq, color=phylogroup)) +
  geom_smooth(se=FALSE)+
  xlab('genome size') +
  ylab('predicted_genes')+scale_color_manual(values=phylogroup.colors)


models %>% unnest(augment) %>% ggplot(aes(predicted_genes, .resid)) +
  geom_line(aes(group=phylogroup),alpha=1/3)+geom_hline(yintercept=0)+
  facet_wrap(~phylogroup)


#=====================================

#===========
# THE BOOTSTRAPPING APPRAOCH
# Uses the resticted dataset, including the assigned habiat ones
# this will analyse only those with both group and sourceclass

#we use bootstrap_group function to extract 25 random subsets, whereby each 95% group conly contains 1 representative at the time
set.seed(1234567890)
bootstrap.tips_2 <- bootstrap_group(metadata, sample.col='sc0.95',n=25)
bootstrap.tips_2 <- bootstrap_group(genome.tbl, sample.col='sc0.95',n=25)

bootstrap.tips_2 <- bootstrap_group(metadata %>% filter(genome %in% sel.genome),
                                    sample.col='sc0.95',n=25)



lengths(bootstrap.tips_2)

#adonis permutations
prm <- 10000




#===============================
# BOOTSTRAP IN TIDYR FORMAT
# models <- genome.phylo.nested %>%
# mutate(mod=purrr::map(data, adonis_model))


#### test


#BOOTSTRAP level tibble containing filtered pantables, distances, and metadata

#Select pan table
pan = panList[[1]]

panz <- list()
metadatz <- list()

for(i in 1:length(bootstrap.tips_2)){
  PRM_data = data.frame(pan[bootstrap.tips_2[[i]],])
  PRM_data <- pan_clean(PRM_data, remove_singletons = TRUE, SD=2)
  panz[[i]] <- PRM_data
  mt <- metadata %>% filter(genome %in% bootstrap.tips_2[[i]])
  rownames(mt) <- mt$genome
  metadatz[[i]] <- mt[bootstrap.tips_2[[i]],]
  }

bootstraps <- tibble(bootstrap=c(1:25),
                     genomes.set = bootstrap.tips_2,
                     pan = panz,
                     metadata = metadatz)


# calculate pairwise distances
ls_dstmn <- list()
for(i in 1:nrow(bootstraps)){
  message('analysing bootsrtap ', i)
  ls_dstmn[[i]] = pan2dist(bootstraps$pan[[i]], dist.method="Manhattan")
}

#add distances to the tibble
bootstraps <- bootstraps %>% mutate(dist = ls_dstmn)

# set up an adonis model
adonis_model <- function(pan,metadata){
  vegan::adonis2(pan ~ group + Completeness + Contamination + sourceClass + Genome_size, data = metadata,permutation=prm)
}

#adonis_model <- function(pan,metadata){
#  vegan::adonis2(pan ~ group + Completeness + Contamination + SS + Genome_size, data = metadata,permutation=prm)
#}


#adonis_model <- function(pan,metadata){
#  vegan::adonis2(pan ~ group + Completeness + Contamination + Genome_size, data = metadata,permutation=prm)
#}


#run over all bootstraps and calculate the adonis models
ls_adonis <- list()
for(i in 1:nrow(bootstraps)){
  message('analysing bootsrtap ', i)
  ls_adonis[[i]]<-adonis_model(bootstraps$dist[[i]], bootstraps$metadata[[i]])
  #mtdt <- bootstraps$metadata[[i]]
  #mtdt[ is.na(mtdt$SS),'SS']<-'unknown'
  #ls_adonis[[i]]<-adonis_model(bootstraps$dist[[i]], mtdt)
  }

#add models to the tibble
bootstraps <- bootstraps %>% mutate(adonis = ls_adonis)

#extract the results (no tidy version out yet)
ls_adonis.res <- list()
for(i in 1:nrow(bootstraps)){
  ls_adonis.res[[i]] = cbind('terms'=rownames(data.frame(bootstraps$adonis[[i]]))[1:nrow(data.frame(bootstraps$adonis[[i]]))-1],data.frame(bootstraps$adonis[[i]])[1:nrow(data.frame(bootstraps$adonis[[i]]))-1,]) %>% mutate(pan_name = i,model='mod1')
  }

#add models to the tibble
bootstraps <- bootstraps %>% mutate(res.adonis = ls_adonis.res)






p1 <- bootstraps %>% tidyr::unnest(res.adonis) %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residuals','Contamination','Completeness','Genome_size','SS','group'))) %>%
  ggplot(aes(x = pan_name, y = R2,fill=terms))+
    geom_bar(stat="identity",position="fill")+
    scale_y_continuous(expand = c(0, 0))+
    theme_bw()+
    theme( panel.background = element_rect(colour = "black", size=0.5))+
    geom_text(aes(label=round(R2*100,2)),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
    scale_fill_brewer(palette = "BuPu")+
    facet_wrap(~ model,scales='free',nrow=1)+
    xlab(label='model')+ylab(label='variance explained (R2)')+
    coord_flip()+
    ggtitle(label='PERMANOVA (adonis) model comparison')
p1
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/bootstrappedPermanova2.pdf',plot=p1, width = 5, height = 3,unit='in')


p2 <- bootstraps %>% tidyr::unnest(res.adonis) %>%
  dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))) %>%
  ggplot(aes(x = terms, y = R2))+
  geom_jitter(shape=21,fill='grey',color='darkgrey', width = 0.25)+
  facet_wrap(~ model,scales='free',nrow=1)+
  xlab(label='terms')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison')+fdb_style()+
  geom_crossbar(data=aggregate(R2 ~ terms, median, data=bootstraps %>% tidyr::unnest(res.adonis) %>%
                                 dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group')))),
                aes(ymin = R2, ymax = R2),
                size=0.5,col="black", width = .75)+
  geom_text(data=aggregate(R2 ~ terms, median, data=bootstraps %>% tidyr::unnest(res.adonis)) %>%
              dplyr::mutate(terms = factor(terms, levels=c('Residual','Contamination','Completeness','Genome_size','sourceClass','group'))),aes(label=round(R2,2)),hjust=-0.1)
p2




ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/bootstrappedPermanova2_g2.pdf',plot=p2, width = 4, height = 3,unit='in')


#==========




#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------







#-------------
#multiple models


adonisdf=c()
adonisdfAIC =c()

for(i in names(panList)){
  message(i)
  pan = panList[[i]]
  PRM_data = data.frame(pan[intersect(metadata$genome,rownames(pan)),])
  dPRM_data = micropan::distManhattan(PRM_data)
  metadf = metadata %>% data.frame()
  rownames(metadf) <- metadata$genome

  #AICc_PERMANOVA(adonis_phylogeny)
  mod1 = vegan::adonis(dPRM_data ~ group + Completeness + Contamination + sourceClass + Genome_size, data = metadf[rownames(PRM_data),],permutation=10000)
  #mod1 = vegan::adonis(dPRM_data ~ group + sourceClass + Genome_size, data = metadata[rownames(PRM_data),],permutation=10000)
  #mod2 = vegan::adonis(dPRM_data ~ group + sourceClass, data = metadata[rownames(PRM_data),],permutation=10000)

    adonisdf = rbind(adonisdf,
                   cbind('terms'=rownames(data.frame(mod1$ aov.tab))[1:nrow(data.frame(mod1$ aov.tab))-1],data.frame(mod1$ aov.tab)[1:nrow(data.frame(mod1$ 	aov.tab))-1,]) %>% mutate(pan_name = i,model='mod1')
  )
  adonisdfAIC = rbind(adonisdfAIC,
                      unlist(c('model'='mod1',AICc_PERMANOVA(mod1),'pan_name'=i)))
}

adonisdfAIC = data.frame(adonisdfAIC)

adonisdfAIC$AIC = as.numeric(as.character(adonisdfAIC$AIC))

p1 <- ggplot(adonisdf, aes(x = pan_name, y = R2,fill=terms))+
  geom_bar(stat="identity",position="fill")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  geom_text(aes(label=round(R2*100,2)),stat='identity',position=position_fill(vjust=0.5),size=2.5,color='white')+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~ model,scales='free',nrow=1)+
  xlab(label='model')+ylab(label='variance explained (R2)')+
  coord_flip()+
  ggtitle(label='PERMANOVA (adonis) model comparison')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_adonis_ANIc0.85.pdf',plot=p1, width = 8, height = 4,unit='in')

p2 <- ggplot(adonisdfAIC, aes(x =  pan_name, y = AIC,fill= k))+
  geom_point()+
  facet_wrap(~ model,scales='free',nrow=1)+
  geom_bar(stat="identity")+
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal()+
  theme( panel.background = element_rect(colour = "black", size=0.5))+
  xlab(label='model')+ylab(label='AIC')+
  coord_flip()+
  ggtitle(label='PERMANOVA model comparison')

p3 <- ggplot(adonisdfAIC, aes(x =  k, y = AIC))+
  geom_point(shape=1,size=0.5) +
  geom_smooth(color='red',size=0.3,alpha=0.3,method='lm')+
  facet_wrap(~ model,scales='free',nrow=1)+
  theme_minimal()+theme( panel.background = element_rect(colour = "black", size=0.5))+ggtitle(label='AIC (adonis) model comparison - k vs AIC')

#---------------


p4 <- adonisdf %>%
  filter(terms != 'Residuals') %>%
  ggplot(aes(terms,pan_name))+
    geom_point(aes(color=-log10(Pr..F.)),size=4)+
    ylab('')+xlab('') +scale_x_discrete(position = "top") +fdb_style()+
      theme(axis.text.x = element_text(angle = 90,hjust=0))+
  viridis::scale_color_viridis(option='E',direction=-1)

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/pan_adonis_ANIc0.85_pval.pdf',plot=p4, width = 3, height = 3,unit='in')

p5 <-adonisdf %>%
  filter(terms != 'Residuals') %>%
  ggplot(aes(F.Model,pan_name))+
  geom_point(aes(color=terms),size=4)


gridExtra::grid.arrange(p4,p1,nrow=1)
#------

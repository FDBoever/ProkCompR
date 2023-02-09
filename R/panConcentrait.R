# FUNCTION FROM ANNOTREE
#' @title Homoplasy slope ratio
#'
#' @description
#' HSR returns the homoplasy slope ratio of binary character states on a
#' phylogenetic tree
#'
#' @details
#' The homoplasy slope (HS) is calculated as follows: HS = ES/(t-3) where t = #
#' taxa and ES = extra steps where ES = 1/CI-1 where CI = consistency index.
#' This can be calculated using phangorn::CI. Then HSR = HS(real)/HS(random)
#' where HS(random) is the average HS of nsims simulated data sets of size t.
#' The probability of either state is equal for all taxa.
#'
#' @param tree  Phylo object
#' @param data  Named vector of states. Names correspond with tip labels of the
#'              tree
#' @param nsims Number of simulations to run for the estimation of a random
#'              homoplasy slope (default: 100)
#' @return Homoplasy slope ratio of the tip state data and tree
#' @export
#' @author Kerrin Mendler <mendlerke@gmail.com>
#' @references Meier R, Kores P, & Darwin S. 1991. Homoplasy slope ratio: a
#'             better measurement of observed homoplasy in cladistic analyses.
#'             Sys. Zool. 40(1):74-88.

HSR <- function(tree, data, nsims=100) {
  #require(phangorn)

  num_taxa_t <- length(data)
  tip_levels <- as.character(unique(data))

  # Homoplasy slope of given data
  data.phyDat <- phangorn::as.phyDat(sapply(data, as.character), type = "USER", levels = tip_levels)
  es_real <-  1/phangorn::CI(tree, data.phyDat) - 1
  hs_real <- es_real/(num_taxa_t - 3)

  # Average homoplasy slope of random data with same # taxa
  random_data <- replicate(nsims,
                           setNames(sample(tip_levels,
                                           num_taxa_t,
                                           replace = TRUE),
                                    tree$tip.label),
                           simplify = FALSE)
  random_data.phyDat <- lapply(random_data,
                               phangorn::as.phyDat, type = "USER", levels = tip_levels)
  es_random <- unlist(lapply(random_data.phyDat, function(tip_states){
    1/phangorn::CI(tree, tip_states) - 1
  }))
  hs_random <- es_random/(num_taxa_t - 3)
  avg_hs_random <- mean(hs_random)

  return(hs_real/avg_hs_random)
}


#pan_HSR
#=================================================#
#wrapper to run HSR on all columns of a pangenome
#' Title
#'
#' @param pan
#' @param tree
#' @param cutoff
#' @param nsims
#'
#' @return
#' @export
#'
#' @examples
#' pan_hsr <- pan_HSR(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75, nsims=100)

pan_HSR <- function(pan, tree, cutoff, nsims=100){
  pan_hsr <- NULL
  for(mod in colnames(pan)){
    message('analysing', mod)
    states <- pan[tip_order(tree),mod]
    names(states) <- rownames(pan[tip_order(tree),])
    states[states<cutoff] <- 0
    states[states>=cutoff] <- 1
    if(sum(states>0)){
      pan_hsr<-rbind(pan_hsr, c('feature'=mod,
                                'hsr'= HSR(tree = lsTrees.94.rooted[[1]],data = states ,nsims = 100)))
    }
  }

  pan_hsr <- pan_hsr %>% data.frame() %>%
    mutate(hsr = as.numeric(as.character(hsr)))

  return(pan_hsr)
}



pan_hsr <- pan_HSR(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75, nsims=100)
pan_hsr <- pan_HSR(pan=OG.pan, tree=lsTrees.94.rooted[[1]], cutoff = 1, nsims=100)


pan_HSR(pan=OG.pan[,1:5], tree=lsTrees.94.rooted[[1]], cutoff = 1, nsims=100)


# trait_patchiness
#=================================================#
#' Title
#' Pathciness is calculated as follows: -ln(CI)/ln(family size)
#' see, Annotree publication (https://academic.oup.com/nar/article/47/9/4442/5432638)
#' @param tree phylogenetic tree object
#' @param data named vector of trait state data (presence/absence as 1/0), names correspond to tree tip labels
#'
#' @return vector with CI, Patchiness and family_size (number of positive tips)
#' @export
#'
#' @examples
#'
trait_patchiness <- function(tree, data) {
  out <- NULL
  #require(phangorn)
  num_taxa_t <- length(data)
  tip_levels <- as.character(unique(data))
  data.phyDat <- phangorn::as.phyDat(sapply(data, as.character), type = "USER", levels = tip_levels)
  ci <- phangorn::CI(tree, data.phyDat)
  #ln in R through log
  patchiness <- -log(ci)/log(length(data[data>0]))
  out <- c('CI'=ci,'patchiness'=patchiness,'family_size'=length(data[data>0]))
  return(out)
}


# pan_patchiness
#=================================================#
#' Title
#'wrapper to calculate patchiness for each column in a feature abundance table
#' @param pan pan table/ feature abundance table
#' @param tree phylogenetic tree object
#' @param cutoff cutoff to assign presence or absence (for example above 1 for genes, or 90% for module completion ratios)
#'
#' @return table with outputs for each feature in the feature abundance table
#' @export
#'
#' @examples
#' pan.patchiness <- pan_patchiness(pan=m.pan, tree=lsTrees.94.rooted[[1]],cutoff=0.75)

pan_patchiness <- function(pan, tree, cutoff){
  out <- NULL
  for(mod in colnames(pan)){
    message('analysing', mod)
    states <- pan[tip_order(tree),mod]
    names(states) <- rownames(pan[tip_order(tree),])
    states[states<cutoff] <- 0
    states[states>=cutoff] <- 1
    if(sum(states>0)){
      out<-rbind(out, c('feature'=mod,
          trait_patchiness(tree = lsTrees.94.rooted[[1]],data = states)))
    }
  }

  out <- out %>% data.frame() %>%
    mutate(ci = as.numeric(as.character(ci))) %>%
    mutate(patchiness = as.numeric(as.character(patchiness))) %>%
    mutate(family_size = as.numeric(as.character(family_size)))

  return(out)
}




#=======================#
#### MAPLE ####

pan.patchiness <- pan_patchiness(pan=m.pan, tree=lsTrees.94.rooted[[1]],cutoff=0.75)
pan.patchiness %>% filter(!is.na(patchiness)) %>% ggplot(aes(patchiness,reorder(feature,patchiness)))+geom_point()

sel.features <- pan.patchiness %>%
  filter(!is.na(patchiness)) %>%
  filter(patchiness > 0) %>% mutate(feature=as.character(feature))%>%
  arrange(patchiness) %>% select(feature) %>% pull() %>% as.character


sel.features <- pan.patchiness %>%
  filter(!is.na(patchiness)) %>%
  filter(patchiness == 0) %>% mutate(feature=as.character(feature))%>%
  arrange(patchiness) %>% select(feature) %>% pull() %>% as.character



trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_Pathciness','.pdf'),plot=combined.plots, width = 30, height = 30,unit='in')



consentrait_hsr_patchiness <- consentrait.out %>%
  #left_join(pan_hsr, by='feature') %>%
  left_join(pan.patchiness, by='feature')


p.patch_meandepth <- consentrait_hsr_patchiness %>%
  filter(mean_depth>0) %>%
  filter(!is.na(mean_depth)) %>%
  filter(patchiness>0) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  ggplot(aes(mean_depth,patchiness))+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F, color='grey')+
  geom_point(shape=21,alpha=0.5,aes(fill=TYPE,size=family_size)) +
  fdb_style()+
  ggrepel::geom_text_repel(aes(label=feature),max.overlaps = 50)+
  scale_size_continuous(range = c(1, 4))+ylab('Patchiness')+xlab('mean phylogenetic depth')


ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_patchiness_v_meandepth','.pdf'),plot=p.patch_meandepth, width = 6, height = 6,unit='in')







consentrait_hsr_patchiness %>%
  ggplot(aes(mean_depth,family_size))+
  geom_point(shape=21)+
  geom_smooth(method='lm') + fdb_style()


sel.func.c <- 'monooxygenase'
sel.func.c <- 'alkane'

sel.ogs <- annotation.tbl %>% filter(grepl(sel.func.c,product,ignore.case=TRUE)) %>% filter(!is.na(OG)) %>% select(OG) %>%
  unique() %>% pull()

conversion_tb <- annotation.tbl %>% filter(grepl(sel.func.c,product,ignore.case=TRUE)) %>% filter(!is.na(OG)) %>% select(OG,product) %>% unique() %>% group_by(OG) %>% summarise(product=paste(unique(product),collapse='; ')) %>% mutate(product=paste(OG,product,sep=' - '))

sel.matrix <- OG.pan[,sel.ogs]
colnames(sel.matrix) <- conversion_tb %>% select(product) %>% pull()

gplots::heatmap.2(as.matrix(t(sel.matrix)),trace='none' ,margins = c(5,20))




#pan_consenTRAIT
#=================================================
#' Title
#'Function using castor's implementation of consentrait
#'This function will run it on an entire feature abundance table
#' Since castor 1.7, we switched to consentrait_depht as get_trait_depth depricated
#' https://github.com/cran/castor/blob/master/R/consentrait_depth.R
#' @param pan pan table/ feature abundance table
#' @param tree phylogenetic tree object
#' @param cutoff cutoff to assign presence or absence (for example above 1 for genes, or 90% for module completion ratios)
#'
#' @return table with outputs for each feature in the feature abundance table
#' @export
#'
#' @examples
pan_consenTRAIT <- function(pan,tree, cutoff){
  consentrait.out <- NULL
  castorlist <- list()
  for(mod in colnames(pan)){
    states <- pan[tip_order(tree),mod]
    names(states) <- rownames(pan[tip_order(tree),])
    states[states<cutoff] <- 0
    states[states>=cutoff] <- 1
    if(sum(states>0)){
      castorlist[[mod]]<-castor::consentrait_depth(tree=tree,
                                                 tip_states=states,
                                                 Npermutations = 1000,
                                                 #count_singletons = F,
                                                 min_fraction = 0.90)
      consentrait.out <- rbind(consentrait.out,
                               c('feature'=mod,
                                 'total_tips'=sum(states),
                                 'mean_depth'=castorlist[[mod]]$mean_depth,
                                 'var_depth'=castorlist[[mod]]$var_depth,
                                 'min_depth'=castorlist[[mod]]$min_depth,
                                 'max_depth'=castorlist[[mod]]$max_depth,
                                 'Npositives'=castorlist[[mod]]$Npositives,
                                 'P'=castorlist[[mod]]$P,
                                 'mean_random_depth'=castorlist[[mod]]$mean_random_depth))
      message(mod,castorlist[[mod]]$mean_depth)
    }

  }


  consentrait.out <- consentrait.out %>% data.frame() %>%
    mutate(total_tips = as.numeric(as.character(total_tips)),
           mean_depth = as.numeric(as.character(mean_depth)),
           var_depth = as.numeric(as.character(var_depth)),
           min_depth = as.numeric(as.character(min_depth)),
           max_depth = as.numeric(as.character(max_depth)),
           Npositives = as.numeric(as.character(Npositives)),
           P = as.numeric(as.character(P)),
           mean_random_depth = as.numeric(as.character(mean_random_depth))
    )

  return(consentrait.out)
}


consentrait.out <- pan_consenTRAIT(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75)
consentrait.out <- pan_consenTRAIT(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75)

consentrait.out.OG <- pan_consenTRAIT(pan=pan_clean(OG.pan[sel.genome,]), tree=lsTrees.94.rooted[[1]], cutoff = 1)
consentrait.out.OG <- pan_consenTRAIT(pan=pan_clean(OG.pan[sel.genome,]), tree=lsTrees.94.rooted[[1]], cutoff = 1)


patchiness.out.OG <- pan_patchiness(pan=pan_clean(OG.pan[sel.genome,]), tree=lsTrees.94.rooted[[1]],cutoff=0.75)



patchiness.out.OG %>% ggplot(aes(patchiness,family_size)) + geom_point()


#if
consentrait.out.OG %>%
  filter(P<0.05) %>%
  left_join(patchiness.out.OG ,by=c('feature'='feature')) %>%
  ggplot(aes(mean_depth,patchiness))+geom_point()+  scale_x_log10()+
  coord_trans(x='reverse')+
  xlab(expression("Mean trait depth ("~Tau[D]~")"))


ANI.hmg.grp.tbl %>%
  filter(!(feature.id %in% CORE.ogs))%>%
  mutate(sign_depleted = ifelse(fdr.p.value_depletion_binary<0.05 |
                                  fdr.value_depletion_rawcounts < 0.05,
                                'lineage_associated','not_associated'),
         sign_enriched = ifelse(fdr.p.value_enriched_binary<0.05 |
                                  fdr.p.value_enriched_rawcounts < 0.05,
                                'lineage_associated','not_associated')) %>%
  tidyr::pivot_longer(cols=c(sign_depleted,sign_enriched),names_to='var',values_to='value') %>%
  left_join(consentrait.out.OG ,by=c('feature.id'='feature')) %>%
  left_join(patchiness.out.OG ,by=c('feature.id'='feature')) %>%
  #mutate(sign = ifelse(P<0.05, 'sign','non-sign')) %>%
  #filter(sign=='sign') %>%
  ggplot(aes(patchiness,var,fill=value)) +
  geom_boxplot() +
  #geom_jitter()  +
  scale_x_log10()+
  coord_trans(x='reverse')+
  xlab(expression("Mean trait depth ("~Tau[D]~")")) +
  scale_fill_manual(values=c('#22B2DA','#555A60')) + facet_wrap(~group, ncol=1)






ANI.hmg.grp.tbl %>%
  filter(!(feature.id %in% CORE.ogs))%>%
  mutate(sign_depleted = ifelse(fdr.p.value_depletion_binary<0.05 |
                                  fdr.value_depletion_rawcounts < 0.05,
                                'lineage_associated','not_associated'),
         sign_enriched = ifelse(fdr.p.value_enriched_binary<0.05 |
                                  fdr.p.value_enriched_rawcounts < 0.05,
                                'lineage_associated','not_associated')) %>%
  tidyr::pivot_longer(cols=c(sign_depleted,sign_enriched),names_to='var',values_to='value') %>%
  left_join(consentrait.out.OG ,by=c('feature.id'='feature')) %>%
  left_join(patchiness.out.OG ,by=c('feature.id'='feature')) %>%
  filter(value=='lineage_associated') %>%
  ggplot(aes(mean_depth,patchiness,color=var)) +
  geom_point(shape=21) +
    scale_x_log10()+
  coord_trans(x='reverse')+
  xlab(expression("Mean trait depth ("~Tau[D]~")")) +
  scale_color_manual(values=c('#22B2DA','#555A60')) + facet_wrap(~group)



p.p2 <- ANI.hmg.grp.tbl %>%
  filter(!(feature.id %in% CORE.ogs))%>%
  mutate(sign_depleted = ifelse(fdr.p.value_depletion_binary<0.05 |
                                  fdr.value_depletion_rawcounts < 0.05,
                                'lineage_associated','not_associated'),
         sign_enriched = ifelse(fdr.p.value_enriched_binary<0.05 |
                                  fdr.p.value_enriched_rawcounts < 0.05,
                                'lineage_associated','not_associated')) %>%
  tidyr::pivot_longer(cols=c(sign_depleted,sign_enriched),names_to='var',values_to='value') %>%
  left_join(consentrait.out.OG ,by=c('feature.id'='feature')) %>%
  left_join(patchiness.out.OG ,by=c('feature.id'='feature')) %>%
  filter(value=='lineage_associated') %>%
  ggplot(aes(mean_depth,group)) +
  geom_point(shape=21,position=position_jitterdodge(),alpha=0.5,aes(color=var)) +
  geom_boxplot(position= position_dodge(width = 0.9),outlier.colour = NA,aes(fill=var))+
  scale_x_log10()+
  coord_trans(x='reverse')+
  xlab(expression("Mean trait depth ("~Tau[D]~")")) +
  scale_color_manual(values=c('#22B2DA','#555A60')) +
  scale_fill_manual(values=c('#22B2DA','#555A60')) + fdb_style()#+ facet_wrap(~group)
#+ facet_wrap(~group)




p.p1 <- ANI.hmg.grp.tbl %>%
  filter(!(feature.id %in% CORE.ogs))%>%
  mutate(sign_depleted = ifelse(fdr.p.value_depletion_binary<0.05 |
                                  fdr.value_depletion_rawcounts < 0.05,
                                'lineage_associated','not_associated'),
         sign_enriched = ifelse(fdr.p.value_enriched_binary<0.05 |
                                  fdr.p.value_enriched_rawcounts < 0.05,
                                'lineage_associated','not_associated')) %>%
  tidyr::pivot_longer(cols=c(sign_depleted,sign_enriched),names_to='var',values_to='value') %>%
  left_join(consentrait.out.OG ,by=c('feature.id'='feature')) %>%
  left_join(patchiness.out.OG ,by=c('feature.id'='feature')) %>%
  filter(value=='lineage_associated') %>%
  ggplot(aes(patchiness,group)) +
  geom_point(shape=21,position=position_jitterdodge(),alpha=0.5,aes(color=var)) +
  geom_boxplot(position= position_dodge(width = 0.9),outlier.colour = NA,aes(fill=var))+
  scale_x_log10()+
  #coord_trans(x='reverse')+
  xlab(expression("Patchiness")) +
  scale_color_manual(values=c('#22B2DA','#555A60')) +
  scale_fill_manual(values=c('#22B2DA','#555A60')) + fdb_style()#+ facet_wrap(~group)


ggpubr::ggarrange(p.p1,p.p2)











ANI.hmg.grp.tbl %>%
  filter(!(feature.id %in% CORE.ogs))%>%
  mutate(sign_depleted = ifelse(fdr.p.value_depletion_binary<0.05 |
                                  fdr.value_depletion_rawcounts < 0.05,
                                'lineage_associated','not_associated'),
         sign_enriched = ifelse(fdr.p.value_enriched_binary<0.05 |
                                  fdr.p.value_enriched_rawcounts < 0.05,
                                'lineage_associated','not_associated')) %>%
  tidyr::pivot_longer(cols=c(sign_depleted,sign_enriched),names_to='var',values_to='value') %>%
  left_join(consentrait.out.OG ,by=c('feature.id'='feature')) %>%
  mutate(sign = ifelse(P<0.05, 'sign','non-sign')) %>%
  filter(sign=='sign') %>%
  ggplot(aes(mean_depth,var,fill=value)) +
  geom_boxplot() +
  #geom_jitter()  +
  scale_x_log10()+
  coord_trans(x='reverse')+
  scale_x_continuous(trans='log10') + xlab(expression("Mean trait depth ("~Tau[D]~")")) +
  scale_fill_manual(values=c('#22B2DA','#555A60')) + facet_wrap(~group, ncol=1)






consentrait.out.OG %>%
  filter(!(feature %in% CORE.ogs))%>%
  mutate(sign = ifelse(P<0.05, 'sign','non-sign')) %>%
  left_join(og.table ,by=c('feature'='feature.id')) %>%
  #filter(sign =='sign') %>%
  #left_join(cog_func, by=c('feature'='COG')) %>%
  ggplot(aes(mean_depth,CORE,fill=sign)) +
  geom_boxplot() +
  #geom_jitter()  +
  scale_x_log10()+
  coord_trans(x='reverse')+
  scale_x_continuous(trans='log10') + xlab(expression("Mean trait depth ("~Tau[D]~")")) +
  scale_fill_manual(values=c('#22B2DA','#555A60'))







consentrait.out <- pan_consenTRAIT(pan=COG.pan, tree=lsTrees.94.rooted[[1]], cutoff = 1)
pan.patchiness <- pan_patchiness(pan=COG.pan, tree=lsTrees.94.rooted[[1]],cutoff=1)

consentrait_hsr_patchiness <- consentrait.out %>%
  #left_join(pan_hsr, by='feature') %>%
  left_join(pan.patchiness, by='feature')


p.patch_meandepth <- consentrait_hsr_patchiness %>%
  filter(mean_depth>0) %>%
  filter(!is.na(mean_depth)) %>%
  filter(patchiness>0) %>%
  left_join(cog_func, by=c('feature'='COG')) %>%
  ggplot(aes(mean_depth,patchiness))+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F, color='grey')+
  geom_point(shape=21,alpha=0.5,aes(fill=func,size=family_size)) +
  fdb_style()+
  scale_size_continuous(range = c(1, 4))+ylab('Patchiness')+xlab('mean phylogenetic depth')

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','COG_pan_patchiness_v_meandepth','.pdf'),plot=p.patch_meandepth, width = 6, height = 6,unit='in')






pan_hsr



consentrait.out %>%  ggplot(aes(x=mean_depth))+geom_histogram(stat='count')

consentrait.out %>%  ggplot(aes(x=mean_depth,total_tips, label=feature))+
  geom_point(shape=21) +
  geom_text_repel() +
  geom_smooth() + fdb_style()


consentrait.out %>% filter(mean_depth > 0.5) %>% left_join(m2n_filtered, by=c('feature'='module'))

consentrait.out %>%
  filter(mean_depth > 0.0) %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>% arrange(desc(total_tips))


consentrait.out %>%
  filter(mean_depth > 0.0) %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>% arrange(desc(total_tips)) %>%
  ggplot(aes(mean_depth,reorder(DESCRIPTION,mean_depth))) + geom_point() +facet_wrap(~TYPE,scales='free_y')

consentrait.out %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='Pathway') %>%
  ggplot(aes(mean_depth,reorder(DESCRIPTION,mean_depth))) +
  geom_point() +
  geom_errorbarh(aes(xmin=mean_depth-var_depth, xmax=mean_depth+var_depth), height =0) +
  fdb_style(aspect.ratio =1.5) + ylab('') + xlab('Phylogenetic depth')




# phylo_traits_plot
#=================================================
#' Title
#'
#' @param pan
#' @param tree
#' @param cutoff
#' @param sf
#'
#' @return
#' @export
#'
#' @examples
phylo_traits_plot <- function(pan, tree, cutoff, sf, description=NULL){
  states <- pan[tip_order(tree),sf]
  names(states) <- rownames(pan[tip_order(tree),])
  states[states<cutoff] <- 0
  states[states>=cutoff] <- 1
  message(states)
  pos_genomes <- names(states[states==1])

  p.tree.annot <- ggtree::ggtree(tree,layout="circular")
  p.tree.annot$data <- p.tree.annot$data %>%
    dplyr::mutate(positive = label %in% pos_genomes)

  p.tree.annot <- p.tree.annot + ggtree::geom_tiplab(data=.%>%dplyr::filter(positive==TRUE),
                                             align = TRUE,
                                             size=0,
                                             linesize=1.25,
                                             linetype = 1,
                                             color='skyblue', offset=0.1)
  if(is.null(description)){
    p.tree.annot <- p.tree.annot +ggtitle(sf)
  }else{
    p.tree.annot <- p.tree.annot + labs(title = sf,
                                        subtitle = description[sf])
  }


  return(p.tree.annot)
}




#' Title
#'
#' @param pan
#' @param tree
#' @param cutoff
#' @param sel.features
#'
#' @return
#' @export
#'
#' @examples
phylo_traits_plots <- function(pan, tree, cutoff, sel.features, description=NULL){
  plotlist<-list()
  for(sf in sel.features){
    trait_plot <- phylo_traits_plot(pan=m.pan,
                                    tree=lsTrees.94.rooted[[1]],
                                    cutoff = 75,
                                    sf = sf,
                                    description=description)
    plotlist[[sf]] <- trait_plot
  }
  return(plotlist)
}


#===========================#

sel.features <- consentrait.out %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='Pathway') %>% select(feature) %>% pull


#single plot
trait_plot <- phylo_traits_plot(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sf = sel.features[2])
trait_plot


#multiple plots
trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features)
gridExtra::grid.arrange(grobs=trait_plots)



description <- as.character(m2n_filtered$DESCRIPTION)
names(description) <- as.character(m2n_filtered$module)

trait_plot <- phylo_traits_plot(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sf = sel.features[2],description=description)
trait_plot

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT','.pdf'),plot=combined.plots, width = 18, height = 14,unit='in')






sel.features <- consentrait.out %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='Pathway') %>% select(feature) %>% pull

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT_pathways','.pdf'),plot=combined.plots, width = 18, height = 14,unit='in')


sel.features <- consentrait.out %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='FuncSet') %>% select(feature) %>% pull

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT_FuncSet','.pdf'),plot=combined.plots, width = 18, height = 14,unit='in')


sel.features <- consentrait.out %>%
  filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='Complex') %>% select(feature) %>% pull

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT_Complex','.pdf'),plot=combined.plots, width = 18, height = 14,unit='in')

sel.features <- consentrait.out %>%
  #filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='Complex') %>% select(feature) %>% pull

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT_Complex_all','.pdf'),plot=combined.plots, width = 18, height = 14,unit='in')




sel.features <- consentrait.out %>%
  #filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='FuncSet') %>% select(feature) %>% pull

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT_FuncSet_all','.pdf'),plot=combined.plots, width = 18, height = 14,unit='in')




sel.features <- consentrait.out %>%
  #filter(mean_depth < 0.5) %>%
  left_join(m2n_filtered, by=c('feature'='module')) %>%
  arrange(desc(total_tips)) %>%
  filter(TYPE=='Pathway') %>% select(feature) %>% pull

trait_plots <- phylo_traits_plots(pan=m.pan, tree=lsTrees.94.rooted[[1]], cutoff = 75,sel.features = sel.features,description=description)
combined.plots <- gridExtra::grid.arrange(grobs=trait_plots)

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_set_of_consenTRAIT_pathways_all','.pdf'),plot=combined.plots, width = 30, height = 30,unit='in')










heatmap.colours <- c("white","grey","seagreen3","darkgreen",
                     "green","brown","tan","red","orange",
                     "pink","magenta","purple","blue","skyblue3",
                     "blue","skyblue2")


heatmap.d <- m.pan[tip_order(lsTrees.94.rooted[[1]]),] %>%
  tibble::rownames_to_column(var='id') %>% data.frame

p <- ggtree(lsTrees.94.rooted[[1]])

ggtree::gheatmap(p, heatmap.d[,1:10], offset = 10,
         colnames_position="top",
         colnames_angle=90, colnames_offset_y = 1,
         hjust=0, font.size=2) +
  viridis::scale_fill_viridis(discrete = TRUE)


pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>%
  gheatmap(m.pan[tip_order(lsTrees.94.rooted[[1]]),], offset=8, width=0.6, colnames=FALSE) %>%
  scale_x_ggtree()
pp + theme(legend.position="right")












#https://ps.uci.edu/scholar/amartiny/files/consentrait.txt
#The script requires two input files as arguments: a Newick Tree and tab delimited text file with names of each taxon in the first column and then 0 or 1 values for each trait in the following columns (and no headers).

library("data.table")
library("adephylo")
library("ape")

#args <- commandArgs(TRUE)
#loading Newick tree (multitree) - replace to read.nexus if using nexus tress
#tree_all = read.tree(args[1],keep.multi = TRUE)
#Loading trait table w. no headers
#table = read.table(args[2], sep = "\t", header=FALSE)



#========================
#
#
#

#presence absence, and repeat sco detection, should select CORE genes

# function to extract core features
core_features <- function(pan){
  pan.presabs<-pan
  pan.presabs[pan.presabs > 0] <- 1
  CORE = apply(t(pan.presabs %>% data.frame)==1,1,all)
  CORE = names(CORE[CORE ==TRUE])
  return(CORE)
}


#=============================================================
cazy.pan

head(cazy.pan)
metadata <- genome.tbl %>% filter(genome %in% sel.genome)

pan <- OG.pan

group <- 'phylogroup'
grps <- unique(metadata[,group]) %>% pull
grps <- grps[!is.na(grps)]
out <- NULL
for(grp in grps){
  message(grp)
  f.genomes <- metadata %>% filter(!!as.name(group) == grp) %>% select(genome) %>% pull
  f.pan <- pan[f.genomes,]
  f.pan <- f.pan[,colSums(f.pan)>0]
  message(ncol(f.pan))
  f.core <- core_features(f.pan)
  message(length(f.core))
  out <- rbind(out,c('group'=grp,
               'pan.size'=ncol(f.pan),
               'core.size'=length(f.core)
               ))
}
out <- data.frame(out) %>%
  mutate(pan.size=as.numeric(as.character(pan.size))) %>%
  mutate(core.size=as.numeric(as.character(core.size)))

out %>% ggplot(aes(pan.size,core.size))+
  geom_point()+geom_smooth(method='lm')+fdb_style()

out %>% ggplot(aes(group,core.size))+geom_col()
out %>% ggplot(aes(group,pan.size))+geom_col()


#diversifying
#=============================================================
pan <-cazy.pan
pan <-OG.pan

tree.list <- lsTrees[[1]]

#implement the multiroot thing
root_tree<-root(lsTrees[[1]],outgroup,resolve.root=T)
tree.list <- list('concat.sco.nuc'=root_tree)


tree.list <- list('concat.sco.nuc'=lsTrees.94.rooted[[1]])


#feautureConservation
#=============================================================
#' Title
#'
#' @param pan
#' @param tree.list
#' @param perc
#'
#' @return
#' @export
#'
#' @examples
feautureConservation <- function(pan, tree.list, perc){
  nodeMatch <- NULL
  pan[pan>1]<-1
  for(tr in 1:length(tree.list)){
    message(tr)
    tree.name <- names(tree.list)[tr]
    tree<-tree.list[[tr]]
    subtree<-subtrees(tree, wait=FALSE)
    cluster_mean<-numeric(length=0)

    for (i in 1:length(subtree)){
      tip_names<-subtree[[i]]$tip.label
      subpan <- pan[tip_names,]
      matched.pan <- subpan[,colSums(subpan)/nrow(subpan)>perc]
      n_matched <- ncol(matched.pan)
      nodeMatch <- rbind(nodeMatch, c('tree'=tree.name,'node'=i,'n'=n_matched))
    }
  }

  ntips <- length(lsTrees[[1]]$tip.label)
  out <- nodeMatch %>% data.frame(stringsAsFactors = FALSE) %>%
    mutate(node=as.numeric(node)) %>%
    mutate(ggnode=node+ntips) %>%
    mutate(n=as.numeric(n))
  return(out)
}




#================================================

nodeMatch <- feautureConservation(OG.pan[sel.genome,], lsTrees.94.rooted[1], perc=0.95)
nodeMatch <- feautureConservation(OG.pan[sel.genome,], lsTrees.94.rooted[1], perc=0.95)

tree = lsTrees.94.rooted[[1]]

internalNodes = data.frame(cbind('node'=1:tree$Nnode,
                                 'nrTips'= castor::count_tips_per_node(tree),
                                 'RED'=castor::get_reds(tree),
                                 'nodeDepths'=castor::get_all_node_depths(tree)),
                           stringsAsFactors=FALSE)


merged <- internalNodes %>% dplyr::left_join(nodeMatch,by='node')


ggplot(merged,aes(RED,n))+geom_point()+geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()
p.b<-ggplot(merged,aes(nodeDepths,n))+geom_point(shape=21,color='darkgrey',fill='grey',alpha=0.7)+geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F,color='forestgreen')+fdb_style()
p.a <- ggplot(merged,aes(nrTips,n))+geom_point(shape=21,color='darkgrey',fill='grey',alpha=0.7)+geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F,color='forestgreen')+fdb_style()
ggplot(merged,aes(nrTips,nodeDepths))+geom_point(shape=21,fill='grey',color='darkgrey',alpha=0.7)+geom_smooth()+fdb_style()

p.c <- gridExtra::grid.arrange(p.b,p.a,nrow=1)
ggsave(paste0('~/DATA/MarbGenomics/Graphs/core0.95_present','phylo','.pdf'),plot=p.c, width = 6, height = 3,unit='in')



merged %>%
  filter(nodeDepths<0.1)%>%
  ggplot(aes(nodeDepths,n))+geom_point(shape=21,fill='grey',color='darkgrey',alpha=0.7)+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F,linetype='dashed',color='red')+
  geom_smooth()+
  fdb_style()

merged %>%
  #filter(nodeDepths<0.2)%>%
  ggplot(aes(nodeDepths,n))+geom_point(shape=21,fill='grey',color='darkgrey',alpha=0.7)+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+
  geom_smooth()+
  fdb_style()


#----
p <- ggtree(root_tree)
p$data <- p$data %>% dplyr::left_join(nodeMatch%>% mutate(node=ggnode),by='node')
p + geom_point(aes(color=n,size=n))

p + geom_point(aes(color=n))

#============================================
phylogroup.nodes <- c(109,112,172,182,115,170,146,169,167,162,163)
phylogroup.names <- paste0('P',seq(1:length(phylogroup.nodes)))
d <- data.frame(node=phylogroup.nodes, type=phylogroup.names)
d

for(grp in unique(phylogroup.offspring$P1)){
  message(grp)
}

phygrps <- phylogroup.offspring$phylogroup %>% as.character() %>% unique()
tr <- tree
d_t <- NULL #ancestral node
p_t <- NULL #all descendent
for(phgr in phygrps){
  pt <- ggtree::ggtree(tr)
  tps <- genome.tbl %>%
    dplyr::filter(phylogroup==phgr) %>%
    dplyr::select(genome) %>%
    dplyr::pull() %>% as.character()

  if(ape::is.monophyletic(tr, tps)){
    message('+')
  }else{
    message('WARNING!! ',phgr, 'is not monophyletic in ',tr.name)
  }
  nd = ape::getMRCA( tr, tps)
  pt.off <- pt %>% tidytree::offspring(nd) %>% select(node) %>% mutate(phylogroup =phgr) %>% data.frame()
  p_t <- rbind(p_t,pt.off)
  d_t <- rbind(d_t,c('node'=nd,'type'=phgr))
}
d_t <- data.frame(d_t)
d_t$node <- as.numeric(as.character(d_t$node))
d_t$type <- as.character(d_t$type)


p.in1 <- merged %>% left_join(p_t, by='node') %>%
  #filter(nodeDepths<0.04)%>%
  filter(phylogroup %in% c('P1','P3','P4','P5','P7','P11'))%>%
  mutate(phylogroup=factor(phylogroup,levels= c('P1','P3','P4','P5','P7','P11'))) %>%
  ggplot(aes(nodeDepths,n))+geom_point(aes(color=phylogroup),size=2)+
  geom_smooth(aes(color=phylogroup),method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+
  geom_smooth(aes(color=phylogroup),method='lm',se=F,linetype='dashed')+
  facet_wrap(~phylogroup,scales='free',ncol=3)+
  scale_color_manual(values=phylogroup.colors)+
  scale_fill_manual(values=phylogroup.colors)+
  fdb_style()+
  theme(legend.position = "none")

ggsave(paste0('~/DATA/MarbGenomics/Graphs/core0.95_present_in_tips_','phylo','.pdf'),plot=p.in1, width = 10, height = 3,unit='in')


merged %>% left_join(p_t, by='node') %>%
  #filter(nodeDepths<0.04)%>%
  filter(phylogroup %in% c('P1','P3','P4','P5','P7','P11'))%>%
  mutate(phylogroup=factor(phylogroup,levels= c('P1','P3','P4','P5','P7','P11'))) %>%
  ggplot(aes(nodeDepths,n))+geom_point(aes(color=phylogroup),size=2)+
  geom_smooth(aes(color=phylogroup),method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+
  geom_smooth(aes(color=phylogroup),method='lm',se=F,linetype='dashed')+
  facet_wrap(~phylogroup,scales='free',ncol=3)+
  scale_color_manual(values=phylogroup.colors)+
  scale_fill_manual(values=phylogroup.colors)+
  fdb_style()+
  theme(legend.position = "none")

















merged %>% left_join(p_t, by='node') %>%
  #filter(nodeDepths<0.2)%>%
  #filter(phylogroup %in% c('P1','P3','P4','P5','P7'))%>%
  ggplot(aes(nodeDepths,n))+geom_point(aes(fill=phylogroup,color=phylogroup),shape=21,alpha=0.7)+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+
  #facet_wrap(~phylogroup)+
  fdb_style()

#




nodeMatch.count <- nodeMatch %>% filter(tree =='t') %>% dplyr::count(subtree)
nodeMatch.count <- nodeMatch.count %>% mutate(node = subtree + ntips)

nodeMatch.count <- nodeMatch.count %>% filter(subtree!=1)

root_tree<-root(lsTrees[[1]],outgroup,resolve.root=T)


p <- ggtree(root_tree)
p$data <- p$data %>% dplyr::left_join(nodeMatch.count,by='node')
p + geom_point(aes(color=n,size=n))

ggtree(subtree[[25]])+geom_tiplab(size=2)

castor::get_subtree_at_node(root_tree,1)

#=======


#================================
#
#
#===============================




#pan from metadata
pan <- mtdat %>% mutate(genome=as.factor(genome))%>%
  dplyr::group_by(genome) %>%
  dplyr::count(SS) %>%
  tidyr::spread(SS,n, fill=0) %>% data.frame()


pan <- pan %>% select(-genome) %>% dplyr::select(which(colSums(.)>5))

rownames(pan)=mtdat$genome
table=pan


table <- cazy.pan
table <- OG.pan
table[table>1]<-1

#The boever way, consider removing core traits, as this is not that informative?

table <- table %>% rownames_to_column(var='genome')


stored_names <- colnames(table)
names(table) <- NULL
tree_all <- lsTrees[1]

maindir <- '~/DATA/MarbGenomics/consenTRAIT/'
subdir <- 'cazy.test/'
subdir <- 'OG.pan/'
subdir <- 'habitat'
out_dir <- paste0(maindir, subdir)
dir.create(maindir)#, showWarnings = FALSE)
dir.create(out_dir)#, showWarnings = FALSE)

#pimpt this so that it removes previous run!


#Starting script
Mean_all<-matrix(nrow=ncol(table)-1,ncol=100)

#FDB implemented nodeMatch, which is counting the matches per node in trees
nodeMatch <- NULL

for (m in 1:length(tree_all)) {#loop through all trees
  tree.name <- names(tree_all)[m]
  tree<-tree_all[[m]]
  #tree <- tree_all
  # testing if table and tree contain the same entries - else drop tips
  z<-subset(tree$tip.label,!(tree$tip.label %in% table[,1]))
  if (length(z) > 0) {
    drop.tip(tree,z)
  }

  #rooting tree with first taxon - change if different root
  #!!!!!!!!!!!!!!!!!!!!!
  #FDB changed this to use our own outgroup
  #root_tree<-root(tree,1,resolve.root=T)

  tree<-root(tree,outgroup,resolve.root=T)

  #replacing negative branch lengths - e.g., from PHYLIP

  tree$edge.length[tree$edge.length <= 0] =  0.00001
  subtree<-subtrees(tree, wait=FALSE)


  cluster_mean<-numeric(length=0)

  # loop through all traits
  for (j in 2:ncol(table)) {
    message(paste0(c("Analyzing",j,"")))
    #Loading trait table
    table_tmp<-table[,c(1,j)]
    colnames(table_tmp)[1]<-"ID";
    colnames(table_tmp)[2]<-"Trait";

    # removing all entries not in tree
    table_tmp2<-data.table(table_tmp)
    setkey(table_tmp2,ID)
    table2<-table_tmp2[intersect(table_tmp2$ID,root_tree$tip.label)]
    setkey(table2,ID)

    #initializing result vectors and file names
    positives<-vector(mode="list",length=0)
    cluster_size<-numeric(length=0)

    #!!!!!!!!!!!! FILE NAMES CHANGES
    cluster_size_file<-paste0(out_dir,"R_cluster_size_",j,".txt")
    cluster_dist<-numeric(length=0)
    cluster_dist_file<-paste0(out_dir,"R_cluster_dist_",j,".txt")

    #initalizing files
    if (m == 1) {
      cat(c("trait","tree","distance","cluster_size"), file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
      cat("\n", file = cluster_size_file, fill = FALSE, labels = NULL,append = TRUE)

      cat(c("trait","tree","distance"), file = cluster_dist_file, sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
      cat("\n", file = cluster_dist_file, fill = FALSE, labels = NULL,append = TRUE)
    }

    #loop through all subtrees and determining if any subtrees have >90% positives
    for (i in 1:length(subtree)){
      tip_names<-subtree[[i]]$tip.label
      if (mean(table2[tip_names][,Trait]) > 0.9 ) {#change this value if you want a new threshold
        match_test<-match(tip_names,positives)
        if (all(is.na(match_test))) {

          #FDB store the node and trait information
          nodeMatch <- rbind(nodeMatch, c('tree'=tree.name, 'subtree'=i, 'trait_nr' =j, 'conserved_in'=mean(table2[tip_names][,Trait])))
          positives<-c(positives,tip_names)

          #!!!! ALTHOUGH I WANT THAT?, we can use any alternatives, and speed up here
          cluster_dist<-distRoot(subtree[[i]],tip_names, method=c("p"))
          cluster_size<-c(cluster_size,mean(cluster_dist))

          # printing to files###
          cat(c(j,m,mean(cluster_dist),length(cluster_dist)), file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
          cat("\n", file = cluster_size_file, fill = FALSE, labels = NULL,append = TRUE)

          cat(j,m,cluster_dist, file = cluster_dist_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
          cat("\n", file = cluster_dist_file, fill = FALSE, labels = NULL,append = TRUE)


          #print(cluster_dist)
        }
        else if (any(is.na(match_test))) {
          print("some NAs - something is weird")
        }
        else {
          #print(tip_names)
          #print("found cluster before")
        }
      }
    }

    ##### find singletons ######
    a<-table2[table2$Trait == 1,][,ID]
    g<-as.character(a)

    singletons_names = setdiff(g,positives)
    if (length(singletons_names) > 0) {
      for (h in 1:length(singletons_names)){
        # weigh singletons with half
        singleton_edges = 0.5*root_tree$edge.length[which.edge(root_tree,singletons_names[h])] #here we use half the distance for singletons
        cluster_size<-c(cluster_size,singleton_edges)

        cat(c(j,m,singleton_edges,1), file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
        cat("\n", file = cluster_size_file, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
      }

    }
    Mean_all[j-1,m] = mean(cluster_size)
  }


}
#output file
write.table(Mean_all,paste0(out_dir,"Mean_all_bootstrap2.txt"), sep = "\t")




Mean_all <- read.delim('~/DATA/MarbGenomics/consenTRAIT/OG.pan/Mean_all_bootstrap2.txt')
hist(Mean_all$V1)
colnames(Mean_all)[1] <- 'depth'
Mean_all <- Mean_all %>% select(depth) %>% mutate(OG = colnames(OG.pan))


Mean_all %>% ggplot(aes(depth)) + geom_histogram()

annotation.tbl %>% filter(genome=='Marinobacter_adhaerens_HP15') %>%
  left_join(Mean_all, by='OG') %>% ggplot(aes(start,depth))+geom_line()+xlim(c(1,200000))


sel.hgm <- ANI.hmg.grp.tbl %>% filter(group=='cl1') %>% filter(level=='sc0.8')

annotation.tbl %>% filter(genome=='Marinobacter_adhaerens_HP15') %>%
  left_join(Mean_all, by='OG') %>%
  left_join(sel.hgm,by=c('OG'='feature.id')) %>%
  ggplot(aes(start,-log10(fdr.p.value_enriched_rawcounts))) +
    geom_line()

annotation.tbl %>% filter(genome=='Marinobacter_adhaerens_HP15') %>%
  left_join(Mean_all, by='OG') %>%
  left_join(sel.hgm,by=c('OG'='feature.id')) %>%
  ggplot(aes(depth,-log10(fdr.p.value_enriched_rawcounts))) +
  geom_point()

deep <- Mean_all %>% filter(depth>0.5) %>% select(OG) %>% pull()
shallow <- Mean_all %>% filter(depth<0.1) %>% select(OG) %>% pull()

annotation.tbl %>% filter(OG %in% deep) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(Name, COG, product, OG) %>% data.frame
annotation.tbl %>% filter(OG %in% shallow) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(Name, COG, product, OG) %>% data.frame


#SELECT THE PROX CLUSTERS
proxclusts <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED','.tsv'),sep='\t')



genom = 'Marinobacter_adhaerens_HP15'
sel.group <- 'cl1'
genom = 'Marinobacter_hydrocarbonoclasticus_ATCC_49840'
sel.group <- 'cl5'
genom = 'Marinobacter_hydrocarbonoclasticus_ATCC_49840'
sel.group <- 'cl5'

genom = 'Marinobacter_sp_Arc7_DN_1'
sel.group <- 'cl3'

genom = 'Marinobacter_salarius_SMR5'
genom = 'Marinobacter_sp_FDB33'
genom = 'Marinobacter_sp_HL_58'
genom = 'Marinobacter_algicola_DG893'
genom = 'Marinobacter_salarius_HK15'

sel.group <- 'cl8'

#==============================

cl <- sorted.ANI.cliques %>%
  left_join(genome.tbl,by='genome') %>%
  select(genome, c0.75:cc0.98) %>%
  mutate(group=paste0('cl',sc0.8)) %>%
  filter(group!='clNA') %>%
  select(genome,group)

sel.col = 'sc0.85'
mtd <- genome.tbl %>% left_join(sorted.ANI.cliques, by='genome') %>%
  left_join(out.colors, by='genome') %>%
  data.frame()
mtd$grp <- mtd[,sel.col]
mtd <- mtd %>% mutate(group=paste0('cl',grp))

cl <- mtd %>%
  filter(group!='clNA') %>%
  select(genome,group) %>% arrange(group)

###### need to rerun, as i needed the below code becaiuse sc0.85 was not consistently named between sets, say the sortedclqie and the hgm
# THIS IS HORRIBLE
#ffs
#extract it from the data
conv <- output.tbl %>%
  filter(sel.col=='sc0.85')%>% select(genome,group,sel.col) %>%
  unique() %>%
  left_join(cl,by='genome') %>% mutate(group=group.y) %>% select(group,group.x)

cl <- cl %>% left_join(conv,by='group') %>% filter(!is.na(group.x)) %>% mutate(group=group.x)
#########


for(i in 1:nrow(cl)){
  genom <- as.character(cl[i,'genome'])
  sel.group <- as.character(cl[i,'group'])
  message('analysing: ',sel.group,' ',genom)

  sel.hgm <- ANI.hmg.grp.tbl %>%
    dplyr::filter(group==sel.group) %>%
    dplyr::filter(level=='sc0.85')

  sig.ogs <- sel.hgm %>%
    tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts), names_to='variable') %>%
    dplyr::filter(value<.05) %>%
    dplyr::select(feature.id) %>%
    dplyr::pull() %>% as.character()

  ltags <- annotation.tbl %>% filter(OG %in% sig.ogs) %>%
    filter(genome == genom) %>%
    select(locus_tag) %>% pull()

  if(length(ltags)>0){
    message('ok')
    clusdetect.sig <- proximity_clusters(annotation.tbl, ltags,distance=5000 )
    clusdetect.sig <- clusdetect.sig  %>% data.frame() %>%
      dplyr::group_by(prox.cluster,seqnames) %>%
      dplyr::summarise(minStart=min(start),maxEnd=max(end))

    sel.proxclust <- proxclusts %>% filter(sel.col=='sc0.8') %>%
      filter(group==sel.group) %>%
      filter(pan=='OGpan')

    clusters_c <- sel.proxclust %>%
      dplyr::group_by(prox.cluster,seqnames) %>%
      dplyr::summarise(minStart=min(start),maxEnd=max(end))

    }

  #cluster SCOs
  ltags <- annotation.tbl %>% filter(OG %in% SCO) %>%
    filter(genome == genom) %>%
    select(locus_tag) %>% pull()

  clusdetect.loc <- proximity_clusters(annotation.tbl, ltags,distance=5000 )
  clusdetect.loc <- clusdetect.loc  %>%
    dplyr::group_by(prox.cluster,seqnames) %>%
    dplyr::summarise(minStart=min(start),maxEnd=max(end))

  outplot <- annotation.tbl %>% filter(genome==genom) %>%
    left_join(Mean_all, by='OG') %>%
    left_join(sel.hgm,by=c('OG'='feature.id')) %>%
    mutate(fdr.p.value_enriched_binary = -log10(fdr.p.value_enriched_binary)) %>%
    mutate(fdr.p.value_enriched_rawcounts = -log10(fdr.p.value_enriched_rawcounts)) %>%
    tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, depth), names_to='variable')%>%
    ggplot(aes(start,value)) + geom_point(aes(color=ifelse(value>1.3,TRUE,FALSE)),size=0.2) +
    scale_color_manual(values=c('black','red'))+
    facet_wrap(~variable,ncol=1,scales='free') + #fdb_style(aspect.ratio=0.5)+
    #add the linaege specific ones
    #geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)+
    #add the sco loci?
    geom_segment(data=clusdetect.loc,aes(x = minStart, y = -2, xend = maxEnd, yend = -2),color='black',size=3)
    #add the sco loci?


  if(length(ltags)>0){
    outplot <- outplot + geom_segment(data=clusdetect.sig,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)
  }

    outplot <- outplot + facet_grid(variable ~ seqnames, scales = "free", space='free_x') +
    ggtitle(genom)+ theme_classic() + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 11, colour = '#000000'),
      axis.text = ggplot2::element_text(size = 10, colour = '#000000'),
      #legend.justification = c(1, 1),
      #legend.key.width = unit(0.25, 'cm'),
      #legend.key.height = unit(0.55, 'cm'),
      #legend.text = ggplot2::element_text(size = 10),
      #legend.title = ggplot2::element_text(size = 11),
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size  = 11,  hjust = 0),
    )+theme(legend.position = 'none')


  ggplot2::ggsave(filename=paste0(outdir,'mapping_',sel.group,'_',genom,'.pdf'),plot=outplot, width = 13, height = 4,unit='in')
  }




#=====



#select right importance.all # see in RF


annotation.tbl %>% filter(genome==genom) %>%
  left_join(Mean_all, by='OG') %>%
  left_join(sel.hgm,by=c('OG'='feature.id')) %>%
  left_join(importance.all %>% filter(group == sel.group), by=c('OG'='feature.id')) %>%
  mutate(fdr.p.value_enriched_binary = -log10(fdr.p.value_enriched_binary)) %>%
  mutate(fdr.p.value_enriched_rawcounts = -log10(fdr.p.value_enriched_rawcounts)) %>%
  #tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, depth,MeanDecreaseGini,MeanDecreaseAccuracy), names_to='variable')%>%
  tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts,MeanDecreaseGini,MeanDecreaseAccuracy), names_to='variable')%>%
  ggplot(aes(start,value)) +
  #geom_bar(stat="identity",aes(fill=ifelse(value>1.3,TRUE,FALSE)),size=0.2) +
  geom_point(aes(color=ifelse(value>1.3,TRUE,FALSE)),size=0.2) +
  scale_color_manual(values=c('black','red'))+
  scale_fill_manual(values=c('black','red'))+
  facet_wrap(~variable,ncol=1,scales='free') + facet_grid(variable ~ seqnames, scales = "free", space='free_x') +
  ggtitle(genom)+ theme_classic() + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    axis.title = ggplot2::element_text(size = 11, colour = '#000000'),
    axis.text = ggplot2::element_text(size = 10, colour = '#000000'),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size  = 11,  hjust = 0),
  )+theme(legend.position = 'none')
#+ #fdb_style(aspect.ratio=0.5)+
  #add the linaege specific ones
  #geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)+
  #add the sco loci?
  #geom_segment(data=clusdetect.loc,aes(x = minStart, y = -2, xend = maxEnd, yend = -2),color='black',size=3)
#add the sco loci?




annotation.tbl %>% filter(genome==genom) %>%
  left_join(Mean_all, by='OG') %>%
  left_join(sel.hgm,by=c('OG'='feature.id')) %>%
  left_join(importance.all %>% filter(group == sel.group), by=c('OG'='feature.id')) %>%
  mutate(fdr.p.value_enriched_binary = -log10(fdr.p.value_enriched_binary)) %>%
  mutate(fdr.p.value_enriched_rawcounts = -log10(fdr.p.value_enriched_rawcounts)) %>%
  #tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, depth,MeanDecreaseGini,MeanDecreaseAccuracy), names_to='variable')%>%
  tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts,MeanDecreaseGini,MeanDecreaseAccuracy), names_to='variable')%>%
  ggplot(aes(start,value)) + geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)#+


annotation.tbl %>% dplyr::filter(genome==genom) %>% mutate(CORE=ifelse(feature.id %in% CORE.ogs))






#=================
# OGs TABLE!!!!!!!!

og.table <- sel.hgm %>%
  dplyr::filter(group == sel.group) %>%
  dplyr::left_join(importance.all %>%
              dplyr::filter(group == sel.group) %>%
              dplyr::select(-group),
            by='feature.id') %>%
  dplyr::left_join(df_counts, by='feature.id') %>%
  dplyr::left_join(shared.df %>% tibble::rownames_to_column(var='feature.id'), by='feature.id') %>%
  dplyr::mutate(CORE=ifelse(feature.id %in% CORE.ogs,'core','accessory'))






annotation.tbl %>% filter(genome==genom) %>%
#  dplyr::mutate(ifelse(product)) %>%
  ggplot()+
  geom_segment(aes(x = start, y = type, xend = end, yend = type,color=type),size=5)+ geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)#+



annotation.tbl %>% filter(genome==genom) %>%
  ggplot()+
  geom_segment(aes(x = start, y = type, xend = end, yend = type,color=type),size=5)



annotation.tbl %>% filter(genome==genom) %>%
  dplyr::mutate(tsp = ifelse(grepl('transposase',product),'tsp',NA)) %>%
  ggplot()+
  geom_segment(aes(x = start, y = tsp, xend = end, yend = tsp,color=tsp),size=5)+ geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)#+



#============================#
#### TESTING THE PREMILINARY IDEA OF SYNTENY COMPASIRON #####
# SCORE BLOCKS
# Level of synteny conservation per block
# gene order conservation scores using a sliding window (8CDS, 1 CDS step) - each profile reference against query


reference <- "Marinobacter_adhaerens_HP15"
query <- 'Marinobacter_sp_CP1'
sliding.window <- 8
sliding.frame <- 1




GOC_scores <- function(reference, query, sliding.window, sliding.frame, chromosome = FALSE){
  start_time <- Sys.time()

  # extract the reference table
  if(chromosome==TRUE){
    longest_contig <- annotation.tbl %>% filter(genome %in% reference) %>% select(seqnames) %>% group_by(seqnames) %>% tally() %>% arrange(desc(n)) %>% select(seqnames) %>% slice(1) %>% pull

    ref.dat <- annotation.tbl %>%
      dplyr::filter(genome==reference) %>%
      dplyr::filter(seqnames==longest_contig) %>%
      dplyr::filter(type=='CDS') %>%
      dplyr::select(genome, seqnames, start, end, OG, type, locus_tag, product) %>%
      dplyr::mutate(pos = 1:nrow(.))
  }else{
    ref.dat <- annotation.tbl %>%
      dplyr::filter(genome==reference) %>%
      dplyr::filter(type=='CDS') %>%
      dplyr::select(genome, seqnames, start, end, OG, type, locus_tag, product) %>%
      dplyr::mutate(pos = 1:nrow(.))
  }

  qr.dat <- annotation.tbl %>%
    dplyr::filter(genome==query)# %>%
    #dplyr::filter(type=='CDS') %>%
    #dplyr::select(genome, seqnames, start, end, OG, type, locus_tag, product)# %>%
    #dplyr::mutate(pos = 1:nrow(.))

  #set last position
  last.pos <- ref.dat[nrow(ref.dat),'pos'] %>% as.numeric() - 8

  out.comb <- c()
  # run over all genes, with increments = sliding.frame
  for(i in seq(1,last.pos,sliding.frame)){
    #if(i%%100==0){
    #  message('analysing ',i)
    #}
    #select genes in current window
    selected.OGs <- ref.dat %>% slice(i:(i+sliding.window)) %>% select(OG) %>% filter(!is.na(OG)) %>% pull()

    #extract corresponding locus_tags form the query genome
    ltags <- qr.dat %>% filter(OG %in% selected.OGs) %>%
      filter(genome == query) %>%
      select(locus_tag) %>% pull()

    if(length(ltags > 0)){
      clusdetect.sig <- proximity_clusters(qr.dat, ltags,distance=1000 )
      top.cluster <- clusdetect.sig %>% select(locus_tag, prox.cluster) %>% group_by(prox.cluster) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(prox.cluster) %>% pull
      out.comb <- rbind(out.comb,
                        c('pos'=i,
                          'ref_locus_tag' =   as.character(ref.dat[i,'locus_tag']),
                          'ref_size' = length(selected.OGs),
                          'query_size' = clusdetect.sig %>% filter(prox.cluster == top.cluster) %>% nrow()
                        )
      )
    }else{
      out.comb <- rbind(out.comb,
                        c('pos'=i,
                          'ref_locus_tag' =   as.character(ref.dat[i,'locus_tag']),
                          'ref_size' = length(selected.OGs),
                          'query_size' = 0
                        )
      )
    }
    #cluster the query locus_tags and extract the longest cluster

    #output some data
  }

  out.comb <- out.comb %>%
    data.frame() %>%
    mutate(pos = as.numeric(as.character(pos)),
           query_size = as.numeric(as.character(query_size)),
           ref_size = as.numeric(as.character(ref_size)),
           score=query_size/ref_size)

  end_time <- Sys.time()
  message(end_time-start_time)
  return(out.comb)
}



GOC.scores.df <-  GOC_scores(reference = "Marinobacter_adhaerens_HP15",
                             query = 'Marinobacter_sp_CP1',
                             sliding.window = 8,
                             sliding.frame = 1)

GOC.scores.df %>%  mutate(score = ifelse(score>1,1,score)) %>% ggplot(aes(pos,1))+geom_point(aes(color=score)) + scale_color_viridis()


# Run over multiple queries
queries <- c('Marinobacter_sp_DSM_26671','Marinobacter_sp_CP1','Marinobacter_sp_EhC06','Marinobacter_sp_N4','Marinobacter_hydrocarbonoclasticus_VT8')

queries <-genome.tbl %>% filter(phylogroup=='P1' & quality_class=='high' & Completeness==100) %>% select(genome) %>% pull()
queries <-queries[1:10]






queries <- genome.tbl %>% filter(TypeStrain==TRUE & quality_class=='high') %>% select(genome) %>% pull()
reference <- "Marinobacter_hydrocarbonoclasticus_VT8"



queries <-genome.tbl %>% filter(phylogroup=='P1' & quality_class=='high' & Completeness==100) %>% select(genome) %>% pull()
tpstr <- genome.tbl %>% filter(TypeStrain==TRUE & quality_class=='high') %>% select(genome) %>% pull()
queries <- unique(c(queries,tpstr))
reference <- "Marinobacter_adhaerens_HP15"
queries <- queries[queries!=reference]


GOC.many.df <-c()
for(qr in queries){
  message(qr)
  GOC.scores.df <- GOC_scores(reference = reference,
                               query = qr,
                               sliding.window = 8,
                               sliding.frame = 1)
  GOC.scores.df <- GOC.scores.df %>% mutate(query_genome = qr)
  GOC.many.df <- rbind(GOC.many.df, GOC.scores.df)
}


write.table(x=GOC.many.df, file='GOC.many.df_HP15_v_p1_tpstr.txt',row.names = FALSE, sep='\t',quote = FALSE)






ref.dat <- annotation.tbl %>%
  dplyr::filter(genome==reference) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::select(genome, seqnames, start, end, OG, type, locus_tag, product,seqnames) %>%
  dplyr::mutate(pos = 1:nrow(.))

GOC.many.df %>% left_join(ref.dat, by=c('ref_locus_tag'='locus_tag')) %>%
mutate(score = ifelse(score>1,1,score)) %>%
  ggplot()+
geom_segment(aes(x = start, y = query_genome, xend = end, yend = query_genome, color=score),size=5)+
  scale_color_viridis()+facet_wrap(~seqnames,ncol=1)


nms.order <- names(sort(COPH[reference,queries]))
nms.order <- GOC.many.df %>% group_by(query_genome) %>% dplyr::summarise(sum=sum(score)) %>% dplyr::arrange(desc(sum)) %>% select(query_genome) %>% pull %>% as.character()
GOC.many.df$query_genome = factor(GOC.many.df$query_genome, levels=rev(nms.order))





#-------------------------------
# PLOTTING COMPOSITE GRAPH
#---------------------------------

p5 <- annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  select(seqnames, genome,start,end,OG,P1,P2,P3,P4,P5,P6,P7,P9,P10,P11) %>%
  tidyr::pivot_longer(cols=c(P1:P11),names_to='phgrp') %>%
  ggplot()+
  geom_segment(aes(x = start, y = phgrp, xend = end, yend = phgrp,color=value),size=5)+
  scale_x_continuous(limits=c(0,4000000))+facet_wrap(~seqnames,ncol=1)+
  scale_color_manual(values=c('grey','darkblue'),na.value=NA)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())




p4 <- annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  mutate(persist.core = ifelse(nrOfSpecies==94,TRUE,FALSE),
         persist.95 = ifelse(nrOfSpecies>94*0.95,TRUE,FALSE),
         persist.90 = ifelse(nrOfSpecies>94*0.9,TRUE,FALSE),
         persist.75 = ifelse(nrOfSpecies>94*0.75,TRUE,FALSE),
         persist.50 = ifelse(nrOfSpecies>94*0.50,TRUE,FALSE),
         persist.25 = ifelse(nrOfSpecies>94*0.25,TRUE,FALSE),
         persist.10 = ifelse(nrOfSpecies>94*0.10,TRUE,FALSE),
         a.rare.10 =  ifelse(nrOfSpecies<94*0.10,TRUE,FALSE)

  ) %>%
  select(seqnames, genome,start,end,OG,persist.core:a.rare.10) %>%
  tidyr::pivot_longer(cols=c(persist.core:a.rare.10),names_to='pan_class') %>%
  dplyr::filter(value==TRUE) %>%
  ggplot()+
  geom_segment(aes(x = start, y = pan_class, xend = end, yend = pan_class,color=value),size=5)+
  scale_x_continuous(limits=c(0,4000000))+facet_wrap(~seqnames,ncol=1)+
  scale_color_manual(values=c('firebrick','skyblue'),na.value=NA)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())



p3 <- annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  geom_segment(aes(x = start, y = CORE, xend = end, yend = CORE,color=CORE),size=5)+
  scale_x_continuous(limits=c(0,4000000))+facet_wrap(~seqnames,ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())


annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  #geom_segment(aes(x = start, y = CORE, xend = end, yend = CORE,color=CORE),size=5)+ geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)+
  geom_point(aes(x=start,y=nrOfGenes)) +
  scale_x_continuous(limits=c(0,4000000))+facet_wrap(~seqnames,ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())

p6<-annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  #geom_segment(aes(x = start, y = CORE, xend = end, yend = CORE,color=CORE),size=5)+ geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)+
  geom_point(aes(x=start,y=nrOfSpecies)) +
  scale_x_continuous(limits=c(0,4000000))+facet_wrap(~seqnames,ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())


p6<-annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  geom_segment(aes(x = start, y = 1, xend = end, yend = 1, color=nrOfSpecies/94),size=5)+
  #geom_segment(data=clusters_c,aes(x = minStart, y = -1, xend = maxEnd, yend = -1),color='red',size=3)+
  #geom_point(aes(x=start,y=nrOfSpecies/94)) +
  scale_x_continuous(limits=c(0,4000000))+facet_wrap(~seqnames,ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())+viridis::scale_color_viridis()


chrsize <- annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  dplyr::filter(seqnames == 'CP001978') %>% arrange(desc(end)) %>% slice(1) %>% select(end) %>% pull()


annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  #geom_segment(aes(x=0,y=1,xend=chrsize+1,yend=1),color='black',size=6)+
  geom_segment(aes(x=1,y=1,xend=chrsize,yend=1),color='grey',size=5)+
  geom_segment(aes(x = start, y = 1, xend = end, yend = 1, color=nrOfSpecies/94),size=5)+
  scale_x_continuous(limits=c(0,gsize))+facet_wrap(~seqnames,ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())+
  viridis::scale_color_viridis()+ theme_void()




#=================================================================#
#### KEY FIGURE ####

p.prox.1 <-annotation.tbl %>% filter(genome==genom) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!is.na(CORE)) %>%
  filter(CORE=='core') %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  geom_segment(aes(x=1,y=1,xend=chrsize ,yend=1),color='grey90',size=5)+
  geom_segment(data=clusters_c,aes(x = minStart,y = 1, xend = maxEnd, yend = 1),
               color=phylogroup.colors['P1'],size=5)+
  scale_x_continuous(limits=c(0,chrsize))+facet_wrap(~seqnames,ncol=1)+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())+
  viridis::scale_color_viridis()+ theme_void()

p.significance.1 <- annotation.tbl %>% filter(genome==genom) %>%
  #left_join(Mean_all, by='OG') %>%
  left_join(sel.hgm,by=c('OG'='feature.id')) %>%
  left_join(importance.all %>% filter(group == sel.group), by=c('OG'='feature.id')) %>%
  mutate(fdr.p.value_enriched_binary = -log10(fdr.p.value_enriched_binary)) %>%
  mutate(fdr.p.value_enriched_rawcounts = -log10(fdr.p.value_enriched_rawcounts)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts,MeanDecreaseAccuracy), names_to='variable')%>%
  ggplot() +
  geom_segment(aes(x=start,xend=end,y=0,yend=value,color=ifelse(value>1.3,TRUE,FALSE)),size=0.15) +
  scale_color_manual(values=c('black',"firebrick"))+
  scale_fill_manual(values=c('black',"firebrick"))+
  facet_grid(variable ~ seqnames, scales = "free", space='free_x') +
  guides(color=guide_legend(title='FDRp<.05'))+
  scale_x_continuous(limits=c(0,chrsize))+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + theme_classic()+
  theme(strip.background=element_blank(),strip.text.x=element_blank())#+ fdb_style(aspect.ratio = 0.2)


ggpubr::ggarrange(p.prox,p.significance,align='v',ncol=1)

p.combined.s <- cowplot::plot_grid(p.prox,p.significance,align='v',axis='rl',ncol=1,rel_heights = c(0.3,0.7),labels = c('a','b'))

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",'Genome_mapping_',sel.group,'_',genom,'.pdf'),
                plot=p.combined.s,
                width = 8, height = 4,unit='in')





#======================#
listedRepresentatives


out.clstall <- c()
for(i in 1:nrow(listedRepresentatives)){
  genom <- as.character(listedRepresentatives[i,'gnm'])
  sel.group <- as.character(listedRepresentatives[i,'sel.grp'])
  message('analysing: ',sel.group,' ',genom)

  selected.contigs <- annotation.tbl %>%
    filter(OG %in% sig.ogs) %>%
    filter(genome == genom) %>%
    select(seqnames) %>%
    unique() %>%
    pull

  if(length(selected.contigs)<4){
  sel.hgm <- ANI.hmg.grp.tbl %>%
    dplyr::filter(group==sel.group) %>%
    dplyr::filter(level=='sc0.85')

  sig.ogs <- sel.hgm %>%
    tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts), names_to='variable') %>%
    dplyr::filter(value<.05) %>%
    dplyr::select(feature.id) %>%
    dplyr::pull() %>% as.character()

  ltags <- annotation.tbl %>% filter(OG %in% sig.ogs) %>%
    filter(genome == genom) %>%
    select(locus_tag) %>% pull()

  if(length(ltags)>0){
    message('ok')
    clusdetect.sig <- proximity_clusters(annotation.tbl, ltags,distance=5000 )
    clusdetect.sig <- clusdetect.sig  %>% data.frame() %>%
      dplyr::group_by(prox.cluster,seqnames) %>%
      dplyr::summarise(minStart=min(start),maxEnd=max(end))

    sel.proxclust <- proxclusts %>% filter(sel.col=='sc0.8') %>%
      filter(group==sel.group) %>%
      filter(pan=='OGpan')

    clusters_c <- sel.proxclust %>%
      dplyr::group_by(prox.cluster,seqnames) %>%
      dplyr::summarise(minStart=min(start),maxEnd=max(end))

  }



  #for each contig

    for(selected.contig in selected.contigs){
      #selct chromosome size
      chrsize <- annotation.tbl %>% filter(genome==genom) %>%
        #left_join(og.table ,by=c('OG'='feature.id')) %>%
        dplyr::filter(seqnames == selected.contig) %>%
        arrange(desc(end)) %>% slice(1) %>% select(end) %>% pull()



      s.col <- col.list[[sel.col]][[sel.group]]

      p.significance <- annotation.tbl %>%
        filter(genome==genom) %>%
        left_join(sel.hgm,by=c('OG'='feature.id')) %>%
        left_join(importance.all %>% filter(group == sel.group), by=c('OG'='feature.id')) %>%
        mutate(fdr.p.value_enriched_binary = -log10(fdr.p.value_enriched_binary)) %>%
        mutate(fdr.p.value_enriched_rawcounts = -log10(fdr.p.value_enriched_rawcounts)) %>%
        dplyr::filter(seqnames == selected.contig) %>%
        tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts,MeanDecreaseAccuracy), names_to='variable')%>%
        ggplot() +
        geom_segment(aes(x=start,xend=end,y=0,yend=value,color=ifelse(value>1.3,TRUE,FALSE)),size=0.15) +
        scale_color_manual(values=c('black',"firebrick"))+
        scale_fill_manual(values=c('black',"firebrick"))+
        facet_grid(variable ~ seqnames, scales = "free", space='free_x') +
        guides(color=guide_legend(title='FDRp<.05'))+
        scale_x_continuous(limits=c(0,chrsize))+
        theme(axis.title.x=element_blank(), axis.text.x = element_blank()) + theme_classic()+
        theme(strip.background=element_blank(),strip.text.x=element_blank())#+ fdb_style(aspect.ratio = 0.2)


      p.prox <- annotation.tbl %>%
        filter(genome==genom) %>%
        left_join(og.table ,by=c('OG'='feature.id')) %>%
        dplyr::filter(seqnames == selected.contig) %>%
        ggplot()+
        geom_segment(aes(x=1,y=1,xend=chrsize,yend=1),color='grey90',size=5)+
        geom_segment(data=clusdetect.sig,aes(x = minStart,y = 1, xend = maxEnd, yend = 1),
                     color=s.col,size=5)+
        scale_x_continuous(limits=c(0,chrsize))+
        theme(axis.title.x=element_blank(), axis.text.x = element_blank())+
        viridis::scale_color_viridis()+ theme_void()+
        labs(subtitle=paste0(selected.contig, ' (', round(chrsize/1000000,2),' Mb) ', nrow(clusdetect.sig),' clusters'))


      p.combined.s <- cowplot::plot_grid(p.prox,p.significance,align='v',axis='rl',ncol=1,rel_heights = c(0.4,0.6),labels = c('a','b'))

      ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",'Genome_mapping_',sel.group,'_',genom,'_',selected.contig,'.pdf'),
                      plot=p.combined.s,
                      width = 8, height = 4,unit='in')


      out.clst <- clusdetect.sig %>% filter(seqnames %in% selected.contig) %>% mutate(genome=genom,group=sel.group,chrsize=chrsize,nsig=nrow(clusdetect.sig)) %>% ungroup()
      out.clstall <- rbind(out.clstall, out.clst)

    }
  }
}




#======================#


p.reference <- out.clstall %>% filter(genome %in% c('Marinobacter_psychrophilus_20041','Marinobacter_sp_LV10R510_8','Marinobacter_adhaerens_HP15','Marinobacter_salarius_HK15','Marinobacter_hydrocarbonoclasticus_VT8')) %>%
  filter(chrsize>2000000) %>%
  ggplot() +
  geom_segment(aes(x=1,y=seqnames,xend=chrsize,yend=seqnames),color='grey90',size=5)+
  geom_segment(aes(x = minStart,y = seqnames, xend = maxEnd, yend = seqnames,color=group),
               size=5)+theme_classic()+scale_color_manual(values=col.list[['sc0.85']])+theme(legend.position = 'none')#+


ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",'Genome_mapping_references',sel.group,'.pdf'),
                plot=p.reference,
                width = 8, height = 2,unit='in')



p.combined.s <- cowplot::plot_grid(p.prox,p.significance,p.reference,align='v',axis='rl',ncol=1,rel_heights = c(0.1,0.5,0.5),labels = c('a','b','c'))

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",'Genome_mapping_references_all',sel.group,'.pdf'),
                plot=p.combined.s,
                width = 9, height = 4,unit='in')















p1 <- GOC.many.df %>% left_join(ref.dat, by=c('ref_locus_tag'='locus_tag')) %>%
  mutate(score = ifelse(score>1,1,score)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  ggplot()+
  geom_segment(aes(x = start, y = query_genome, xend = end, yend = query_genome, color=score),size=5)+
  viridis::scale_color_viridis()+
  facet_wrap(~seqnames,ncol=1)+
  scale_x_continuous(limits=c(0,4000000))

p2 <- annotation.tbl %>% filter(genome==genom) %>%
  left_join(Mean_all, by='OG') %>%
  left_join(sel.hgm,by=c('OG'='feature.id')) %>%
  left_join(importance.all %>% filter(group == sel.group), by=c('OG'='feature.id')) %>%
  mutate(fdr.p.value_enriched_binary = -log10(fdr.p.value_enriched_binary)) %>%
  mutate(fdr.p.value_enriched_rawcounts = -log10(fdr.p.value_enriched_rawcounts)) %>%
  dplyr::filter(seqnames == 'CP001978') %>%
  #tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, depth,MeanDecreaseGini,MeanDecreaseAccuracy), names_to='variable')%>%
  tidyr::pivot_longer(cols=c(fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts,MeanDecreaseGini,MeanDecreaseAccuracy), names_to='variable')%>%
  ggplot(aes(start,value)) +
  #geom_bar(stat="identity",aes(fill=ifelse(value>1.3,TRUE,FALSE)),size=0.2) +
  geom_point(aes(color=ifelse(value>1.3,TRUE,FALSE)),size=0.2) +
  scale_color_manual(values=c('black','red'))+
  scale_fill_manual(values=c('black','red'))+
  #facet_wrap(~variable,ncol=1,scales='free') +
  facet_grid(variable ~ seqnames, scales = "free", space='free_x') +
  guides(color=guide_legend(title='FDRp<.05'))+
  scale_x_continuous(limits=c(0,4000000))+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank())



cowplot::plot_grid(p3,p4,p2,p5,p1,align='v',axis='rl',ncol=1,rel_heights = c(0.5,0.5,1,0.5,1))

cowplot::plot_grid(p3,p4,p2,p5,p1,p6,align='v',axis='rl',ncol=1,rel_heights = c(0.5,0.5,1,0.5,1))

cowplot::plot_grid(p3,p4,p2,p6,p1,align='v',axis='rl',ncol=1,rel_heights = c(0.5,0.5,1,0.5,1))




#write.table(x=GOC.many.df, file='GOC.many.df_VT8_v_p5.txt',row.names = FALSE, sep='\t',quote = FALSE)
#write.table(x=GOC.many.df, file='GOC.many.df_VT8_v_types.txt',row.names = FALSE, sep='\t',quote = FALSE)


GOC.many.df.2 <- read.delim('GOC.many.df_VT8_v_p5.txt')
GOC.combined <- rbind(GOC.many.df, GOC.many.df.2)

nms.order <- GOC.combined %>% mutate(score = ifelse(score>1,1,score)) %>% group_by(query_genome) %>% dplyr::summarise(sum=sum(score)) %>% dplyr::arrange(desc(sum)) %>% select(query_genome) %>% pull %>% as.character()
#nms.order <- names(sort(COPH[reference,queries]))
GOC.combined$query_genome = factor(GOC.combined$query_genome, levels=nms.order)

GOC.combined %>% left_join(ref.dat, by=c('ref_locus_tag'='locus_tag')) %>%
  mutate(score = ifelse(score>1,1,score)) %>%
  #dplyr::filter(seqnames == 'NC_008740') %>%
  ggplot()+
  geom_segment(aes(x = start, y = query_genome, xend = end, yend = query_genome, color=score),size=5)+
  scale_color_viridis()+facet_wrap(~seqnames,ncol=1)









clusdetect.sig <- clusdetect.sig  %>% data.frame() %>%
  dplyr::group_by(prox.cluster,seqnames) %>%
  dplyr::summarise(minStart=min(start),maxEnd=max(end))









clusdetect.sig <- proximity_clusters(annotation.tbl, ltags,distance=5000 )
clusdetect.sig <- clusdetect.sig  %>% data.frame() %>%
  dplyr::group_by(prox.cluster,seqnames) %>%
  dplyr::summarise(minStart=min(start),maxEnd=max(end))


#===============
annotation.tbl %>% filter(genome %in% sel.genome) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(genome=='Marinobacter_adhaerens_HP15') %>%
  #dplyr::filter(seqnames == 'CP001978') %>%
  mutate(persist.core = ifelse(nrOfSpecies==94,TRUE,FALSE),
         persist.95 = ifelse(nrOfSpecies>94*0.95,TRUE,FALSE),
         persist.90 = ifelse(nrOfSpecies>94*0.9,TRUE,FALSE),
         persist.75 = ifelse(nrOfSpecies>94*0.75,TRUE,FALSE),
         persist.50 = ifelse(nrOfSpecies>94*0.50,TRUE,FALSE),
         persist.25 = ifelse(nrOfSpecies>94*0.25,TRUE,FALSE),
         persist.10 = ifelse(nrOfSpecies>94*0.10,TRUE,FALSE),
         a.rare.10 =  ifelse(nrOfSpecies<94*0.10,TRUE,FALSE)

  ) %>% filter(persist.core ==TRUE)


annotation.tbl %>% filter(genome %in% sel.genome) %>%
  left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(genome=='Marinobacter_adhaerens_HP15') %>%
  #dplyr::filter(seqnames == 'CP001978') %>%
  mutate(persist.core = ifelse(nrOfSpecies==94,TRUE,FALSE),
         persist.95 = ifelse(nrOfSpecies>94*0.95,TRUE,FALSE),
         persist.90 = ifelse(nrOfSpecies>94*0.9,TRUE,FALSE),
         persist.75 = ifelse(nrOfSpecies>94*0.75,TRUE,FALSE),
         persist.50 = ifelse(nrOfSpecies>94*0.50,TRUE,FALSE),
         persist.25 = ifelse(nrOfSpecies>94*0.25,TRUE,FALSE),
         persist.10 = ifelse(nrOfSpecies>94*0.10,TRUE,FALSE),
         a.rare.10 =  ifelse(nrOfSpecies<94*0.10,TRUE,FALSE)

  ) %>% filter(persist.core ==TRUE) %>%
  select(locus_tag,product,gene,COG, Name, cazy_domain) %>%
  data.frame() %>% arrange(product) %>%
  filter(grepl('methionine',product,ignore.case=TRUE))







annotation.tbl %>%
  #left_join(og.table ,by=c('OG'='feature.id')) %>%
  filter(!grepl('Marinobacter',genome)) %>%
  select(genome, locus_tag,product,gene,COG, Name, cazy_domain) %>%
  data.frame() %>% arrange(product) %>%
  filter(grepl('methionine',product,ignore.case=TRUE)) %>% filter(grepl('met',Name))

annotation.tbl %>%
  mutate(mb = ifelse(grepl('Marinobacter',genome),T,F)) %>%
  select(locus_tag,product,gene,COG, Name, cazy_domain,mb) %>%
  data.frame() %>% arrange(product) %>%
  filter(grepl('methionine',product,ignore.case=TRUE)) %>% group_by(product,mb) %>% tally() %>% arrange(product) #%>% filter(grepl('met',Name))

annotation.tbl %>%
  filter(!grepl('Marinobacter',genome)) %>%
  select(locus_tag,product,gene,COG, Name, cazy_domain) %>%
  data.frame() %>% arrange(product) %>%
  filter(grepl('methionine',product,ignore.case=TRUE)) %>% group_by(product) %>% tally() #%>% filter(grepl('met',Name))



annotation.tbl %>%
  mutate(mb = ifelse(grepl('Marinobacter',genome),T,F)) %>%
  select(locus_tag,product,gene,COG, Name, cazy_domain,mb) %>%
  data.frame() %>% arrange(product) %>%
  filter(grepl('methionine',product,ignore.case=TRUE)) %>% group_by(product,mb) %>% tally() %>% arrange(product) %>%
  ggplot(aes(n,product,color=mb))+geom_point()#%>% filter(grepl('met',Name))





#===============fg


#---------------




#sel.proxclust %>% ggplot()+
#  geom_segment(aes(x = start, y = seqnames, xend = end, yend = seqnames, colour = "segment"))


clusters_c %>% ggplot()+ geom_segment(aes(x = minStart, y = 1, xend = maxEnd, yend = 1, colour = "segment"),size=20)


annotation.tbl %>% filter(genome=='Marinobacter_hydrocarbonoclasticus_VT8') %>%
  left_join(Mean_all, by='OG') %>% ggplot(aes(start,depth))+geom_point(size=0.2)+facet_wrap(~seqnames,ncol=1)

annotation.tbl %>% filter(genome=='Marinobacter_algicola_DG893') %>%
  group_by(seqnames) %>% slice(1:3)%>% ungroup()%>%
  left_join(Mean_all, by='OG') %>% ggplot(aes(start,depth))+geom_point(size=0.2)+facet_wrap(~seqnames,ncol=1)


s.gnms <- cl %>% filter(group.x=='cl11') %>% select(genome) %>% pull
sel.tbl <- annotation.tbl%>% filter(genome %in% c('Marinobacter_adhaerens_HP15','Marinobacter_sp_CP1'))
sel.tbl <- annotation.tbl%>% filter(genome %in% c('Marinobacter_hydrocarbonoclasticus_ATCC_49840','Marinobacter_hydrocarbonoclasticus_VT8'))
sel.tbl <- annotation.tbl%>% filter(genome %in% s.gnms)
sel.tbl <- annotation.tbl%>% filter(genome %in% c('Marinobacter_salarius_HK15','Marinobacter_salarius_SMR5'))
s.gnms = c('Marinobacter_hydrocarbonoclasticus_ATCC_49840','Marinobacter_hydrocarbonoclasticus_VT8',
           'Marinobacter_salarius_HK15','Marinobacter_salarius_SMR5',
           'Marinobacter_adhaerens_HP15','Marinobacter_sp_CP1')
sel.tbl <- annotation.tbl%>% filter(genome %in% s.gnms)



sel.tbl %>% mutate(OG=factor(OG,levels = unique(sel.tbl$OG))) %>%
  mutate(genome=factor(genome, levels=s.gnms)) %>%
  ggplot(aes(start,genome))+geom_point(aes(color=OG)) +
  scale_color_viridis(discrete=TRUE)+
  theme(legend.position='none') + theme_minimal()





#=====

sig.ogs


pcor <- cor(OG.pan[,sig.ogs],method='spearman')

gplots::heatmap.2(as.matrix(pcor),trace='none',col=rev(seriation::bluered(16)),symm=T)

corrplot::corrplot(as.matrix(pcor))
#======



nodeMatch <- nodeMatch %>% data.frame(stringsAsFactors = FALSE) %>%
  mutate(subtree=as.numeric(subtree)) %>%
  mutate(trait_nr=as.numeric(trait_nr))

ntips <- length(lsTrees[[1]]$tip.label)

nodeMatch.count <- nodeMatch %>% filter(tree =='concat.sco.nuc') %>% dplyr::count(subtree)
nodeMatch.count <- nodeMatch.count %>% mutate(node = subtree + ntips)

nodeMatch.count <- nodeMatch.count %>% filter(subtree!=1)


nodeMatch %>% ggplot(aes(trait_nr)) + geom_histogram()



write.table(nodeMatch.count,paste0(out_dir,"nodeMatch.count.txt"), sep = "\t")
write.table(nodeMatch,paste0(out_dir,"nodeMatch.txt"), sep = "\t")

nodeMatch <- read.delim('~/DATA/MarbGenomics/consenTRAIT/OG.pan/nodeMatch.txt')
nodeMatch.count <- read.delim('~/DATA/MarbGenomics/consenTRAIT/OG.pan/nodeMatch.count.txt')


root_tree<-root(sco.106.concat.nuc.tree,outgroup,resolve.root=T)


p <- ggtree(root_tree)
p$data <- p$data %>% dplyr::left_join(nodeMatch.count,by='node')

p1 <- p + geom_point(aes(color=n,size=n))

ggsave('~/DATA/MarbGenomics/Graphs/OG_pan_concentrait_per_node.pdf',plot=p1, width = 10, height = 10,unit='in')


#p

p$data %>% filter(isTip==FALSE)
p$data %>% filter(label=='Root')



rownames(Mean_all) <- stored_names[2:length(stored_names)]
outdf <- Mean_all %>% data.frame %>% rownames_to_column(var='feature.id')

outdf %>% arrange(desc(X1)) %>% head(20)
hist(outdf$X1)
outdf[reorder(outdf$X1),]

plot(outdf[,2],outdf[,3])

plot(outdf[,2],outdf[,4])


plot(outdf[,2],outdf[,7])


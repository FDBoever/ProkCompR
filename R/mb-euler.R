#library(eulerr)
#library(ggplot2)
#library(dplyr)

#DEPENDENCIES TO BE FIXED
#col.list initiated elsewhere
#col.list needs to be loaded load colors.R or something sendible
#col.list[['sc0.85']]

#
cog_func_single <- cog_func %>%
  dplyr::group_by(COG) %>%
  dplyr::summarise(func=paste(unique(func),collapse=';'),
                   name=paste(unique(name),collapse=';')) %>%
  dplyr::ungroup()


#
df.pfam_per_locus <- df.pfam.all %>%
  select(seq_id, hmm_acc, hmm_name) %>%
  group_by(seq_id) %>%
  dplyr::summarise(hmm_acc=paste(hmm_acc,collapse=';'),
                   hmm_name=paste(hmm_name,collapse=';'))





# function that converts color to lower alpha value and returns it
color_alpha<-function(color){
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - 50) * 255 / 100)
  return(t.col)
}

#visualise colors
#col_pal <- c('grey93',wesanderson::wes_palettes[[12]],wesanderson::wes_palettes[[5]])
#plot(1:length(col_pal), 1:length(col_pal), col = col_pal, cex = 5, pch=19)



remove_small_groups = TRUE

selcol='sc0.85'
metadata <-genome.tbl %>%
  dplyr::left_join(sorted.cliques, by='genome') %>%
  dplyr::left_join(out.colors, by='genome') %>%
  dplyr::filter(!is.na(!!as.name(selcol))) %>%
  data.frame()
metadata$grp <- metadata[,selcol]

metadata <- metadata %>%
  dplyr::mutate(group=paste0('cl',grp)) %>%
  dplyr::select(-grp)

larger_groups <- metadata %>%
  dplyr::group_by(group) %>%
  dplyr::tally() %>%
  dplyr::filter(n>2) %>%
  dplyr::select(group) %>%
  dplyr::pull()

if(remove_small_groups == TRUE){
        metadata <- metadata %>%
          dplyr::filter(group %in% larger_groups)
}



### Initialise ####
#====================================================================#
#select grouping level, here set to 'sc.85' for ANI > 0.85, can be sc.0.7, sc0.8 etc...
# set pan.name to OGpan, but can be any of the pan tables KOpan, COGpan etc.
sel.col = 'sc0.85'
pan.name = "OGpan"

##data loading####
#Load data from HGM test analysis
#ANI.hmg.grp.tbl <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'hgm_allgroups.tsv'))
ANI.hmg.grp.tbl <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'TRUEhgm_allgroups.tsv')) %>%
  as_tibble()

#Load data from RF
#load(paste0("~/DATA/MarbGenomics/",'multiRF',pan.name,sel.col,'.RData'))
load(paste0("~/DATA/MarbGenomics/",'multiRF',pan.name,sel.col,'TRUE.RData'))

#extract importance values
importance.all <- NULL
importance.wide <-  rownames(RF_perGroup[[1]]$importance)  %>% data.frame()
for(grp in names(RF_perGroup)){
  set <- RF_perGroup[[grp]]$importance %>% data.frame() %>%
    tibble::rownames_to_column(var='feature.id') %>%
    dplyr::mutate(group=grp)

  importance.all<-rbind(importance.all, set)
  s <- set[,c('MeanDecreaseAccuracy','MeanDecreaseGini')]
  colnames(s) <- paste(grp,c('MDA','MDG'),sep='_')
  importance.wide <- cbind(importance.wide, s)
}

rf.all <- importance.all %>%
  mutate(comb = paste0(feature.id,'_',group))


#combine HGM and RF datasets
HGM_RF_COMB <- ANI.hmg.grp.tbl %>%
  filter(level=='sc0.85') %>%
  mutate(comb = paste0(feature.id,'_',group)) %>%
  left_join( rf.all,by='comb')


### Euler plots ####
#============================================================================#
# depends on eulerr library - library(eulerr)

# Set color palette
euler_palette <- c('grey93',color_alpha(wesanderson::wes_palettes[[12]][2]),wesanderson::wes_palettes[[12]][2],color_alpha(wesanderson::wes_palettes[[12]][1]),wesanderson::wes_palettes[[12]][1])

euler_palette <- c('grey93', #RF
                   color_alpha(rgb(174,99,49,max=255)), #EnrichBin
                   rgb(174,99,49,max=255), #Enrich_raw
                   color_alpha(rgb(125,155,166,max=255)), # Deple_bin
                   rgb(125,155,166,max=255)) # deple_raw

# Identify groups
groups <- ANI.hmg.grp.tbl %>%
  dplyr::filter(level=='sc0.85') %>%
  dplyr::select(group) %>%
  unique() %>%
  dplyr::pull() %>%
  as.character() %>%
  sort

#for each group plot Euler
plotlist <- list()
for(sel.grp in groups){
  message(sel.grp)
  top100RF <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
    head(100) %>%
    dplyr::select(feature.id.x) %>%
    dplyr::pull() %>%
    as.character()

  out.mat <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::mutate(top100RF = ifelse(feature.id.x %in% top100RF,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Deple_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Deple_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::select(top100RF,Enrich_bin,Enrich_raw,Deple_bin,Deple_raw) %>%
    dplyr::filter(rowSums(.)!=0) %>%
    data.frame() %>%
    as.matrix()

  fit3 <- eulerr::euler(out.mat, #shape = "ellipse",
  )

  tmp.plot <- plot(fit3,
       quantities = TRUE,
       fill = euler_palette,
       lty = c(1,0,0,0,0),
       labels = NULL,main = sel.grp)

  plotlist[[sel.grp]]<- tmp.plot
}


#combine and save
p.comb <- gridExtra::grid.arrange(grobs=plotlist)
ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/",pan.name,"_plots_comparison_top250RF",'_combined',".pdf"),
                plot=p.comb,
                width = 9,
                height = 9,
                unit='in')






#==================================================================
# Using 250 most important varaibles (MeanDecreaseAccuracy) as identified by Random Forest

plotlist <- list()
for(sel.grp in groups){
  message(sel.grp)
  top250RF <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
    head(250) %>%
    dplyr::select(feature.id.x) %>%
    dplyr::pull() %>%
    as.character()

  out.mat <- HGM_RF_COMB %>%
    filter(group.x==sel.grp) %>%
    dplyr::mutate(top250RF = ifelse(feature.id.x %in% top250RF,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Deple_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Deple_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::select(top250RF,Enrich_bin,Enrich_raw,Deple_bin,Deple_raw) %>%
    dplyr::filter(rowSums(.)!=0) %>%
    data.frame() %>%
    as.matrix()

  fit3 <- euler(out.mat, #shape = "ellipse",
  )


  tmp.plot <- plot(fit3,
                   quantities = TRUE,
                   fill = euler_palette,
                   lty = c(1,0,0,0,0),
                   labels = NULL,main = sel.grp)

  #pdf( paste0("~/DATA/MarbGenomics/Graphs/",pan.name,"_plots_comparison",sel.grp,".pdf"), width=3, height=3)
  #print(tmp.plot)
  #dev.off()

  plotlist[[sel.grp]]<- tmp.plot
}



p.comb <- gridExtra::grid.arrange(grobs=plotlist)
ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/",pan.name,"_plots_comparison_top250RF",'_combined',".pdf"),
                plot=p.comb,
                width = 9,
                height = 9,
                unit='in')




#### Ranked Feature RF plots ####
#==========================================================================#
### CODE WILL GENERATE RANKED FEATURE PLOT with MEAN DECREASE ACCURACY ####
selcol = 'sc0.85'
colcol = 'cc0.85'

#format metadata
metdad <-genome.tbl %>%
  dplyr::left_join(sorted.cliques,by=c('genome'='genome')) %>%
  dplyr::left_join(out.colors, by=c('genome'='genome')) %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(!is.na(!!as.name(selcol)))%>%
  data.frame()
rownames(metdad) <- metdad$genome
metdad$group <- paste0('cl',metdad[,selcol])

#extract colors
colors <- metdad %>%
  dplyr::select(paste0('col.',gsub('sc','cc',selcol)),colcol, group) %>%
  unique() %>%
  data.frame()

grid.colors <- as.character(colors[,1])
names(grid.colors) <- as.character(colors[,'group'])
grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
grid.colors <- c(grid.colors,'clNA'='white')


# Select groups
groups <- HGM_RF_COMB %>%
  dplyr::select(group.x) %>%
  unique() %>%
  dplyr::pull() %>%
  as.character() %>%
  sort()

#rank per group
RF_ranked <- NULL
for(sel.grp in groups){
  tmp <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
    dplyr::slice_max(order_by = MeanDecreaseAccuracy, n = 250) %>%
    dplyr::mutate(rank=1:nrow(.))

  RF_ranked <- rbind(RF_ranked, tmp)
}

#for labeling the lines
tophits<- RF_ranked %>%
  dplyr::group_by(group.x) %>%
  dplyr::slice_max(order_by = MeanDecreaseAccuracy, n = 1)

#make line pplot
p.impor <- RF_ranked %>%
  ggplot2::ggplot(ggplot2::aes(rank,MeanDecreaseAccuracy))+
  ggplot2::geom_line(aes(color=group.x),size=1)+
  ggplot2::xlab('feature rank')+
  ggplot2::scale_color_manual(values=grid.colors)+
  fdb_style()

#add labels
p.impor <- p.impor +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = group.x), data = tophits,
    fontface ="plain", color = "black"
  )

#save pdf
ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/",pan.name,"_varImportance",'_overall_ranked',".pdf"),
                plot=p.impor,
                width = 5,
                height = 5,
                unit='in')






#=========================================================================#
# ==== same plot for Hyper Geometric tests results

HGM_ranked <- NULL
for(sel.grp in groups){
  tmp <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    arrange(desc(-log10(fdr.p.value_enriched_binary))) %>%
    dplyr::slice_max(order_by = -log10(fdr.p.value_enriched_binary), n = 1000) %>%
    mutate(rank=1:nrow(.)) %>% mutate(FDR = fdr.p.value_enriched_binary) %>% mutate(type = 'enriched_binary')

  HGM_ranked <- rbind(HGM_ranked, tmp)
}

for(sel.grp in groups){
  tmp <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    arrange(desc(-log10(fdr.p.value_enriched_rawcounts))) %>%
    dplyr::slice_max(order_by = -log10(fdr.p.value_enriched_rawcounts), n = 1000) %>%
    mutate(rank=1:nrow(.)) %>% mutate(FDR = fdr.p.value_enriched_rawcounts) %>% mutate(type = 'enriched_raw')

  HGM_ranked <- rbind(HGM_ranked, tmp)
}

for(sel.grp in groups){
  tmp <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    arrange(desc(-log10(fdr.p.value_depletion_binary))) %>%
    dplyr::slice_max(order_by = -log10(fdr.p.value_depletion_binary), n = 1000) %>%
    mutate(rank=1:nrow(.)) %>% mutate(FDR = fdr.p.value_depletion_binary) %>% mutate(type = 'depletion_binary')

  HGM_ranked <- rbind(HGM_ranked, tmp)
}

for(sel.grp in groups){
  tmp <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    arrange(desc(-log10(fdr.value_depletion_rawcounts))) %>%
    dplyr::slice_max(order_by = -log10(fdr.value_depletion_rawcounts), n = 1000) %>%
    mutate(rank=1:nrow(.)) %>% mutate(FDR = fdr.value_depletion_rawcounts) %>% mutate(type = 'depletion_raw')

  HGM_ranked <- rbind(HGM_ranked, tmp)
}


#make line plot
p.pval <- HGM_ranked %>%
  ggplot2::ggplot(aes(rank,-log10(FDR)))+
  ggplot2::geom_hline(yintercept = -log10(0.05),color='grey80',linetype=1) +
  ggplot2::geom_hline(yintercept = -log10(0.01),color='grey90',linetype=2) +
  ggplot2::geom_hline(yintercept = -log10(0.001),color='grey70',linetype=3) +
  ggplot2::geom_line(aes(color=group.x),size=1)+
  ggplot2::scale_color_manual(values=grid.colors)+
  ggplot2::xlim(c(0,1000))+
  ggplot2::xlab('feature rank')+
  ggplot2::ylab('-log10(FDR)')+
  facet_wrap(~type,ncol=2,scales='free_y')+
  fdb_style()

#save
ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/",pan.name,"_enriched_binary",'_overall_ranked',".pdf"),
                plot=p.pval,
                width = 7,
                height = 7,
                unit='in')



#------------------------------------------------------------------------#
#HGM_RF_COMB  %>%
#  ggplot2::ggplot(aes(z.score_enriched_rawcounts,-log10(fdr.p.value_enriched_rawcounts))) + geom_point(aes(color=group.x)) +
#  scale_color_manual(values=col.list[['sc0.85']])


#------------------------------------------------------------------------------------#
####  DETECT CO-LOCALISATION OF GENE SETS ####
####  detect and visualise ####
# for each genome in each group, we identify lineage-associated gene loci and look for groups of co-localised genes in on the chromosome (max distance of 3000 nucleotides, groupsize>=3)
# data is summarised per co-localised set of genes ('proximity cluster') and exported as output file
# proximity clusters are visualised per genome and annotated with singificance levels and Random Forest output
# we assessed the distribution of each proximity clusters in all genomes and visualise the output along a phylogenetic tree

listedRepresentatives <- genome.tbl %>%
  dplyr::left_join(sorted.cliques, by='genome') %>%
  dplyr::left_join(out.colors, by='genome') %>%
  dplyr::filter(!is.na(sc0.85)) %>%
  dplyr::select(genome, sc0.85, Completeness, quality_class, TypeStrain) %>%
  dplyr::mutate(sel.grp = paste0('cl',sc0.85), gnm = genome) %>%
  dplyr::select(sel.grp,gnm) %>%
  dplyr::arrange(sel.grp) %>%
  data.frame()

if(remove_small_groups == TRUE){
  listedRepresentatives <- listedRepresentatives %>% filter(sel.grp %in% larger_groups)
}



#--------------------------------------------------------------------#
#initialise output
prox.summary.out <- c()
prox.annotated.out <- c()

#i=1

# Big loop that runs over all representatives per group, looks for gene clusters, and project this on all marinobacter
for(i in 1:nrow(listedRepresentatives)){
  sel.grp = unname(listedRepresentatives[i,'sel.grp'])
  gnm = unname(listedRepresentatives[i,'gnm'])
  message('=== analysing ',gnm)

  #extract top 250 OGs ranked according to variable importance of RF
  RF.250 <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
    dplyr::slice(1:250) %>%
    dplyr::select(feature.id.x) %>%
    dplyr::pull() %>%
    as.character()

  #combine HGM and RF data
  out.mat <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
    dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(RFtop250 = ifelse(feature.id.x %in% RF.250,TRUE,FALSE))


  # select features that are enriched (HGM_FDR<.05) in either raw and binary counts
  selected.features <-
    out.mat %>%
    dplyr::filter(Enrich_raw ==TRUE | Enrich_bin ==TRUE)  %>%
    dplyr::select(feature.id.x) %>%
    dplyr::pull() %>%
    as.character()

  #report
  message(length(selected.features), ' selected features')

  #select locus_tags from reference genome corresponding to the selected OGs
  sel.loci <- annotation.tbl %>%
    dplyr::filter(genome %in% sel.genome) %>%
    dplyr::filter(OG %in% selected.features) %>%
    dplyr::filter(type=='CDS') %>% dplyr::filter(!is.na(OG))%>%
    dplyr::mutate(OG = coalesce(OG,locus_tag)) %>%
    dplyr::filter(genome==gnm) %>%
    dplyr::select(locus_tag) %>% pull()

  #of those enriched features, which ones are present?
  present.features <- annotation.tbl %>%
    dplyr::filter(genome %in% sel.genome) %>%
    dplyr::filter(OG %in% selected.features) %>%
    dplyr::filter(type=='CDS') %>% dplyr::filter(!is.na(OG))%>%
    dplyr::mutate(OG = coalesce(OG,locus_tag)) %>%
    dplyr::filter(genome==gnm) %>%
    dplyr::select(OG) %>% pull() %>% unique() %>% length()

  #analyse gene proximity
  prox.clust <- proximity_clusters(annotation.tbl, sel.loci,distance=3000 )


  if(!is.null(prox.clust)){
    # extract the largest clusters formed by >3 clusters
    largest.clusts_all <- prox.clust %>%
      dplyr::group_by() %>%
      dplyr::count(prox.cluster) %>%
      dplyr::filter(n>3) %>%
      dplyr::arrange(desc(n)) %>%
      dplyr::select(prox.cluster) %>%
      dplyr::pull()

    if(length(largest.clusts_all)>0){
      {# # plot them all?
      # for(sel.clust in largest.clusts_all){
      #   selectedOG <- prox.clust %>% filter(prox.cluster == sel.clust) %>% select(OG) %>% pull() %>% as.character()
      #
      #   clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
      #   clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
      #   clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
      #   clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)
      #
      #   #clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=FALSE, OGs=mapOG)
      #   mapOGs <- detect_sco_in_region(clusdetect.sorted)
      #   mapOG = mapOGs[1]
      #   if(is.null(mapOG)){
      #     message('consider alternative')
      #   }
      #
      #   # -----
      #   #Supply more than one OG to OGs to organise them around center (if not all have an SCO in this region)
      #   clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
      #
      #   #use function to add dummies
      #   clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
      #
      #   clusdetect.center$genome <- factor(clusdetect.center$genome,levels = tip.order)
      #
      #   focusOG=selectedOG
      #
      #   orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
      #   og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()
      #
      #
      #
      #   out.p <- clusdetect.center %>%
      #     mutate(direction=ifelse(strand=="+",1,-1)) %>%
      #     mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
      #     mutate(OG = factor(OG, levels=og.oder)) %>%
      #     ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=Name,fill=OG)) + #color=SELCOG
      #     gggenes::geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) +
      #     gggenes::geom_gene_label(align = "left")+
      #     viridis::scale_fill_viridis(discrete = TRUE)+
      #     #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
      #     theme(legend.position = "none")
      #
      #   out.p
      #
      #   ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/operons/prox.",gnm,'_',sel.clust,"_testing4",'clusterdistribution',".pdf"),
      #                   plot=out.p,
      #                   width = 7,
      #                   height = 14,
      #                   unit='in')
      # }
}
      #fill and recalculate proximity clusters
      filtered <- prox.clust %>% filter(prox.cluster %in% largest.clusts_all)
      filled <- fillgaps2(in.cluster=filtered)
      prox.clust2 <- proximity_clusters(annotation.tbl, filled$locus_tag,distance=3000 )

      #the alternative is to plot is same scale, start all the clusters at 0 position
      prox.clust.origin <- c()
      for(prx.c in unique(prox.clust2$prox.cluster)){
        t.start <- prox.clust2 %>% filter(prox.cluster == prx.c) %>% slice(1) %>% select(start) %>% pull()
        tmp <- prox.clust2 %>% filter(prox.cluster == prx.c) %>% mutate(start=start - t.start, end=end-t.start)
        prox.clust.origin = rbind(prox.clust.origin, tmp)
      }

      #summarise proximity clusters for output
      prox.summary <- prox.clust.origin %>%
        dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
        dplyr::mutate(enriched = ifelse(OG %in% selected.features,TRUE,FALSE)) %>%
        dplyr::left_join(annotation.q.per.OG %>% select(OG, annotation_quality), by='OG') %>%
        dplyr::group_by(genome, seqnames, prox.cluster) %>%
        dplyr::summarise(
                  level=sel.col,
                  group=sel.grp,
                  n_genes=n(),
                  mean_length = mean(width),
                  sd_length=sd(width),
                  cluster_length=max(end),
                  known_with_pfam = sum(annotation_quality=='known with pfam'),
                  known_without_pfam = sum(annotation_quality=='known without pfam'),
                  unknown_with_pfam = sum(annotation_quality=='unknown with pfam'),
                  unknown_without_pfam = sum(annotation_quality=='unknown without pfam'),
                  #n_hypothetical = sum(grepl('hypothetical',product)),
                  #nopfam=sum(is.na(pfam_acc)),
                  #n_poor = sum((grepl('hypothetical',product)& is.na(pfam_acc))),
                  #perc_hypotheical = n_hypothetical/n_genes,
                  n_transposase = sum(grepl('transposase',product)),
                  n_enriched = sum(enriched),
                  perc_of_total = n_enriched/length(sel.loci),
                  total_enriched = length(sel.loci),
        ) %>%
        dplyr::ungroup()

      #save summary
      write.table(prox.summary,
                  file=paste0("~/DATA/MarbGenomics/Graphs/operons2/prox.summary.",gnm,'_',sel.clust,"_final.txt"),
                  quote=FALSE,
                  sep='\t',
                  row.names = FALSE)

      #KO_annotated_clean <- KO_annotated %>% unique() %>%
      #  dplyr::filter(genome==gnm) %>%
      #  dplyr::group_by(locus_tag) %>%
      #  dplyr::mutate(KO=as.character(KO),
      #         ko_annotation=as.character(ko_annotation),
      #         MOD_F = as.character(MOD_F),
      #         TYPE=as.character(TYPE),
      #         DESCRIPTION=as.character(DESCRIPTION)) %>%
      #  dplyr::summarise(KO=paste0(unique(KO),collapse=';'),
      #            ko_annotation=paste0(unique(ko_annotation),collapse=';'),
      #            MOD_F=paste0(unique(MOD_F),collapse=';'),
      #            TYPE=paste0(unique(TYPE),collapse=';'),
      #            DESCRIPTION=paste0(unique(DESCRIPTION),collapse=';')) %>% ungroup()

      prox.annotated <- prox.clust2 %>%
        dplyr::left_join(ko2ann, by=c('kfm_domain'='ko')) %>%
        dplyr::left_join(cog_func_single, by='COG') %>%
        dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>%
        dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
        dplyr::select(-score,-note,-rpt_family,-rpt_type,-rpt_unit_seq)

      if(!is.null(prokka2ncbi.annotated.list[[gnm]])){
        message('adding NCBI accessions and annotations')
        prox.annotated <- prox.annotated %>% left_join(prokka2ncbi.annotated.list[[gnm]], by=c('locus_tag'='prokka_locus_tag'))
      }else{
        prox.annotated <- prox.annotated %>% mutate(
          protein_id=NA,
          evalue=NA,
          non.redundant_refseq=NA,
          ncbi_product=NA,
          ncbi_symbol=NA,
          ncbi_GeneID=NA,
          ncbi_locus_tag=NA,
          ncbi_old_locus_tag=NA
        )
        }


      #save summary
      write.table(prox.annotated,
                  file=paste0("~/DATA/MarbGenomics/Graphs/operons2/prox.annotated.",gnm,'_',sel.clust,"_final.txt"),
                  quote=FALSE,
                  sep='\t',
                  row.names = FALSE)

      v.cls <- prox.clust.origin %>%
        dplyr::select(prox.cluster) %>%
        dplyr::pull() %>%
        unique()

      # Make blocks of 10 gene clusters for plotting purposes,
      plotlist <- list()
      plotlist2 <- list()
      io <- 0

      #s=split(v.cls, ceiling(seq_along(v.cls)/10))[[1]]
      for(s in split(v.cls, ceiling(seq_along(v.cls)/10))){
        io<-io+1
        largest.clusts <- s
        f.clust.origin <- prox.clust.origin %>% dplyr::filter(prox.cluster %in% largest.clusts)
        n_in_clusters <- f.clust.origin %>% nrow()

        out.p <-  f.clust.origin %>%
          #prox.clust.origin %>%
          dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
          dplyr::mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
          dplyr::mutate(sel.loci = ifelse(locus_tag %in% sel.loci,1,0.5)) %>%
          dplyr::left_join(out.mat,by=c('OG'='feature.id.x'))%>%
          dplyr::mutate(Enrich_raw = ifelse(Enrich_raw==TRUE,TRUE,NA)) %>%
          dplyr::mutate(Enrich_bin = ifelse(Enrich_bin==TRUE,TRUE,NA)) %>%
          dplyr::mutate(middle=(start+end)/2) %>%
          ggplot2::ggplot(aes(x=middle, xmin = start, xmax = end, y = seqnames,forward=direction,label=Name)) +
          gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),aes(fill=transposase)) +
          gggenes::geom_gene_label(align = "left")+
          geom_text(data=.%>%filter(!is.na(Name)), aes(label=locus_tag),angle=45,hjust=-.2,size=2)+
          ggplot2::scale_fill_manual(values=c('grey83','grey20'))+
          ggplot2::theme(legend.position = "none")+
          ggplot2::facet_wrap(~ prox.cluster, ncol = 1, scales='free_y')+
          {if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_bin == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_bin == TRUE),aes(x = start, y = '0binary', xend = end, yend = '0binary', colour = -log10(fdr.p.value_enriched_binary)),size=1.5)} +
          {if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_raw == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_raw == TRUE), aes(x = start, y = '00raw', xend = end, yend = '00raw', colour = log10(fdr.p.value_enriched_rawcounts)),size=1.5)} +
          scale_color_gradientn(colours=c(adjustcolor('dodgerblue',alpha=0.1),'dodgerblue')) +
          ggnewscale::new_scale_color()+#
          scale_alpha(range=c(0.5,1))+ ggtitle(paste0(gnm,'')) +
          geom_segment(aes(x = start, y = '0RF', xend = end, yend = '0RF', colour = log10(MeanDecreaseAccuracy)),size=1.5) +
          scale_color_gradientn(colours=c(adjustcolor('purple',alpha=0.1),'purple')) +
          fdb_style(aspect.ratio=0.5) +
          ggplot2::theme_grey() +
          ggplot2::theme(panel.background = ggplot2::element_blank(),
                         panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 0.5),
                         panel.grid.minor.y = ggplot2::element_blank(),
                         panel.grid.minor.x = ggplot2::element_blank(),
                         panel.grid.major.x = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
                         axis.ticks.x = ggplot2::element_line(colour = "grey20", size=0.5),
                         strip.text = ggplot2::element_text(size  = 9,  hjust = 0),
                         strip.background = ggplot2::element_blank()
          )
        plotlist2[[io]] <- out.p
      }

      if(length(plotlist2)>0){
        for(p in 1:length(plotlist2)){
          message('saving ',p)
          ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.tight_o_",sel.col,'_',sel.grp,'_',gnm,'_',p,'_proximity_clusters2',".pdf"),
                          plot=plotlist2[[p]],
                          width = 12,
                          height = 8,
                          unit='in')
        }
      }


      #storing in bigger objects
      prox.summary.out <- rbind(prox.summary.out,prox.summary)
      prox.annotated.out <- rbind(prox.annotated.out,prox.annotated)


      #-------------------------------------#
       #plot along phylogeny
  #     largest.clusts_all <- prox.clust2 %>% select(prox.cluster) %>% unique() %>% pull
  #     cluster.out.list <-list()
  #     for(sel.clust in largest.clusts_all){
  #       selectedOG <- prox.clust2 %>% dplyr::filter(prox.cluster == sel.clust) %>% dplyr::filter(!is.na(OG)) %>% dplyr::select(OG) %>% pull() %>% as.character()
  #       if(length(selectedOG)>0){
  #         clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
  #         clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
  #         clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
  #         clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)
  #
  #         mapOGs <- detect_sco_in_region(clusdetect.sorted)
  #         mapOG = mapOGs[1]
  #         if(is.null(mapOG)){
  #           message('consider alternative')
  #         }
  #
  #         clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
  #         clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
  #         clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
  #         focusOG=selectedOG
  #
  #         orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
  #         og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()
  #
  #         clusdetect.center <- clusdetect.center %>%
  #           dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
  #           dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
  #           dplyr::mutate(OG = factor(OG, levels=og.oder))
  #
  #         cluster.out.list[[sel.clust]] <- clusdetect.center
  #       }
  #     }
  #
  #     v.cls <- prox.clust.origin %>%
  #       select(prox.cluster) %>% pull() %>% unique()
  #
  #     o=0
  #     for(s in split(1:length(largest.clusts_all), ceiling(seq_along(1:length(largest.clusts_all))/10))){
  #       o=o+1
  #       v.options <- rep(c('viridis','plasma','magma','cividis','inferno','mako'),10)
  #       p = ggtree::ggtree(sco.106.concat.tree) + ggtree::geom_tiplab(align=TRUE)
  #
  #       cluster.out.s <- cluster.out.list[s]
  #
  #       for(i in 1:length(names(cluster.out.s))){
  #         message(v.options[i])
  #         p <-    ggtree::facet_plot(p,
  #                            panel=names(cluster.out.s)[i],
  #                            mapping = aes(x = start, xend = end, yend=y,color=OG),
  #                            data=cluster.out.s[[i]] %>%
  #                              mutate(id=genome) %>%
  #                              select(id, everything()),
  #                            geom=geom_segment,size=1.5)+
  #           theme(legend.position = "none")+
  #           viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
  #           ggnewscale::new_scale_color()
  #        }
  #
  #       ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.phylo_22_26_.",sel.col,'_',sel.grp,'_',gnm,'_',o,'_proximity_clusters',".pdf"),
  #                       plot=p,
  #                       width = 12,
  #                       height = 8,
  #                       unit='in')
  #     }
     }
   }
}





#====================================================================#
prox.s.annotated <- read.delim(file='/Users/sa01fd/DATA/MarbGenomics/Graphs/operons2/prox.annotated.Marinobacter_algicola_DG893_cluster3.txt')
prox.s.annotated <- read.delim(file='/Users/sa01fd/DATA/MarbGenomics/Graphs/operons2/prox.annotated.Marinobacter_adhaerens_HP15_cluster3.txt')
prox.s.annotated <- read.delim(file='/Users/sa01fd/DATA/MarbGenomics/Graphs/operons2/prox.annotated.Marinobacter_hydrocarbonoclasticus_ATCC_49840_cluster3.txt')
prox.s.annotated <- read.delim(file='/Users/sa01fd/DATA/MarbGenomics/Graphs/operons2/prox.annotated.Marinobacter_psychrophilus_20041_cluster3.txt')
prox.s.annotated <- read.delim(file='/Users/sa01fd/DATA/MarbGenomics/Graphs/operons2/prox.annotated.Marinobacter_sp_LV10R510_8_cluster3.txt')



#select clusters
prox.clust.origin <- c()
for(prx.c in unique(prox.s.annotated$prox.cluster)){
  t.start <- prox.s.annotated %>% filter(prox.cluster == prx.c) %>% slice(1) %>% select(start) %>% pull()
  tmp <- prox.s.annotated %>% filter(prox.cluster == prx.c) %>% mutate(start=start - t.start, end=end-t.start)
  prox.clust.origin = rbind(prox.clust.origin, tmp)
}

v.cls <- prox.clust.origin %>% select(prox.cluster) %>% unique() %>% pull() %>% as.character()
s <- split(v.cls, ceiling(seq_along(v.cls)/10))[[1]]

#algicola
s <- c('cluster3','cluster6','cluster22','cluster23','cluster26','cluster27','cluster45')

#adhaerens
s <- c('cluster13','cluster10','cluster18','cluster29','cluster27')

#ATCC
s <- c('cluster14','cluster17','cluster19','cluster5','cluster4','cluster22')

#20041
s <- c('cluster21','cluster15','cluster11','cluster1','cluster3','cluster19')

#LVR
s <- c('cluster10','cluster15','cluster16','cluster26','cluster28','cluster37')



for(s in split(v.cls, ceiling(seq_along(v.cls)/10))){
  io<-io+1
  largest.clusts <- s
  f.clust.origin <- prox.clust.origin %>% dplyr::filter(prox.cluster %in% largest.clusts)
  n_in_clusters <- f.clust.origin %>% nrow()

  out.p <-  f.clust.origin %>%
    dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
    dplyr::mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
    dplyr::mutate(sel.loci = ifelse(locus_tag %in% sel.loci,1,0.5)) %>%
    dplyr::mutate(Enrich_raw = ifelse(Enrich_raw==TRUE,TRUE,NA)) %>%
    dplyr::mutate(Enrich_bin = ifelse(Enrich_bin==TRUE,TRUE,NA)) %>%
    dplyr::mutate(middle=(start+end)/2) %>%
    dplyr::mutate(ko_short = gsub("\\;.*","",ko_annotation)) %>%
    ggplot2::ggplot(aes(x=middle, xmin = start, xmax = end, y = seqnames,forward=direction,label=ncbi_old_locus_tag)) +
    gggenes::geom_gene_arrow(arrowhead_width = unit(2, "mm"),aes(fill=transposase)) +
    gggenes::geom_gene_label(align = "left")+
    ggplot2::geom_text(data=.%>%filter(!is.na(ko_short)), aes(label=ko_short,x=middle),angle=90,hjust=-.5,size=3)+
    ggplot2::scale_fill_manual(values=c('grey83','grey20'))+
    ggplot2::theme(legend.position = "none")+
    ggplot2::facet_wrap(~ prox.cluster,
                        ncol = 1,
                        scales='free_y',
                        strip.position="right")+
    {if( f.clust.origin %>% filter(Enrich_bin == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_bin == TRUE),aes(x = start, y = '0binary', xend = end, yend = '0binary', colour = -log10(fdr.p.value_enriched_binary)),size=1.5)} +
    {if( f.clust.origin %>% filter(Enrich_raw == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_raw == TRUE), aes(x = start, y = '00raw', xend = end, yend = '00raw', colour = log10(fdr.p.value_enriched_rawcounts)),size=1.5)} +
    scale_color_gradientn(colours=c(adjustcolor('forestgreen',alpha=0.1),'forestgreen')) +
    ggnewscale::new_scale_color()+#
    scale_alpha(range=c(0.5,1))+
    ggtitle(paste0(gnm,'')) +
    geom_segment(aes(x = start, y = '0RF', xend = end, yend = '0RF', colour = log10(MeanDecreaseAccuracy)),size=1.5) +
    scale_color_gradientn(colours=c(adjustcolor('purple',alpha=0.1),'purple'),na.value = 'white') +
    fdb_style(aspect.ratio=0.5) +
    ggplot2::theme_grey() +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 0.5),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
                   axis.ticks.x = ggplot2::element_line(colour = "grey20", size=0.5),
                   strip.text = ggplot2::element_text(size  = 9,  hjust = 0),
                   strip.background = ggplot2::element_blank()
    )
  plotlist2[[io]] <- out.p
}


ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.LVR5_selection2",".pdf"),
                plot=out.p,
                width = 12,
                height = 8,
                unit='in')







annotation.tbl %>%
  dplyr::filter(locus_tag %in% prokka_loci) %>%
  dplyr::left_join(annotation.quality.tbl %>%
                     select(locus_tag, ncbi_old_locus_tag,ko_annotation, pfam_acc, pfam_name, n_known_pfam, annotation_quality),
                   by='locus_tag') %>%
  mutate(ko_short = gsub("\\;.*","",ko_annotation)) %>%

  #prox.clust.origin %>%
  dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
  dplyr::mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
  dplyr::mutate(sel.loci = ifelse(locus_tag %in% sel.loci,1,0.5)) %>%
  dplyr::left_join(out.mat,by=c('OG'='feature.id.x'))%>%
  dplyr::mutate(Enrich_raw = ifelse(Enrich_raw==TRUE,TRUE,NA)) %>%
  dplyr::mutate(Enrich_bin = ifelse(Enrich_bin==TRUE,TRUE,NA)) %>%
  dplyr::mutate(middle=(start+end)/2) %>%
  ggplot2::ggplot(aes(x=middle, xmin = start, xmax = end, y = seqnames,forward=direction,label=ncbi_old_locus_tag)) +
  gggenes::geom_gene_arrow(#arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(2, "mm"),aes(fill=annotation_quality)) +
  #gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),aes(fill=transposase)) +
  gggenes::geom_gene_label(align = "left")+
  geom_text(data=.%>%filter(!is.na(ko_short)), aes(label=ko_short),angle=45,hjust=-.5)+
  #ggplot2::scale_fill_manual(values=c('grey83','grey20'))+
  ggplot2::theme(legend.position = "none")+
  #ggplot2::facet_wrap(~ prox.cluster, ncol = 1, scales='free_y')+
  #{if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_bin == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_bin == TRUE),aes(x = start, y = '0binary', xend = end, yend = '0binary', colour = -log10(fdr.p.value_enriched_binary)),size=1.5)} +
  #{if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_raw == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_raw == TRUE), aes(x = start, y = '00raw', xend = end, yend = '00raw', colour = log10(fdr.p.value_enriched_rawcounts)),size=1.5)} +
  scale_color_gradientn(colours=c(adjustcolor('dodgerblue',alpha=0.1),'dodgerblue')) +
  ggnewscale::new_scale_color()+#
  scale_alpha(range=c(0.5,1))+ ggtitle(paste0(ref.genome,'')) +
  #geom_segment(aes(x = start, y = '0RF', xend = end, yend = '0RF', colour = log10(MeanDecreaseAccuracy)),size=1.5) +
  scale_color_gradientn(colours=c(adjustcolor('purple',alpha=0.1),'purple')) +
  fdb_style(aspect.ratio=0.5) +
  ggplot2::theme_grey() +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 0.5),
                 panel.grid.minor.y = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
                 axis.ticks.x = ggplot2::element_line(colour = "grey20", size=0.5),
                 strip.text = ggplot2::element_text(size  = 9,  hjust = 0),
                 strip.background = ggplot2::element_blank()
  )



#=======================================================================#


prox.summary.out %>% arrange(desc(perc_hypotheical)) %>% ggplot(aes(n_genes,cluster_length))+geom_point(shape=21,aes(color=perc_hypotheical)) + fdb_style()

prox.summary.out %>% arrange(desc(n_genes)) %>% ggplot(aes(n_genes))+geom_histogram() + fdb_style()

prox.summary.out %>%
  arrange(desc(n_genes)) %>%
  ggplot(aes(x=1:nrow(prox.summary.out),y=n_genes)) +
  geom_point(shape=21,aes(color=perc_hypotheical))+xlab('clusters sorted by ngenes')+ylab('n genes in cluster')








gnm <- 'Marinobacter_hydrocarbonoclasticus_ATCC_49840'
gnm <- 'Marinobacter_adhaerens_HP15'

sel.clusters <- paste0('cluster',5:11)
#c('cluster17','cluster6','cluster14','cluster4','cluster5')
 #c('cluster1','cluster2','cluster3')

prox.s <- prox.annotated.out %>%
  filter(genome==gnm)%>%
  #filter(prox.cluster %in% sel.clusters)%>%
  dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
  dplyr::mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
  dplyr::mutate(sel.loci = ifelse(locus_tag %in% sel.loci,1,0.5))

prox.clust.origin <- c()
for(prx.c in unique(prox.s$prox.cluster)){
  t.start <- prox.s %>% filter(prox.cluster == prx.c) %>% slice(1) %>% select(start) %>% pull()
  tmp <- prox.s %>% filter(prox.cluster == prx.c) %>% mutate(start=start - t.start, end=end-t.start)
  prox.clust.origin = rbind(prox.clust.origin, tmp)
}


p <- prox.clust.origin %>%
  dplyr::mutate(middle=(start+end)/2) %>%
  ggplot2::ggplot(aes(x=middle, xmin = start, xmax = end, y = prox.cluster,forward=direction,label=ncbi_old_locus_tag)) +
  gggenes::geom_gene_arrow(#arrowhead_height = unit(3, "mm"),
                           arrowhead_width = unit(2, "mm"),
                           aes(fill=-log10(fdr.p.value_enriched_binary))) +
  gggenes::geom_gene_label(align = "left",color='white')+
  scale_fill_gradient(low='white',high="#AE6331") + fdb_style(aspect.ratio=1)+
  geom_text(data=.%>%filter(!is.na(Name)), aes(label=Name),angle=45,hjust=-.9,size=2)+
  ggplot2::theme(legend.position = "none")+
  ggplot2::facet_wrap(~ seqnames,
                      ncol = 1,
                      scales ='free_y',
                      strip.position = "right")

ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/test_gen_aesthetics",".pdf"),
                       plot=p,
                       width = 8,
                       height = 4,
                       unit='in')




#  {if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_bin == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_bin == TRUE),aes(x = start, y = '0binary', xend = end, yend = '0binary', colour = -log10(fdr.p.value_enriched_binary)),size=1.5)} +
#  {if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_raw == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_raw == TRUE), aes(x = start, y = '00raw', xend = end, yend = '00raw', colour = log10(fdr.p.value_enriched_rawcounts)),size=1.5)} +
#  scale_color_gradientn(colours=c(adjustcolor('dodgerblue',alpha=0.1),'dodgerblue')) +
#  ggnewscale::new_scale_color()+#
#  scale_alpha(range=c(0.5,1))+ ggtitle(paste0(gnm,'')) +
#  geom_segment(aes(x = start, y = '0RF', xend = end, yend = '0RF', colour = log10(MeanDecreaseAccuracy)),size=1.5) +
#  scale_color_gradientn(colours=c(adjustcolor('purple',alpha=0.1),'purple')) +
#  fdb_style(aspect.ratio=0.5) +
#  ggplot2::theme_grey() +
#  ggplot2::theme(panel.background = ggplot2::element_blank(),
#                 panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 0.5),
#                 panel.grid.minor.y = ggplot2::element_blank(),
#                 panel.grid.minor.x = ggplot2::element_blank(),
#                 panel.grid.major.x = ggplot2::element_blank(),
#                 axis.ticks.y = ggplot2::element_blank(),
#                 axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
#                 axis.ticks.x = ggplot2::element_line(colour = "grey20", size=0.5),
#                 strip.text = ggplot2::element_text(size  = 9,  hjust = 0),
#                 strip.background = ggplot2::element_blank()
#  )







prox.tight_o_sc0.85_cl3_Marinobacter_hydrocarbonoclasticus_ATCC_49840_2_proximity_clusters2

#==============================================================#
### clustering in softcore ####

softcore.ogs <- og.table %>%
  dplyr::mutate(softcore = ifelse(nrOfSpecies>=94*0.95,TRUE,FALSE)) %>%
  dplyr::filter(softcore==TRUE) %>%
  dplyr::select(feature.id) %>%
  dplyr::pull()


#select locus_tags from reference genome corresponding to the selected OGs
gnm <- 'Marinobacter_adhaerens_HP15'
gnm <- 'Marinobacter_hydrocarbonoclasticus_ATCC_49840'


sel.loci <- annotation.tbl %>%
  dplyr::filter(OG %in% softcore.ogs) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter(!is.na(OG))%>%
  dplyr::mutate(OG = coalesce(OG,locus_tag)) %>%
  dplyr::filter(genome==gnm) %>%
  dplyr::select(locus_tag) %>% pull()


#analyse gene proximity
prox.clust <- proximity_clusters(annotation.tbl, sel.loci,distance=3000 )


if(!is.null(prox.clust)){
  # extract the largest clusters formed by >3 clusters
  largest.clusts_all <- prox.clust %>%
    dplyr::group_by() %>%
    dplyr::count(prox.cluster) %>%
    dplyr::filter(n>3) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::select(prox.cluster) %>%
    dplyr::pull()

  if(length(largest.clusts_all)>0){
    #fill and recalculate proximity clusters
    filtered <- prox.clust %>% dplyr::filter(prox.cluster %in% largest.clusts_all)
    filled <- fillgaps2(in.cluster=filtered)
    prox.clust2 <- proximity_clusters(annotation.tbl, filled$locus_tag,distance=3000 )

    #the alternative is to plot is same scale, start all the clusters at 0 position
    prox.clust.origin <- c()
    for(prx.c in unique(prox.clust2$prox.cluster)){
      t.start <- prox.clust2 %>%
        dplyr::filter(prox.cluster == prx.c) %>%
        dplyr::dplyr::slice(1) %>%
        dplyr::select(start) %>%
        dplyr::pull()
      tmp <- prox.clust2 %>% dplyr::filter(prox.cluster == prx.c) %>% dplyr::mutate(start=start - t.start, end=end-t.start)
      prox.clust.origin = rbind(prox.clust.origin, tmp)
    }

    #summarise proximity clusters for output
    prox.summary <- prox.clust.origin %>%
      dplyr::left_join(df.pfam_per_locus %>%
                         dplyr::select(-genome), by=c('locus_tag'='seq_id')) %>%
      dplyr::mutate(enriched = ifelse(OG %in% selected.features,TRUE,FALSE)) %>%
      dplyr::group_by(genome, seqnames, prox.cluster) %>%
      dplyr::summarise(
        level=sel.col,
        group=sel.grp,
        n_genes=n(),
        mean_length = mean(width),
        sd_length=sd(width),
        cluster_length=max(end),
        n_hypothetical = sum(grepl('hypothetical',product)),
        nopfam=sum(is.na(pfam_acc)),
        n_poor = sum((grepl('hypothetical',product)& is.na(pfam_acc))),
        perc_hypotheical = n_hypothetical/n_genes,
        n_transposase = sum(grepl('transposase',product)),
        n_enriched = sum(enriched),
        perc_of_total = n_enriched/length(sel.loci),
        total_enriched = length(sel.loci),
      ) %>%
      dplyr::ungroup()

    #save summary
    write.table(prox.summary,
                file=paste0("~/DATA/MarbGenomics/Graphs/operons/soft_core.prox.summary.",gnm,'_',sel.clust,".txt"),
                quote=FALSE,
                sep='\t',
                row.names = FALSE)

    prox.annotated <- prox.clust2 %>%
      dplyr::left_join(ko2ann, by=c('kfm_domain'='ko')) %>%
      dplyr::left_join(cog_func_single, by='COG') %>%
      dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>%
      dplyr::left_join(df.pfam_per_locus %>%
                         dplyr::select(-genome), by=c('locus_tag'='seq_id')) %>%
      dplyr::select(-score,-note,-rpt_family,-rpt_type,-rpt_unit_seq)


    if(!is.null(prokka2ncbi.annotated.list[[gnm]])){
      message('adding NCBI accessions and annotations')
      prox.annotated <- prox.annotated %>%
        dplyr::left_join(prokka2ncbi.annotated.list[[gnm]], by=c('locus_tag'='prokka_locus_tag'))
    }


    #save summary
    write.table(prox.annotated,
                file=paste0("~/DATA/MarbGenomics/Graphs/operons/soft_core.prox.annotated.",gnm,'_',sel.clust,".txt"),
                quote=FALSE,
                sep='\t',
                row.names = FALSE)

    v.cls <- prox.clust.origin %>%
      dplyr::select(prox.cluster) %>%
      dplyr::pull() %>%
      unique()

    # Make blocks of 10 gene clusters for plotting purposes,
    plotlist <- list()
    plotlist2 <- list()
    io <- 0

    for(s in split(v.cls, ceiling(seq_along(v.cls)/10))){
      io<-io+1
      largest.clusts <- s
      f.clust.origin <- prox.clust.origin %>% dplyr::filter(prox.cluster %in% largest.clusts)
      n_in_clusters <- f.clust.origin %>% nrow()

      out.p <-  f.clust.origin %>%
        dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
        dplyr::mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
        dplyr::mutate(sel.loci = ifelse(locus_tag %in% sel.loci,1,0.5)) %>%
        dplyr::left_join(out.mat,by=c('OG'='feature.id.x'))%>%
        dplyr::mutate(Enrich_raw = ifelse(Enrich_raw==TRUE,TRUE,NA)) %>%
        dplyr::mutate(Enrich_bin = ifelse(Enrich_bin==TRUE,TRUE,NA)) %>%
        dplyr::mutate(middle=(start+end)/2) %>%
        ggplot2::ggplot(aes(x=middle, xmin = start, xmax = end, y = seqnames,forward=direction,label=Name)) +
        gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),aes(fill=transposase)) +
        gggenes::geom_gene_label(align = "left")+
        geom_text(data=.%>%filter(!is.na(Name)), aes(label=locus_tag),angle=45,hjust=-.2,size=2)+
        ggplot2::scale_fill_manual(values=c('grey83','grey20'))+
        ggplot2::theme(legend.position = "none")+
        ggplot2::facet_wrap(~ prox.cluster, ncol = 1, scales='free_y')+
        #{if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_bin == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_bin == TRUE),aes(x = start, y = '0binary', xend = end, yend = '0binary', colour = -log10(fdr.p.value_enriched_binary)),size=1.5)} +
        #{if( f.clust.origin %>% dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>% filter(Enrich_raw == TRUE) %>% nrow > 0) geom_segment(data= . %>% filter(Enrich_raw == TRUE), aes(x = start, y = '00raw', xend = end, yend = '00raw', colour = log10(fdr.p.value_enriched_rawcounts)),size=1.5)} +
        scale_color_gradientn(colours=c(adjustcolor('dodgerblue',alpha=0.1),'dodgerblue')) +
        ggnewscale::new_scale_color()+#
        scale_alpha(range=c(0.5,1))+ ggtitle(paste0(gnm,'')) +
        #geom_segment(aes(x = start, y = '0RF', xend = end, yend = '0RF', colour = log10(MeanDecreaseAccuracy)),size=1.5) +
        scale_color_gradientn(colours=c(adjustcolor('purple',alpha=0.1),'purple')) +
        fdb_style(aspect.ratio=0.5) +
        ggplot2::theme_grey() +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 0.5),
                       panel.grid.minor.y = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
                       axis.ticks.x = ggplot2::element_line(colour = "grey20", size=0.5),
                       strip.text = ggplot2::element_text(size  = 9,  hjust = 0),
                       strip.background = ggplot2::element_blank()
        )
      plotlist2[[io]] <- out.p
    }

    if(length(plotlist2)>0){
      for(p in 1:length(plotlist2)){
        message('saving ',p)
        ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/softcore.tight_o_",sel.col,'_',sel.grp,'_',gnm,'_',p,'_proximity_clusters',".pdf"),
                        plot=plotlist2[[p]],
                        width = 12,
                        height = 8,
                        unit='in')
      }
    }
    #prox.summary.out <- rbind(prox.summary.out,prox.summary)



    #-------------------------------------#
    #plot along phylogeny
        largest.clusts_all <- prox.clust2 %>% select(prox.cluster) %>% unique() %>% pull
        cluster.out.list <-list()
        for(sel.clust in largest.clusts_all){
          selectedOG <- prox.clust2 %>% dplyr::filter(prox.cluster == sel.clust) %>% dplyr::filter(!is.na(OG)) %>% dplyr::select(OG) %>% pull() %>% as.character()
          if(length(selectedOG)>0){
            clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
            clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
            clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
            clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)

            mapOGs <- detect_sco_in_region(clusdetect.sorted)
            mapOG = mapOGs[1]
            if(is.null(mapOG)){
              message('consider alternative')
            }

            clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
            clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
            clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
            focusOG=selectedOG

            orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
            og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()

            clusdetect.center <- clusdetect.center %>%
              dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
              dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
              dplyr::mutate(OG = factor(OG, levels=og.oder))

            cluster.out.list[[sel.clust]] <- clusdetect.center
          }
        }

        v.cls <- prox.clust.origin %>%
          select(prox.cluster) %>% pull() %>% unique()

        o=0
        for(s in split(1:length(largest.clusts_all), ceiling(seq_along(1:length(largest.clusts_all))/10))){
          o=o+1
          v.options <- rep(c('viridis','plasma','magma','cividis','inferno','mako'),10)
          p = ggtree::ggtree(sco.106.concat.tree) + ggtree::geom_tiplab(align=TRUE)

          cluster.out.s <- cluster.out.list[s]

          for(i in 1:length(names(cluster.out.s))){
            message(v.options[i])
            p <-    ggtree::facet_plot(p,
                               panel=names(cluster.out.s)[i],
                               mapping = aes(x = start, xend = end, yend=y,color=OG),
                               data=cluster.out.s[[i]] %>%
                                 mutate(id=genome) %>%
                                 select(id, everything()),
                               geom=geom_segment,size=1.5)+
              theme(legend.position = "none")+
              viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
              ggnewscale::new_scale_color()
           }

          ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/softcore.phylo_22_26_.",sel.col,'_',sel.grp,'_',gnm,'_',o,'_proximity_clusters',".pdf"),
                          plot=p,
                          width = 12,
                          height = 8,
                          unit='in')
        }


  }
}





#=================================================================================#
#### FIGURE ####
# # DSTRIBUTION OF CLUSTERS
#=================================================================================#

sel.genome <- sco.106.concat.tree$tip.label

prokka_loci <- prokka2ncbi.annotated.list[[ref.genome]] %>%
  dplyr::filter(grepl(paste(publ_loci,collapse='|'),ncbi_old_locus_tag)) %>%
  dplyr::select(prokka_locus_tag) %>%
  dplyr::pull() %>% as.character()


associated_OGs <- annotation.tbl %>%
  filter(locus_tag %in% prokka_loci) %>%
  filter(!is.na(OG)) %>%
  select(OG) %>% pull()
selectedOG <- associated_OGs



prokka2ncbi.tbl <- tibble(genome=names(prokka2ncbi.annotated.list),prokka2ncbi.annotated.list) %>% unnest(prokka2ncbi.annotated.list)

compiled.clusters.df <- read.delim('~/DATA/MarbGenomics/grouped_loci.txt')


#select clusters >1 lenght
compiled.clusters <- compiled.clusters.df %>%
  group_by(Group) %>%
  tally() %>%
  filter(n>1) %>%
  select(Group) %>%
  pull() %>%
  as.character()

c.c <- compiled.clusters[46]
c.c <- 'Benzoate_degradation'

all.out.compiled <- c()
for(c.c in compiled.clusters){
  publ_loci <- compiled.clusters.df %>% dplyr::filter(Group == c.c) %>% select(Gene_names) %>% pull() %>% as.character()
  publ_loci <- gsub(' ','',publ_loci)
  ncbi.selected <- prokka2ncbi.tbl %>%
    filter(grepl(paste(publ_loci,collapse='|'),ncbi_old_locus_tag))

  prokka_loci <- ncbi.selected %>%
    select(prokka_locus_tag) %>% pull() %>% as.character()



  ref.genome <- ncbi.selected %>% select(genome) %>% unique() %>% pull() %>% as.character()

  #nice vis--#
  annotation.tbl %>%
    dplyr::filter(locus_tag %in% prokka_loci) %>%
    dplyr::left_join(annotation.quality.tbl %>%
                       select(locus_tag, ncbi_old_locus_tag,ko_annotation, pfam_acc, pfam_name, n_known_pfam, annotation_quality),
                     by='locus_tag') %>%
    mutate(ko_short = gsub("\\;.*","",ko_annotation)) %>%

    #prox.clust.origin %>%
    dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
    dplyr::mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
    dplyr::mutate(sel.loci = ifelse(locus_tag %in% sel.loci,1,0.5)) %>%
    dplyr::left_join(out.mat,by=c('OG'='feature.id.x'))%>%
    dplyr::mutate(Enrich_raw = ifelse(Enrich_raw==TRUE,TRUE,NA)) %>%
    dplyr::mutate(Enrich_bin = ifelse(Enrich_bin==TRUE,TRUE,NA)) %>%
    dplyr::mutate(middle=(start+end)/2) %>%
    ggplot2::ggplot(aes(x=middle, xmin = start, xmax = end, y = seqnames,forward=direction,label=ncbi_old_locus_tag)) +
    gggenes::geom_gene_arrow(arrowhead_width = unit(2, "mm"),
                             aes(fill=annotation_quality)) +
    gggenes::geom_gene_label(align = "left")+
    geom_text(data=.%>%filter(!is.na(ko_short)),
              aes(label=ko_short),
              angle=45,hjust=-.5)+
    ggplot2::theme(legend.position = "none")+
    scale_color_gradientn(colours=c(adjustcolor('dodgerblue',alpha=0.1),'dodgerblue')) +
    ggnewscale::new_scale_color()+#
    scale_alpha(range=c(0.5,1))+ ggtitle(paste0(ref.genome,'')) +
    scale_color_gradientn(colours=c(adjustcolor('purple',alpha=0.1),'purple')) +
    fdb_style(aspect.ratio=0.5) +
    ggplot2::theme_grey() +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 0.5),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
                   axis.ticks.x = ggplot2::element_line(colour = "grey20", size=0.5),
                   strip.text = ggplot2::element_text(size  = 9,  hjust = 0),
                   strip.background = ggplot2::element_blank()
    )

  #-----------#

  associated_OGs <- annotation.tbl %>%
    filter(locus_tag %in% prokka_loci) %>%
    filter(!is.na(OG)) %>%
    select(OG) %>% pull()

  selectedOG <- associated_OGs


  if(length(selectedOG) > 1){

    clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
    clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
    clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
    clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)
    #

    if(!is.null(clusdetect.sorted)){
      mapOGs <- detect_sco_in_region(clusdetect.sorted)
      mapOG = mapOGs[1]
      if(is.null(mapOG)){
        message('consider alternative')
      }
      #
      clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
      clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
      clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
      focusOG=selectedOG
      #
      orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
      og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()
      #
      clusdetect.center <- clusdetect.center %>%
        dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
        dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
        dplyr::mutate(OG = factor(OG, levels=og.oder))


      p = ggtree::ggtree(sco.106.concat.tree) + ggtree::geom_tiplab(align=TRUE)

      #g = ggtree(tree) + geom_tiplab(align=TRUE)
      p <- ggtree::facet_plot(p,
                              panel=c.c,
                              mapping = aes(x = start, xend = end, yend=y,color=OG),
                              data=clusdetect.center %>%
                                mutate(id=genome) %>%
                                select(id, everything()),
                              geom=geom_segment,size=1.5)+
        theme(legend.position = "none")+
        viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
        ggnewscale::new_scale_color()

      p

      ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/proximity_clusters_unk_",c.c,".pdf"),
                      plot=p,
                      width = 12,
                      height = 8,
                      unit='in')

      all.out.compiled <- rbind(all.out.compiled, clusdetect.center %>%
                                  mutate(id=genome) %>%
                                  select(id, everything()) %>%
                                  mutate(c=c.c))
    }
  }
}


#--------------------------------------------------------------#









clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)
#
mapOGs <- detect_sco_in_region(clusdetect.sorted)
mapOG = mapOGs[1]
if(is.null(mapOG)){
message('consider alternative')
}
#
clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
 clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
focusOG=selectedOG
#
orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()
#
clusdetect.center <- clusdetect.center %>%
dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
dplyr::mutate(OG = factor(OG, levels=og.oder))


p = ggtree::ggtree(sco.106.concat.tree) + ggtree::geom_tiplab(align=TRUE)


p <-    ggtree::facet_plot(p,
                            panel=clustername,
                            mapping = aes(x = start, xend = end, yend=y,color=OG),
                            data=clusdetect.center %>%
                              mutate(id=genome) %>%
                              select(id, everything()),
                            geom=geom_segment,size=1.5)+
theme(legend.position = "none")+
viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
ggnewscale::new_scale_color()

p

#       ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.phylo_22_26_.",sel.col,'_',sel.grp,'_',gnm,'_',o,'_proximity_clusters',".pdf"),
#                       plot=p,
#                       width = 12,
#                       height = 8,
#                       unit='in')



#======

### Analyse summary out ####
prox.stat.genome <- prox.summary.out %>%
  ungroup() %>%
  group_by(genome,level,group) %>%
  summarise(n_enriched=sum(n_enriched),
            total_enriched=total_enriched,
            perc_enriched=n_enriched/total_enriched,
            not_assigned = total_enriched - n_enriched ,
            n_genes = sum(n_genes), #number of genes in clusters
            total_length=sum(cluster_length), #total lenght of gene clusters in bp
            n_transposase=sum(n_transposase), #total number of transposases in clusters
            known_with_pfam = sum(known_with_pfam),
            perc_known_with_pfam = known_with_pfam/n_genes,
            known_without_pfam = sum(known_with_pfam),
            perc_known_without_pfam = known_without_pfam/n_genes,
            unknown_with_pfam = sum(unknown_with_pfam),
            perc_unknown_with_pfam = unknown_with_pfam/n_genes,
            unknown_without_pfam = sum(unknown_with_pfam),
            perc_unknown_without_pfam = unknown_without_pfam/n_genes
            ) %>%
  unique() %>% ungroup()



write.table(prox.summary.out,
            file=paste0("~/DATA/MarbGenomics/Graphs/operons/prox.summary.out2.",'_',".txt"),
            quote=FALSE,
            sep='\t',
            row.names = FALSE)

write.table(prox.stat.genome,
            file=paste0("~/DATA/MarbGenomics/Graphs/operons/prox.stat.genome.2",'_',".txt"),
            quote=FALSE,
            sep='\t',
            row.names = FALSE)






prox.stat.genome %>%
  tidyr::pivot_longer(cols=c(n_enriched,not_assigned),names_to='enriched',values_to='n') %>%
  ggplot(aes(x=genome,y=n,fill=enriched)) +
  geom_bar(position='stack',stat='identity')

tmp.tree <- ape::drop.tip(sco.106.concat.tree,tip=setdiff(sco.106.concat.tree$tip.label, unique(prox.summary.out$genome)))



library(ggstance)


p <- ggtree::ggtree(tmp.tree)
p$data <- p$data %>% left_join(prox.stat.genome %>% select(genome,group), by=c('label'='genome'))
p <- p + geom_tippoint(aes(color=group))+ scale_color_manual(values=col.list[[2]]) +
  ggtree::geom_tiplab(aes(color=group),align=TRUE)+ theme_tree2()


p <- ggtree::facet_plot(p,
                   panel='sig loci',
                   mapping = aes(x = n, fill=enriched),
                   data= prox.stat.genome %>% mutate(id=genome) %>% select(id, everything()) %>%
                     tidyr::pivot_longer(cols=c(n_enriched,not_assigned),names_to='enriched',values_to='n') %>%
                     mutate(enriched=ifelse(enriched=='n_enriched','xn_assigned','not_assigned')),
                   geom=geom_barh,stat='identity')+
  theme(legend.position = "none")+
  ggnewscale::new_scale_color() + scale_fill_manual(values=c('grey','black'))

p <- ggtree::facet_plot(p,
                        panel='total_length(Mb)',
                        mapping = aes(x = total_length/1000),
                        data= prox.stat.genome %>% mutate(id=genome) %>% select(id, everything()),
                        geom=geom_barh,stat='identity')+
  theme(legend.position = "none")+
  ggnewscale::new_scale_color() + scale_fill_manual(values=c('grey','black'))

p <- ggtree::facet_plot(p,
                        panel='perc_hypothetical',
                        mapping = aes(x = perc_hypothetical),
                        data= prox.stat.genome %>% mutate(id=genome) %>% select(id, everything()),
                        geom=geom_barh,stat='identity')+
  theme(legend.position = "none")+
  viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
  ggnewscale::new_scale_color() + scale_fill_manual(values=c('grey','black'))


p <- ggtree::facet_plot(p,
                        panel='n_transposase',
                        mapping = aes(x = n_transposase),
                        data= prox.stat.genome %>% mutate(id=genome) %>% select(id, everything()),
                        geom=geom_barh,stat='identity')+
  theme(legend.position = "none")+
  ggnewscale::new_scale_color()+ scale_fill_manual(values=c('grey','black')) + ggnewscale::new_scale_color()



ggtree::facet_plot(p + ggnewscale::new_scale_fill(),
                   panel='perc_hypothetical',
                   mapping = aes(x = value,fill=name),
                   data= prox.stat.genome %>% mutate(id=genome) %>% select(id, everything()) %>%
                     mutate(annotated = n_genes- n_hypothetical,hyp_with_pfam=n_hypothetical-n_poor, hyp_without_pfam=n_poor) %>%
                     tidyr::pivot_longer(cols=c(annotated,hyp_with_pfam,hyp_without_pfam)),
                   geom=geom_barh,stat='identity')+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete=TRUE,option=v.options[2],na.value='grey') +
  ggnewscale::new_scale_fill()





###Figure caption ####
# characteristics of co-localised gene clusters
# number of lineage specific loci that group into gene clusters (black) (unassigned loci, grey) (overall average is 47.8%
# total lenght of gene clusters in bp
# percentage of hypothetical proteins in geneclusters
# number of transposases in gene clusters

ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.phylo_summarised2.",sel.col,".pdf"),
                plot=p,
                width = 6,
                height = 5,
                unit='in')






p.hist <- prox.summary.out %>% arrange(desc(n_genes)) %>% ggplot(aes(x=n_genes)) + geom_histogram(binwidth = 1) + scale_y_continuous(expand=c(0,0)) + fdb_style(aspect.ratio=0.5)


ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.phylo_nGenes_histogram.",sel.col,".pdf"),
                plot=p.hist,
                width = 4,
                height = 4,
                unit='in')




prox.stat.genome %>% ggplot(aes(n_enriched, total_length)) + geom_point(aes(color=group)) + fdb_style()


prox.stat.genome %>% ggplot(aes(group, total_length)) + geom_jitter(aes(color=group)) + fdb_style()

prox.stat.genome %>% ggplot(aes(group, perc_enriched)) + geom_jitter(aes(color=group)) + fdb_style()

prox.stat.genome %>% ungroup() %>% group_by(level) %>% summarise(mean=mean(perc_enriched))


p.out <- prox.stat.genome %>%
  select(group, n_enriched, perc_enriched, total_length, n_transposase) %>%
  tidyr::pivot_longer(cols=-group) %>%
  mutate(name=factor(name,levels=c('n_enriched','perc_enriched', 'total_length', 'n_transposase'))) %>%
  mutate(group=factor(group,levels=rev(c('cl1','cl10','cl12','cl3','cl4','cl5','cl8'))))%>%
  ggplot(aes(group, value)) +
  geom_bar(data=.%>% group_by(group,name) %>% summarise(mean=mean(value)),
           aes(x=group, y=mean,fill=group),stat='identity',alpha=0.6) +
  geom_errorbar(data=.%>% group_by(group,name) %>% summarise(mean=mean(value),sd=sd(value)),
                aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd,color=group),
                width=0.2)+
  geom_jitter(aes(color=group)) +
  fdb_style()+coord_flip()+
  scale_color_manual(values=grid.colors) +
  scale_fill_manual(values=grid.colors) + facet_wrap(~name,nrow=1,scales = 'free_x')


ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.stat.genome_grouped.",sel.col,".pdf"),
                plot=p.out,
                width = 8,
                height = 5,
                unit='in')



prox.stat.genome %>%
  select(group, n_enriched, perc_enriched, total_length, n_transposase) %>%
  tidyr::pivot_longer(cols=-group) %>%
  mutate(name=factor(name,levels=c('n_enriched','perc_enriched', 'total_length', 'n_transposase'))) %>%
  mutate(group=factor(group,levels=rev(c('cl1','cl10','cl12','cl3','cl4','cl5','cl8'))))%>%
  ggplot(aes(group, value)) +
  geom_bar(data=.%>% group_by(group,name) %>% summarise(mean=mean(value)),
           aes(x=group, y=mean),color='black',stat='identity',alpha=0.6) +
  geom_errorbar(data=.%>% group_by(group,name) %>% summarise(mean=mean(value),sd=sd(value)),
                aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd),
                width=0.2)+
  geom_jitter(color='grey20') +
  #fdb_style()+
  coord_flip()+
  facet_wrap(~name,nrow=1,scales = 'free_x')



prox.stat.genome %>%
  select(group, n_enriched, perc_enriched, total_length, n_transposase) %>%
  tidyr::pivot_longer(cols=-group) %>%
  mutate(name=factor(name,levels=c('n_enriched','perc_enriched', 'total_length', 'n_transposase'))) %>%
  mutate(group=factor(group,levels=rev(c('cl1','cl10','cl12','cl3','cl4','cl5','cl8'))))%>%
  ggplot(aes(group, value)) +
  geom_bar(data=.%>% group_by(group,name) %>% summarise(mean=mean(value)),
           aes(x=group, y=mean),color='grey30',fill='grey',stat='identity',width=0.7) +
  geom_errorbar(data=.%>% group_by(group,name) %>% summarise(mean=mean(value),sd=sd(value)),
                aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd),
                width=0.2)+
  #geom_jitter(color='grey20',shape=21) +
  scale_y_continuous(expand=c(0,0))+
  coord_flip()+
  facet_wrap(~name,nrow=1,scales = 'free_x') +
  fdb_style()
  #theme_bw()+
  #theme(
  #  strip.background=element_blank(),
  #  strip.text.x = element_blank(),
  #  panel.grid.major = element_blank(),
  #  panel.grid.minor = element_blank(),
  #  axis.title = ggplot2::element_text(size = 11, colour = '#000000'),
  #  axis.text = ggplot2::element_text(size = 10, colour = '#000000'),
  #
  #  legend.justification = c(1, 1),
  #  legend.key.width = unit(0.25, 'cm'),
  #  legend.key.height = unit(0.55, 'cm'),
  #  legend.text = ggplot2::element_text(size = 10),
  #  legend.title = ggplot2::element_text(size = 11),
  #  panel.background = ggplot2::element_blank(),
  # )
  #)










#annotation quality derived from OG
prox.annotated.out %>% left_join(annotation.q.per.OG, by='OG') %>%
  filter(!is.na(annotation_quality)) %>%
  group_by(genome,annotation_quality) %>% tally() %>%
  ggplot2::ggplot(aes(x=genome,y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="stack", stat="identity") +
  #geom_point(shape=21) +
  coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  geom_hline(yintercept = 0)+
  fdb_style(aspect.ratio=0.5)+
  ylab('No of significant gene families') + xlab(' ')


#annotation quality derived from OG
prox.annotated.out %>% left_join(annotation.q.per.OG, by='OG') %>%
  filter(!is.na(annotation_quality)) %>%
  group_by(group.x,annotation_quality) %>% tally() %>%
  ggplot2::ggplot(aes(x=group.x,y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="stack", stat="identity") +
  #geom_point(shape=21) +
  coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  geom_hline(yintercept = 0)+
  fdb_style(aspect.ratio=0.5)+
  ylab('No of significant gene families') + xlab(' ')








sel.loci <- paste0('HJKODCBM_033',15:45)
sel.loci <- paste0('HJKODCBM_02',226:303)
sel.loci <- paste0('HJKODCBM_02',395:415)

loci.list <- list(paste0('HJKODCBM_033',15:45),
                  paste0('HJKODCBM_02',226:303),
                  paste0('HJKODCBM_02',395:415))
names(loci.list) <- c('T6SS','flagellum','pilus')



loci.list <- list(c('MLBPJMHG_04007','MLBPJMHG_04008','MLBPJMHG_04009','MLBPJMHG_04010'),
                  c('MLBPJMHG_01527','MLBPJMHG_01528','MLBPJMHG_01529','MLBPJMHG_01530'),
                  c('MLBPJMHG_03085','MLBPJMHG_03086','MLBPJMHG_03087','MLBPJMHG_03088'),
                  paste0('IFOEGCPI_036',84:90),
                  paste0('IFOEGCPI_037',10:20),
                  paste0('NHGOFBLF_003',35:46))

names(loci.list) <- c('ADI','pta_acka','pyr_dehyr','nar','nor','nuo')

annotation.tbl %>% filter(locus_tag %in% sel.loci)




cluster.out.list <-list()
  for(sn in names(loci.list)){
    sel.loci = loci.list[[sn]]
    prox.clust <- proximity_clusters(annotation.tbl, sel.loci,distance=3000 )
    sel.clust <- 'cluster1'
    selectedOG <- prox.clust %>% filter(prox.cluster == sel.clust) %>% select(OG) %>% pull() %>% as.character()
    clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
    clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
    clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
    clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)

    mapOGs <- detect_sco_in_region(clusdetect.sorted)
    mapOG = mapOGs[1]
    if(is.null(mapOG)){
      message('consider alternative')
    }

    clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
    clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
    clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
    focusOG=selectedOG

    orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()

    clusdetect.center <- clusdetect.center %>%
      dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
      dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
      dplyr::mutate(OG = factor(OG, levels=og.oder))

    cluster.out.list[[sn]] <- clusdetect.center
  }





#=====#
#alternative using proxclusters
prox.clust=prox.clust2

largest.clusts_all <- prox.clust2 %>% select(prox.cluster) %>% unique() %>% pull
cluster.out.list <-list()
for(sel.clust in largest.clusts_all){
  selectedOG <- prox.clust %>% dplyr::filter(prox.cluster == sel.clust) %>% dplyr::filter(!is.na(OG)) %>% dplyr::select(OG) %>% pull() %>% as.character()
  if(length(selectedOG)>0){
    clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
    clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
    clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
    clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)

    mapOGs <- detect_sco_in_region(clusdetect.sorted)
    mapOG = mapOGs[1]
    if(is.null(mapOG)){
      message('consider alternative')
    }

    clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
    clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
    clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
    focusOG=selectedOG

    orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()

    clusdetect.center <- clusdetect.center %>%
      dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
      dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
      dplyr::mutate(OG = factor(OG, levels=og.oder))

    cluster.out.list[[sel.clust]] <- clusdetect.center
  }
}




v.options <- rep(c('viridis','plasma','magma','cividis','inferno','mako'),10)
p = ggtree(sco.106.concat.tree) + geom_tiplab(align=TRUE)

for(i in 1:length(names(cluster.out.list))){
  message(v.options[i])
  p <-    facet_plot(p,
                     panel=names(cluster.out.list)[i],
                     mapping = aes(x = start, xend = end, yend=y,color=OG),
                     data=cluster.out.list[[i]] %>%
                       mutate(id=genome) %>%
                       select(id, everything()),
                     geom=geom_segment,size=1.5)+
    theme(legend.position = "none")+
    viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
    ggnewscale::new_scale_color()
  #for some reason need to print in order for colors to work
  #print(p)
}


ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.phylo_22_22226_.",sel.col,'_',sel.grp,'_',gnm,'_proximity_clusters',".pdf"),
                plot=p,
                width = 12,
                height = 8,
                unit='in')











  out.p <- clusdetect.center %>%
    dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
    dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
    dplyr::mutate(OG = factor(OG, levels=og.oder)) %>%
    ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=Name,fill=OG)) + #color=SELCOG
    gggenes::geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) +
    gggenes::geom_gene_label(align = "left")+
    #viridis::scale_fill_viridis(discrete = TRUE)+
    #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
    theme(legend.position = "none")

  out.p
#----------------------#

  clusdetect.center %>% #filter(genome=='Marinobacter_sp_FDB33') %>%
    dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
    dplyr::mutate(OG = factor(OG, levels=og.oder)) %>%
    ggplot(aes(x = start, xend = end, y = genome, yend=genome, forward=direction,color=OG)) + #color=SELCOG
    geom_segment(size=2) +
    # gggenes::geom_gene_label(align = "left")+
    #ggplot2::facet_wrap(~genome,ncol = 2) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )




plot.data <- clusdetect.center %>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(id=genome) %>%
  select(id, everything()) %>% mutate(gene=OG)%>% data.frame()

g = ggtree(sco.106.concat.tree) + geom_tiplab(align=TRUE)
p <- facet_plot(g,
                panel='alignment',
                mapping = aes(x = start, xend = end, yend=y,color=OG),
                data=plot.data,
                geom=geom_segment,size=1.5)+
  theme(legend.position = "none")+
  #scale_fill_manual(values=c('blue','red','black','white'))+
  theme(strip.text=element_blank(),panel.spacing=unit(0, 'cm'))
p






#=================================================================================#
##### THe reverse #####
#=================================================================================#
# aim is to detect absent operons, operons detected on depleted genes, but in other genomes...
# sel.grp='cl5'
# gnm = 'Marinobacter_vinifirmus_FB1'

for(i in 1:nrow(listedRepresentatives)){
  sel.grp <- listedRepresentatives[i,1]
  gnm <- listedRepresentatives[i,2]

  top100RF <- HGM_RF_COMB %>%
    filter(group.x==sel.grp) %>%
    arrange(desc(MeanDecreaseAccuracy)) %>%
    head(100) %>%
    select(feature.id.x) %>%
    pull() %>%
    as.character()

  selected.features <- HGM_RF_COMB %>%
    filter(group.x==sel.grp) %>%
    mutate(top100RF = ifelse(feature.id.x %in% top100RF,TRUE,FALSE)) %>%
    mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    mutate(Deple_bin = ifelse(fdr.p.value_depletion_binary<0.1,TRUE,FALSE)) %>%
    mutate(Deple_raw = ifelse(fdr.value_depletion_rawcounts<0.1,TRUE,FALSE)) %>%
    filter(Deple_bin ==TRUE | Deple_raw ==TRUE | top100RF == TRUE) %>%
    select(feature.id.x) %>%
    pull() %>%
    as.character()

  RF.250 <- HGM_RF_COMB %>%
    filter(group.x==sel.grp) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
    slice(1:250) %>%
    select(feature.id.x) %>%
    pull() %>%
    as.character()

  #combine HGM and RF data
  out.mat <- HGM_RF_COMB %>%
    filter(group.x==sel.grp) %>%
    select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
    mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    mutate(RFtop250 = ifelse(feature.id.x %in% RF.250,TRUE,FALSE))

  message(length(selected.features))

  if(length(selected.features)!=0){
    gnm2= gnm


    for(gnm in c('Marinobacter_hydrocarbonoclasticus_ATCC_49840','Marinobacter_adhaerens_HP15')){
      #look for loci in other genome
      sel.loci <- annotation.tbl %>%
        filter(genome %in% gnm) %>%
        filter(OG %in% selected.features) %>%
        filter(type=='CDS') %>% dplyr::filter(!is.na(OG))%>%
        mutate(OG = coalesce(OG,locus_tag)) %>%
        dplyr::select(locus_tag) %>% pull()

      #analyse gene proximity
      prox.clust <- proximity_clusters(annotation.tbl, sel.loci,distance=3000)
      filled.clust <- fillgaps2(in.cluster=prox.clust)
      prox.clust <- proximity_clusters(annotation.tbl, filled.clust$locus_tag,distance=3000)

      prox.annotated <- prox.clust %>%
        dplyr::left_join(ko2ann, by=c('kfm_domain'='ko')) %>%
        dplyr::left_join(cog_func_single, by='COG') %>%
       dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>%
        dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
        dplyr::select(-score,-note,-rpt_family,-rpt_type,-rpt_unit_seq)

      if(!is.null(prokka2ncbi.annotated.list[[gnm]])){
        message('adding NCBI accessions and annotations')
        prox.annotated <- prox.annotated %>% left_join(prokka2ncbi.annotated.list[[gnm]], by=c('locus_tag'='prokka_locus_tag'))
      }

      #save summary
      write.table(prox.annotated,
                  file=paste0("~/DATA/MarbGenomics/Graphs/operons/DEPLETED.prox.annotated.",gnm2,'_',gnm,".txt"),
                  quote=FALSE,
                  sep='\t',
                  row.names = FALSE)
    }
    }

}


if(!is.null(prox.clust)){
  # extract the largest clusters formed by >3 clusters
  largest.clusts_all <- prox.clust %>%
    group_by() %>%
    dplyr::count(prox.cluster) %>%
    filter(n>3) %>%
    arrange(desc(n)) %>%
    select(prox.cluster) %>% pull()

  if(length(largest.clusts_all)>0){
    filtered <- prox.clust %>% filter(prox.cluster %in% largest.clusts_all)
    filled <- fillgaps2(in.cluster=filtered)
    prox.clust2 <- proximity_clusters(annotation.tbl, filled$locus_tag,distance=3000 )

    #the alternative is to plot is same scale, start all the clusters at 0 position
    prox.clust.origin <- c()
    for(prx.c in unique(prox.clust2$prox.cluster)){
      t.start <- prox.clust2 %>% filter(prox.cluster == prx.c) %>% slice(1) %>% select(start) %>% pull()
      tmp <- prox.clust2 %>% filter(prox.cluster == prx.c) %>% mutate(start=start - t.start, end=end-t.start)
      prox.clust.origin = rbind(prox.clust.origin, tmp)
    }

# plot them all?
for(sel.clust in largest.clusts_all){
  selectedOG <- prox.clust %>% filter(prox.cluster == sel.clust) %>% select(OG) %>% pull() %>% as.character()

  clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
  clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
  clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
  clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)

  #clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=FALSE, OGs=mapOG)
  mapOGs <- detect_sco_in_region(clusdetect.sorted)
  mapOG = mapOGs[1]
  if(is.null(mapOG)){
    message('consider alternative')
  }

  # -----
  #Supply more than one OG to OGs to organise them around center (if not all have an SCO in this region)
  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)

  #use function to add dummies
  clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)

  clusdetect.center$genome <- factor(clusdetect.center$genome,levels = tip.order)

  focusOG=selectedOG

  orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
  og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()

  out.p <- clusdetect.center %>%
    mutate(direction=ifelse(strand=="+",1,-1)) %>%
    mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
    mutate(OG = factor(OG, levels=og.oder)) %>%
    ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=Name,fill=OG)) + #color=SELCOG
    gggenes::geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) +
    gggenes::geom_gene_label(align = "left")+
    viridis::scale_fill_viridis(discrete = TRUE)+
    #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
    theme(legend.position = "none")

  out.p

  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/operons/",'depleted',gnm2,'_',sel.clust,"_testing4",'clusterdistribution',".pdf"),
                  plot=out.p,
                  width = 7,
                  height = 14,
                  unit='in')
  }
}
}


prox.clust %>% filter(prox.cluster %in% largest.clusts_all) %>%
  left_join(KO_annotated %>% select(-genome) %>% unique(), by='locus_tag') %>%
  select(genome, locus_tag, product, Name, eC_number, prox.cluster, KO, module,ko_annotation,MOD_F, TYPE, DESCRIPTION) %>% data.frame





#===================================================≠≠#
###?????####
#pfam per locus

#===========

#
all250topRF <- rf.all %>% group_by(group) %>% arrange(desc(MeanDecreaseAccuracy)) %>% slice(1:250) %>% select(feature.id) %>% pull()


allSig_hgm <- HGM_RF_COMB %>%
  dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts,fdr.p.value_depletion_binary,fdr.value_depletion_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
  dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::select(feature.id.x, group.x, Enrich_bin,Depl_bin,Depl_raw, Enrich_raw) %>%
  dplyr::mutate(signEnriched = ifelse(Enrich_bin == TRUE | Enrich_raw==TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(signDepleted = ifelse(Depl_bin == TRUE | Depl_raw==TRUE, TRUE, FALSE)) %>%
  dplyr::filter(signEnriched==TRUE | signDepleted==TRUE) %>% select(feature.id.x) %>% pull() %>% unique() %>% as.character()


allSig_combined <- unique(c(all250topRF,allSig_hgm))

allSig_hgm %>%length()
all250topRF %>%length()
allSig_combined %>%length()




col_fun = circlize::colorRamp2(c(0.6,0.8, 1), c("deepskyblue2", "black", "yellow"))

metadata <- genome.tbl %>% filter(genome %in% sco.94.concat.tree.rooted$tip.label)# %>% data.frame()


tip.order <- rev(tip_order(sco.94.concat.tree.rooted))

plotdat <-genome.tbl %>%
  left_join(sorted.cliques,by=c('genome'='genome')) %>%
  left_join(out.colors, by=c('genome'='genome')) %>%
  filter(genome %in% sco.94.concat.tree.rooted$tip.label) %>% data.frame()

rownames(plotdat) = plotdat$genome

#metadata <- metadata %>% mutate(phylogroup=phylogroup.x)
ph2gen <- metadata$phylogroup
names(ph2gen) <- metadata$genome

type2gen <- metadata$TypeStrain
names(type2gen) <- metadata$genome



th2 <- ComplexHeatmap::HeatmapAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                         phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                         ANI0.8 = paste0('cl',plotdat[rev(tip.order),'sc0.8']),
                                         ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),
                                         ANI0.9 = paste0('cl',plotdat[rev(tip.order),'sc0.9']),
                                         ANI0.95 = paste0('cl',plotdat[rev(tip.order),'sc0.95']),
                                         ANI0.98 = paste0('cl',plotdat[rev(tip.order),'sc0.98']),

                                         col=list(
                                           phylogroup = phylogroup.colors,
                                           ANI0.8 = col.list[[1]],
                                           ANI0.85 = col.list[[2]],
                                           ANI0.9 = col.list[[3]],
                                           ANI0.95 = col.list[[4]],
                                           ANI0.98 = col.list[[5]],
                                           type=c('TRUE'='black','FALSE'='white')#,
                                         )
)

th3 <- ComplexHeatmap::rowAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                     phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                     ANI0.8 = paste0('cl',plotdat[rev(tip.order),'sc0.8']),
                                     ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),
                                     ANI0.9 = paste0('cl',plotdat[rev(tip.order),'sc0.9']),
                                     ANI0.95 = paste0('cl',plotdat[rev(tip.order),'sc0.95']),
                                     ANI0.98 = paste0('cl',plotdat[rev(tip.order),'sc0.98']),

                                     col=list(
                                       phylogroup = phylogroup.colors,
                                       ANI0.8 = col.list[[1]],
                                       ANI0.85 = col.list[[2]],
                                       ANI0.9 = col.list[[3]],
                                       ANI0.95 = col.list[[4]],
                                       ANI0.98 = col.list[[5]],
                                       type=c('TRUE'='black','FALSE'='white')#,
                                     )
)

th4 <- ComplexHeatmap::rowAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                     phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                     ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),

                                     col=list(
                                       phylogroup = phylogroup.colors,
                                       ANI0.85 = col.list[[2]],
                                       type=c('TRUE'='black','FALSE'='white')#,
                                     )
)


th <- ComplexHeatmap::HeatmapAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                        ANI0.98 = paste0('cl',plotdat[rev(tip.order),'sc0.98']),
                                        ANI0.95 = paste0('cl',plotdat[rev(tip.order),'sc0.95']),
                                        ANI0.9 = paste0('cl',plotdat[rev(tip.order),'sc0.9']),
                                        ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),
                                        ANI0.8 = paste0('cl',plotdat[rev(tip.order),'sc0.8']),
                                        phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                        col=list(
                                          phylogroup = phylogroup.colors,
                                          ANI0.8 = col.list[[1]],
                                          ANI0.85 = col.list[[2]],
                                          ANI0.9 = col.list[[3]],
                                          ANI0.95 = col.list[[4]],
                                          ANI0.98 = col.list[[5]],
                                          type=c('TRUE'='black','FALSE'='white')#,
                                        )
)

rh <- ComplexHeatmap::rowAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                    phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                    col=list(phylogroup = phylogroup.colors,
                                             type=c('TRUE'='black','FALSE'='white')#,
                                    )
)

# For AAI
col_fun2 = circlize::colorRamp2(c(0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1), viridis::viridis(8))
col_fun = circlize::colorRamp2(c(0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1), colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(8))



allSig_combined

pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_allsig_hgm_3.pdf",sep=''), width=8, height=4)
ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_hgm]),
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()



pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_allsig_hgm_4.pdf",sep=''), width=8, height=4)
ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_combined]),
                        left_annotation= th4,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()








pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_AAI_upsidedown.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(AAI.sel),
                        top_annotation = th2,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE
)
dev.off()

pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_AAI.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(AAI.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE
)
dev.off()

pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_AAI_viridis.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(AAI.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()


# For ANIb
col_fun2 = circlize::colorRamp2(c(0.7,0.75,0.8,0.85,0.9,0.95, 1), viridis::viridis(7))
col_fun = circlize::colorRamp2(c(0.7,0.75,0.8,0.85,0.9,0.95, 1), colorRampPalette(brewer.pal(9, "Blues"))(7))


pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_ANIb.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(ANIb.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()

pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_ANIb_viridis.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(ANIb.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()




pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_GRR_viridis.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(GRR.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()


#================================
pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_full.pdf",sep=''), width=14, height=4)

ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), ]),
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_allsig_hgm.pdf",sep=''), width=14, height=4)

ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_hgm]),
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()


allSig_combined

pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_allsig_hgm_3.pdf",sep=''), width=8, height=4)
ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_hgm]),
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()




pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam__heatmap_v.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pan_clean(pfam.pan[rev(tip.order), sigpfam ])),
                        left_annotation= rh,
                        #col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()

ComplexHeatmap::Heatmap(as.matrix(pan_clean(pfam.pan[rev(tip.order), sig.genes ])),
                        left_annotation= rh,
                        #col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)



pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('Cob|Dct|amylase|TPR|Ank|Transposase',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)

sig.genes[grepl('Transposase',sig.genes)]


pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_HTH_hmg_heatmap.pdf",sep=''), width=7, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('HTH',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()



pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_allegaartje_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('Nap|rhodopsin|lactamas|Cyt|toxin|MCP|resist|Phage|transposase|Tnp',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_DUF_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('DUF',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10), c(colorRampPalette(brewer.pal(9, "Blues"))(10))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_NAD_ATP_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('ATP|NAD',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('Gly|aldo|argin|alg|aldedh|pqq|malt|Mtase',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(seq(0,20,by=1), c(colorRampPalette(brewer.pal(9, "Blues"))(21))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)



groups <- ANI.hmg.grp.tbl %>%
  filter(level=='sc0.85') %>% select(group) %>% unique() %>% pull %>% as.character() %>% sort

for(sel.grp in groups){
  message(sel.grp)
  ssp <- ANI.hmg.grp.tbl %>%
    filter(level=='sc0.85') %>%
    filter(fdr.p.value_depletion_binary<0.05) %>%
    filter(group==sel.grp) %>% select(feature.id) %>%
    pull %>% as.character()


  pdf( paste0("~/DATA/MarbGenomics/Graphs/","pfam_hmg_heatmap_de",sel.grp,".pdf"), width=0.2*length(ssp), height=4)
  print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(ssp)]),
                                left_annotation= th3,
                                col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                                show_row_names = FALSE,
                                #show_column_names = FALSE,
                                cluster_rows = FALSE,
                                cluster_columns = TRUE))
  dev.off()
}



sel.grp='cl1'

mda <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(MeanDecreaseAccuracy) %>% pull()
names(mda) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

penrbin <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_enriched_binary) %>% pull() %>%
  names(penrbin) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

penrbin <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_enriched_binary) %>% pull()
names(penrbin) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

penrraw <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_enriched_rawcounts) %>% pull()
names(penrraw) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

pdepbin <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_depletion_binary) %>% pull()
names(pdepbin) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

pdepraw <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.value_depletion_rawcounts) %>% pull()
names(pdepraw) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()


sel.features <- sig.genes

ch <- ComplexHeatmap::HeatmapAnnotation(MeanDecrAcc = mda[sel.features],
                                        penrbin=-log10(penrbin[sel.features]),
                                        penrraw=-log10(penrraw[sel.features]),
                                        pdepbin=-log10(pdepbin[sel.features]),
                                        pdepraw=-log10(pdepraw[sel.features]),
                                        col=list(
                                          mda = circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          penrraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          pdepbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          pdepraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")))
)


print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order),sel.features ]),
                              left_annotation= th3,
                              top_annotation= ch ,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



mda1 <- HGM_RF_COMB %>% filter(group.x=='cl1') %>% select(MeanDecreaseAccuracy) %>% pull()
mda2 <- HGM_RF_COMB %>% filter(group.x=='cl2') %>% select(MeanDecreaseAccuracy) %>% pull()
mda3 <- HGM_RF_COMB %>% filter(group.x=='cl3') %>% select(MeanDecreaseAccuracy) %>% pull()
mda4 <- HGM_RF_COMB %>% filter(group.x=='cl4') %>% select(MeanDecreaseAccuracy) %>% pull()
mda5 <- HGM_RF_COMB %>% filter(group.x=='cl5') %>% select(MeanDecreaseAccuracy) %>% pull()
mda6 <- HGM_RF_COMB %>% filter(group.x=='cl6') %>% select(MeanDecreaseAccuracy) %>% pull()
mda7 <- HGM_RF_COMB %>% filter(group.x=='cl7') %>% select(MeanDecreaseAccuracy) %>% pull()
mda8 <- HGM_RF_COMB %>% filter(group.x=='cl8') %>% select(MeanDecreaseAccuracy) %>% pull()
mda9 <- HGM_RF_COMB %>% filter(group.x=='cl9') %>% select(MeanDecreaseAccuracy) %>% pull()
mda10 <- HGM_RF_COMB %>% filter(group.x=='cl10') %>% select(MeanDecreaseAccuracy) %>% pull()
mda11 <- HGM_RF_COMB %>% filter(group.x=='cl11') %>% select(MeanDecreaseAccuracy) %>% pull()

names(mda1) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda2) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda3) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda4) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda5) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda6) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda7) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda8) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda9) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda10) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda11) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()



ch <- ComplexHeatmap::HeatmapAnnotation(cl1=mda1[sel.features],
                                        cl2=mda2[sel.features],
                                        cl3=mda3[sel.features],
                                        cl4=mda4[sel.features],
                                        cl5=mda5[sel.features],
                                        cl6=mda6[sel.features],
                                        cl7=mda7[sel.features],
                                        cl8=mda8[sel.features],
                                        cl9=mda9[sel.features],
                                        cl10=mda10[sel.features],
                                        cl11=mda11[sel.features]
                                        #col=list(
                                        #  mda = circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  penrraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  pdepbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  pdepraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")))
)


print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order),sel.features ]),
                              left_annotation= th3,
                              top_annotation= ch ,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), grepl('RCC|NLR|LRR',colnames(pfam.pan))]),
                              left_annotation= rh,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))


print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), grepl('Arm-DNA-bind_1|recombinase|resolvase',colnames(pfam.pan),ignore.case=T)]),
                              left_annotation= rh,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), grepl('Arm-DNA-bind_1|recombinase|resolvase',colnames(pfam.pan),ignore.case=T)]),
                              left_annotation= rh,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))




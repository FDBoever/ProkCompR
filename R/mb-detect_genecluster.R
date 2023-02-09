#   DETECT GENE CLUSTERS IN GIVEN SET OF LOCUS TAGS?
#-----------------------------------------------------------------------------#
#   UTILITY FUNCTIONS   ####
#-----------------------------------------------------------------------------#

#' proximity_clusters is a function to group locus_tags into clusters based on gene proximity
#'
#' @param an.tbl annotation table object
#' @param ltags vector containing locus tags
#' @param distance intergenic distance bewteen two locus_tags
#'
#' @return
#' @export
#'
#' @examples
#' proximity_clusters(an.tbl, ltags , distance=3000)

# proximity_clusters
#=================================================================================#
proximity_clusters <- function(an.tbl, ltags, distance = 3000){
  f.tbl <- an.tbl %>% filter(locus_tag %in% ltags)

  contigs <- f.tbl %>%
    dplyr::select(seqnames) %>%
    dplyr::pull() %>%
    as.character() %>% unique()

  tracker <- 0 #stores the cluster number
  output <- NULL
  for(contig in contigs){
    #subset table
    sub.tbl <-  f.tbl %>% dplyr::filter(seqnames == contig )

    glued.tbl <- sub.tbl %>%
      dplyr::mutate(diff = start - lag(end,default=first(end))) %>%
      dplyr::mutate(glued = ifelse(diff < distance, TRUE, FALSE))# %>%
      #dplyr::select(seqnames, locus_tag,start, end , diff, glued)

    glued.tbl[1,"glued"] <- FALSE

    #logic, for each false, make new group
    clusters <- c()
    for(o in glued.tbl$glued){
      if(!o){
        tracker <- tracker+1
        clusters <- c(clusters,paste0('cluster',tracker))
      }else{
        clusters <- c(clusters,paste0('cluster',tracker))
      }
    }

    glued.tbl <- glued.tbl %>%
      dplyr::mutate(prox.cluster = clusters)

    output <- rbind(output, glued.tbl)

  }
  return(output)
}



proximity_clusters(an.tbl, ltags,distance=3000 )





# ANALYSIS
#=================================================================================#

#selected.og Can be calculated for example in Gene enrichment script
selected.og

selected.og <- hmg.grp.tbl %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::filter(group == 'cl14') %>%
  dplyr::filter(fdr.p.value_enriched_binary < 0.05) %>% dplyr::select(feature.id) %>% dplyr::pull()


#=================================================================================#
#=================================================================================#
#=================================================================================#

outdir = '~/DATA/MarbGenomics/Graphs/'

for(sel.col in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  #Select relevant metadata
  #sel.col <- 'sc0.85'
  mtd <- genome.tbl %>% left_join(sorted.ANI.cliques, by='genome') %>%
    left_join(out.colors, by='genome') %>%
    data.frame()
  mtd$grp <- mtd[,sel.col]
  mtd <- mtd %>% mutate(group=paste0('cl',grp))

  #set color scheme
  colors <- mtd %>% select(paste0('col.',sel.col), group) %>% unique()
  colors <-colors[order(colors$group),]
  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,2])
  clique.colors <- grid.colors

  # select representative genomes
  # GENERATE A LITTLE TABLE WITH CLIQUE REPRESENTATIVES
  representatives <- mtd %>% select(genome,group,TypeStrain,Completeness) %>%
    dplyr::arrange(genome,desc(Completeness),desc(TypeStrain),group_by=group) %>%
    group_by(group) %>% slice(1) %>% select(genome, group)

  clique2representative <- representatives$genome
  names(clique2representative)  <- representatives$group

  for(sel.group in names(clique2representative)[names(clique2representative)!='clNA']){
    sel.genome <- unname(clique2representative[sel.group])

    sel.features <- ANI.hmg.grp.tbl %>%
      filter(level==sel.col) %>%
      filter(group==sel.group) %>%
      filter(fdr.p.value_enriched_binary<0.05) %>%
      select(feature.id) %>%
      pull() %>%
      as.character()

    if(length(sel.features)!=0){
      hgm.bin.enriched <- annotation.tbl %>%
        dplyr::filter(genome==sel.genome) %>%
        dplyr::filter(OG %in% selected.og) %>%
        dplyr::select(locus_tag) %>%
        dplyr::pull()

      out.prox <- proximity_clusters(annotation.tbl, hgm.bin.enriched,distance=3000 )

      #Select all >3 large proximity clusters
      prox.large <- out.prox %>% dplyr::group_by(prox.cluster) %>%
        select(prox.cluster) %>%
        dplyr::count() %>%
        dplyr::filter(n>3) %>%
        dplyr::select(prox.cluster) %>%
        dplyr::pull() %>% as.character()

      out.tbl <- NULL
      for(prox.clust in prox.large){
        message(prox.clust)
        cluster.loci <- out.prox %>%
          filter(prox.cluster == prox.clust) %>%
          select(locus_tag) %>%
          pull()

        sel.annotations <- annotation.tbl %>%
          filter(locus_tag %in% cluster.loci) %>%
          mutate(prox.cluster = prox.clust)

        selected.summary <- sel.annotations %>%
          select(locus_tag, product, OG, COG, gene, kfm_domain, cazy_domain,tcdb_domain)

        print(selected.summary)

        gaps.filled <- fillgaps(in.cluster=sel.annotations) %>%
          mutate(hit = ifelse(locus_tag %in% sel.annotations$locus_tag, TRUE, FALSE)) %>%
          mutate(prox.cluster = prox.clust) %>%
          dplyr::mutate(direction=ifelse(strand=='+',1,-1))

        out.tbl <- rbind(out.tbl, gaps.filled)
      }

      p.plt <- ggplot2::ggplot(out.tbl, aes(xmin = start, xmax = end, y = seqnames,forward=direction,label=Name,fill=prox.cluster,alpha=as.numeric(hit))) +
        gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
        gggenes::geom_gene_label(align = "left")+
        facet_wrap(~ prox.cluster, scales = "free", ncol = 1) +
        theme_genes()+
        theme(legend.position = "none")+scale_alpha(range=c(0.5,1))+ggtitle(paste0(sel.genome,' - ',sel.col, ' - ', sel.group ))

      ggplot2::ggsave(filename=paste0(outdir,'proximity_cluster_',sel.genome,sel.col,'.pdf'),plot=p.plt, width = 10, height = 25,unit='in')
    }
  }
}









#how many do form clusters?
out.prox %>% select(prox.cluster) %>%count() %>% filter(freq>3) %>% select(freq) %>% sum()

cluster.loci <- out.prox %>% filter(prox.cluster == "cluster82") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster110") %>% select(locus_tag) %>% pull() #alginate operon
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster121") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster142") %>% select(locus_tag) %>% pull() #maltose utilisation!

cluster.loci <- out.prox %>% filter(prox.cluster == "cluster16") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster19") %>% select(locus_tag) %>% pull() # --> NUO WHOLE LOT
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster23") %>% select(locus_tag) %>% pull() #maltose utilisation!

cluster.loci <- out.prox %>% filter(prox.cluster == "cluster49") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster50") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster57") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster73") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster74") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster75") %>% select(locus_tag) %>% pull() #maltose utilisation!
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster9") %>% select(locus_tag) %>% pull() #maltose utilisation!


annotation.tbl %>% filter(locus_tag %in% cluster.loci) %>% select(locus_tag, product, COG, gene, kfm_domain, cazy_domain,tcdb_domain) %>% data.frame




#------------------------#

selected.og <- hmg.grp.tbl %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::filter(group == 'cl7') %>%
  dplyr::filter(fdr.p.value_enriched_binary < 0.01) %>% dplyr::select(feature.id) %>% dplyr::pull()

hgm.bin.enriched.cl17 <- annotation.tbl %>%
  dplyr::filter(genome=="Marinobacter_hydrocarbonoclasticus_ATCC_49840") %>%
  dplyr::filter(OG %in% selected.og) %>%
  dplyr::select(locus_tag) %>%
  dplyr::pull()


out.prox <- proximity_clusters(an.tbl, hgm.bin.enriched.cl17,distance=3000 )
out.prox %>% select(prox.cluster) %>%count()

cluster.loci <- out.prox %>% filter(prox.cluster == "cluster56") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster84") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster97") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster103") %>% select(locus_tag) %>% pull() #

annotation.tbl %>% filter(locus_tag %in% cluster.loci) %>% select(locus_tag, product, COG, gene, kfm_domain, cazy_domain,tcdb_domain) %>% data.frame


#------------------------#

selected.og <- hmg.grp.tbl %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::filter(group == 'cl1') %>%
  dplyr::filter(fdr.p.value_enriched_binary < 0.01) %>% dplyr::select(feature.id) %>% dplyr::pull()

hgm.bin.enriched.cl1 <- annotation.tbl %>%
  dplyr::filter(genome=="Marinobacter_adhaerens_HP15") %>%
  dplyr::filter(OG %in% selected.og) %>%
  dplyr::select(locus_tag) %>%
  dplyr::pull()


out.prox <- proximity_clusters(an.tbl, hgm.bin.enriched.cl1,distance=3000 )
out.prox %>% select(prox.cluster) %>%count()

cluster.loci <- out.prox %>% filter(prox.cluster == "cluster108") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster109") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster110") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster34") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster35") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster65") %>% select(locus_tag) %>% pull() #

annotation.tbl %>% filter(locus_tag %in% cluster.loci) %>% select(locus_tag, product, COG, gene, kfm_domain, cazy_domain,tcdb_domain) %>% data.frame


#------------------------#

selected.og <- hmg.grp.tbl %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::filter(group == 'cl3') %>%
  dplyr::filter(fdr.p.value_enriched_binary < 0.05) %>% dplyr::select(feature.id) %>% dplyr::pull()

hgm.bin.enriched.cl3 <- annotation.tbl %>%
  dplyr::filter(genome=="Marinobacter_psychrophilus_20041") %>%
  dplyr::filter(OG %in% selected.og) %>%
  dplyr::select(locus_tag) %>%
  dplyr::pull()


out.prox <- proximity_clusters(an.tbl, hgm.bin.enriched.cl3,distance=3000 )
out.prox %>% select(prox.cluster) %>%count()

cluster.loci <- out.prox %>% filter(prox.cluster == "cluster104") %>% select(locus_tag) %>% pull() # phosphonate operon, C-P lyase operon phn
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster111") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster46") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster62") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster67") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster84") %>% select(locus_tag) %>% pull() #
cluster.loci <- out.prox %>% filter(prox.cluster == "cluster89") %>% select(locus_tag) %>% pull() #

annotation.tbl %>% filter(locus_tag %in% cluster.loci) %>% select(locus_tag, product, COG, gene, kfm_domain, cazy_domain,tcdb_domain) %>% data.frame







# detect_genecluster
#=================================================================================#


detect_genecluster <- function()

# precede
# follow
# distance to nearest

#look for nitrate in product name, and just check if these are in close proximity "clusters" or not :)
# EXAMPLE, EXPLORATION

ltags <- annotation.tbl %>%
  dplyr::filter(genome=="Marinobacter_hydrocarbonoclasticus_ATCC_49840") %>%
  dplyr::filter(grepl("nitrate",product,ignore.case=TRUE)) %>%
  dplyr::select(locus_tag) %>%
  dplyr::pull()

ltags <- annotation.tbl %>%
  dplyr::filter(genome=="Marinobacter_hydrocarbonoclasticus_105B") %>%
  dplyr::filter(grepl("nitrate",product,ignore.case=TRUE)) %>%
  dplyr::select(locus_tag) %>%
  dplyr::pull()

# Explore this
# How many contigs (clusters sadly cant cross multiple contigs)

annotation.tbl %>%
  dplyr::filter( locus_tag %in% ltags) %>%
  dplyr::group_by(seqnames) %>%
  dplyr::count()


# two approaches?
# 1 could be looking at the locus_tag00000 nr and look for close neighbours. say within boundary of x nrs we define it as close proximity
# 2 could be looking at the start and end positions themselves, say within 5000 nucleotides new gene must occur

# Let me explore this in another way

annotation.tbl %>%
  filter(is.na(OG)) %>%
  filter(genome == 'Marinobacter_sp_FDB33') %>%
  select(seqnames, locus_tag, product, type, OG) %>%
  filter(type == 'CDS') %>%
  group_by(seqnames) %>%
  count()

annotation.tbl %>%
  dplyr::filter(is.na(OG)) %>%
  dplyr::filter(genome == 'Marinobacter_sp_FDB33') %>%
  #select(seqnames, locus_tag, product, type, OG) %>%
  dplyr::filter(type == 'CDS') %>%
  dplyr::filter(seqnames =='scaffold15.1')

#Section of genes, just dumped there!
# not homologous to other marinoabacter, yet in a nice cluster, interesting!

annotation.tbl %>%
  dplyr::filter(is.na(OG)) %>%
  dplyr::filter(genome == 'Marinobacter_algicola_DG893') %>%
  dplyr::select(seqnames, locus_tag, product, type, OG) %>%
  dplyr::filter(type == 'CDS') %>%
  dplyr::group_by(seqnames) %>%
  dplyr::filter(seqnames =='NZ_ABCP01000002')






#-------------

#TEST - within genome

#Function input
an.tbl <- annotation.tbl
ltags

#


# proximity_clusters is a function to group locus_tags into clusters based on gene proximity




#-------------------------------------------------------------#






#------------------------------------------------------------------
# tp.ht = hmm.gene %>% arrange(sequence_evalue) %>% slice(1)
# topfeauture <- tp.ht %>% select(domain_name) %>% pull
# starttophit <- tp.ht %>% select(ali_from) %>% pull
# endtophit <-tp.ht %>% select(ali_to) %>% pull
#
# hmm.remove_overlaps <-hmm.gene %>%
#   mutate(tophit = ifelse(domain_name == topfeauture, TRUE, FALSE)) %>%
#   mutate(overlap = ifelse(ali_from <= endtophit & ali_to >= starttophit,TRUE,FALSE)) %>%
#   filter(overlap == TRUE & tophit == TRUE | overlap == FALSE & tophit == FALSE) #%>% select(-tophit,-overlap)
# return(hmm.remove_overlaps)


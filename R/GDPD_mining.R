#
#   LOCUS CENTERED - PARALLEL COORDINATES PLOT and more
#
# to do, need to pass annotation.tbl in functions
# instead of refering to it hard coded...
#


#=======================================================
# Extract synteny blocks...
#n()==1 used to select only single size orgogroups

genoA <- annotation.tbl %>% filter(genome %in% c('Marinobacter_algicola_DG893')) %>% select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_sp_FDB33')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_sp_MCTG268')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())

genoA <- annotation.tbl %>% filter(genome %in% c('Marinobacter_sp_EhN04')) %>% select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_sp_EhC06')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())

genoA <- annotation.tbl %>% filter(genome %in% c('Marinobacter_hydrocarbonoclasticus_ATCC_49840')) %>% select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_adhaerens_HP15')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())

genoA <- annotation.tbl %>% filter(genome %in% c('Marinobacter_hydrocarbonoclasticus_ATCC_49840')) %>% select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_hydrocarbonoclasticus_VT8')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())

genoA <- annotation.tbl %>% filter(genome %in% c('Marinobacter_hydrocarbonoclasticus_ATCC_49840')) %>% select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_psychrophilus_20041')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())
genoB <- annotation.tbl %>% filter(genome %in% c('Marinobacter_hydrocarbonoclasticus_NI9')) %>%  select(OG) %>% group_by(OG) %>% filter(n() == 1) %>% ungroup() %>% dplyr::mutate(n = 1:n())


genoSyn <- genoA %>% dplyr::left_join(genoB, by=c('OG'='OG'))

ggplot(genoSyn,aes(n.x,n.y))+geom_point(shape=21,size=0.5)+fdb_style()






#=================================================================================#
# Functions
#=================================================================================#


#' proximity_clusters is a function to group locus_tags into clusters based on gene proximity
#'
#' @param an.tbl annotation table object
#' @param ltags vector containing locus tags
#' @param distance intergenic distance bewteen two locus_tags
#'
#' @return identified gene clusters
#' @export
#'
#' @examples
#' proximity_clusters(an.tbl, ltags , distance=3000)

# proximity_clusters
#=================================================================================#
proximity_clusters <- function(an.tbl, ltags, distance = 3000){
  f.tbl <- an.tbl %>%
    dplyr::filter(locus_tag %in% ltags)

  contigs <- f.tbl %>%
    dplyr::select(seqnames) %>%
    dplyr::pull() %>%
    as.character() %>%
    unique()

  tracker <- 0 #stores the cluster number
  output <- NULL
  for(contig in contigs){
    #subset table
    sub.tbl <-  f.tbl %>%
      dplyr::filter(seqnames == contig )

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



#find_nearby_genes
#======================================================================#
#' find nearby genes in a single genome, based on proximity cluster of OGs
#' function is handy when you want to find synteny anchors that surround a highly variable region
#' For example, when 10% of genomes have it, regions outside this cluster can be highly conserved.
#'
#' @param OGs
#' @param genome
#' @param ann.tbl
#' @param before
#' @param after
#'
#' @return
#' @export
#'
#' @examples
#' nearby_genes.clust <- find_nearby_genes(OGs = selectedOG, genome ="Marinobacter_algicola_DG893",ann.tbl=annotation.tbl,before=20000, after=20000)
#' nearby_genes.clust$OG

find_nearby_genes <- function(OGs, genome, ann.tbl, before=10000, after=10000){
  genom<-genome

  #Select the locus_tags in the specified genome
  ltags <- ann.tbl %>%
    dplyr::filter(genome==genom) %>%
    dplyr::filter(OG %in% OGs) %>%
    dplyr::select(locus_tag) %>%
    dplyr::pull()

  #detect proximity clusters
  prox.clust <- proximity_clusters(ann.tbl, ltags,distance=3000 )

  #filter largest cluster
  largest.clust <- prox.clust %>%
    dplyr::group_by() %>%
    dplyr::count(prox.cluster) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::select(prox.cluster) %>%
    dplyr::pull() %>% head(n=1)

  prox.clust <- prox.clust %>% filter(prox.cluster==largest.clust)

  #obtain region before and after
  minstart <- prox.clust %>%
    dplyr::select(start) %>%
    min()
  maxend <- prox.clust %>%
    dplyr::select(end) %>%
    max()

  contig <- prox.clust %>%
    dplyr::select(seqnames) %>%
    dplyr::pull() %>%
    as.character() %>%
    unique()

  output <- ann.tbl %>%
    dplyr::filter(genome == genom) %>%
    dplyr::filter(seqnames == contig) %>%
    dplyr::filter(start > minstart - before) %>%
    dplyr::filter(end < maxend + after)

  return(output)
}



#find_nearby_loci
#======================================================================#
#' find nearby genes in a single genome, based on gene proximity to specified locus_tag
#' returns gene neigbourhood in specific genome, on contig of the specified locus_tag
#' mainly helpful when trying to understand what genes are up and downstream of the specified locus
#' @param locus
#' @param genome
#' @param ann.tbl
#' @param before
#' @param after
#'
#' @return
#' @export
#'
#' @examples
find_nearby_loci <- function(locus, genome, ann.tbl, before=10000, after=10000 ){
  genom <- genome
  selected_locus <- ann.tbl %>%
    dplyr::filter(locus_tag ==locus)

  minstart = selected_locus$start
  minend = selected_locus$end
  sqname = selected_locus$seqnames

  output <- annotation.tbl %>%
    dplyr::filter(genome == genom) %>%
    dplyr::filter(seqnames==sqname) %>%
    dplyr::filter(start > minstart - before) %>%
    dplyr::filter(end < minend + after)

  return(output)
}



#fillgaps
#======================================================================#
#' Function that fills in the gaps (genes that were not detected but within were detected in the genomic regions)
#'
#' @param in.cluster
#'
#' @return
#' @export
#'
#' @examples
#' otpt <- fillgaps(in.cluster=output)

fillgaps <- function(in.cluster){
  otpt <- NULL
  for(genom in unique(as.character(in.cluster$genome))){
    sb <- in.cluster %>% filter(genome == genom)
    for(sqnm in unique(as.character(sb$seqnames))){
      minstart <- sb %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::select(start) %>% min()
      maxend <- sb %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::select(end) %>% max()

      filled <- annotation.tbl %>%
        dplyr::filter(genome==genom) %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::filter(start>=minstart & end<=maxend)

      otpt <- rbind(otpt,filled)
    }
  }
  return(otpt)
}



#fillgaps
#======================================================================#
#' Title
#'
#' @param in.cluster
#'
#' @return
#' @export
#'
#' @examples
fillgaps2 <- function(in.cluster){
  otpt <- NULL
  for(genom in unique(as.character(in.cluster$genome))){
    sb <- in.cluster %>%
      dplyr::filter(genome == genom)

    for(sqnm in unique(as.character(sb$seqnames))){
      zb <- sb %>%
        dplyr::filter(seqnames == sqnm)

      for(prx in unique(as.character(zb$prox.cluster))){
      minstart <- zb %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::filter(prox.cluster == prx) %>%
        dplyr::select(start) %>% min()

      maxend <- zb %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::filter(prox.cluster == prx) %>%
        dplyr::select(end) %>%
        max()

      filled <- annotation.tbl %>%
        dplyr::filter(genome==genom) %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::filter(start>=minstart & end<=maxend)

      otpt <- rbind(otpt,filled)
      }
    }
  }
  return(otpt)
}

#expandregion
#======================================================================#
#' Funciton to expand regions, in both directions as specified
#'
#' @param in.cluster
#' @param before
#' @param after
#'
#' @return
#' @export
#'
#' @examples
#' otpt <- expandregion(in.cluster = otpt, before=10000,after=10000)

expandregion <- function(in.cluster, before, after){
  otpt <- NULL
  for(genom in unique(as.character(in.cluster$genome))){
    sb <- in.cluster %>%
      dplyr::filter(genome == genom)

    for(sqnm in unique(as.character(sb$seqnames))){
      sb <- sb %>%
        dplyr::filter(seqnames == sqnm)

      dirctn <- sb %>%
        dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
        dplyr::select(direction) %>%
        sum()

      if(dirctn >0){
        minstart <- sb %>%
          dplyr::filter(seqnames == sqnm) %>%
          dplyr::select(start) %>% min()

        minstart = minstart - before
        maxend <- sb %>%
          dplyr::filter(seqnames == sqnm) %>%
          dplyr::select(end) %>% max()

        maxend = maxend + after
      }else{
        minstart <- sb %>%
          dplyr::filter(seqnames == sqnm) %>%
          dplyr::select(end) %>%
          min()

        minstart = minstart - after

        maxend <- sb %>%
          dplyr::filter(seqnames == sqnm) %>%
          dplyr::select(start) %>%
          max()
        maxend = maxend + before
      }

      added <- annotation.tbl %>%
        dplyr::filter(genome==genom) %>%
        dplyr::filter(seqnames == sqnm) %>%
        filter(start>=minstart & end<=maxend)

      otpt <- rbind(otpt,added)
    }
  }
  return(otpt)
}


#sortregions
#=================================================================================#
#' Function that sets start to 0, and use majority rule direction (most on + strand, reverse visual)
#'
#' @param in.cluster
#' @param majority TRUE when you want to use majority rule
#' @param OGs specify OG strand to define direction
#'
#' @return
#' @export
#'
#' @examples
#' otpt <- sortregions(in.cluster=otpt)

sortregions <- function(in.cluster,majority=TRUE,OGs=NULL){
  otpt <- NULL
  #for each genome
  #genom = unique(as.character(in.cluster$genome))[1]
  for(genom in unique(as.character(in.cluster$genome))){
    sb <- in.cluster %>%
      dplyr::filter(genome == genom)

    #for each contig
    #sqnm <-  unique(as.character(sb$seqnames))[1]
    for(sqnm in unique(as.character(sb$seqnames))){
      sb <- sb %>%
        dplyr::filter(seqnames == sqnm)

      if(!is.null(OGs)){
        #dirctn <- sb %>% mutate(direction=ifelse(strand=="+",1,-1)) %>% select(direction) %>% sum()
        dirctn <- sb %>%
          dplyr::filter(OG == OGs) %>%
          dplyr::select(strand) %>%
          dplyr::pull() %>%
          dplyr::head(n=1) %>%
          dplyr::as.character()

        dirctn <- ifelse(dirctn=='+',1,-1)
      }
      if(majority==TRUE){
        dirctn <- sb %>%
          dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
          dplyr::select(direction) %>%
          sum()
      }

      if(dirctn <0){
        maxend <- sb %>%
          dplyr::filter(seqnames == sqnm) %>%
          dplyr::select(end) %>%
          max()

        sb <- sb %>%
          #dplyr::mutate(start=maxend-start) %>%
          #dplyr::mutate(end=maxend-end) %>%
          #recently added to make sure the genes are in correct orrientation aftr 'rotating'
          dplyr::mutate(end2=maxend-start) %>%
          dplyr::mutate(start2=maxend-end) %>%
          mutate(start=start2,
                 end = end2) %>%
          select(-end2, -start2) %>%
        mutate(strand=ifelse(strand=='-','+','-'))
      }

      #General, all trimm to 0 start
      minstart <- sb %>%
        dplyr::filter(seqnames == sqnm) %>%
        dplyr::select(start) %>% min()

      sb <- sb %>%
        dplyr::mutate(start=start-minstart) %>%
        dplyr::mutate(end=end-minstart)

      otpt <- rbind(otpt,sb)
    }
  }
  return(otpt)
}


#detect_sco_in_region
#=================================================================================#
#' Detect possible positions to center on, detect SCOs in genomic region
#'
#' @param in.cluster
#'
#' @return
#' @export
#'
#' @examples
#' mapOGs <- detect_sco_in_region(clusdetect.sorted)

detect_sco_in_region <- function(in.cluster){
  output = NULL
  nr_of_genomes =  length(unique(as.character(in.cluster$genome)))
  putSCO.OG <- in.cluster %>%
    dplyr::group_by(OG) %>%
    dplyr::select(OG) %>%
    dplyr::count(sort=TRUE) %>%
    dplyr::filter(n==nr_of_genomes) %>%
    dplyr::select(OG) %>%
    dplyr::pull()

  for(og in putSCO.OG){
    abundancecheck <- in.cluster %>%
      dplyr::filter(OG == og) %>%
      dplyr::group_by(OG,genome) %>%
      dplyr::count(sort=TRUE)

    if(nrow(abundancecheck) == nr_of_genomes){
      #abundancecheck is king
      message(og)
      output = c(output, og)
    }else{
      message('not present in all')
      #non shared amongst all...
    }
  }
  return(output)
}


#centerregions
#=================================================================================#
#' Function that centers regions at the specified OG (needs to be present in all!)
#' CONSIDER using our function detect_sco_in_region() to explore whether you have a valid SCO in the region
#' if you dont have a SCO in this specific region, you could supply a vector to OGs (instead of a single value)
#' in case of a vector, the function will use all orthologs detected in given genome, and center based on distance in nucleotide
#' @param in.cluster
#' @param OGs either one (center on this one), or a vector of multiple OGs (center)
#'
#' @return
#' @export
#'
#' @examples
#' clusdetect.center <- centerregions(in.cluster=clusdetect.expanded, OGs=mapOG)
#
#     CHECK WHETHER THIS IS A VALID SCO!
#

centerregions <- function(in.cluster, OGs){
  otpt <- NULL
  for(genom in unique(as.character(in.cluster$genome))){
    sb <- in.cluster %>%
      dplyr::filter(genome == genom)

    for(sqnm in unique(as.character(sb$seqnames))){
      sb <- sb %>%
        dplyr::filter(seqnames == sqnm)

      #if you only supply one, the direction of the gene is assesed and start position taken as 0
      if(length(OGs)==1){
        # could include a check whether this is a SCO, as we need it in all genomes for this procedure
        direction <- sb %>%
          dplyr::filter(OG==OGs) %>%
          dplyr::select(strand) %>%
          dplyr::pull() %>%
          head(n=1) %>%
          as.character()

        if(length(direction)>0){
          if(direction == '-'){
            minstart <- sb %>%
              dplyr::filter(OG==OGs) %>%
              dplyr::select(start) %>%
              dplyr::pull() %>% head(n=1)

          }else{
            minstart <- sb %>%
              dplyr::filter(OG==OGs) %>%
              dplyr::select(end) %>%
              dplyr::pull() %>%
              head(n=1)
          }
        }else{
          minstart=0
        }

      #if you have more than one OG selected, it will look for the harmonic center of the regions
      }else{
        #not sure if direction matters, need to play more
        #dirctn <- sb %>% mutate(direction=ifelse(strand=="+",1,-1)) %>% select(direction) %>% sum()
        mstart <- sb %>%
          dplyr::filter(OG%in%OGs) %>%
          dplyr::select(start) %>%
          dplyr::pull() %>%
          min()

        mend <- sb %>%
          dplyr::filter(OG%in%OGs) %>%
          dplyr::select(end) %>%
          dplyr::pull() %>%
          max()

        minstart = round(((mend-mstart)/2),0)
      }

      #in both cases, we continue
      sb <- sb %>%
        dplyr::mutate(start=start-minstart) %>%
        dplyr::mutate(end=end-minstart)

      otpt <- rbind(otpt,sb)
    }
  }
  return(otpt)
}


#multigenedetect
#=================================================================================#
#' detect proximity based gene clusters in sets of Orthologous genes in multiple genomes
#'
#' @param OGs
#' @param genomes
#' @param largest.cluster
#'
#' @return
#' @export
#'
#' @examples
multigenedetect <- function(OGs, genomes,largest.cluster=TRUE){
  output <- NULL
  for(genom in genomes){
    message(genom)
    ltags <- annotation.tbl %>%
      dplyr::filter(OG %in% OGs) %>%
      dplyr::filter(genome == genom) %>%
      dplyr::select(locus_tag) %>%
      dplyr::pull()

    if(length(ltags)>0){
      prox.clust <- proximity_clusters(annotation.tbl, ltags,distance=3000 )

      largest.clust <- prox.clust %>%
        dplyr::group_by() %>%
        dplyr::count(prox.cluster) %>%
        dplyr::filter(n == max(n)) %>%
        dplyr::select(prox.cluster) %>%
        dplyr::pull() %>%
        head(n=1)

      #questioning the top_n, as somethimes it might be relevant to look into them both
      filtered <- prox.clust %>%
        dplyr::filter(prox.cluster == largest.clust)

      output<- rbind(output,filtered )
    }
  }
  return(output)
}





#genome_dummy_genes
#=================================================================================#
#' Function to add dummy values for additional genomes, without orthologs
#' this is sometimes requied, for example when plotting operons against phylogeny
#' @param in.cluster gene cluster object
#' @param genomes vector of genonomes (for example, all genomes in phylogeny)
#'
#' @return gene cluster object with dummies where no OGs were identified
#' @export
#'
#' @examples
genome_dummy_genes <- function(in.cluster, genomes){

  #genomes in incluster
  pg <- in.cluster %>%
    dplyr::select(genome) %>%
    unique() %>%
    dplyr::pull() %>%
    as.character()


  for(gn in setdiff(genomes, pg)){

    tmp <- annotation.tbl %>%
      dplyr::filter(genome==gn) %>%
      dplyr::slicee(1) %>%
      dplyr::mutate(start=NA,
                    end=NA,
                    width=NA,
                    source=NA,
                    type=NA,
                    score=NA,
                    phase=NA,
                    ID=NA,
                    Name =NA,
                    COG=NA,
                    gene=NA,
                    locus_tag=NA,
                    product=NA,
                    note=NA,
                    eC_number=NA,
                    rpt_family=NA,
                    OG=NA,
                    kfm_evalue=NA,
                    kfm_domain=NA,
                    cazy_evalue=NA,
                    cazy_domain=NA,
                    tcdb_evalue=NA,
                    tcdb_domain=NA)

    in.cluster <- rbind(in.cluster, tmp)
  }
  return(in.cluster)
}




#=================================================================================#













#some T6SS for carl
selected_locus <- annotation.tbl %>%
  filter(genome=='Marinobacter_algicola_DG893') %>%
  filter(grepl('vgrg',product,ignore.case=T)) %>%
  select(locus_tag) %>%
  pull()


out <- NULL
for(i in selected_locus){
  tmp <- find_nearby_loci(locus=i,
                          genome='Marinobacter_algicola_DG893',
                          ann.tbl=annotation.tbl,
                          before=20000,
                          after=20000)
  out <- rbind(out, tmp)
  }



annotated.out <- out %>% left_join(KO_annotated %>% select(-genome), by='locus_tag') %>% select(-eC_number, -kfm_evalue, -kfm_domain) %>% unique()

ag1<-  aggregate(KO~locus_tag, data = annotated.out, paste0, collapse=";")
ag2 <-  aggregate(MOD_F~locus_tag, data = annotated.out, paste0, collapse=";")
ag3 <-  aggregate(TYPE ~locus_tag, data = annotated.out, paste0, collapse=";")
ag4 <- aggregate(DESCRIPTION ~locus_tag, data = annotated.out, paste0, collapse=";")

annotated.out <- out %>% left_join(ag1, by="locus_tag") %>%
  left_join(ag2, by="locus_tag") %>%
  left_join(ag3, by="locus_tag") %>%
  left_join(ag4, by="locus_tag") %>%
  select(-eC_number, -kfm_evalue, -kfm_domain)

write.table(annotated.out, '~/DATA/MarbGenomics/EImodules_DG893.tsv',sep='\t',quote = FALSE, row.names = FALSE)
write.table(annotated.out %>% filter(seqnames=="NZ_ABCP01000034"), '~/DATA/MarbGenomics/T6SS_DG893.tsv',sep='\t',quote = FALSE, row.names = FALSE)




p1 <- annotated.out%>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  ggplot(aes(xmin = start, xmax = end, y = strand, fill = KO,forward=direction,label=Name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  facet_wrap(~ seqnames, scales = "free", ncol = 1) + geom_gene_label(align = "left") + theme(legend.position = 'none')+
  scale_fill_manual(values=colorRampPalette(brewer.pal(8, "Set2"))(43), na.value="white")

ggsave('~/DATA/MarbGenomics/Graphs/EImodules_DG893.pdf',plot=p1, width = 15, height = 7,unit='cm')


p1 <- annotated.out%>% filter(seqnames=="NZ_ABCP01000034")%>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  ggplot(aes(xmin = start, xmax = end, y = strand, fill = KO,forward=direction,label=Name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  facet_wrap(~ seqnames, scales = "free", ncol = 1) + geom_gene_label(align = "left") + theme(legend.position = 'none')+
  scale_fill_manual(values=colorRampPalette(brewer.pal(8, "Set2"))(43), na.value="white")

ggsave('~/DATA/MarbGenomics/Graphs/T6SS_DG893.pdf',plot=p1, width = 15, height = 3,unit='cm')





out2 <- NULL
for(i in selected_locus[3]){
  tmp <- find_nearby_loci(locus=i,
                          genome='Marinobacter_algicola_DG893',
                          ann.tbl=annotation.tbl,
                          before=20000,
                          after=50000)
  out2 <- rbind(out2, tmp)
}




#=========================

selected_locus <- KO_annotated %>%
  filter(grepl('urea',ko_annotation ,ignore.case=T)) %>%
  group_by(ko_annotation ,DESCRIPTION) %>%
  filter(grepl('VT8',genome)) %>% select(locus_tag) %>%
  pull()


out <- NULL
for(i in selected_locus){
  tmp <- find_nearby_loci(locus=i,
                          genome='Marinobacter_hydrocarbonoclasticus_VT8',
                          ann.tbl=annotation.tbl,
                          before=2000,
                          after=2000)
  out <- rbind(out, tmp)
}

out <- fillgaps(out)


annotated.out <-
  out %>%
  left_join(KO_annotated %>%
  select(-genome), by='locus_tag') %>%
  select(-eC_number, -kfm_evalue, -kfm_domain) %>%
  unique()


ag1<-  aggregate(KO~locus_tag, data = annotated.out, paste0, collapse=";")
ag2 <-  aggregate(MOD_F~locus_tag, data = annotated.out, paste0, collapse=";")
ag3 <-  aggregate(TYPE ~locus_tag, data = annotated.out, paste0, collapse=";")
ag4 <- aggregate(DESCRIPTION ~locus_tag, data = annotated.out, paste0, collapse=";")

annotated.out <- out %>% left_join(ag1, by="locus_tag") %>%
  left_join(ag2, by="locus_tag") %>%
  left_join(ag3, by="locus_tag") %>%
  left_join(ag4, by="locus_tag") %>%
  select(-eC_number, -kfm_evalue, -kfm_domain)


annotated.out %>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  ggplot(aes(xmin = start, xmax = end, y = strand, fill = KO,forward=direction,label=Name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  facet_wrap(~ seqnames, scales = "free", ncol = 1) + geom_gene_label(align = "left") + theme(legend.position = 'none')+
  scale_fill_manual(values=colorRampPalette(brewer.pal(8, "Set2"))(43), na.value="white")



selectedOG <- annotated.out %>% select(OG) %>% pull()


#RUN
clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=1000,after=1000)
#clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded)
#clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)

mapOGs <- detect_sco_in_region(clusdetect.sorted)
mapOG = mapOGs[1]
if(length(mapOG>0)){
  clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=FALSE, OGs=mapOG)
  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)
}else{

}



annotated.out <-
  clusdetect.expanded %>%
  left_join(KO_annotated %>%
              select(-genome), by='locus_tag') %>%
  select(-eC_number, -kfm_evalue, -kfm_domain) %>%
  unique()


annotated.out %>%
  filter(grepl('VT8',genome))%>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  ggplot(aes(xmin = start, xmax = end, y = strand,forward=direction,label=Name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  facet_wrap(~genome, scales = "free", ncol = 1)  + theme(legend.position = "none")



+ geom_gene_label(align = "left") + theme(legend.position = 'none')+
  scale_fill_manual(values=colorRampPalette(brewer.pal(8, "Set2"))(43), na.value="white")










#===========

#GVD
# can we detect nearest conserved region at both ends? and merge that as input for cluster detection

nearby_genes.clust <- find_nearby_loci(locus=selected_locus,
                                       genome='Marinobacter_algicola_DG893',
                                       ann.tbl=annotation.tbl,
                                       before=10000,
                                       after=10000)

selectedOG <- nearby_genes.clust$OG
clusdetect.loc <- multigenedetect(OG=selectedOG, genomes='Marinobacter_algicola_DG893')



selected_locus <- annotation.tbl %>% filter(COG=='COG0722')%>% filter(genome=='Marinobacter_algicola_DG893') %>% select(locus_tag) %>% pull

nearby_genes.clust <- find_nearby_loci(locus=selected_locus[1],
                                       genome='Marinobacter_algicola_DG893',
                                       ann.tbl=annotation.tbl,
                                       before=10000,
                                       after=10000)
selectedOG <- nearby_genes.clust %>% select(OG) %>% pull()
outname <- 'COG0722_1'

nearby_genes.clust <- find_nearby_loci(locus=selected_locus[2],
                                       genome='Marinobacter_algicola_DG893',
                                       ann.tbl=annotation.tbl,
                                       before=10000,
                                       after=10000)
selectedOG <- nearby_genes.clust %>% select(OG) %>% pull()
outname <- 'COG0722_2'
focusOG <- selected_locus[2]

#COG0147
selected_locus <- annotation.tbl %>% filter(COG=='COG0147')%>% filter(genome=='Marinobacter_algicola_DG893') %>% select(locus_tag) %>% pull

nearby_genes.clust <- find_nearby_loci(locus=selected_locus[1],
                                       genome='Marinobacter_algicola_DG893',
                                       ann.tbl=annotation.tbl,
                                       before=10000,
                                       after=10000)
selectedOG <- nearby_genes.clust %>% select(OG) %>% pull()
outname <- 'COG0147_1'
focusOG <- selected_locus[1]

nearby_genes.clust <- find_nearby_loci(locus=selected_locus[2],
                                       genome='Marinobacter_algicola_DG893',
                                       ann.tbl=annotation.tbl,
                                       before=10000,
                                       after=10000)
selectedOG <- nearby_genes.clust %>% select(OG) %>% pull()
outname <- 'COG0147_2'
focusOG <- selected_locus[2]




#=================================================================================#
# PIPELINE

#SET GENOMES
selgenomes <- tree$tip.label
selgenomes <- tr.tree$tip.label
tree=tr.tree

#SET OGs
selectedOG <- sig.genes2
selectedOG <- c('OG0002720','OG0002559','OG0002698','OG0002906',
                'OG0002618','OG0002821','OG0001780','OG0002430',
                'OG0001877','OG0001932','OG0002209','OG0002331',
                'OG0001842','OG0001842','OG0002232','OG0001733',
                'OG0002067')

#RUN
clusdetect.og <- multigenedetect(OG=selectedOG, genomes=selgenomes)
clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=40000,after=40000)
#clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded)
#clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)

mapOGs <- detect_sco_in_region(clusdetect.sorted)
mapOG = mapOGs[1]
if(length(mapOG>0)){
  clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=FALSE, OGs=mapOG)
  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)
  }else{

}

#=================================================================================#

selectedOG <- sco %>% as.character()
outdir <- '~/DATA/MarbGenomics/Graphs/'
outname <- 'hyperconserved'


#=========
possibleSCOs = setdiff(sco, unique(plot.data$OG))
annotation.tbl %>% filter(OG  %in% possibleSCOs) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(locus_tag, product,Name,COG, OG)

#nr_of_genomes =  length(unique(as.character(in.cluster$genome)))
#putSCO.OG <- in.cluster %>% group_by(OG) %>% select(OG) %>% dplyr::count(sort=TRUE) %>% filter(n==nr_of_genomes) %>% select(OG) %>% pull()


selected_locus = 'NHGOFBLF_02968'
outname <- 'murD'

#-----
possibleSCOs = setdiff(possibleSCOs, unique(plot.data$OG))
annotation.tbl %>% filter(OG  %in% possibleSCOs) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(locus_tag, product,Name,COG, OG)

outname <- 'hslU'
selected_locus = 'NHGOFBLF_03018'



#----- look for phenol degradation
# BASED ON ONE HIT!? -- one gene retrieved, passed to selected_locus
# In this case, inspection shows that there is only one hit!
# OG0003741
selected_locus <- annotation.tbl %>% filter(grepl('catechol',product,ignore.case=T)) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(locus_tag) %>% pull()
outdir <- '~/DATA/MarbGenomics/Graphs/'
outname <- 'catechol'

#of alginate lyase test
selected_locus= 'NHGOFBLF_01686'
outdir <- '~/DATA/MarbGenomics/Graphs/'
outname <- 'alginate'




#==========

#we can use find_nearby_loci, to inspect the neighbourhood of a specific locus_tag
nearby_genes.clust <- find_nearby_loci(locus=selected_locus,
                                       genome='Marinobacter_algicola_DG893',
                                       ann.tbl=annotation.tbl,
                                       before=10000,
                                       after=10000)

selectedOG <- nearby_genes.clust$OG

#---- FOR UNCONSERVED REGIONS, SAY LOSS/AQUISITION of LARGE OPERON, alginate, glycogen etc...
# the idea is to iteratively optimise the detection
# first, we obtain gene cluster in the organism that has it
clusdetect.loc <- multigenedetect(OG=selectedOG, genomes='Marinobacter_algicola_DG893')
clusdetect.wider <- expandregion(in.cluster = clusdetect.loc, before=10000,after=10000)

selectedOG <- clusdetect.wider$OG

#remove potential duplicates, and NAs
selectedOG <- unique(selectedOG[!is.na(selectedOG)])


clusdetect.og <- multigenedetect(OG=selectedOG, genomes=selgenomes)
clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=20000,after=20000)
clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)

mapOGs <- detect_sco_in_region(clusdetect.sorted)
mapOG = mapOGs[1]

if(is.null(mapOG)){
  #clusterregions detects lenght of OGs, if 1, use that (in case of SCO), if more, use all.
  message('-- no SCO detected, average center')
  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)

}else{
  message('SCO detected, center on first SCO!')
  clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=FALSE, OGs = mapOG)
  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)
}

clusdetect.center$genome <- factor(clusdetect.center$genome,levels = rev(tree$tip.label))
focusOG=selectedOG

out.p <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, color = SELOG,forward=direction,label=Name,fill=OG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
  scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
    theme(legend.position = "none")

#out.p
outfile <- paste0(outname,'_ordered.pdf')
outpath <- paste0(outdir,outfile)

ggsave(outpath,plot=out.p, width = 65, height = 65,unit='cm')


#focusOG = selectedOG
#plot.data <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
#  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>% mutate(id=genome) %>%
#  select(id, everything()) %>% mutate(gene=OG)%>% data.frame()

#rownames(plot.data) <- plot.data$genome


#p <- ggtree(tree, branch.length='none') +
#  geom_tiplab() +
#  ggtree::geom_facet(
#        mapping = aes(xmin = start, xmax = end),
#        data = plot.data,
#        geom = geom_gene_arrow,
#        panel = 'Alignment'
#             #arrowhead_height = unit(3, "mm"),
#             # = unit(2, "mm"),
#             ) +
#  theme(legend.position = "none")+
#  #scale_fill_brewer(palette = "Set3") +
#  #scale_x_continuous(expand=c(0,0)) +
#  theme(strip.text=element_blank(),
#        panel.spacing=unit(0, 'cm'))
#out.p <- facet_widths(p, widths=c(1,4))

#ggsave('~/DATA/MarbGenomics/Graphs/catechol_phylo.pdf',plot=out.p, width = 65, height = 65,unit='cm')


# Alternative to above masked , idially we find a solution to using geom_gene_arrow in facet_plot!
#----- using geom_motif, not ideal, need to specify one common gene, else region gets dropped
addOGs = selectedOG
#focusOG = selectedOG[2]

#ANCHOR DETECTION
#
anchors <- clusdetect.center %>% dplyr::group_by(genome,seqnames) %>% dplyr::slice(which.min(abs(start))) %>% ungroup() %>% select(locus_tag) %>% pull()



plot.data <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(locus_tag%in%anchors,'anchor',ifelse(OG%in%focusOG,'focal',ifelse(OG%in%addOGs,'t6ss','z')))) %>%
  mutate(id=genome) %>%
  select(id, everything()) %>% mutate(gene=OG)%>% data.frame()

g = ggtree(tree) + geom_tiplab(align=TRUE)
p <- facet_plot(g,
                panel='alignment',
                mapping = aes(xmin = start,
                              xmax = end,
                              y = genome,
                              color = OG,
                              forward=direction,
  ),
                data=plot.data,
                on='anchor',
                geom=geom_motif,
                arrowhead_height = unit(3, "mm"),
                arrowhead_width = unit(2, "mm"))+
  theme(legend.position = "none")+scale_fill_manual(values=c('blue','red','black','white'))+
  theme(strip.text=element_blank(),panel.spacing=unit(0, 'cm'))

gt = ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
#gtable_show_layout(gt) # will show you the layout - very handy function
#gt # see plot layout in table format
#gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[5] = 0.25*gt$widths[5] # in this case it was colmun 7 - reduce the width by a half
#out.p <- grid::grid.draw(gt) # plot with grid draw

outfile <- paste0(outname,'_phylo.pdf')
outpath <- paste0(outdir,outfile)
ggsave(outpath,plot=gt, width = 65, height = 65,unit='cm')














annotation.tbl %>%
  filter(genome=='Marinobacter_algicola_DG893') %>%
  slice(1) %>%
  mutate(start=NA, end=NA,width=NA,source=NA)






#-----  look for NUO
# This is an example where the region in absent in most genomes,
# first detect the region, and extract relevant OGs


#selgenomes <- metadata %>% filter(group %in% c('cl14'))  %>% select(genome) %>% pull()


selected_gene_name <- 'nuo2'
selectedOG = annotation.tbl %>% filter(grepl('nuo',Name,ignore.case=T)) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(OG) %>% pull()
selectedOG = annotation.tbl %>% filter(grepl('nuo',Name,ignore.case=T)) %>% filter(genome=='Marinobacter_sp_FDB33') %>% select(OG) %>% pull()
nearby_genes.clust <- find_nearby_genes(OGs = selectedOG, genome ="Marinobacter_sp_FDB33",ann.tbl=annotation.tbl,before=5000, after=5000)


selected_gene_name <- 'urea'
selectedOG = annotation.tbl %>% filter(grepl('urea',product,ignore.case=T)) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(OG, seqnames, start,end,locus_tag,Name, product) %>% select(OG) %>% pull()


selected_gene_name <- 'flg'
selectedOG = annotation.tbl %>% filter(grepl('flg',product,ignore.case=T)) %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(OG, seqnames, start,end,locus_tag,Name, product)   %>% select(OG) %>% pull()

# Use function to retrieve nearby genes
#nearby_genes.clust <- find_nearby_genes(OGs = selectedOG, genome ="Marinobacter_algicola_DG893",ann.tbl=annotation.tbl,before=20000, after=20000)
nearby_genes.clust <- find_nearby_genes(OGs = selectedOG, genome ="Marinobacter_algicola_DG893",ann.tbl=annotation.tbl,before=0, after=0)

#include the larger region into selectedOGs
selectedOG <- nearby_genes.clust$OG

mapOG = "OG0003273"

clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)


clusdetect.og <-

clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)
#clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)


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

#out.p <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
out.p <- clusdetect.center %>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=Name,fill=OG)) + #color=SELCOG
  geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
  scale_fill_viridis(discrete = TRUE)+
  #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
  theme(legend.position = "none")

out.p

ggsave(paste0('~/DATA/MarbGenomics/Graphs/NUO_region_phylo_sorted',selected_gene_name,'.pdf'),plot=out.p, width = 20, height = 35,unit='cm')





selected.loci_all <- KO_annotated %>% filter(grepl('nuo',ko_annotation,ignore.case=T)) %>% select(locus_tag) %>% pull()
clusdetect.og <- annotation.tbl %>% filter(locus_tag %in% selected.loci_all)
clusdetect.og <- clusdetect.og %>% filter(genome=='Marinobacter_sp_FDB33')

clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded,majority = TRUE)
#clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)


clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)

out.p <- clusdetect.filled %>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  #mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=Name,fill=OG)) + #color=SELCOG
  geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
  scale_fill_viridis(discrete = TRUE)+
  #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
  theme(legend.position = "none")
ggsave(paste0('~/DATA/MarbGenomics/Graphs/kegg_based_test_region_phylo_sorted',selected_gene_name,'.pdf'),plot=out.p, width = 20, height = 35,unit='cm')




#BCBNPALP_05983

# GDPD in HP15
#     MLBPJMHG_03489
#
# OG::: OG0002720
#

annotation.tbl %>% filter(locus_tag %in% c('MLBPJMHG_03481','MLBPJMHG_03482','MLBPJMHG_03483','MLBPJMHG_03484','MLBPJMHG_03485','MLBPJMHG_03486','MLBPJMHG_03487','MLBPJMHG_03488','MLBPJMHG_03489','MLBPJMHG_03490','MLBPJMHG_03491','MLBPJMHG_03492','MLBPJMHG_03493','MLBPJMHG_03494','MLBPJMHG_03495','MLBPJMHG_03496','MLBPJMHG_03497','MLBPJMHG_03498','MLBPJMHG_03499','MLBPJMHG_03500','MLBPJMHG_03501')) %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain) %>% select(OG) %>%pull()

selectedOG <- annotation.tbl %>% filter(locus_tag %in% c('MLBPJMHG_03481','MLBPJMHG_03482','MLBPJMHG_03483','MLBPJMHG_03484','MLBPJMHG_03485','MLBPJMHG_03486','MLBPJMHG_03487','MLBPJMHG_03488','MLBPJMHG_03489','MLBPJMHG_03490','MLBPJMHG_03491','MLBPJMHG_03492','MLBPJMHG_03493','MLBPJMHG_03494','MLBPJMHG_03495','MLBPJMHG_03496','MLBPJMHG_03497','MLBPJMHG_03498','MLBPJMHG_03499','MLBPJMHG_03500','MLBPJMHG_03501','MLBPJMHG_03502','MLBPJMHG_03503','MLBPJMHG_03504','MLBPJMHG_03505')) %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain) %>% select(OG) %>%pull()


selectedOG <- annotation.tbl %>% filter(locus_tag %in% c('MLBPJMHG_03488','MLBPJMHG_03489','MLBPJMHG_03490','MLBPJMHG_03491','MLBPJMHG_03492')) %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain) %>% select(OG) %>%pull()

#selectedOG <- annotation.tbl %>% filter(locus_tag %in% c('MLBPJMHG_03489')) %>%
#  select(locus_tag, genome, product, OG, COG, kfm_domain) %>% select(OG) %>%pull()

annotation.tbl %>% filter(OG %in% selectedOG) %>% filter(genome=='Marinobacter_algicola_DG893') %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain)

annotation.tbl %>% filter(OG == 'OG0002720') %>% filter(genome=='Marinobacter_sp_FDB33') %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain)

annotation.tbl %>% filter(kfm_domain == 'K01126') %>% filter(genome=='Marinobacter_sp_FDB33') %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain)

annotation.tbl %>% filter(kfm_domain == 'K01851') %>% filter(genome=='Marinobacter_sp_FDB33') %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain)

annotation.tbl %>% filter(kfm_domain == 'K04781') %>% filter(genome=='Marinobacter_sp_FDB33') %>%
  select(locus_tag, genome, product, OG, COG, kfm_domain)


annotation.tbl %>% filter(OG %in% selectedOG) %>% filter(genome %in% c('Marinobacter_algicola_DG893','Marinobacter_adhaerens_HP15')) %>%
  ggplot(aes(x=start, y=seqnames))+geom_point()






# QUICK CODE RUN, IF FUNCTIONS ARE SET UP
selgenomes <- c('Marinobacter_flavimaris_LMG_23834','Marinobacter_hydrocarbonoclasticus_VT8','Marinobacter_sp_DSM_26671','Marinobacter_salarius_HK15','Marinobacter_adhaerens_HP15','Marinobacter_algicola_DG893','Marinobacter_sp_C18','Marinobacter_sp_FDB33')
selgenomes <-as.character(metadata$genome)

selgenomes <- as.character(metadata$genome[1:30])
selgenomes <- metadata %>% filter(group %in% c('cl14','cl7','cl1','cl10','cl19'))  %>% select(genome) %>% pull()

# WTF IS THIS ONE FFS
#selgenomes <-  c('Marinobacter_flavimaris_LMG_23834','Marinobacter_hydrocarbonoclasticus_VT8','Marinobacter_sp_DSM_26671','Marinobacter_salarius_HK15','Marinobacter_adhaerens_HP15','Marinobacter_algicola_DG893','Marinobacter_sp_C18','Marinobacter_sp_FDB33','Marinobacter_sp_ELB17')

selgenomes <-  c('Marinobacter_flavimaris_LMG_23834','Marinobacter_hydrocarbonoclasticus_VT8','Marinobacter_sp_DSM_26671','Marinobacter_salarius_HK15','Marinobacter_adhaerens_HP15','Marinobacter_algicola_DG893','Marinobacter_sp_C18','Marinobacter_sp_FDB33','Marinobacter_adhaerens_strain_PBVC038','Marinobacter_flavimaris_KCTC_12185','Marinobacter_halophilus_JCM_30472','Marinobacter_lutaoensis_T5054')

selgenomes <-  c('Marinobacter_lutaoensis_T5054','Marinobacter_sp_ELB17')
mapOG = 'OG0002720'


#detect, fill, expand and sort
clusdetect.og <- multigenedetect(OG=selectedOG, genomes=selgenomes)
clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=40000,after=40000)
clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded)
clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)


#Visualise
#clusdetect.sorted %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
#  ggplot(aes(xmin = start, xmax = end, y = genome, fill = OG,forward=direction,label=Name)) +
#  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")

#clusdetect.sorted %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
#  mutate(SELOG = ifelse(OG%in%mapOG,TRUE,FALSE)) %>%
#  ggplot(aes(xmin = start, xmax = end, y = genome, fill = OG,forward=direction,label=Name,color=SELOG)) +
#  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+scale_color_manual(values=c('white','red'))+fdb_style(aspect.ratio=0.5) +
#  theme(legend.position = "none")

#theme(legend.position="bottom",
#      legend.spacing.x = unit(0, 'cm')

clusdetect.center$genome <- factor(clusdetect.center$genome,levels = tree$tip.label)

clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,TRUE,FALSE)) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = OG,forward=direction,label=Name,color=SELOG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+scale_color_manual(values=c('white','red'))+fdb_style(aspect.ratio=0.5) +
  theme(legend.position = "none")

levels(clusdetect.center$genome)
tree$tip.label

library(gggenes)

output %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  ggplot(aes(xmin = start, xmax = end, y = strand, fill = OG,forward=direction,label=Name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  facet_wrap(~ genome, scales = "free", ncol = 1) + geom_gene_label(align = "left")




clusdetect.sorted$genome <- factor(clusdetect.sorted$genome,levels = tree$tip.label)

clusdetect.sorted %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,TRUE,FALSE)) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = OG,forward=direction,label=Name,color=SELOG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+scale_color_manual(values=c('white','red'))+fdb_style(aspect.ratio=0.5) +
  theme(legend.position = "none")



selectedOG <- c('OG0002720','OG0002559','OG0002698','OG0002906',
                'OG0002618','OG0002821','OG0001780','OG0002430',
                'OG0001877','OG0001932','OG0002209','OG0002331',
                'OG0001842','OG0001842','OG0002232','OG0001733',
                'OG0002067')

selgenomes <-  c('Marinobacter_lutaoensis_T5054','Marinobacter_sp_ELB17')
mapOG = 'OG0002906'
selgenomes <- as.character(metadata$genome[1:30])
selgenomes <- tree$tip.label


#detect, fill, expand and sort
clusdetect.og <- multigenedetect(OG=selectedOG, genomes=selgenomes)
clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=40000,after=40000)
clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=TRUE)



#detect possible anchors
mapOGs <- detect_sco_in_region(clusdetect.sorted)
mapOG = mapOGs[1]

clusdetect.sorted <- sortregions(in.cluster=clusdetect.expanded, majority=FALSE, OGs=mapOGs)
clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)
clusdetect.center$genome <- factor(clusdetect.center$genome,levels = rev(tree$tip.label))

out.p <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,TRUE,FALSE)) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = OG,forward=direction,label=Name,color=SELOG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+scale_color_manual(values=c('white','red'))+fdb_style(aspect.ratio=1) +
  theme(legend.position = "none")


ggsave('~/DATA/MarbGenomics/Graphs/GDPD_region_phylo_sorted.pdf',plot=out.p, width = 65, height = 65,unit='cm')


focusOG = c('OG0002720')
addOGs = clusdetect.center %>% filter(genome=='Marinobacter_hydrocarbonoclasticus_ATCC_49840') %>% select(locus_tag, product, kfm_domain, OG) %>% data.frame %>% tail(20) %>% select(OG) %>% pull()

plot.data <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal',ifelse(OG%in%addOGs,'t6ss','z')))) %>% mutate(id=genome) %>%
  select(id, everything()) %>% mutate(gene=OG)%>% data.frame()

g = ggtree(tree) + geom_tiplab(align=TRUE)
p <- facet_plot(g,
                panel='alignment',
                mapping = aes(xmin = start,
                              xmax = end,
                              y = genome,
                              color = OG,
                              forward=direction,
                              #label=Name,
                              fill=SELOG),
                data=plot.data,
                on='anchor',
                geom=geom_motif)+
  theme(legend.position = "none")+scale_fill_manual(values=c('blue','red','black','white'))+
  theme(strip.text=element_blank(),panel.spacing=unit(0, 'cm'))

gt = ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
#gtable_show_layout(gt) # will show you the layout - very handy function
#gt # see plot layout in table format
#gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[5] = 0.25*gt$widths[5] # in this case it was colmun 7 - reduce the width by a half
#out.p <- grid::grid.draw(gt) # plot with grid draw
ggsave('~/DATA/MarbGenomics/Graphs/GDPD_phylo.pdf',plot=gt, width = 65, height = 65,unit='cm')









focusOG = c('OG0002720')
addOGs = clusdetect.center %>% filter(genome=='Marinobacter_hydrocarbonoclasticus_ATCC_49840') %>% select(locus_tag, product, kfm_domain, OG) %>% data.frame %>% tail(20) %>% select(OG) %>% pull()

out.p <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal',ifelse(OG%in%addOGs,'t6ss','z')))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = OG,forward=direction,label=Name,color=SELOG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
  scale_color_manual(values=c('blue','red','black','white'))+fdb_style(aspect.ratio=1) +
  theme(legend.position = "none")



out.p <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal',ifelse(OG%in%addOGs,'t6ss','z')))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = SELOG,forward=direction,label=Name,color=OG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
  scale_fill_manual(values=c('blue','red','black','white'))+fdb_style(aspect.ratio=1) +
  theme(legend.position = "none")

ggsave('~/DATA/MarbGenomics/Graphs/GDPD_region_phylo_sorted_2.pdf',plot=out.p, width = 65, height = 65,unit='cm')


out.p <- clusdetect.center %>% mutate(direction=ifelse(strand=="+",1,-1)) %>%
  mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal',ifelse(OG%in%addOGs,'t6ss','z')))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = kfm_domain,forward=direction,label=Name,color=SELOG)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
  scale_colour_manual(values=c('blue','red','black','white'))+fdb_style(aspect.ratio=1) +
  theme(legend.position = "none")

ggsave('~/DATA/MarbGenomics/Graphs/GDPD_region_phylo_sorted_4.pdf',plot=out.p, width = 65, height = 65,unit='cm')








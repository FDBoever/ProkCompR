#### FUNCTIONS ####
#function to vonvert ncbi locus_tags to prokka locus_tags
ltag_ncbi2prokka <- function(prokka2ncbi.tbl, loci){
  message('translating ',length(loci), ' loci from ncbi to prokka')
  ncbi.selected <- prokka2ncbi.tbl %>%
    filter(grepl(paste(loci,collapse='|'),ncbi_old_locus_tag))

  prokka_loci <- ncbi.selected %>%
    select(prokka_locus_tag) %>% pull() %>% as.character()
  message('revovered ',length(prokka_loci),'/',length(loci), ' prokka loci')

  return(prokka_loci)
}

#example
#publ_loci <- c('D777_03329','D777_03330','D777_03331','D777_03332','D777_03333','D777_03334','D777_03335','D777_03336','D777_03337')
#ltag_ncbi2prokka(prokka2ncbi.tbl, loci=publ_loci)

#function to


#function to vonvert prokka locus_tags prokka to ncbi locus_tags
ltag_prokka2ncbi <- function(prokka2ncbi.tbl, loci){
  message('translating ',length(loci), ' loci from prokka to ncbi')
  ncbi.selected <- prokka2ncbi.tbl %>%
    filter(grepl(paste(loci,collapse='|'),prokka_locus_tag))

  ncbi_loci <- ncbi.selected %>%
    select(ncbi_old_locus_tag) %>% pull() %>% as.character()
  message('revovered ',length(ncbi_loci),'/',length(loci), ' ncbi loci')

  return(ncbi_loci)
}

#example
#publ_loci <- c('D777_03329','D777_03330','D777_03331','D777_03332','D777_03333','D777_03334','D777_03335','D777_03336','D777_03337')
#prokka_loci <- ltag_ncbi2prokka(prokka2ncbi.tbl, loci=publ_loci)
#ncbi_loci <- ltag_prokka2ncbi(prokka2ncbi.tbl, loci=prokka_loci)

#==========================≠#

final.clusters.df <- read.delim('~/DATA/MarbGenomics/final_cluster_selection.txt')

#select clusters >1 lenght
compiled.clusters <- final.clusters.df %>%
  dplyr::group_by(c_function) %>%
  dplyr::tally() %>%
  dplyr::filter(n>1) %>%
  dplyr::select(c_function) %>%
  dplyr::filter(c_function != '') %>% #!!! fix earlier on
  dplyr::pull() %>%
  as.character()

#outstat_old is the non-KEGG based additions
#outstat_old <- outstat

#initialise
fromKEGG=FALSE
outstat <- c()

c.c <- compiled.clusters[22]
c.c <- "methionine transport"


all.out.compiled <- c()
for(c.c in compiled.clusters){
 prokka_loci <- final.clusters.df %>%
    dplyr::filter(c_function == c.c) %>%
    dplyr::select(locus_tag) %>%
    dplyr::pull() %>%
    as.character()

  associated_OGs <- annotation.tbl %>%
    filter(locus_tag %in% prokka_loci) %>%
    filter(!is.na(OG)) %>%
    select(OG) %>%
    unique() %>%
    pull()

  if(fromKEGG == TRUE){
    # Experimental, derive it from KO numbers, needs checking!
    associated_KOs <-final.clusters.df %>%
      dplyr::filter(c_function == c.c) %>%
      dplyr::select(kfm_domain) %>%
      dplyr::filter(!is.na(kfm_domain)) %>%
      dplyr::pull() %>% as.character()

    #
    OGfromKEGG <- annotation.tbl %>%
      dplyr::filter(grepl(paste(associated_KOs,collapse='|'),kfm_domain)) %>%
      dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
      dplyr::left_join(ko2ann, by=c('kfm_domain'='ko')) %>%
      dplyr::select(OG,kfm_domain,ko_annotation, product,hmm_acc) %>%
      unique() %>% select(OG) %>% unique() %>% dplyr::pull() %>% as.character()

    associated_OGs <- c(associated_OGs, OGfromKEGG) %>% unique()
  }




  selectedOG <- associated_OGs


  if(length(selectedOG) > 1){

    clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
    clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
    #clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=0,after=0)
    clusdetect.sorted <- sortregions(in.cluster=clusdetect.filled,majority = TRUE)
    #
    clusdetect.center <- clusdetect.sorted
    #??????????????????????
    #mapOGs <- detect_sco_in_region(clusdetect.sorted)
    #mapOG = mapOGs[1]
    #if(is.null(mapOG)){
    #  message('consider alternative')
    #  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)
    #}else{
    #  clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=mapOG)
    #}


    #clusdetect.center <- centerregions(in.cluster=clusdetect.sorted, OGs=selectedOG)

    #make function, but now here, to clean up low certainty, at least 2 genes long
    cleaned.g <- clusdetect.center %>% group_by(genome) %>% tally() %>% filter(n>1) %>% select(genome) %>% pull()
    clusdetect.center <- clusdetect.center %>% filter(genome %in% cleaned.g)

    clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = sel.genome)
    clusdetect.center$genome <- factor(clusdetect.center$genome,levels = sco.106.concat.tree$tip.label)
    focusOG=selectedOG

    #get order of OGs from the longest detected one
    orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()
    #


    tmp.outstat <- clusdetect.center %>%
      filter(!is.na(locus_tag)) %>%
      group_by(genome) %>% tally() %>%
      mutate(perc=n/max(n), func=c.c)
    outstat <- rbind(outstat, tmp.outstat)



    if(nrow(clusdetect.center %>% filter(!is.na(locus_tag)))>1){
      clusdetect.center <- clusdetect.center %>%
        dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
        dplyr::mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
        dplyr::mutate(OG = factor(OG, levels=og.oder))


      p = ggtree::ggtree(sco.106.concat.tree) + ggtree::geom_tiplab(align=TRUE)


      p <-    ggtree::facet_plot(p,
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



      ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/test_final_",gsub('\\/','',gsub(' ','_',gsub('\\+','',c.c))),'_fromkegg_',fromKEGG,".pdf"),
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



#-----------#
#all.out.compiled <- c()
#write.table(all.out.compiled)
write.table(all.out.compiled,
            file='~/DATA/MarbGenomics/Graphs/all.out.compuled.txt',
            sep='\t',
            quote = FALSE)

#-----------------------------------------------#



#=============================================================#
### visualise a selection ####

cl.names <- unique(all.out.compiled$c)

cl.names <- c("alginate",
              "glycogen",
              "purine",
              "vibrioferrin",
              "petrobactin",
              "phenol hydroxylase",
              "nap operon",
              "betacarotene_PR",
              "nitrate reductase",
              "nitric oxide reductase",
              "nitrous oxide reductase")

v.options <- rep(c('viridis','plasma','magma','cividis','inferno','mako'),10)
p = ggtree(sco.106.concat.tree) + geom_tiplab(align=TRUE,size=0)

for(i in 1:length(cl.names)){
  message(v.options[i])
  p <-    facet_plot(p,
                     panel=cl.names[i],
                     mapping = aes(x = start, xend = end, yend=y,color=OG),
                     data=all.out.compiled %>% filter(c==cl.names[i]),
                     geom=geom_segment,size=1.5)+
    theme(legend.position = "none")+
    viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
    ggnewscale::new_scale_color()
  #for some reason need to print in order for colors to work
  #print(p)
}

p <-p + xlim_tree(6) + coord_cartesian(clip='off')

ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.grgr",'_proximity_clusters2',".pdf"),
                plot=p,
                width = 12,
                height = 8,
                unit='in')

#=============================================================#



cl.names <- unique(all.out.compiled$c)


#

categories.l <- final.clusters.df %>% select(c_category) %>% unique() %>% pull() %>% as.character()

for(o in 1:length(categories.l)){
  cl.names <- final.clusters.df %>% filter(c_category == categories.l[o]) %>%
    select(c_function) %>%
    unique() %>% pull() %>%
    as.character()
  generic.name <- categories.l[o]

  p = ggtree(sco.106.concat.tree) + geom_tiplab(align=TRUE,size=0)

  lenghts.s <- c()
  for(i in 1:length(cl.names)){
    message(v.options[i])
    c.og <- all.out.compiled %>% filter(c==cl.names[i]) %>%
      filter(!is.na(OG)) %>%
      select(OG) %>% group_by(OG) %>%
      tally() %>% arrange(desc(n)) %>%
      select(OG) %>% slice(1) %>% pull() %>%
      as.character()
    s.og.clust <- all.out.compiled %>% filter(c==cl.names[i])


    #get order of OGs from the longest detected one
    orderfromg <- s.og.clust %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    og.oder <- s.og.clust %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull() %>% as.character()
    #dplyr::mutate(OG = factor(OG, levels=og.oder))


    s.t <- s.og.clust %>% arrange(start) %>% select(start) %>% slice(1) %>% pull()
    e.t <- s.og.clust %>% arrange(desc(end)) %>% select(start) %>% slice(1) %>% pull()
    l.l <- e.t - s.t

    lenghts.s <- c(lenghts.s, l.l)
    p <-
      facet_plot(p,
                 panel=cl.names[i],
                 mapping = aes(xmin = start, xmax = end, fill=OG),
                 data=s.og.clust %>%
                   dplyr::mutate(orientation=ifelse(strand=="+",1,-1),
                                 strand = ifelse(strand=='+','forward','reverse')) %>%
                   dplyr::mutate(OG = factor(OG, levels=og.oder)),
                 geom=geom_motif,
                 on=c.og,
                 label='gene',
                 align='centre',
                 arrowhead_width = unit(2, "mm"))+
      theme(legend.position = "none")+
      viridis::scale_fill_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
      ggnewscale::new_scale_fill()
  }


  #p


  lenghts.s <- lenghts.s/sum(lenghts.s)


  p <- facet_widths(p,width=c(0.25,lenghts.s))

  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.grgr_",generic.name,'_proximity_clusters3',".pdf"),
                  plot=p,
                  width = 18,
                  height = 20,
                  unit='in')
}



generic.name <- 'phosphoryl'
cl.names <- c("cytochrome c",
              "cytochrome o",
              "cytochrome bd",
              "Rnf",
              "atpase_ii",
              "Na+/H+ antiporter",'na_h_antiporter')



generic.name <- 'transp'
cl.names <- c("arabinose",
              "arginine_ornithine",
              "branched_chain_aa",
              "C4_dicarbox",
              "dicitrate",
              'fructose PTS',
              'general L-amino acid',
              'glucose mannose transport system',
              'glycerol transport',
              'glucose mannose transport system',
              'glycine betaine_proline',
              "hydroxyproline transporter",
              "iron transport",
              "iron_complex",
              "maltose uptake",
              "methionine transport",
              "microcin",
              "osmoprotectant transport system",
              "phenylacetate_operon_importer",
              "phosphonate uptake",
              "polar amino acid",
              "putative thiamine transport",
              "sulfate transport",
              "zink")


generic.name <- 'motility'
cl.names <- c("flagellum",
              "chemotaxis",
              "tight adherence",
              "rare_pilus")


generic.name <- 'final_selection'

cl.names <- c("alginate",
              "pel gene cluster",
              "glycogen",
              "general L-amino acid" ,"hydroxyproline transporter","putative thiamine transport",
              "putitaive maltose operon","maltose uptake","arabinose","fructose PTS",
              'glycerol transport',
              "purine",
              "putative methylmalonyl pathway",
              "phenylacetate","phenylacetate_operon_importer",
              'phosphonate utilisation',
              "osmoprotectant transport system",
              "NADH-quinone oxidoreductase",
              "Na+/H+ antiporter",
              "vibrioferrin",
              "petrobactin",
              "phenol hydroxylase",
              "betacarotene_PR",
              "tellurite reistance",
              "cytochrome o",
              "atpase_ii",
              "nap operon",
              "nitrate reductase",
              "nitric oxide reductase",
              "nitrous oxide reductase",
              "pecorrin_B12",
              "cobalamin_1",
              "cobalamin_2",
              "dmdA gene cluster",
              "dmsD_formate dehdrogenase",
              "sulfate transport",
              "tight adherence",
              "flagellum",
              "chemotaxis")




p = ggtree(sco.106.concat.tree) + geom_tiplab(align=TRUE,size=0)

  lenghts.s <- c()
  for(i in 1:length(cl.names)){
    c.og <- all.out.compiled %>% filter(c==cl.names[i]) %>%
      filter(!is.na(OG)) %>%
      select(OG) %>% group_by(OG) %>%
      tally() %>% arrange(desc(n)) %>%
      select(OG) %>% slice(1) %>% pull() %>%
      as.character()
    s.og.clust <- all.out.compiled %>% filter(c==cl.names[i])

    #get order of OGs from the longest detected one
    orderfromg <- s.og.clust %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    og.oder <- s.og.clust %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull() %>% as.character()
    #dplyr::mutate(OG = factor(OG, levels=og.oder))

    s.t <- s.og.clust %>% arrange(start) %>% select(start) %>% slice(1) %>% pull()
    e.t <- s.og.clust %>% arrange(desc(end)) %>% select(start) %>% slice(1) %>% pull()
    l.l <- e.t - s.t
    message(i,l.l)

    lenghts.s <- c(lenghts.s, l.l)
    p <-
      facet_plot(p,
                 panel=cl.names[i],
                 mapping = aes(xmin = start, xmax = end, fill=OG),
                 data=s.og.clust %>%
                   dplyr::mutate(orientation=ifelse(strand=="+",1,-1),
                                 strand = ifelse(strand=='+','forward','reverse')) %>%
                   dplyr::mutate(OG = factor(OG, levels=og.oder)),
                 geom=geom_motif,
                 on=c.og,
                 label='gene',
                 align='centre',
                 arrowhead_width = unit(2.5, "mm"),
                 arrow_body_height=grid::unit(2.5,'mm'))+
      theme(legend.position = "none")+
      viridis::scale_fill_viridis(discrete=TRUE,option='viridis',na.value='grey') +
      ggnewscale::new_scale_fill()
  }


  #p


  lenghts.s <- lenghts.s/sum(lenghts.s)


  p <- facet_widths(p,width=c(0.25,lenghts.s))

  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.",generic.name,'_proximity_clusters3',".pdf"),
                  plot=p,
                  width = 1000,
                  height = 5000,
                  unit='mm')




#=========================
# setting up a annotated phylogenetic tree



#assign backbone
p = ggtree(sco.106.concat.tree) +
    geom_tiplab(align=TRUE,size=0,hjust=-1,linetype='dashed',linesize=0.3)

#attach metadata (and set NA genome characteristics if genome quality is too low)
p$data <- p$data %>%
  left_join(genome.tbl, by=c('label'='genome')) %>%
  mutate(Completeness = ifelse(quality_class =='low',NA, Completeness))%>%
  mutate(Contamination = ifelse(quality_class =='low',NA, Contamination))%>%
  mutate(Genome_size = ifelse(quality_class =='low',NA, Genome_size))


#build layers
p <- p +
    ggnewscale::new_scale_fill()  +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=phylogroup),
                            color = NA, size=0,width=0.1)+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=phylogroup.colors)+
    ggnewscale::new_scale_fill()+
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=SS),
                          color = NA,
                          size=0,width=0.1,offset = 0.1)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=habitat.colors)+
  ggnewscale::new_scale_fill()#+
  # ggtreeExtra::geom_fruit(geom=geom_tile,
  #                         mapping=aes(fill=Genome_size),
  #                         color = NA,
  #                         size=0,width=0.1,offset = 0.1)+
  # scale_fill_distiller(palette='Greens',direction=1)+
  # ggnewscale::new_scale_fill() +
  # ggtreeExtra::geom_fruit(geom=geom_tile,
  #                         mapping=aes(fill=Completeness),
  #                         color = NA,
  #                         size=0,width=0.1,offset = 0.1)+
  # scale_fill_distiller(palette='Blues',direction=1)+
  # ggnewscale::new_scale_fill() +
  # ggtreeExtra::geom_fruit(geom=geom_tile,
  #                         mapping=aes(fill=Contamination),
  #                         color = NA,
  #                         size=0,width=0.1,offset = 0.1)+
  # scale_fill_distiller(palette='Purples',direction=1)


#-------------------------------------#
# ---- USING SEGMENTS INSTEAD OF GGGENES
#--------------------------------------#


  lenghts.s <- c()
  for(i in 1:length(cl.names)){
    #message(v.options[i])
    c.og <- all.out.compiled %>% filter(c==cl.names[i]) %>%
      filter(!is.na(OG)) %>%
      select(OG) %>% group_by(OG) %>%
      tally() %>% arrange(desc(n)) %>%
      select(OG) %>% slice(1) %>% pull() %>%
      as.character()
    s.og.clust <- all.out.compiled %>% filter(c==cl.names[i])

    #get order of OGs from the longest detected one
    orderfromg <- s.og.clust %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    og.oder <- s.og.clust %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull() %>% as.character()
    #dplyr::mutate(OG = factor(OG, levels=og.oder))

    s.t <- s.og.clust %>% arrange(start) %>% select(start) %>% slice(1) %>% pull()
    e.t <- s.og.clust %>% arrange(desc(end)) %>% select(start) %>% slice(1) %>% pull()
    l.l <- e.t - s.t
    message(i,l.l)

    lenghts.s <- c(lenghts.s, l.l)
    p <-
      facet_plot(p,
                 panel=cl.names[i],
                 mapping = aes(x = start, xend = end,yend=y, color=OG),
                 data=s.og.clust %>%
                   dplyr::mutate(orientation=ifelse(strand=="+",1,-1),
                                 strand = ifelse(strand=='+','forward','reverse')) %>%
                   dplyr::mutate(OG = factor(OG, levels=og.oder)),
                 geom=geom_segment,size=1.4)+
      ggplot2::theme_grey() +
      ggplot2::scale_y_continuous(minor_breaks = seq(0,106,1), expand=c(0.005,0.005))+
      ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.border=ggplot2::element_rect(colour='black',fill=NA),
            panel.spacing.x = unit(0.5,'mm'),
            legend.position = "none",
            strip.text.x = ggplot2::element_text(size=6,angle=90,hjust=0),
            strip.background = ggplot2::element_blank())+
      viridis::scale_color_viridis(discrete=TRUE,option='marko',na.value='grey') +
      ggnewscale::new_scale_color()


  }


  #p


  lenghts.s <- lenghts.s/sum(lenghts.s)


  p <- facet_widths(p,width=c(0.25,lenghts.s))

  #p <- p + theme_classic()

  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.segments3_annotated4.",generic.name,'_proximity_clusters3',".pdf"),
                  plot=p,
                  width = 360,
                  height = 200,
                  unit='mm')

  #p.rot <- p + coord_flip()
  #ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/prox.segments_flipped.",generic.name,'_proximity_clusters3',".pdf"),
  #                plot=p.rot,
  #                width = 180,
  #                height = 360,
  #                unit='mm')

#---------------------------------------------------------------#




outstat %>% ggplot(aes(genome,func))+
  geom_point(aes(color=perc)) +
  coord_flip()

outstat %>%
  tidyr::pivot_wider(id_cols=genome, names_from=func,values_from=perc)




matrix.stat <- outstat %>% tidyr::pivot_wider(id_cols=genome, names_from=func,values_from=perc) %>% data.frame()
rownames(matrix.stat)=matrix.stat$genome
matrix.stat <- matrix.stat[,2:ncol(matrix.stat)]
matrix.stat[is.na(matrix.stat)] <- 0


h.clust <- hclust(dist(t(matrix.stat)))
plot(h.clust)

pl.outstat<- outstat %>%
  mutate(ID=genome) %>%
  select(ID, everything())

pl.outstat$func <- factor(pl.outstat$func,levels=unique(outstat$func)[h.clust$order])

pl.outstat <- pl.outstat %>% left_join(final.clusters.df %>% select(c_function, c_category) %>% unique(),
                                       by=c('func'='c_function'))

#normalise by original cluster lenght (else it takes the lognest cluster, and that is not correct in many cases where some clusters are returned that are longer, like tandem repeats, we scale back to the original)
n_genes_original <- final.clusters.df %>%
  group_by(c_function) %>%
  tally() %>%
  mutate(original_n=n) %>%
  select(-n)

#
pl.outstat <- pl.outstat %>% left_join(n_genes_original %>% unique(),
                                       by=c('func'='c_function')) %>%
  mutate(dummy_n =ifelse(n>original_n,original_n,n)) %>%
  mutate(perc=dummy_n/original_n)



nrofcat <- pl.outstat$c_category %>% unique() %>% as.character() %>% length()

cat_colors <- viridis::viridis(nrofcat+10)[1:nrofcat]

ppp <- pl.outstat %>%
  mutate(genome=factor(pl.outstat$genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  dplyr::filter(perc>0.33) %>%
  ggplot(aes(genome,func))+geom_point(aes(color=c_category,alpha=perc)) +
  coord_flip() +
  facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(legend.position='none',
        panel.border=element_blank(),
        panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  scale_color_manual(values=cat_colors)


ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/test_final_fig_part1",".pdf"),
                plot=ppp,
                width = 10,
                height = 10,
                unit='in')


ppp <- pl.outstat %>%
  mutate(genome=factor(pl.outstat$genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  dplyr::filter(perc>0.33) %>%
  ggplot(aes(genome,func))+geom_point(aes(color=c_category,alpha=perc)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(legend.position='none',
        panel.border=element_blank(),
        panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  scale_color_manual(values=cat_colors)


ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/test_final_fig_part1_nonfacetted",".pdf"),
                plot=ppp,
                width = 10,
                height = 10,
                unit='in')






annotation.tbl %>% filter(grepl('toxin',product,ignore.case=TRUE)) %>%
  group_by(genome,product) %>% tally() %>%
  arrange(desc(n)) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  ggplot(aes(genome,product))+geom_point(aes(color=n)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(#legend.position='none',
        panel.border=element_blank(),
        panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  viridis::scale_color_viridis(option='magma')
  #scale_color_manual(values=cat_colors)



annotation.tbl %>% filter(grepl('alkane',product,ignore.case=TRUE)) %>%
  group_by(genome,product) %>% tally() %>%
  arrange(desc(n)) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  ggplot(aes(genome,product))+geom_point(aes(color=n)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(#legend.position='none',
    panel.border=element_blank(),
    panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  viridis::scale_color_viridis(option='magma')
#scale_color_manual(values=cat_colors)


alkane
cytochrome
oxidoreductase

annotation.tbl %>%
  dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
  filter(grepl('atp',hmm_name,ignore.case=TRUE)) %>%
  group_by(genome,hmm_name) %>% tally() %>%
  arrange(desc(n)) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  ggplot(aes(genome,hmm_name))+geom_point(aes(color=n)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(#legend.position='none',
    panel.border=element_blank(),
    panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  viridis::scale_color_viridis(option='magma')
#scale_color_manual(values=cat_colors)


annotation.tbl %>%
  dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
  filter(grepl('atp-synt',hmm_name,ignore.case=TRUE)) %>%
  group_by(genome,product) %>% tally() %>%
  arrange(desc(n)) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  ggplot(aes(genome,product))+geom_point(aes(color=n)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(#legend.position='none',
    panel.border=element_blank(),
    panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  viridis::scale_color_viridis(option='magma')
#scale_color_manual(values=cat_colors)




annotation.tbl %>%
  dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id')) %>%
  filter(grepl('p450',hmm_name,ignore.case=TRUE)) %>%
  group_by(genome) %>% tally() %>%
  arrange(desc(n)) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  ggplot(aes(genome,'GGDEF'))+geom_point(aes(color=n)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(#legend.position='none',
    panel.border=element_blank(),
    panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)#+
  viridis::scale_color_viridis(option='magma')
#scale_color_manual(values=cat_colors)






#+
#coord_flip()




pl.outstat %>%
  mutate(genome=factor(pl.outstat$genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  dplyr::filter(perc>0.33) %>%
  ggplot(aes(genome,func))+geom_point(aes(color=c_category,alpha=perc)) +
  facet_grid(c_category~.,scales='free_y',space='free_y')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(legend.position='none',
        panel.border=element_blank(),
        panel.spacing.y=unit(.1,'line'))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(values=cat_colors)





p = ggtree::ggtree(sco.106.concat.tree,layout = 'circular') + ggtree::geom_tiplab(align=TRUE,size=0)
p = ggtree::ggtree(sco.106.concat.tree) + ggtree::geom_tiplab(align=TRUE,size=0)

p + ggtreeExtra::geom_fruit(data=pl.outstat,
                            geom=geom_tile,
                            mapping=aes(y=ID,x=func, alpha=perc), pwidth=0.4)


p + ggtreeExtra::geom_fruit(data=pl.outstat,
                            geom=geom_point,
                            mapping=aes(y=ID,x=func,color = c_category, alpha=perc), pwidth=0.4)






p <-    ggtree::facet_plot(p,
                           panel=c.c,
                           mapping = aes(x = func),
                           data=outstat %>%
                             mutate(id=genome) %>%
                             select(id, everything()),
                           geom=geom_point,size=1.5)+
  theme(legend.position = "none")+
  viridis::scale_color_viridis(discrete=TRUE,option=v.options[i],na.value='grey') +
  ggnewscale::new_scale_color()

p








#========================================#

pan.name <- 'pfam'

#Load data from HGM test analysis
ANI.hmg.grp.pfam <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'hgm_allgroups.tsv'))
ANI.hmg.grp.pfam <- ANI.hmg.grp.pfam %>% as_tibble()
#ANI.hmg.grp.pfam$group <- factor(ANI.hmg.grp.pfam$group,levels=paste0('cl',1:length(unique(ANI.hmg.grp.pfam$group))))

#Load data from RF
load(paste0("~/DATA/MarbGenomics/",'multiRF',pan.name,sel.col,'.RData'))

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

HGM_RF_COMB_pfam <- ANI.hmg.grp.pfam %>%
  filter(level=='sc0.85') %>%
  mutate(comb = paste0(feature.id,'_',group)) %>%
  left_join( rf.all,by='comb')





#
selected.domains<- ANI.hmg.grp.pfam %>%
  group_by(group) %>%
  filter(level=='sc0.85') %>%
  arrange(fdr.p.value_enriched_rawcounts) %>%
  select(feature.id,fdr.p.value_enriched_rawcounts) %>%
  slice(1:20) %>%
  unique() %>%
  select(feature.id) %>%
  pull() %>%
  as.character()


pfam.plot.tbl <- pfam.pan[,selected.domains] %>% tibble::rownames_to_column(var='genome')%>% tidyr::pivot_longer(cols=-genome,names_to='feature.id',values_to='n')


sorted.domains <- pfam.pan[,selected.domains] %>% colSums() %>% sort(decreasing = TRUE) %>% names()
sorted.domains[1:20]


pfam.plot.tbl %>%
  filter(feature.id %in% sorted.domains[20:30]) %>%
  filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  mutate(n=ifelse(n==0,NA,n)) %>%
  ggplot(aes(genome,feature.id))+geom_point(aes(color=n)) +
  coord_flip() +
  #facet_grid(.~c_category,scales='free_x',space='free_x')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1))+
  theme(#legend.position='none',
    panel.border=element_blank(),
    panel.spacing.x=unit(0,'line'))+
  scale_y_discrete(limits=rev)+
  viridis::scale_color_viridis(option='magma',na.value=NA)




#------------------------------------------------------------------------------------------------#

# mess with ncbi to prokka and reverse, to be able to screen easily things like grouped_loci.txt

#### RUN CODE ####

#load precompiles table with all the intersting marinobacter functions
compiled.clusters.df <- read.delim('~/DATA/MarbGenomics/grouped_loci.txt')


#prokka2ncbi! from blastback.R
#prokka2ncbi.tbl <- tibble(genome=names(prokka2ncbi.annotated.list),prokka2ncbi.annotated.list) %>% unnest(prokka2ncbi.annotated.list)

prokka2ncbi.tbl %>% head()






publ_loci <- c('D777_03329','D777_03330','D777_03331','D777_03332','D777_03333','D777_03334','D777_03335','D777_03336','D777_03337')

#translate to prokka
prokka_loci <- ltag_ncbi2prokka(prokka2ncbi.tbl, loci=publ_loci)

#extract associated OGs
associated_OGs <- annotation.tbl %>%
  filter(locus_tag %in% prokka_loci) %>%
  filter(!is.na(OG)) %>%
  select(OG) %>% unique() %>% pull()


#look up in specified genome
annotation.tbl %>%
  dplyr::filter(OG %in% associated_OGs) %>%
  dplyr::filter(genome=='Marinobacter_algicola_DG893') %>%
  dplyr::select(locus_tag,Name, product,OG) %>%



clst <- multigenedetect(OG=associated_OGs, genomes='Marinobacter_algicola_DG893')

largest.cluster <- clst %>%
   dplyr::group_by(prox.cluster) %>%
   dplyr::tally() %>%
   arrange(desc(n)) %>%
   dplyr::slice(1) %>%
   select(prox.cluster) %>% pull()

clst <- clst %>% dplyr::filter(prox.cluster ==largest.cluster)

ref.prokka.loci <- clst %>% select(locus_tag) %>% pull()



#translate back to ncbi locus tags
ncbi_loci <- ltag_prokka2ncbi(prokka2ncbi.tbl, loci=ref.prokka.loci)


paste(ncbi_loci, collapse=' or ')



#






ref.genome <- ncbi.selected %>% select(genome) %>% unique() %>% pull() %>% as.character()


associated_OGs














#select clusters >1 lenght
compiled.clusters <- compiled.clusters.df %>%
  group_by(Group) %>%
  tally() %>%
  filter(n>1) %>%
  select(Group) %>%
  pull() %>%
  as.character()


c.c <- "urea_amidolyase"

publ_loci <- compiled.clusters.df %>%
  filter(Group %in% c.c) %>%
  select(Gene_names) %>%
  pull() %>%
  as.character()








c.c <- compiled.clusters[11]
for(c.c in compiled.clusters){
  publ_loci <- compiled.clusters.df %>% filter(Group %in% c.c) %>% select(Gene_names) %>% pull() %>% as.character()

  ncbi.selected <- prokka2ncbi.tbl %>%
    filter(grepl(paste(publ_loci,collapse='|'),ncbi_old_locus_tag))

  prokka_loci <- ncbi.selected %>%
    select(prokka_locus_tag) %>% pull() %>% as.character()

  ref.genome <- ncbi.selected %>% select(genome) %>% unique() %>% pull() %>% as.character()

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



    ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/proximity_clusters",c.c,".pdf"),
                    plot=p,
                    width = 12,
                    height = 8,
                    unit='in')

  }

}



















lactate utilisation




annotation.tbl %>% filter(kfm_domain == 'K00925') %>% select(Name, COG, gene,product, OG) %>% unique()

Fermentation	Pyruvate oxidation		porA	pyruvate ferredoxin oxidoreductase alpha subunit	K00169.hmm	K00169
Fermentation	Alcohol utilization		adh	alcohol dehydrogenase	K00001.hmm	K00001
Fermentation	Lactate utilization		ldh	L-lactate dehydrogenase	K00016.hmm	K00016
Fermentation	Acetogenesis		acdA	acetyl coenzyme A synthetase (ADP forming), alpha domain	TIGR02717.hmm	K01905
Fermentation	Acetogenesis		ack	acetate kinase	TIGR00016.hmm	K00925
Fermentation	Acetogenesis		pta	phosphate acetyltransferase	TIGR00651.hmm	K00625
Fermentation	Acetate to acetyl-CoA		acs	acetyl-CoA synthetase	TIGR02188.hmm	K01895
Fermentation	Pyruvate <=> acetyl-CoA + formate 		pflD	formate C-acetyltransferase	K00656.hmm	K00656


annotation.tbl %>% filter(kfm_domain == 'K00169') %>% select(Name, COG, gene,product, OG) %>% unique()
annotation.tbl %>% filter(kfm_domain == 'K00001') %>% select(Name, COG, gene,product, OG) %>% unique()
annotation.tbl %>% filter(kfm_domain == 'K00016') %>% select(Name, COG, gene,product, OG) %>% unique()
annotation.tbl %>% filter(kfm_domain == 'K01905') %>% select(Name, COG, gene,product, OG) %>% unique()
annotation.tbl %>% filter(kfm_domain == 'K00925') %>% select(Name, COG, gene,product, OG) %>% unique()
annotation.tbl %>% filter(kfm_domain == 'K00625') %>% select(Name, COG, gene,product, OG) %>% unique()



annotation.tbl %>% filter(OG %in% c('OG0002400','OG0002501'))





ko2mod <-







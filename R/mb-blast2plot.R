
i='NKSG1_ATPase.results.out'

balst.res <- read.delim(paste0('~/DATA/MarbGenomics/',i),header=FALSE)


colnames(balst.res) =c('qseqid','sseqid','slen', 'qstart', 'qend', 'length', 'mismatch', 'gapopen', 'gaps', 'evalue', 'bitscore',  'qcovs', 'qcovhsp', 'sseq')
balst.res <- balst.res %>% as_tibble()

selected.genome <- balst.res %>%
  left_join(annotation.tbl, by=c('sseqid'='locus_tag')) %>%
  filter(evalue==0)



Hahella_chejuensis_KCTC_2396
annotation.tbl %>% filter(genome=='Hahella_chejuensis_KCTC_2396') %>%
  filter(grepl('atp',gene)) %>% select(locus_tag, gene, product, COG, OG) %>% data.frame()

annotation.tbl %>% filter(genome=='Marinobacter_santoriniensis_NKSG1') %>%
  filter(grepl('atp',gene)) %>% select(locus_tag, gene, product, COG, OG) %>% data.frame()



selected.genome <- 'Hahella_chejuensis_KCTC_2396'
i='KCTC_2396_ATPase.results.out'
sel.loci <- c('BCBNPALP_00845',
              'BCBNPALP_00846',
              'BCBNPALP_00847',
              'BCBNPALP_00848',
              'BCBNPALP_00849',
              'BCBNPALP_00850',
              'BCBNPALP_00851',
              'BCBNPALP_00852')

selected.genome <- 'Marinobacter_santoriniensis_NKSG1'
i='NKSG1_ATPase.results.out'

sel.loci <- c('NDLGLMCA_00106',
'NDLGLMCA_00107',
'NDLGLMCA_00108',
'NDLGLMCA_00109',
'NDLGLMCA_00110',
'NDLGLMCA_00111',
'NDLGLMCA_00112',
'NDLGLMCA_00113')

i='NKSG1_ATPase2.results.out'

sel.loci <- c('NDLGLMCA_00851',
              'NDLGLMCA_00852',
              'NDLGLMCA_00853',
              'NDLGLMCA_00854',
              'NDLGLMCA_00855',
              'NDLGLMCA_00856',
              'NDLGLMCA_00857',
              'NDLGLMCA_00858',
              'NDLGLMCA_00859')

i='NKSG1_ATPase3.results.out'

sel.loci <- c('NDLGLMCA_02910',
              'NDLGLMCA_02911',
              'NDLGLMCA_02912',
              'NDLGLMCA_02913',
              'NDLGLMCA_02914',
              'NDLGLMCA_02915',
              'NDLGLMCA_02916',
              'NDLGLMCA_02917')

 NDLGLMCA_02917
#=
# blast 2 plot


balst.res <- read.delim('~/DATA/MarbGenomics/antismash_algicola_cluster.results.out',header=FALSE)

antismash_sid_SM19_cluster.results.out

i = 'antismash_sid_vt8.results.out'
selected.genome = 'Marinobacter_hydrocarbonoclasticus_VT8'

i = 'phenol_Hiroaki.results.out'
i = 'vib_alg.results.out'
selected.genome = 'Marinobacter_algicola_DG893'


i = 'antismash_sid_SM19_cluster.results.out'
selected.genome = 'Marinobacter_lipolyticus_SM19'



i = 'antismash_nrps_d15_cluster.results.out'
selected.genome = 'Marinobacter_nanhaiticus_D15_8W'

i = 'D15W8_marinobactin.results.out'
selected.genome = 'Marinobacter_nanhaiticus_D15_8W'



tip.order = tip_order(lsTrees.106.rooted[[1]])


sel.loci <- c('HJKODCBM_00534','HJKODCBM_00535','HJKODCBM_00536','HJKODCBM_00537','HJKODCBM_00538','HJKODCBM_00539')
i = 'FDB33_pvdQ'
selected.genome = 'Marinobacter_sp_FDB33'

sel.loci <- c('EHILCICM_03866','EHILCICM_03867','EHILCICM_03868','EHILCICM_03869')
i = 'D15_W8_pvdQ'
selected.genome = 'Marinobacter_nanhaiticus_D15_8W'


i = 'ELB17_xanthoferrin.results.out'
selected.genome = 'Marinobacter_sp_ELB17'


i = 'glycolate_operon.results.out'
capsular_operon.results.out
i = 'flp_operon.results.out'

i='NKSG1_ATPase.results.out'

Glycolate_HP15_2.results.out


i = 'Glycolate_HP15_2.results.out'
selected.genome = 'Marinobacter_adhaerens_HP15'

i = 'maltose_redo.results.out'
i = 'perchlorate.results.out'
i='nitrite_red_algicola.results.out'

i='pel_operon.results.out'


balst.res <- read.delim(paste0('~/DATA/MarbGenomics/',i),header=FALSE)


colnames(balst.res) =c('qseqid','sseqid','slen', 'qstart', 'qend', 'length', 'mismatch', 'gapopen', 'gaps', 'evalue', 'bitscore',  'qcovs', 'qcovhsp', 'sseq')
balst.res <- balst.res %>% as_tibble()

selected.genome <- balst.res %>%
  left_join(annotation.tbl, by=c('sseqid'='locus_tag')) %>%
  filter(evalue==0) %>%
  select(sseqid, evalue, bitscore, genome, gene, product, COG, OG) %>% group_by(genome) %>%
  select(sseqid) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull()
message('selected genome:', selected.genome)

sel.loci <- balst.res %>%
  left_join(annotation.tbl, by=c('sseqid'='locus_tag')) %>%
  #filter(genome=='Marinobacter_algicola_DG893') %>%
  filter(genome==selected.genome) %>%
  filter(evalue==0) %>%
  select(sseqid, evalue, bitscore, genome, gene, product, COG, OG) %>%
  select(sseqid) %>%
  pull() %>%
  unique() %>%
  as.character()
message('# loci:', length(sel.loci))


#sel.loci=c("MEIJJOJH_02097","MEIJJOJH_02098","MEIJJOJH_02099","MEIJJOJH_02100","MEIJJOJH_02101")
#sel.loci <- sel.loci[2:(length(sel.loci)-1)]

#sel.loci = c(sel.loci, 'NHGOFBLF_00090','NHGOFBLF_00089','NHGOFBLF_00088','NHGOFBLF_00103','NHGOFBLF_00104','NHGOFBLF_00105')

#the loop

prox.clust <- proximity_clusters(annotation.tbl, sel.loci,distance=3000 )

  #extract clusters formed by >3 clusters
largest.clusts_all <- prox.clust %>%
    #group_by(prox.cluster) %>%
    dplyr::count(prox.cluster) %>%
    filter(n>1) %>%
    arrange(desc(n)) %>%
    select(prox.cluster) %>% pull()


if(length(largest.clusts_all)>0){
  selectedOG <- prox.clust %>% filter(prox.cluster == largest.clusts_all[1]) %>% select(OG) %>% pull() %>% as.character()

  clusdetect.og <- multigenedetect(OG=selectedOG, genomes=sel.genome)
  clusdetect.filled <- fillgaps(in.cluster=clusdetect.og)
  clusdetect.expanded <- expandregion(in.cluster = clusdetect.filled, before=1500,after=1500)
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
    clusdetect.center <- genome_dummy_genes(in.cluster = clusdetect.center, genomes = tip.order)

    clusdetect.center$genome <- factor(clusdetect.center$genome,levels = tip.order)

    focusOG=selectedOG

    #orderfromg <- clusdetect.center %>% group_by(genome) %>% tally() %>% arrange(desc(n)) %>% slice(1) %>% select(genome) %>% pull %>% as.character()
    #
    orderfromg <- selected.genome
    og.oder <- clusdetect.center %>% filter(genome==orderfromg) %>% select(OG)  %>% unique() %>% pull()



    out.p <- clusdetect.center %>%
      mutate(direction=ifelse(strand=="+",1,-1)) %>%
      mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
      mutate(OG = factor(OG, levels=og.oder)) %>%
      mutate(genome=factor(genome,levels=(tip.order))) %>%
      ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=Name,fill=OG)) + #color=SELCOG
      geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label(align = "left")+
      scale_fill_viridis(discrete = TRUE)+
      #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +
      theme(legend.position = "none")

    out.p
}


    ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/operons/",i,'_',"_testing4",'clusterdistribution',".pdf"),
                    plot=out.p,
                    width = 7,
                    height = 14,
                    unit='in')



clusdetect.center %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(locus_tag, OG, product,tcdb_domain) %>% left_join(df.pfam.all , by = c('locus_tag'='seq_id')) %>% select(locus_tag, OG, product,tcdb_domain,hmm_acc, hmm_name) %>% data.frame()
clusdetect.center %>% filter(genome=='Marinobacter_sp_ELB17') %>% select(locus_tag, OG, product,tcdb_domain) %>% left_join(df.pfam.all , by = c('locus_tag'='seq_id')) %>% select(locus_tag, OG, product,tcdb_domain,hmm_acc, hmm_name) %>% data.frame()


sel.loci<-


i = 'manual_vibrioferrin'

vibrioferrin_annotated <- rbind(
  c('locus_tag'="NHGOFBLF_00090",'gene_name'='pvuE','class'='Siderophore mediated iron acquisition ABC transporters'),
  c('locus_tag'="NHGOFBLF_00091",'gene_name'='pvuD','class'='Siderophore mediated iron acquisition ABC transporters'),
  c('locus_tag'="NHGOFBLF_00092",'gene_name'='pvuC','class'='Siderophore mediated iron acquisition ABC transporters'),
  c('locus_tag'="NHGOFBLF_00093",'gene_name'='pvuB','class'='Siderophore mediated iron acquisition ABC transporters'),
  c('locus_tag'="NHGOFBLF_00094",'gene_name'='pvuA','class'='TonB-dependend siderophore outer membrane receptors'),
  c('locus_tag'="NHGOFBLF_00095",'gene_name'='pvsE','class'='Vibrioferrin biosynthesis'),
  c('locus_tag'="NHGOFBLF_00096",'gene_name'='pvsD','class'='Vibrioferrin biosynthesis'),
  c('locus_tag'="NHGOFBLF_00097",'gene_name'='pvsC','class'='Vibrioferrin biosynthesis'),
  c('locus_tag'="NHGOFBLF_00098",'gene_name'='pvsB','class'='Vibrioferrin biosynthesis'),
  c('locus_tag'="NHGOFBLF_00099",'gene_name'='pvsA','class'='Vibrioferrin biosynthesis'),
  c('locus_tag'="NHGOFBLF_00100",'gene_name'='pvsX','class'='Vibrioferrin biosynthesis')
                                ) %>% data.frame(stringAsFactors=FALSE)

vibrioferrin_annotated <- vibrioferrin_annotated %>%
  left_join(annotation.tbl, by='locus_tag') %>%
  select(locus_tag, gene_name, class, OG)

sel.loci <- vibrioferrin_annotated$locus_tag %>% as.character()

#ok, ELB17 has it but it has been split, by contig
sel.loci <- vibrioferrin_annotated %>% filter(grepl('Vibrioferrin',class)) %>% select(locus_tag) %>% pull()


#=======

# or from genbank accession numbers. obsolete
sel.loci <- vibrioferrin_annotated$locus_tag %>% as.character()

i = 'DG893_vibrio_amin2012.results.out'
selected.genome = 'Marinobacter_algicola_DG893'





  out.p <- clusdetect.center %>%
    left_join(vibrioferrin_annotated, by='OG') %>%
    mutate(direction=ifelse(strand=="+",1,-1)) %>%
    mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
    mutate(OG = factor(OG, levels=og.oder)) %>%
    mutate(genome=factor(genome,levels=tip.order)) %>%
    ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=gene_name,fill=class)) + #color=SELCOG
    geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label()+#align = "left")+
    scale_fill_manual(values=c('skyblue','skyblue3','salmon'))+
    ggplot2::theme_grey() + ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "darkgrey", size = 0.5),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
      axis.ticks.x = ggplot2::element_line(colour = "grey20", size = 0.5),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )+ theme(legend.position = "bottom")#+scale_fill_viridis(discrete = TRUE)#+
    #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +


  out.p


  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/operons/",i,'_',"_testing6",'clusterdistribution',".pdf"),
                  plot=out.p,
                  width = 8,
                  height = 17,
                  unit='in')



  i = 'VT8_petro_amin2012.results.out'
  selected.genome = 'Marinobacter_hydrocarbonoclasticus_VT8'

  vibrioferrin_annotated <- rbind(
    c('locus_tag'="IFOEGCPI_04145",'gene_name'='fatB','class'='TonB-dependend siderophore outer membrane receptors'),
    c('locus_tag'="IFOEGCPI_04144",'gene_name'='fpuA','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="IFOEGCPI_04143",'gene_name'='fpuB','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="IFOEGCPI_04142",'gene_name'='fbuC','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="IFOEGCPI_04141",'gene_name'='fbuD','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="IFOEGCPI_04140",'gene_name'='asbA','class'='Petrobactin biosynthesis'),
    c('locus_tag'="IFOEGCPI_04139",'gene_name'='asbB','class'='Petrobactin biosynthesis'),
    c('locus_tag'="IFOEGCPI_04138",'gene_name'='asbC','class'='Petrobactin biosynthesis'),
    c('locus_tag'="IFOEGCPI_04137",'gene_name'='asbD','class'='Petrobactin biosynthesis'),
    c('locus_tag'="IFOEGCPI_04136",'gene_name'='asbE','class'='Petrobactin biosynthesis'),
    c('locus_tag'="IFOEGCPI_04135",'gene_name'='asbF','class'='Petrobactin biosynthesis')
  ) %>% data.frame(stringAsFactors=FALSE)

  vibrioferrin_annotated <- vibrioferrin_annotated %>%
    left_join(annotation.tbl, by='locus_tag') %>%
    select(locus_tag, gene_name, class, OG)

  out.p <- clusdetect.center %>%
    left_join(vibrioferrin_annotated, by='OG') %>%
    mutate(direction=ifelse(strand=="+",1,-1)) %>%
    mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
    mutate(OG = factor(OG, levels=og.oder)) %>%
    mutate(genome=factor(genome,levels=tip.order)) %>%
    ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=gene_name,fill=class)) + #color=SELCOG
    geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label()+#align = "left")+
    scale_fill_manual(values=c('salmon','skyblue','skyblue3'))+
    ggplot2::theme_grey() + ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "darkgrey", size = 0.5),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
      axis.ticks.x = ggplot2::element_line(colour = "grey20", size = 0.5),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )+ theme(legend.position = "bottom")#+scale_fill_viridis(discrete = TRUE)#+
  #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +


  out.p


  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/operons/",i,'_',"_testing6",'clusterdistribution',".pdf"),
                  plot=out.p,
                  width = 8,
                  height = 17,
                  unit='in')




  i = 'DG893_hydrox_amin2012.results.out'
  selected.genome = 'Marinobacter_algicola_DG893'





  vibrioferrin_annotated <- rbind(
    c('locus_tag'="NHGOFBLF_00375",'gene_name'='fhuB','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="NHGOFBLF_00376",'gene_name'='fhuC','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="NHGOFBLF_00377",'gene_name'='fhuD','class'='Siderophore mediated iron acquisition ABC transporters'),
    c('locus_tag'="NHGOFBLF_00378",'gene_name'='fhux','class'='Siderophore mediated iron acquisition ABC transporters')
  ) %>% data.frame(stringAsFactors=FALSE)

  vibrioferrin_annotated <- vibrioferrin_annotated %>%
    left_join(annotation.tbl, by='locus_tag') %>%
    select(locus_tag, gene_name, class, OG)

  out.p <- clusdetect.center %>%
    left_join(vibrioferrin_annotated, by='OG') %>%
    mutate(direction=ifelse(strand=="+",1,-1)) %>%
    mutate(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal','z'))) %>%
    mutate(OG = factor(OG, levels=og.oder)) %>%
    mutate(genome=factor(genome,levels=tip.order)) %>%
    ggplot(aes(xmin = start, xmax = end, y = genome, forward=direction,label=gene_name,fill=gene_name)) + #color=SELCOG
    geom_gene_arrow(arrowhead_height = unit(2.7, "mm"), arrowhead_width = unit(2, "mm")) + geom_gene_label()+#align = "left")+
   # scale_fill_manual(values=c('salmon','skyblue','skyblue3'))+
    ggplot2::theme_grey() + ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "darkgrey", size = 0.5),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
      axis.ticks.x = ggplot2::element_line(colour = "grey20", size = 0.5),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )+ theme(legend.position = "bottom")+
  scale_fill_viridis(discrete = TRUE)#+
  #scale_color_manual(values=c('blue','red','white'))+fdb_style(aspect.ratio=1) +


  out.p


  ggplot2::ggsave(filename= paste0("~/DATA/MarbGenomics/Graphs/operons/",i,'_',"_testing6",'clusterdistribution',".pdf"),
                  plot=out.p,
                  width = 8,
                  height = 17,
                  unit='in')


#==== pvdQ


  annotation.tbl %>% filter(grepl('pvdQ',Name))



  testAnn <- genome_dummy_genes(in.cluster = annotation.tbl %>% filter(grepl('pvdQ',Name)), genomes = tip.order)

  testAnn%>% group_by(genome, gene) %>% tally %>% mutate(n=ifelse(grepl('pvdQ',gene),n,0)) %>% mutate(genome=factor(genome,levels=(tip.order))) %>% ggplot(aes(x=gene,y=genome))+geom_point(aes(size=n))




  #============


dginc <- loadTree('~/DATA/MarbGenomics/DG_inclusive.fa.treefile',clean_genome_names = TRUE)

  dginc <- midpoint.root(dginc)
ggtree(dginc)+geom_tiplab(size=3) + geom_text2(aes(label=label, subset = isTip))

ggtree(dginc) +geom_tiplab(size=3) + geom_text2(aes(subset=!isTip, label=label),hjust=1,vjust=0.5,size=3)

andginc <- data.frame(rbind(c('AY258106.1','Marinobacter sp. DG870'),
c('AY258107.1', 'Marinobacter sp. DG879'),
c('AY258110.2', 'Marinobacter algicola DG893'),
c('AY258112.1', 'Marinobacter sp. DG979'),
c('AY258113.2', 'Marinobacter sp. DG980'),
c('AY258116.2', 'Marinobacter algicola strain DG1136'),
c('EU052743.1', 'Marinobacter sp. MH125a')))
colnames(andginc)=c('accession','strain')

dginc$tip.label[dginc$tip.label == 'AY258106_1'] <- 'Marinobacter sp. DG870'
dginc$tip.label[dginc$tip.label == 'AY258107_1'] <- 'Marinobacter sp. DG879'
dginc$tip.label[dginc$tip.label == 'AY258110_2'] <- 'Marinobacter algicola DG893'
dginc$tip.label[dginc$tip.label == 'AY258112_1'] <- 'Marinobacter sp. DG979'
dginc$tip.label[dginc$tip.label == 'AY258113_2'] <- 'Marinobacter sp. DG980'
dginc$tip.label[dginc$tip.label == 'AY258116_2'] <- 'Marinobacter algicola strain DG1136'
dginc$tip.label[dginc$tip.label == 'EU052743_1'] <- 'Marinobacter sp. MH125a'

ggtree(dginc) +geom_tiplab(size=3,aes(color=ifelse(label %in% andginc$strain,'a','b'))) +
  geom_text2(aes(subset=!isTip, label=label),hjust=1,vjust=0.5,size=3) +
  theme(legend.position = 'none')+
  scale_color_manual(values=c('firebrick','black'))+
  xlim(c(0,1))

p <- ggtree(dginc)
p$data <- p$data %>% left_join(genome.tbl, by=c('label'='genome'))
p +geom_tiplab(size=3) +
  #geom_text2(aes(subset=!isTip, label=label),hjust=1,vjust=0.5,size=3) +
  theme(legend.position = 'none')+
  geom_point(aes(color=phylogroup))+
  scale_color_manual(values=phylogroup.colors)+
  xlim(c(0,1))



dginc <- loadTree('~/DATA/MarbGenomics/dgall_inc.fa.treefile',clean_genome_names = TRUE)
dginc <- midpoint.root(dginc)

p <- ggtree(dginc)
p$data <- p$data %>% left_join(genome.tbl, by=c('label'='genome'))
p +geom_tiplab(size=3) +
  #geom_text2(aes(subset=!isTip, label=label),hjust=1,vjust=0.5,size=3) +
  theme(legend.position = 'none')+
  scale_color_manual(values=c(phylogroup.colors,'black'))+
  geom_point(aes(color=phylogroup))+
  #scale_fill_manual(values=c(phylogroup.colors,'black'))+
  xlim(c(0,0.3))

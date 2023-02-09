


# Absense of ANR transcriptional regulator in many strains of P7

annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('Transcriptional activator protein Anr',product,ignore.case=TRUE)) %>%
  filter(phylogroup=='P7') %>%
  group_by(genome) %>%
  tally() %>%
  arrange(n) %>%
  data.frame()




#ldhA and lldD absent in P7.1, but present in P7.2 (W6S, HL_58, halotolerans and confluentis)
annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('D-lactate dehydrogenase',product,ignore.case=TRUE)) %>%
  select(genome, product, Name, phylogroup) %>%
  #filter(phylogroup=='P7') %>%
  group_by(genome) %>%
  tally() %>%
  arrange(n) %>%
  data.frame()

annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('adha',Name,ignore.case=TRUE)) %>%
  select(genome, product, Name, phylogroup) %>%
  #filter(phylogroup=='P7') %>%
  group_by(genome) %>%
  tally() %>%
  arrange(n) %>%
  data.frame()



annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('acetaldehyde',product,ignore.case=TRUE)) %>%
  select(genome, product, Name, phylogroup,kfm_domain) %>%
  filter(phylogroup=='P7') %>%
  group_by(genome) %>%
  tally() %>%
  arrange(n) %>%
  data.frame()



annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('oxygenase',product,ignore.case=TRUE)) %>%
  select(genome, product, Name, phylogroup,kfm_domain) %>%
  filter(phylogroup=='P7') %>%
  group_by(genome) %>%
  tally() %>%
  arrange(n) %>%
  data.frame()




annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('anaer',product,ignore.case=TRUE)) %>%
  select(genome, product, Name, phylogroup,kfm_domain) %>%
  #filter(phylogroup=='P7') %>%
  group_by(product) %>%
  tally()# %>%
  #arrange(n) %>%
  #data.frame()











annotation.tbl %>%
  left_join(genome.tbl %>%
              select(genome, phylogroup),by='genome') %>%
  filter(grepl('lactate dehydrogenase',product,ignore.case=TRUE)) %>%
  select(genome, product, Name, phylogroup) %>%
  filter(phylogroup=='P7') %>%
  group_by(genome) %>%
  tally() %>%
  arrange(n) %>%
  data.frame()




#' function that writes fasta file from a two column dataframe with names and sequences
#'
#' @param data a two column dataframe with names and sequences as columns
#' @param filename name of the output file
#'
#' @return no return, it outputs a file in the specified location instrad
#' @export
#'
#' @examples
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,1], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,2]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


#' wrapper function to plot values from dataframe using ggplot2
#'
#' @param df dataframe
#' @param value column of interest
#' @param xlab label for the x axis
#' @param binwidth binwidth of the geom_histogram
#' @param div can be set to divide the values by given number (default = 1) for example to scale to whole numbers for Megabases, set div to 1000000
#'
#' @return gglot2 object
#' @export
#'
#' @examples
#'
plot.distibution<- function(df, value, xlab, binwidth, div=1){
  df[,value]=df[,value]/div
  p.pg <- df %>%
    ggplot2::ggplot(aes_string(x = value)) +
    ggplot2::geom_histogram( aes(y = ..density..), binwidth = binwidth, color = 'white', fill = "grey30", alpha = .3) +
    ggplot2::geom_density(alpha = .3, fill = "antiquewhite3",color ='antiquewhite3') +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::xlab(xlab) +
    ggplot2::xlim(min(df[,value]),max(df[,value])) +
    ggplot2::geom_vline(aes(xintercept=base::mean(df[,value]%>% dplyr::pull())),color = "grey50", linetype = "dashed")+
    ggplot2::geom_vline(aes(xintercept=stats::median(df[,value]%>% dplyr::pull())),color = "grey50")+
    ggplot2::theme_classic()

  xa <- ""
  p.pg.bx <- df %>%
    ggplot2::ggplot(aes_string(x=1, y = value)) +
    ggplot2::geom_boxplot(fill = "antiquewhite3",color="antiquewhite3", alpha=0.5) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylim(min(df[,value]%>% dplyr::pull()),max(df[,value]%>% dplyr::pull())) +
    ggplot2::theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  p.pg_b <- cowplot::plot_grid(p.pg.bx, p.pg,
                               ncol = 1, rel_heights = c(0.3, 2),
                               align = 'v', axis = 'lr')

  return(p.pg_b)
}



#' Title
#' function to check if genne family is shared by all selected genomes (in larger pangenome)
#' @param pan pan/abundance table
#' @param selected_genomes selected genomes
#'
#' @return named vector boolean, TRUE if shared by all. names are feature names
#' @export
#'
#' @examples
#' selected.genomes <- c("Marinobacter_bohaiensis_T17","Marinobacter_fonticola_CS412","Marinobacter_nanhaiticus_D15_8W")
#' is.shared(pan, selected.genomes)
is.shared <- function(pan, selected_genomes){
  pan.presabs <- pan[selected_genomes,]
  pan.presabs[pan.presabs>1]<-1
  out <- colSums(pan.presabs) == length(selected_genomes)
  return(out)
}



#========================================#
## partitioning of annotations ####
# ---
# ---
# ---

# collect pfam domains

# knowns if annotated with at least one Pfam domain
# unknowns if no pfam domain detected or if entirely annotated with Pfam domains of unknown function (DUF)

# known with Pfam
# known without Pfam (excluding DUF)
# unknown with DUF
# unknown without DUF


# summarise PFAM annotation quality per locus
pfam.ann.q <- df.pfam.all %>%
  dplyr::select(genome,seq_id, hmm_acc, hmm_name) %>%
  dplyr::mutate(unknown_pfam_function = ifelse(grepl('DUF',hmm_name),TRUE,FALSE)) %>%
  dplyr::group_by(genome,seq_id) %>%
  dplyr::summarise(n_known_pfam = sum(!unknown_pfam_function),
                   n_DUF_pfam = sum(unknown_pfam_function)) %>%
  dplyr::mutate(pfam_annotation_q = ifelse(n_known_pfam>0, 'known', 'unknown')) %>%
  dplyr::ungroup() %>%
  dplyr::select(-genome)

# construct tibble from prokka2ncbi dataset
prokka2ncbi.annotated.tbl <- tibble::tibble(genome=names(prokka2ncbi.annotated.list),
                                            df = prokka2ncbi.annotated.list) %>%
                              tidyr::unnest(df)




# annotation quality per locus_tag
annotation.quality.tbl <-
  annotation.tbl %>% filter(type == 'CDS') %>%
  filter(genome %in% lsTrees.94[[1]]$tip.label) %>%
  mutate(kfm_domain = ifelse(kfm_domain=='NA',NA,kfm_domain)) %>%
  dplyr::left_join(ko2ann, by=c('kfm_domain'='ko')) %>%
  dplyr::left_join(cog_func_single, by='COG') %>%
  dplyr::left_join(df.pfam_per_locus %>%
                     mutate(pfam_acc=hmm_acc,
                            pfam_name=hmm_name,
                            pfam_type=type) %>%
                     select(seq_id, pfam_acc, pfam_name, pfam_type)
                   , by=c('locus_tag'='seq_id')) %>%
  dplyr::left_join(pfam.ann.q, by=c('locus_tag'='seq_id')) %>%
  dplyr::mutate(pfam_annotation_q =ifelse(is.na(pfam_annotation_q),'none',pfam_annotation_q)) %>%
  dplyr::select(-score,-note,-rpt_family,-rpt_type,-rpt_unit_seq) %>%
  dplyr::left_join(prokka2ncbi.annotated.tbl %>%
                     dplyr::select(-genome),
                   by=c('locus_tag'='prokka_locus_tag')) %>%
  dplyr::select(genome, locus_tag,ncbi_old_locus_tag,OG, Name, COG, gene, product, kfm_domain,ko_annotation,pfam_acc, pfam_name,n_known_pfam,n_DUF_pfam,pfam_annotation_q) %>%
  dplyr::mutate(knownKO = ifelse(!is.na(kfm_domain) & !grepl('uncharacterized protein',ko_annotation), TRUE, FALSE),
                knownCOG = ifelse(!is.na(COG), TRUE, FALSE),
                KOorCOG = ifelse(knownKO|knownCOG, TRUE, FALSE),
                knownPFAM = ifelse(pfam_annotation_q == 'known', TRUE, FALSE)
                ) %>%
  dplyr::mutate(annotation_quality = ifelse(KOorCOG==TRUE & knownPFAM == TRUE,'known with pfam',
                                     ifelse(KOorCOG==TRUE & knownPFAM == FALSE,'known without pfam',
                                     ifelse(KOorCOG==FALSE & knownPFAM == TRUE,'unknown with pfam','unknown without pfam')))) %>%
  dplyr::mutate(annotation_quality_strict = ifelse(knownCOG==TRUE & knownPFAM == TRUE,'known with pfam',
                                     ifelse(knownCOG==TRUE & knownPFAM == FALSE,'known without pfam',
                                     ifelse(knownCOG==FALSE & knownPFAM == TRUE,'unknown with pfam','unknown without pfam'))))



#save as txt
#write.table(annotation.quality.tbl, file=paste0(outdir,'annotation_quality_per_locus_tag.txt'),
#            sep='\t',quote = FALSE, row.names = FALSE)

#annotation.quality.tbl <- read.delim(file=paste0('~/DATA/MarbGenomics/Graphs/','annotation_quality_per_locus_tag.txt')) %>%
#  tibble::tibble()



# annotation quality per OG
# uses majority rule to determine annotation quality of each of the OrthoFinder OGs
# if equal, take the conservative option of 'unknown>known'
annotation.q.per.OG <- annotation.quality.tbl %>%
  dplyr::filter(!is.na(OG)) %>%
  dplyr::group_by(OG,annotation_quality) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate( quality_homogeneity= n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(OG) %>%
  dplyr::arrange(desc(n),desc(annotation_quality)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()


annotation.q.per.unique <- annotation.quality.tbl %>%
  filter(is.na(OG)) %>%
  select(locus_tag,annotation_quality) %>%
  mutate(OG=locus_tag) %>%
  select(OG, annotation_quality,-locus_tag) %>%
  mutate(n=1,quality_homogeneity=NA)


#save as txt
#write.table(annotation.q.per.OG, file=paste0('~/DATA/MarbGenomics/Graphs/','annotation_quality_per_OG.txt'),
#            sep='\t',quote = FALSE, row.names = FALSE)

#annotation.q.per.OG <- read.delim(file=paste0('~/DATA/MarbGenomics/Graphs/','annotation_quality_per_OG.txt')) %>%
#  tibble::tibble()



# for OGs with n>1 we derive it from the majority (conservative appraoch)
# while for singletons we need alternative
dummy.softcore <- annotation.q.per.OG %>% filter(OG %in% CORE.ogs) %>% mutate(partition='SoftCore')

# Visualise per pangenome partition
p.summary <- rbind(annotation.q.per.OG,annotation.q.per.unique) %>%
 # filter(OG %in% c(CORE.ogs,accessory.ogs))
  dplyr::mutate(partition = ifelse(OG %in% CORE.ogs, 'Core',
                                   ifelse(OG %in% softcore.ogs, 'SoftCore',
                                          ifelse(OG %in% accessory.ogs, 'Accessory', 'Unique')))) %>%
  dplyr::bind_rows(dummy.softcore) %>%
  dplyr::group_by(partition,annotation_quality) %>% tally() %>%
  dplyr::mutate(partition=factor(partition, levels=rev(c('Core', 'SoftCore','Accessory','Unique')))) %>%
  ggplot2::ggplot(aes(x=partition,y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="fill", stat="identity") + coord_flip()+
  ggplot2::geom_text(data = . %>%
                            dplyr::mutate(p = n / sum(n)) %>%
                            dplyr::ungroup(),
                           aes(y = p, label = scales::percent(p)),
                           position = position_stack(vjust = 0.5),
                           show.legend = FALSE,
                           color='white',
                     size=3)+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  fdb_style(aspect.ratio=0.2)+
  ylab('Percentage of gene families') + xlab('pangenome partition')


ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/annotation_quality_pangenomepartitions.pdf',plot=p.summary, width = 8, height = 4,unit='in')

#plot for overall pangenome annotation qualities
p.overall_ann.q <- rbind(annotation.q.per.OG,annotation.q.per.unique) %>%
  dplyr::group_by(annotation_quality) %>% tally() %>%
  ggplot2::ggplot(aes(x='all',y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="stack", stat="identity") + coord_flip()+
  ggplot2::geom_text(data = . %>%
                       dplyr::mutate(p = n / sum(n)) %>%
                       dplyr::ungroup(),
                     aes(y = p, label = scales::percent(p)),
                     position = position_stack(vjust = 0.5),
                     show.legend = FALSE,
                     color='white',
                     size=3)+
  scale_fill_manual(values=wesanderson::wes_palette(names(wesanderson::wes_palettes[13])))+
  fdb_style(aspect.ratio=0.075)+
  ylab('No of gene families') + xlab('pangenome partition')

#plot for overall paritioning of pangenome
p.overall_partitioning <- rbind(annotation.q.per.OG,annotation.q.per.unique) %>%
  dplyr::mutate(partition = ifelse(OG %in% CORE.ogs, 'Core',
                                   ifelse(OG %in% softcore.ogs, 'SoftCore',
                                          ifelse(OG %in% accessory.ogs, 'Accessory', 'Unique')))) %>%
  dplyr::group_by(partition) %>% tally() %>%
  dplyr::mutate(partition=factor(partition, levels=rev(c('Core','SoftCore', 'Accessory','Unique')))) %>%
  ggplot2::ggplot(aes(x='all',y=n, fill=partition)) +
  ggplot2::geom_bar(position="stack", stat="identity") + coord_flip()+
  ggplot2::geom_text(data = . %>%
                       dplyr::mutate(p = n / sum(n)) %>%
                       dplyr::ungroup(),
                     aes(y = p, label = scales::percent(p)),
                     position = position_stack(vjust = 0.5),
                     show.legend = FALSE,
                     color='white',
                     size=3)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(name='Greys',n=9)[5:9])+
  fdb_style(aspect.ratio=0.075)+
  ylab('No of gene families') + xlab('pangenome partition')

p.comb.q <- ggpubr::ggarrange(p.overall_ann.q,p.summary,p.overall_partitioning,align = 'v',ncol = 1,labels = c('a','b','c'))


ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/annotation_quality_pangenomepartitions_summary.pdf',plot=p.comb.q, width = 8, height = 6,unit='in')


#=================================================================≠##

# aim to construct a classical OG freq vs genome number plot annotated with annotation quality



tst.p <- annotation.quality.tbl %>%
  select(genome,locus_tag,OG) %>%
  #mutate(OG=ifelse(is.na(OG),locus_tag,OG)) %>% mutate(OG = as.character(OG))%>%
  left_join(rbind(annotation.q.per.OG,annotation.q.per.unique) %>% mutate(OG = as.character(OG)), by='OG')  %>%
  select(genome,OG,annotation_quality) %>%
  unique() %>%
  group_by(OG,annotation_quality) %>%
  summarise(n_genomes = n()) %>%
  ungroup() %>%
  group_by(annotation_quality, n_genomes) %>%
  tally() %>% filter(n_genomes>1) %>%

ggplot2::ggplot(aes(x=n_genomes,y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="stack", stat="identity") + #coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  ggplot2::ylab("# of gene families")+
  ggplot2::xlab("# of genomes ")+
  fdb_style(aspect.ratio=0.33)+scale_y_continuous(expand = c(0, 0))#,trans = scales::log10_trans())


ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/annotation_quality_pangenomepartitions_summary_rg.pdf',plot=tst.p, width = 8, height = 6,unit='in')


p.comb.q2 <- ggpubr::ggarrange(heights = c(1,0.5,0.5,0.5),tst.p,p.overall_ann.q,p.summary,p.overall_partitioning,align = 'v',ncol = 1,labels = c('a','b','c','d'))
ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/annotation_quality_pangenomepartitions_summary23.pdf',plot=p.comb.q2, width = 8, height = 8,unit='in')



p.l2 <- annotation.quality.tbl %>%
  select(genome,locus_tag,OG) %>%
  #mutate(OG=ifelse(is.na(OG),locus_tag,OG)) %>% mutate(OG = as.character(OG))%>%
  left_join(rbind(annotation.q.per.OG,annotation.q.per.unique) %>% mutate(OG = as.character(OG)), by='OG')  %>%
  select(genome,OG,annotation_quality) %>%
  unique() %>%
  group_by(OG,annotation_quality) %>%
  summarise(n_genomes = n()) %>%
  ungroup() %>%
  group_by(annotation_quality, n_genomes) %>%
  tally() %>% #filter(n_genomes>1) %>%

  ggplot2::ggplot(aes(x=n_genomes,y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position='fill', stat="identity") + #coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  ggplot2::ylab("# of gene families")+
  ggplot2::xlab("# of genomes ")+
  fdb_style(aspect.ratio=0.15)+scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(limits = c(0,94))
  #theme(legend.position = 'none')
#,trans = scales::log10_trans())



p.l1 <- annotation.quality.tbl %>%
  select(genome,locus_tag,OG) %>%
  #mutate(OG=ifelse(is.na(OG),locus_tag,OG)) %>% mutate(OG = as.character(OG))%>%
  left_join(rbind(annotation.q.per.OG,annotation.q.per.unique) %>% mutate(OG = as.character(OG)), by='OG')  %>%
  select(genome,OG,annotation_quality) %>%
  unique() %>%
  group_by(OG) %>%
  summarise(n_genomes = n()) %>%
  ungroup() %>%
  group_by(n_genomes) %>%
  tally() %>% #filter(n_genomes>1) %>%
  ggplot2::ggplot(aes(x=n_genomes,y=n)) +
  ggplot2::geom_bar(position="stack", stat="identity") + #coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  ggplot2::ylab("# of gene families")+
  ggplot2::xlab("# of genomes ")+
  fdb_style(aspect.ratio=0.3)+
  scale_y_continuous(expand = c(0, 0),trans = scales::log10_trans())+
  scale_x_continuous(limits = c(0,94))
  #theme(legend.position = 'none')


p.comb.q2 <- ggpubr::ggarrange(heights = c(1,0.5,0.5,0.5),
                               p.l1+theme(legend.position = 'none'),
                               p.l2+theme(legend.position = 'none'),
                               p.summary+theme(legend.position = 'none'),
                               p.overall_partitioning+theme(legend.position = 'none'),
                               align = 'v',
                               ncol = 1,
                               axis = 'lr',
                               labels = c('a','b','c','d'))

#p.comb.q2 <- ggpubr::ggarrange(heights = c(1,0.5,0.5,0.5),p.l1,p.l2,align = 'v',ncol = 1,labels = c('a','b','c','d'),axis = 'lr')
ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/annotation_quality_pangenomepartitions_summary3.pdf',plot=p.comb.q2, width = 8, height = 8,unit='in')


#==================================================================#
#### significant lineage associated
#


#Enriched, when either of the raw counts or binary counts hits are significant

p.sig <- HGM_RF_COMB %>%
  dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts,fdr.p.value_depletion_binary,fdr.value_depletion_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
  dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::left_join(annotation.q.per.OG,by=c("feature.id.x"="OG")) %>%
  dplyr::select(feature.id.x, group.x, Enrich_bin,Depl_bin,Depl_raw, Enrich_raw,annotation_quality) %>%
  dplyr::mutate(signEnriched = ifelse(Enrich_bin == TRUE | Enrich_raw==TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(signDepleted = ifelse(Depl_bin == TRUE | Depl_raw==TRUE, TRUE, FALSE)) %>%
  dplyr::filter(signEnriched==TRUE | signDepleted==TRUE) %>%
  dplyr::mutate(direction=ifelse(signEnriched==TRUE,'More','Less')) %>%
  dplyr::group_by(annotation_quality, group.x,direction) %>% tally() %>%
  dplyr::mutate(n=ifelse(direction=='Less',-n,n)) %>%
  #dplyr::mutate(group.x=factor(group.x, levels=rev(paste0('cl',1:12)))) %>%
  dplyr::mutate(group=factor(group.x,levels=rev(c('cl1','cl10','cl12','cl3','cl4','cl5','cl8'))))%>%

  ggplot2::ggplot(aes(x=group,y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="stack", stat="identity") + coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  geom_hline(yintercept = 0)+
  fdb_style(aspect.ratio=0.5)+
  ylab('No of significant gene families') + xlab(' ')


ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/significant_0.85_annotation_quality_pangenomepartitions.pdf',plot=p.sig, width = 6, height = 4,unit='in')



# of accessory gene content, what is lineage associated?

p.sig2 <- HGM_RF_COMB %>%
  dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts,fdr.p.value_depletion_binary,fdr.value_depletion_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
  dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::left_join(annotation.q.per.OG,by=c("feature.id.x"="OG")) %>%
  dplyr::select(feature.id.x, group.x, Enrich_bin,Depl_bin,Depl_raw, Enrich_raw,annotation_quality) %>%
  dplyr::mutate(signEnriched = ifelse(Enrich_bin == TRUE | Enrich_raw==TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(signDepleted = ifelse(Depl_bin == TRUE | Depl_raw==TRUE, TRUE, FALSE)) %>%
  group_by(feature.id.x) %>% summarise(signEnriched=any(signEnriched),
                                                          signDepleted=any(signDepleted))%>%
  group_by(signEnriched,signDepleted) %>% tally() %>%
  ggplot2::ggplot(aes(x='accessory',y=n, fill=interaction(signEnriched,signDepleted))) +
  ggplot2::geom_bar(position="stack", stat="identity") + coord_flip()+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise1'))+
  geom_hline(yintercept = 0)+
  fdb_style(aspect.ratio=0.05)+
  ylab('No of significant gene families') + xlab(' ') + geom_text(
    aes(label = stat(y)),
    stat = 'summary', fun = sum,# vjust = -1
  )


p.comb.q2 <- ggpubr::ggarrange(heights = c(0.25,0.25,0.5),p.sig2,p.sig2,p.sig,align = 'v',ncol = 1,labels = c('a','b','c','d'))
ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/sign_annotation_quality_pangenomepartitions_summary2.pdf',plot=p.comb.q2, width = 8, height = 8,unit='in')





HGM_RF_COMB %>%
  dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts,fdr.p.value_depletion_binary,fdr.value_depletion_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
  dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_bin = ifelse(fdr.p.value_depletion_binary<0.05,TRUE,FALSE)) %>%
  dplyr::mutate(Depl_raw = ifelse(fdr.value_depletion_rawcounts<0.05,TRUE,FALSE)) %>%
  dplyr::left_join(annotation.q.per.OG,by=c("feature.id.x"="OG")) %>%
  dplyr::select(feature.id.x, group.x, Enrich_bin,Depl_bin,Depl_raw, Enrich_raw,annotation_quality) %>%
  dplyr::mutate(signEnriched = ifelse(Enrich_bin == TRUE | Enrich_raw==TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(signDepleted = ifelse(Depl_bin == TRUE | Depl_raw==TRUE, TRUE, FALSE)) %>%
  group_by(feature.id.x,annotation_quality) %>% summarise(signEnriched=any(signEnriched),
                                                          signDepleted=any(signDepleted))%>%
  group_by(signEnriched,signDepleted) %>% tally() %>%
  group_by(signEnriched,signDepleted) %>%
  summarise(freq=n/sum(n))















#=================================================================##


p.summary <- annotation.q.per.OG %>%
  ggplot2::ggplot(aes(x='all',y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="fill", stat="identity") + coord_flip()+
  #ggplot2::geom_text(data = . %>%
  #                     dplyr::mutate(p = n / sum(n)) %>%
  #                     dplyr::ungroup(),
  #                   aes(y = p, label = scales::percent(p)),
  #                   position = position_stack(vjust = 0.5),
  #                   show.legend = FALSE,
  #                   color='white',
  #                   size=3)+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  fdb_style(aspect.ratio=0.2)+
  ylab('Percentage of gene families') + xlab('pangenome partition')

annotation.q.per.OG %>%
  ggplot2::ggplot(aes(x='all',y=n, fill=annotation_quality)) +
  ggplot2::geom_bar(position="fill", stat="identity") + coord_flip()+
  ggplot2::geom_text(data = . %>% ungroup() %>%
                       dplyr::mutate(p = n / sum(n)) %>%
                       dplyr::ungroup(),
                     aes(y = p, label = scales::percent(p)),
                     #position = position_stack(vjust = 0.5),
                     show.legend = FALSE,
                     color='white',
                     size=3)+
  scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))+
  fdb_style(aspect.ratio=0.2)+
  ylab('Percentage of gene families') + xlab('pangenome partition')







annotation.quality.tbl %>%
  filter(!is.na(OG)) %>%
  group_by(OG,annotation_quality) %>%
  summarise(n=n()) %>%
  mutate( quality_homogeneity= n/sum(n)) %>%
  ungroup() %>%
  group_by(OG) %>%
  arrange(desc(n)) %>%
  slice(1) %>% arrange(quality_homogeneity) %>% ggplot(aes(x=quality_homogeneity,y=n)) + geom_point()







annotation.quality.tbl%>%
  filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  mutate(genome = factor(genome, levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  #group_by(knownKO,knownCOG,KOorCOG,knownPFAM) %>% tally()
  group_by(genome, annotation_quality) %>% tally() %>%
  ggplot(aes(x=genome,y=n, fill=annotation_quality)) +
  geom_bar(position="stack", stat="identity") + coord_flip()+ scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))


annotation.quality.tbl%>%
  filter(genome %in% lsTrees.94.rooted[[1]]$tip.label) %>%
  mutate(genome = factor(genome, levels=tip_order(lsTrees.94.rooted[[1]]))) %>%
  #group_by(knownKO,knownCOG,KOorCOG,knownPFAM) %>% tally()
  group_by(genome, annotation_quality) %>% tally() %>%
  ggplot(aes(x=genome,y=n, fill=annotation_quality)) +
  geom_bar(position="stack", stat="identity") + coord_flip()+ scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))


annotation.quality.tbl%>%
  filter(genome %in% lsTrees.94.rooted[[1]]$tip.label) %>%
  mutate(genome = factor(genome, levels=tip_order(lsTrees.94.rooted[[1]]))) %>%
  #group_by(knownKO,knownCOG,KOorCOG,knownPFAM) %>% tally()
  group_by(genome, annotation_quality_strict ) %>% tally() %>%
  ggplot(aes(x=genome,y=n, fill=annotation_quality_strict )) +
  geom_bar(position="stack", stat="identity") + coord_flip()+ scale_fill_manual(values=wesanderson::wes_palette('Moonrise2'))




#===================≠#
#=====================#








#======
plotdata

genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  ggplot(aes(Contamination, genome)) + geom_point(aes(color=quality_class))




genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  #dplyr::filter(quality_class != "low") %>%
  ggplot(aes(Genome_size, Contamination)) +
  geom_point(aes(color=quality_class,shape=ifelse(genus=='Marinobacter','Marinobacter','non-Marinobacter'))) +
               fdb_style()



genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  #dplyr::filter(quality_class != "low") %>%
  ggplot(aes(Genome_size, predicted_genes)) +
  geom_point(aes(color=quality_class)) +
  geom_smooth(aes(group=quality_class, color=quality_class),method='lm')+
  fdb_style()


genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  dplyr::filter(quality_class != "low") %>%
  ggplot(aes(N50_scaffolds, N50_contigs)) + geom_point(aes(color=quality_class))


genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  select(genome, quality_class, Completeness, Contamination, N50_contigs) %>%
  tidyr::pivot_longer(cols=c(Completeness, Contamination, N50_contigs)) %>%
  ggplot(aes(value,genome)) +
    geom_point(aes(color=quality_class)) +
    facet_wrap(~name,ncol=3,scale='free_x')




# Plot genome statistics next to the genome names ~ ordered based on phylogenetic tips
d.p <- genome.tbl %>%
  dplyr::filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  dplyr::mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  dplyr::select(genome, quality_class, Genome_size,Coding_density, predicted_genes,Mean_scaffold_length_bp,scaffolds, GC, Completeness, Contamination, N50_scaffolds) %>%
  dplyr::mutate(Genome_size=Genome_size/1000000) %>%
  dplyr::mutate(N50_scaffolds=N50_scaffolds/1000000) %>%
  dplyr::mutate(Mean_scaffold_length_bp=Mean_scaffold_length_bp/1000000) %>%
  dplyr::mutate(predicted_genes=predicted_genes/1000) %>%
  tidyr::pivot_longer(cols=c(Genome_size, predicted_genes,Coding_density,Mean_scaffold_length_bp, scaffolds, GC, Completeness, Contamination, N50_scaffolds)) %>%
  dplyr::mutate(name = factor(name, levels=c('Genome_size','predicted_genes','GC','Coding_density','Completeness','Contamination','scaffolds','N50_scaffolds','Mean_scaffold_length_bp'))) %>%
  ggplot2::ggplot(ggplot2::aes(value,genome)) +
    ggplot2::geom_point(ggplot2::aes(color=quality_class)) +
    ggplot2::scale_x_continuous(position = 'top')+
    ggplot2::facet_wrap(~name,nrow=1,scale='free_x')+
    fdb_style(aspect.ratio=8)+
    ggplot2::theme(panel.spacing.x = unit(0.5, "lines"),
          #strip.background=element_rect(fill="red"),
          panel.background = element_rect(fill = "lightgrey",
                                          colour = "lightgrey",
                                          size = 0.5, linetype = "solid"),
          axis.text.y = ggplot2::element_text(size =6 , colour = '#000000'),)

d.p
ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/Quality_all_r.pdf',plot=d.p, width = 15, height = 20,unit='in')




#--------------- along phylogeny
p = ggtree(sco.106.concat.tree)+ geom_tiplab(size=2)+theme_tree2()
for(sc in c('Genome_size','predicted_genes','GC','Coding_density','Completeness','Contamination','scaffolds','N50_scaffolds','Mean_scaffold_length_bp')){
  p <- facet_plot(p,
             panel=sc,
             mapping = aes(x = !!as.name(sc), y=y, color= quality_class),
             data=genome.tbl %>%
               mutate(id=genome) %>%
               select(id, everything()) %>%
               dplyr::mutate(Genome_size=Genome_size/1000000) %>%
               dplyr::mutate(N50_scaffolds=N50_scaffolds/1000000) %>%
               dplyr::mutate(Mean_scaffold_length_bp=Mean_scaffold_length_bp/1000000) %>%
               dplyr::mutate(predicted_genes=predicted_genes/1000),
             geom=geom_point,size=1.5)#+ ggnewscale::new_scale()#theme(legend.position = "none")
 # p <- p + ggnewscale::new_scale('color')#theme(legend.position = "none")
  }
p

#p+scale_color_manual(values=phylogroup.colors)

p <- p+ggplot2::theme(panel.spacing.x = unit(0.5, "lines"),
                 panel.background = element_rect(fill = "lightgrey",
                                                 colour = "lightgrey",
                                                 size = 0.5, linetype = "solid"),
                 axis.text.y = ggplot2::element_text(size =6 , colour = '#000000'),)+
  ggplot2::scale_x_continuous(position = 'top')

ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/Quality_all_phylogenetuc.pdf',plot=p, width = 15, height = 8,unit='in')





d <- genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  mutate(genome=factor(genome,levels=tip_order(lsTrees.106.rooted[[1]]))) %>%
  select(genome,TypeStrain, quality_class, Completeness, Contamination, N50_contigs) %>%
  gridExtra::tableGrob()

ggsave('~/DATA/MarbGenomics/Graphs/Quality_allr.pdf',plot=d, width = 15, height = 20,unit='in')




#======================

genome.tbl %>%   filter(quality_class!='low') %>%
ggplot(aes(phylogroup, GC)) + geom_point()

genome.tbl %>%
  arrange(desc(GC)) %>%
  filter(quality_class!='low') %>%
  filter(genome %in% sel.genome) %>%
  select(genome, phylogroup, GC) %>%
  group_by(phylogroup) %>%
  dplyr::summarise(mean=mean(GC), sd=sd(GC))


genome.tbl %>%
  arrange(desc(Coding_density)) %>%
  filter(quality_class!='low') %>%
  filter(genome %in% sel.genome) %>%
  select(genome, phylogroup, Coding_density) %>%
  group_by(phylogroup) %>%
  dplyr::summarise(mean=mean(Coding_density), sd=sd(Coding_density))


genome.tbl %>%
  arrange(desc(Genome_size)) %>%
  filter(quality_class!='low') %>%
  filter(genome %in% sel.genome) %>%
  select(genome, phylogroup, Genome_size) %>%
  group_by(phylogroup) %>%
  dplyr::summarise(mean=mean(Genome_size), sd=sd(Genome_size))



#=============

annotation.tbl %>% filter(grepl('dsr',Name)) %>% group_by(Name,product) %>% tally %>% data.frame
annotation.tbl %>% filter(grepl('cox',Name)) %>% group_by(Name,product) %>% tally %>% data.frame
annotation.tbl %>% filter(grepl('sqr',Name)) %>% group_by(Name,product) %>% tally %>% data.frame

annotation.tbl %>% filter(grepl('tau',Name,ignore.case=TRUE)) %>% group_by(Name,product) %>% tally %>% data.frame


annotation.tbl %>% filter(grepl('sulfite',product,ignore.case=T)) %>% group_by(Name,product) %>% tally %>% data.frame
annotation.tbl %>% filter(grepl('cox',product,ignore.case=T)) %>% group_by(Name,product) %>% tally %>% data.frame



KO_annotated %>% filter(grepl('sulf',DESCRIPTION ,ignore.case=T)) %>% group_by(ko_annotation ,DESCRIPTION) %>% tally()

KO_annotated %>% filter(grepl('urea',DESCRIPTION ,ignore.case=T)) %>% group_by(ko_annotation ,DESCRIPTION) %>% tally()


KO_annotated %>% filter(grepl('dsr',ko_annotation ,ignore.case=T)) %>% group_by(ko_annotation ,DESCRIPTION) %>% tally()


#========================================================================

#set pan
pan = OG.pan[sel.genome,]
pan.presabs <- pan
pan.presabs[pan.presabs>1]<-1

#Runs over all phylogroups, and checks if feautres are shared by all members of groups
shared.df <- NULL
for(phylogrp in phylogroup.names){
  message(phylogrp)
  selected_genomes <- genome.tbl %>%
    dplyr::filter(phylogroup == phylogrp) %>%
    dplyr::select(genome) %>%
    dplyr::pull()

  #uses is.shared (custom function)
  shared <- is.shared(pan, selected_genomes)
  shared.df <-  cbind(shared.df, shared)
}

colnames(shared.df) <- phylogroup.names
shared.df <- data.frame(shared.df)

#select those unique in P1
selUnique <- shared.df[shared.df$'P3'==TRUE & rowSums(shared.df)==1,]
s.ogs <- rownames(selUnique)

annotation.tbl %>% select(genome, Name, COG, gene, product,OG) %>% filter(genome=='Marinobacter_adhaerens_HP15') %>% filter(OG %in% s.ogs)
annotation.tbl %>% select(genome, Name, COG, gene, product,OG) %>% filter(genome=='Marinobacter_adhaerens_HP15') %>% filter(OG %in% s.ogs)






df_counts <- data.frame(
  cbind('feature.id' = colnames(pan),
        'nrOfGenes' = colSums(pan),
        'nrOfSpecies' = colSums(pan.presabs),
        'shared'=rowSums(shared.df)),stringsAsFactors = FALSE) %>%
  mutate(nrOfGenes = as.numeric(as.character(nrOfGenes))) %>%
  mutate(nrOfSpecies = as.numeric(as.character(nrOfSpecies))) %>%
  mutate(shared = as.numeric(as.character(shared)))


p.pan_distr <- df_counts %>%
  ggplot(aes(nrOfSpecies, nrOfGenes)) +
    geom_point(size=0.5, shape=21,aes(color=shared)) +
    fdb_style(aspect.ratio=0.5)+scale_color_distiller(palette='BrBG',aesthetics = 'color')+xlab('number of genomes') + ylab('number of genes')

ggsave('~/DATA/MarbGenomics/Graphs/pan_genome_distribution_shared.pdf',plot=p.pan_distr, width = 5, height = 3,unit='in')











#========================================================================
#Supplementary Figure
plotdata <- genome.tbl %>%
  dplyr::filter(grepl('Marinobacter',genome)) %>%
  dplyr::filter(quality_class != "low")


#use plot.distribution to visualise several genome characteristics
p.d.pg <- plot.distibution(df = plotdata,value = 'predicted_genes', binwidth= 100, xlab='# predicted genes')
p.d.gs <- plot.distibution(df = plotdata,value = 'Genome_size', xlab='Genome size (Mbp)', binwidth= 0.1, div=1000000 )
p.d.gc <- plot.distibution(df = plotdata,value = 'GC', binwidth= 0.25, xlab='%GC content')
p.d.cd <- plot.distibution(df = plotdata,value = 'Coding_density', binwidth= 0.2, xlab='Coding denisty (%)')
p.d.contigs <- plot.distibution(df = plotdata,value = 'contigs', xlab='# contigs', binwidth= 10 )

p.combined <- cowplot::plot_grid(p.d.pg, p.d.gs, p.d.gc, p.d.cd, p.d.contigs,
                   ncol = 1, align = 'v', axis = 'lr')

ggsave('~/DATA/MarbGenomics/Graphs/allMB_combined_genome_char.pdf',plot=p.combined, width = 4, height = 8,unit='in')


#========================================================================


p.pg_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size, y= predicted_genes)) +
  ggplot2::geom_point(shape=21,alpha=.5,aes(size=Coding_density)) +
  ggplot2::geom_smooth(method='lm') +
  ggplot2::facet_wrap(~quality_class) +
  fdb_style()

p.pg_v_gs

p.pg_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size/1000000, y= predicted_genes))+
  ggplot2::geom_smooth(method='lm',color='grey20')+
  ggplot2::geom_point(shape=21,alpha=.5,size=2 ,fill='grey')+
  ggplot2::xlab('Genome size (Mbp)') +
  ggplot2::ylab('# predicted genes') +
  fdb_style()

p.pg_v_gs


p.pg_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size/1000000, y= predicted_genes))+
  ggplot2::geom_smooth(method='lm',color='grey20')+
  ggplot2::geom_point(aes(shape=quality_class,fill=quality_class),alpha=.5,size=2)+
  ggplot2::xlab('Genome size (Mbp)') +
  ggplot2::ylab('# predicted genes') +
  ggplot2::scale_shape_manual(values=c(21,23))+
  ggplot2::scale_fill_manual(values=c('grey80','grey30'))+
  fdb_style()

p.pg_v_gs

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pg_v_gs.pdf',plot=p.pg_v_gs, width = 6, height = 2.5,unit='in')



p.cont_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size/1000000, y= contigs))+
  ggplot2::geom_smooth(method='lm',color='grey20')+
  ggplot2::geom_point(aes(shape=quality_class,fill=quality_class),alpha=.5,size=2)+
  ggplot2::xlab('Genome size (Mbp)') +
  ggplot2::ylab('# contigs') +
  ggplot2::scale_shape_manual(values=c(21,23))+
  ggplot2::scale_fill_manual(values=c('grey80','grey30'))+
  fdb_style()
p.cont_v_gs


ggsave('~/DATA/MarbGenomics/Graphs/allMB_cont_v_gs.pdf',plot=p.cont_v_gs, width = 6, height = 2.5,unit='in')


p.cd_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size/1000000, y= Coding_density))+
  ggplot2::geom_smooth(method='lm',color='grey20')+
  ggplot2::geom_point(aes(shape=quality_class,fill=quality_class),alpha=.5,size=2)+
  ggplot2::xlab('Genome size (Mbp)') +
  ggplot2::ylab('Coding density (%)') +
  ggplot2::scale_shape_manual(values=c(21,23))+
  ggplot2::scale_fill_manual(values=c('grey80','grey30'))+
  fdb_style()
p.cd_v_gs


ggsave('~/DATA/MarbGenomics/Graphs/allMB_cd_v_gs.pdf',plot=p.cd_v_gs, width = 6, height = 2.5,unit='in')


p.gc_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size/1000000, y= GC))+
  ggplot2::geom_smooth(method='lm',color='grey20')+
  ggplot2::geom_point(aes(shape=quality_class,fill=quality_class),alpha=.5,size=2)+
  ggplot2::xlab('Genome size (Mbp)') +
  ggplot2::ylab('%GC') +
  ggplot2::scale_shape_manual(values=c(21,23))+
  ggplot2::scale_fill_manual(values=c('grey80','grey30'))+
  fdb_style()
p.gc_v_gs


ggsave('~/DATA/MarbGenomics/Graphs/allMB_gc_v_gs.pdf',plot=p.gc_v_gs, width = 6, height = 2.5,unit='in')

p.ctn_v_gs <- plotdata %>%
  ggplot2::ggplot(aes(x = Genome_size/1000000, y= Contamination))+
  ggplot2::geom_smooth(method='lm',color='grey20')+
  ggplot2::geom_point(aes(shape=quality_class,fill=quality_class),alpha=.5,size=2)+
  ggplot2::xlab('Genome size (Mbp)') +
  ggplot2::ylab('Contamination') +
  ggplot2::scale_shape_manual(values=c(21,23))+
  ggplot2::scale_fill_manual(values=c('grey80','grey30'))+
  fdb_style()
p.ctn_v_gs


ggsave('~/DATA/MarbGenomics/Graphs/allMB_cnt_v_gs.pdf',plot=p.ctn_v_gs, width = 6, height = 2.5,unit='in')


p.combined2 <- cowplot::plot_grid(p.pg_v_gs,
                                  p.gc_v_gs,
                                  p.cd_v_gs,
                                  p.cont_v_gs,
                                  p.ctn_v_gs,
                                 ncol = 1, align = 'v', axis = 'lr')

ggsave('~/DATA/MarbGenomics/Graphs/allMB_combined2_genome_char.pdf',plot=p.combined2, width = 4, height = 12,unit='in')

outfig <- cowplot::plot_grid(p.combined,
                             p.combined2,
                   ncol = 2, align = 'h', axis = 'lr',rel_widths = c(0.8, 1.4))

ggsave('~/DATA/MarbGenomics/Graphs/allMB_outfig_char.pdf',plot=outfig, width = 9, height = 12,unit='in')



#====================================================================
# PAN GENOME ANALYSIS
# set sel.genome to select the genomes you want to incclude

#Selection of genomes
# genus Marinobacter, not manually removed, not low quality
sel.genome <- genome.tbl %>%
  dplyr::filter(genus=="Marinobacter") %>%
  dplyr::filter( manualy_removed ==FALSE)%>%
  dplyr::filter(quality_class!='low') %>%
  dplyr::select(genome) %>%
  dplyr::pull()

length(sel.genome)

#Extract in high-medium quality genomes
#extract singletons, not assigned to OG
pan.size <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  nrow()

#Assigned to OG
pan.assigned <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter(!is.na(OG)) %>%
  nrow()

pan.singletons <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter(is.na(OG)) %>%
  nrow()

pan.OGs <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter(!is.na(OG)) %>%
  dplyr::select(OG) %>%
  dplyr::pull() %>%
  base::unique() %>%
  base::length()

#return basic output statistics
if((pan.assigned + pan.singletons)/pan.size==1){
  #proportions
  message('total CDS: ',pan.size)
  message('total OGs: ',pan.OGs)
  message('nr of OG assigned CDS: ',pan.assigned)
  message('nr of singletons: ',pan.singletons)
  message('percentage assigned: ',round((pan.assigned/pan.size)*100,2),'%')
  message('percentage singletons: ',round((pan.singletons/pan.size)*100,2),'%')

}else{
  message('warning! something went wrong, you should check')
}


#====================================================================================
# PAN GENOMICS
#===== Build pangenome table including singletons! stored in OG.full.pan
# we coalesce OG and locus tag, if OG is NA (=singleton) then it takes the respective locus tag

#ASSIGN LOCUS TAGS TO SINGLETONS!
pan.tbl <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::mutate(OG = coalesce(OG,locus_tag))

#use tbl2pan function to convert it to a pan table!
OG.full.pan <- tbl2pan(tbl = pan.tbl %>% select(locus_tag, genome, OG, product),
                       feature.col='OG')

# we remove the ones with rowSums 0 (none)
OG.full.pan = OG.full.pan[ rowSums(OG.full.pan)!=0, ]
#sanity check
dim(OG.full.pan)



#### MARKER SET DETECTION #####
#-------------------------------------------------------
#Detect single copy orthologs!
SCO = apply(t(OG.full.pan %>% data.frame)==1,1,all)
SCO = names(SCO[SCO ==TRUE])
pan.scos <- base::length(SCO)



# AS we now detected 741 SCOs we need to extract these and re-run phylogeny within marinobacter!
og.sequences <- tblAnnotation %>% filter(OG %in% SCO) %>% filter(genome %in% sel.genome) %>%
  select(genome,OG,locus_tag,aa_sequence,nuc_sequence) %>% collect()

#write each SCO down
for(i in 1:length(SCO)){
  message('writing ',SCO[i])
  #save faa
  sel.og.seq <- og.sequences %>%
    dplyr::filter(OG == SCO[i]) %>%
    dplyr::select(genome,aa_sequence) %>%
    data.frame()

  f.name = paste0("~/DATA/MarbGenomics/MB_SCO/aa/",SCO[i],'.faa')
  writeFasta(data=sel.og.seq, filename=f.name)

  #save fna
  sel.og.seq <- og.sequences %>%
    dplyr::filter(OG == SCO[i]) %>%
    dplyr::select(genome,nuc_sequence) %>%
    data.frame()

  f.name = paste0("~/DATA/MarbGenomics/MB_SCO/nuc/",SCO[i],'.fna')
  writeFasta(data=sel.og.seq, filename=f.name)
}


#Explore more the SCO's what are they? and export a Supplementary table lis of Marinobacter SCO's
# Ideally, one should report all the phylogenetics stats of it as well, dispersion, discordance etc...
annotation.tbl %>%
  dplyr::filter(OG %in% SCO) %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::select(genome,OG,locus_tag,product)

#Look for consistency in annotation
consistent.sco <- annotation.tbl %>%
  dplyr::filter(OG %in% SCO)  %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::arrange(OG,genome) %>%
  dplyr::select(OG,product) %>%
  dplyr::group_by(OG) %>%
  dplyr::count(product) %>%
  dplyr::filter(n==94)
write.table(consistent.sco,
           paste0('~/DATA/MarbGenomics/','Supplementary_MB_SCO_consistently_annotated','.tsv'),
           sep='\t',
           quote=FALSE)

con.SCO <- consistent.sco %>% select(OG) %>% pull
writeLines(con.SCO, "~/DATA/MarbGenomics/','MB_SCO_consistent_annotation.txt")

bashcopy <- paste0('cp ','SCO_Protein_alignments_mafft_genafpair/',con.SCO,'.faa ','conSCO_Protein_alignments_mafft_genafpair/',con.SCO,'.faa')
writeLines(bashcopy, paste0("~/DATA/MarbGenomics/","bashcopy_run_consistent_annot.txt"))


#look for ribosomal proteins only
ribosomal.sco <- consistent.sco %>%
  dplyr::filter(grepl('30S|50S',product,ignore.case=TRUE)) %>%
  dplyr::filter(!grepl('methyltransferase',product,ignore.case=TRUE)) %>%
  dplyr::arrange(product)%>%
  data.frame()

write.table(ribosomal.sco,
            paste0('~/DATA/MarbGenomics/','Supplementary_MB_SCO_ribosomal_proteins','.tsv'),
            sep='\t',
            quote=FALSE)

rib.SCO <- ribosomal.sco %>% select(OG) %>% pull
writeLines(rib.SCO, "~/DATA/MarbGenomics/','MB_SCO_ribosomal.txt")

bashcopy <- paste0('cp ','SCO_Protein_alignments_mafft_genafpair/',rib.SCO,'.faa ','RibSCO_Protein_alignments_mafft_genafpair/',rib.SCO,'.faa')
writeLines(bashcopy, paste0("~/DATA/MarbGenomics/","bashcopy_run_ribosomal.txt"))




#=======
# after PHIPACK, we can continue filtering
# RUN CODE IN mb-rPHI.R

length(setdiff(SCO,SigPHIOfs$OG))
passed_phipack = setdiff(SCO,SigPHIOfs$OG)
bashcopy <- paste0('cp ','SCO_Protein_alignments_mafft_genafpair2/',passed_phipack,'.faa ','SCO_Protein_alignments_mafft_genafpair/',passed_phipack,'.faa')
writeLines(bashcopy, paste0("~/DATA/MarbGenomics/","bashcopy_filterPHI.txt"))




#======================



#presence absence, and repeat sco detection, should select CORE genes
pan.presabs<-OG.full.pan
pan.presabs[pan.presabs > 0] <- 1
CORE = apply(t(pan.presabs %>% data.frame)==1,1,all)
CORE = names(CORE[CORE ==TRUE])
pan.CORE <-  base::length(CORE)

#store output measures
pan.core.all <- annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>%
  dplyr::filter(OG %in% CORE) %>% nrow()

pan.sco.all <- annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>%
  dplyr::filter(OG %in% SCO) %>% nrow()

#report basic numbers
if(pan.scos * length(sel.genome)==pan.sco.all){
  message('total CORE CDS: ',pan.core.all)
  message('total CORE OGs: ',pan.CORE)
  message('total SCO CDS: ',pan.sco.all)
  message('total SCO OGs: ',pan.scos)
  message('percentage CORE: ',round((pan.core.all/pan.size)*100,2),'%')
  message('percentage SCO: ',round((pan.sco.all/pan.size)*100,2),'%')

  }else{
  message('warning! something went wrong, you should check')
}

#====================================================================

pan = OG.full.pan


pan2represented <- function(pan) {
  out <- data.frame(table(apply(t(pan)>0,1,sum)))
  colnames(out) <- c('genomes','orthogroups')
  return(out)
}

#Examples
OFpan.represented = pan2represented(pan)

p.repr <- OFpan.represented %>%
  ggplot2::ggplot(aes(x = as.numeric(genomes), y = orthogroups)) +
  ggplot2::geom_col(fill='grey30',color='white',size=0.2) +
  #ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::ylab("# of gene families")+
  ggplot2::xlab("# of genomes ")+
  fdb_style(aspect.ratio=0.5)+scale_y_continuous(expand = c(0, 0),trans = scales::log10_trans())

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pangenom_repr_log10.pdf',plot=p.repr, width = 4, height = 3,unit='in')

p.repr <- OFpan.represented %>%
  ggplot2::ggplot(aes(x = as.numeric(genomes), y = orthogroups)) +
  ggplot2::geom_col(fill='grey30',color='white',size=0.2) +
  ggplot2::scale_y_log10(expand = c(0, 0)) +
  ggplot2::ylab("# of gene families")+
  ggplot2::xlab("# of genomes ")+
  fdb_style(aspect.ratio=2)+
  ggplot2::coord_flip()

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pangenom_repr.pdf',plot=p.repr, width = 4, height = 3,unit='in')




OFpan.represented2 = pan2represented(pan2)

p.repr2 <- OFpan.represented2 %>%
  ggplot2::ggplot(aes(x = as.numeric(genomes), y = orthogroups)) +
  ggplot2::geom_col(fill='grey30',color='white',size=0.2) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::ylab("# of gene families")+
  ggplot2::xlab("# of genomes ")+
  fdb_style(aspect.ratio=0.5)

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pangenom_repr_pan2.pdf',plot=p.repr2, width = 4, height = 3,unit='in')


#==========
# Setting up variations of the pangenome, to one exclude singletons, and to only retain HQ genome

# pan2 = removal of singletons
pan2 = pan[,which(colSums(pan)!=1)]

#select hq genomes
sel.genome.hq <- genome.tbl %>%
  dplyr::filter(genus=="Marinobacter") %>%
  dplyr::filter( manualy_removed ==FALSE)%>%
  dplyr::filter(quality_class=='high') %>%
  dplyr::select(genome) %>%
  pull()

#pan3=hq genomes
pan3 = pan[intersect(sel.genome.hq, rownames(pan)), ]
pan3 = pan3[,which(colSums(pan3)!=0)]

dim(pan)
dim(pan2)
dim(pan3)
#==========


#sp <- vegan::specaccum(pan, 'random', permutations=100)
#summary(sp)
#plot(sp, ci.type='poly', col='darkred', lwd=2, ci.lty=0, ci.col='darkred', xlab='Genomes', ylab='Proteins', main='Gene accumulation plot')
#boxplot(sp, col='white', add=TRUE,cex=0.5,pch = 21)
#mods <- vegan::fitspecaccum(sp, "arrh")


#====== pan genome openness
#calculate HEAPS LAW
heaps_law.est1 <- micropan::heaps(pan,n.perm=500)
heaps_law.est2 <- micropan::heaps(pan2,n.perm=500)
heaps_law.est3 <- micropan::heaps(pan3,n.perm=500)

# calculate pangenome rarefaction curves
rar.tbl <- micropan::rarefaction(pan, n.perm = 500)
rar.tbl2 <- micropan::rarefaction(pan2, n.perm = 500)
rar.tbl3 <- micropan::rarefaction(pan3, n.perm = 500)

# cacultating the mean of each permutation to have a line in the plot
rar.mean.tbl <-rar.tbl %>% rsample::gather(key = "Permutation", value = "Clusters", -Genome) %>%
  group_by(Genome) %>%
  dplyr::summarise(mean = mean(Clusters), n = n()) %>%
  ungroup() %>%
  mutate(Clusters=mean)

rar.mean.tbl2 <-rar.tbl2 %>% rsample::gather(key = "Permutation", value = "Clusters", -Genome) %>%
  group_by(Genome) %>%
  dplyr::summarise(mean = mean(Clusters), n = n()) %>%
  ungroup() %>%
  mutate(Clusters=mean)

rar.mean.tbl3 <-rar.tbl3 %>% rsample::gather(key = "Permutation", value = "Clusters", -Genome) %>%
  group_by(Genome) %>%
  dplyr::summarise(mean = mean(Clusters), n = n()) %>%
  ungroup() %>%
  mutate(Clusters=mean)

#tidy data
rar.tbl_l <- rar.tbl %>% rsample::gather(key = "Permutation", value = "Clusters", -Genome) %>% mutate(type='pan-genome')
rar.tbl2_l <- rar.tbl2 %>% rsample::gather(key = "Permutation", value = "Clusters", -Genome) %>% mutate(type='woSing')
rar.tbl3_l <- rar.tbl3 %>% rsample::gather(key = "Permutation", value = "Clusters", -Genome) %>% mutate(type='hq')

#merge data
rar.all.tbl <- rbind(rar.tbl_l,rar.tbl2_l,rar.tbl3_l)

rar.all.tbl %>%
  ggplot(aes(x = Genome, y = Clusters)) +
  geom_line(alpha=0.2,aes(group = Permutation)) +
  geom_line(data=rar.mean.tbl2, aes(x = Genome, y = Clusters,),color='red',size=1.5) +
  fdb_style(aspect.ratio = 0.5)

rar.mean.tbl <-rar.mean.tbl %>% mutate(type='pan-genome')
rar.mean.tbl2 <-rar.mean.tbl2 %>% mutate(type='woSing')
rar.mean.tbl3 <-rar.mean.tbl3 %>% mutate(type='hq')
rar.all.mean.tbl <- rbind(rar.mean.tbl2,rar.mean.tbl,rar.mean.tbl3)

p.singleplot <- rar.all.tbl %>% ggplot(aes(x=Genome, y=Clusters)) +
  geom_line(aes(colour=type,  group=interaction(Permutation, type)),alpha=0.2) +
  geom_line(data=rar.all.mean.tbl, aes(x = Genome, y = Clusters, group=type),color='black',size=0.5) +
  scale_color_manual(values=c('grey50','grey70','grey80'))+
  fdb_style(aspect.ratio = 0.5)
ggsave('~/DATA/MarbGenomics/Graphs/allMB_pangenomCurve.pdf',plot=p.singleplot, width = 4, height = 3,unit='in')


p.facetted <- rar.all.tbl %>% ggplot(aes(x=Genome, y=Clusters)) +
  geom_line(aes(colour=type,  group=interaction(Permutation, type)),alpha=0.2) +
  geom_line(data=rar.all.mean.tbl, aes(x = Genome, y = Clusters, group=type),color='black',size=0.5) +
  scale_color_manual(values=c('grey50','grey70','grey80'))+
  fdb_style(aspect.ratio = 0.5) + facet_wrap(~type)
ggsave('~/DATA/MarbGenomics/Graphs/allMB_pangenomCurve_facet.pdf',plot=p.facetted, width = 4, height = 3,unit='in')


#==========================================
# CORE GENOME CURVE
rmat <- micropan::rarefact(pan, what = 'coregenome', n.perm = 500)

# tidy the rmat data and rbind with the pangenome curve
rar.tbl_cor <- rmat %>% data.frame() %>%
  rownames_to_column(var='Genome') %>%
  pivot_longer(permut_1:permut_500) %>%
  mutate(Permutation=name) %>%
  mutate(Clusters =as.numeric(as.character(value ))) %>%
  select(-value,-name) %>%
  mutate(type='corall') %>%
  select(Genome, Permutation,Clusters,type)

rar.tbl_cor_l <- rar.tbl_cor %>% gather(key = "Permutation", value = "Clusters", -Genome) %>% mutate(type='cor')

#calculate mean of the permutaitons to fit the central line in the plots
rar.mean.cor <-rar.tbl_cor %>%
  group_by(Genome) %>%
  dplyr::summarise(mean = mean(Clusters), n = n()) %>%
  ungroup() %>%
  mutate(Clusters=mean) %>% mutate(type='corall')

#combine data with pan-genome curves, visualise all together
rar.all.mean.tbl2 <- rbind(rar.mean.cor%>% mutate(cor=TRUE),rar.all.mean.tbl%>% mutate(cor=FALSE)) %>% mutate(Genome=as.numeric(as.character(Genome)))

#Vislualise all together
p.cmb <- rar.tbl_cor %>% mutate(cor=TRUE) %>% rbind(.,rar.all.tbl %>% mutate(cor=FALSE)) %>%
  mutate(Clusters =as.numeric(as.character(Clusters ))) %>%
  mutate(Genome =as.numeric(as.character(Genome ))) %>%
  ggplot(aes(x=Genome, y=Clusters)) +
  geom_line(aes(colour=type,  group=interaction(Permutation, type)),alpha=0.2) +
  geom_line(data=rar.all.mean.tbl2, aes(x = Genome, y = Clusters, group=type),color='black',size=0.5) +
  ggtitle(paste0("Heap's Law: intercept=",round(heaps_law.est1[1]),' alpha=',round(heaps_law.est1[2],3))) +
  ggplot2::expand_limits(x = 0, y = 0)+
  fdb_style(aspect.ratio = 0.5)+facet_wrap(~cor,ncol=1,scale='free')+scale_color_manual(values=c('grey60','grey60','grey70','grey80'))

ggsave('~/DATA/MarbGenomics/Graphs/allMB_combined_rarefraction.pdf',plot=p.cmb, width = 5, height = 5,unit='in')


core.n.perm <- rar.tbl_cor %>% filter(Genome==length(sel.genome)) %>% select(Clusters) %>% pull() %>% mean()
pan.n <- rar.tbl_l %>% filter(Genome==length(sel.genome)) %>% select(Clusters) %>% pull() %>% mean()
pan2.n <- rar.tbl2_l%>% filter(Genome==length(sel.genome)) %>% select(Clusters) %>% pull() %>% mean()
pan3.n <- rar.tbl3_l%>% filter(Genome==length(sel.genome.hq)) %>% select(Clusters) %>% pull() %>% mean()

message('nr of predicted core families: ',core.n.perm)
message('pan all nr of families: ',pan.n)
message('pan no singletons nr of families: ',pan2.n)
message('pan HQ nr of families: ',pan3.n)

message(paste0("Heap's Law: intercept=",round(heaps_law.est1[1]),' alpha=',round(heaps_law.est1[2],3)))
message(paste0("Heap's Law: intercept=",round(heaps_law.est2[1]),' alpha=',round(heaps_law.est2[2],3)))
message(paste0("Heap's Law: intercept=",round(heaps_law.est3[1]),' alpha=',round(heaps_law.est3[2],3)))




#==========================================================================
# Binomial mixture models fitted on pan table
#if BIC is last row, we need to increase the K.range
binmix.lst <- binomixEstimate(pan, K.range = 3:25)

binmix.lst$BIC.tbl  %>% slice(which.min(BIC))
ncomp <- binmix.lst$BIC.tbl  %>% slice(which.min(BIC)) %>% select(K.range) %>% pull()

print(binmix.lst$BIC.tbl) # minimum BIC at 3 components
## Not run:
# The pan-genome gene distribution as a pie-chart library(ggplot2)
# prepare
binmix.comb1 <- binmix.lst$Mix.tbl %>%
  filter(Components == ncomp) %>%
  mutate(class='1.pan')
binmix.comb2 <- binmix.lst$Mix.tbl %>%
  filter(Components == ncomp) %>%
  mutate(Mixing.proportion = Mixing.proportion * Detection.prob)%>%
  mutate(class='2.avg')

binmix.comb <- rbind(binmix.comb1, binmix.comb2)

p.pan.dist <- binmix.comb %>%
  ggplot(aes(x = class, y = Mixing.proportion, fill = Detection.prob)) +
  geom_bar(position='fill',color='white',stat = 'identity',size=0.5) +
  geom_text(aes(label = round(Detection.prob,3)),
            position = position_fill(vjust = .5),size=2)+
  labs(x = "", y = "", title = "Pan-genome gene distribution") +
  scale_fill_gradientn(colors = c("forestgreen",'white',"purple")) +
  fdb_style(aspect.ratio=3)+ theme(axis.text.x = element_text(angle = 90))
p.pan.dist

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pangnome_gene_distribution.pdf',plot=p.pan.dist, width = 3, height = 5,unit='in')


#average genome contributions
total <- sum(binmix.comb2$Mixing.proportion)
core.prop <- binmix.comb2[binmix.comb2$Detection.prob ==1,'Mixing.proportion'] %>% pull()
#core genome as fraction of average genome
core.prop/total

softcore.prop <- binmix.comb2[binmix.comb2$Detection.prob > 0.95,'Mixing.proportion'] %>% pull() %>% sum()

#soft core genome as fraction of average genome
softcore.prop/total

singleton.prop <- binmix.comb2[binmix.comb2$Detection.prob < 1/94,'Mixing.proportion'] %>% pull() %>% sum()
#singleton faction of average genome
singleton.prop/total

#singleton in total pan
singleton.inall.prop <- binmix.comb1[binmix.comb1$Detection.prob < 1/94,'Mixing.proportion'] %>% pull() %>% sum()
singleton.inall.prop



#============================================================
# functional enrichment analysis in pangenome sets
# This is not trivial, so think and implement in a smart way

CORE.ogs <- CORE  %>% unique()

singletons_loci <- annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>% dplyr::filter(is.na(OG))%>%
  mutate(OG = coalesce(OG,locus_tag)) %>% select(locus_tag) %>% pull()  %>% unique()

accessory.ogs <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter(!is.na(OG)) %>%
  dplyr::filter(!(OG %in% CORE.ogs)) %>% select(OG) %>% pull() %>% unique()

length(CORE.ogs)
length(singletons_loci)
length(accessory.ogs)

df <- data.frame(c('3.CORE'=length(CORE.ogs),
  '2.Accessory'=length(accessory.ogs),
  '1.singletons'=length(singletons_loci)))
colnames(df) =c('n')

df %>% rownames_to_column(var='type') %>%
  ggplot(aes(x='',y=n,fill=type)) +
  geom_bar(position='fill',color='white',stat = 'identity',size=0.5) + coord_flip() +
  fdb_style()
df<- df %>% mutate(prop=n/(sum(df$n)))

p.pan.classes <- df  %>% rownames_to_column(var='type') %>%
  ggplot2::ggplot(aes(x='',y=n,fill=type)) +
  ggplot2::geom_bar(color='white',stat = 'identity',size=0.5) +
  geom_text(aes(label = round(prop,3)),
            position = position_stack(vjust = .5),size=4) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values=c('grey70','grey50','grey30'))+
  fdb_style(aspect.ratio=0.2)

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pan_c.pdf',plot=p.pan.classes, width = 5, height = 1,unit='in')


#============================================================
#plot the accesory pan genome here
gplots::heatmap.2(as.matrix(pan.presabs[sel.genome, accessory.ogs]),trace='none')
gplots::heatmap.2(as.matrix(remove_rare(pan.presabs[sel.genome, accessory.ogs]),fraction=0.1),trace='none', col=colorRampPalette(c("white", "black"))(n = 2),dendrogram=c("row"))


AccesoryDF = t(AccesoryDF)
#force row order so that it matches the order of leafs in rep_tree_d


#============================================#
#SORT pangenome along tree
#can’t have any branch lengths of zero or downstream commands will collapse those nodes…
tree2 = drop.tip(lsTrees.rooted[[1]], tip =setdiff(lsTrees.rooted[[1]]$tip.label, sel.genome))
tree2$edge.length[which(tree2$edge.length == 0)] <- 0.00001
rep_tree_um <- chronopl(tree2,lambda = 0.1,tol = 0)
rep_tree_d <- as.dendrogram(as.hclust.phylo(rep_tree_um))

AccesoryDF = remove_rare(pan.presabs[sel.genome, accessory.ogs])

clade_order <- order.dendrogram(rep_tree_d)
clade_name <- labels(rep_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(AccesoryDF))
combined_ordered_matrix <- AccesoryDF[new_order,]

heatmap.2(as.matrix(combined_ordered_matrix),trace="none",col = colorRampPalette(c("white", "black"))(n = 2),Rowv= rep_tree_d,labCol = FALSE,margins = c(2, 20),dendrogram="row")



#===============


unique_count <- annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  group_by(genome) %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::filter(is.na(OG)) %>%
  dplyr::tally()



unique_count %>% ggplot(aes(x=n)) + geom_histogram()



p <- ggtree::ggtree(sco.94.concat.tree.rooted) + ggtree::geom_tiplab(align=TRUE,linetype='dashed',linesize=0.3, size=0)

cct <- cazyCat.pan %>%
  tibble::rownames_to_column(var='id')  %>%
  tidyr::pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT') %>%
  dplyr::filter(id %in% sco.94.concat.tree.rooted$tip.label)


metadata <- genome.tbl %>% filter(genome %in% sco.94.concat.tree.rooted$tip.label)# %>% data.frame()
#rownames(metadata) <- metadata$genome
metadata <- metadata %>% mutate(id=genome) %>% select(-genome) %>% select(ID, everything())


habitat.colors =c(Oil = "#FFD422",
                  Photic ="#43997A",
                  Saline_soil = 'grey20' ,
                  Phototroph ="#658E67",
                  Polar ="#5D6795",
                  Terrestial ="#B85F49",
                  Hydrothermal_Vent="#E41A1C",
                  Sediment ="#A35390",
                  salt_lake ='purple',
                  Lake   ='orange',
                  Saltern = 'grey20',
                  Deep  ='black'   ,
                  Tidal_flat ='brown',
                  Other    ='grey',
                  'blue')

unique(metadata$SS) %>% as.character()

habitat.colors  <- c(
  'Phototroph' ="#658E67",
  'Other'='grey',
  'Polar' ="#5D6795",
  'Sediment' ="#A35390",
  #''='blue',
  'Photic' ="#43997A",
  'Saltern' = 'grey20',
  'Oil' = "#FFD422",
  'salt_lake' ='purple',
  'Terrestial' ="#B85F49",
  'Hydrothermal_Vent'="#E41A1C",
  'Lake'   ='orange',
  'Deep'  ='black',
  'Saline_soil' ='brown2',
  'Tidal_flat' ='brown')





p1 <- p %<+% metadata
p1$data <- p1$data %>% mutate(SS = as.character(SS))%>% mutate(SS = ifelse(nchar(SS)<1,NA,SS))

p1<-p1 + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = "white", offset = 0.04,size = 0.02,width=0.03)+
  #theme(legend.position='none')+
  scale_fill_manual(values=phylogroup.colors)+ggnewscale::new_scale_fill()+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=SS),
                          color = "white",
                          offset = 0.04,
                          size = 0.02,
                          width=0.03)+
                          #theme(legend.position='none')+
  scale_fill_manual(values=habitat.colors)+
  ggnewscale::new_scale_fill()
p1

#library(ggstance)
#install.packages('ggstance')
unique.p <- ggtree::facet_plot(p1, panel = '# unique CDS', data = unique_count,
                                  geom = ggstance::geom_barh,
                                  mapping = aes(x = n,),
                                  stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])
unique.p



ggtree::facet_plot(p1, panel = 'CAZY', data = cct,
                   geom = ggstance::geom_barh,
                   mapping = aes(x = value, fill = cazyCAT),
                   stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])


p.facet1<- ggtree::facet_plot(p, panel = 'size', data = metadata,
                   #geom = ggstance::geom_barh,
                   geom = geom_point,
                   mapping = aes(x = Genome_size/1000000,color='Genome_size'),
                   stat='identity' )  #+ theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])


p.facet2<- ggtree::facet_plot(p.facet1, panel = 'GC', data = metadata,
                              geom = geom_point,
                              mapping = aes(x = GC,color='GC'),
                              stat='identity' )  #+ theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])

p.facet3<- ggtree::facet_plot(p.facet2, panel = 'coding_density', data = metadata,
                              geom = geom_point,
                              mapping = aes(x = Coding_density,color='Coding_density'),
                              stat='identity' )  #+ theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])

p.facet4<- ggtree::facet_plot(p.facet3, panel = 'piBias', data = metadata,
                              geom = geom_point,
                              mapping = aes(x = b,color='pi'),
                              stat='identity' )  + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_color_manual(values=ochRe::ochre_palettes[['emu_woman_paired']])



d=data.frame(y=1:94,
             .panel=c('piBias'))
p.facet4 <- p.facet4 + geom_hline(data=d, aes(yintercept=y),size=.5, color='lightgrey') + theme(legend.position = 'none')
p.facet4



ggsave('~/DATA/MarbGenomics/Graphs/p.genome_stats_inMarb.pdf',p.facet4 ,width=7,height = 10,units='in')





#namatjira_qual, namatjira_div mccrea parliament tasmania nolan_ned winmar olsen_qual williams_pilbara healthy_reef, dead_reef galah, lorikeet, jumping_frog emu_woman_paired

#=============


annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  #dplyr::filter((OG %in% CORE.ogs)) %>% #select(COG) %>%
  select(OG, genome) %>%
  #dplyr::left_join(cog_func, by='COG') %>%
  group_by(genome) %>% dplyr::count() %>% plot()


pan.name = "OGpan"
pan.name = "COGpan"
pan.name = "COGpan"
pan.name = "CAZYpan"

sel.col = 'sc0.85'

#Load data
ANI.hmg.grp.tbl <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'hgm_allgroups.tsv'))
ANI.hmg.grp.tbl <- ANI.hmg.grp.tbl %>% as_tibble()
ANI.hmg.grp.tbl$group <- factor(ANI.hmg.grp.tbl$group,levels=paste0('cl',1:length(unique(ANI.hmg.grp.tbl$group))))



s4.metadata %>% filter(sc0.85=='1')
s4.metadata %>% filter(sc0.85=='8')


h_gen <- rbind(c('cl1','Marinobacter_adhaerens_HP15'),
               c('cl2','Marinobacter_guineae_M3B'),
               c('cl3','Marinobacter_hydrocarbonoclasticus_ATCC_49840'),
               c('cl4','Marinobacter_shengliensis_SL013A34A2'),
               c('cl5','Marinobacter_vinifirmus_FB1'),
               c('cl6','Marinobacter_pelagius_114J'),
               c('cl7','Marinobacter_vulgaris_F01'),
               c('cl8','Marinobacter_algicola_DG893'),
               c('cl9','Marinobacter_sp_R17'),
               c('cl10','Marinobacter_sp_LV10R510_8'),
               c('cl11','Marinobacter_piscensis_Abdou3'),
               c('cl12','Marinobacter_psychrophilus_20041')
               )

grp <- 'cl3'
sh_genom <- h_gen[h_gen[,1]==grp,2]
sh_genom

sig_enriched_binary <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.p.value_enriched_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_enriched_raw <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.p.value_enriched_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_depleted_binary <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.p.value_depletion_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_depleted_raw <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.value_depletion_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()


#Cazy-based
annotation.tbl %>%
  dplyr::filter(grepl(paste(unique(c(sig_enriched_binary, sig_enriched_raw)),collapse='|'),cazy_domain)) %>%
  filter(genome==sh_genom) %>%
  select(locus_tag,Name, COG, gene, product, eC_number,OG,kfm_domain,cazy_domain,tcdb_domain) %>%
  left_join(cog_func, by='COG') %>% arrange(gene) %>% data.frame()


#OG BASED
annotation.tbl %>%
  dplyr::filter(OG %in% c(sig_enriched_binary, sig_enriched_raw)) %>%
  filter(genome==sh_genom) %>%
  select(locus_tag,Name, COG, gene, product, eC_number,OG,kfm_domain,cazy_domain,tcdb_domain) %>%
  left_join(cog_func, by='COG') %>% arrange(gene) %>% data.frame()



#OG BASED
annotation.tbl %>%
  dplyr::filter(grepl(paste(unique(c(sig_enriched_binary, sig_enriched_raw)),collapse='|'),cazy_domain)) %>%
  filter(genome==sh_genom) %>%
  select(locus_tag,Name, COG, gene, product, eC_number,OG,kfm_domain,cazy_domain,tcdb_domain) %>%
  left_join(cog_func, by='COG') %>% arrange(gene) %>% data.frame()




# COG specific
cog_func %>% filter(COG %in% c(sig_enriched_binary, sig_enriched_raw))

# COG specific
annotation.tbl %>%
  dplyr::filter(COG %in% c(sig_enriched_binary, sig_enriched_raw)) %>%
  filter(genome==sh_genom) %>%
  select(locus_tag,Name, COG, gene, product, eC_number,OG,kfm_domain,cazy_domain,tcdb_domain) %>%
  left_join(cog_func, by='COG') %>% arrange(gene) %>% data.frame()








annotation.tbl %>% filter(COG %in% c(sig_enriched_binary, sig_enriched_raw)) %>%
  filter(genome=='Mairnobacter_adhaerens_HP15')






CORE.ogs

cog_func
cog_desc


COG_core <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter((OG %in% CORE.ogs)) %>% select(COG) %>%
  dplyr::left_join(cog_func, by='COG') %>% group_by(func) %>% dplyr::count()%>%
  dplyr::mutate(n_core=n) %>%
  dplyr::select(-n)

COG_acc <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter((OG %in% accessory.ogs)) %>% select(COG) %>%
  dplyr::left_join(cog_func, by='COG') %>% group_by(func) %>% dplyr::count()%>%
  dplyr::mutate(n_acc=n) %>%
  dplyr::select(-n)

COG_sing <- annotation.tbl %>%
  dplyr::filter(genome %in% sel.genome) %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter((locus_tag %in% singletons_loci )) %>% select(COG) %>%
  dplyr::left_join(cog_func, by='COG') %>% group_by(func) %>% dplyr::count() %>%
  dplyr::mutate(n_singc=n) %>%
  dplyr::select(-n)

COG.pan.df <- COG_core %>% dplyr::left_join(COG_acc, by='func' )%>%
  dplyr::left_join(COG_sing, by='func') %>% data.frame()


COG.pan.df$func <- as.character(COG.pan.df$func)

COG.pan.df[COG.pan.df$func=='',]
COG.pan.df[is.na(COG.pan.df$func),]

COG.pan.df[is.na(COG.pan.df$func),'func'] <- ''
COG.pan.df[COG.pan.df$func=='','func'] <- c('unassigned','unassigned')

COG.pan.df <- COG.pan.df %>% group_by(func)%>%
  summarise_each(funs(sum)) %>% data.frame()

#tidy up , this is marinobacter specific, needs a function
COG.pan.df[is.na(COG.pan.df$n_singc),'n_singc'] <- c(0,0)
COG.pan.df


#ok now we should be able to build contigency tables!

pan.class.df <- COG.pan.df
rownames(pan.class.df) <- pan.class.df$func
pan.class.df<- pan.class.df[,2:ncol(pan.class.df)]
pan.class.df <- t(pan.class.df)

rw.nm <- 'n_acc'

df <- pan.class.df

#filter the categories, to those that make sence (get rid of unassigned ones?)
#or low characterized
df <- df[,c('J','K','L','D','V','T','M','N','U','O','C','G','E','F','H','I','P','Q','R','S')]

outout <- NULL
for(rw.nm in rownames(df)){

nfeatures.ingroup <- sum(df[rownames(df) == rw.nm, ]) # total features in the selected group
nfeatures.outgroup <- sum(df[rownames(df) != rw.nm, ])	# total features in the non-selected genomes

df.met <- data.frame(cbind('type'=rownames(df)))

  output <- NULL
  for(feature in 1:ncol(df)){
    feature.id <- colnames(df)[feature]
    df.met$feature <- df[,feature]

    # Binary version
    pval <- phyper(
      q = nrow(subset(df.met, type == rw.nm & feature > 0)) - 1,
      m = sum(df.met$feature > 0),
      n = nrow(subset(df.met,feature == 0)),
      k = nrow(subset(df.met, type == rw.nm)),
      lower.tail=FALSE
    )

    score <- -log10(pval)

    pval2 <- phyper(
      q = sum(subset(df.met, type == rw.nm)$feature) - 1,
      m = nfeatures.ingroup,
      n = nfeatures.outgroup,
      k = sum(df.met$feature),
      lower.tail=FALSE
    )

    #Test for depletion binary version
    pval_dep <- phyper(
      q = nrow(subset(df.met, type == rw.nm & feature > 0)),
      m = sum(metadata$feature > 0),
      n = nrow(subset(df.met,feature == 0)),
      k = nrow(subset(df.met, type == rw.nm)),
      lower.tail=TRUE
    )

    score_dep <- -log10(pval_dep)

    #Test for depletion Raw counts version
    pval2_dep <- phyper(q = sum(subset(df.met, type == rw.nm)$feature),
                        m = nfeatures.ingroup, n = nfeatures.outgroup, k = sum(df.met$feature),lower.tail=TRUE)
    score2_dep<- -log10(pval2_dep)

    res <- data.frame(
      feature.id = feature.id,
      score_enriched_binary = score,
      p.value_enriched_binary = pval,
      score_depletion_binary=score_dep,
      p.value_depletion_binary=pval_dep,
      score_enriched_rawcounts = -log10(pval2),
      p.value_enriched_rawcounts = pval2,
      score_depletion_rawcounts=score2_dep,
      p.value_depletion_rawcounts=pval2_dep
    )

    metadata$feature <- NULL
    output <- rbind(output,res)
  }

  output<-data.frame(
    feature.id=output$feature.id,
    group = rw.nm,
    score_enriched_binary=output$score_enriched_binary,
    z.score_enriched_binary= (output$score_enriched_binary - mean(output$score_enriched_binary)) / sd(output$score_enriched_binary),
    p.value_enriched_binary=output$p.value_enriched_binary,
    fdr.p.value_enriched_binary = p.adjust(output$p.value_enriched_binary, method = 'fdr'),

    score_depletion_binary=output$score_depletion_binary,
    z.score_depletion_binary=(output$score_depletion_binary - mean(output$score_depletion_binary)) / sd(output$score_depletion_binary),
    p.value_depletion_binary=output$p.value_depletion_binary,
    fdr.p.value_depletion_binary=p.adjust(output$p.value_depletion_binary, method = 'fdr'),

    score_enriched_rawcounts=output$score_enriched_rawcounts,
    z.score_enriched_rawcounts=(output$score_enriched_rawcounts - mean(output$score_enriched_rawcounts)) / sd(output$score_enriched_rawcounts),
    p.value_enriched_rawcounts=output$p.value_enriched_rawcounts,
    fdr.p.value_enriched_rawcounts =p.adjust(output$p.value_enriched_rawcounts, method = 'fdr'),

    score_depletion_rawcounts=output$score_depletion_rawcounts,
    z.score_depletion_rawcounts=(output$score_depletion_rawcounts - mean(output$score_depletion_rawcounts)) / sd(output$score_depletion_rawcounts),
    p.value_depletion_rawcounts=output$p.value_depletion_rawcounts,
    fdr.value_depletion_rawcounts=p.adjust(output$p.value_depletion_rawcounts, method = 'fdr')
  )

  outout <- rbind(outout,output)

}

outout %>% filter(fdr.p.value_enriched_rawcounts<0.05)

df_plot <- df %>% data.frame() %>% rownames_to_column(var='type') %>% melt(.,id="type") %>% mutate(comb = paste0(type,variable))
df_sign <- outout %>% mutate(comb=paste0(group,feature.id))
df_plot <- df_plot %>% left_join(df_sign,by='comb') %>% mutate(p.value = fdr.p.value_enriched_rawcounts) %>% select(type,variable,value,comb,p.value)

df_plot$stars <- cut(df_plot$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

#calculate proportions per group
df1 <- df_plot %>%  filter(type=='n_core') %>%
  mutate(per=value/sum(value))
df2 <- df_plot %>%  filter(type=='n_acc') %>%
  mutate(per=value/sum(value))
df3 <- df_plot %>%  filter(type=='n_singc') %>%
  mutate(per=value/sum(value))

#merge together
df_plot <- rbind(df1,df2,df3)
df_plot$type <- factor(df_plot$type,levels=rev(c('n_core','n_acc','n_singc')))

#plot
p.hyperg <- df_plot %>%
  ggplot2::ggplot(aes(x=variable, y=type, fill=per))+
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "#2C7BB6", na.value = NA)+
  ggplot2::geom_text(aes(label=stars), color="black", size=5) +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::labs(y=NULL, x=NULL) +
  fdb_style(aspect.ratio =0.15)

ggsave('~/DATA/MarbGenomics/Graphs/allMB_COGCAT_hyperg_pan.pdf',plot=p.hyperg, width = 8, height = 2,unit='in')



#=======================================================================
# Percemtahe of COG categoeies in the core, accessory and singleton gene families

COG.pan.df

c1 <- COG.pan.df %>% filter(!(func %in% c('R','S','unassigned')))%>%select(-func)%>%colSums()
c2 <- COG.pan.df %>% filter(func %in% c('R','S'))%>%select(-func)%>%colSums()
c3 <- COG.pan.df %>% filter(func %in% 'unassigned')%>%select(-func)%>%colSums()

out.df <- rbind('3.annotated'=c1,'2.poorly'=c2,'1.unknown'=c3)

out.df <- out.df %>% data.frame() %>% rownames_to_column(var='type') %>% melt(.,id='type')
out.df$variable <- factor(out.df$variable,rev(c('n_core','n_acc','n_singc')))

percentData <- out.df %>% dplyr::group_by(variable) %>%
  dplyr::mutate(ratio=scales::percent(value/sum(value)))

p.k <- out.df%>%
  ggplot2::ggplot(aes(x=variable,y=value,fill=type)) +
  ggplot2::geom_bar(position='fill',color='white',stat = 'identity',size=0.5) +
  geom_text(data=percentData, aes(y=value,label=ratio),
            position=position_fill(vjust=0.5),color='white',size=3)+
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values=c('grey70','grey50','grey30'))+
  fdb_style(aspect.ratio=0.2)

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pan_COGCATannotationrate.pdf',plot=p.k, width = 6, height = 2,unit='in')


#=======================================================================


#=======================================================================


fluid.all <- micropan::fluidity(pan.matrix=as.matrix(OG.pan),n.sim = 1000)

tracker=c('phylogroup'='all', fluid.all['Mean'] , fluid.all['Std'])

for(ph in levels(genome.tbl$phylogroup)){
  selected_genomes <- genome.tbl %>% filter(phylogroup==ph) %>% select(genome) %>% pull
  fluid <- micropan::fluidity(pan.matrix=as.matrix(OG.pan[selected_genomes,]),n.sim = 1000)
  tracker <- rbind(tracker, c('phylogroup'=ph, fluid['Mean'] , fluid['Std']))

}


tracker <- tracker %>% data.frame
tracker %>%
  mutate(Mean=as.numeric(as.character(Mean)), Std=as.numeric(as.character(Std))) %>%
  ggplot2::ggplot(ggplot2::aes(x=phylogroup,y=Mean, fill = phylogroup)) +
  ggplot2::geom_bar(stat='identity') +
  ggplot2::coord_flip()+
  ggplot2::geom_errorbar(ggplot2::aes(ymin=Mean-Std, ymax=Mean+Std), width=.2,
                position=ggplot2::position_dodge(.9)) + ggplot2::scale_fill_manual(values=c('grey',phylogroup.colors)) + fdb_style(aspect.ratio=1) + ggplot2::scale_y_continuous(expand=c(0,0))





imngs[,4:13]

gtdb <- read.delim('/Users/sa01fd/DATA/MarbGenomics/gtdb/mb.bac120.summary.tsv')
gtdb %>% tibble::tibble() %>% select(user_genome, classification)
head(gtdb)


gtdb %>% tibble::tibble() %>%
  filter(user_genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  select(user_genome, classification) %>%
  group_by(classification) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(n>1) %>%
  ggplot(aes(reorder(classification,n),n)) +
  geom_col() +
  coord_flip()








chao.pansize <- chao(panmat)

fluid <- fluidity(panmat)
geneWeights(panmat, type = c("shell", "cloud"))


ppca <- panPca(panmat)

ggplot(ppca$Evar.tbl) +geom_col(aes(x = Component, y = Explained.variance)) # Plotting scores

ggplot(ppca$Scores.tbl) +
  geom_text(aes(x = PC1, y = PC2, label = GID.tag)) # Plotting loadings
ggplot(ppca$Loadings.tbl) +
  geom_text(aes(x = PC1, y = PC2, label = Cluster))



#=======================================================================
# Explore singletons before ignoring them further



#identify distance to nearest tip
nearest_tip_dist <- COPH %>% data.frame() %>%
  tibble::rownames_to_column(var='genome') %>%
  tidyr::pivot_longer(cols = -genome, names_to='genome_B',values_to='COPH') %>%
  group_by(genome) %>%
  filter(genome != genome_B) %>%
  arrange(COPH) %>% slice(1)


singleton_counts <- annotation.tbl %>% filter(locus_tag %in% singletons_loci) %>% group_by(genome) %>% tally()

singleton_counts %>% left_join(nearest_tip_dist, by='genome') %>%
  left_join(genome.tbl, by='genome') %>%
  ggplot(aes(n,COPH)) + geom_point(aes(size=Contamination),shape=21,fill='grey') +
  fdb_style()+ylab('cophenetic distance to nearest tip')


singleton_counts %>% left_join(nearest_tip_dist, by='genome') %>%
  left_join(genome.tbl, by='genome') %>%
  ggplot(aes(n,COPH)) + geom_point(shape=21,fill='grey') +
  fdb_style()+ylab('cophenetic distance to nearest tip')

nearest_tip_dist



singletons.df <- annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>% dplyr::filter(is.na(OG))%>%
  mutate(OG = coalesce(OG,locus_tag)) %>%
  dplyr::group_by(genome) %>%
  dplyr::count() %>%
  dplyr::mutate(singletons = n) %>%
  dplyr::left_join(genome.tbl,by='genome') %>%
  dplyr::left_join(nearest_tip_dist, by='genome')

p.sing.1 <- singletons.df %>%
  ggplot2::ggplot(aes(x=Genome_size/1000000, singletons))+
  geom_smooth(method='lm',color='black',size=.5)+
  geom_point(shape=21,alpha=0.7,fill='grey50',size=2)+
  ggpubr::stat_cor() + fdb_style() +
  xlab('Genome size (Mbp)') +
  ylab('# of singletons')

p.sing.2 <- singletons.df %>%
  ggplot2::ggplot(aes(x=Contamination, singletons))+
  geom_smooth(method='lm',color='black',size=.5)+
  geom_point(shape=21,alpha=0.7,fill='grey50',size=2)+
  ggpubr::stat_cor() + fdb_style() +  ylab('# of singletons')

p.sing.3 <- singletons.df %>%
  ggplot2::ggplot(aes(x=COPH, singletons))+
  geom_smooth(method='lm',color='black',size=.5)+
  geom_point(shape=21,alpha=0.7,fill='grey50',size=2)+
  ggpubr::stat_cor() + fdb_style() +  ylab('# of singletons')


p.sing.comb <- gridExtra::grid.arrange(p.sing.1,p.sing.2,p.sing.3,ncol=3)
ggsave('~/DATA/MarbGenomics/Graphs/allMB_pan_singletons_cor_2.pdf',plot=p.sing.comb, width = 8, height = 3,unit='in')

p.singletons <- singletons.df %>%
  ggplot2::ggplot(aes(x='', singletons))+
   geom_hline(yintercept = mean(singletons.df$singletons),color='grey',linetype='dashed') +
  geom_hline(yintercept = median(singletons.df$singletons),color='grey') +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=21,alpha=0.7,fill='grey50',size=2)+
  fdb_style(aspect.ratio = 4)+ggtitle(paste0('mean=',round(mean(singletons.df$singletons),2),
                                             ' median=',round(median(singletons.df$singletons),2),
                                             ' sd=',round(sd(singletons.df$singletons),2)))+
  theme(plot.title = element_text(size = 10))+ylab('# of singletons')

ggsave('~/DATA/MarbGenomics/Graphs/allMB_pan_singletons.pdf',plot=p.singletons, width = 6, height = 4,unit='in')

#====
# genome specific singletons, cluster into blocks of genes

outdir = '~/DATA/MarbGenomics/Graphs/'
library(gggenes)


stat.out <- NULL
global.out <- NULL
for(gnm in sel.genome){
  message('analysing ',gnm)
  #Select singleton loci
  gen.spec.loci <- annotation.tbl %>%
    filter(genome %in% sel.genome) %>%
    filter(type=='CDS') %>% dplyr::filter(is.na(OG))%>%
    mutate(OG = coalesce(OG,locus_tag)) %>%
    dplyr::filter(genome==gnm) %>%
    dplyr::select(locus_tag) %>% pull()

  #analyse gene proximity
  prox.clust <- proximity_clusters(annotation.tbl, gen.spec.loci,distance=3000 )

  #extract >2 clusters
  largest.clusts <- prox.clust %>%
    group_by() %>%
    dplyr::count(prox.cluster) %>% filter(n>2) %>% select(prox.cluster) %>% pull()
  filtered <- prox.clust %>% filter(prox.cluster %in% largest.clusts)
  n_in_clusters <- filtered %>% nrow()
  #fill gaps and assign to prox
  filled <- fillgaps2(in.cluster=filtered)
  prox.clust <- proximity_clusters(annotation.tbl, filled$locus_tag,distance=3000 )

  if(!is.null(prox.clust)){
    out.p <- prox.clust %>% mutate(direction=ifelse(strand=="+",1,-1)) %>% mutate(transposase=grepl('transposase',product,ignore.case=TRUE)) %>%
      mutate(singleton = ifelse(locus_tag %in% gen.spec.loci,1,0.5)) %>%
      ggplot(aes(xmin = start, xmax = end, y = seqnames, fill = product,forward=direction,label=Name,alpha=singleton,color=transposase)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
      geom_gene_label(align = "left")+fdb_style(aspect.ratio=0.5) +
      theme(legend.position = "none")+facet_wrap(~ prox.cluster, scales = "free", ncol = 1) + geom_gene_label(align = "left") +
      theme_genes()+
      scale_alpha(range=c(0.5,1))+ ggtitle(paste0(gnm,' singletons'))+theme(legend.position = "none")

    out.tb <- prox.clust %>%  mutate(direction=ifelse(strand=="+",1,-1)) %>%
      mutate(singleton = ifelse(locus_tag %in% gen.spec.loci,1,0.5)) %>%
      select(seqnames, start, end, strand, locus_tag, genome, product, COG, kfm_domain, cazy_domain,prox.cluster,singleton) %>% data.frame

    ggsave(paste0(outdir,'Singletons_',gnm,'.pdf'),plot=out.p, width = 6, height = length(unique(prox.clust$prox.cluster)),unit='in')
    global.out <- rbind(global.out, out.tb)


    stat.out <- rbind(stat.out, c('genome' = gnm,
      'n_loci' = length(gen.spec.loci),
      'n_clustered' = n_in_clusters,
      'prop_clust' = n_in_clusters/length(gen.spec.loci),
      'n_clust' = length(unique(prox.clust$prox.cluster)))
    )

  }
}


global.out <- global.out %>%
  dplyr::left_join(ko2ann, by=c('kfm_domain'='ko')) %>%
  dplyr::left_join(cog_func_single, by='COG') %>%
  #dplyr::left_join(out.mat,by=c('OG'='feature.id.x')) %>%
  dplyr::left_join(df.pfam_per_locus %>% select(-genome), by=c('locus_tag'='seq_id'))



write.table(global.out,file=paste0(outdir,'Singletons_all','.tsv'),sep='\t',row.names = FALSE, quote = FALSE)

stat.out %>% data.frame() %>%
  mutate(prop_clust=as.numeric(as.character(prop_clust))) %>%
  mutate(n_loci=as.numeric(as.character(n_loci))) %>%
  ggplot(aes(x=n_loci,y=prop_clust))+
  geom_smooth(method='lm',color='black',size=.5)+
  geom_point(shape=21,alpha=0.7,fill='grey50',size=2)+
  ggpubr::stat_cor() + fdb_style() +
  xlab('# singleton loci') +
  ylab('proportion that forms clusters')


plotmeta <- metadata %>% left_join(out.colors, by=c('genome'='genome')) %>% data.frame()
plotmeta$group <- paste0('cl',plotmeta[,colcol])
plotmeta$group <- factor(plotmeta$group,levels = c(paste0('cl',1:(length(unique(plotmeta$group))-1)),'clNA'))

colcol = 'sc0.8'
colors <- plotmeta  %>% select(paste0('col.',colcol),colcol, group) %>% unique() %>% data.frame()
grid.colors <- as.character(colors[,1])
names(grid.colors) <- as.character(colors[,'group'])

p.sing1 <- stat.out %>% data.frame() %>%
  mutate(prop_clust=as.numeric(as.character(prop_clust))) %>%
  mutate(n_loci=as.numeric(as.character(n_loci))) %>%
  left_join(plotmeta,by='genome')%>%
  ggplot(aes(x=group,y=n_loci))+
  geom_boxplot(aes(fill=group),alpha=0.3,outlier.size=-1)+
  geom_jitter(aes(fill=group),height=0,width=0.2,size=2,shape=21,alpha=0.7) +
  ylab('# of singletons') +
  xlab('ANI cliques') +
  scale_fill_manual(values=grid.colors) + fdb_style(aspect.ratio=0.5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

p.sing2 <- stat.out %>% data.frame() %>%
  mutate(prop_clust=as.numeric(as.character(prop_clust))) %>%
  mutate(n_loci=as.numeric(as.character(n_loci))) %>%
  left_join(plotmeta,by='genome')%>%
  ggplot(aes(x=group,y=prop_clust))+
  geom_boxplot(aes(fill=group),alpha=0.3,outlier.size=-1)+
  geom_jitter(aes(fill=group),height=0,width=0.2,size=2,shape=21,alpha=0.7) +
  ylab('prop. clustering genes') +
  xlab('ANI cliques') +
  scale_fill_manual(values=grid.colors) + fdb_style(aspect.ratio=0.5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

p.comb <- gridExtra::grid.arrange(p.sing1, p.sing2, ncol=1)
ggplot2::ggsave(filename=paste0(outdir,'Singleton_ANIb_cliques',colcol,'.pdf'),plot=p.comb, width = 3, height = 5,unit='in')






stat.out %>% data.frame() %>%
  mutate(prop_clust=as.numeric(as.character(prop_clust))) %>%
  mutate(n_loci=as.numeric(as.character(n_loci))) %>% filter(n_loci>200)



metadata

global.out %>% dplyr::group_by(product) %>% dplyr::count() %>% arrange(desc(n)) %>% data.frame()

annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>% dplyr::filter(is.na(OG))%>%
  mutate(OG = coalesce(OG,locus_tag))%>% dplyr::group_by(product) %>%
  dplyr::count() %>%
  arrange(desc(n)) %>%
  data.frame() %>% dplyr::filter(grepl('transport',product,ignore.case=TRUE))

annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>%
  dplyr::filter(is.na(OG))%>%
  mutate(OG = coalesce(OG,locus_tag))%>%
  dplyr::group_by(product) %>%
  dplyr::count() %>%
  arrange(desc(n)) %>%
  data.frame() %>% dplyr::filter(grepl('transposase',product,ignore.case=TRUE))

#within cluster singletons
w_sing <- global.out %>% dplyr::filter(singleton==1.0)
#within cluster non-singletons
w_non_sing <- global.out %>% dplyr::filter(singleton==0.5)

w_sing %>% dplyr::group_by(product) %>% dplyr::count() %>% arrange(desc(n)) %>% dplyr::filter(grepl('transposase',product,ignore.case=TRUE)) #%>% select(n) %>% sum()
w_non_sing %>% dplyr::group_by(product) %>% dplyr::count() %>% arrange(desc(n)) %>% dplyr::filter(grepl('transposase',product,ignore.case=TRUE)) #%>% select(n) %>% sum()



#==========

#============================
#
#Evolution of gene repertoire relatedness with time.

# First look at dividing time into units of similarity
# Hierarchical units of similarity
# non-overlapping distance classes etc

#nr_of_non_marinobacter_genomes
n_non_MB <-
  genome.tbl %>% filter(manualy_removed==FALSE) %>%
  group_by(genus) %>%
  dplyr::count() %>%
  filter(genus!='Marinobacter') %>% select(n) %>% pull() %>% sum()

#nnr_of_non_marinobacter_genomes
n_non_MB_nLQ <-
nr_of_non_marinobacter_genomes_nonLOW <-
  genome.tbl %>% filter(manualy_removed==FALSE) %>%
  filter(quality_class !='low') %>%
  group_by(genus) %>%
  dplyr::count() %>%
  filter(genus!='Marinobacter') %>% select(n) %>% pull() %>% sum()



#USING THE ENTIRE DATASET
#========================================
#---- note, we lose Talminaduibacer
# no singleton FALSE!, tp predict all other possible non group comparisons
multi.clique.ANIb <- multiClique(distance=ANIb[sel.genome,sel.genome], cutoffs = c(0.7, 0.75, 0.80, 0.85,0.90,0.95,0.98),no_singletons = FALSE)

#subset, ignore 0.7 when ingroup
multi.clique.ANI <- multi.clique.ANIb[,c('genome','c0.7','c0.75','c0.8', 'c0.85', 'c0.9', 'c0.95', 'c0.98')]
sorted.ANI.cliques <- sortCliques(multi.clique.ANI, tree = tree)

s4.metadata <- genome.tbl %>%
  dplyr::filter(manualy_removed==FALSE) %>%
  dplyr::filter(quality_class !='low') %>%
  left_join(sorted.ANI.cliques, by='genome') #%>%
  #left_join(out.colors, by='genome')

#extract some values to report in text
n_cliques <- s4.metadata %>% filter(genus=='Marinobacter') %>% group_by(c0.95) %>% dplyr::count() %>% nrow()
n_singleton_cliques <- s4.metadata %>% filter(genus=='Marinobacter') %>% group_by(c0.95) %>% dplyr::count() %>% filter(n==1) %>% nrow()
n_cliques
n_singleton_cliques/n_cliques

s4.metadata %>% filter(genus=='Marinobacter') %>% group_by(c0.95) %>% dplyr::count() %>% ggplot(aes_string(x='c0.95','n')) + geom_point()
s4.metadata %>% filter(genus=='Marinobacter') %>% group_by(c0.95) %>% dplyr::count() %>% ggplot(aes_string(x='as.factor(c0.95)')) + geom_histogram(stat='count')
s4.metadata %>% filter(genus=='Marinobacter') %>% group_by(c0.95) %>% dplyr::count() %>% ggplot(aes_string(x='n')) + geom_histogram(stat='count')
s4.metadata %>% filter(genus=='Marinobacter') %>% group_by(c0.95) %>% dplyr::count() %>% filter(n>1) %>% ggplot(aes_string(x='n')) + geom_histogram(stat='count')


sel.genome <- s4.metadata %>%
  dplyr::select(genome) %>%
  dplyr::pull() %>%
  base::unique()


dim(ANIb)
dim(AAI)
dim(COPH)
dim(Ortho.Fraction)
dim(pan.ggr)
dim(distMan)

sim.list <- list(
    'COPH' = 1-COPH[sel.genome,sel.genome],
    'ANIb' = ANIb[sel.genome,sel.genome],
    'AAI' = AAI[sel.genome,sel.genome],
    'OF' =  Ortho.Fraction[sel.genome,sel.genome],
    'GGR' = pan.ggr[sel.genome,sel.genome]
    )#,
    #'Manh' = distMan)
names(sim.list)[1]



mantal.1 = vegan::mantel(sim.list[['COPH']], sim.list[['ANIb']],method='pearson',permutations=999)
mantal.1 = vegan::mantel(sim.list[['COPH']], sim.list[['AAI']],method='pearson',permutations=999)
mantal.1 = vegan::mantel(sim.list[['COPH']], sim.list[['OF']],method='pearson',permutations=999)
mantal.1 = vegan::mantel(sim.list[['COPH']], sim.list[['GGR']],method='pearson',permutations=999)
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['ANIb']],method='pearson',permutations=999)

permutation.data <- vegan::permustats(mantal.1)


summary(mantal.1)
densityplot(perm)
qqmath(perm)



#Genome characteristics distances

gs.v <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(Genome_size) %>% pull
names(gs.v) <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(genome) %>% pull
gs.dist = dist(gs.v, method = "euclidean")

pg.v <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(predicted_genes) %>% pull
names(pg.v) <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(genome) %>% pull
pg.dist = dist(pg.v, method = "euclidean")

cd.v <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(Coding_density) %>% pull
names(cd.v) <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(genome) %>% pull
cd.dist = dist(cd.v, method = "euclidean")

b.v <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(b) %>% pull
names(b.v) <- genome.tbl %>% filter(genome %in% sel.genome) %>% select(genome) %>% pull
b.dist = dist(b.v, method = "euclidean")

#
cumulative.v <- cbind(gs.v,pg.v,cd.v,b.v,b.v)
scale.v = scale(cumulative.v, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.cumv = dist(scale.v, method = "euclidean")


df.distances = data.frame(cbind(
  #'COPH'=c(as.matrix(sim.list[['COPH']])),
  'GGR'=c(as.matrix(sim.list[['GGR']])),
  'ANIb'=c(as.matrix(sim.list[['ANIb']])),
  'Gsizedist'=c(as.matrix(gs.dist)),
  'PGdist'=c(as.matrix(pg.dist[sel.genome,sel.genome])),
  'GCdist'=c(as.matrix(gc.dist[sel.genome,sel.genome])),
  'CDdist'=c(as.matrix(cd.dist[sel.genome,sel.genome])),
  'PIdist'=c(as.matrix(b.dist[sel.genome,sel.genome])),
  'cmltv'=c(as.matrix(dist.cumv[sel.genome,sel.genome]))
  ))

df.distances <- df.distances %>% filter(ANIb!=1)
cop.vs.all = melt(df.distances,id='ANIb')
p3.chardist <- ggplot(cop.vs.all,aes(x= ANIb,y= value)) +
  geom_hex(bins = 35, colour = NA) +
  scale_fill_viridis(option="magma",trans = 'sqrt', name = 'Frequency')+
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 7.5),
    legend.text.align = 1,
    legend.title = element_text(size = 9))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11))+
  facet_wrap(~variable,scales='free_y',strip.position = "left",nrow=1)+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 11))+ylab('')+xlab('ANI similarity')

p3.chardist
ggplot2::ggsave(filename=paste0(outdir,'genome_characteristics_distances','.pdf'),plot=p3.chardist, width = 15, height = 4,unit='in')


#Pearson
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['ANIb']],method='pearson',permutations=999)
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['AAI']],method='pearson',permutations=999)
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['COPH']],method='pearson',permutations=999)

#Spearman
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['ANIb']],method='spearman',permutations=999)
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['AAI']],method='spearman',permutations=999)
mantal.1 = vegan::mantel(sim.list[['GGR']], sim.list[['COPH']],method='spearman',permutations=999)


#======== Global ANALYSIS OF RELATIONSHIPS
# wihtin marinobacter



mantal.1 = vegan::mantel(sim.list[['COPH']], sim.list[['ANIb']],method='pearson',permutations=999)


df.distances = data.frame(cbind(
  'COPH'=c(as.matrix(sim.list[['COPH']])),
  'ANIb'=c(as.matrix(sim.list[['ANIb']])),
  'AAI'=c(as.matrix(sim.list[['AAI']])),
  'OF'=c(as.matrix(sim.list[['OF']])),
  'GGR'=c(as.matrix(sim.list[['GGR']]))))

#remove itself
df.distances <- df.distances %>% filter(ANIb!=1)

cop.vs.all = melt(df.distances,id='COPH')

ggplot(cop.vs.all,aes(x= 1-COPH,y= value)) +
  geom_hex(bins = 35, colour = NA) +
  scale_fill_viridis(option="magma",trans = 'sqrt', name = 'Frequency')+
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 7.5),
    legend.text.align = 1,
    legend.title = element_text(size = 9))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11))+
  #geom_smooth(linetype='dashed',color='grey',size=0.5)+
  geom_smooth(method='lm',linetype='dashed',color='grey',size=0.5)+
  geom_smooth(method='lm',formula=y~poly(x,2),color='grey40',size=0.5)+
  #geom_smooth(method='lm',formula=y~poly(x,3),color='black',size=0.5)+
  facet_wrap(~variable,scales='free_y',strip.position = "left",nrow=1)+
  #           labeller = as_labeller(c(Manhattan = "Manhattan distance", Jaccard = "Jaccard distance",GRR='GRR')))+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 11))+ylab('')+xlab('Cophenetic distance')


ggcorrplot::ggcorrplot(M)

M<-cor(df.distances)
corrplot::corrplot(cor_mat, method="circle")

cor_mat <- cor(df.distances,method='pearson')
cor_pmat  <- ggcorrplot::cor_pmat(df.distances)
corrplot:corrplot(M)
ggcorrplot::ggcorrplot(M, p.mat = cor_pmat, hc.order = TRUE, type = 'lower')





#======== ANALYSIS OF INTRA- VS INTER GROUP VARIATION

grp.name <- "c0.95"

grps <- c('genus','c0.75','c0.8','c0.85','c0.9','c0.95','c0.9','c0.95')

outout <- NULL
for(grp.name in grps){
  out <- NULL
  dist.name <- names(sim.list)[1]
  df.dist <- sim.list[[dist.name]]
  df.dist <- df.dist[sel.genome,sel.genome]
  out <- interaDistances(df.dist = df.dist,
                             var_name = names(sim.list)[1],
                             metadata = s4.metadata,
                             grp = grp.name)

  for(dist.name in names(sim.list)[2:length(names(sim.list))]){
    df.dist <- sim.list[[dist.name]]
    df.dist <- df.dist[sel.genome,sel.genome]
    intDist <- interaDistances(df.dist = df.dist,
                               var_name = dist.name,
                               metadata = s4.metadata,
                               grp = grp.name)
    out <- cbind(out, dist.name = intDist$value)
  }
  #out <- out %>% dplyr::mutate(dist.name = factor(out$dist.name, levels=names(sim.list)))
  colnames(out)[colnames(out) %in% c('value','dist.name')] = names(sim.list)
  outout <- rbind(outout, out %>% mutate(grp.name=grp.name))

}


outout %>%
  #dplyr::filter(group.x !='clNA') %>%
  #dplyr::filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>%
  ggplot2::ggplot(aes(ANIb,GGR))+
  ggplot2::geom_point(aes(color=groupin.2),shape=21) +
  ggplot2::geom_smooth(aes(color=groupin.2),method='lm')+
  fdb_style()

outout %>%
  dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>%
  ggplot2::ggplot(aes(1-COPH,ANIb))+
  ggplot2::geom_smooth(aes(color=grouping),method='lm')+
  facet_wrap(~grp.name,nrow=1)+
  ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
  fdb_style()




#===================
# Extract some numbers and tables
#ANI 85 groups with more than 1 type strain
s4.metadata %>% filter(sc0.85==12) %>% select(Organism, strain, TypeStrain, genome)
s4.metadata %>% filter(sc0.85==1) %>% select(Organism, strain, TypeStrain, genome)

s4.metadata %>% group_by(phylogroup,sc0.85, TypeStrain) %>% tally() %>% data.frame %>% filter(TypeStrain=='TRUE')
#======
# Generates mutliple plots,
# x-axis is 1-cophenetic distance
# grouping is used to subset the data according to inter- intra groups
x.s <- "1-COPH"
col.s <- 'grouping'

#double coloured hexbins
p.list <- list()
for(y.s in c('ANIb','AAI','OF','GGR')){
  p.tmp <- outout %>% dplyr::filter(itself==FALSE) %>%
    filter(grepl('Marinobacter',genome)) %>% # !!!!!!!! CONDISDER
    filter(grepl('Marinobacter',genome_b)) %>% # !!!!!!!! CONDISDER
    filter(grp.name!='genus') %>%
    dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>% {
      ggplot2::ggplot(.) +
        geom_hex(data = dplyr::filter(., grouping == 'Inter'),aes_string(x.s,y.s), bins = 35) +
        scale_fill_gradientn(colors = brewer.pal(9,"Blues")[3:9],trans = 'sqrt')+
        ggnewscale::new_scale("fill")+
        geom_hex(data = dplyr::filter(., grouping == 'Intra'),aes_string(x.s,y.s), bins = 35) +
        scale_fill_gradientn(colors = brewer.pal(9,"Oranges")[3:9],trans = 'sqrt')+
        ggplot2::geom_smooth(data=.,aes_string(x.s, y.s, color=col.s),method='lm')+
        #ggplot2::geom_smooth(data=.,aes_string(x.s, y.s, color=col.s), method = 'nls', formula = y ~ a*exp(b *-x), se = FALSE, start = list(a=1,b=0.01), linetype = 2)+
        facet_wrap(~grp.name,nrow=1)+
        ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
        fdb_style() + theme( legend.justification = c(1, 1),
                             legend.key.width = unit(0.12, 'cm'),
                             legend.key.height = unit(0.25, 'cm'),
                             legend.text = ggplot2::element_text(size = 5),
                             legend.title = ggplot2::element_text(size = 5))
      }
  p.list[[y.s]] <- p.tmp

  }

p.comb <- gridExtra::grid.arrange(grobs = p.list, ncol = 1)
ggplot2::ggsave(filename=paste0(outdir,'divergence_all_',y.s,'vs_',x.s,'.pdf'),plot=p.comb, width = 6, height = 6,unit='in')

#==================
# For each grouping level, store the edges (pair-wise comparison) that is within a clique

Intra.c0.7 <- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='c0.7') %>% filter(grouping=='Intra') %>% select(concatenated) %>% pull()
Intra.c0.75 <- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='c0.75') %>% filter(grouping=='Intra') %>% select(concatenated) %>% pull()
Intra.c0.8 <- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='c0.8') %>% filter(grouping=='Intra') %>% select(concatenated) %>% pull()
Intra.c0.85 <- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='c0.85') %>% filter(grouping=='Intra') %>% select(concatenated) %>% pull()
Intra.c0.9 <- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='c0.9') %>% filter(grouping=='Intra') %>% select(concatenated) %>% pull()
Intra.c0.95 <- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='c0.95') %>% filter(grouping=='Intra') %>% select(concatenated) %>% pull()
Inter.genus <-outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  filter(grp.name=='genus') %>% filter(grouping=='Inter') %>% select(concatenated) %>% pull()


y.trans <- 'log'
y.trans <- NULL

#y.s = 'GGR'
p.list <- list()
for(y.s in c('ANIb','AAI','OF','GGR')){
  if(y.trans=='log'){
    y.s <- paste0('log(',y.s,')')
  }
  p.tmp<- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
    dplyr::filter(itself==FALSE) %>% filter(grepl('Marinobacter',genome)) %>%
    dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>%
    dplyr::arrange(desc(groupin.2)) %>% {
      ggplot(.) +
        ggplot2::geom_point(data=filter(.,groupin.2 == 'Inter_genus'),aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.75') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.8)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.8') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.85)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.85') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.9)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.9') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.95)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.95') %>% filter(grouping=='Intra') ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        geom_smooth(data=. ,aes_string(x.s,y.s),linetype='dashed',color='grey',size=0.5)+
        geom_smooth(data=. ,aes_string(x.s,y.s),method='lm',linetype='dashed',color='grey',size=0.5)+
        geom_smooth(data=. ,aes_string(x.s,y.s),method='lm',formula=y~poly(x,2),color='grey40',size=0.5)+
        geom_smooth(data=. ,aes_string(x.s,y.s),method='lm',formula=y~poly(x,3),color='black',size=0.5)+
        scale_color_viridis(discrete = TRUE)+
        fdb_style()
    }
  p.list[[y.s]] <- p.tmp
}

p.comb <- gridExtra::grid.arrange(grobs = p.list, ncol = 1)
p.comb
ggplot2::ggsave(filename=paste0(outdir,'divergence_',y.s,'vs_',x.s,'.pdf'),plot=p.comb, width = 6, height = 8,unit='in')


#=====
# Compare lm models (linear, vs polinomial vs non linear)
#nlm_2 <- nls(bone~a*(1-exp(-c*age)),start=list(a=0.7,c=-2),data=data)
#nlm_2 <- nls(bone~a*(1-exp(-c*age)),start=list(a=120,c=0.064),data=data)
#summary(nlm_2)

#data$predictions <- predict(nlm_2) #assign the values from the model to a column in data sheet
#ggplot(data,aes(x=age,y=bone))+
#  geom_line(aes(y=predictions), colour="blue")+
#  geom_point()+theme_few()+ylab("Jaw Bone Length") + xlab("Age")



#

#
#y.s = 'GGR'
p.list <- list()
for(y.s in c('ANIb','AAI','OF','GGR')){
  p.tmp<- outout %>% mutate(concatenated = paste0(genome,'_', genome_b)) %>%
    dplyr::filter(itself==FALSE) %>% filter(grepl('Marinobacter',genome)) %>%
    dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>%
    dplyr::arrange(desc(groupin.2)) %>% {
      ggplot(.) +
        #ggplot2::geom_point(data=filter(.,groupin.2 == 'Inter_genus'),aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.75') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.8)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.8') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.85)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.85') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.9)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.9') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.95)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        ggplot2::geom_point(data= . %>% filter(grp.name=='c0.95') %>% filter(grouping=='Intra') ,aes_string(x.s,y.s, color='groupin.2'),shape=21,size=0.3) +
        #geom_smooth(data=. ,aes_string(x.s,y.s),linetype='dashed')+
        #scale_colour_brewer()+
        scale_color_viridis(discrete = TRUE)+
        fdb_style()
    }
  p.list[[y.s]] <- p.tmp
}

p.comb <- gridExtra::grid.arrange(grobs = p.list, ncol = 1)
p.comb
ggplot2::ggsave(filename=paste0(outdir,'divergence_onlyMB_',y.s,'vs_',x.s,'.pdf'),plot=p.comb, width = 6, height = 8,unit='in')








outout %>%  mutate(concatenated = paste0(genome,'_', genome_b)) %>%
  dplyr::filter(itself==FALSE) %>% filter(grepl('Marinobacter',genome)) %>%
  dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>% {
    ggplot(.) +
      ggplot2::geom_point(data= . %>% filter(grp.name=='c0.95') %>% filter(grouping=='Intra') ,aes_string(x.s,y.s, color='groupin.2'),shape=21) +
      ggplot2::geom_point(data= . %>% filter(grp.name=='c0.9') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.95)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21) +
      ggplot2::geom_point(data= . %>% filter(grp.name=='c0.85') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.9)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21) +
      ggplot2::geom_point(data= . %>% filter(grp.name=='c0.8') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.85)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21) +
      ggplot2::geom_point(data= . %>% filter(grp.name=='c0.75') %>% filter(grouping=='Intra') %>% filter(!(concatenated %in% Intra.c0.8)) ,aes_string(x.s,y.s, color='groupin.2'),shape=21) +
      geom_smooth(data=.%>%filter(!(concatenated %in% Inter.genus)) ,aes_string(x.s,y.s),linetype='dashed')+
      fdb_style()
  }



  #no clue
outout %>%
  dplyr::filter(group.x !='clNA') %>%
  dplyr::filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::mutate(groupin.2 = paste0(grouping,'_',grp.name)) %>%
  dplyr::filter(groupin.2 != 'Inter_genus') %>%
  dplyr::filter(groupin.2 != 'Intra_genus') %>%
  dplyr::filter(groupin.2 != 'Intra_c0.8') %>%
  dplyr::filter(groupin.2 != 'Intra_c0.85') %>%
  dplyr::filter(groupin.2 != 'Intra_c0.9') %>%
  ggplot2::ggplot(aes(1-COPH,GGR))+
  ggplot2::geom_point(aes(color=groupin.2),shape=21) +
  ggplot2::geom_smooth(aes(color=groupin.2),method='lm')+
  fdb_style()



#=======

outout <- NULL
for(grp.name in grps){
  out <- NULL
  for(dist.name in names(sim.list)){
    df.dist <- sim.list[[dist.name]]
    df.dist <- df.dist[sel.genome,sel.genome]
    intDist <- interaDistances(df.dist = df.dist,
                               var_name = dist.name,
                               metadata = s4.metadata,
                               grp = grp.name)
    out <- rbind(out, intDist %>% mutate(dist.name=dist.name))
    out <- out %>% dplyr::mutate(dist.name = factor(out$dist.name, levels=names(sim.list)))
  }
  p.violin <- out %>%
    dplyr::filter(group.x !='clNA') %>%
    dplyr::filter(group.y !='clNA') %>%
    dplyr::filter(itself==FALSE) %>%
    dplyr::filter(!is.na(grouping)) %>%
    ggplot2::ggplot(aes(grouping, value)) +
    ggplot2::geom_violin(trim=FALSE,aes(fill=grouping,color=grouping),size=0.25,alpha=0.4)+
    ggplot2::geom_boxplot(aes(fill=grouping),alpha=0.7,outlier.shape=NA,width=.2,size=0.5)+
    ggplot2::ylab('simmilarity') +
    ggplot2::xlab(grp.name) +
    ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
    ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
    facet_wrap(~dist.name,nrow=1,scales='free')+
    fdb_style(aspect.ratio=2)
  ggplot2::ggsave(filename=paste0(outdir,'intraViolin_distances',grp.name,'.pdf'),plot=p.violin, width = 6, height = 4,unit='in')


  p.hist <- out %>%
    filter(group.x !='clNA') %>%
    filter(group.y !='clNA') %>%
    dplyr::filter(itself==FALSE) %>%
    dplyr::filter(!is.na(grouping)) %>% ggplot(aes(value))  +
    geom_histogram(binwidth = 0.01,aes(fill=grouping), color = 'white', alpha = .5, position="identity")+
    geom_density(aes(y=0.01 * ..count.., color=grouping, fill=grouping),alpha=.3) +
    scale_y_continuous(expand=c(0,0)) +
    ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
    ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
    facet_wrap(~dist.name,ncol=1,scales='free')+
    fdb_style(aspect.ratio=0.5) +
    theme(strip.placement = "outside")
  ggplot2::ggsave(filename=paste0(outdir,'intraHist_distances',grp.name,'.pdf'),plot=p.hist, width = 5, height = 8,unit='in')

  outout <- rbind(outout, out %>% mutate(grp.name=grp.name))
}

outout <- outout %>% dplyr::mutate(dist.name = factor(outout$dist.name, levels=names(sim.list)))


p.violin <- outout %>%
  dplyr::filter(group.x !='clNA') %>%
  dplyr::filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot2::ggplot(aes(grouping, value)) +
  ggplot2::geom_violin(trim=FALSE,aes(fill=grouping,color=grouping),size=0.25,alpha=0.4)+
  ggplot2::geom_boxplot(aes(fill=grouping),alpha=0.7,outlier.shape=NA,width=.2,size=0.5)+
  ggplot2::ylab('simmilarity') +
  ggplot2::xlab(grp.name) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
  ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
  facet_grid(dist.name~grp.name,scales='free')+
  fdb_style(aspect.ratio=2)
ggplot2::ggsave(filename=paste0(outdir,'ALL_intraViolin_distances',grp.name,'.pdf'),plot=p.violin, width = 6, height = 6,unit='in')


p.hist <- outout %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>% ggplot(aes(value))  +
  geom_histogram(binwidth = 0.01,aes(fill=grouping), color = 'white', alpha = .5, position="identity")+
  geom_density(aes(y=0.01 * ..count.., color=grouping, fill=grouping),alpha=.3) +
  scale_y_continuous(expand=c(0,0)) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
  ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
  facet_grid(grp.name~dist.name,scales='free')+
  fdb_style(aspect.ratio=0.5) +
  theme(strip.placement = "outside")

ggplot2::ggsave(filename=paste0(outdir,'ALL_intraHist_distances',grp.name,'.pdf'),plot=p.hist, width = 5, height = 8,unit='in')




intraGenus <- interaDistances(df.dist = df.dist,
                              var_name = dist.name,
                              metadata = s4.metadata,
                              grp = grp.name)
intraGenus%>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(group.x, value)) +
  geom_boxplot(aes(color=group.x),outlier.shape=NULL) +
  geom_jitter(aes(color=group.x, fill=group.x),shape=21,width=0.2,alpha=0.5) + facet_wrap(~grouping)+
  fdb_style()

intraGenus

intraGenus %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(as.factor(group.x), value)) +
  geom_boxplot(aes(color=group.x),outlier.shape=NULL) +
  geom_jitter(aes(color=group.x, fill=group.x),shape=21,width=0.2,alpha=0.5) + facet_wrap(~grouping)+
  fdb_style()

intraGenus %>%
  dplyr::filter(group.x !='clNA') %>%
  dplyr::filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot2::ggplot(aes(grouping, value,fill=grouping)) +
  ggplot2::geom_jitter(aes(color=grouping),shape=21,alpha=0.5,size=1)+
  ggplot2::geom_boxplot(aes(fill=grouping),alpha=0.7,outlier.shape=NA,width=.2)+
  ggplot2::ylab(dist.name) +
  ggplot2::xlab(grp.name) +
  fdb_style(aspect.ratio=3)

intraGenus %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot2::ggplot(aes(grouping, value)) +
  ggplot2::geom_violin(trim=FALSE,aes(fill=grouping),alpha=0.4)+
  ggplot2::geom_boxplot(aes(fill=grouping),alpha=0.7,outlier.shape=NA,width=.2)+
  ggplot2::ylab(dist.name) +
  ggplot2::xlab(grp.name) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00')) +
  ggplot2::facet_wrap(~group.x, nrow=1)+
  fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


intraGenus %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  dplyr::filter(group.x=='Marinobacter') %>%
  ggplot2::ggplot(aes(grouping, value)) +
  ggplot2::geom_jitter(aes(color=grouping),shape=21,alpha=0.5,size=1)+
  ggplot2::geom_violin(trim=FALSE,aes(fill=grouping),alpha=0.4)+
  ggplot2::geom_boxplot(aes(fill=grouping),alpha=0.7,outlier.shape=NA,width=.2)+
  #stat_summary(aes(fill=grouping),alpha=0.6,fun.data="mean_sdl", geom="crossbar", width=0.1 )+
  ggplot2::ylab(dist.name) +
  ggplot2::xlab(grp.name) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00')) +
  ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00')) +
  ggplot2::facet_wrap(~group.x, nrow=1)+
  fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


intraGenus %>%
  dplyr::filter(group.x !='clNA') %>%
  dplyr::filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot2::ggplot(aes(grouping, value)) +
  ggplot2::geom_violin(trim=FALSE,aes(fill=grouping),alpha=0.4)+
  ggplot2::geom_boxplot(aes(fill=grouping),alpha=0.7,outlier.shape=NA,width=.2)+
  ggplot2::ylab(dist.name) +
  ggplot2::xlab(grp.name) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
  fdb_style(aspect.ratio=2)






intraGenus %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(value)) +
  geom_histogram(aes(fill=grouping),alpha=0.7, position="identity")+
  fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))



#allowed count, instead of density
intraGenus %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>% ggplot(aes(value))  +
  geom_histogram(binwidth = 0.01,aes(fill=grouping), color = 'white', alpha = .5, position="identity")+
  geom_density(aes(y=0.01 * ..count.., color=grouping, fill=grouping),alpha=.3) +
  scale_y_continuous(expand=c(0,0)) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
  ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
  fdb_style(aspect.ratio=0.5)

#Desnity instead
intraGenus %>%
  filter(group.x !='clNA') %>%
  filter(group.y !='clNA') %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(value))  +
  geom_histogram( aes(y = ..density.., fill=grouping), color = 'white', alpha = .5, position="identity") +
  geom_density(aes(fill=grouping,color=grouping),alpha = .3)+
  scale_y_continuous(expand=c(0,0)) +
  ggplot2::scale_fill_manual(values=c('#56B4E9','#E69F00'))+
  ggplot2::scale_color_manual(values=c('#56B4E9','#E69F00'))+
  fdb_style(aspect.ratio=0.5)# +



#===========

#testDF = taxData %>% unique()

testDF=  genome.tbl %>% left_join(sorted.ANI.cliques) %>%
  filter(genome %in% sco.106.concat.tree$tip.label) %>%
  select(sc0.75:sc0.98,genome) %>%
  data.frame() %>%
  unique()

#estDF= s4.metadata %>% select(sc0.7:sc0.98,genome) %>% data.frame() %>% unique()
testDF$size = 2
testDF $pathString <- paste(testDF[,1],
                            testDF[,2],
                            testDF[,3],
                            testDF[,4],
                            testDF[,5],
                            testDF[,6],
                            testDF[,7],
                           testDF[,8],
                            sep = "/")
testDF <- testDF %>% left_join(s4.metadata,by='genome')

testDF $pathString <- gsub('\\/NA','',testDF $pathString)

population <- data.tree::as.Node(testDF)
edges <- data.tree::ToDataFrameNetwork(population, "name")

mygraph <- graph_from_data_frame( edges )

# PLot final ! #tidy a litte as this is quite a thing
#------------------------------------------------

p.circpak <- ggraph(mygraph, layout = 'circlepack') +
  geom_node_circle(data= . %>%
                     mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                     left_join(genome.tbl,by='genome') %>%
                     mutate(phylogroup = as.character(phylogroup)) %>%
                     mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(is.na(grp)),
                   aes(fill = as.factor(depth))#%>% filter(genome%in%sco.94.concat.tree$tip.label) %>% select(.ggraph.index,genome,phylogroup,name)
  ) + scale_fill_manual(values=RColorBrewer::brewer.pal(n=9,name='Greys')[2:7]) + ggnewscale::new_scale_fill() +
  geom_node_circle(data= . %>%
                     mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                     left_join(genome.tbl,by='genome') %>%
                     mutate(phylogroup = as.character(phylogroup)) %>%
                     mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(!is.na(grp)),
                   aes(fill = grp )#%>% filter(genome%in%sco.94.concat.tree$tip.label) %>% select(.ggraph.index,genome,phylogroup,name)
  ) +
  scale_fill_manual(values=phylogroup.colors,na.value='white') + coord_fixed()+
  geom_node_text(data= . %>%
                   mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                   left_join(genome.tbl,by='genome') %>%
                   mutate(phylogroup = as.character(phylogroup)) %>%
                   mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(TypeStrain==TRUE) %>% mutate(TS = 'T') %>% mutate(Tamil = grepl('Tamilnadui',genome)),
                 aes(label=TS,color=Tamil),
                 size=3) +
  scale_color_manual(values=c('black','red'))+
  theme_void()


ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",'ANI_circlepacked_final_annotated_phylogroup','.pdf'),plot=p.circpak, width = 7, height = 7,unit='in')






testDF= genome.tbl %>% left_join(sorted.AAI.cliques) %>%
  filter(genome %in% sco.106.concat.tree$tip.label) %>%
  filter(grepl('Marinobacter|Tamilnadui',genome)) %>%
  select(sc0.75:sc0.98,genome) %>%
  data.frame() %>%
  unique()

#estDF= s4.metadata %>% select(sc0.7:sc0.98,genome) %>% data.frame() %>% unique()
testDF$size = 2
testDF $pathString <- paste(testDF[,1],
                            testDF[,2],
                            testDF[,3],
                            testDF[,4],
                            testDF[,5],
                            testDF[,6],
                            testDF[,7],
                            testDF[,8],
                            sep = "/")
testDF <- testDF %>% left_join(s4.metadata,by='genome')

testDF $pathString <- gsub('\\/NA','',testDF $pathString)

population <- data.tree::as.Node(testDF)
edges <- data.tree::ToDataFrameNetwork(population, "name")

mygraph <- graph_from_data_frame( edges )



# PLot final ! #tidy a litte as this is quite a thing
#------------------------------------------------

p.circpak <- ggraph(mygraph, layout = 'circlepack') +
  geom_node_circle(data= . %>%
                     mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                     left_join(genome.tbl,by='genome') %>%
                     mutate(phylogroup = as.character(phylogroup)) %>%
                     mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(is.na(grp)),
                   aes(fill = as.factor(depth))#%>% filter(genome%in%sco.94.concat.tree$tip.label) %>% select(.ggraph.index,genome,phylogroup,name)
  ) + scale_fill_manual(values=RColorBrewer::brewer.pal(n=9,name='Greys')[2:7]) + ggnewscale::new_scale_fill() +
  geom_node_circle(data= . %>%
                     mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                     left_join(genome.tbl,by='genome') %>%
                     mutate(phylogroup = as.character(phylogroup)) %>%
                     mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(!is.na(grp)),
                   aes(fill = grp )#%>% filter(genome%in%sco.94.concat.tree$tip.label) %>% select(.ggraph.index,genome,phylogroup,name)
  ) +
  scale_fill_manual(values=phylogroup.colors,na.value='white') + coord_fixed()+
  geom_node_text(data= . %>%
                   mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                   left_join(genome.tbl,by='genome') %>%
                   mutate(phylogroup = as.character(phylogroup)) %>%
                   mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(TypeStrain==TRUE) %>% mutate(TS = 'T') %>% mutate(Tamil = grepl('Tamilnadui',genome)),
                 aes(label=TS,color=Tamil),
                 size=3) +
  scale_color_manual(values=c('black','red'))+
  theme_void()


ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",'AAI_circlepacked_final_annotated_phylogroup','.pdf'),plot=p.circpak, width = 7, height = 7,unit='in')





















p.circpak <- ggraph(mygraph, layout = 'linear') +
  geom_node_circle(data= . %>%
                     mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                     left_join(genome.tbl,by='genome') %>%
                     mutate(phylogroup = as.character(phylogroup)) %>%
                     mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(is.na(grp)),
                   aes(fill = as.factor(depth))#%>% filter(genome%in%sco.94.concat.tree$tip.label) %>% select(.ggraph.index,genome,phylogroup,name)
  ) + scale_fill_manual(values=RColorBrewer::brewer.pal(n=9,name='Greys')[2:7]) + ggnewscale::new_scale_fill() +
  geom_node_circle(data= . %>%
                     mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                     left_join(genome.tbl,by='genome') %>%
                     mutate(phylogroup = as.character(phylogroup)) %>%
                     mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(!is.na(grp)),
                   aes(fill = grp )#%>% filter(genome%in%sco.94.concat.tree$tip.label) %>% select(.ggraph.index,genome,phylogroup,name)
  ) +
  scale_fill_manual(values=phylogroup.colors,na.value='white') + coord_fixed()+
  geom_node_text(data= . %>%
                   mutate(genome = sub(".*/", "", gsub("\\/2","",name))) %>%
                   left_join(genome.tbl,by='genome') %>%
                   mutate(phylogroup = as.character(phylogroup)) %>%
                   mutate(grp = ifelse(grepl('P',phylogroup),as.character(phylogroup),NA)) %>% filter(TypeStrain==TRUE) %>% mutate(TS = 'T'),
                 aes(label=TS),size=3) + theme_void()





ggraph(mygraph, layout = 'linear') + geom_edge_arc()

+ geom_node_circle() +
  geom_node_circle(aes(fill = category )) + scale_fill_manual(values=phylogroup.colors)





ggraph(mygraph, layout = 'circlepack') +
  geom_node_circle(aes(fill = as.factor(depth),alpha = as.factor(depth), color = as.factor(depth) )) +
  scale_fill_manual(values=c("0" = "white", "1" = magma(7)[1], "2" = magma(7)[2], "3" = magma(7)[3], "4"=magma(7)[4],"5"=magma(7)[5],"6"=magma(7)[6],"6"=magma(7)[7])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black",'5'="black",'7'='black') ) +
  theme_void() +
  theme(legend.position="FALSE")


l <- ggraph(mygraph, layout = 'partition', circular = TRUE)
l + geom_node_arc_bar(aes(fill = depth)) +
  coord_fixed()

l + geom_edge_diagonal() +
  geom_node_point(aes(colour = depth)) +
  coord_fixed()

graph_tree <- ggraph(mygraph, 'tree') +
  geom_edge_diagonal()
ggplot2::ggsave(filename=paste0(outdir,'ANI_finescale_tree_mb_only_with singletons',grp.name,'.pdf'),plot=graph_tree, width = 5, height = 3,unit='in')


ggraph(mygraph, layout='dendrogram', circular=TRUE) +
  geom_edge_diagonal() +
  theme_void() +
  theme(legend.position="none")

ggraph(mygraph, layout='dendrogram', circular=FALSE) +
  geom_edge_diagonal() +
  theme_void() +
  theme(legend.position="none")

ggraph(mygraph, 'treemap', weight = 'size') +
  geom_node_tile(aes(fill = depth), size = 0.25) +
  theme_void() +
  theme(legend.position="none")

ggraph(mygraph, 'partition', circular = TRUE) +
  geom_node_arc_bar(aes(fill = depth), size = 0.25) +
  theme_void() +
  theme(legend.position="none")

p.graph_edge_link <- ggraph(mygraph) +
  geom_edge_link() +
  geom_node_point(size=.2) +
  theme_void() +
  theme(legend.position="none")
ggplot2::ggsave(filename=paste0(outdir,'ANI_finescale_tree_levels_mb_only_singletons',grp.name,'.pdf'),plot=p.graph_edge_link, width = 5, height = 3,unit='in')

p.graph_edge_link2 <- ggraph(mygraph, layout = 'dendrogram') +
  geom_edge_diagonal() +
  geom_node_point(size=.2) +
  theme_void() +
  theme(legend.position="none")+
  geom_node_text(aes( label=name, filter=leaf) , angle=90 , hjust=1, nudge_y = -0.04)+
  geom_node_point(aes(filter=leaf, size=value, color=group) , alpha=0.6)

ggplot2::ggsave(filename=paste0(outdir,'ANI_finescale_tree_levels_mb_only_singletons2',grp.name,'.pdf'),plot=p.graph_edge_link2, width = 5, height = 3,unit='in')






# create a vertices data.frame. One line per object of our hierarchy
vertices = data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) ,
  value = runif(length(unique(c(as.character(edges$from), as.character(edges$to)))))
)
# Let's add a column with the group of each name. It will be useful later to color points
vertices$group = edges$from[ match( vertices$name, edges$to ) ]


#Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
#calculate the ANGLE of the labels
vertices$id=NA
myleaves=which(is.na( match(vertices$name, edges$from) ))
nleaves=length(myleaves)
vertices$id[ myleaves ] = seq(1:nleaves)
#vertices$angle= 90 - 360 * vertices$id / nleaves

#vertices$name[ grepl('Marinobacter',vertices$name) ]
basename(as.character(vertices$name[ grepl('Marinobacter',vertices$name) ]))

vertices <- vertices %>% mutate(name2=basename(as.character(name))) %>%
  left_join(s5.metadata,by=c('name2'='genome'))  %>% mutate(genome=name2)


# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
#vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
#vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Create a graph object
mygraph <- graph_from_data_frame( edges, vertices=vertices )

#mycolor <- colormap(colormap = colormaps$viridis, nshades = 6, format = "hex", alpha = 1, reverse = FALSE)[sample(c(1:6), 10, replace=TRUE)]

# Make the plot
p.annotated.ggraph <- ggraph(mygraph, layout = 'dendrogram') +
  geom_edge_diagonal(colour="grey") +
  scale_edge_colour_distiller(palette = "RdPu") +
  geom_node_text(aes( filter = leaf, label=genome, colour=as.factor(sc0.8)), size=2.7, alpha=1,angle=90) +
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=as.factor(sc0.8), size=GC, alpha=0.2)) +
  scale_size_continuous( range = c(0.1,7) ) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-10, 10))
p.annotated.ggraph

ggplot2::ggsave(filename=paste0(outdir,'ANI_finescale_ggraph_annotated_mb_only_singletons2','.pdf'),plot=p.annotated.ggraph, width = 10, height = 5,unit='in')














#======


MCR.pan <- read.table("~/DATA/MarbGenomics/MCR.pan.tsv",sep='\t')




#==============
# DMDA DDDL
lc <- read.delim('/Users/sa01fd/DATA/MarbGenomics/Graphs/operons/Marinobacter_algicola_DG893_cluster52_clusterannotated.txt')
Marinobacter_hydrocarbonoclasticus_ATCC_49840_cluster91_clusterannotated

lc %>% filter(prox.cluster=='cluster45') %>% left_join(ko2ann, by=c('kfm_domain'='ko'))

s.g <- annotation.tbl %>%
  filter(gene == 'dmdA') %>%
  left_join(genome.tbl, by='genome') %>%
  select(genome, gene, product, phylogroup, host) %>%
  data.frame()

s.g$genome

s.g1 <- annotation.tbl %>%
  filter(gene == 'dddL') %>%
  left_join(genome.tbl, by='genome') %>%
  select(genome, gene, product, phylogroup, host) %>%
  data.frame()

s.g1$genome



intersect(s.g$genome,s.g1$genome)

gplots::venn(list('dmdA'=s.g$genome,'dddL'=s.g1$genome))

IOFGMAMH_02740:IOFGMAMH_02770



lc <- read.delim('/Users/sa01fd/DATA/MarbGenomics/Graphs/operons/Marinobacter_hydrocarbonoclasticus_ATCC_49840_cluster91_clusterannotated.txt')

lc %>% filter(prox.cluster=='cluster93') %>% left_join(ko2ann, by=c('kfm_domain'='ko'))

annotation.tbl %>% filter(locus_tag %in% paste0('IOFGMAMH_027',40:70)) %>%
  dplyr::left_join(cog_func, by='COG') %>%
  mutate(direction=ifelse(strand=="+",1,-1)) %>%
  #(SELOG = ifelse(OG%in%mapOG,'anchor',ifelse(OG%in%focusOG,'focal',ifelse(OG%in%addOGs,'t6ss','z')))) %>%
  ggplot(aes(xmin = start, xmax = end, y = genome, fill = func,forward=direction,label=Name)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  geom_gene_label(align = "left")+
  scale_color_manual(values=c('blue','red','black','white'))+fdb_style(aspect.ratio=1) + scale_fill_viridis(discrete = TRUE)



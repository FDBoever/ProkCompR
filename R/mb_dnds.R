#==== FDB_ DSDN

#install.packages(c("reshape2", "seqinr"), repos='http://cran.us.r-project.org')
library(seqinr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

infile = args[1]
groupfile = args[2]
group_type = args[3]
outfile = args[4]

groups = read.table(groupfile, sep='\t', row.names='locus', header=TRUE)

is_in_diff_group <- function(s1, s2){
  result <- rep(0L, length(s1))
  for (i in seq_along(s1)){
    result[i] = groups[as.character(s1[i]), group_type] != groups[as.character(s2[i]), group_type]
  }
  return(result)
}



infile = '/Users/sa01fd/DATA/OrthoFinderMarinobacter/SCO_pal2nal_asFasta/OG0000643.pal2nal.fasta'
filelist <- list.files(path='/Users/sa01fd/DATA/OrthoFinderMarinobacter/SCO_pal2nal_asFasta/',patter='.fasta',full.names = TRUE)

out.all <- NULL
for(infile in filelist){
  prot = strsplit(basename(infile), split='.', fixed=TRUE)[[1]][[1]]
  message(prot)
  align_nt = seqinr::read.alignment(infile, format='fasta')
  kaks(align_nt) -> result
  omega = as.matrix(result$ka / result$ks)
  ka = as.matrix(result$ka)
  ks = as.matrix(result$ks)

  # turn into upper triangular matrix
  omega[upper.tri(omega)] <- NA
  ka[upper.tri(ka)] <- NA
  ks[upper.tri(ks)] <- NA

  # turn into long format
  omega = melt(omega)
  ka = melt(ka)
  ks = melt(ks)

  # remove NA rows
  #omega = na.omit(omega)
  #ka = na.omit(ka)
  #ks = na.omit(ks)

  colnames(omega) = c('locus1', 'locus2', 'dnds')
  colnames(ka) = c('locus1', 'locus2', 'ka')
  colnames(ks) = c('locus1', 'locus2', 'ks')


  omega = merge(omega, ka, by=c('locus1', 'locus2'))
  omega = merge(omega, ks, by=c('locus1', 'locus2'))
  omega = na.omit(omega)

  #omega = transform(omega, diff_group=is_in_diff_group(locus1, locus2))
  omega = transform(omega, prot=prot)

  out.all <- rbind(out.all, omega)
}


loc2genome <- annotation.tbl %>% select(genome, locus_tag)

out.ann.all <- out.all %>%
  dplyr::left_join(loc2genome, by=c('locus1'='locus_tag')) %>%
  dplyr::rename(genome1 = genome)

out.ann.all <- out.ann.all %>%
  dplyr::left_join(loc2genome, by=c('locus2'='locus_tag')) %>%
  dplyr::rename(genome2 = genome)


out.ann.all %>% filter(dnds>=0) %>% ggplot(aes(dnds,ks)) +geom_point()

out.ann.all %>% filter(dnds>1) %>% ggplot(aes(dnds,ks)) +geom_point()
out.ann.all %>% filter(dnds>1) %>% filter(dnds != Inf) %>% select(genome1, genome2, dnds, prot)



out.ann.all %>%
  filter(!is.na(genome1)) %>%
  filter(!is.na(genome2))


#inter an intra group analysis
sorted.ANI.cliques[sorted.ANI.cliques$genome=='Marinobacter_sp_DSM_26671','sc0.85']
groupdf <- sorted.ANI.cliques %>% select(genome,'sc0.85')
colnames(groupdf) = c('genome','group')

grouped.ann.all <- out.ann.all %>%
  filter(!is.na(genome1)) %>%
  filter(!is.na(genome2)) %>%
  left_join(groupdf, by=c('genome1'='genome')) %>%
  left_join(groupdf, by=c('genome2'='genome')) %>%
  filter(!is.na(group.x)) %>% filter(!is.na(group.y)) %>%
  mutate(comp = ifelse(group.x == group.y,'intra','inter'))


grouped.ann.all %>% filter(prot=='OG0000643') %>%
  ggplot(aes(comp,ka))+geom_violin()+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
  #facet_wrap(~group.x,nrow=1)+
  fdb_style()

grouped.ann.all %>% filter(prot=='OG0000643') %>% filter(ks<5)  %>%
  ggplot(aes(comp,ks))+geom_violin()+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
 # facet_wrap(~group.x,nrow=1)+
  fdb_style()

grouped.ann.all %>% filter(prot=='OG0000643') %>% filter(dnds != Inf)   %>%
  ggplot(aes(comp,dnds))+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+geom_boxplot()+
  # facet_wrap(~group.x,nrow=1)+
  fdb_style()

#grouped.ann.all %>%  filter(ka<2) %>%
 # ggplot(aes(comp,ka))+geom_violin()+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
  #facet_wrap(~group.x,nrow=1)+
  #fdb_style()


lowinterka.og <- grouped.ann.all %>%  filter(ka<0.01) %>% filter(comp=='inter') %>% select(prot) %>% unique() %>% pull %>% as.character()
grouped.ann.all %>% filter(ka<2) %>% filter(prot %in% lowinterka.og) %>%
  ggplot(aes(comp,ka))+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
  facet_wrap(~group.x,nrow=1)+
  fdb_style()

#medians per gene
med.genes <- grouped.ann.all %>% group_by(comp,prot) %>%filter(ks<2) %>% filter(ka<5) %>%
  dplyr::summarise(dnds_med = median(dnds),ks_med=median(ks),ka_med=median(ka))

med.genes %>%
  ggplot(aes(comp,ka_med)) +geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+  fdb_style()
med.genes %>%
  ggplot(aes(comp,ks_med)) +geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+  fdb_style()
med.genes %>%
  ggplot(aes(comp,dnds_med)) +geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+  fdb_style()


median_strain_comp <- out.ann.all %>%
  filter(!is.na(genome1)) %>%
  filter(!is.na(genome2)) %>%
  filter(dnds != Inf) %>%
  filter(ks<2) %>% filter(ka<5)%>%
  #filter(dnds>0) %>%
  group_by(genome1, genome2) %>%
  dplyr::summarise(dnds_med = median(dnds),ks_med=median(ks),ka_med=median(ka))%>%
  filter(genome1!=genome2)%>%
  left_join(groupdf, by=c('genome1'='genome')) %>%
  left_join(groupdf, by=c('genome2'='genome')) %>%
  filter(!is.na(group.x)) %>% filter(!is.na(group.y)) %>%
  mutate(comp = ifelse(group.x == group.y,'intra','inter'))

median_strain_comp %>% ggplot(aes(comp,ks_med))+geom_violin()+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
  #facet_wrap(~group.x,nrow=1)+
  fdb_style()
median_strain_comp %>% ggplot(aes(comp,ka_med))+geom_violin()+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
  #facet_wrap(~group.x,nrow=1)+
  fdb_style()
median_strain_comp %>% ggplot(aes(comp,dnds_med))+geom_violin()+geom_jitter(width=.1,shape=21,size=1,fill='grey',alpha=.5)+
  #facet_wrap(~group.x,nrow=1)+
  fdb_style()







median_strain_comp <- out.ann.all %>%
  filter(!is.na(genome1)) %>%
  filter(!is.na(genome2)) %>%
  filter(!is.na(dnds))%>%
  #filter(dnds <5) %>%
  filter(dnds != Inf) %>%
  #filter(dnds>0) %>%
  filter(genome1!=genome2) %>%
  group_by(genome1, genome2) %>%
  dplyr::summarise(dnds_med = median(dnds)) %>%
  ungroup() %>%pivot_wider(names_from = genome2, values_from = dnds_med) %>% data.frame()


rownames(median_strain_comp) <- median_strain_comp$genome1
median_strain_comp <- median_strain_comp[,2:ncol(median_strain_comp)]
m <- as.matrix(median_strain_comp)
m[upper.tri(m)] <- t(m)[upper.tri(m)]
gplots::heatmap.2(m,trace='none',margin=c(10,10))


dnds.pan <- m



out.ann.all %>% #filter(dnds>1) %>%
  filter(dnds != Inf) %>%
  select(genome1, genome2, dnds, prot) %>%
  ggplot(aes(prot,dnds)) +geom_point()



selection <- out.ann.all %>% filter(dnds>1) %>% filter(dnds != Inf) %>% select(prot) %>% unique() %>% pull()
annotation.tbl %>% filter(OG %in% selection) %>% select(OG, product) %>% unique() %>% pull()

ggplot(aes(dnds,ks)) +geom_point()

#write.csv(omega, outfile, row.names=FALSE)



out.ann.all %>%
  #filter(dnds>1) %>%
  ggplot(aes(dnds,ks)) +
  geom_point()



#===========

align_nt = seqinr::read.alignment(infile, format='fasta')

prot = strsplit(basename(infile), split='.', fixed=TRUE)[[1]][[1]]
message(prot)

# PERCENT SEQUENCE IDENTITY!!!
align_nt = bio3d::read.fasta(infile)
ide.mat <- bio3d::seqidentity(align_nt)
bio3d::plot.dmat(ide.mat,
          main="Sequence Identity", xlab="Structure No.",
          ylab="Structure No.")

#Prokka pseudogenes


prokka_pseudogenes = read.delim('~/DATA/MarbGenomics/Prokka_Pseudogenes.txt' , header=FALSE)
colnames(prokka_pseudogenes) = c('genome','c1','call','product','seqname','position','c2','c3')


prokka_pseudogenes %>% select(genome, call, product, seqname,position) %>%
  group_by(genome) %>%
  dplyr::count(call) %>%
  filter(n<100) %>%
  ggplot(aes(n))+
    geom_histogram()+
    fdb_style()


prokka_pseudogenes %>% select(genome, call, product, seqname,position) %>%
  group_by(genome) %>%
  dplyr::count(call) %>%
  dplyr::left_join(metadata, by=c('genome'='genome')) %>%
  select(genome, predicted_genes, group, n) %>%
  ggplot(aes(x=group,y=n/predicted_genes,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  fdb_style(aspect.ratio=0.5)


prokka_pseudogenes %>%
  dplyr::select(genome, call, product, seqname,position) %>%
  dplyr::group_by(genome) %>%
  dplyr::count(call) %>%
  dplyr::left_join(metadata, by=c('genome'='genome')) %>%
  dplyr::select(genome, Genome_size, group, n) %>%
  filter(Genome_size<7000000) %>%
  ggplot2::ggplot(aes(x=Genome_size/1000000,y=n))+
  geom_smooth(method='lm',color="red")+
  geom_point(size=2,shape=21,fill='grey50',alpha=0.5)+
  xlab('Genome size (Mbp)')+
  ylab('nr of putative pseudogenes')+
  fdb_style(aspect.ratio=1)



#   SOME CAZY SPECIFIC? FUNCTIONS
#
# obtain cazy.pan from hmm scripts



#-----------------------------------------------------------------------------#
#   UTILITY FUNCTIONS   ####
#-----------------------------------------------------------------------------#


# CAZYpan2CAZYCAT
#=================================================================================#
#' condense CAZY pan into its broad categories, GH, GT, PL, CE, CBM
#'
#' @param CAZYpan a CAZY pan table (genomes, features)
#'
#' @return
#' @export
#'
#' @examples
#' cazyCat.pan <-CAZYpan2CAZYCAT(cazy.pan)

CAZYpan2CAZYCAT <- function(CAZYpan){
  cazyGroups = c()
  for( i in c('GH','GT','PL','CE','CBM')){
    cazyGroups <- cbind(cazyGroups , rowSums(CAZYpan[,grepl(i,colnames(CAZYpan))]))
  }
  colnames(cazyGroups) <- c('GH','GT','PL','CE','CBM')
  cazyGroups <- data.frame(cazyGroups)
  cazyGroups$name <- rownames(cazyGroups)
  cazyGroups <- cazyGroups[,1:(ncol(cazyGroups)-1)]
  return(cazyGroups)
}

#=================================================================================#
# Some code that displays it function, also visualisable as such

cazyCat.pan <-CAZYpan2CAZYCAT(cazy.pan)

# SOM INTERESTING PAN MINGLING GOING ON HERE

# SHOW GH SPECIFICALLY
#cazyCat.pan = cazyCat.pan[rownames(metadata),]
cazyCat.pan %>% rownames_to_column(var='genome') %>%
  left_join(metadata,by=c('genome'='genome')) %>%
  ggplot(aes(x=group,y=GH,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  fdb_style(aspect.ratio=0.5)

# SHOW all categories
cazyCat.pan %>% rownames_to_column(var='genome') %>%
  left_join(metadata,by=c('genome'='genome')) %>%
  pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT') %>%
  ggplot(aes(x=group,y=value,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  facet_wrap(~cazyCAT,scales='free_y',ncol=1)+
  fdb_style(aspect.ratio=0.3)

# SHOW all bloody domains
cazyCat.pan %>% rownames_to_column(var='genome') %>%
  left_join(metadata,by=c('genome'='genome')) %>%
  pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT') %>%
  select(genome,cazyCAT, group, value) %>%
  group_by(genome) %>%
  mutate(total = sum(value)) %>% distinct(genome, .keep_all=TRUE)%>%
  ungroup() %>%
  ggplot(aes(x=group,y=value,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  ylab('total CAZY domains detected')+
  fdb_style(aspect.ratio=0.5)


plot.dat <- cazyCat.pan %>% rownames_to_column(var='id') %>%
  left_join(genome.tbl,by=c('id'='genome')) %>%
  pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT')


plot.dat %>%
  ggplot(aes(x=id,y=value,fill=cazyCAT)) +
  geom_col() #+
#  geom_jitter(shape=21,size=2,width = 0.2)



p <- ggtree(sco.94.concat.tree.rooted) + geom_tiplab(align=TRUE,linetype='dashed',linesize=0.3, size=0)

cct <- cazyCat.pan %>% rownames_to_column(var='id')  %>%
  pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT') %>%
  filter(id %in% sco.94.concat.tree.rooted$tip.label)


metadata <- genome.tbl %>% filter(genome %in% sco.94.concat.tree.rooted$tip.label)# %>% data.frame()
#rownames(metadata) <- metadata$genome
metadata <- metadata %>% mutate(ID=genome) %>% select(-genome) %>% select(ID, everything())

p1 <- p %<+% metadata
p1<-p1 + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
             mapping=aes(fill=phylogroup),
             color = "white", offset = 0.04,size = 0.02,width=0.03)+theme(legend.position='none')+
  scale_fill_manual(values=phylogroup.colors)+ggnewscale::new_scale_fill()


#library(ggstance)
#install.packages('ggstance')
p.cazy.tree <- ggtree::facet_plot(p1, panel = 'CAZY', data = cct,
           geom = ggstance::geom_barh,
           mapping = aes(x = value, fill = cazyCAT),
           stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])
p.cazy.tree
ggsave('~/DATA/MarbGenomics/Graphs/p.cazy.tree.pdf',p.cazy.tree ,width=5,height = 10,units='in')


p1 <- p1 %<+% metadata

ggtree::facet_plot(p1, panel = 's', data = metadata %>% as_tibble(),
                  geom = geom_point,
                  mapping = aes( x= Completeness)) +
  theme_tree2()#+ #theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])

p2 <- ggtree::facet_plot(p1, panel = 'Dot', data = metadata,
                   geom = geom_point,
                   mapping = aes( x= GC)) +
  theme_tree2()

p2 + ggtree::xlim_expand(c(50, 70), 'Dot') +
  scale_x_continuous(name=c("Distance Dot Units")) +
  theme_bw()






selected.var <- c('Completeness',
                  'Contamination',
                  'Genome_size',
                  'ambiguous_bases',
                  'scaffolds',
                  'N50_scaffolds',
                  'Mean_contig_length_bp',
                  'Longest_scaffold_bp',
                  'GC',
                  'GC_std',
                  'Coding_density',
                  'Genome_size',
                  'predicted_genes')

p2 <- p1
for(var in selected.var){
  p2 <- p2 + geom_fruit(
    data=metadata,
    geom=geom_point,
    mapping = aes_string(
      y='ID',
      x=var#,
      #group=phylogroup,
      #fill=phylogroup,
    ),
    #size=2,
    #shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  )

}
p2


mapping_fun <- cazy.pan %>%
  rownames_to_column(var='genome') %>%
  filter(genome %in% sel.genome) %>%
  pivot_longer(cols=-genome, names_to='var', values_to='value')


p15 <- p + ggnewscale::new_scale_fill() +
  geom_fruit(data=mapping_fun, geom=geom_tile,
             mapping=aes(y=genome, x=var, alpha=value),
             color = "grey30", offset = 0.04,size = 3)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))

ggsave('~/DATA/MarbGenomics/Graphs/p.cazy.tree_crazy.pdf',p15 ,width=6,height = 10,units='in')





p + ggnewscale::new_scale_fill() +
  geom_fruit(data=cct, geom=geom_tile,
             mapping=aes(y=id, x=cazyCAT, alpha=value, fill=cazyCAT),
             color = "grey30", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))








+
  geom_fruit(data=dat3, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38,
             orientation="y",
             stat="identity",
  )







#===============
facet_plot(p1, panel="dot", data=cct, geom=geom_point, aes(x=value), color="firebrick")
facet_plot(p1, panel="dot", data=metadata, geom=geom_point, aes(x=Contamination), color="firebrick")

test <- cct %>% select(id) %>% unique() %>% left_join(metadata, by=c('id'='ID')) %>% select(id,Contamination)
test <- test %>% left_join(p$data, by=c('id'= 'label')) %>% mutate(label=id)
facet_plot(p1, panel="dot", data=test, geom=geom_point, aes(x=Contamination), color="firebrick")



facet_plot(p1, panel="Genome_size", data= CheckM3_annotated, geom= geom_point,
           aes(x= Genome_size), color='firebrick',stat='identity') + theme_tree2()

#--------------------------------------------------------------------------------------------
#	ANNOTATED TREE
#---------------------------------------------------------------------------------------------
library(ggstance)

CheckM3_annotated = test %>% data.frame()
rownames(CheckM3_annotated) = test$id
CheckM3_annotated  = CheckM3_annotated[sco.94.concat.tree.rooted$tip.label,]
p <- ggtree(sco.94.concat.tree.rooted,branch.length='none') %<+% CheckM3_annotated

facet_plot(p, panel="SNP", geom=geom_point, mapping=aes(x= GC),data= CheckM3_annotated, color="firebrick")  + theme_tree2()# %>%

#facet_plot("BAR", CheckM2_annotated, geom_segment, aes(x=0, xend=dummy_bar_value, y=y, yend=y)) + theme_tree2()


p <- ggtree(tree2)


p <- ggtree(sco.94.concat.tree.rooted,branch.length='none')
p1 <- p %<+% CheckM3_annotated[,c(1:10)]


p2 <- facet_plot(p1, panel="Genome_size", data= CheckM3_annotated, geom= geom_point,
                 aes(x= Genome_size), color='firebrick',stat='identity')

p2 = p2 + geom_hline(data=d,aes(yintercept=y))

p3  = facet_plot(p2, panel="predicted_genes", data= CheckM3_annotated, geom= geom_point,
                 aes(x= predicted_genes), color='firebrick',stat='identity') + theme_tree2()
d = data.frame(y=1:57,panel='predicted_genes')
p3 = p3 + geom_hline(data=d,aes(yintercept=y))

p4  = facet_plot(p3, panel="GCcontent", data= CheckM3_annotated, geom= geom_point,
                 aes(x= GC), color='firebrick',stat='identity') + theme_tree2()
d = data.frame(y=1:57,panel='GCcontent')
p4 = p4 + geom_hline(data=d,aes(yintercept=y))

p5  = facet_plot(p4, panel="Coding_density", data= CheckM3_annotated, geom= geom_point,
                 aes(x= Coding_density), color='firebrick',stat='identity') + theme_tree2()
d = data.frame(y=1:57,panel='Coding_density')
p5 = p5 + geom_hline(data=d,aes(yintercept=y))

p6  = facet_plot(p5, panel="Contamination", data= CheckM3_annotated, geom= geom_point,
                 aes(x= Contamination), color='firebrick',stat='identity') + theme_tree2()
d = data.frame(y=1:57,panel='Contamination')
p6 = p6 + geom_hline(data=d,aes(yintercept=y))

p7  = facet_plot(p6, panel="Completeness", data= CheckM3_annotated, geom= geom_point,
                 aes(x= Completeness), color='firebrick',stat='identity') + theme_tree2()
d = data.frame(y=1:57,panel='Completeness')
p7 = p7 + geom_hline(data=d,aes(yintercept=y))

p7 + theme(axis.text.x = element_text(angle = 90, hjust = 1))




#install.packages('ggstance')
ggtree::facet_plot(p, panel = 'Stacked Barplot', data = cct,
                   geom = ggstance::geom_barh,
                   mapping = aes(x = value, fill = cazyCAT),
                   stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['dead_reef']])
#install.packages('ggstance')
ggtree::facet_plot(p, panel = 'Stacked Barplot', data = cct,
                   geom = ggstance::geom_barh,
                   mapping = aes(x = value, fill = cazyCAT),
                   stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['namatjira_qual']][2:6])


ggtree::facet_plot(p, panel = 'Stacked Barplot', data = cct,
                   geom = ggstance::geom_barh,
                   mapping = aes(x = value, fill = cazyCAT),
                   stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=rev(ochRe::ochre_palettes[['dead_reef']][2:6]))


ggtree::facet_plot(p, panel = 'Stacked Barplot', data = cct,
                   geom = ggstance::geom_barh,
                   mapping = aes(x = value, fill = cazyCAT),
                   stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=rev(ochRe::ochre_palettes[['nolan_ned']]))







#=================================================================================#

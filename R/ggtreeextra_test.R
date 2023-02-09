library(ggplot2)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)


metadata = df.chkm
subtree = drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(metadata))])






p <- ggtree(subtree, layout="fan", open.angle=10, size=0.5) + geom_tiplab(align=TRUE,size=0,linetype='dashed',linesize=0.3)
p
p1 <- p %<+% metadata

p2 <- p1 #+
      #geom_point(shape=21, mapping=aes(fill= group, size= 2))


p3 = p2 + new_scale_fill() +
         geom_fruit(geom=geom_tile,
                    mapping=aes(fill=group),
                    color = "white", offset = 0.04,size = 0.02,width=0.03)+
         scale_alpha_continuous(range=c(0, 1),
                             guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))

p4 = p3	+ new_scale_fill() +
         geom_fruit(geom=geom_tile,
                    mapping=aes(fill=GC),
                    color = "white", offset = 0.06,size = 0.02,width=0.03)+
        scale_fill_gradientn(colors=c('forestgreen','white','purple'))

p5 = p4	+ new_scale_fill() +
         geom_fruit(geom=geom_tile,
                    mapping=aes(fill= Coding_density),
                    color = "white", offset = 0.06,size = 0.02,width=0.03)+
        scale_fill_gradientn(colors=c('blue','white','red'))

p6 = p5 +
      new_scale_fill() +
      geom_fruit(
          geom=geom_bar,
          mapping=aes(y= label, x= selgene, fill=group),
          pwidth=1, #0.3
          offset=0.05,
          stat="identity",
          orientation="y", # the orientation of axis.
          axis.params=list(
                          axis="x", # add axis text of the layer.
                          text.angle=-45, # the text size of axis.
                          hjust=0  # adust the horizontal position of text of axis.
                      ),
          grid.params=list() # add the grid line of the external bar plot.
      )



p7=  p6   +       geom_treescale(fontsize=2, linesize=0.3) +
         theme(legend.position=c(0.93, 0.5),
               legend.background=element_rect(fill=NA),
               legend.title=element_text(size=6.5),
               legend.text=element_text(size=4.5),
               legend.spacing.y = unit(0.02, "cm"),
             )

p7

#p7 + layout_rectangular()
p4 + layout_rectangular()













#-------------------- <- #


p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data=dat2,
    geom=geom_tile,
    mapping=aes(y=ID, x=Pos, fill=Type),
    offset=0.08,   # The distance between layers, default is 0.03 of x range of tree.
    pwidth=0.25, # width of the layer, default is 0.2 of x range of tree.
    axis.params=list(
      axis="x",  # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      hjust=0  # adust the horizontal position of text of axis.
    )
  ) +
  scale_fill_manual(
    values=c("#339933", "#dfac03"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  )





#=========================================================================
#----- COLOR THIS ONE BASED ON THE OUT.COLORS DATASET


plotdat <-genome.tbl %>%
  left_join(sorted.cliques,by=c('genome'='genome')) %>%
  left_join(out.colors, by=c('genome'='genome')) %>%
  mutate(phylogroup=phylogroup.x) %>%
  #filter(genome %in% sel.genome) %>%
  data.frame()
plotdat$group <- paste0('cl',plotdat[,colcol])
plotdat$group <- factor(plotdat$group,levels = c(paste0('cl',1:(length(unique(plotdat$group))-1)),'clNA'))


col.list <- list()
for(colcol in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  plotdat <-genome.tbl %>%
    left_join(sorted.cliques,by=c('genome'='genome')) %>%
    left_join(out.colors, by=c('genome'='genome')) %>%
    mutate(phylogroup=phylogroup.x) %>%
    filter(genome %in% sel.genome) %>%
    data.frame()
  plotdat$group <- paste0('cl',plotdat[,colcol])
  plotdat$group <- factor(plotdat$group,levels = c(paste0('cl',1:(length(unique(plotdat$group))-1)),'clNA'))

  colors <- plotdat%>% select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>%
    filter(group!='clNA') %>%
    unique() %>% data.frame()

  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,'group'])
  grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
  grid.colors <- c(grid.colors,'clNA'='white')
  col.list[[colcol]]<-grid.colors
}



mtd <- plotdat
rownames(mtd)<-plotdat$genome

#without outgroup
drop.outgroup = TRUE


#with outgroup
drop.outgroup=FALSE

if(drop.outgroup==TRUE){
  c.offset=0.07
  c.size=0.2 #if you want separating lines
  c.width=0.05
}else{
  c.offset=0.01
  c.size=0.4
  c.width=0.2
}

clique.colors <- grid.colors

circ <- ggtree(tr.tree, layout="fan", open.angle=10) +
  geom_tiplab2(align=TRUE,size=3,linetype='dashed',linesize=0.3)


#the inclusive annotated tree from phylogeny.R
#circ <- p.o #from phylogeny.r

#circ <-  tree_labeled

#
#0.03
c.offset=0.06
c.size=0.3 #if you want separating lines
c.width=0.07

c.offset=0.1
c.size=0.5 #if you want separating lines
c.width=0.07

c.offset2=0.05
c.size2=0.5 #if you want separating lines
c.width2=0.07



circ <- ggtree(tr.tree, layout="fan", open.angle=10) +
  geom_tiplab2(align=TRUE,size=2,hjust=-1,linetype='dashed',linesize=0.3)

circ$data <- circ$data %>% left_join(mtd, by=c('label'='genome'))




p1 <- circ
#p1 <-tree_labeled
p2 <- p1 + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=TypeStrain),
                          color = "white", shape=21,size=1.5)+
  scale_fill_manual(values=c('white','black'), na.value='white')+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=GC),
                          color = "white", offset = c.offset, size=c.size,width=c.width)+
   scale_fill_distiller(palette='Greens',direction=1)+
#  scale_fill_distiller(palette='Oranges',direction=1)+
  #scale_fill_viridis()+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=Coding_density),
                          color = "white", offset = c.offset, size=c.size,width=c.width)+
  scale_fill_distiller(palette='Purples',direction=1)+
  #scale_fill_viridis(option = 'cividis')+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=b),
                          color = "white", offset = c.offset, size=c.size,width=c.width)+
  scale_fill_distiller(palette='Blues',direction=1)+
  #scale_fill_viridis(option = 'magma')+
  #scale_fill_manual(values=phylogroup.colors)+
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = "white", offset = c.offset+0.05,size=c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=phylogroup.colors)+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.8)),
                          color = "white", offset = c.offset2,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[1]])+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.85)),
                          color = "white", offset = c.offset2,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[2]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.9)),
                          color = "white", offset = c.offset2,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[3]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.95)),
                          color = "white", offset = c.offset2,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=col.list[[4]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.98)),
                          color = "white", offset = c.offset2,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')

p2 <- p2 + geom_treescale()
#p2 + geom_tiplab(aes(angle=angle),size=2 ,align=TRUE,hjust=-1)


open_tree(p2+geom_tiplab2(),angle=90)
p2



ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/full_annotated_tree_test.pdf',plot=p2, width = 15, height = 10,unit='in')



ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_sortedCliques_coloured_ch.pdf',plot=p2, width = 15, height = 15,unit='in')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_sortedCliques_coloured_smaller_ch.pdf',plot=p2, width = 7, height = 7,unit='in')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_sortedCliques_coloured_smaller_open_ch.pdf',plot=open_tree(p2,angle=180), width = 8, height = 8,unit='in')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_sortedCliques_coloured_smaller_open90_ch.pdf',plot=open_tree(p2,angle=90), width = 8, height = 8,unit='in')

p.named <- p1 + geom_tiplab(align=T,size=2)
p.named



#=======
c.offset=0.05
c.size=0.5 #if you want separating lines
c.width=0.05



circ <- ggtree(tr.tree, layout="fan", open.angle=10)
circ$data <- circ$data %>% left_join(mtd, by=c('label'='genome'))


circ + geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dashed',linesize=0.3)


p1 <- circ

border.color <- 'lightgrey'

p1 +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=TypeStrain),
                          color = "white", shape=21,size=2)+
  scale_fill_manual(values=c('white','black'), na.value='white')+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=Coding_density),
                          color = "white", offset = c.offset, size=c.size,width=c.width)+
  scale_fill_distiller(palette='Purples',direction=1)+
  #scale_fill_viridis(option = 'cividis')+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=b),
                          color = "white", offset = c.offset, size=c.size,width=c.width)+
  scale_fill_distiller(palette='Blues',direction=1)+
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=GC),
                          color = "white", offset = c.offset, size=c.size,width=c.width)+
  scale_fill_distiller(palette='Greens',direction=1)+
  ggnewscale::new_scale_fill() +

  #---phylogroup
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = border.color, offset = 0.1,size=c.size,width=c.width)+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=phylogroup.colors)+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.8)),
                          color = border.color, offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[1]])+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.85)),
                          color = border.color, offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[2]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.9)),
                          color = border.color, offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[3]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.95)),
                          color = border.color, offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=col.list[[4]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.98)),
                          color = border.color, offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+
  ggnewscale::new_scale_fill()



#=======================
#### ATTEMPT FRESH ####

border.color <- 'lightgrey'

circ <- ggtree(tr.tree, layout="fan", open.angle=90)
circ$data <- circ$data %>% left_join(mtd, by=c('label'='genome'))
circ <- rotate_tree(circ,90)

#circ + geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dashed',linesize=0.3)
#circ

p2 <- circ +
  #typestrain
  #ggtreeExtra::geom_fruit(geom=geom_point,
  #                        mapping=aes(fill=TypeStrain),
  #                        color = "white", shape=21,size=2)+
  #scale_fill_manual(values=c('white','black'), na.value='white')+
  #ggnewscale::new_scale_fill()+
  ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                         mapping=aes(fill=TypeStrain),
                                                         color = border.color, offset = 0.1,size=c.size,width=0.1,
                                                         axis.params = c(axis='x',text='Type strain',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=c('grey','black'))+

  #----others
  ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                         mapping=aes(fill=quality_class),
                                                         color = border.color, offset = 0.1,size=c.size,width=0.1,
                                                         axis.params = c(axis='x',text='quality class',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=c('grey','black','red'))+

  #Habitat
  ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                         mapping=aes(fill=SS),
                                                         color = border.color, offset = 0.1,size=c.size,width=0.1,
                                                         axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=habitat.colors) +
  #----absolute lattidude
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=abs(lat)),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

  #phylogroup
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = border.color, offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=phylogroup.colors)+

  #ANI
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.8)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[1]])+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.85)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[2]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.9)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[3]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.95)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=col.list[[4]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.98)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+


  #Genome characteristics

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=b),
                          color = "white", offset = 0.15,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='piBias',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Blues',direction=1)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=Coding_density),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% Coding density',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Purples',direction=1)+
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=GC),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greens',direction=1)+



  #bars
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(data=k.df,
                          geom=geom_bar,
                          mapping=aes(y=genome,x=n,fill=phylogroup),
                          pwidth=0.35,orientation='y', stat='identity', offset=0.2,
                          axis.params = c(axis='x',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_manual(values=phylogroup.colors)





ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_sortedCliques_coloured_smaller_open90_cblablah.pdf',plot=p2, width = 12, height =12,unit='in')








#============= VERSION 1

p2 <- circ +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=TypeStrain),
                          color = "white", shape=21,size=2)+
  scale_fill_manual(values=c('white','black'), na.value='white')+
  ggnewscale::new_scale_fill()+

  #Habitat
  ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                         mapping=aes(fill=SS),
                                                         color = 'white', offset = 0.1,size=0.5,width=0.1,
                                                         axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=habitat.colors) +
  #----absolute latidude
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=abs(lat)),
                          color = "white", offset = 0.1,size=0.5,width=0.1,
                          axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

  #phylogroup
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = border.color, offset = 0.15,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=phylogroup.colors)+

  #ANI
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.8)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[1]])+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.85)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[2]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.9)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[3]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.95)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=col.list[[4]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.98)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+


  #Genome characteristics

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=b),
                          color = "white", offset = 0.15,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='piBias (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=Coding_density),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% Coding density (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=GC),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% GC-content (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  #bars
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(data=k.df,
                          geom=geom_bar,
                          mapping=aes(y=genome,x=n,fill=phylogroup),
                          pwidth=0.35,orientation='y', stat='identity', offset=0.1,
                          axis.params = c(axis='x',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_manual(values=phylogroup.colors) + geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dotted',linesize=0.1,color='lightgrey')


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/Figure_phylogeny_open_version1.pdf',plot=p2, width = 12, height =12,unit='in')


#============= VERSION 2

p2 <- circ +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=TypeStrain),
                          color = "white", shape=21,size=2)+
  scale_fill_manual(values=c('white','black'), na.value='white')+
  ggnewscale::new_scale_fill()+

  #Habitat
  ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                         mapping=aes(fill=SS),
                                                         color = 'white', offset = 0.1,size=0.5,width=0.1,
                                                         axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=habitat.colors) +
  #----absolute latidude
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=abs(lat)),
                          color = "white", offset = 0.1,size=0.5,width=0.1,
                          axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

  #Genome characteristics

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=b),
                          color = "white", offset = 0.15,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='piBias (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=Coding_density),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% Coding density (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=GC),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% GC-content (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+

  #phylogroup
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = border.color, offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=phylogroup.colors)+

  #ANI
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.8)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[1]])+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.85)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[2]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.9)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[3]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.95)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=col.list[[4]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.98)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+



  #bars
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(data=k.df,
                          geom=geom_bar,
                          mapping=aes(y=genome,x=n,fill=phylogroup),
                          pwidth=0.35,orientation='y', stat='identity', offset=0.1,
                          axis.params = c(axis='x',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_manual(values=phylogroup.colors) + geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dotted',linesize=0.1,color='lightgrey')


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/Figure_phylogeny_open_version2.pdf',plot=p2, width = 12, height =12,unit='in')


#============= VERSION 3

p2 <- circ +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=TypeStrain),
                          color = "white", shape=21,size=2)+
  scale_fill_manual(values=c('white','black'), na.value='white')+
  ggnewscale::new_scale_fill()+


  #Genome characteristics

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=b),
                          color = "white", offset = 0.15,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='piBias (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=Coding_density),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% Coding density (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=GC),
                          color = "white", offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='% GC-content (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='Greys',direction=1)+

  #phylogroup
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=phylogroup),
                          color = border.color, offset = 0.1,size=c.size,width=0.1,
                          axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=phylogroup.colors)+

  #ANI
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.8)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[1]])+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.85)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[2]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.9)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[3]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.95)),
                          color = border.color, offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=col.list[[4]])+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=paste0('cl',sc0.98)),
                          color = border.color,  offset = 0.075,size=0.1,width=0.075,
                          axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+

  #Habitat
  ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                         mapping=aes(fill=SS),
                                                         color = 'white', offset = 0.1,size=0.5,width=0.1,
                                                         axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_discrete(na.value = 'grey')+
  scale_fill_manual(values=habitat.colors) +
  #----absolute latidude
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=abs(lat)),
                          color = "white", offset = 0.1,size=0.5,width=0.1,
                          axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

  #bars
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(data=k.df,
                          geom=geom_bar,
                          mapping=aes(y=genome,x=n,fill=phylogroup),
                          pwidth=0.35,orientation='y', stat='identity', offset=0.1,
                          axis.params = c(axis='x',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
  scale_fill_manual(values=phylogroup.colors) + geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dotted',linesize=0.1,color='lightgrey')

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/Figure_phylogeny_open_version3.pdf',plot=p2, width = 12, height =12,unit='in')









c.offset=0.1
c.size=0 #if you want separating lines
c.width=0.07
spacer <- 0.05

# more of the tree
#c.offset=0.1 * 0.75
#c.size=0 #if you want separating lines
#c.width=0.07 * 0.75
#spacer <- 0.05 * 0.75



{

  #============= VERSION 1

  p2 <- circ +
    ggtreeExtra::geom_fruit(geom=geom_point,
                            mapping=aes(fill=TypeStrain),
                            color = "white", shape=21,size=2)+
    scale_fill_manual(values=c('white','black'), na.value='white')+
    ggnewscale::new_scale_fill()+

    #Habitat
    ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                           mapping=aes(fill=SS),
                                                           color = 'white', offset = c.offset,size=0.5,width=c.width,
                                                           axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'grey')+
    scale_fill_manual(values=habitat.colors) +
    #----absolute latidude
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=abs(lat)),
                            color = "white", offset = c.offset ,size=0.5,width=c.width,
                            axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

    #phylogroup
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=phylogroup),
                            color = border.color, offset = c.offset +spacer ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'grey')+
    scale_fill_manual(values=phylogroup.colors)+

    #ANI
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.8)),
                            color = border.color, offset= c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[1]])+
    ggnewscale::new_scale_fill() +

    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.85)),
                            color = border.color,  offset= c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[2]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.9)),
                            color = border.color, offset=   c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[3]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.95)),
                            color = border.color, offset = c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
    scale_fill_manual(values=col.list[[4]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.98)),
                            color = border.color,  offset = c.offset, size=c.size,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+


    #Genome characteristics

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=b),
                            color = "white", offset = c.offset +spacer ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='piBias (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=Coding_density),
                            color = "white",  offset = c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='% Coding density (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill()  +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=GC),
                            color = "white",  offset = c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='% GC-content (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    #bars
    geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dotted',linesize=0.1,color='lightgrey')


  ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/Figure_phylogeny_open_version1_nobar.pdf',plot=p2, width = 12, height =12,unit='in')


  #============= VERSION 2


  p2 <- circ +
    ggtreeExtra::geom_fruit(geom=geom_point,
                            mapping=aes(fill=TypeStrain),
                            color = "white", shape=21,size=2)+
    scale_fill_manual(values=c('white','black'), na.value='white')+
    ggnewscale::new_scale_fill()+

    #Habitat
    ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                           mapping=aes(fill=SS),
                                                           color = 'white', offset = c.offset,size=0.5,width=c.width,
                                                           axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'grey')+
    scale_fill_manual(values=habitat.colors) +
    #----absolute latidude
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=abs(lat)),
                            color = "white", offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

    #Genome characteristics

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=b),
                            color = "white", offset = c.offset + spacer,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='piBias (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=Coding_density),
                            color = "white", offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='% Coding density (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill()  +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=GC),
                            color = "white", offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='% GC-content (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+

    #phylogroup
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=phylogroup),
                            color = border.color, offset = c.offset + spacer ,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'grey')+
    scale_fill_manual(values=phylogroup.colors)+

    #ANI
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.8)),
                            color = border.color, offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[1]])+
    ggnewscale::new_scale_fill() +

    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.85)),
                            color = border.color,  offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[2]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.9)),
                            color = border.color,  offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[3]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.95)),
                            color = border.color, offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
    scale_fill_manual(values=col.list[[4]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.98)),
                            color = border.color,  offset = c.offset,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+



    #bars
   geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dotted',linesize=0.1,color='lightgrey')


  ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/Figure_phylogeny_open_version2_nobar.pdf',plot=p2, width = 12, height =12,unit='in')


  #============= VERSION 3



  border.color <- 'lightgrey'

  circ <- ggtree(tr.tree, layout="fan", open.angle=90,size=0.3)
  circ$data <- circ$data %>% left_join(mtd, by=c('label'='genome'))
  circ <- rotate_tree(circ,90)



  p2 <- circ +
    ggtreeExtra::geom_fruit(geom=geom_point,
                            mapping=aes(fill=TypeStrain),
                            shape=21,size=1.5,color='white')+
    scale_fill_manual(values=c(NA,'black'),na.value=NA)+
    ggnewscale::new_scale_fill()+

    #ggtreeExtra::geom_fruit(geom=geom_star,
    #                        mapping=aes(fill=TypeStrain),
    #                        starshape=26,
    #                        color = NA,
    #                        starstroke=0,
    #                        size=3,
    #                        offset=spacer*2)+
    #scale_fill_manual(values=c('white','black'), na.value=NA)+
    #ggnewscale::new_scale_fill()+


    #Genome characteristics

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=b),
                            color = "white", offset =  c.offset, size=c.size,width=c.width,
                            axis.params = c(axis='x',text='piBias (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=Coding_density),
                            color = "white",  offset =  c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='% Coding density (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill()  +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=GC),
                            color = "white",  offset =  c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='% GC-content (123)',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+
    ggnewscale::new_scale_fill()  +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=Genome_size),
                            color = "white",  offset =  c.offset ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='Genome size',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='Greys',direction=1)+

    #phylogroup
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=phylogroup),
                            color = border.color,  offset =  c.offset + spacer ,size=c.size,width=c.width,
                            axis.params = c(axis='x',text='Phylogroup',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'grey')+
    scale_fill_manual(values=phylogroup.colors)+

    #ANI
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.8)),
                            color = border.color,  offset =  c.offset ,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.8',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[1]])+
    ggnewscale::new_scale_fill() +

    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.85)),
                            color = border.color,   offset =  c.offset ,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.85',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[2]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.9)),
                            color = border.color,  offset =  c.offset ,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.9',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[3]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.95)),
                            color = border.color,  offset =  c.offset ,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.95',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
    scale_fill_manual(values=col.list[[4]])+

    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=paste0('cl',sc0.98)),
                            color = border.color,  offset =  c.offset ,size=0.1,width=c.width,
                            axis.params = c(axis='x',text='ANI>0.98',text.angle=0, hjust=0,text.size=3,fontface='bold'))+

    scale_fill_discrete(na.value = 'white')+
    scale_fill_manual(values=col.list[[5]]) + theme(legend.position='left')+

    #Habitat
    ggnewscale::new_scale_fill()+  ggtreeExtra::geom_fruit(geom=geom_tile,
                                                           mapping=aes(fill=SS),
                                                           color = 'white',  offset =  c.offset +spacer ,size=0.5,width=c.width,
                                                           axis.params = c(axis='x',text='Environment',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_discrete(na.value = 'grey')+
    scale_fill_manual(values=habitat.colors) +
    #----absolute latidude
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_tile,
                            mapping=aes(fill=abs(lat)),
                            color = "white", offset =  c.offset ,size=0.5,width=c.width,
                            axis.params = c(axis='x',text='Absolute latitude',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_distiller(palette='RdBu',direction=1,na.value = 'lightgrey') +

    #host-associated?
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(geom=geom_star,
                            mapping=aes(fill=ifelse(host=='',NA,TRUE)),
                            starshape=26,
                            color = NA,
                            starstroke=0,
                            offset =  c.offset + spacer + spacer ,
                            size=3.5,
                            axis.params = c(axis='x',text='host-associated',text.angle=0, hjust=0,text.size=3,fontface='bold'))+
    scale_fill_manual(values=c('black',NA),na.value=NA)+

    #bars
    geom_tiplab2(aes(label=strain),align=TRUE,size=0,hjust=-1,linetype='dotted',linesize=0.1,color='lightgrey')


  ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/Figure_phylogeny_open_version3_nobar_impsrosvalsw.pdf',plot=p2, width = 11, height =11,unit='in')



}




















Genome_size


circ +  #bars
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(geom=geom_col,
                          mapping=aes(x=label,y=Genome_size,fill=phylogroup),
                          pwidth=0.35,orientation='y',position=ggtreeExtra::position_stackx(),
                          axis.params = c(axis='y',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))#+






k.df <- annotation.tbl %>% filter(grepl('transpos',product,ignore.case = TRUE)) %>% group_by(genome) %>% tally()

circ +  #bars
  ggnewscale::new_scale_fill()  +
  ggtreeExtra::geom_fruit(data=k.df,
                          geom=geom_bar,
                          mapping=aes(y=genome,x=n,fill=phylogroup),
                          pwidth=0.35,orientation='y', stat='identity',
                          axis.params = c(axis='x',text='%GC',text.angle=0, hjust=0,text.size=3,fontface='bold'))













#+


  #--clades

  geom_cladelab(data=cladedat,mapping=aes(node=nodeid,label=label,color=label),
                align=T,
                fontsize=3,
                barsize=1,offset.text=0.1,offset=0.4) +
  scale_color_manual(values=phylogroup.colors)



cladedat <- data.frame(nodeid=phylogroup.nodes, label=paste0('P',1:11))




#=================

p2 + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_bar, mapping=aes(
    x=GC, y=genome#,# fill=phylogroup
  ), stat="identity",
  orientation="y", offset=0.2, pwidth=0.5, axis.params=list(
    axis = "x", text.angle = -45, hjust = 0,
    vjust = 0.5, nbreak = 4
  )
  )+ scale_fill_manual(values=phylogroup.colors)





p.cazy.tree



ggtree::facet_plot(p1, panel = 'CAZY', data = cct,
                   geom = ggstance::geom_barh,
                   mapping = aes(x = value, fill = cazyCAT),
                   stat='identity' ) + theme_tree2()+fdb_style(aspect.ratio = 3)+scale_fill_manual(values=ochRe::ochre_palettes[['tasmania']])
p.cazy.tree


#----
#with geom point

p1 + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=as.factor(c0.75)),
                          color = "white", shape=21,size=2)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=as.factor(c0.8)),
                          color = "white", shape=21,size=2)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=as.factor(c0.85)),
                          color = "white", shape=21,size=2)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=as.factor(c0.9)),
                          color = "white", shape=21,size=2)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=as.factor(c0.95)),
                          color = "white", shape=21,size=2)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_point,
                          mapping=aes(fill=as.factor(c0.98)),
                          color = "white", shape=21,size=2)




#----
# with geom_tile

if(drop.outgroup==TRUE){
  c.offset=0.07
  #c.size=0 #turn off line width
  #add line width, spacing the blocks
  c.size=0.2 #if you want separating lines
  c.width=0.05
}else{
  c.offset=0.01
  c.size=0.4
  c.width=0.2
}
#p2 <- p1 + ggnewscale::new_scale_fill() +
p1 + #ggnewscale::new_scale_fill() +
  #ggtreeExtra::geom_fruit(geom=geom_tile,
  #           mapping=aes(fill=as.factor(c0.7)),
  #           color = "white", size=c.size,width=c.width)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=as.factor(c0.75)),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=as.factor(c0.8)),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=as.factor(c0.85)),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=as.factor(c0.9)),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  #scale_fill_gradientn(colors=c('forestgreen','white','purple'))
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=as.factor(c0.95)),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  #scale_fill_gradientn(colors=c('forestgreen','white','purple'))
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=as.factor(c0.98)),
                          color = "white", offset = c.offset,size = c.size,width=c.width)
#scale_fill_gradientn(colors=c('forestgreen','white','purple'))
p2

p2  + ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=c85),
    color = "white", offset = 0.04,size = 0.02,width=0.03)+scale_fill_gradientn(colors=c('black','red'))+
  #ggnewscale::new_scale_fill() +
  #ggtreeExtra::geom_fruit(geom=geom_tile,
  #           mapping=aes(fill=c0.85),
  #           color = "white", offset = 0.04,size = 0.02,width=0.03)+
  #scale_fill_gradientn(colors=c('forestgreen','white','purple'))
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=c90),
    color = "white", offset = 0.04,size = 0.02,width=0.03)+scale_fill_gradientn(colors=c('black','red'))+
  #scale_fill_gradientn(colors=c('forestgreen','white','purple'))
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=c95),
    color = "white", offset = 0.04,size = 0.02,width=0.03)+scale_fill_gradientn(colors=c('black','red'))+
  #scale_fill_gradientn(colors=c('forestgreen','white','purple'))
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=c98),
    color = "white", offset = 0.04,size = 0.02,width=0.03)+scale_fill_gradientn(colors=c('black','red'))





#====
#REVERSE, LIKE A ANVIO PLOT


pan = pan_clean(OG.pan[tip.order, ])
#presabs
pan[pan>1]<-1


rev.dist <- vegan::vegdist(t(as.matrix(pan)),method='jaccard')
rev.dist[is.na(rev.dist)] <- 1

rev.clust <- hclust(rev.dist)

p.r <- ggtree(rev.clust,layout='circular')
p.heat <- gheatmap(p.r, t(pan))
p.opened <- open_tree(p.heat,angle=90)

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/replicate2_anvio_open90_ch.pdf',plot=p.opened, width = 8, height = 8,unit='in')




#, offset=5, width=0.5, font.size=3, colnames_angle=-45, hjust=0)


+
  scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"),
                    values=c("steelblue", "firebrick", "darkgreen"), name="genotype")



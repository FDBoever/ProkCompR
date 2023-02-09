#ASTRAL

#tree
concat.rps <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins/concat.treefile', clean_genome_names=TRUE)
concat.rps.nuc <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins_pal2nal/concat.treefile', clean_genome_names=TRUE)
concat.sco <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile', clean_genome_names=TRUE)
concat.sco.nuc <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
concat.rDNA <- loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees//concat.treefile', clean_genome_names=TRUE)
rDNA16S <- loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees/16S.treefile', clean_genome_names=TRUE)

# SCO
#java -jar astral.5.7.5.jar -i /Users/sa01fd/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/loci.treefile -o /Users/sa01fd/DATA/MarbGenomics/all_SCO_genafpear_trimai_astral.tree 2>/Users/sa01fd/DATA/MarbGenomics/all_SCO_genafpear_trimai_astral.tree.log
#java -jar astral.5.7.5.jar -i /Users/sa01fd/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/loci.treefile -o /Users/sa01fd/DATA/MarbGenomics/all_SCO_genafpear_trimai_astral.tree 2>/Users/sa01fd/DATA/MarbGenomics/all_SCO_genafpear_trimai_astral.tree.log
#java -jar astral.5.7.5.jar -i ~/DATA/MarbGenomics/RibSCO_Protein_alignments_mafft_genafpair_pal2nal_trimAl/loci.treefile -o /Users/sa01fd/DATA/MarbGenomics/RibSco_MB_SCO_genafpear_trimai_astral_pal2nal.tree 2>/Users/sa01fd/DATA/MarbGenomics/RibSco_MB_SCO_genafpear_trimai_astral_pal2nal.tree.log

#Idialy find a solution for this, as this really is odd...
concat.rps$tip.label[concat.rps$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.rps.nuc$tip.label[concat.rps.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.sco$tip.label[concat.sco$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.sco.nuc$tip.label[concat.sco.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.rDNA$tip.label[concat.rDNA$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
rDNA16S$tip.label[rDNA16S$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"


#========================================================================
#   Phylogenetic of the Marinobacter paper
#=======================================================================

#on GTDB, following genomes are cassified as "Marinobacter A"
Marinobacter_A <- c("Marinobacter_sp_UBA2698",
                    "Marinobacter_sp_YJ_S3_2",
                    "Marinobacter_sp_R17",
                    "Marinobacter_bohaiensis_T17",
                    "Marinobacter_nanhaiticus_D15_8W",
                    "Marinobacter_fonticola_CS412")


#Read 106 genome trees (containing outgroup genera)
#----------------------------------------------------------------#
#based on amino acid alignments
rib.106.concat.tree = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins/concat.treefile', clean_genome_names=TRUE)
sco.106.concat.tree = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile', clean_genome_names=TRUE)
astral.106.tree = loadTree(inPath='/Users/sa01fd/DATA/MarbGenomics/all_SCO_genafpear_trimai_astral.tree', clean_genome_names=TRUE)
astral.106.tree$edge.length[is.na(astral.106.tree$edge.length)] <- 0

#based on nucleotide alignments
rib.106.concat.nuc.tree <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins_pal2nal/concat.treefile', clean_genome_names=TRUE)
sco.106.concat.nuc.tree <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
astral.106.nuc.tree = loadTree(inPath='/Users/sa01fd/DATA/MarbGenomics/SCO_filtered_mafft_genafpair_pal2nal_trimAl_ASTRAL.tree', clean_genome_names=TRUE)
astral.106.nuc.tree$edge.length[is.na(astral.106.nuc.tree$edge.length)] <- 0

#concat.rps$tip.label[concat.rps$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
rib.106.concat.tree$tip.label[rib.106.concat.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
sco.106.concat.tree$tip.label[sco.106.concat.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
astral.106.tree$tip.label[astral.106.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
rib.106.concat.nuc.tree$tip.label[rib.106.concat.nuc.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
sco.106.concat.nuc.tree$tip.label[sco.106.concat.nuc.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
astral.106.nuc.tree$tip.label[astral.106.nuc.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"

#concatenated rDNA operon
concat.rDNA <- loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees//concat.treefile', clean_genome_names=TRUE)

#16S rDNA gene
rDNA16S <- loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees/16S.treefile', clean_genome_names=TRUE)



# Read 94 genome species trees (containing ingroup only)
#----------------------------------------------------------------#
rib.94.concat.tree = loadTree(inPath='~/DATA/MarbGenomics/RibMB_SCO_Protein_alignments_mafft_genafpair/concat.treefile', clean_genome_names=TRUE)
sco.94.concat.tree = loadTree(inPath='~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair_trimAl/concat.treefile', clean_genome_names=TRUE)
astral.94.tree = loadTree(inPath='~/DATA/MarbGenomics/ASTRAL.94.MBSCO.trimal.tree', clean_genome_names=TRUE)
astral.94.tree$edge.length[is.na(astral.94.tree$edge.length)] <- 0
astral.94.tree$edge.length <- astral.94.tree$edge.length + .0001

rib.94.concat.nuc.tree = loadTree(inPath='~/DATA/MarbGenomics/RibSCO_Protein_alignments_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
sco.94.concat.nuc.tree = loadTree(inPath='~/DATA/MarbGenomics/MB_con_SCO_Protein_alignments_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
astral.94.nuc.tree= loadTree(inPath='~/DATA/MarbGenomics/con_MB_SCO_genafpear_trimai_astral_pal2nal.tree', clean_genome_names=TRUE)
astral.94.nuc.tree$edge.length[is.na(astral.94.nuc.tree$edge.length)] <- 0
astral.94.nuc.tree$edge.length <- astral.94.nuc.tree$edge.length + .0001


rib.94.concat.tree$tip.label[rib.94.concat.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
sco.94.concat.tree$tip.label[sco.94.concat.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
astral.94.tree$tip.label[astral.94.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
rib.94.concat.nuc.tree$tip.label[rib.94.concat.nuc.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
sco.94.concat.nuc.tree$tip.label[sco.94.concat.nuc.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
astral.94.nuc.tree$tip.label[astral.94.nuc.tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"





#root trees
#----------------------------------------------------------------#
rib.94.concat.tree.rooted <- ape::root(rib.94.concat.tree, outgroup = Marinobacter_A)
sco.94.concat.tree.rooted <- ape::root(sco.94.concat.tree, outgroup = Marinobacter_A)
astral.94.tree.rooted <- ape::root(astral.94.tree, outgroup = Marinobacter_A)
rib.94.concat.nuc.tree.rooted <- ape::root(rib.94.concat.nuc.tree, outgroup = Marinobacter_A)
sco.94.concat.nuc.tree.rooted <- ape::root(sco.94.concat.nuc.tree, outgroup = Marinobacter_A)
astral.94.nuc.tree.rooted <- ape::root(astral.94.nuc.tree, outgroup = Marinobacter_A)


#-- read loci trees (leave unrooted?)
rib.94.loci.trees <- ape::read.tree('~/DATA/MarbGenomics/RibMB_SCO_Protein_alignments_mafft_genafpair/loci.treefile')
sco.94.loci.trees <- ape::read.tree('~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair_trimAl/loci.treefile')

rib.94.nuc.loci.trees <- ape::read.tree('~/DATA/MarbGenomics/RibSCO_Protein_alignments_mafft_genafpair_pal2nal_trimAl/loci.treefile')
sco.94.nuc.loci.trees <- ape::read.tree('~/DATA/MarbGenomics/MB_con_SCO_Protein_alignments_mafft_genafpair_pal2nal_trimAl/loci.treefile')

sco.106.loci.trees <- ape::read.tree('~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_pal2nal_trimAl/loci.treefile')



# Change pelagius !!!!!!!!!!




#sco.94.loci.trees.rooted <- multi.root(sco.94.loci.trees,Marinobacter_A)
#from concord.log, extract a table... did this manually, but should be easy to load automaticlly

trees.meta <- read.delim('~/DATA/MarbGenomics/concord_all_trees.txt')
trees.meta <- trees.meta %>% mutate(name=gsub('\\.faa','',Name))

sum(trees.meta$Sites)
sum(trees.meta$Infor)
sum(trees.meta$Infor)/sum(trees.meta$Sites)


#list of trees
#----------------------------------------------------------------#

#for 94 genome trees
lsTrees.94 =
  list('sco.94.concat'=sco.94.concat.tree,
       'rib.94.concat'=rib.94.concat.tree,
       'astral.94.tree'=astral.94.tree,
       'sco.94.nuc.concat'=sco.94.concat.nuc.tree,
       'rib.94.nuc.concat'=rib.94.concat.nuc.tree,
       'astral.94.nuc.tree'=astral.94.nuc.tree
  )
class(lsTrees.94) = 'multiPhylo'
lsTrees.94.rooted <- multi.root(lsTrees.94, Marinobacter_A)

#for 106 genome trees
lsTrees.106 =
  list('sco.106.concat'=sco.106.concat.tree,
       'rib.106.concat'=rib.106.concat.tree,
       'astral.106.tree'=astral.106.tree,
       'sco.106.nuc.concat'=sco.106.concat.nuc.tree,
       'rib.106.nuc.concat'=rib.106.concat.nuc.tree,
       'astral.106.nuc'=astral.106.nuc.tree
  )
class(lsTrees.106) = 'multiPhylo'
lsTrees.106.rooted <- multi.root(lsTrees.106, outgroup)


#Combine both datasets
lsTrees.rooted <- c(lsTrees.94.rooted, lsTrees.106.rooted)
lsTree.rooted <- lsTrees.94.rooted

#Assess monophyely in all
mono.tdtree <- robustMonophyly(lsTrees.rooted)

ggtree::ggtree(mono.tdtree, branch.length="none") + geom_tiplab(size=2) +
  geom_point(data=.%>%dplyr::filter(isTip==FALSE),aes(color=support)) +
  geom_label(data=. %>% filter(isTip==FALSE),aes(label=node,color=support),size=3)


ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=.%>%dplyr::filter(isTip==FALSE)%>%filter(support>0.4),aes(color=support,size=1.5))
ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=.%>%dplyr::filter(isTip==FALSE)%>%filter(support>0.6),aes(color=support,size=1.5))
ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=.%>%dplyr::filter(isTip==FALSE)%>%filter(support==1),aes(color=support,size=1.5))
ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=.%>%dplyr::filter(isTip==FALSE),aes(color=support))
ggtree(mono.tdtree) + geom_hilight(mapping=aes(subset = node %in% c(147,110, 117)))
ggtree(mono.tdtree) + geom_hilight(mapping=aes(subset = node %in% c(105,112,172,182,115,170,146,169,167,162,163)))
ggtree(mono.tdtree) + geom_hilight(mapping=aes(subset = node %in% c(105,112,172,182,115,146,169,167,162,163)))



#================================================================#
# visualise unrooted trees as part of Figure 1
#----------------------------------------------------------------#

plot(sco.94.concat.tree,
     type = "unrooted",
     no.margin = TRUE
)
tiplabels(pch=21, col='black', bg='white',cex=1)


selected.tree <- sco.94.concat.tree
selected.tree <-lsTrees.94[["astral.94.tree"]]

tmp <- genome.tbl %>% select(genome, phylogroup) %>% filter(genome %in% selected.tree$tip.label) %>% data.frame()
colors.p <- phylogroup.colors[tmp$phylogroup]
names(colors.p) <- tmp$genome
colors.p <- colors.p[selected.tree$tip.label]

ape::plot.phylo(selected.tree,
     type = "unrooted",
     no.margin = TRUE,show.tip.label = FALSE
)
tiplabels(pch=21, col='black', bg=colors.p,cex=1)



grDevices::pdf(paste0(outdir,'/unrooted_trees_94','.pdf'), height = 5, width =10)
par(mfrow=c(1,3))
for(i in c("sco.94.concat","rib.94.concat" ,"astral.94.tree")){
  selected.tree <-lsTrees.94[[i]]

  tmp <- genome.tbl %>% select(genome, phylogroup) %>% filter(genome %in% selected.tree$tip.label) %>% data.frame()
  colors.p <- phylogroup.colors[tmp$phylogroup]
  names(colors.p) <- tmp$genome
  colors.p <- colors.p[selected.tree$tip.label]

  ape::plot.phylo(selected.tree,
                  type = "unrooted",
                  no.margin = TRUE,show.tip.label = FALSE
  )
  #tiplabels(pch=21, col='black', bg=colors.p,cex=1)
  tiplabels(pch=16, col=colors.p,cex=2)
  add.scale.bar(cex=0.7, font=2,col='black')

}
dev.off()





#=================================================================
#MANUAL INSPECTION RESULTED IN THESE NODES
phylogroup.nodes <- c(109,112,172,182,115,170,146,169,167,162,163)
#phylogroup.nodes <- c(105,112,172,182,115,170,146,169,167,162,163)

phylogroup.names <- paste0('P',seq(1:length(phylogroup.nodes)))
d <- data.frame(node=phylogroup.nodes, type=phylogroup.names)

#Assigning colors is importand here
phylogroup.colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(phylogroup.nodes))
names(phylogroup.colors) <- phylogroup.names
phylogroup.colors[3]<-"#87ceeb" #skyblue


#Assign base tree, on which the searches have been done
base.tree <- sco.94.concat.tree.rooted
p.base <- ggtree(base.tree)

# STORE OFFSPRING GENOMES PER NODE
phylogroup.offspring <- NULL
for(j in 1:dim(d)[1]){
  n <- d[j,'node']
  t <- as.character(d[j,'type'])

  offspring <- p.base$data %>% tidytree::offspring(n) %>%
    dplyr::filter(isTip==TRUE) %>% select(label) %>% pull

  phylogroup.offspring <- rbind(phylogroup.offspring,cbind('phylogroup'=t,'genome'=offspring))
  }
phylogroup.offspring <- phylogroup.offspring %>% data.frame()

#==================


#===========================================================
#
# annotate genome table!
# genome.tbl <- genome.tbl %>% select(-phylogroup.x, -phylogroup.y)
# genome.tbl <- genome.tbl %>% select(-phylogroup)
genome.tbl <- genome.tbl %>% left_join(phylogroup.offspring, by='genome') %>%
  mutate(phylogroup = factor(phylogroup,levels=as.character(unique(phylogroup.offspring$phylogroup))))

genome.tbl %>% select(genome,phylogroup)


#===========================================================
# ANNOTATE ALL TREES USING THESE PHYLOGROUPS
# DETECT THE ANCESTRAL NODE, AND COLOR THOSE

#example
tree_plot_ann(tr=lsTrees.106.rooted[[1]], phylogroup='phylogroup',layout='rectangular',nrow=1,extensive=TRUE, tipsize=2.5,nodes=TRUE)


tree_plot_ann <- function(tr,phylogroup='phylogroup',layout='circular',branch.length="",nrow=1, extensive=FALSE, tipsize=3, nodes=FALSE,highlight=TRUE,tr.name='tree',extendto=NULL,clade_label_align=NULL){
  phygrps <- genome.tbl %>% select(phylogroup) %>% unique() %>% pull %>% as.character() %>% sort()
  if(branch.length == 'none'){
    p.tmp <- ggtree(tr,branch.length='none',layout=layout)
  }else{
    p.tmp <- ggtree(tr,layout=layout)
  }

  d_t <- NULL
  for(phgr in phygrps){
    tps <- genome.tbl %>% filter(phylogroup==phgr) %>% select(genome) %>% pull %>% as.character()
    if(is.monophyletic(tr, tps)){

    }else{
      message('WARNING!! ',phgr, 'is not monophyletic in ')
    }
    nd = getMRCA( tr, tps)
    d_t <- rbind(d_t,c('node'=nd,'type'=phgr))
  }
  d_t <- data.frame(d_t)
  d_t$node <- as.numeric(as.character(d_t$node))
  d_t$type <- as.character(d_t$type)

  if(highlight==TRUE){
    if(is.null(extendto)){
      p <- p.tmp +
        geom_hilight(data=d_t, aes(node=node, fill=type),alpha=0.5)+
        scale_fill_manual(values =phylogroup.colors)
    }else{
      p <- p.tmp +
        geom_hilight(data=d_t, aes(node=node, fill=type),alpha=0.5,extendto=extendto)+
        scale_fill_manual(values =phylogroup.colors)
    }

  }else{
    p <- p.tmp
  }


  for(j in 1:dim(d_t)[1]){
    if(!is.null(clade_label_align)){
      p <- p + geom_cladelabel(node=d_t[j,'node'],
                               label=as.character(d_t[j,'type']),
                               barsize=1.5,
                               color=phylogroup.colors[as.character(d_t[j,'type'])])
    }else{
      p <- p + geom_cladelabel(node=d_t[j,'node'],
                               label=as.character(d_t[j,'type']),
                               barsize=1.5,
                               color=phylogroup.colors[as.character(d_t[j,'type'])],align=TRUE)
    }

  }

  if(extensive==TRUE){
    p <- p +
      geom_tiplab(size=tipsize) +
      geom_treescale()
  }

  if(nodes==TRUE){
    maxsupport <- p$data %>%filter(isTip==FALSE) %>% select(label) %>% filter(label!="") %>% pull() %>% as.numeric() %>% max()
    if(maxsupport>1){
      maxsupport=100
    }else{
      maxsupport=1
    }

    sf <- p$data %>%filter(isTip==FALSE)%>% filter(label!=maxsupport) %>% mutate(label=as.numeric(label))
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_point(data=sf,
                 aes(fill=label),
                 shape=21)+
      geom_text(data=sf,
                aes(label=label),
                size=tipsize, hjust=1.2,vjust=-0.5)+
      scale_fill_distiller(palette='Greys',direction = 1) +
      geom_treescale()
  }
  p <- p + ggtitle(tr.name) + theme(legend.position = "none")
  return(p)
}

#' Title
#'
#' @param lsTrees
#' @param phylogroup
#' @param layout
#' @param branch.length
#' @param nrow
#' @param extensive
#' @param tipsize
#' @param nodes
#'
#' @return
#' @export
#'
#' @examples
multitree_plot_ann <- function(lsTrees,phylogroup='phylogroup',layout='circular',branch.length="",nrow=1, extensive=FALSE, tipsize=3, nodes=FALSE,highlight=TRUE){
  p.list <- NULL
  for(tr.name in names(lsTrees)){
    p.list[[tr.name]] <- tree_plot_ann(lsTrees[[tr.name]],
                                            phylogroup=phylogroup,
                                            layout=layout,
                                            branch.length=branch.length,
                                            nrow=nrow,
                                            extensive=extensive,
                                            tipsize=tipsize,
                                            nodes=nodes,
                                            highlight=highlight,
                                            tr.name=tr.name)

  }

  p.comb <- grid.arrange(grobs=p.list,nrow=nrow)
  return(p.comb)
}

#==== 94 trees
p.comb_circular <- multitree_plot_ann(lsTrees=lsTrees.94.rooted, phylogroup='phylogroup',layout='circular',nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_speciestrees_circular.pdf',plot=p.comb_circular, width = 10, height = 8,unit='in')

p.comb_rectangular <- multitree_plot_ann(lsTrees=lsTrees.94.rooted, phylogroup='phylogroup',layout='rectangular',nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_speciestrees_rect.pdf',plot=p.comb_rectangular, width = 7, height = 10,unit='in')

p.comb_circular <- multitree_plot_ann(lsTrees=lsTrees.94.rooted, phylogroup='phylogroup',layout='circular',branch.length="none",nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_speciestrees_circular_no_brlengths.pdf',plot=p.comb_circular, width = 10, height = 10,unit='in')

p.comb_rectangular <- multitree_plot_ann(lsTrees=lsTrees.94.rooted, phylogroup='phylogroup',layout='rectangular',branch.length="none",nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_speciestrees_rect_no_brlengths.pdf',plot=p.comb_rectangular, width = 7, height = 10,unit='in')

#==== 106 trees
p.comb_106_rectangular <- multitree_plot_ann(lsTrees=lsTrees.106.rooted, phylogroup='phylogroup',layout='rectangular',nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_106_speciestrees_rect_no_brlengths.pdf',plot=p.comb_106_rectangular, width = 7, height = 8,unit='in')

p.comb_106_circular <- multitree_plot_ann(lsTrees=lsTrees.106.rooted, phylogroup='phylogroup',layout='circular',nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_106_speciestrees_circular.pdf',plot=p.comb_106_circular, width = 10, height = 10,unit='in')

p.comb_106_circular <- multitree_plot_ann(lsTrees=lsTrees.106.rooted, phylogroup='phylogroup',layout='circular',nrow=2,branch.length="none")
ggsave('~/DATA/MarbGenomics/Graphs/annotated_106_speciestrees_circular_no_brlengths.pdf',plot=p.comb_106_circular, width = 10, height = 10,unit='in')


#== using extensive=TRUE
p.comb_106_ext <- multitree_plot_ann(lsTrees=lsTrees.106.rooted, phylogroup='phylogroup',layout='rectangular',nrow=2,extensive=TRUE, tipsize=2,nodes=TRUE)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_106_speciestrees_circular_ext.pdf',plot=p.comb_106_ext, width = 13, height = 18,unit='in')

p.comb_106_ext <- multitree_plot_ann(lsTrees=lsTrees.106.rooted[1], phylogroup='phylogroup',layout='rectangular',nrow=1,extensive=TRUE, tipsize=2.5,nodes=TRUE)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_106_speciestrees_circular_ex1t.pdf',plot=p.comb_106_ext, width = 7, height = 10,unit='in')

p.comb_106_ext <- multitree_plot_ann(lsTrees=lsTrees.106.rooted[6], phylogroup='phylogroup',layout='rectangular',nrow=1,extensive=TRUE, tipsize=2.5,nodes=TRUE)
ggsave('~/DATA/MarbGenomics/Graphs/annotated_106_speciestrees_circular_ex16.pdf',plot=p.comb_106_ext, width = 7, height = 10,unit='in')


#=========================================================
# RETRUNS A GGPLOT OBJECT
tree_plot_ann(tr=lsTrees.106.rooted[[1]], phylogroup='phylogroup',layout='rectangular',nrow=1,extensive=TRUE, tipsize=2.5,nodes=TRUE)

#=========================================================


tree_plot_ann




base.tree <- sco.94.concat.tree.rooted
base.tree <- lsTrees.106.rooted[[1]]
p.base <- ggtree(base.tree, layout = 'circular')



metadata <- genome.tbl %>% filter(genome %in% base.tree$tip.label)
rownames(metadata) <- metadata$genome
p1 <- p.base %<+% metadata

p2 <- p1 #+
#geom_point(shape=21, mapping=aes(fill= group, size= 2))

habitat.colors = c(
  'Hydrothermal_Vent'="#E41A1C",
  'Lake'="#F37912",
  'Oil'="#FFD422",
  'Photic'="#43997A",
  'Phototroph'="#658E67",
  'Polar'="#5D6795",
  'Sediment'="#A35390",
  'Terrestial'="#B85F49",
  #''='white',
  'Other'='grey',
  'Tidal_flat'='brown',
  'Deep'='black',
  'salt_lake'='purple',
  'Saltern'='orange',
  'grey','grey20'
)


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

p2$data <- p2$data %>% mutate(SS = as.character(SS))%>% mutate(SS = ifelse(nchar(SS)<1,NA,SS))


p3 = p2 + new_scale_fill() +
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=phylogroup),
             color = "white", offset = 0.04,size = 0.02,width=0.03)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))+
  scale_fill_manual(values =phylogroup.colors)+
  new_scale_fill() +
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=SS),
             color = "white", offset = 0.04,size = 0.02,width=0.03)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))+
  scale_fill_manual(values =habitat.colors)

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

p7=  p5   +       geom_treescale(fontsize=2, linesize=0.3) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )

p7








#====================================================================================
e <- e %>% mutate(node=ID+1)
RED = castor::get_reds(sco.94.concat.tree.rooted)
nodeDepths = castor::get_all_node_depths(sco.94.concat.tree.rooted)

e$RED <- RED[2:(length(RED)-1)]
e$nodeDepths = nodeDepths<- RED[2:(length(RED)-1)]





p <- ggtree(sco.94.concat.tree.rooted)
p <- tree_plot_ann(tr=lsTrees.94.rooted[[1]], phylogroup='phylogroup',layout='rectangular',nrow=1)

p <- tree_plot_ann(tr=lsTrees.106.rooted[[1]], phylogroup='phylogroup',layout='rectangular',nrow=1)

anndt <- genome.tbl %>% left_join(tblAssembly, by = c("Accession"="BioSampleAccn"))

p$data <- p$data %>% left_join(e, by='node') %>%
  mutate(branchsupport=paste(bootstrap,gCF,sCF,sep='|')) %>%
  left_join(anndt, by = c("label"="genome")) %>%
  mutate(label = paste0( AssemblyAccession, " - ", Organism.x, ' ', strain))

p.o <- p + geom_text(data = p$data %>% filter(isTip != TRUE), aes(label = branchsupport), hjust = 1, vjust = -0.4, size = 2) +
  geom_tiplab(align=TRUE,size=2.5) + geom_treescale()

p.o

ggplot2::ggsave(filename=paste0(outdir,'tet333NOT_nuc_nodelabels','.pdf'),plot=p.o, width = 7, height = 12,unit='in')

p.o <- p + geom_text(data = p$data %>% filter(isTip != TRUE), aes(label = branchsupport), hjust = 1, vjust = -0.4, size = 2) +
  geom_tiplab(size=2.5) + geom_treescale()

p.o

ggplot2::ggsave(filename=paste0(outdir,'tet323NOT_nuc_nodelabels','.pdf'),plot=p.o, width = 7, height = 12,unit='in')

p.o$data <- p.o$data  %>% select(-Completeness)

plotdata <- genome.tbl %>%
  mutate(id=genome) %>%
  select(-genome) %>%
  select(id, everything()) %>%
  data.frame()



metadata <- genome.tbl %>% filter(genome %in% lsTrees.106.rooted[[1]]$tip.label) %>%
  select(-Completeness,-strain)%>%
  left_join(anndt, by = c("genome"="genome")) %>%
  mutate(genome = paste0( AssemblyAccession, " - ", Organism.x, ' ', strain))




#========
#Circular version
# still need to make more trees, say all trees, and annotate them better by including things like bootstrap

p <- tree_plot_ann(tr=lsTrees.106.rooted[[1]], phylogroup='phylogroup',layout='circular',nrow=1,extendto=1.22)

anndt <- genome.tbl %>% left_join(tblAssembly, by = c("Accession"="BioSampleAccn"))

#p$data <- p$data %>% left_join(e, by='node') %>%
#  mutate(branchsupport=paste(bootstrap,gCF,sCF,sep='|')) %>%
#  left_join(anndt, by = c("label"="genome")) %>%
#  mutate(label = paste0(Organism.x, ' ', strain))

p.o <- p + geom_tiplab(size=1.7) + geom_treescale()

p.o

ggplot2::ggsave(filename=paste0(outdir,'circular_aa_nodelabels','.pdf'),plot=p.o, width = 7, height = 10,unit='in')


#========

# %>% data.frame()
#rownames(metadata) <- metadata$genome
metadata <- metadata %>% mutate(ID=genome) %>% select(-genome) %>% select(ID, everything())

ggtree::facet_plot(p.o, panel = 's', data = plotdata %>% as_tibble(),
                   geom = geom_point,
                   mapping = aes( x= Completeness)) +
  theme_tree2()


p <- tree_plot_ann(tr=lsTrees.106.rooted[[1]], phylogroup='phylogroup',layout='rectangular',nrow=1)

p_gs <- ggtree::facet_plot(p, panel = 's', data = plotdata %>% as_tibble(),
                   geom = ggstance::geom_barh,
                   mapping = aes( x= Genome_size),stat='identity')

p_gs2 <- ggtree::facet_plot(p_gs, panel = 'g', data = plotdata %>% as_tibble(),
                           geom = geom_point,
                           mapping = aes( x= Coding_density)) +
  theme_tree2()
p_gs2


cogcat.plot <- COGCATpan %>% t() %>% data.frame ()%>% rownames_to_column(var = 'id') %>% pivot_longer(cols = c(-genome))


ggtree::facet_plot(p +   ggnewscale::new_scale_fill(), panel = 's', data = cogcat.plot %>% as_tibble(),
                   geom = ggstance::geom_barh,
                   mapping = aes( x= value,fill=name),stat='identity')


facet_plot(p.o, panel = ' Barplot', data = plotdata,
           geom = geom_point,
           mapping = aes(x = Completeness),
           stat='identity' )


p <- tree_plot_ann(tr=lsTrees.106.rooted[[1]], phylogroup='phylogroup',layout='circular',nrow=1)

p$data <- p$data %>% left_join(e, by='node') %>%
  mutate(branchsupport=paste(bootstrap,gCF,sCF,sep='|')) %>%
  left_join(anndt, by = c("label"="genome")) %>%
  mutate(label = paste0( AssemblyAccession, " - ", Organism.x, ' ', strain))

p.o <- p + geom_text(data = p$data %>% filter(isTip != TRUE), aes(label = branchsupport), hjust = 1, vjust = -0.4, size = 2) +
  geom_tiplab(align=TRUE, size=2.5) + geom_treescale()
p.o




p$data %>%
  ggplot(aes(RED, branchlength))+
  geom_point(aes(fill=gCF),shape=21,size=2) +
  viridis::scale_fill_viridis(option='viridis',direction = -1, discrete = FALSE)+
  fdb_style()


p$data %>%
  ggplot(aes(branchlength, gCF))+
  geom_point(aes(fill=RED),shape=21,size=2) +
  viridis::scale_fill_viridis(option='viridis',direction = -1, discrete = FALSE)+
  fdb_style()

testm <- glm(gCF~RED+branchlength+RED:branchlength, data=p$data)
summary(testm)


p + geom_text(data = p$data %>% filter(isTip != TRUE), aes(label = branchsupport), hjust = 1, vjust = -0.4, size = 3) + geom_tiplab(size=3) + geom_point(aes(color=gCF))
p + geom_point(aes(color=gCF,size=gCF))

p + geom_text(data = p$data %>% filter(isTip != TRUE), aes(label = branchsupport), hjust = 1, vjust = -0.4, size = 3) + geom_tiplab(size=3) + geom_point(aes(color=RED))


#===rib.94.concat.tree
z <- rib.94.concat.tree.rooted
x <- sco.94.concat.tree.rooted
y <- astral.94.tree.rooted

#y = nj.Man
p1 <- ggtree(x, layout='rectangular')
p2 <- ggtree(y)

d1 <- p1$data
d2 <- fortify(y)
d3 <- fortify(z)
d2$x <- d2$x + max(d1$x) + 1
d3$x <- d3$x + max(d2$x) + 1

dd = bind_rows(d1, d2, d3) %>%
  filter(!is.na(label))

p1 + geom_tree(data = d2) + geom_tree(data = d3) + geom_tiplab(data=d3) +
  geom_line(aes(x, y, group=label, color=node < 15), data=dd, alpha=.3)



#-------------------#
#z <- rib.94.concat.tree.rooted
y <- astral.94.tree.rooted
x <- astral.94.tree.rooted
y <- sco.94.concat.tree.rooted

x <- ape::ladderize(nj.GGR)#midpoint(nj.GGR) #%>% ape::ladderize()
x.name <- 'nj.GGR'

x <- ape::ladderize(nj.ANIb)#midpoint(nj.GGR) #%>% ape::ladderize()
x.name <- 'nj.ANIb'

x <- ape::ladderize(nj.AAI)#midpoint(nj.GGR) #%>% ape::ladderize()
x.name <- 'nj.AAI'


p1 <- ggtree(y, layout='rectangular') +
  geom_hilight(
    mapping=aes(subset = node %in% c(38, 48, 58, 36),
                node = node,
                fill = as.factor(node)
    )
  ) +
  labs(fill = "clades for tree in left" )

p2 <- ggtree(x)

d1 <- p1$data
d2 <- p2$data

#IF YOU WANT PANGENOME
#d2 <- d2 %>% mutate(branch.length=branch.length/20)
#d2 <- d2 %>% mutate(x=x/10)

# for astral
#d2 <- d2 %>% mutate(branch.length=branch.length/15)
#d2 <- d2 %>% mutate(x=x/15)

# for nj.ggr
d2 <- d2 %>% mutate(branch.length=branch.length*3)
d2 <- d2 %>% mutate(x=x*3)


## reverse x-axis and
## set offset to make the tree in the right hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2) +
  ggnewscale::new_scale_fill() +
  labs(fill = "clades for tree in right" )

dd <- bind_rows(d1, d2) %>%
  filter(!is.na(label)) %>% filter(isTip==TRUE) %>% left_join(genome.tbl, by=c('label'='genome'))

po <- pp + geom_line(aes(x, y, group=label,  color=phylogroup), data=dd) +
  geom_point(data = dd,aes(x, y, color=phylogroup), size=2)+
  ggtitle(paste0(x.name,' vs ML species tree'))+
  coord_flip()+scale_fill_manual(values=phylogroup.colors)+
  scale_color_manual(values=phylogroup.colors, na.value="lightgrey")+
  theme(legend.position='none')

po
ggplot2::ggsave(filename=paste0(outdir,'/',x.name,'_vs_species_tree','.pdf'), plot= po, width = 8, height = 2.75,unit='in')








#=====================================
#======



#================ SUPPORT CONCORDANCE STUFF
p4 <- ggtree(x, layout='circular',branch.length = 'none')

p4 <- tree_plot_ann(tr=lsTrees.94.rooted[[1]], phylogroup='phylogroup',layout='rectangular',nrow=1,branch.length = "none",highlight = FALSE)
#p4 <- ggtree(lsTrees.94.rooted[[1]], layout='rectangular',branch.length = 'none')

p4$data <- p4$data %>% left_join(genome.tbl, by=c('label'='genome')) %>%
  left_join(e, by=c('node'='ID'))
#p4 <- open_tree(p4, 180)
p.gCF <- p4 + aes(color=gCF) + scale_color_gradient(low = 'red',high = 'navy' ,limits = c(0,100))+theme(legend.position="bottom", legend.box = "horizontal") + ggtitle('gCF') + theme(legend.key.width = unit(0.55, 'cm'), legend.key.height = unit(0.25, 'cm'))   #+ geom_point(data = p4$data %>% filter(gCF<100), aes(color=gCF),size=1)
p.sCF <- p4 + aes(color=sCF) + scale_color_gradient(low = 'red',high = 'navy',limits = c(0,100) )+theme(legend.position="bottom", legend.box = "horizontal") + ggtitle('sCF') + theme(legend.key.width = unit(0.55, 'cm'), legend.key.height = unit(0.25, 'cm'))#+ geom_point(data = p4$data %>% filter(sCF<100),aes(color=sCF),size=1)
p.bootstrap <-p4 + aes(color=bootstrap) + scale_color_gradient(low = 'red',high = 'navy',limits = c(0,100) )+theme(legend.position="bottom", legend.box = "horizontal") + ggtitle('bootstrap') + theme(legend.key.width = unit(0.55, 'cm'), legend.key.height = unit(0.25, 'cm'))#+ geom_point(data = p4$data %>% filter(bootstrap<100),aes(color=bootstrap),size=1)

p.comparisonOfSupport <- gridExtra::grid.arrange(p.bootstrap, p.gCF,p.sCF, nrow=1)
ggplot2::ggsave(filename=paste0(outdir,'comparisonsOfSupport_2_noncolor','.pdf'), plot= p.comparisonOfSupport, width = 2.5, height = 4,unit='in')

#==================

concord.p1 <- e %>% ggplot2::ggplot(aes(x = branchlength, y = gCF)) +
  ggplot2::geom_point(aes(fill = bootstrap),shape=21,size=2) +
  scale_fill_gradient(low = 'red',high = 'navy' ) +
  ggplot2::ylim(0, 100) +
  geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), size = 0.5,color="red", se=F) +
  fdb_style()

concord.p2 <- e %>% ggplot2::ggplot(aes(x = branchlength, y = bootstrap)) +
  ggplot2::geom_point(aes(fill = bootstrap),shape=21,size=2) +
  scale_fill_gradient(low = 'red',high = 'navy') +
  ggplot2::ylim(0, 100) +
  fdb_style()

concord.p3 <- e %>% ggplot2::ggplot(aes(x = gCF, y = sCF)) +
  ggplot2::geom_point(aes(fill = bootstrap),shape=21,size=2) +
  scale_fill_gradient(low = 'red',high = 'navy') +
  ggplot2::xlim(0, 100) +
  ggplot2::ylim(0, 100) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  ggplot2::geom_smooth(method='lm', size = 0.5, color='red')+
  fdb_style()

concord.p4 <- e %>% ggplot2::ggplot(aes(x = branchlength, y = sCF)) +
  ggplot2::geom_point(aes(fill = bootstrap),shape=21,size=2) +
  scale_fill_gradient(low = 'red',high = 'navy') +
  ggplot2::ylim(0, 100) +
  fdb_style()


p.comb <- gridExtra::grid.arrange(concord.p3,concord.p2,concord.p1,concord.p4, nrow=1)
ggplot2::ggsave(filename=paste0(outdir,'concordance_factors_reds','.pdf'),plot=p.comb, width = 13, height = 4,unit='in')


gridExtra::grid.arrange( p.comparisonOfSupport, p.comb, ncol=1)




e %>% dplyr::mutate(gcf_bin =ifelse(gCF>50,TRUE,FALSE)) %>%
  ggplot2::ggplot(aes(x = gcf_bin, y = branchlength)) +
  ggplot2::geom_jitter(aes(fill = sCF),shape=21,size=2,width=0.2) +
  scale_fill_gradient(low = 'red',high = 'navy' ) +
  fdb_style(aspect.ratio=1.5)

e %>% dplyr::mutate(gcf_bin =ifelse(gCF>50,TRUE,FALSE)) %>%
  ggplot2::ggplot(aes(x = gcf_bin, y = RED)) +
  ggplot2::geom_jitter(aes(fill = sCF),shape=21,size=2,width=0.2) +
  scale_fill_gradient(low = 'red',high = 'navy' ) +
  fdb_style(aspect.ratio=1.5)





#========================================
p.astral <- ggtree(astral.94.tree.rooted, layout='rectangular',branch.length = 'none')
#p.astral <- open_tree(p.astral, 180)

#p.astral <- ggtree(astral.94.tree.rooted)
p.astral$data <- p.astral$data %>%
  mutate(support = label) %>%
  mutate(support=ifelse(isTip==TRUE,NA,support)) %>%
  mutate(support = as.numeric(as.character(support)))

p.astral <- p.astral + aes(color=support) +
  #geom_point(data=p.astral$data %>% filter(isTip==FALSE) %>% filter(support <1),aes(color=support)) +
  scale_color_gradient(low = 'red',high = 'navy') +
  #scale_color_distiller(palette='Spectral', direction = 1)+
  theme(legend.position="bottom", legend.box = "horizontal") +
  ggtitle('ASTRAL') +
  theme(legend.key.width = unit(0.55, 'cm'), legend.key.height = unit(0.25, 'cm'))

p.comparisonOfSupport <- gridExtra::grid.arrange(p.bootstrap, p.gCF,p.sCF, p.astral, nrow=1)
ggplot2::ggsave(filename=paste0(outdir,'comparisonsOfSupport_plusastral_rect_reds','.pdf'), plot= p.comparisonOfSupport, width = 4, height = 4,unit='in')



#======================================================================
# Extract numbers to report in the writing

non100boot <- e %>% filter(bootstrap!=100) %>%nrow()
alltips <- e %>% nrow()

#percentage of 100bootsrtap nodes
(alltips-non100boot)/alltips


posterior.prob <- as.numeric(astral.94.tree$node.label)
posterior.prob <- posterior.prob[!is.na(posterior.prob)]

length(posterior.prob[posterior.prob==1])/length(posterior.prob)


#how many are of short branch lenght, or RED
e %>% filter(bootstrap!=100) %>% filter(RED<0.98)
e %>% filter(bootstrap!=100) %>% filter(RED>0.98) %>% pull(RED) %>% mean()


e %>% ggplot(aes(RED,branchlength)) + geom_point()+fdb_style()


#======================================================
# DENSITREE
#======================================================
tiplabel.colors.df <- genome.tbl %>%
  dplyr::filter(genome %in% sco.94.concat.tree.rooted$tip.label) %>%
  dplyr::select(genome,phylogroup) %>%
  dplyr::mutate(tipcol = phylogroup.colors[phylogroup])

tiplabel.colors <- tiplabel.colors.df$tipcol
names(tiplabel.colors) <- tiplabel.colors.df$genome


names(sco.94.loci.trees) <- trees.meta$name

species.trees <- list('sco.94.concat.tree.rooted'=sco.94.concat.tree.rooted,
                      'rib.94.concat.tree.rooted'=rib.94.concat.tree.rooted,
                      'astral.94.tree'=astral.94.tree.rooted)

class(species.trees) = 'multiPhylo'

all.94.trees <- c(species.trees, sco.94.loci.trees)

length(all.94.trees)

# FOR LOCI 106

densi.trees <- sco.106.loci.trees
for(t.name in names(densi.trees)){
  densi.trees[[t.name]] <- phangorn::midpoint(densi.trees[[t.name]])
  densi.trees[[t.name]]$edge.length  <- NULL

}

#densiTree phangron, sorted along ML species tree!
grDevices::pdf(paste0(outdir,'densiTree_loci_trees.106','.pdf'), height = 7.5, width =20)
phangorn::densiTree(densi.trees,type="cladogram",
                    consensus = lsTrees.106.rooted[[1]],
                    col="navyblue",
                    direction="downwards",
                    width=1,
                    scaleX = TRUE,alpha=0.02)
dev.off()


# FOR LOCI
densi.trees <- sco.94.loci.trees
for(t.name in names(densi.trees)){
  densi.trees[[t.name]] <- phangorn::midpoint(densi.trees[[t.name]])
  densi.trees[[t.name]]$edge.length  <- NULL

}


#tiplabel.colors.df <- genome.tbl %>% filter(genome %in% sco.94.concat.tree.rooted$tip.label) %>% select(genome,phylogroup) %>% mutate(tipcol = phylogroup.colors[phylogroup])
#tiplabel.colors <- tiplabel.colors.df$tipcol
#names(tiplabel.colors) <- tiplabel.colors.df$genome

#densiTree phangron, sorted along ML species tree!
grDevices::pdf(paste0(outdir,'densiTree_loci_trees_94','.pdf'), height = 7.5, width =15)
phangorn::densiTree(densi.trees,type="cladogram",
                    consensus = sco.94.concat.tree.rooted,
                    col="black",
                    direction="downwards",
                    width=1,
                    scaleX = TRUE,alpha=0.02,
                    tip.color = tiplabel.colors[tip_order(sco.94.concat.tree.rooted)])
dev.off()

# FOR Species trees
densi.trees <- species.trees
densi.trees <- lsTrees.rooted


for(t.name in names(densi.trees)){
  densi.trees[[t.name]]$edge.length  <- NULL
}

grDevices::pdf(paste0(outdir,'densiTree_species_trees','.pdf'), height = 7.5, width =15)
phangorn::densiTree(densi.trees,type="cladogram",
                    consensus = sco.94.concat.tree.rooted,
                    col="red",
                    direction="downwards",
                    width=1,
                    scaleX = TRUE,alpha=0.5)
dev.off()

#≠===========
densi.trees <- c(nj.AAI,nj.ANIb,nj.GGR,nj.Man)

for(t.name in names(densi.trees)){
  densi.trees[[t.name]]$edge.length  <- NULL
}

grDevices::pdf(paste0(outdir,'densiTree_ANIetc_trees','.pdf'), height = 7.5, width =20)
phangorn::densiTree(densi.trees,type="cladogram",
                    consensus = sco.94.concat.tree.rooted,
                    col="red",
                    direction="downwards",
                    width=1,
                    scaleX = TRUE,alpha=0.5)
dev.off()

#=============

densi.trees <- sco.94.loci.trees[names(sco.94.loci.trees[intersect(rib.SCO, names(sco.94.loci.trees))])]
for(t.name in names(densi.trees)){
  densi.trees[[t.name]] <- phangorn::midpoint(densi.trees[[t.name]])
  densi.trees[[t.name]]$edge.length  <- NULL
}

#densiTree phangron, sorted along ML species tree!
grDevices::pdf(paste0(outdir,'densiTree_rib.loci_trees','.pdf'), height = 7.5, width =15)
phangorn::densiTree(densi.trees,type="cladogram",
                    consensus = sco.94.concat.tree.rooted,
                    col="navyblue",
                    direction="downwards",
                    width=1,
                    scaleX = TRUE,alpha=0.1)
dev.off()



#densiTree phangron, sorted along ML species tree!
grDevices::pdf(paste0(outdir,'densiTree_rib.loci_trees_up','.pdf'), height = 7.5, width =15)
phangorn::densiTree(densi.trees,type="cladogram",
                    consensus = sco.94.concat.tree.rooted,
                    col="navyblue",
                    direction="upwards",
                    width=1,
                    scaleX = TRUE,alpha=0.1)
dev.off()



#============================================================================================================
#not normalised


##############################################################
#====== EXCLUSIVITY
#exclusivity
minimum_outgroup
maximim_ingroup


x <- sco.94.concat.tree.rooted
p1 <- ggtree(x,branch.length='none')

coph.p <- cophenetic(x)

internal <- p1$data %>% filter(isTip == FALSE)  %>% filter(label!="Root")
exlusivity <- NULL
for(node in internal$node){
  message('analysing', node)

  ingroup <- p1$data %>% tidytree::offspring(node) %>% filter(isTip==TRUE) %>% select(label) %>% pull
  max.ingroup <- max(coph.p[ingroup,ingroup])
  min.outgroup <- min(coph.p[ingroup,setdiff(rownames(coph.p),ingroup)])
  exl <- c('node'=node, 'exl' = min.outgroup - max.ingroup)
  exlusivity <- rbind(exlusivity, exl)
}

exlusivity <- exlusivity %>% data.frame()
p1$data <- p1$data %>% left_join(exlusivity, by='node')

p1 + aes(color=exl) +
  #geom_point(data=p.astral$data %>% filter(isTip==FALSE) %>% filter(support <1),aes(color=support)) +
  #scale_color_gradient(low = 'navyblue',high = 'deepskyblue') +
  #scale_color_distiller(palette='Reds', direction = 1)+
  colorspace::scale_colour_continuous_diverging(palette = "Blue-Red 2",
                                                l1 = 30,
                                                l2 = 70,
                                                p1 = .7,
                                                p2 = 1.5)+
  theme(legend.position="bottom", legend.box = "horizontal") +
  ggtitle('exclusivity') +
  theme(legend.key.width = unit(0.55, 'cm'), legend.key.height = unit(0.25, 'cm'))


p1$data %>% ggplot(aes(branch.length,exl)) + geom_point()+fdb_style()

#=========



x <- sco.94.concat.tree.rooted
p1 <- ggtree(x,branch.length='none')

coph.p <- cophenetic(x)

internal <- p1$data %>% filter(isTip == FALSE)  %>% filter(label!="Root")

grps <- genome.tbl %>%
  dplyr::filter(!is.na(phylogroup)) %>%
  dplyr::select(phylogroup) %>%
  dplyr::pull() %>%
  unique() %>%
  as.character()

grps <- genome.tbl %>% select(-phylogroup) %>%
  left_join(sorted.cliques,by='genome') %>%
  dplyr::filter(!is.na(phylogroup)) %>%
  dplyr::select(sc0.8) %>%
  dplyr::pull() %>%
  unique() %>%
  as.character() %>% paste0('cl',.)
#is an input

sel.col <- 'phylogroup'




distD <- coph.p
distD <- 1-ANIb[sel.genome,sel.genome]
distD <- 1-pan.ggr[sel.genome,sel.genome]

exclusion.out <- NULL
for(sel.col in c('phylogroup','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  mtdt <- genome.tbl %>%
    select(-phylogroup) %>%
    left_join(sorted.cliques,by='genome') %>%
    dplyr::filter(!is.na(phylogroup)) %>% data.frame()

  mtdt <- mtdt[!is.na(mtdt[,sel.col]),]

  grps <- mtdt %>%
    dplyr::select(sel.col) %>%
    dplyr::pull() %>%
    unique() %>%
    as.character() #%>% #paste0('cl',.)

  exlusivity <- NULL
  for(grp in grps){
    message('analysing', grp)
    ingroup <- mtdt[mtdt[,sel.col]==grp,'genome']
    max.ingroup <- max(distD[ingroup,ingroup])
    min.outgroup <- min(distD[ingroup,setdiff(rownames(distD),ingroup)])
    exl <- c('grp'=grp, 'exl' = min.outgroup - max.ingroup)
    exlusivity <- rbind(exlusivity, exl)
  }

  exlusivity <- exlusivity  %>% data.frame() %>% mutate(level=sel.col)
  exclusion.out <- rbind(exclusion.out, exlusivity)

}
exclusion.out <- exclusion.out %>% mutate(exl = as.numeric(as.character(exl)))


#p1$data <- p1$data %>% left_join(exlusivity, by='node')
#p1$data %>% ggplot(aes(branch.length,exl)) + geom_point()+fdb_style()

exclusion.out %>% ggplot(aes(x=exl,y=level)) + geom_point()+geom_vline(xintercept=0) + fdb_style()








sel.col = 'sc0.8'

mtdt <- genome.tbl %>%
  select(-phylogroup) %>%
  left_join(sorted.cliques,by='genome') %>%
  dplyr::filter(!is.na(phylogroup)) %>% data.frame()

mtdt <- mtdt[!is.na(mtdt[,sel.col]),]

if(grepl('sc',sel.col)){
  mtdt[,sel.col] <- paste0('cl',mtdt[,sel.col])
}

rownames(mtdt) <- mtdt$genome

metadata

grps <- mtdt %>%
  dplyr::select(sel.col) %>%
  dplyr::pull() %>%
  unique() %>%
  as.character() #%>% #paste0('cl',.)

cazy.dist <- pan2dist(OG.pan[mtdt$genome,], dist.method="Manhattan")

pan2pcoa(dist.pan = cazy.dist,
         metadata=mtdt,
         color=sel.col,
         title_suffix=sel.col,
         color_scale = unname(as.character(unique(out.colors[,'col.cc0.8']))))

#cazy.dist <- pan2dist(OG.pan[sel.genome,], dist.method="Manhattan")
#pan2pcoa(dist.pan = cazy.dist,
#         metadata=mtdt,
#         color="phylogroup",
#         title_suffix="CAZY") + scale_fill_manual(values=phylogroup.colors)









##############################################################
#====== PREPARING CONFIG FILES

# preparing input files for angst analysis, one per gene tree

out.path = '/Users/sa01fd/DATA/MarbGenomics/angst/'
ultrametric = 'True'
species.tree = '/Users/sa01fd/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile'
genes.dir = '/Users/sa01fd/DATA/MarbGenomics/loci_boottrees'
penalties.path = '/Users/sa01fd/Genomics/angst/example/penalty.file'
config.dir <- '/Users/sa01fd/DATA/MarbGenomics/angst_configs/'


#   complicated now, i read locally, but prepare config files to work on CLIMB
#   on CLIMB
#   /home/ubuntu/fred/angst/
out.path = '/home/ubuntu/fred/angst/output/'
ultrametric = 'True'
species.tree = '/home/ubuntu/fred/angst/angst.species.treefile'
genes.dir = '/Users/sa01fd/DATA/MarbGenomics/loci_boottrees'
genes.dir.out = '/home/ubuntu/fred/angst/loci_boottrees'
penalties.path = '/home/ubuntu/fred/angst/example/penalty.file'
config.dir <- '/Users/sa01fd/DATA/MarbGenomics/angst_configs/'


files <- list.files(genes.dir,pattern='.tree',full.names = TRUE)
for(file in files){
  gene.name = gsub('.tree','',basename(file))
  output.path <- paste0(out.path, gene.name)
  out.lines <-c(paste0('species=',species.tree),
                paste0('gene=',paste0(genes.dir.out,'/',basename(file))),
                paste0('penalties=',penalties.path),
                paste0('output=',output.path),
                paste0('ultrametric=',ultrametric))
  file.out = paste0(config.dir,gene.name,'.input')
  writeLines(out.lines,file.out)
  }


##############################################################
#====== LOADING THE OUTPUT FILES
# COPY FILES BACK FROM CLIMB

#ANGsT2tree
#===============================================================
# Function to map ANGST to tree !
#' Title
#'
#' @param angst.intree
#' @param angst.events
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'

angst.intree
angst.events= paste0(dir,'/AnGST.events')
dir
ANGsT2tree <- function(angst.intree, angst.events, dir){
  gene.product <- annotation.tbl %>%
    dplyr::filter(OG==basename(dir)) %>%
    dplyr::select(product) %>%
    dplyr::count(product) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::slice(1) %>%
    dplyr::select(product) %>%
    dplyr::pull()
  angst.events <- base::readLines(angst.events)

  p <- ggtree(angst.intree, branch.length = "none")
  #p$data <- p$data %>% left_join(names.conv, by=c('label'='angst'))
  #p <- ggtree(angst.intree, branch.length = "none")

  #loading ANGsT data
  #cur <- extant gene copy
  cur <- angst.events[grepl('\\[cur\\]',angst.events)] %>% gsub('\\[cur\\]: ','',.)
  cur <- cur %>% table() %>% data.frame() %>% `colnames<-`(c("tip", "cur"))
  #brn <- gene family birth (say hirizontal aq)
  brn <- angst.events[grepl('\\[brn\\]',angst.events)] %>% gsub('\\[brn\\]: ','',.) %>% strsplit(.,split='-')
  #spc <- speciation
  spc <- angst.events[grepl('\\[spc\\]',angst.events)] %>% gsub('\\[spc\\]: ','',.) %>% strsplit(.,split='-')
  #los <- gene loss
  los <- angst.events[grepl('\\[los\\]',angst.events)] %>% gsub('\\[los\\]: ','',.) %>% strsplit(.,split='-')
  #dup <- gene duplication
  dup <- angst.events[grepl('\\[dup\\]',angst.events)] %>% gsub('\\[dup\\]: ','',.) %>% strsplit(.,split='-')
  #hgt <- horizontal transfer
  hgt <- angst.events[grepl('\\[hgt\\]',angst.events)] %>% gsub('\\[hgt\\]: ','',.) %>% strsplit(.,split=' --> ')

  #assessing node in phylogeny associated with events
  #BRN
  mrca_nodes <- NULL
  for(i in 1:length(brn)){
    mrca_node <- p$data %>% tidytree::MRCA(brn[[i]]) %>% select(node) %>% pull()
    mrca_nodes <- c(mrca_nodes, mrca_node)
  }
  brn_nodes = mrca_nodes

  #LOS
  mrca_nodes <- NULL
  for(i in 1:length(los)){
    mrca_node <- p$data %>% tidytree::MRCA(los[[i]]) %>% select(node) %>% pull()
    mrca_nodes <- c(mrca_nodes, mrca_node)
  }
  los_nodes = mrca_nodes

  #DUP
  if(length(dup)!=0){
    mrca_nodes <- NULL
    for(i in 1:length(dup)){
      mrca_node <- p$data %>% tidytree::MRCA(dup[[i]]) %>% select(node) %>% pull()
      mrca_nodes <- c(mrca_nodes, mrca_node)
    }
    dup_nodes = mrca_nodes
  }else{dup_nodes<-NULL}


  mrca_nodes <- NULL
  for(i in 1:length(spc)){
    mrca_node <- p$data %>% tidytree::MRCA(spc[[i]]) %>% select(node) %>% pull()
    mrca_nodes <- c(mrca_nodes, mrca_node)
  }
  spc_nodes = mrca_nodes

  # Assessing nodes from and to for HGT
  hgt.df<-NULL
  for(i in hgt){
    hgt.df<- rbind(hgt.df , c('from'=i[1],'to'=i[2]))
  }

  hgt.df %>% as_tibble()
  f_nodes <- NULL
  t_nodes <- NULL
  for(i in 1:nrow(hgt.df)){
    from <- hgt.df[i,'from']
    mrca_from <- p$data %>% tidytree::MRCA(unname(unlist(strsplit(from, split='-')))) %>% select(node) %>% pull()
    f_nodes <- c(f_nodes, mrca_from)
    to <- hgt.df[i,'to']
    mrca_to <- p$data %>% tidytree::MRCA(unname(unlist(strsplit(to, split='-')))) %>% select(node) %>% pull()
    t_nodes <- c(t_nodes, mrca_to)
  }

  hgt_nodes <- cbind('from'=f_nodes,'to'=t_nodes)

  #count the arrivals of hgt
  count.hgt.to <- data.frame(table(hgt_nodes[,'to']))
  colnames(count.hgt.to) = c('node','hgt.to.count')
  count.hgt.to <- count.hgt.to %>% mutate(node=as.numeric(as.character(node)))

  p$data <- p$data %>%
    dplyr::mutate(los = ifelse(node %in% los_nodes, 1, NA)) %>%
    dplyr::mutate(dup = ifelse(node %in% dup_nodes, 1, NA)) %>%
    dplyr::mutate(brn = ifelse(node %in% brn_nodes, 1, NA)) %>%
    dplyr::mutate(spc = ifelse(node %in% spc_nodes, 1, NA)) %>%
    dplyr::left_join(cur, by=c('label'='tip')) %>%
    dplyr::left_join(count.hgt.to, by='node')

  p2 <- p
  #hgt_nodes <- data.frame(hgt_nodes)
  #for(i in 1:nrow(hgt_nodes)){
  #  p2 <- p2 +
  #    ggtree::geom_taxalink(unname(hgt_nodes[i,1]),
  #                          color="blue3",
  #                          unname(hgt_nodes[i,2]),
  #                          curvature = 0.5)
  #

  p2 <- p2 +
    ggtree::geom_point(aes(size=los),color='black') +
    ggtree::geom_tippoint(aes(color=cur)) +
    ggtree::geom_point(aes(size=brn),color='red') +
    ggtree::geom_point(aes(size=hgt.to.count),color='blue')

  if(!is.null(dup_nodes)){
    p2 <- p2 + ggtree::geom_point(aes(size=dup),color='green')
  }

  p2 <- p2 + ggtitle(paste0('ANGsT_',gene,' ',gene.product))

  ggplot2::ggsave(filename=paste0('/Users/sa01fd/DATA/MarbGenomics/Graphs/angst/','ANGsT_',gene,'.pdf'),plot=p2, width = 5, height = 10,unit='in')
  return(p)
}

#===============================================================
#loading in the angst output files
# make tree dataframe based on species tree
# convert the
angst.intree <- sco.94.concat.tree.rooted
angst.intree$tip.label <-  gsub('_','',gsub('\\.','',angst.intree$tip.label))
names.conv <- data.frame(cbind('original'=sco.94.concat.tree.rooted$tip.label,'angst'=angst.intree$tip.label))


angst.events <- readLines(paste0('/Users/sa01fd/DATA/MarbGenomics/angst/output/','OG0000246','/AnGST.events'))

alldir <- list.dirs(path = "/Users/sa01fd/DATA/MarbGenomics/angst/output", full.names = TRUE)
alldir <- alldir[2:length(alldir)]
length(alldir)

#======================================================================
#     THIS IS WHERE ALL ANALYSIS IS LOADED AND MERGED INTO BIG TABLE
#======================================================================
allgene.p <- NULL
for(dir in alldir){
  gene <- basename(dir)
  message('loading',gene)
  p.m <- ANGsT2tree(angst.intree,paste0(dir,'/AnGST.events'), dir)
  allgene.p <- rbind(allgene.p, p.m$data %>% mutate(OG=gene))
  # PLOT THE LINKS FOR GIVEN ANALYSIS
}

save(allgene.p, file='~/DATA/MarbGenomics/allgene_angst.RData')



#allgene.p

allgene.p %>% filter(!is.na(los)) %>%
  group_by(node,OG) %>%
  summarise(sum=sum(los))

#========================
#allgene.p

load(file='~/DATA/MarbGenomics/allgene_angst.RData')


filter(OG %in% sig.genes) %>%

  sig.genes2 <- annotation.tbl %>% filter(grepl("nuo",Name)) %>% filter(!is.na(OG)) %>% select(OG) %>% unique() %>% pull() %>% as.character()


#========================
#culculate totals of all events
totals <- allgene.p %>%
  #FILTER OR NOT , make function!
  filter(OG %in% sig.genes2) %>%
  mutate(los = ifelse(!is.na(los),los,0)) %>%
  mutate(dup = ifelse(!is.na(dup),dup,0)) %>%
  mutate(brn = ifelse(!is.na(spc),brn,0)) %>%
  mutate(spc = ifelse(!is.na(spc),spc,0)) %>%
  mutate(cur = ifelse(!is.na(cur),cur,0)) %>%
  mutate(hgt.to.count = ifelse(!is.na(hgt.to.count),hgt.to.count,0)) %>%
  group_by(node) %>%
  dplyr::summarise(los = sum(los),
                   dup = sum(dup),
                   brn = sum(brn),
                   spc = sum(spc),
                   cur = sum(cur),
                   hgt.to.count = sum(hgt.to.count))

p <- ggtree(angst.intree, branch.length = "none")
p$data <- p$data %>%
  dplyr::left_join(totals, by='node')

#p +   ggtree::geom_point(aes(size=los),color='black') + geom_text(aes(label=los),color='white',size=4)# +
#p +   ggtree::geom_point(aes(size=dup),color='black') + geom_text(aes(label=dup),color='white',size=4)# +
#p +   ggtree::geom_point(aes(size=spc),color='black') + geom_text(aes(label=spc),color='white',size=4)# +
#p +   ggtree::geom_point(aes(size=cur),color='black') + geom_text(aes(label=cur),color='white',size=4)# +
#p +   ggtree::geom_point(aes(size=hgt.to.count),color='black') + geom_text(aes(label=hgt.to.count),color='white',size=4)# +

#if one wants to ignore the events at the tips, we can use p2,
#probably this is not desirable
p2 <- p
p2$data <- p2$data %>% mutate(hgt.to.count = ifelse(isTip==TRUE,NA,hgt.to.count))


p.hgt <- p + aes(color=hgt.to.count) +ggtree::geom_point(data = p$data %>% filter(isTip==FALSE), aes(size=hgt.to.count,color=hgt.to.count)) +
  scale_colour_gradient(low = "black", high = "green") +
  scale_size(range=c(0.2,3))+ggtitle('Intra-HGT')

p.los <- p + aes(color=los) +
  ggtree::geom_point(data = p$data %>% filter(isTip==FALSE), aes(size=los,color=los)) +
  scale_colour_gradient(low = "navyblue", high = "yellow")+
  scale_size(range=c(0.2,3))+ggtitle('Losses')

p.dup <- p + aes(color=dup) +
  ggtree::geom_point(data = p$data %>% filter(isTip==FALSE), aes(size=dup,color=dup)) +
  scale_colour_gradient(low = "blue", high = "red")+
  scale_size(range=c(0.2,3))+ggtitle('Duplications')

p.spc <- p + aes(color=spc) +
  ggtree::geom_point(data = p$data %>% filter(isTip==FALSE), aes(size=spc,color=spc)) +
  scale_colour_gradient(low = "deepskyblue", high = "deeppink")+
  scale_size(range=c(0.2,3))+ggtitle('Origination')

p.cur <- p + aes(color=cur) +
  ggtree::geom_point(data = p$data %>% filter(isTip==FALSE), aes(size=cur,color=cur)) +
  scale_colour_gradient(low = "navyblue", high = "yellow")+
  scale_size(range=c(0.2,3))

p.brn <- p + aes(color=brn) +
  ggtree::geom_point(data = p$data %>% filter(isTip==FALSE), aes(size=brn,color=brn)) +
  scale_colour_gradient(low = "navyblue", high = "yellow")+
  scale_size(range=c(0.2,3))

p.angst.comb <- gridExtra::grid.arrange(p.los, p.dup,p.hgt, nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/angst_sub2.pdf',plot=p.angst.comb, width = 7, height = 7,unit='in')

p.angst.comb <- gridExtra::grid.arrange(p.los, p.dup,p.hgt,p.spc, nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/angst_sub3.pdf',plot=p.angst.comb, width = 13, height = 7,unit='in')



#========================================================
p.dup + geom_text(aes(label=node),hjust=1.2,vjust=-0.2)


sel.og  <- allgene.p %>% filter(node==147) %>% filter(!is.na(dup)) %>% select(OG) %>% pull()
clusdetect.loc <- multigenedetect(OG=sel.og, genomes='Marinobacter_algicola_DG893')

annotation.tbl %>%
  filter(genome=='Marinobacter_algicola_DG893') %>%
  filter(OG %in% sel.og) %>%
  select(OG, COG, gene, product, eC_number) %>%
  dplyr::arrange(OG) %>% data.frame()


annotation.tbl %>%
  filter(genome=='Marinobacter_adhaerens_HP15') %>%
  filter(OG %in% sel.og) %>%
  select(OG, COG, gene, product, eC_number) %>%
  dplyr::arrange(OG) %>% data.frame()


sel.og  <- allgene.p %>% filter(node==172) %>% filter(!is.na(dup)) %>% select(OG) %>% pull()
p$data %>% tidytree::offspring(110)

clusdetect.loc <- multigenedetect(OG=sel.og, genomes='Marinobacter_adhaerens_HP15')

Marinobacter_antarcticusCGMCC110835
annotation.tbl %>%
  filter(genome=='Marinobacter_adhaerens_HP15') %>%
  filter(OG %in% sel.og) %>%
  select(OG, COG, gene, product, eC_number) %>%
  dplyr::arrange(OG) %>% data.frame()





sel.og  <- allgene.p %>% filter(node==117) %>% filter(!is.na(dup)) %>% select(OG) %>% pull()
clusdetect.loc <- multigenedetect(OG=sel.og, genomes='Marinobacter_hydrocarbonoclasticus_VT8')
clusdetect.loc %>% select(product, OG, prox.cluster, locus_tag)


annotation.tbl %>%
  filter(genome=='Marinobacter_hydrocarbonoclasticus_VT8') %>%
  filter(OG %in% sel.og) %>%
  select(OG, COG, gene, product, eC_number) %>%
  dplyr::arrange(OG) %>% data.frame()




ribosomal.ogs <- annotation.tbl %>%
  filter(grepl('ribosomal protein',product)) %>%
  select(OG) %>% filter(!is.na(OG)) %>% unique() %>% pull
ribosomal.ogs




plot.angst.counts<- function(sel.ogs, highlighted.ogs){
  p.count <- allgene.p %>% filter(OG %in% sel.ogs) %>%
    mutate(los=ifelse(!is.na(los),los,0)) %>%
    mutate(dup=ifelse(!is.na(dup),dup,0)) %>%
    mutate(brn =ifelse(!is.na(brn ),brn ,0)) %>%
    mutate(spc=ifelse(!is.na(spc),spc,0)) %>%
    mutate(cur=ifelse(!is.na( cur), cur,0)) %>%
    mutate(hgt.to.count=ifelse(!is.na(hgt.to.count ),hgt.to.count ,0)) %>%
    group_by(OG) %>% dplyr::summarise(los=sum(los),
                                      dup=sum(dup),
                                      #brn=sum(brn),
                                      #spc=sum(spc),
                                      #cur=sum(cur),
                                      hgt.to.count=sum(hgt.to.count)) %>%
    #tidyr::pivot_longer(cols=c(los,dup,brn,spc, cur, hgt.to.count), names_to='event') %>%
    tidyr::pivot_longer(cols=c(los,dup, hgt.to.count), names_to='event') %>%
    mutate(sub = ifelse(OG %in% highlighted.ogs,TRUE,FALSE)) %>%
    ggplot(aes(event,value)) + geom_jitter(alpha=0.5,shape=21,width=0.3) +geom_boxplot(outlier.shape = NA,alpha=0.4,width=0.5)+
    ggtitle(paste0("n=",length(sel.ogs)))+scale_fill_manual(values=c('black','orange'))+
    facet_wrap(~sub,nrow=1)+
    fdb_style(aspect.ratio=2.5)+
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    return(p.count)
}

p.SCO <- plot.angst.counts(SCO, rib.SCO)

ggsave('~/DATA/MarbGenomics/Graphs/angst_SCO2.pdf',plot=p.SCO, width = 5, height = 4,unit='in')



#===================================================================================
#
#
# RERUN THIS ANALYSIS USING NORMALISED RF VALUES
# NOW WE ARE OBSCURING THIS I FEEL WITH REAL VALUES

RFoulds <- phytools::multiRF(all.94.trees)

colnames(RFoulds) = names_of_trees
rownames(RFoulds) = names_of_trees

#normalised
#The normalized Robinson-Foulds distance is derived by dividing d(T_1, T_2) by the maximal possible distance i(T_1) + i(T_2). If both trees are unrooted and binary this value is 2n-6
RFoulds <- as.matrix(phangorn::RF.dist(all.94.trees, normalize = TRUE, check.labels = TRUE))
hist(as.matrix(RFoulds))

save(RFoulds, file='~/DATA/MarbGenomics/RFoulds_94.RData')



#PCoA
MDS = cmdscale(RFoulds, eig = TRUE, x.ret=TRUE)
MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
MDSdata = data.frame(MDS$points)
MDSdata$Name = rownames(MDSdata)
MDSdata$Name[1:3]= c('sco.concat','rib.concat','astral.concat')
MDSdata$type = as.numeric(grepl('concat',MDSdata$Name))
MDSdata$type = ifelse(MDSdata$type == 1, 'species-tree','gene-tree')
colnames(MDSdata) = c("x",'y','Name','type')

p <- ggplot(MDSdata,aes(x=x,y=y,label=Name)) +
  theme_classic() +
  xlab(paste('PCoA1 (', MDS.var.perc[1], '%',')',sep='')) +
  ylab(paste('PCoA2 (', MDS.var.perc[2], '%',')',sep='')) +
  labs(fill = "Genus") +
  ggtitle(label='Robinsons Foulds')+
  geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  #geom_hex()+
  geom_point(aes(fill=type),shape = 21, size = 1.5,alpha=.7) +
  scale_fill_manual(values=c('black','red'))+
  geom_text_repel(data=MDSdata[1:3,], aes(label=Name))+
  fdb_style()


p
ggsave('~/DATA/MarbGenomics/Graphs/RFoulds_647+3speciestrees.pdf',plot=p, width = 4, height = 4,unit='in')


#NO convergence after 500 iterations, seems odd?
MDSdata <- pan2metaMDS(RFoulds, trymax=500)



#===================================================================================
# CORRELATION WITH RED
p <- p2
RED = castor::get_reds(angst.intree)
REDdf <- data.frame(cbind(RED, 'node' = (1:length(RED)) + Ntip(angst.intree)))
p$data <- p$data %>% left_join(REDdf,by='node')

p$data %>% filter(RED!=0) %>%
  ggplot(aes(RED, los)) +
  geom_point() +
  geom_smooth(method='lm',color='red')+
  fdb_style()

p$data %>% pivot_longer(cols=c(los,dup,spc,hgt.to.count),names_to='event')%>%
  filter(RED!=0) %>%
  ggplot(aes(RED, value)) +
  geom_point() +  ggpubr::stat_cor()+
  geom_smooth(method='lm',color='red') + facet_wrap(~event,nrow=1,scales = 'free_y') +
  fdb_style()



p$data %>% pivot_longer(cols=c(los,dup,spc,branch.length),names_to='event')  %>%#%>% filter(RED!=0) %>%
  ggplot(aes(value,hgt.to.count)) +
  geom_point(aes(color=RED)) +  ggpubr::stat_cor() +
  geom_smooth(method='lm',color='red') +
  facet_wrap(~event,nrow=1,scales = 'free') +
  fdb_style()

# INCLUDING THE TIPS, HOW DOES IT CHANGE
p$data %>% pivot_longer(cols=c(los,dup,spc,branch.length),names_to='event')  %>%
  ggplot(aes(value,hgt.to.count)) +
  geom_point(aes(color=RED)) +  ggpubr::stat_cor() +
  geom_smooth(method='lm',color='red') +
  facet_wrap(~event,nrow=1,scales = 'free') +
  fdb_style()


#====
# TIP ONLY
p$data %>% filter(isTip==TRUE) %>% pivot_longer(cols=c(los,dup,cur,branch.length),names_to='event')  %>%
  ggplot(aes(value,hgt.to.count)) +
  geom_point(aes(color=RED)) +  ggpubr::stat_cor() +
  geom_smooth(method='lm',color='red') +
  facet_wrap(~event,nrow=1,scales = 'free') +
  fdb_style()

angst.metdat <- genome.tbl %>% mutate(angs.name = gsub('_','',genome))

#merge with genome table
p$data %>% filter(isTip==TRUE) %>% left_join(angst.metdat, by = c('label'='angs.name')) %>%
  pivot_longer(cols=c(Completeness,Contamination,Genome_size,GC),names_to='event')  %>%
  ggplot(aes(value,hgt.to.count))+geom_point(aes(color=RED)) +  ggpubr::stat_cor() +
  geom_smooth(method='lm',color='red') +
  facet_wrap(~event,nrow=1,scales = 'free') +
  fdb_style()

p$data %>% filter(isTip==TRUE) %>% left_join(angst.metdat, by = c('label'='angs.name')) %>%
  pivot_longer(cols=c(Completeness,Contamination,Genome_size,GC),names_to='event')  %>%
  ggplot(aes(hgt.to.count, lat)) + geom_point(aes(color=RED))+
  fdb_style()

#can we sum up events? like since ancestor?
# or is that a recicoulous proposition?
# the outgroup would have very little changes, which may not be really so?






p$data %>% pivot_longer(cols=c(los,hgt.to.count,spc,branch.length),names_to='event')%>% filter(RED!=0) %>%
  ggplot(aes(value,dup)) +
  geom_point(aes(color=RED)) +  ggpubr::stat_cor() +
  geom_smooth(method='lm',color='red') +
  facet_wrap(~event,nrow=1,scales = 'free') +
  fdb_style()


corrdf <- p$data %>% filter(isTip==FALSE) %>% filter(RED!=0) %>% select(branch.length,los,dup,spc,hgt.to.count,RED)
GGally::ggpairs(corrdf)






p2 <- p

p2 <- p2 + ggtree::geom_point(aes(size=los),color='blue')+
  #ggnewscale::new_scale() +
  ggtree::geom_tippoint(aes(color=cur)) +
  #ggnewscale::new_scale() +
  ggtree::geom_point(aes(size=brn),color='red') +
  ggtree::geom_point(aes(size=dup),color='green')
#p + ggtree::geom_point(aes(size=spc),color='red')# +
p2




p$data %>% tidytree::offspring(118)






#map that to the tree!

spc <- cur %>% table() %>% data.frame() %>% `colnames<-`(c("tip", "cur"))







# We check in the the overall tree with outgroup (hahella and co) to see what is a propper outgroup in the Marinobcter only phylogenies
# this is needed for ANGST


#this Marinobacter_A clade is clearly monophyletic even in out sub tree
ape::is.monophyletic(sco.concat.tree, tips = Marinobacter_A)
ape::is.monophyletic(rib.concat.tree, tips = Marinobacter_A)

#we this assume that this is a good clade to root the tree on
sco.concat.tree.rooted <- ape::root(sco.concat.tree, outgroup = Marinobacter_A)
rib.concat.tree.rooted <- ape::root(rib.concat.tree, outgroup = Marinobacter_A)

ape::write.tree(rib.concat.tree.rooted, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/concat.rooted.treefile' )
ggtree::ggtree(sco.concat.tree.rooted) + geom_tiplab()

ggtree::ggtree(rib.concat.tree.rooted) + geom_tiplab()

rib.concat.tree




#Bloody hell
# IMPLEMENT AN ANGST MAKER
#needs to be rooted

#IN SPECIES TREE, NO underscores, spaces or dots
angst.species.tree <- sco.94.concat.tree.rooted
angst.species.tree$tip.label <- gsub('_','',angst.species.tree$tip.label)

#ape::is.rooted(angst.species.tree)
#ape::is.rooted(angst.species.tree)
#angst.species.tree <- ape::root.phylo(angst.species.tree, outgroup = gsub('_','',Marinobacter_A))

#boostraps are ignored, so I remove them
angst.species.tree$node.label <- NULL

#all edge lenghts need to be >0, so I add rediculously low value to everything
angst.species.tree$edge.length <- angst.species.tree$edge.length + 0.0001


angst.species.tree <- ape::read.tree('~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile')

#This makes a single-number vector that represents the length of the consensus$edge.length
#n <- length(angst.species.tree$edge.length)

#need to reduce the number of elements by 1... We also add a very small amount to make sure no branch length is zero.
#angst.species.tree$edge.length <- angst.species.tree$edge.length + .0000001
#print(angst.species.tree$edge.length)
#how does it look now?  It should not have a zero for the last row.  Plot to confirm it works.
#plotTree(angst.species.tree, fsize=0.6)

#using chronopl to make it ultrametric
# to reduce ANGst from over estimating the HGT events!
angst.species.tree.ultra <- ape::chronos(angst.species.tree,lambda=0)
#angst.species.tree.ultra <- ape::root.phylo(angst.species.tree.ultra, outgroup = gsub('_','',Marinobacter_A))

#write downt the file
ape::write.tree(angst.species.tree.ultra, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile' )

is.ultrametric(read.tree('~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile'))



#ape::write.tree(angst.species.tree, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile' )



angst.species.tree$tip.label

angst.species.tree$tip.label <- 1:length(angst.species.tree$tip.label)
ape::write.tree(angst.species.tree, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile' )

is.rooted(read.tree('~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile'))



# THIS IS NOW JUST BECAUSE I WANT TO TEST IT WITH THE SCO
# IDEALLY, WE GENERATE FASTA FILES FOR ALIGNMENTS, THAT HAVE APPROPORATE ANGST FORMATTED HEADERS!
# BEING GSUBBED underscores, and additional numbering _1, _2 per genome
#in gene trees, need a identifier, _1 etc
angst.loci.trees <- list()
#class(angst.loci.trees) = 'multiPhylo'

for(gene.tree.name in names(sco.loci.trees)){
  gene.tree <- sco.loci.trees[[gene.tree.name]]
  #print(gene.tree)
  gene.tree$tip.label <- paste0(gsub('_','',gene.tree$tip.label),'_1')
  angst.loci.trees[[gene.tree.name]] <- gene.tree
}
class(angst.loci.trees) = 'multiPhylo'
ape::write.tree(angst.loci.trees, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.genes.treefile' )

ape::write.tree(angst.loci.trees[[1]], '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.gene.treefile' )



gene.tree <- sco.loci.trees[[5]]
gene.tree <- rib.concat.tree



#gene.tree$tip.label <- 1:length(gene.tree$tip.label)
gene.tree$tip.label <- paste0(gsub('_','',gene.tree$tip.label),'_1')
gene.tree$edge.length <- gene.tree$edge.length + 0.0001
ape::write.tree(gene.tree, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.gene.treefile' )

gene.tree <- sco.loci.trees[[5]]






gene.tree <- rtree(10)
gene.tree$tip.label <- 1:length(gene.tree$tip.label)
gene.tree$tip.label <- paste0(gsub('_','',gene.tree$tip.label),'_1')
angst.species.tree <- ape::root.phylo(angst.species.tree, outgroup = gsub('_','',Marinobacter_A))

ape::write.tree(gene.tree, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.gene.treefile' )






angst.species.tree <- rtree(10)

angst.species.tree$tip.label <- 1:length(angst.species.tree$tip.label)
angst.species.tree <- ape::root.phylo(angst.species.tree, outgroup = 10)
angst.species.tree <- ape::root(angst.species.tree, outgroup = 10)
angst.species.tree$edge.length <- angst.species.tree$edge.length + 0.0001

ape::write.tree(angst.species.tree, '~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/angst.species.treefile' )




#==========



#PopCOGent
popcognet <- read.table('~/DATA/MarbGenomics/summary.res.txt')
popcognet <- read.table('~/DATA/MarbGenomics/summary_out.txt')

sel.genome = sco.94.concat.tree$tip.label

unique(popcognet$V1 )

popcognet$V1 = gsub('-','_',gsub('\\.','_', gsub('.tsv','',popcognet$V1)))
popcognet$V2 = gsub('-','_',gsub('\\.','_', gsub('.tsv','',popcognet$V2)))

popcognet <- popcognet %>% filter(V1 %in% sel.genome) %>% filter(V2 %in% sel.genome)
lenght.bias_matrix <- popcognet %>%
  dplyr::select(V1,V2,V7) %>%
  tidyr::pivot_wider(names_from=V2, values_from=V7) %>%
  data.frame()

rownames(lenght.bias_matrix) <- lenght.bias_matrix$V1

lenght.bias_matrix <- lenght.bias_matrix[,2:ncol(lenght.bias_matrix)]

head(lenght.bias_matrix)

df2 <- data.frame(cbind(
      'ANIb' = c(as.matrix(ANIb[rownames(lenght.bias_matrix),colnames(lenght.bias_matrix)])),
      'AAI' = c(as.matrix(AAI[rownames(lenght.bias_matrix),colnames(lenght.bias_matrix)])),
      'length.bias' = c(as.matrix(lenght.bias_matrix[rownames(lenght.bias_matrix),colnames(lenght.bias_matrix)])),
      'GRR' = c(as.matrix(pan.ggr[rownames(lenght.bias_matrix),colnames(lenght.bias_matrix)])),
      'COPH' = c(as.matrix(COPH[rownames(lenght.bias_matrix),colnames(lenght.bias_matrix)])),
      'nj.grr' = c(as.matrix(cophenetic(nj.GGR)[rownames(lenght.bias_matrix),colnames(lenght.bias_matrix)]))
))

df2$ANIb <- as.numeric(as.character(df2$ANIb))
df2$AAI <- as.numeric(as.character(df2$AAI))
df2$length.bias <- as.numeric(as.character(df2$length.bias))


popcognet <- popcognet %>%
  dplyr::mutate(comparison = paste(V1,V2,sep='_'))

ani.long <- ANIb %>%
  tibble::rownames_to_column(var='V2') %>%
  tidyr::pivot_longer(cols=-V2, names_to='V1') %>%
  dplyr::mutate(comparison=paste(V1,V2,sep='_')) %>%
  dplyr::mutate(ANIb = value) %>%
  dplyr::select(comparison,ANIb)

geo.long <- geo.dist %>%
  data.frame() %>%
  tibble::rownames_to_column(var='V2') %>%
  tidyr::pivot_longer(cols=-V2, names_to='V1') %>%
  dplyr::mutate(comparison=paste(V1,V2,sep='_')) %>%
  dplyr::mutate(geodist = value) %>%
  dplyr::select(comparison,geodist)


popcognet <- popcognet %>%
  dplyr::left_join(ani.long, by='comparison') %>%
  dplyr::left_join(geo.long, by='comparison')



p.pop1 <- popcognet %>%
  filter(grepl('Marinobacter',V1)) %>%
  filter(grepl('Marinobacter',V2)) %>%
  mutate(length.bias=V7)%>%
  dplyr::filter(!is.na(length.bias)) %>%
  ggplot2::ggplot(aes(x=ANIb, y=length.bias,label)) +
  ggplot2::geom_point() +
  fdb_style(aspect.ratio=0.5) +
  ggrepel::geom_text_repel(data = . %>% filter(length.bias>80), aes(label = comparison),size=2) +
  ggrepel::geom_text_repel(data = . %>% filter(length.bias>20) %>% filter(ANIb<0.8), aes(label = comparison),size=2)

ggplot2::ggsave(filename=paste0(outdir,'length.bias_v_ANIb','.pdf'),plot=p.pop1, width = 8, height = 5,unit='in')


p.pop1 <-popcognet %>%
  filter(grepl('Marinobacter',V1)) %>%
  filter(grepl('Marinobacter',V2)) %>%
  mutate(length.bias=V7)%>%
  filter(ANIb>0.95)%>%
  dplyr::filter(!is.na(length.bias)) %>%
  ggplot2::ggplot(aes(x=ANIb, y=length.bias,label)) +
  ggplot2::geom_point(aes(color=geodist)) +
  viridis::scale_color_viridis()+
  fdb_style(aspect.ratio=0.5) +
  ggrepel::geom_text_repel(data = . %>% filter(length.bias>80), aes(label = comparison),size=2) +
  ggrepel::geom_text_repel(data = . %>% filter(ANIb>0.975), aes(label = comparison),size=2)

ggplot2::ggsave(filename=paste0(outdir,'length.bias_v_ANIb_zoomed','.pdf'),plot=p.pop1, width = 8, height = 5,unit='in')




p.pop2 <-popcognet %>%
  filter(grepl('Marinobacter',V1)) %>%
  filter(grepl('Marinobacter',V2)) %>%
  mutate(length.bias=V7)%>%
  dplyr::filter(!is.na(length.bias)) %>%
  ggplot2::ggplot(aes(x=geodist/1000, y=length.bias)) +
  ggplot2::xlab('Geographical distance (km)')+
  ggplot2::geom_point(aes(fill=ANIb),shape=21,size=2) +
  viridis::scale_fill_viridis(option='inferno',discrete = FALSE)+
  fdb_style()

p.pop3 <- popcognet %>%
  filter(grepl('Marinobacter',V1)) %>%
  filter(grepl('Marinobacter',V2)) %>%
  mutate(length.bias=V7)%>%
  dplyr::filter(!is.na(geodist)) %>%
  ggplot2::ggplot(aes(x=ANIb, y=length.bias)) +
  #ggplot2::xlab('Geographical distance (km)')+
  ggplot2::geom_point(aes(fill=geodist),shape=21,size=2) +
  viridis::scale_fill_viridis(option='viridis',discrete = FALSE)+
  fdb_style()

p.pop4 <- popcognet %>%
  filter(grepl('Marinobacter',V1)) %>%
  filter(grepl('Marinobacter',V2)) %>%
  mutate(length.bias=V7)%>%
  dplyr::filter(!is.na(length.bias)) %>%
  ggplot2::ggplot(aes(x=geodist/1000, y=ANIb)) +
  ggplot2::xlab('Geographical distance (km)')+
  ggplot2::geom_point(aes(fill=length.bias),shape=21,size=2) +
  viridis::scale_fill_viridis(option='inferno',discrete = FALSE)+
  fdb_style()

p.comb <- gridExtra::grid.arrange(p.pop2,p.pop3,p.pop4, nrow=1)
ggplot2::ggsave(filename=paste0(outdir,'lenghtbias_ani_geodist','.pdf'),plot=p.comb, width = 10, height = 4,unit='in')








modeldat <- df2 %>%
  dplyr::filter(ANIb>0.95) %>%
  dplyr::filter(!is.na(length.bias)) %>%
  dplyr::filter(ANIb<1)

exp.model <-lm(length.bias ~ exp(ANIb),modeldat )
exp.model.df <- data.frame(x = modeldat$ANIb,
                           y = exp(fitted(exp.model)))

modeldat_below <- df2 %>%
  dplyr::filter(ANIb<0.95) %>%
  dplyr::filter(!is.na(length.bias))



##############=== insets
#plot different cutoffs

for(cutoff in c(.8,.85,.9,.95,.98)){
  p.1 <- df2 %>% filter(ANIb>cutoff) %>%
    dplyr::filter(!is.na(length.bias)) %>%
    dplyr::filter(ANIb<1)%>%
    ggplot2::ggplot(aes(x=ANIb, y=length.bias)) +
    ggplot2::geom_point(color='grey',shape=21) +
    ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),color='forestgreen',se=F) +
    #ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),method='lm') +
    fdb_style(aspect.ratio=0.5)
  ggplot2::ggsave(filename=paste0(outdir,'lenght_bias_v_ani_inset',cutoff,'.pdf'),plot=p.1, width = 4, height = 4,unit='in')


  p.2 <- df2 %>% filter(ANIb>cutoff) %>%
    dplyr::filter(!is.na(length.bias)) %>%
    dplyr::filter(ANIb<1)%>%
    ggplot2::ggplot(aes(x=COPH, y=length.bias)) +
    ggplot2::geom_point(color='grey',shape=21) +
    ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=length.bias),color='forestgreen',se=F) +
    #ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),method='lm') +
    fdb_style(aspect.ratio=0.5)+ ylab('Length bias')
  ggplot2::ggsave(filename=paste0(outdir,'lenght_bias_v_coph_inset',cutoff,'.pdf'),plot=p.2, width = 4, height = 4,unit='in')


  p.3 <- df2 %>% filter(ANIb>cutoff) %>%
    dplyr::filter(!is.na(length.bias)) %>%
    dplyr::filter(ANIb<1)%>%
    ggplot2::ggplot(aes(x=COPH, y=nj.grr)) +
    ggplot2::geom_point(color='grey',shape=21) +
    ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=nj.grr),color='forestgreen',se=F) +
    #ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),method='lm') +
    fdb_style(aspect.ratio=0.5)+xlab('COPH') + ylab('gene-content tree distance')
  ggplot2::ggsave(filename=paste0(outdir,'geneconetent_v_species_inset',cutoff,'.pdf'),plot=p.3, width = 4, height = 4,unit='in')

  p.4 <- gridExtra::grid.arrange(p.1,p.2,p.3,ncol=1)
  ggplot2::ggsave(filename=paste0(outdir,'lb_comb_inset',cutoff,'.pdf'),plot=p.4, width = 4, height = 5,unit='in')
}




df2 %>%
  ggplot2::ggplot(aes(x=COPH, y=GRR)) +
  ggplot2::geom_point(color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=GRR),color='forestgreen',se=F) +

  #ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),method='lm') +
  fdb_style(aspect.ratio=1)+ ylab('Length bias')





####==================================#
#### FIGURE 2 ####
###===================================#


a.diff <- max(df2$GRR) - min(df2$GRR)
b.diff <-max(df2[!is.na(df2$length.bias),]$length.bias) - min(df2[!is.na(df2$length.bias),]$length.bias)
a.min <- min(df2$GRR)
b.min <- min(df2[!is.na(df2$length.bias),]$length.bias)


p1 <- df2 %>%
  ggplot2::ggplot(aes(x=COPH, y=GRR)) +
  ggplot2::geom_point(aes(x=COPH, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=ANIb)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=GRR),color='black',se=F,size=1.5) +
  #  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=GRR),method='lm',color='black',se=F,size=0.5) +
  ggplot2::geom_smooth(data=modeldat_below,aes(x=COPH, y=GRR),method='lm',color='red',se=F,size=0.5) +
  ggplot2::geom_smooth(data=modeldat_below,aes(x=COPH, y=GRR),color='red',se=F,linetype='dashed',size=0.5) +
  #ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),method='lm') +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                                             name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(limits=c(0.7,1))



p2 <- df2 %>% filter(ANIb>0.95) %>%
  ggplot2::ggplot(aes(x=COPH, y=GRR)) +
  ggplot2::geom_point(aes(x=COPH, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=ANIb)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=GRR),color='black',se=F,size=1.5) +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                          name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(limits=c(0.7,1))

ggplot2::ggsave(filename=paste0(outdir,'allrangeGRR_Coph','.pdf'),plot=p1, width = 4, height = 4,unit='in')
ggplot2::ggsave(filename=paste0(outdir,'inset_GRR_Coph','.pdf'),plot=p2, width = 4, height = 4,unit='in')


#make inset
p3 <- cowplot::ggdraw() +
  cowplot::draw_plot(p1) +
  cowplot::draw_plot(x=0.35,y=0.4, width=0.65,height=0.65)

ggplot2::ggsave(filename=paste0(outdir,'inclusive_inset_GRR_Coph_fr','.pdf'),plot=p3, width = 5, height = 5,unit='in')


p4 <- cowplot::ggdraw() +
  cowplot::draw_plot(p1) +
  cowplot::draw_plot(p2+fdb_style(aspect.ratio=0.6)+theme(legend.position = 'none'),x=0.35,y=0.4, width=0.5,height=0.5)

ggplot2::ggsave(filename=paste0(outdir,'inclusive_inset_GRR_Coph_fr_aspect_ratio_0.5nl','.pdf'),plot=p4, width = 5, height = 5,unit='in')


#--------------#
# vs ANI

p1 <- df2 %>%
  ggplot2::ggplot(aes(x=ANIb, y=GRR)) +
  ggplot2::geom_point(aes(x=ANIb, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=GRR),color='black',se=F,size=1.5) +
  #  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=GRR),method='lm',color='black',se=F,size=0.5) +
  ggplot2::geom_smooth(data=modeldat_below,aes(x=ANIb, y=GRR),method='lm',color='red',se=F,size=0.5,linetype='dashed') +
  ggplot2::geom_smooth(data=modeldat_below,aes(x=ANIb, y=GRR),color='red',se=F,size=0.5) +
  #ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=length.bias),method='lm') +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                          name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))
#+ viridis::scale_color_viridis(limits=c(0.7,1))



p2 <- df2 %>% filter(ANIb>0.95) %>%
  ggplot2::ggplot(aes(x=ANIb, y=GRR)) +
  ggplot2::geom_point(aes(x=ANIb, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=ANIb, y=GRR),color='black',se=F,size=1.5) +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                          name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))

ggplot2::ggsave(filename=paste0(outdir,'allrangeGRR_ANIb','.pdf'),plot=p1, width = 4, height = 4,unit='in')
ggplot2::ggsave(filename=paste0(outdir,'inset_GRR_ANIb','.pdf'),plot=p2, width = 4, height = 4,unit='in')


#make inset
p3 <- cowplot::ggdraw() +
  cowplot::draw_plot(p1) +
  cowplot::draw_plot(p2, x=0.35,y=0.4, width=0.65,height=0.65)

ggplot2::ggsave(filename=paste0(outdir,'inclusive_inset_GRR_ANIb_fr','.pdf'),plot=p3, width = 5, height = 5,unit='in')

p4 <- cowplot::ggdraw() +
  cowplot::draw_plot(p1) +
  cowplot::draw_plot(p2+fdb_style(aspect.ratio=0.6)+theme(legend.position = 'none'),x=0.35,y=0.5, width=0.5,height=0.5)

ggplot2::ggsave(filename=paste0(outdir,'inclusive_inset_GRR_ANIb_fr_aspect_ratio_0.5nl','.pdf'),plot=p4, width = 5, height = 5,unit='in')


#------------------ <- #
#AAI

p1 <- df2 %>%
  ggplot2::ggplot(aes(x=AAI, y=GRR)) +
  ggplot2::geom_point(aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=AAI, y=GRR),color='black',se=F,size=1.5) +
  ggplot2::geom_smooth(data=modeldat_below,aes(x=AAI, y=GRR),method='lm',color='red',se=F,size=0.5,linetype='dashed') +
  ggplot2::geom_smooth(data=modeldat_below,aes(x=AAI, y=GRR),color='red',se=F,size=0.5) +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                          name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))
#+ viridis::scale_color_viridis(limits=c(0.7,1))



p2 <- df2 %>% filter(ANIb>0.95) %>%
  ggplot2::ggplot(aes(x=AAI, y=GRR)) +
  ggplot2::geom_point(aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=AAI, y=GRR),color='black',se=F,size=1.5) +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                          name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))

ggplot2::ggsave(filename=paste0(outdir,'allrangeGRR_AAI','.pdf'),plot=p1, width = 4, height = 4,unit='in')
ggplot2::ggsave(filename=paste0(outdir,'inset_GRR_AAI','.pdf'),plot=p2, width = 4, height = 4,unit='in')


#make inset
p3 <- cowplot::ggdraw() +
  cowplot::draw_plot(p1) +
  cowplot::draw_plot(p2, x=0.35,y=0.4, width=0.65,height=0.65)

ggplot2::ggsave(filename=paste0(outdir,'inclusive_inset_GRR_AAI_fr','.pdf'),plot=p3, width = 5, height = 5,unit='in')

p4 <- cowplot::ggdraw() +
  cowplot::draw_plot(p1) +
  cowplot::draw_plot(p2+fdb_style(aspect.ratio=0.6)+theme(legend.position = 'none'),x=0.35,y=0.5, width=0.5,height=0.5)

ggplot2::ggsave(filename=paste0(outdir,'inclusive_inset_GRR_AAI_fr_aspect_ratio_0.5nl','.pdf'),plot=p4, width = 5, height = 5,unit='in')




df2 %>%
  ggplot2::ggplot(aes(x=AAI, y=COPH)) +
  #ggplot2::geom_point(aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(aes(x=AAI, y=COPH),color='grey',se=F) +
  ggplot2::geom_smooth(data=modeldat,aes(x=AAI, y=COPH),color='black',se=F,size=1.5) +
  fdb_style(aspect.ratio=1) +
  viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))


df2 %>%
  ggplot2::ggplot(aes(x=AAI, y=ANIb)) +
  #ggplot2::geom_point(aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(color='grey',se=F) +
  ggplot2::geom_smooth(data= df2 %>% filter(ANIb>0.95), color='black',se=F,size=1.5) +
  fdb_style(aspect.ratio=1) +
  viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))


df2 %>% tidyr::pivot_longer(cols=-COPH) %>%
  ggplot2::ggplot(aes(x=COPH, y=value)) +
  #ggplot2::geom_point(aes(x=AAI, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(color='grey',se=F) +
  #ggplot2::geom_smooth(data= tidyr::pivot_longer(cols=-COPH) %>% filter(ANIb>0.95), color='black',se=F,size=1.5) +
  fdb_style(aspect.ratio=1) +
  viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))+facet_wrap(~name,scales='free_y')



popcognet <-
  ggplot2::ggplot(aes(x=COPH, y=value)) +
  ggplot2::geom_point(shape=21,aes(color=COPH)) +
  ggplot2::geom_smooth(color='grey',se=F) +
  fdb_style(aspect.ratio=1) +
  viridis::scale_color_viridis(direction = -1,limits=c(0,0.8))+facet_wrap(~name,scales='free_y')







modeldat2 <- df2 %>%
  dplyr::filter(ANIb<0.95) %>%
  dplyr::filter(ANIb>0.85) %>%
  dplyr::filter(!is.na(length.bias)) %>%
  dplyr::filter(ANIb<1)
modeldat3 <- df2 %>%
  dplyr::filter(ANIb<0.85) %>%
  dplyr::filter(ANIb>0.8) %>%
  dplyr::filter(!is.na(length.bias)) %>%
  dplyr::filter(ANIb<1)
modeldat4 <- df2 %>%
  dplyr::filter(ANIb<0.8) %>%
  dplyr::filter(!is.na(length.bias)) %>%
  dplyr::filter(ANIb<1)

plot.p <- df2 %>%
  dplyr::filter(!is.na(length.bias)) %>%
  dplyr::filter(ANIb<1)%>%
  tidyr::pivot_longer(cols=c(ANIb, length.bias, GRR, nj.grr)) %>%
  ggplot2::ggplot(aes(x=COPH, y=value)) +
  ggplot2::geom_point(color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat %>%  tidyr::pivot_longer(cols=c(ANIb, length.bias, GRR, nj.grr)) ,
                       aes(x=COPH, y=value),color='red',se=F) +
  ggplot2::geom_smooth(data=modeldat2 %>%  tidyr::pivot_longer(cols=c(ANIb, length.bias, GRR, nj.grr)) ,
                       aes(x=COPH, y=value),color='blue',se=F,method='lm') +
  ggplot2::geom_smooth(data=modeldat3 %>%  tidyr::pivot_longer(cols=c(ANIb, length.bias, GRR, nj.grr)) ,
                       aes(x=COPH, y=value),color='forestgreen',se=F,method='lm') +
  ggplot2::geom_smooth(data=modeldat4 %>%  tidyr::pivot_longer(cols=c(ANIb, length.bias, GRR, nj.grr)) ,
                       aes(x=COPH, y=value),color='purple',se=F,method='lm') +
  fdb_style(aspect.ratio=0.5)+xlab('species tree distance') +facet_wrap(~name,ncol=1,scales='free')
ggplot2::ggsave(filename=paste0(outdir,'recombination_at_COPH','.pdf'),plot=plot.p, width = 4, height = 10,unit='in')


plot.p2 <- df2 %>%
  dplyr::filter(!is.na(length.bias)) %>%
  dplyr::filter(ANIb<1)%>%
  tidyr::pivot_longer(cols=c(COPH, length.bias, GRR, nj.grr)) %>%
  ggplot2::ggplot(aes(x=ANIb, y=value)) +
  ggplot2::geom_point(color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat %>%  tidyr::pivot_longer(cols=c(COPH, length.bias, GRR, nj.grr)) ,
                       aes(x=ANIb, y=value),color='red',se=F) +
  ggplot2::geom_smooth(data=modeldat2 %>%  tidyr::pivot_longer(cols=c(COPH, length.bias, GRR, nj.grr)) ,
                       aes(x=ANIb, y=value),color='blue',se=F,method='lm') +
  ggplot2::geom_smooth(data=modeldat3 %>%  tidyr::pivot_longer(cols=c(COPH, length.bias, GRR, nj.grr)) ,
                       aes(x=ANIb, y=value),color='forestgreen',se=F,method='lm') +
  ggplot2::geom_smooth(data=modeldat4 %>%  tidyr::pivot_longer(cols=c(COPH, length.bias, GRR, nj.grr)) ,
                       aes(x=ANIb, y=value),color='purple',se=F,method='lm') +
  fdb_style(aspect.ratio=0.5)+xlab('Average Nucletodie Identity') +facet_wrap(~name,ncol=2,scales='free')
ggplot2::ggsave(filename=paste0(outdir,'recombination_at_ANIb','.pdf'),plot=plot.p2, width = 4, height = 10,unit='in')







#================================================================#

p.1 <- df2 %>%
  ggplot2::ggplot(aes(x=COPH, y=ANIb)) +
  ggplot2::geom_point(shape=21,aes(color=ANIb)) +
  ggplot2::geom_smooth(aes(x=COPH, y=ANIb),color='red',se=F,linetype='dashed',size=0.5) +
  fdb_style(aspect.ratio=1)+ xlab('COPH')+ ylab('ANIb') + viridis::scale_color_viridis(limits=c(0.7,1))

p.2 <- df2 %>%
  ggplot2::ggplot(aes(x=COPH, y=AAI)) +
  ggplot2::geom_point(shape=21,aes(color=ANIb)) +
  ggplot2::geom_smooth(aes(x=COPH, y=AAI),color='red',se=F,linetype='dashed',size=0.5) +
  fdb_style(aspect.ratio=1)+ xlab('COPH') + ylab('AAI') + viridis::scale_color_viridis(limits=c(0.7,1))

p.3 <- df2 %>%
  ggplot2::ggplot(aes(x=ANIb, y=AAI)) +
  ggplot2::geom_point(shape=21,aes(color=ANIb)) +
  ggplot2::geom_smooth(aes(x=ANIb, y=AAI),color='red',se=F,linetype='dashed',size=0.5) +
  fdb_style(aspect.ratio=1)+ xlab('ANIb') + ylab('AAI') + viridis::scale_color_viridis(limits=c(0.7,1))

p.4 <- gridExtra::grid.arrange(p.1,p.2,p.3,ncol=3)
ggplot2::ggsave(filename=paste0(outdir,'/S8','.pdf'),plot=p.4, width = 10, height = 6,unit='in')




#================================================================#
# df2 %>%
#   ggplot2::ggplot(aes(x=1-ANIb, y=GRR)) +
#   ggplot2::geom_point(shape=21,aes(color=ANIb)) +
#   ggplot2::geom_smooth(data=. %>% filter(ANIb>=0.80),aes(x=1-ANIb, y=GRR),color='black',se=F,size=1.5) +
#   ggplot2::geom_smooth(data=. %>% filter(ANIb>=0.80),method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F) +
#   ggplot2::geom_smooth(data=. %>% filter(ANIb>=0.80),method='lm',se=F,linetype='dashed') +
#   ggplot2::geom_smooth(data=modeldat_below %>% filter(ANIb<0.80),aes(x=1-ANIb, y=GRR),method='lm',color='red',se=F,size=0.5) +
#   ggplot2::geom_smooth(data=modeldat_below %>% filter(ANIb<0.80),aes(x=1-ANIb, y=GRR),color='red',se=F,linetype='dashed',size=0.5) +
#   fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(limits=c(0.7,1))
#
# df2 %>%
#   ggplot2::ggplot(aes(x=1-ANIb, y=GRR)) +
#   ggplot2::geom_point(shape=21,aes(color=ANIb)) +
#   ggplot2::geom_smooth(data=. %>% filter(ANIb>=0.80),method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F) +
#   ggplot2::geom_smooth(data=modeldat_below %>% filter(ANIb<0.80),aes(x=1-ANIb, y=GRR),method='lm',color='red',se=F,size=0.5) +
#   ggplot2::geom_smooth(data=modeldat_below %>% filter(ANIb<0.80),aes(x=1-ANIb, y=GRR),color='red',se=F,linetype='dashed',size=0.5) +
#   fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(limits=c(0.7,1))

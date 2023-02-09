#----------------------------------------
#
#
#
#-----------------------------------------


pan = OG.pan
data = PfamData
data = CAZYl
data = t(domains.wide)
data = MAPLE

data.frame(pan[intersect(metadata$genome,rownames(pan)),])

pan = pan[intersect(metadata$genome,rownames(pan)),]

metadf = metadf[rownames(PRM_data),]


dPRM_data = distManhattan(PRM_data)
#adonis2(dPRM_data ~ group2 * SS, data = metadata[rownames(PRM_data),])

#SIMPER ANALYSIS
simper.group = simper(data.frame(pan), metadf$group, permutations=100)

adonis_phylogeny = vegan::adonis(dPRM_data ~ group, data = metadata[rownames(PRM_data),])
adonis_phylogeny

anosim_phylogeny = vegan::anosim(dPRM_data, metadata[rownames(PRM_data),]$group)
anosim_phylogeny # take a look at results


dispersion<-vegan::betadisper(dPRM_data, group= metadata[rownames(PRM_data),]$group)
vegan::permutest(dispersion, pairwise=T, permutation=1000)
PCoA1=dispersion$vectors[,1]
PCoA2=dispersion$vectors[,2]
PCoA.plot<-cbind(PCoA1, PCoA2, metadata[rownames(PRM_data),])

phylo.pcoa.p <-ggplot(PCoA.plot, aes(PCoA1, PCoA2))+
  geom_point(position=position_jitter(.1),aes(color=group))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill= group))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+ theme(legend.position = "none")
phylo.pcoa.p

average_dispersion_from_centroid = cbind(data.frame(dispersion$distances),metadata)
phylo.disp.p = ggplot(average_dispersion_from_centroid,aes(x= group2, y=dispersion.distances,fill= group2, color= group2))+geom_boxplot(alpha=0.7)+geom_point()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "none")




SIMPER = vegan::simper(PRM_data, metadata$group, permutations=100)








phyloMDS<-metaMDS(dPRM_data, k=2, trymax=35, autotransform=TRUE) ##k is the number of dimensions
phyloMDS ##metaMDS takes eaither a distance matrix or your community matrix (then requires method for 'distance=')
stressplot(phyloMDS)

NMDS1 <- phyloMDS$points[,1] ##also found using scores(phyloMDS)
NMDS2 <- phyloMDS$points[,2]

mds.plot<-cbind(NMDS1, NMDS2, metadata[rownames(PRM_data),])

phylo.p<-ggplot(mds.plot, aes(NMDS1, NMDS2))+
  geom_point(position=position_jitter(.1),aes(color=group2))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill= group2))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+ theme(legend.position = "none")
phylo.p

multiplot(phylo.p, phylo.pcoa.p ,Habitat.p, Habitat.pcoa.p ,cols=2)
adonis_phylogeny
adonis_habitat


multiplot(phylo.disp.p, habitat.disp.p,cols=2)



#-------
if(!require(fpc)) install.packages("fpc", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(factoextra)) install.packages("factoextra", repos = "http://cran.us.r-project.org",dependencies=TRUE)


#-----------------------------------------------------------------------------
#
# 	ProkCompR 	UTILITY FUNCTIONS
#
#-----------------------------------------------------------------------------
#	Loading phylogenetic tree

loadTree <- function(inPath='', clean_genome_names=TRUE){
	if(inPath==''){
		print("!Error: You'd need to specify a directory using inPath")
	}else{
		tree = ape::read.tree(inPath)
		if(clean_genome_names==TRUE){
			tree$tip.label = gsub('-','_',gsub('\\.','_',tree$tip.label))
		}
	return(tree)
	}
}

tree = loadTree(inPath='~/DATA/MarbGenomics/SCO_genefpair_trimAl/concat.treefile', clean_genome_names=TRUE)
#Annoying, but this below codeline is needed for the Marinobacter comparative genomics study
tree$tip.label[tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"



#-----------------------------------------------
excluded <- genome.tbl %>% filter(!(genome %in% tree$tip.label))
excluded %>% dplyr::select(genome, quality_class,manualy_removed)

metadata <- genome.tbl %>% filter(genome %in% tree$tip.label)

metadata %>% dplyr::count(genus)
genome.tbl %>% filter(manualy_removed == FALSE) %>% dplyr::select(genome, quality_class,manualy_removed)
genome.tbl %>% filter(manualy_removed == FALSE)  %>% filter(!(genome %in% tree$tip.label)) %>% dplyr::select(genome, quality_class,manualy_removed)


#------- build phylogroup
#	Inspired by Levy et al 2017
#

#optimalNclust

COPH = cophenetic(tree)




library(cluster)
library(factoextra)


COPHmb = COPH[grepl("Marinobacter|Tamil",rownames(COPH)),grepl("Marinobacter|Tamil",rownames(COPH))]


fviz_nbclust(COPHmb, kmeans, method = "silhouette",k.max=25)+
  labs(subtitle = "Silhouette method")

#Either 9 ot 20, peak at arround 0.6
#interesting, what should one do here?

pam_out = fpc::pamk(COPHmb,krange=18)
#pam_out = pamk(COPHmb,krange=9)




pam_out$pamobject$clustering

pam_cl = data.frame(pam_out$pamobject$clustering)
colnames(pam_cl) = c('cluster')




#-----------------------------------------------


subtree = drop.tip(tree,tree$tip.label[!grepl('Marinobacter|Tamil', tree$tip.label)])

p <- ggtree(subtree, layout = "circular", open.angle = 50,branch.length="none", size=0.5)
p2 <- open_tree(p, 180)

p <- ggtree(subtree, size=0.5)
p2 <-p


gheatmap(p2, pam_cl,offset=-1, width=0.2)+
  scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdBu')),name='PAM cluster')+
  theme(
		aspect.ratio = 1,
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
        #legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.key.width = unit(0.25, 'cm'),
        legend.key.height = unit(0.55, 'cm'),
        legend.text = element_text(size = 10),
        #legend.text.align = 1,
        legend.title = element_text(size = 11))



MDS = cmdscale(subDistMan, eig = TRUE)
MDS = cmdscale(subDistJac, eig = TRUE)

MDSdata = data.frame(MDS$points)
#MDSdata$nrofgenes =  str_count(rownames(MDSdata) , pattern = "_")
#MDSdata$nrofgenes =  MDSdata$nrofgenes + 1
MDSdata$Name = rownames(MDSdata)
MDSdata$Genus = gsub( "_.*$", "", rownames(MDSdata) )
colnames(MDSdata) = c("x",'y','Name',"Genus")
MDSdata$PAMcluster = pam_cl[gsub('-','_',gsub('\\.','_', rownames(MDSdata))),]

ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(PAMcluster),label=Name)) +
	theme_classic() +
	xlab('Dimension 1') +
	ylab('Dimension 2') +
	labs(fill = "PAM cluster") +
	ggtitle(label='Manhattan distance')+
	geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_point(shape = 21, size = 2) +
	#scale_fill_brewer(palette = "Set1",direction=-1)+
	theme(
		aspect.ratio = 1,
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(size = 11, colour = '#000000'),
		axis.text = element_text(size = 10, colour = '#000000'),
		#legend.position = c(1, 1),
		legend.justification = c(1, 1),
		legend.key.width = unit(0.25, 'cm'),
		legend.key.height = unit(0.55, 'cm'),
		legend.text = element_text(size = 10),
		#legend.text.align = 1,
		legend.title = element_text(size = 11))


ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(PAMcluster),label=Name)) +
	theme_classic() +
	xlab('Dimension 1') +
	ylab('Dimension 2') +
	labs(fill = "PAM cluster") +
	ggtitle(label='Manhattan distance')+
	geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_point(shape = 21, size = 2) +
	stat_ellipse(aes(color= as.factor(PAMcluster))) +
	#scale_fill_brewer(palette = "BuBG",direction=-1)+
	theme(
		aspect.ratio = 1,
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(size = 11, colour = '#000000'),
		axis.text = element_text(size = 10, colour = '#000000'),
		#legend.position = c(1, 1),
		legend.justification = c(1, 1),
		legend.key.width = unit(0.25, 'cm'),
		legend.key.height = unit(0.55, 'cm'),
		legend.text = element_text(size = 10),
		#legend.text.align = 1,
		legend.title = element_text(size = 11))



CheckM[rownames(pam_cl),]

rownames(pam_cl)
pam_cl[]

#-----


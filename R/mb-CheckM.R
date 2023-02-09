#------------------------------------------------------------#
#	  Generates genome.tbl
#   Identify sufficient quality genomes based on CheckM
#	    Marinobacter comparative genomics
#	    Frederik De Boever, 26-nov-2020
#------------------------------------------------------------#

# loadCheckM
#=================================================================================#
#' Title
#'
#' @param inPath
#' @param clean_genome_names
#'
#' @return
#' @export
#'
#' @examples

loadCheckM <- function(inPath='', clean_genome_names=TRUE){
	if(inPath==''){
		print("You'd need to specify a directory using inPath")
	}else{
		CheckM = read.table(inPath,header=FALSE,sep="\t",skip=1)

		colnames(CheckM) = c(
			'BinId',
			'Marker_lineage',
			'genomes',
			'markers',
			'marker_sets',
			'Completeness',
			'Contamination',
			'Strain_heterogeneity',
			'Genome_size',
			'ambiguous_bases',
			'scaffolds',
			'contigs',
			'N50_scaffolds',
			'N50_contigs',
			'Mean_scaffold_length_bp',
			'Mean_contig_length_bp',
			'Longest_scaffold_bp',
			'Longest_contig_bp',
			'GC',
			'GC_std',
			'Coding_density',
			'Translation_table',
			'predicted_genes',
			'c0',
			'c1',
			'c2',
			'c3',
			'c4',
			'c5plus'
			)

		if(clean_genome_names==TRUE){
			CheckM$BinId = gsub('-','_',gsub('\\.','_', CheckM$BinId))
		}

		rownames(CheckM) = CheckM$BinId
		CheckM$genome <- rownames(CheckM)

		return(CheckM)
	}
}

# loadAdditionalData
#=================================================================================#
#' function to load additional data
#'
#' @param inPath
#' @param name_col select the column to be used as genome identifier
#' @param clean_genome_names
#'
#' @return dataframe with genome identifiers as rownames
#' @export
#'
#' @examples
#' mergedEnv <- loadAdditionalData(inPath='~/DATA/MarbGenomics/mergedEnv.txt',name_col='genome_name_dotted',clean_genome_names=TRUE)

loadAdditionalData <- function(inPath='', name_col ,clean_genome_names=TRUE){
  if(inPath==''){
    print("You'd need to specify a directory using inPath")
  }else{
    meta.df = read.delim(inPath)
    if(clean_genome_names==TRUE){
      meta.df[,name_col] = gsub('-','_',gsub('\\.','_', meta.df[,name_col]))
    }
    rownames(meta.df) <- meta.df[,name_col]
    meta.df$genome <- rownames(meta.df)
    meta.df <- meta.df %>%
      dplyr::select(-as.name(name_col))
    return(meta.df)
  }
}


###RUN CODE####
#=================================================================================#
# LOAD DATA WITH LOADING FUNCTIONS
#=================================================================================#

#Load CheckM data
CheckM = loadCheckM(
  inPath= '~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/checkm_out/qa2.txt',
  clean_genome_names=TRUE)

#Load metadata
mergedEnv <- loadAdditionalData(
  inPath='~/DATA/MarbGenomics/mergedEnv.txt',
  name_col='genome_name_dotted',
  clean_genome_names=TRUE)

#specific to Marinobacter case, select these variables
mergedEnv <- mergedEnv %>%
  dplyr::select(genome,
                Accession,
                Organism,
                strain,
                TypeStrain,
                lifestyle,
                SS,
                sourceClass,
                lat,
                lon,
                geo_loc_name,
                host,
                isolation_source,
                origin.Continent,
                colony_morphology.Colony.shape,
                colony_morphology.Colony.color
  )


#=================================================================================#
### Generate Genome tibble ####
# makeGenomeTable
#=================================================================================#

#merge loaded checkm and metadata tables into genome tibble
genome.tbl <- CheckM %>%
    dplyr::select(-BinId) %>%
    dplyr::left_join(mergedEnv, by = 'genome') %>%
    dplyr::select(genome, Accession, Organism, strain, TypeStrain, everything()) %>%
    tibble::as_tibble()

#add a column that denotes the genus
# this is derived from genome identifier, so likely specific to Marinobacter, as I used the full name in initial files.
genome.tbl <- genome.tbl %>%
  dplyr::mutate(genus = gsub( "_.*$", "", genome))


### QUALITY ASSESMENT OF GENOMES ####
#quality assignment (low, medium, high)
# generates a column quality_class that specifies the classes low medium and high
genome.tbl <- genome.tbl %>% dplyr::mutate(quality_class = ifelse(Completeness >=98 &
                                                      Contamination < 5 &
                                                      N50_scaffolds > 50000, ifelse(Completeness >=99 &
                                                                                    Contamination < 1 &
                                                                                    N50_scaffolds > 50000, 'high', 'medium'), 'low')) %>%
                            dplyr::mutate(quality_class = factor(quality_class, levels = c('low', 'medium', 'high')))


#----------------------------------------#
# Manual inspection removed duplicates
# CheckM results used to determine the lowest quality genomes of a set of duplicates
# Or removed, because of bad quality
manual_removal <- c(
  'Marinobacter_sp_HK15',
  'Marinobacter_hydrocarbonoclasticus_ASM1536',
  'Marinobacter_guineae_strain_M3B',
  'Marinobacter_dagiaonensis_CGMCC_1_9167',
  'Marinobacter_sp_DSM_26291',
  'Marinobacter_confluentis_HJM_18',
  'Marinobacter_lipoliticus_BF04_CF_4',
  'Marinobacter_lutaoensis_KAZ22',
  'Marinobacter_sp_UBA4489',
  'Marinobacter_persicus_KCTC_23561',
  'Marinobacter_flavimaris_LMG_23834',
  'Marinobacter_flavimaris',
  'Marinobacter_salarius_N1',
  'Marinobacter_salexigens_strain_HJR7',
  'Oleiphilus_sp__HI0009',
  'Marinobacter_similis__A3d10',
  'Marinobacter_sp_PWS21'
)


#Optional, store this information into the genome.tbl, retaining the genomes but allowing filtering
genome.tbl <- genome.tbl %>%
  dplyr::mutate(manualy_removed = genome %in% manual_removal)




#=================================================================================#
#Rest of the script is optional, and exploratory, patterns in quality?
#------#




#Annotate the phylogenetic cluster (need to move functions here I'd assume)

#pam_cl <- pam_cl %>% rownames_to_column(var='genome') %>% mutate(group=paste0('cl',cluster))
#genome.tbl <- genome.tbl %>% dplyr::left_join(pam_cl, by='genome')



#======


ProtAna=  read.delim('~/DATA/MarbGenomics/pI_final_shrunk.csv',sep=' ',header=TRUE)
pI = ProtAna
colnames(pI) = c('nr','genome','locus_tag','pI','molecular_weigth')
head(pI)

#if(clean_genome_names==TRUE){
pI$genome = gsub('-','_',gsub('\\.','_', gsub('\\.faa','', pI$genome)))
#}


pI$acidic <- ifelse(pI$pI > 7.5 ,"basic", "acidic")

g1 = ggplot(pI, aes(x= pI, color= "black", fill= "black")) +
  geom_histogram(aes(y=..density..), position="identity", alpha=1,binwidth = 0.1)+ scale_y_continuous(expand = c(0, 0))+scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black"))+
  labs(title= strBias ,x="class of pI", y = "Density")+
  theme_classic()+theme(legend.position="none")

##################################################
# plot visialisation per proteome


myplots <- list()
myplots2 <- list()
geno= c()
bias = c()

i = 1
for (organism in unique(pI$genome)){
  pIsub =subset(pI, genome ==organism)
  Nbasic = length(which(pIsub$acidic == 'basic'))
  Nacidic = length(which(pIsub$acidic == 'acidic'))
  pIBias = (Nbasic-Nacidic)/(Nbasic+Nacidic)*100
  strBias = paste(organism,":",round(pIBias,2),"%",sep=" ")

  bias = c(bias, pIBias)
  geno = c(geno, organism)
  p1 <- ggplot(pIsub, aes(x= pI, color= "black", fill= "black")) + geom_histogram(aes(y=..density..), position="identity", alpha=1,binwidth = 0.1) + scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black"))+ labs(title= strBias, x="class of pI", y = "Density")+ theme_classic()+theme(legend.position="none",plot.title = element_text(size=6))
  g2 = ggplot(pIsub, aes(x=pI,y= log(length),color = acidic)) + geom_point()+scale_color_manual(values=c("grey", "black")) + scale_fill_manual(values=c("grey", "black"))+ labs(title= strBias,x="pI", y = "log Length")+ theme_classic()+theme(legend.position="none",plot.title = element_text(size=6))
  myplots[[i]] <- p1
  myplots2[[i]] <- g2
  i = i + 1
}
#stores all the calculated b-values in a new data.frame
bias_data = data.frame(cbind(genome=geno,b=bias))
bias_data$b <- as.numeric(as.character(bias_data$b))

genome.tbl <- genome.tbl %>%
  dplyr::left_join(bias_data,by='genome')


#=================================================================================#
#
#=================================================================================#
# Explore
genome_all.tbl <- genome.tbl

genome.tbl <- genome.tbl %>% dplyr::filter(manualy_removed==FALSE)

#strains per typestrain, isolate, MAG within Marinobacter
genome.tbl %>% dplyr::filter(genus=='Marinobacter') %>% select(TypeStrain) %>% count()

#Howmany per class
genome.tbl %>% select(quality_class) %>% count()

#----
# WAIT TO DROP THE LOW QUALITY GENOMES, AS WE DO NOT DISCARD THEM HERE

#count the number of genomes per quality_class
genome.tbl %>% dplyr::select(quality_class) %>%
  dplyr::count()

genome.tbl %>%  dplyr::count(quality_class, genus , sort = TRUE)


#visualise some stuff
p.genome_quality <- genome.tbl %>%
  ggplot2::ggplot(ggplot2::aes(quality_class, Genome_size/1000000))+
  ggplot2::geom_boxplot(outlier.shape=NA,color='grey80')+
  ggplot2::geom_jitter(width=.2,aes(fill=genus,size=Contamination),shape=21)+
  ggplot2::ylab('Genome size (Mbp)') +
  ggplot2::xlab('Genome Quality') +
  fdb_style()+ggplot2::scale_fill_brewer(palette='Dark2')

ggsave('~/DATA/MarbGenomics/Graphs/Genome_size_genome_quality.pdf',plot=p.genome_quality, width = 4, height = 4,unit='in')

p.genome_quality <- genome.tbl %>% filter(genus=='Marinobacter') %>% #filter(quality_class!='low')%>%
  ggplot(aes(Genome_size/1000000,predicted_genes, group=quality_class))+
  geom_smooth(method='lm',aes(color=quality_class,fill=quality_class),se=FALSE,size=0.5,fullrange=TRUE)+
  geom_point(shape=21,size=2,aes(fill=quality_class))+
  xlab('Genome size (Mbp)') +
  ylab('nr of predicted genes')+
  ggtitle('Marinobacter only, quality')+
  scale_color_manual(values=c('firebrick','orange','forestgreen'))+
  scale_fill_manual(values=c('firebrick','orange','forestgreen'))+
  fdb_style()

ggsave('~/DATA/MarbGenomics/Graphs/qenome_quality_Marinobacter_only.pdf',plot=p.genome_quality, width = 4, height = 4,unit='in')

genome.tbl %>% filter(genus=='Marinobacter') %>%
  ggplot(aes(Genome_size/1000000,GC, group=quality_class))+
  geom_smooth(method='lm',aes(color=quality_class))+
  geom_point(shape=21,size=2,aes(fill=quality_class))+
  xlab('Genome size (Mbp)') +
  ylab('%GC-content')+
  facet_wrap(~quality_class)+
  scale_color_manual(values=c('firebrick','orange','forestgreen'))+
  scale_fill_manual(values=c('firebrick','orange','forestgreen'))+
  fdb_style()

p.genome_quality <- genome.tbl %>% filter(quality_class != 'low') %>%  filter(genus=='Marinobacter') %>%
  ggplot(aes(Genome_size/1000000,Contamination, group=quality_class))+
  geom_smooth(method='lm',aes(color=quality_class),size=0.5)+
  geom_point(shape=21,size=2,aes(fill=quality_class))+
  xlab('Genome size (Mbp)') +
  ylab('Contamination')+
  facet_wrap(~quality_class)+
  scale_color_manual(values=c('orange','forestgreen'))+
  scale_fill_manual(values=c('orange','forestgreen'))+
  ggtitle('Marinobacter only, quality')+
  fdb_style()

ggsave('~/DATA/MarbGenomics/Graphs/qenome_quality_size_v_contamination_Marinobacter_only.pdf',plot=p.genome_quality, width = 6, height = 3,unit='in')


#---------------------------------------------------------------------------------------#

# REMOVAL OF LOW QUALITY GENOMES!

#---------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------#
# OLDER VISUALISATION

#kernel density plots
g.c1 <- genome.tbl %>% filter(genus=='Marinobacter') #%>% filter(quality_class!='low') %>%
  ggplot(aes(x=GC)) +
		xlim(c(round(min(genome.tbl$GC)) - 2, round(max(genome.tbl$GC)) + 2)) +
		geom_density(fill='grey')+
		xlab('GC-content (%)') +
		ylab('Density') +
		#labs(fill = "Genus") +
		geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
	  facet_wrap(~quality_class)+fdb_style()

g.c2 <-
	ggplot(df.chkm, aes(x= Coding_density)) +
		xlim(c(round(min(df.chkm$Coding_density)) - 2, round(max(df.chkm$Coding_density)) + 2)) +
		geom_density(fill='grey')+
		xlab('Coding density (%)') +
		ylab('Density') +
		#labs(fill = "Genus") +
		geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
		theme_classic() +
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

g.c3 <-
	ggplot(df.chkm, aes(x= Genome_size/1000000)) +
		xlim(c(round(min(df.chkm$Genome_size /1000000)) - 0.5, round(max(df.chkm$Genome_size /1000000)) + 0.5)) +
		geom_density(fill='grey')+
		xlab('Genome size (bp)') +
		ylab('Density') +
		#labs(fill = "Genus") +
		geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
		theme_classic() +
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



grid.arrange(g.c3 , g.c1, g.c2, ncol=3)


#

g.c4 <-
	ggplot(df.chkm, aes(x= group, y = GC, fill=group)) +
		#xlim(c(round(min(df.chkm$GC)) - 2, round(max(df.chkm$GC)) + 2)) +
		geom_hline(yintercept = mean(df.chkm$GC), size = 0.25, colour = '#bdbdbd') +
		geom_boxplot(alpha=0.3,outlier.size=-1)+
		geom_jitter(height=0,width=0.2,size=2,shape=21,alpha=0.7) +
		ylab('GC-content (%)') +
		xlab('clade') +
		theme_classic() +
		theme(
			aspect.ratio = 1/2,
			plot.title = element_text(hjust = 0.5),
			axis.title = element_text(size = 11, colour = '#000000'),
	        axis.text = element_text(size = 10, colour = '#000000'),
	        #legend.position = c(1, 1),
			legend.justification = c(1, 1),
			legend.key.width = unit(0.25, 'cm'),
    	    legend.key.height = unit(0.55, 'cm'),
        	legend.text = element_text(size = 10),
        	#legend.text.align = 1,
        	legend.position = "none",
        	axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

g.c5 <-
	ggplot(df.chkm, aes(x= group, y = Coding_density, fill=group)) +
		#xlim(c(round(min(df.chkm$GC)) - 2, round(max(df.chkm$GC)) + 2)) +
		geom_hline(yintercept = mean(df.chkm$Coding_density), size = 0.25, colour = '#bdbdbd') +
		geom_boxplot(alpha=0.3,outlier.size=-1)+
		geom_jitter(height=0,width=0.2,size=2,shape=21,alpha=0.7) +
		ylab('Coding density (%)') +
		xlab('clade') +
		theme_classic() +
		theme(
			aspect.ratio = 1/2,
			plot.title = element_text(hjust = 0.5),
			axis.title = element_text(size = 11, colour = '#000000'),
	        axis.text = element_text(size = 10, colour = '#000000'),
	        #legend.position = c(1, 1),
			legend.justification = c(1, 1),
			legend.key.width = unit(0.25, 'cm'),
    	    legend.key.height = unit(0.55, 'cm'),
        	legend.text = element_text(size = 10),
        	#legend.text.align = 1,
        	legend.position = "none",
        	axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

g.c6 <-
	ggplot(df.chkm, aes(x= group, y = Genome_size/1000000, fill=group)) +
		#xlim(c(round(min(df.chkm$GC)) - 2, round(max(df.chkm$GC)) + 2)) +
		geom_hline(yintercept = (mean(df.chkm$Genome_size)/1000000), size = 0.25, colour = '#bdbdbd') +
		geom_boxplot(alpha=0.3,outlier.size=-1)+
		geom_jitter(height=0,width=0.2,size=2,shape=21,alpha=0.7) +
		ylab('Genome size (Mbp)') +
		xlab('clade') +
		theme_classic() +

		theme(
			aspect.ratio = 1/2,
			plot.title = element_text(hjust = 0.5),
			axis.title = element_text(size = 11, colour = '#000000'),
	        axis.text = element_text(size = 10, colour = '#000000'),
	        #legend.position = c(1, 1),
			legend.justification = c(1, 1),
			legend.key.width = unit(0.25, 'cm'),
    	    legend.key.height = unit(0.55, 'cm'),
        	legend.text = element_text(size = 10),
        	#legend.text.align = 1,
        	legend.position = "none",
        	axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

grid.arrange(g.c6, g.c4, g.c5, ncol=1)


#-----------------------------------------------------------------------------#

ggplot(df.chkm,aes(Genome_size, Coding_density,fill= group)) +
	xlab('Genome size (bp)') +
	ylab('Coding density (%)') +
	geom_point(shape = 21, size = 2) +
	stat_ellipse(aes(color=group),level=0.75) +
	theme_classic() +
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

#----------------#

g.cr1 <-
ggplot(df.chkm,aes(Genome_size/1000000, predicted_genes)) +
	xlab('Genome size (Mbp)') +
	ylab('Detected ORFs') +
	geom_point(shape = 21, size = 2,fill='grey') +
	theme_classic() +
	geom_smooth(method='lm',size=0.5,color='red')+
	ggpubr::stat_cor(method = "pearson") +
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

g.cr2 <-
ggplot(df.chkm,aes(Genome_size/1000000, GC)) +
	xlab('Genome size (Mbp)') +
	ylab('GC-content (%)') +
	geom_point(shape = 21, size = 2,fill='grey') +
	theme_classic() +
	geom_smooth(method='lm',size=0.5,color='red')+
	ggpubr::stat_cor(method = "pearson") +
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


g.cr3 <-
ggplot(df.chkm,aes(Genome_size/1000000, Coding_density)) +
	xlab('Genome size (Mbp)') +
	ylab('Coding density (%)') +
	geom_point(shape = 21, size = 2,fill='grey') +
	theme_classic() +
	geom_smooth(method='lm',size=0.5,color='red')+
	ggpubr::stat_cor(method = "pearson") +
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


grid.arrange(g.cr1, g.cr2, g.cr3, ncol=3)


#----------------#



#---------------------------#


#----------------------------#
# Phylogenetic signals!?
# 	low p-values, suggesting non-random distribution of trait values, that is, fit well to phylogeny, indicating that evolutionary histroy have had a relevant impacrt on the trait value. and that one needs to correct for phylogeny if one wants to perform stratistical comparisons (say of habitat classes)

x = metadata[tree$tip.label,'group']
names(x) = tree$tip.label

phytools::phylosig(tree,x,method="lambda",test=TRUE)





#Need a model for discrete variables
SS.phylosig <- geiger::fitDiscrete(tree, x)
summary(SS.phylosig)
plot(SS.phylosig,main="fitDiscrete model=\"ARD\"",show.zeros=FALSE,
     cex.traits=0.9)

phytools::phylosig(tree, x, method="K", test=TRUE, nsim=10000)
phytools::phylosig(tree, x, method="lambda", test=TRUE, nsim=10000)



#--------------------#
set.seed(1234)

sbtree = drop.tip(tree, tree$tip.label[!(tree$tip.label %in% df.chkm$BinId)] )
x = df.chkm[sbtree$tip.label,'GC']
names(x) = sbtree$tip.label

phylosig(sbtree, x, method="K", test=TRUE, nsim=10000)
phylosig(sbtree, x, method="lambda", test=TRUE, nsim=10000)


x = df.chkm[sbtree$tip.label,'Coding_density']
names(x) = sbtree$tip.label
phylosig(sbtree, x, method="K", test=TRUE, nsim=10000)
phylosig(sbtree, x, method="lambda", test=TRUE, nsim=10000)


x = df.chkm[sbtree$tip.label,'Genome_size']
names(x) = sbtree$tip.label
phylosig(sbtree, x, method="K", test=TRUE, nsim=10000)
phylosig(sbtree, x, method="lambda", test=TRUE, nsim=10000)

#----------------------------#

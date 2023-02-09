phi = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/PHI_output.txt')
resultsDir <- "~/DATA/MarbGenomics/Ortho/OrthoFinder/Results_Oct23_2"

#---------------------------------
#
#   OrthoFinder2Pan
#
#----------------------------------
# ---- partially obsolete? as we now use loadOGs() function in AnnotateTable!

#	load_orthogroups
#	splits the comma separated paralogs, and saves them in list of lists
load_orthogroups <- function(inPath) {
	dfOrthogroup <- read.table(inPath , header=T, sep = "\t", row.names = 1, stringsAsFactors=F)
	lsOrthogroup <- strsplit( as.matrix(dfOrthogroup),", ")

	#tidyverse::enframe(lsOrthogroup) %>% # creates the 'value' as a `list` column
	#  mutate(value = map(value, as.character)) %>% # change to single type
  # unnest(col=c(value))

	#dim(lsOrthogroup) <- dim(dfOrthogroup)
	names(lsOrthogroup) <- rownames(dfOrthogroup)
	return(lsOrthogroup)
}

#Example
lsOrthogroup = load_orthogroups(inPath = '~/DATA/MarbGenomics/Ortho/OrthoFinder/Results_Oct23_3/Orthogroups/Orthogroups.tsv')
lsOrthogroup["OG0001369", ]
lsOrthogroup[["OG0001369", "Marinobacter_sp_LZ_8"]]



##-------------------------------------------------
#	select genes in OGs
#	default, genome is set to all, showing all genes in OG per genome
#	you can specify the genome to obtain those genes in the genome of choise
#	currently does not work with more than one?

genes_in_OG <- function(listOG, OG=NULL, genome=NULL) {
	if(is.null(genome)){
		out = listOG[OG,]
	}else{
		out = listOG[OG, genome]
	}
	return(out)
}

#Examples
genes_in_OG(listOG = lsOrthogroup, OG = 'OG0001369', genome = 'Marinobacter_sp_LZ_8')
genes_in_OG(listOG = lsOrthogroup, OG = 'OG0001369')

#or
lsOrthogroup["OG0001369", ]
lsOrthogroup[["OG0001369", "Marinobacter_sp_LZ_8"]]


##-------------------------------------------------
#	convert to pan-genome table, with sizes of orthogroups per genome
#

og2pan <- function(listOG, clean_genome_names=TRUE) {
	out <- apply(listOG,2,sapply,length)
	if(clean_genome_names==TRUE){
	  colnames(out) = gsub('-','_',gsub('\\.','_', colnames(out)))
	}
	return(t(out))
}

#Examples
OFpan = og2pan(lsOrthogroup, clean_genome_names=TRUE)
#colnames(pan) = gsub('-','_',gsub('\\.','_', colnames(pan)))

#-----


##-------------------------------------------------
#	genomes represented per orthogroup
#

pan2represented <- function(pan) {
  out <- data.frame(table(apply(t(pan)>0,1,sum)))
  colnames(out) <- c('genomes','orthogroups')
  return(out)
}


#Examples
OFpan.represented = pan2represented(OG.pan)

ggplot(OFpan.represented , aes(x = as.numeric(genomes), y = orthogroups)) +
					geom_col(fill='#1380A1',color='white') +
					scale_y_continuous(expand = c(0, 0)) +
					xlab("") +
					ylab("no. orthogroups")+
					xlab("no. genomes represented")+
					theme_bw() +
					theme(panel.border = element_blank(),
						panel.grid.major = element_blank(),
						panel.grid.minor = element_blank(),
						axis.line = element_line(colour = "black"),
						axis.title = element_text(size=11),
						axis.text = element_text(size=11),
						legend.position = "none")


#----------


subpan = OFpan[,grepl('Marinobacter|Tamilnaduibacter',colnames(OFpan))]
subpan = subpan[rowSums(subpan)>0,]



##-------------------------------------------------
#	Detect Single Copy Orthologs


pan2sco <- function(pan) {
	pan = t(pan)
  SCO = apply(pan==1,1,all)
	SCO = names(SCO[SCO ==TRUE])
	return(t(SCO))
}

#Examples
sco = pan2sco(OG.pan)
length(sco)

#----
subsco = pan2sco(subpan)




genes_in_OG(listOG = lsOrthogroup, OG = 'OG0001172')
genes_in_OG(listOG = lsOrthogroup, OG = sco)





#------------------------

library(micropan) # readFasta and writeFasta
library(phangorn) # tree inference

#loading fasta files



load_faa_sequences <- function(inPath, suffix = 'faa') {
	fastaFiles <- dir(inPath, pattern = paste('*.', suffix ,sep=''), full.names = T)
	names(fastaFiles) <- basename(fastaFiles)

		lapply(fastaFiles, function(fastaFile){
	  		dfGenomeSeq <- data.frame(readFasta(fastaFile))
  			print(paste('...loading', fastaFile))
  			row.names(dfGenomeSeq) <- sub("^([^ ]+).*", "\\1", dfGenomeSeq$Header, perl=T)
  			return(dfGenomeSeq)
			}) -> allSequences

	return(allSequences)
}


allSequences = load_faa_sequences(inPath="~/DATA/MarbGenomics/Ortho", suffix = 'faa')

#show head of sequences in genome 1
head(allSequences[[1]])



gOG = genes_in_OG(listOG = lsOrthogroup, OG = 'OG0000018')
gOG = genes_in_OG(listOG = lsOrthogroup, OG = 'OG0001172')



annotate_ogs <- function(alSeq , OG , suffix = 'faa'){
	conv = names(alSeq)
	names(conv) = paste(gsub("-",".", gsub(paste('.', suffix ,sep=''),'', conv)), paste('.', suffix ,sep=''),sep='')
	lsAnn = list()

	for(selOG in OG){
		gOG = genes_in_OG(listOG = lsOrthogroup, OG = selOG)
		print(selOG)
		tmpDF=c()
		for(genome in names(gOG)){
			fname = as.character(conv[paste(genome,'.faa',sep='')])
			if(nrow(alSeq[[fname]][as.character(unlist(gOG[genome])),])>0){
					tmpDF = rbind(tmpDF, cbind('file'= fname,'genome'= genome, alSeq[[fname]][as.character(unlist(gOG[genome])),]))
				}
			}
	tmpDF = cbind(tmpDF,colsplit(tmpDF $Header," ",c("locus_tag","product")))
	lsAnn[[selOG]] = tmpDF
	}

	return(lsAnn)
}

#lsAnn = annotate_ogs(allSequences, OG = 'OG0000018',suffix='faa')
sco_ann = annotate_ogs(allSequences, OG = sco ,suffix='faa')


lapply(names(sco_ann), function(x) unique(sco_ann[[x]][,'product']))

allhits =lapply(names(sco_ann), function(x) unique(sco_ann[[x]][,'product']))
names(allhits) = names(sco_ann)
 allhits[grepl('ribosomal protein',allhits)]

ogsHit = names( allhits[grepl('ribosomal protein',allhits)] )


#--------------------------------------------

#----------------------
#	PAN stats
#	total number of Orthologous groups

number_orthologroups = nrow(pan)
number_genomes = ncol(pan)
number_of_genes = sum(rowSums(pan))
singletons = length(rowSums(pan)[rowSums(pan)==1])
non_singletons = length(rowSums(pan)[rowSums(pan)>1])
number_of_species_specific_orthogroups = length(rowSums(pan.presabs)[rowSums(pan.presabs)==1])
number_of_genes_in_species_specific_orthogroups = sum(rowSums(pan[species_specific_orthogroups,]))
number_shared_orthogroups_min_2_genomes = length(rowSums(pan.presabs)[rowSums(pan.presabs)>1])
number_of_orthogroups_present_in_all = length(rowSums(pan.presabs)[rowSums(pan.presabs)==ncol(pan.presabs)])
orthogroups_present_in_all = names(rowSums(pan.presabs)[rowSums(pan.presabs)==ncol(pan.presabs)])
orthogroups_of_size_nspecies = names(rowSums(pan)[rowSums(pan)==ncol(pan)])
single_copy_orthologs = intersect(orthogroups_of_size_nspecies, orthogroups_present_in_all)
number_of_single_copy_orthologs = length(intersect(orthogroups_of_size_nspecies, orthogroups_present_in_all))
species_specific_orthogroups = names(rowSums(pan.presabs)[rowSums(pan.presabs)==1])
mean_orthogroup_size = mean(rowSums(pan))




paste(
  c(	number_genomes,' genomes, ',
     number_of_genes,' genes were grouped into ',
     number_orthologroups,' orthologous groups',
     ' of which ', singletons,' were singletons (',round((100 * singletons /number_orthologroups),2),'%)'
  ),
  collapse='')

#-------------------------






sco_ann [['OG0001206']]
sco_ann [['OG0001207']]
OG0001261
OG0001366
OG0001448


sco_ann[1]

phihits = phi[phi$ALL_SIG == TRUE,]$OG
phi_ann = annotate_ogs(allSequences, OG = phihits ,suffix='faa')

phi_ann[[phihits[1]]][,c(2,3)]
phi_ann[[phihits[2]]][,c(2,3)]
phi_ann[[phihits[3]]][,c(2,3)]
phi_ann[[phihits[4]]][,c(2,3)]
phi_ann[[phihits[5]]][,c(2,3)]



lapply(names(phi_ann), function(x) unique(phi_ann[[x]][,'Header']))
lapply(names(sco_ann), function(x) unique(sco_ann[[x]][,'Header']))

nonPHI_SCO = names(sco_ann)[!names(sco_ann) %in% names(phi_ann)]
lapply(nonPHI_SCO, function(x) unique(sco_ann[[x]][,'Header']))
names()


lapply(nonPHI_SCO, function(x) unique(sco_ann[[x]][,'Header']))



#This is interesting, we can look at the length of each sequence, and say that a gene can not be longer than x% of the shortest sequence?
lengthVariation = lapply(sco, function(x) nchar(sco_ann[[x]][,'Sequence']))
names(lengthVariation) = sco

dflengthVariation = data.frame(
	cbind(
		'OG'= sco,
		'mean'= lapply(names(lengthVariation), function(x) mean(lengthVariation[[x]])),
		'sd'= lapply(names(lengthVariation), function(x) sd(lengthVariation[[x]])),
		'se'=lapply(names(lengthVariation), function(x) sd(lengthVariation[[x]])/sqrt(length(lengthVariation[[x]])))
		)
	)

head(dflengthVariation)


dflengthVariation[dflengthVariation$sd>10,]


var_too_long = unlist(dflengthVariation[dflengthVariation$sd>10,]$OG)


sco[!sco %in% c(var_too_long, phihits) ]
venn(list('to_long'=var_too_long, 'PHI'=phihits))


filtered_sco = sco[!sco %in% c(var_too_long, phihits) ]

filtered_sco_ann = annotate_ogs(allSequences, OG = filtered_sco ,suffix='faa')
lapply(names(filtered_sco_ann), function(x) unique(filtered_sco_ann[[x]][,'Header']))



length(var_too_long)


var_too_long_ann = annotate_ogs(allSequences, OG = intersect(var_too_long, phihits) ,suffix='faa')
lapply(names(var_too_long_ann), function(x) unique(var_too_long_ann[[x]][,'Header']))
lapply(names(var_too_long_ann), function(x) unique(var_too_long_ann[[x]][,'Sequence']))


filtered_sco_ann

#--------------------------------

#--------------------
#kinky function to copy files around! subsetting to new sets of orthologs

subset_fastaFiles <- function(input_dir , out_dir, subset , suffix = 'fa'){
	files = list.files(input_dir,pattern=paste('*.',suffix,sep=''),full.names=TRUE)
	files_filter = files[grepl(paste(subset,collapse='|'),files)]
	dir.create(out_dir)
	file.copy(files_filter, out_dir)
}


#--------------------------------

dir.create("~/DATA/MarbGenomics/filtered_alignments/")

subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_protein_alignments_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/SCO_filtered_full/",
					subset = filtered_sco,
					suffix = 'fa')

subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_protein_alignments_mafft_genafpair_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/SCO_filtered_full_mafft_genafpair/",
					subset = filtered_sco,
					suffix = 'fa')

subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_pal2nal_asFasta_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/SCO_filtered_full_pal2nal/",
					subset = filtered_sco,
					suffix = 'fa')


subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_mafft_genafpair_pal2nal_asFasta_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/SCO_filtered_full_mafft_genafpair_pal2nal/",
					subset = filtered_sco,
					suffix = 'fa')

#ogsHit for Ribosomal Proteins

subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_protein_alignments_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/RibosomalProteins/",
					subset = ogsHit,
					suffix = 'fa')

subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_protein_alignments_mafft_genafpair_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/RibosomalProteins_mafft_genafpair/",
					subset = ogsHit,
					suffix = 'fa')

subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_pal2nal_asFasta_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/RibosomalProteins_pal2nal/",
					subset = ogsHit,
					suffix = 'fa')


subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/SCO_mafft_genafpair_pal2nal_asFasta_GenomeHeader/",
					out_dir="~/DATA/MarbGenomics/filtered_alignments/RibosomalProteins_mafft_enafpair_pal2nal/",
					subset = ogsHit,
					suffix = 'fa')


#--------------

#!!!!!!!!
#extract the sequences from results OF?



subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/Ortho/OrthoFinder/Results_Oct23_3/Orthogroup_Sequences/",
					out_dir="~/DATA/MarbGenomics/SCO_MB/",
					subset = subsco,
					suffix = 'fa')

scosub_ann = annotate_ogs(allSequences, OG = subsco ,suffix='faa')



lengthVariation = lapply(subsco, function(x) nchar(scosub_ann[[x]][,'Sequence']))
names(lengthVariation) = subsco

dflengthVariation = data.frame(
	cbind(
		'OG'= subsco,
		'mean'= lapply(names(lengthVariation), function(x) mean(lengthVariation[[x]])),
		'sd'= lapply(names(lengthVariation), function(x) sd(lengthVariation[[x]])),
		'se'=lapply(names(lengthVariation), function(x) sd(lengthVariation[[x]])/sqrt(length(lengthVariation[[x]])))
		)
	)

head(dflengthVariation)


dflengthVariation[dflengthVariation$sd>10,]


var_too_long = unlist(dflengthVariation[dflengthVariation$sd>10,]$OG)



subset_fastaFiles(	input_dir="~/DATA/MarbGenomics/Ortho/OrthoFinder/Results_Oct23_3/Orthogroup_Sequences/",
					out_dir="~/DATA/MarbGenomics/SCO_filtered_MB/",
					subset = subsco[!subsco %in% var_too_long],
					suffix = 'fa')




#--------------------------------------------------
#	Evolutionary distances derived from pan-genome
#	GGR
#	-----------
#	gen repertoire relatedness index (GRR).
#	Between two genomes defined as the number of common gene families (the intersection) divided by the number of genes in the smallest genome.
#	It is close to 100% if the gene repertoires are very similar (or one is a subset of the orhter) and lower otherwise
##-------------------------------------------------
#	convert to pan-genome table, with sizes of orthogroups per genome


pan2ggr <- function(pan) {
	pan = t(pan)
  g1tr = c()
	for(g1 in colnames(pan)){
		pan.presabs = pan
		pan.presabs[pan.presabs > 0] <- 1
		g2tr = c()
		for(g2 in colnames(pan)){
			sub= pan.presabs[,c(g1, g2)]
			sub = sub[rowSums(sub)>0,]

			n.shared = length(rowSums(sub)[rowSums(sub) ==2])
			ng1 = unname(colSums(sub)[1])
			ng2 = unname(colSums(sub)[2])
			smallest=sort(c(ng1,ng2))[1]

			ggrval = n.shared/smallest
			g2tr = c(g2tr, ggrval)
		}
		g1tr = rbind(g1tr, g2tr)
	}
	colnames(g1tr) = colnames(pan)
	rownames(g1tr) = colnames(pan)

	g1tr = g1tr[order(row.names(g1tr)),]
	g1tr = g1tr[,order(colnames(g1tr))]

	return(g1tr)
}

#Examples
pan.ggr = pan2ggr(OG.pan)
head(pan.ggr)
pan.ggr = pan2ggr(OG.pan)

#---------

MDS = cmdscale(pan.ggr, eig = TRUE)
MDSdata = data.frame(MDS$points)
MDSdata$Name = rownames(MDSdata)
MDSdata$Genus = gsub( "_.*$", "", rownames(MDSdata) )
colnames(MDSdata) = c("x",'y','Name',"Genus")


ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(Genus),label=Name)) +
	theme_classic() +
	xlab('Dimension 1') +
	ylab('Dimension 2') +
	labs(fill = "Genus") +
	ggtitle(label='GRR')+
	geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_point(shape = 21, size = 2) +
	theme(
		aspect.ratio = 1,
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
        legend.justification = c(1, 1),
        legend.key.width = unit(0.25, 'cm'),
        legend.key.height = unit(0.55, 'cm'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))



#------------
#DISTANCES

concat_tree = read.tree("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile")
COPH = cophenetic(concat_tree)
COPH = COPH[ order(row.names(COPH)), ]
COPH = COPH[ , order(colnames(COPH))]



distJac <- distJaccard(t(OFpan))
distJac= as.matrix(distJac)
distJac = distJac[ order(row.names(distJac)), ]
distJac = distJac[ , order(colnames(distJac))]
subDistJac = distJac[grepl('Marinobacter|Tamil',rownames(distJac)),grepl('Marinobacter|Tamil',colnames(distJac))]


distMan <- distManhattan(t(OFpan))
distMan = as.matrix(distMan)
distMan = distMan[ order(row.names(distMan)), ]
distMan = distMan[ , order(colnames(distMan))]
subDistMan = distMan[grepl('Marinobacter|Tamil',rownames(distMan)),grepl('Marinobacter|Tamil',colnames(distMan))]



pan[,grepl('Marinobacter|Tamul',colnames(pan))]

df.distances = data.frame(cbind('cophenetic'=c(COPH),'Manhattan'=c(distMan),'Jaccard'=c(distJac),'GRR'=c(pan.ggr)))
cop.vs.all = melt(df.distances,id='cophenetic')


ggplot(cop.vs.all,aes(x= cophenetic,y= value)) +
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
	facet_wrap(~variable,scales='free_y',strip.position = "left",
                labeller = as_labeller(c(Manhattan = "Manhattan distance", Jaccard = "Jaccard distance",GRR='GRR')))+
     theme(strip.background = element_blank(),
           strip.placement = "outside",
           strip.text = element_text(size = 11))+ylab('')+xlab('Cophenetic distance')



#-------- Single plot
#ggplot(data.frame(df.distances), aes(x= cophenetic,y= Manhattan)) +
#	geom_hex(bins = 35, colour = NA) +
#	scale_fill_viridis(option="magma",trans = 'sqrt', name = 'Frequency')+
#	theme_bw() +
#	theme(
#		aspect.ratio = 1,
#		axis.title = element_text(size = 11, colour = '#000000'),
#        axis.text = element_text(size = 10, colour = '#000000'),
#        legend.justification = c(1, 1),
#        legend.key.width = unit(0.25, 'cm'),
#       legend.key.height = unit(0.55, 'cm'),
#        legend.text = element_text(size = 7.5),
#       legend.text.align = 1,
#        legend.title = element_text(size = 9))+
#	theme(panel.border = element_blank(),
#		panel.grid.major = element_blank(),
#		panel.grid.minor = element_blank(),
#		axis.line = element_line(colour = "black"),
#		axis.title = element_text(size=11),
#		axis.text = element_text(size=11),
#		)



MDS = cmdscale(distJac, eig = TRUE)
MDSdata = data.frame(MDS$points)
MDSdata$Name = rownames(MDSdata)
MDSdata$Genus = gsub( "_.*$", "", rownames(MDSdata) )
colnames(MDSdata) = c("x",'y','Name',"Genus")

ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(Genus),label=Name)) +
	theme_classic() +
	xlab('Dimension 1') +
	ylab('Dimension 2') +
	labs(fill = "Genus") +
	ggtitle(label='Jaccard distance')+
	geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
	geom_point(shape = 21, size = 2) +
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




MDS = cmdscale(distMan, eig = TRUE)
MDSdata = data.frame(MDS$points)
#MDSdata$nrofgenes =  str_count(rownames(MDSdata) , pattern = "_")
#MDSdata$nrofgenes =  MDSdata$nrofgenes + 1
MDSdata$Name = rownames(MDSdata)
MDSdata$Genus = gsub( "_.*$", "", rownames(MDSdata) )
colnames(MDSdata) = c("x",'y','Name',"Genus")



ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(Genus),label=Name)) +
	theme_classic() +
	xlab('Dimension 1') +
	ylab('Dimension 2') +
	labs(fill = "Genus") +
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




#---------

MDS = cmdscale(subDistMan, eig = TRUE)
MDSdata = data.frame(MDS$points)
#MDSdata$nrofgenes =  str_count(rownames(MDSdata) , pattern = "_")
#MDSdata$nrofgenes =  MDSdata$nrofgenes + 1
MDSdata$Name = rownames(MDSdata)
MDSdata$Genus = gsub( "_.*$", "", rownames(MDSdata) )
colnames(MDSdata) = c("x",'y','Name',"Genus")

  ggplot(MDSdata,aes(x=x,y=y,fill=as.factor(Genus),label=Name)) +
	theme_classic() +
	xlab('Dimension 1') +
	ylab('Dimension 2') +
	labs(fill = "Genus") +
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


#---------------------------------------------

library(Rtsne) # Load package
set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(as.matrix(t(pan))) # Run TSNE
#plot(tsne_out$Y,col=iris_unique$Species,asp=1) # Plot the result

plot(tsne_out$Y,col=as.factor(gsub( "_.*$", "", colnames(pan) )),asp=1) # Plot the result

tsne_model_1 = Rtsne(as.matrix(t(pan)), check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)


d_tsne_1 = as.data.frame(tsne_model_1$Y)

ggplot(d_tsne_1, aes(x=V1, y=V2,fill=as.factor(gsub( "_.*$", "", colnames(pan) )))) +
  	theme_classic() +
	xlab('tSNE 1') +
	ylab('tSNE 2') +
	labs(fill = "Genus") +
	ggtitle(label='pan-genome tSNE')+
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



#--------------------------------------------------------------
# Filtering Pan tables


pan_table = t(panOF)
dim(pan_table)

RFdf = data.frame(t(pan),'group'=as.factor(gsub( "_.*$", "", colnames(pan) )))
dim(RFdf)

PAN_nonzero_counts <- apply(RFdf, 1, function(y) sum(length(which(y > 0))))
hist(PAN_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of PCs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

pan_table_rare_removed <- remove_rare(table= t(pan), cutoff_pro=0.05)
dim(pan_table_rare_removed)

pan_table_rare_removed_norm <- sweep(pan_table_rare_removed, 2, colSums(pan_table_rare_removed) , '/')*100
pan_table_scaled <- scale(pan_table_rare_removed_norm, center = TRUE, scale = TRUE)
pan_table_asinh_mean_centred <- scale( asinh(t(pan)), center=TRUE, scale=FALSE)



#---------------------------------------------------------------
# RANDOM FOREST CLASSIFIER



RFdf = data.frame(t(pan),'group'=as.factor(gsub( "_.*$", "", rownames(t(pan)) )))


#RF <- randomForest(group ~ ., data= RFdf)
#RF <- randomForest(group ~ ., data= RFdf, importance=T, proximity=T,ntree=1500,keep.forest=F)

#TRAIN THE RF (look for code in RandomForest_Brieuc_modification.R)

RF <- randomForest(group ~ ., data= RFdf,importance=TRUE ,proximity=TRUE, ntree=30000, strata= group)




#RF_classify <- randomForest( x= RFdf[,1:(ncol(RFdf)-1)] , y= as.character(RFdf[ , ncol(RFdf)]) , ntree=501, importance=TRUE, proximities=TRUE )

#RF_state_classify_sig <- rf.significance( x= RF_classify ,  xdata= RFdf[,1:(ncol(RFdf)-1)] , nperm=1000 , ntree=501 )


RF_classify_imp <- as.data.frame( RF $importance )
RF_classify_imp$features <- rownames( RF_classify_imp )
RF_classify_imp_sorted <- arrange( RF_classify_imp  , desc(MeanDecreaseAccuracy)  )

barplot(RF_classify_imp_sorted $MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

barplot(RF_classify_imp_sorted[1:30,"MeanDecreaseAccuracy"], names.arg= RF_classify_imp_sorted[1:30,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=2, main="Classification RF")

#EggNog[RF_classify_imp_sorted[1:30,"features"],"eggNOG_HMM"]


heatmap.2(t(pan[RF_classify_imp_sorted[1:100,"features"],]),trace="none")




scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)


heatmap.2(t(pan[RF_classify_imp_sorted[1:50,"features"],]), col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)









#-------------
#-------------
#-------------
#-------------
library(kdetrees)
library(ape)
library("seqRFLP")
library(msa)

# TREES

library(ape)
library(phytools)
library(ggtree)

concat_tree = read.tree('~/DATA/MarbGenomics/SCO_genefpair_trimAl/concat.treefile')

p <- ggtree(concat_tree ) + geom_tiplab()

loci_trees = read.tree('~/DATA/MarbGenomics/SCO_genefpair_trimAl/loci.treefile')
p <- ggtree(loci_trees[2] ) + geom_tiplab()


#class(treelist) <- "multiPhylo"

#kdeRes =  kdetrees(loci_trees)


#RobinsonF_trees = multiRF(loci_trees)
#rownames(RobinsonF_trees) = names(loci_trees)
#colnames(RobinsonF_trees) = names(loci_trees)

CMD  = (cmdscale(RobinsonF_trees, eig = TRUE))

dfCMD = data.frame(CMD$points)
colnames(dfCMD) = c('x','y')
dfCMD$names = rownames(dfCMD)

ggplot(data= dfCMD, aes(x=x, y=y, label=names))+geom_point()+geom_text()
ggplot(data= dfCMD, aes(x=x, y=y))+geom_point()

#-----

concat_tree = read.tree("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile")
concat_tree = read.tree("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile")
COPH = cophenetic(concat_tree)

p <- ggtree(concat_tree) + geom_tiplab(size=2)

ggtree(concat_tree) +
	geom_text2(size=2,hjust=0.5,vjust=0.2, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+
	geom_tiplab(size=2)

ggtree(concat_tree) +
	geom_text2(size=2, hjust=1.2,vjust=-0.26, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+
	geom_tiplab(size=2)
	#-----


#-------------

outgroup = concat_tree$tip.label[!grepl("Marinobacter",concat_tree$tip.label)]
outgroup = outgroup[! outgroup %in% 'Tamilnaduibacter_salinus_Mi_7']

concat_tree.rooted = root(concat_tree, outgroup = outgroup, resolve.root = TRUE)
p <- ggtree(concat_tree.rooted ) + geom_tiplab()

#--------------


subtree = drop.tip(concat_tree,outgroup)

ggtree(subtree) +
	geom_text2(size=2, hjust=1.2,vjust=-0.26, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+
	geom_tiplab(size=2)



COPH = cophenetic(subtree)
kms2 = kmeans(COPH, 2)
kms3 = kmeans(COPH, 3)
kms4 = kmeans(COPH, 4)
kms5 = kmeans(COPH, 5)
kms6 = kmeans(COPH, 6)
kms7 = kmeans(COPH, 7)
kms8 = kmeans(COPH, 8)
kms9 = kmeans(COPH, 9)
kms10 = kmeans(COPH, 10)



CMD  = (cmdscale(COPH, eig = TRUE))

dfCMD = data.frame(CMD$points)
colnames(dfCMD) = c('x','y')
dfCMD$names = rownames(dfCMD)

ggplot(data= dfCMD, aes(x=x, y=y, label=names))+geom_point()+geom_text()
ggplot(data= dfCMD, aes(x=x, y=y))+geom_point()


dfCMD$c8 = kms8$cluster[dfCMD$names]

ggplot(data= dfCMD, aes(x=x, y=y,color=c8))+geom_point()




p.sub = ggtree(subtree)

p.sub$data  = p.sub$data %>% mutate(c2 = kms2$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c3 = kms3$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c4 = kms4$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c5 = kms5$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c6 = kms6$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c7 = kms7$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c8 = kms8$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c9 = kms9$cluster[label])
p.sub$data  = p.sub$data %>% mutate(c10 = kms10$cluster[label])


p.sub +
	geom_text2(size=2, hjust=1.2,vjust=-0.26, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+
	geom_tiplab(size=2,aes(color= c10)) +	scale_colour_viridis(direction = -1)

p.sub +
	geom_text2(size=2, hjust=1.2,vjust=-0.26, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+
	geom_tiplab(size=2,aes(color= c2)) +	scale_colour_viridis(direction = -1)






####------------------

library(gridExtra)

MDSdata$c10 = kms10$cluster[MDSdata$Name]

cols=colorRampPalette(colors = brewer.pal(n = 9,name = "Set1"))(length(unique(MDSdata$c10)))
p.mds =
	ggplot(MDSdata,aes(x=x,y=y,color=as.factor(c10),label=Name))+geom_point() + scale_color_manual(values=cols)+ylab("Dimension 2")+xlab("Dimension 1")+
    theme(axis.line = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour =   "#D9D9D9"),
          panel.grid.minor = element_line(colour = "#D9D9D9"),
          panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
          axis.ticks = element_line(colour = "black",size = 1),
          axis.text.x = element_text(family = "Arial",face = "plain",size =12,colour="black"),
          axis.text.y = element_text(family = "Arial",face="plain",size=12,colour="black"),
          axis.title.x = element_text(family = "Arial",face="plain",size = 14,colour = "black"),
          axis.title.y = element_text(family = "Arial",face="plain",size=14,colour="black"),
          legend.background = element_blank(),legend.key.size = unit(2,"point"),
          legend.title=element_blank(),legend.key = element_blank(),
          legend.text = element_text(size=12,family = "Arial",face = "plain",colour = "black"),legend.position ="right")


p.subann =
	p.sub +
		geom_text2(size=2, hjust=1.2,vjust=-0.26, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+
		geom_tiplab(size=2,aes(color= as.factor(c10))) +scale_color_manual(values=cols)+ theme(legend.position = "none")

grid.arrange(p.subann,p.mds,ncol=2)




#-----------------------


ggtree(concat_tree)$data








library(castor)

#we calculated the relative evolutionary vivergence (RED)of each node in a rooted phylogenetic tree.
#The RED of the nodes is a measure of its relative placement bewteen the root and the node's descending tips (Parks et al., 2018)
#The root's RED is always 0 and the RED pf each tip is 1, and the RED of each node is between 0 and 1

RED = get_reds(concat_tree)

#-------
#Can I caltulate a % sequence identity? which all remain decently close? eg >90% minimum?



outgroup = concat_tree$tip.label[!grepl("Marinobacter",concat_tree$tip.label)]
outgroup = outgroup[! outgroup %in% 'Tamilnaduibacter_salinus_Mi_7']

concat_tree.rooted = root(concat_tree, outgroup = outgroup, resolve.root = TRUE)
p <- ggtree(concat_tree.rooted ) + geom_tiplab()

p

#-----
#we calculated the relative evolutionary vivergence (RED)of each node in a rooted phylogenetic tree.
#The RED of the nodes is a measure of its relative placement bewteen the root and the node's descending tips (Parks et al., 2018)
#The root's RED is always 0 and the RED pf each tip is 1, and the RED of each node is between 0 and 1
RED = get_reds(concat_tree.rooted)

#Given a rooted phylogenetic tree, calculate the phylogenetic depth of each node (mean distance to its descending tips).
nodeDepths = get_all_node_depths(concat_tree.rooted)

#Calculate phylogenetic ("patristic") distances between all pairs of tips or nodes in the tree, or among a subset of tips/nodes requested.
pairwiseDist = get_all_pairwise_distances(concat_tree.rooted)

dfNodeStats = data.frame(cbind('RED'=RED,'nodeDepths'=nodeDepths))

dfNodeStats %>% ggplot(.,aes(RED, nodeDepths)) + geom_point() + geom_smooth()
dfNodeStats %>% ggplot(.,aes(RED, nodeDepths)) + geom_smooth(method='lm') + geom_point()





CMD  = (cmdscale(pairwiseDist, eig = TRUE))

dfCMD = data.frame(CMD$points)
colnames(dfCMD) = c('x','y')
dfCMD$names = rownames(dfCMD)

ggplot(data= dfCMD, aes(x=x, y=y, label=names))+geom_point()+geom_text()
ggplot(data= dfCMD, aes(x=x, y=y))+geom_point()














library(Biostrings)
library(ggmsa)
library(ggplot2)
protein_sequences <- system.file('~/DATA/MarbGenomics/SCO_protein_alignments_GenomeHeader/OG0001217.fa')
#ggmsa(protein_sequences, 164, 213, color = "Chemistry_AA")
#ggmsa(x , 164, 213, color = "Chemistry_AA")



x <- readAAStringSet('~/DATA/MarbGenomics/SCO_protein_alignments_GenomeHeader/OG0001306.fa')
d <- as.dist(stringDist(x, method = "hamming")/width(x)[1])
library(ape)
tree <- bionj(d)
library(ggtree)
p <- ggtree(tree ) + geom_tiplab()

data = tidy_msa(x)
p + geom_facet(geom = geom_msa, data = data,  panel = 'msa',
               font = NULL, color = "Chemistry_AA") +
    xlim_tree(1)







#ggmsa(x , 164, 213, color = "Chemistry_AA")
#ggmsa(x, color = "Chemistry_AA")





gOG = genes_in_OG(listOG = lsOrthogroup, OG = phihits[1])



phi_ann = annotate_ogs(allSequences, OG = phihits ,suffix='faa')


tmpDF[,c(1,2,3)]

head(tmpDF)
tmpDF


#------------------------
#Now extract the sequences for a specific orthogroup (OG0000008 in this example):


grpID <- "OG0001369"

# get sequences for an orthogroup for each species and combine with rbind
do.call( rbind, lapply(names(allSequences), function(spc){
  allSequences[[spc]][lsOrthogroup[[grpID,sub("\\.faa","", spc)]], ]
})) -> orthoSequences

#------------------------
#Save the sequences to a fasta file:

# I want to have a separate directory for the orthogroup fasta files
grpFastasDir <- "orthofinder/grpFastas"
dir.create(grpFastasDir)

writeFasta(orthoSequences, file.path(grpFastasDir, "OG0000018.faa"))









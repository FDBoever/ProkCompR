
#----------------------------------------------------------------------------------------
#output is composed of several colums
#for more information on the specific columns refer to PHIPACK paper/manual
~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output
#	OG: 	ORTHOGROUP
#	NSS: 	permutation based p-value for NSS
#	MAX: 	permutation based p-value for MAX
#	PHIp: 	permutation based p-value for PHI
#	PHIn: 	non-permutation based p-value for NSS
#	nallSites: 	number of sites
#	nInfSites: 	number of informative sites
#	estDiversity: 	estimated diversity (%) calculated as (pairwise deletion - ignoring missing/ambig)
#	npolyMorphSites: n polymorphic sites
#	nTaxa: 	number of taxa(/sequences) in alignment

#----------------------------------------------------------------------------------------

#libraries
library(gplots)
library(ggplot2)

#Loading data and clean/set rownames
phi = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/PHI_output.txt')
phi = read.delim('~/DATA/MarbGenomics/MB_SCO/PHI_output.txt')


OFs = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/Ortho/OrthoFinder/Results_Sep30/Orthogroups/Orthogroups.tsv')
rownames(OFs) = OFs$Orthogroup
phi$OG = gsub(' ','',phi$OG)


#cutoff can be 0.1 as used in literature (cite)
#cutoff = 0.05 #S6 S18 L32 L35
cutoff=0.1

PHISigList = list()
PHISigList[['NSS']] = as.character(na.omit(phi[phi$NSS< cutoff,]$OG))
PHISigList[['MAX']] = as.character(na.omit(phi[phi$MAX < cutoff,]$OG))
PHISigList[['PHI']] = as.character(na.omit(phi[phi$PHIp < cutoff,]$OG))

venPHISIG = venn(PHISigList)
intersectPHIPACK = attr(venPHISIG,"intersections")$`NSS:MAX:PHI`

phi$NSS_SIG = ifelse(phi$OG %in% PHISigList[['NSS']],TRUE,FALSE)s
phi$MAX_SIG = ifelse(phi$OG %in% PHISigList[['MAX']],TRUE,FALSE)
phi$PHI_SIG = ifelse(phi$OG %in% PHISigList[['PHI']],TRUE,FALSE)
phi$ALL_SIG = ifelse(phi$OG %in% intersectPHIPACK,TRUE,FALSE)
phi$estDiv = as.numeric(gsub('% ','',as.character(phi$estDiversity)))/100
phi$annotation =  as.character(OFs[phi$OG,2])

phi[phi$ALL_SIG == TRUE & grepl('ribosomal',phi$annotation,ignore.case=TRUE),]
phi[phi$ALL_SIG == TRUE,]

table(phi$ALL_SIG)

#--- percentages, nr of significant hits vs the entire set
table(phi$ALL_SIG)[2] / (table(phi$ALL_SIG)[1] + table(phi$ALL_SIG)[2] )
table(phi$NSS_SIG)[2] / (table(phi$NSS_SIG)[1] + table(phi$NSS_SIG)[2] )
table(phi$MAX_SIG)[2] / (table(phi$MAX_SIG)[1] + table(phi$MAX_SIG)[2] )
table(phi$PHI_SIG)[2] / (table(phi$PHI_SIG)[1] + table(phi$PHI_SIG)[2] )

#--- look for specific annotations, say ribosomal in Marinobacter case
SigPHIOfs[grepl('ribosomal',SigPHIOfs[,2],ignore.case=T),]

SigPHIOfs = phi[phi$ALL_SIG == TRUE,]



SigPHIOfs

#---------
# Visualisation
plot(venPHISIG)

ggplot(phi, aes(x= nallSites, fill= ALL_SIG)) +
	geom_histogram(color= 'white',position="identity",bins = 45)+
	scale_y_continuous( expand = c(0, 0)) +
	theme_classic()

ggplot(phi, aes(x= nInfSites, fill= ALL_SIG)) +
	geom_histogram(color= 'white',position="identity",bins = 45)+
	scale_y_continuous( expand = c(0, 0)) +
	theme_classic()

ggplot(phi, aes(x= npolyMorphSites, fill= ALL_SIG)) +
	geom_histogram(color= 'white',position="identity",bins = 45)+
	scale_y_continuous( expand = c(0, 0)) +
	theme_classic()

ggplot(phi, aes(x= estDiv, fill= ALL_SIG)) +
	geom_histogram(color= 'white',position="identity",bins = 45)+
	scale_y_continuous( expand = c(0, 0)) +
	theme_classic()

#--- SOME RESULT STATEMENTS
#x% of all x core genome genes were identified to have significant evidence for recombination (p<0.05) that was detected by at least one of the three tests
#notably,\more informative sites and higher nucelotide diversutt
#-----




library(seqinr)
library(chopper)

fastafiles = list.files(path = "~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_pal2nal_asFasta_GenomeHeader",pattern='.fa',full.names=TRUE)

#for(fst in fastafiles){
#	print(paste('converting',fst,sep=' '))
#	fas2phy(file = fst, format='sequential')
#}

#---
concatenateAlignments(	pattern='.fa',
						input.format='fasta',
						path = '~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_pal2nal_asFasta_GenomeHeader',
						output="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_pal2nal_concatenated_alignment.fasta")



startPosition = 1
partitions = c()
for(alnFile in fastafiles){
	OGname = basename(alnFile)
	aln = read.alignment(alnFile, format='fasta')
	lnght = nchar(aln[[3]][[1]])
	partitions = c(partitions ,paste(c('DNA, ', OGname,' = ',startPosition,'-', lnght+(startPosition-1)),collapse=''))
	startPosition=lnght+startPosition
}

writeLines(partitions,"~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_pal2nal_concatenated_alignment.partitions")



#IF DNA
paste('DNA')

#IF AA
paste('WAG')






fastafiles


#---------------------------







#fastafiles[grepl(paste(SigPHIOfs$OG,collapse='|'), fastafiles)]


dir.create("~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_asFasta_GenomeHeader")
dir.create("~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_asFasta_GenomeHeader")



# Select PHI filtered files and stor in new directory
file.copy(fastafiles[!grepl(paste(SigPHIOfs$OG,collapse='|'), fastafiles)], "~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_asFasta_GenomeHeader")


concatenateAlignments(	pattern='.fa',
						input.format='fasta',
						path = '~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_asFasta_GenomeHeader',
						output="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_concatenated_alignment.fasta")



fastafiles = list.files(path = "~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_asFasta_GenomeHeader",pattern='.fa',full.names=TRUE)

startPosition = 1
partitions = c()
for(alnFile in fastafiles){
	OGname = basename(alnFile)
	aln = read.alignment(alnFile, format='fasta')
	lnght = nchar(aln[[3]][[1]])
	partitions = c(partitions ,paste(c('DNA, ', OGname,' = ',startPosition,'-', lnght+(startPosition-1)),collapse=''))
	startPosition=lnght+startPosition
}

writeLines(partitions,"~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/SCO_noPHI_pal2nal_concatenated_alignment.partitions")








######---- OLD BABBLE

v = as.character(na.omit(phi[phi$PHIp<0.05,]$OG))
v = gsub(' ','',v)

OFs[OFs$Orthogroup %in% v ,c(1,2,3)]

#-----

SigPHIOfs = OFs[OFs$Orthogroup %in% v ,c(1,2,3)]
SigPHIOfs[grepl('ribosomal',SigPHIOfs[,2],ignore.case=T),]

OFs[OFs$Orthogroup == as.character(na.omit(phi[phi$PHIp<0.05,]$OG))[1]  ,c(1,2,3)]
OFs[OFs$Orthogroup == 'OG0000678'  ,c(1,2,3)]
OFs[OFs$Orthogroup == 'OG0000699'  ,c(1,2,3)]


OFs[OFs$Orthogroup == as.character(na.omit(phi[phi$PHIp<0.05,]$OG))[2]  ,c(1,2,3)]
OFs[OFs$Orthogroup == as.character(na.omit(phi[phi$PHIp<0.05,]$OG))[3]  ,c(1,2,3)]


OFs[OFs$Orthogroup == as.character(na.omit(phi[phi$PHIp<0.05,]$OG))[3]  ,c(1,2,3)]

as.character(OFs[grepl('ribosomal', OFs[,2],ignore.case=T),1])

allribosomal = phi[phi$OG %in% as.character(OFs[grepl('ribosomal', OFs[,2],ignore.case=T),1]),]

plot(allribosomal$PHIp)





OFs[OFs$Orthogroup == as.character(na.omit(phi[phi$PHIp<0.05,]$OG))[3]  ,c(1,2,3)]

phi[phi$NSS<0.05 & phi$MAX<0.05 & phi$PHIp<0.05,]$OG

phi[phi$NSS<0.05 & phi$MAX<0.05 & phi$PHIp<0.05,]$OG

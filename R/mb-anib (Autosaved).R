ANI16S = read.delim('~/DATA/MarbGenomics/ANIb_16Sseq_output/ANIb_percentage_identity.tab')
rownames(ANI16S)= ANI16S[,1]
ANI16S = ANI16S[,c(2:(ncol(ANI16S)))]


ANIb = read.delim('~/DATA/MarbGenomics/ANIb_output/ANIb_percentage_identity.tab')
rownames(ANIb)= ANIb[,1]
ANIb = ANIb[,c(2:(ncol(ANIb)))]

head(ANIb)


mergedEnv = read.delim('~/DATA/MarbGenomics/mergedEnv.txt')

#-----------------
# Reformatting the names so it is consistent for each data type
tree = read.tree('~/DATA/MarbGenomics/SCO_genefpair_trimAl/concat.treefile')

tree$tip.label[tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
tree$tip.label = gsub('-','_',gsub('\\.','_',tree$tip.label))
rownames(ANIb) = gsub('-','_',gsub('\\.','_', rownames(ANIb)))
colnames(ANIb) = gsub('-','_',gsub('\\.','_', colnames(ANIb)))

rownames(ANI16S) = gsub('-','_',gsub('\\.','_', rownames(ANI16S)))
colnames(ANI16S) = gsub('-','_',gsub('\\.','_', colnames(ANI16S)))

colnames(pan) = gsub('-','_',gsub('\\.','_', colnames(pan)))

rownames(mergedEnv) = gsub('-','_',gsub('\\.','_', mergedEnv$genome_name_dotted))

#------------------------------

ANIb[tree$tip.label, tree$tip.label]

ANIb = ANIb[tree$tip.label, tree$tip.label]
ANI16S = ANI16S[intersect(rownames(ANI16S), tree$tip.label), intersect(rownames(ANI16S), tree$tip.label)]



#-------------------------------


library(phangorn)

#1-ANIb as distance matrix, using Neighbor-joining algorithm implemented in phangorn.

nj.ANIb = NJ(1-as.matrix(ANIb))
plot(nj.ANIb)




#countTable2vector =

ctv = unname(unlist(c(ANIb)))
hist(ctv)

ggplot() + aes(ctv)+ geom_histogram(colour="black", fill="white")
ggplot() + aes(ctv)+ geom_density(fill='grey')

CMD  = (cmdscale(ANIb, eig = TRUE))
dfCMD = data.frame(CMD$points)
colnames(dfCMD) = c('x','y')
dfCMD$names = rownames(dfCMD)
ggplot(data= dfCMD, aes(x=x, y=y, label=names))+geom_point()+geom_text()
ggplot(data= dfCMD, aes(x=x, y=y))+geom_point()

#dfCMD$c8 = kms8$cluster[dfCMD$names]

#ggplot(data= dfCMD, aes(x=x, y=y,color=as.factor(c8)))+geom_point()+scale_color_manual(values=cols)
#ggplot(data= dfCMD, aes(x=x, y=y,color=as.factor(c8)))+geom_point()+scale_color_manual(values=cols)

#---------------------------

#ANIb bewteen gebines were calculated using pyANI, with default parameters.
#Cliques are a complete graph wghere each vertex represents a genome and every vcertex is linked with evry other vertex, links are formed between genomes at a specific ANI threshold such as >= 95% (See Parks, species cliques)


library(igraph)
#full network isolates
#sim16s2 =  sim16s[!grepl('Marinobacter_sp_FDB33|Alteromonas_sp_FDB36',rownames(sim16s)),!grepl('Marinobacter_sp_FDB33|Alteromonas_sp_FDB36',colnames(sim16s))]
#those from bioassay
#sim16s2 =  sim16s[rownames(trait.binairy),rownames(trait.binairy)]



g <- graph_from_adjacency_matrix(as.matrix(ANI16S), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <0.99))
plot(g2)
wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)



g <- graph_from_adjacency_matrix(as.matrix(ANIb), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <0.95))
plot(g2)
wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)


cliques = membership(wc)


#--------------------------------------------



cliques[cliques==12]

sortedCliqueSizes = sort(table(cliques),decreasing=TRUE)
ggplot(as.data.frame(sortedCliqueSizes), aes(x=cliques, y = Freq)) +
	scale_y_continuous(expand = c(0, 0)) +
	geom_bar(stat="identity",fill='#999999',color='white')+
	scale_color_manual(values='lightgrey')+
	ylab("Genomes")+
	xlab("ANI (95%) clique")+
	theme_bw() +
	theme(panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"),
		axis.title = element_text(size=11),
		axis.text = element_text(size=11),
		legend.position = "none")

ggplot(as.data.frame(sortedCliqueSizes[sortedCliqueSizes>1]), aes(x=cliques, y = Freq)) +
	scale_y_continuous(expand = c(0, 0)) +
	geom_bar(stat="identity",fill='#999999',color='white')+
	scale_color_manual(values='lightgrey')+
	ylab("Genomes")+
	xlab("ANI (95%) clique")+
	theme_bw() +
	theme(panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"),
		axis.title = element_text(size=11),
		axis.text = element_text(size=11),
		legend.position = "none")



#---------------------------------------

#none singleton cliques

non_signleton_cliques = sortedCliqueSizes[sortedCliqueSizes>1]
singleton_cliques = sortedCliqueSizes[sortedCliqueSizes==1]

length(sortedCliqueSizes)
length(non_signleton_cliques)
length(singleton_cliques)


#---------------------------------------

for(clq in names(non_signleton_cliques)){
	print('--------------------------')
	print(paste('---', clq))
	print(names(cliques[cliques == clq]))
	print('---')
	print(length(names(cliques[cliques == clq])))

	print(mergedEnv[names(cliques[cliques == clq]),c('TypeStrain', 'geo_loc_name','SS')])
	print('is monophyletic?')
	print(is.monophyletic(tree, tips= names(cliques[cliques == clq])))
}





#------------------------------------------------




#Those genomes in clique 5
names(cliques[cliques==5])

#Check if names correspond to the tree
names(cliques[cliques==5]) %in% tree$tip.label
gsub('\\.','_',names(cliques[cliques==5])) %in% tree$tip.label
gsub('\\.','-',names(cliques[cliques==5])) %in% tree$tip.label

#Check if clique is monophyletic in phylogeny
names(cliques[cliques==5]) %in% tree$tip.label



g <- graph_from_adjacency_matrix(as.matrix(ANIb), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <0.94))
plot(g2)

wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)

#-----------------------

g <- graph_from_adjacency_matrix(as.matrix(ANIb), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <0.75))
plot(g2)

wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)

#------------------------












par(mfrow=c(2,5))
#par(mfrow=c(1,1))

g5 <- delete.edges(g, which(E(g)$weight <0.99))
E(g5)$weight = 0.5
LO = layout_with_fr(g5)
V(g5)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g5))))]
#plot(g5, layout=LO, vertex.size=5, vertex.label="")
plot(g5, vertex.size=5, vertex.label="")

wc <- cluster_walktrap(g5)
title(paste(length(unique(membership(wc))), 0.99))

g3 <- delete.edges(g, which(E(g)$weight <0.98))
E(g3)$weight = 0.5
LO = layout_with_fr(g3)
V(g3)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g3))))]
#plot(g3, layout=LO, vertex.size=5, vertex.label="")
plot(g3, vertex.size=5, vertex.label="")

wc <- cluster_walktrap(g3)
title(paste(length(unique(membership(wc))), 0.98))


g4 <- delete.edges(g, which(E(g)$weight <0.97))
E(g4)$weight = 0.5
LO = layout_with_fr(g4)
V(g4)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g4))))]
#plot(g4, layout=LO, vertex.size=5, vertex.label="")
plot(g4, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g4)
title(paste(length(unique(membership(wc))), 0.97))


g5 <- delete.edges(g, which(E(g)$weight <0.96))
E(g5)$weight = 0.5
LO = layout_with_fr(g5)
V(g5)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g5))))]
#plot(g5, layout=LO, vertex.size=5, vertex.label="")
plot(g5, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g5)
title(paste(length(unique(membership(wc))), 0.96))

g3 <- delete.edges(g, which(E(g)$weight <0.95))
E(g3)$weight = 0.5
LO = layout_with_fr(g3)
V(g3)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g3))))]
#plot(g3, layout=LO, vertex.size=5, vertex.label="")
plot(g3, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g3)
title(paste(length(unique(membership(wc))), 0.95))


#----

g5 <- delete.edges(g, which(E(g)$weight <0.925))
E(g5)$weight = 0.5
LO = layout_with_fr(g5)
V(g5)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g5))))]
plot(g5, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g5)
title(paste(length(unique(membership(wc))), 0.925))

g3 <- delete.edges(g, which(E(g)$weight <0.9))
E(g3)$weight = 0.5
LO = layout_with_fr(g3)
V(g3)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g3))))]
plot(g3, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g3)
title(paste(length(unique(membership(wc))), 0.9))


g4 <- delete.edges(g, which(E(g)$weight <0.85))
E(g4)$weight = 0.5
LO = layout_with_fr(g4)
V(g4)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g4))))]
plot(g4, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g4)
title(paste(length(unique(membership(wc))), 0.85))


g5 <- delete.edges(g, which(E(g)$weight <0.80))
E(g5)$weight = 0.5
LO = layout_with_fr(g5)
V(g5)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g5))))]
plot(g5, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g5)
title(paste(length(unique(membership(wc))), 0.80))

g3 <- delete.edges(g, which(E(g)$weight <0.75))
E(g3)$weight = 0.5
LO = layout_with_fr(g3)
V(g3)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g3))))]
plot(g3, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g3)
title(paste(length(unique(membership(wc))), 0.75))





#------------------------------------------------------------------

ggtree(tree)+geom_tiplab(size=3)


#---------------------------------------
par(mfrow=c(1,1))

#	Select only Marinobacter, and the odd one Tamilnaduibacter
#	Can I find evidence that this is a problem?
filteredANIb = ANIb[grepl('Marinobacter|Tamilnaduibacter',rownames(ANIb)),grepl('Marinobacter|Tamilnaduibacter',colnames(ANIb))]

Tamilna_vs_Marb = filteredANIb[grepl('Tamilnaduibacter',rownames(filteredANIb)),!grepl('Tamilnaduibacter',colnames(filteredANIb))]

g <- graph_from_adjacency_matrix(as.matrix(filteredANIb), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <0.75))
plot(g2)
wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2, vertex.label="", vertex.size=5)
plot(wc, g2,  vertex.size=5, layout=LO)


#-------------------------------------------------------------------

filteredANIb = ANIb[grepl('Marinobacter|Tamilnaduibacter',rownames(ANIb)),grepl('Marinobacter|Tamilnaduibacter',colnames(ANIb))]

Tamilna_vs_Marb = filteredANIb[grepl('Tamilnaduibacter',rownames(filteredANIb)),!grepl('Tamilnaduibacter',colnames(filteredANIb))]

g <- graph_from_adjacency_matrix(as.matrix(filteredANIb), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <0.76))
plot(g2)
wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2, vertex.label="", vertex.size=5)
plot(wc, g2,  vertex.size=5)



#-------------------------------------------------------------------







sim16s[sim16s>95]
msim16s = melt(sim16s)


ggplot(msim16s, aes(x=value)) + geom_histogram(color="black", fill="goldenrod",binwidth=1) +scale_y_continuous(expand = c(0, 0))+ theme_classic()

longFMT16S = data.frame(col=colnames(sim16s)[col(sim16s)], row=rownames(sim16s)[row(sim16s)], dist=c(as.matrix(sim16s)))

selectedGenus = longFMT16S[grepl('Vibrio', longFMT16S$row) & grepl('Vibrio', longFMT16S$col) & longFMT16S$row != longFMT16S$col,]
ggplot(selectedGenus, aes(x=dist)) + geom_histogram(color="black", fill="goldenrod",binwidth=1) +scale_y_continuous(expand = c(0, 0))+ theme_classic()

plot(sim16s[,1],sim16s[,2])

notItself16S = longFMT16S[longFMT16S$row != longFMT16S$col,]
notItself16S

notItself16S[!notItself16S$dist>=99,]

dim(longFMT16S)
sqrt(nrow(longFMT16S[longFMT16S $dist<40,]))


longFMT16S[longFMT16S$row == 'Marinobacter_sp_FDB33' & grepl('Marinobacter', longFMT16S$col),]
longFMT16S[longFMT16S$row == 'Vibrio_sp_IJ49' & grepl('Vibrio', longFMT16S$col),]

notItself16S$GenusCol =gsub('_.*$','', as.character(notItself16S$col))
notItself16S$GenusRow =gsub('_.*$','', as.character(notItself16S$row))

notItself16S$sameGenus = ifelse(notItself16S$GenusCol==notItself16S$GenusRow , TRUE, FALSE)

ggplot(notItself16S, aes(x=dist)) + geom_histogram(color="black", aes(fill= sameGenus),binwidth=1) +scale_y_continuous(expand = c(0, 0))+ theme_classic()+scale_fill_manual(values=c('darkolivegreen4 ','goldenrod'))

TaxTable = read.table('~/DATA/MarinobacterGenomics/miscl/TaxTable.txt')

TaxTable = TaxTable[TaxTable$Genus_name %in%notItself16S$GenusCol,]
rownames(TaxTable)= TaxTable$Genus_name
TaxTable[notItself16S$GenusCol,]

notItself16S$FamilyCol = TaxTable[notItself16S$GenusCol,'Family']
notItself16S$FamilyRow = TaxTable[notItself16S$GenusRow,'Family']
notItself16S$sameFamily = ifelse(notItself16S$FamilyCol ==notItself16S$FamilyRow , TRUE, FALSE)

notItself16S$OrderCol = TaxTable[notItself16S$GenusCol,'Order']
notItself16S$OrderRow = TaxTable[notItself16S$GenusRow,'Order']
notItself16S$sameOrder = ifelse(notItself16S$OrderCol ==notItself16S$OrderRow , TRUE, FALSE)

notItself16S$ClassCol = TaxTable[notItself16S$GenusCol,'Class']
notItself16S$ClassRow = TaxTable[notItself16S$GenusRow,'Class']
notItself16S$sameClass = ifelse(notItself16S$ClassCol ==notItself16S$ClassRow , TRUE, FALSE)

notItself16S$PhylumCol = TaxTable[notItself16S$GenusCol,'Phylum']
notItself16S$PhylumRow = TaxTable[notItself16S$GenusRow,'Phylum']
notItself16S$samePhylum = ifelse(notItself16S$PhylumCol ==notItself16S$PhylumRow , TRUE, FALSE)


notItself16S$taxGrouping = ifelse(notItself16S$sameGenus==TRUE,'Genus',ifelse(notItself16S$sameFamily==TRUE,'Family',ifelse(notItself16S$sameOrder==TRUE,'Order',ifelse(notItself16S$sameClass==TRUE,'Class',ifelse(notItself16S$samePhylum==TRUE,'Phylum','nope')))))

ggplot(notItself16S, aes(x=dist,fill= taxGrouping)) + geom_histogram(alpha=0.5,color="black",binwidth=1,position="identity") +scale_y_continuous(expand = c(0, 0))+ theme_classic()+scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(6))

ggplot(notItself16S, aes(y=dist,x= reorder(taxGrouping,-dist))) +geom_jitter(aes(color= taxGrouping))+ geom_boxplot(color="black", aes(fill= taxGrouping)) +scale_y_continuous(expand = c(0, 0))+ theme_classic()


notItself16S$taxGrouping = ordered(notItself16S$taxGrouping, levels = c("Genus", "Family", "Order",'Class',"Phylum",'nope'))

ggplot(notItself16S, aes(y=dist,x= taxGrouping )) +geom_jitter(aes(color= taxGrouping))+ geom_boxplot(color="black", aes(fill= taxGrouping)) +scale_y_continuous(expand = c(0, 0))+ theme_classic()+scale_fill_manual(values=gray.colors(6, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))+scale_color_manual(values=gray.colors(6, start = 0.4, end = 0.8, gamma = 2.2, alpha = NULL))

ggplot(notItself16S, aes(y=dist,x= taxGrouping )) +geom_jitter(aes(color= taxGrouping),alpha=0.2)+ geom_boxplot(color="black", aes(fill= taxGrouping)) +scale_y_continuous(expand = c(0, 0))+ theme_classic()+scale_fill_manual(values= colorRampPalette(brewer.pal(9, "Set1"))(6))+scale_color_manual(values= colorRampPalette(brewer.pal(9, "Set1"))(6))




library(igraph)
#full network isolates
sim16s2 =  sim16s[!grepl('Marinobacter_sp_FDB33|Alteromonas_sp_FDB36',rownames(sim16s)),!grepl('Marinobacter_sp_FDB33|Alteromonas_sp_FDB36',colnames(sim16s))]
#those from bioassay
sim16s2 =  sim16s[rownames(trait.binairy),rownames(trait.binairy)]


g <- graph_from_adjacency_matrix(as.matrix(sim16s2), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight <95))
plot(g2)

wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)



par(mfrow=c(2,4))
g5 <- delete.edges(g, which(E(g)$weight <92.5))
E(g5)$weight = 0.5
LO = layout_with_fr(g5)
V(g5)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g5))))]
plot(g5, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g5)
title(paste(length(unique(membership(wc))),92.5))

E(g2)$weight = 0.5
LO = layout_with_fr(g2)
V(g2)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g2))))]
plot(g2, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g2)
title(paste(length(unique(membership(wc))),95))

g3 <- delete.edges(g, which(E(g)$weight <97.5))
E(g3)$weight = 0.5
LO = layout_with_fr(g3)
V(g3)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g3))))]
plot(g3, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g3)
title(paste(length(unique(membership(wc))),97.5))


g4 <- delete.edges(g, which(E(g)$weight <99))
E(g4)$weight = 0.5
LO = layout_with_fr(g4)
V(g4)$color <- colorRampPalette(brewer.pal(9, "Set1"))(38)[as.factor(gsub('_.*$','', names(V(g4))))]
plot(g4, layout=LO, vertex.size=5, vertex.label="")
wc <- cluster_walktrap(g4)
title(paste(length(unique(membership(wc))),99))

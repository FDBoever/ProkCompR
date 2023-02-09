# PRECISION RECALL
#=================================================================================#

#usage
multi.phylo.cluster <- multi.phylo.cluster %>% rownames_to_column(var='taxa')
#MonoPhy::AssessMonophyly(tree= lsTrees.ingroup[[6]],taxonomy=multi.phylo.cluster)

taxonomy= multi.phylo.cluster[,c('taxa','concat.sco.nuc')]
tree = lsTrees.ingroup[[1]]
#tree = lsTrees.ingroup[[3]]

check.monoph(tree,taxonomy)

mdt <- taxonomy
rownames(mdt) <- mdt$taxa

p <- ggtree(tree)
p$data <- p$data %>% left_join(mdt, by=c('label'='taxa'))
p+geom_tiplab(size=2)+geom_point(data=p$data[p$data$isTip==TRUE,],aes(x=x,y=y,fill=concat.sco.nuc),shape=21)


p+geom_point(data=p$data[p$data$isTip==TRUE,],aes(x=x,y=y,fill=concat.sco.nuc),shape=21)

#=================================================================================#
#=================================================================================#

internalNodes$stNode = as.numeric(internalNodes$stNode)
p = ggtree(tr)
p$data <- p$data %>% dplyr::left_join(internalNodes, by=c('node'='stNode'))
q <- p + geom_point(data=p$data[p$data$isTip==FALSE,],aes( x,y,fill=RED,size=nodeDepths),shape=21) +
  scale_fill_gradientn(colours = brewer.pal(11, 'RdBu'))
q

internalNodes %>% ggplot(aes(nodeDepths,RED))+
  geom_smooth(method='lm',size=.5,color='black')+
  geom_point(aes(fill=RED),shape=21,size=2,alpha=.7)+
  scale_fill_gradientn(colours = brewer.pal(11, 'RdBu'))+
  fdb_style()


#=================================================================================#

#precision recall mining of positions in tree
Edges = tree$edge
Nodes = unique(Edges[,1])

#Calculating subtrees for each internal node
st = ape::subtrees(tree)
#naming the nodes according to parent tree
names(st)= Nodes

#---	STEP 1 ------
#	to speed up the process, we first store all descending tips in a list per node number
lsNodes2Tips = list()
for(i in Nodes){
  lsNodes2Tips[[as.character(i)]] = st[[as.character(i)]]$tip.label
}

#---	STEP 2 ------
#	generate a count matrix (nodes x taxa) with number of descending nodes that contain with name

#First we generate a list containing all unique taxonomic names (irrespective of RANK)
#allTaxa = unique(unlist(str_split(SINAclassification$lca_tax_slv,';')))
allTaxa = as.character(unique(taxonomy[,2]))

DFdescTipName = c()
for(i in names(lsNodes2Tips)){
  descTipName = c()
  for(taxon in allTaxa){
    descTipName = c(descTipName , table(as.character(taxonomy[taxonomy$taxa %in% st[[as.character(i)]]$tip.label,][,2] == taxon))['TRUE'])

  }
  DFdescTipName = rbind(DFdescTipName, descTipName )
}
rownames(DFdescTipName)=names(lsNodes2Tips)
colnames(DFdescTipName) = allTaxa
head(DFdescTipName)

#---	STEP 3 ------
#	total number of tips under each node
tipsPerNode = lengths(lsNodes2Tips)

#---	STEP 4 ------
#	count tips in total tree that have taxon
tipsInTree = c()
for(taxon in allTaxa){
  tipsInTree = c(tipsInTree , table(taxonomy[,2]==taxon)['TRUE'])
}
names(tipsInTree) = allTaxa

#---	STEP 5 ------
#	Caulculate Precision & Recall
# 	Precision defined as the fraction of informative tips with a given name under a given node out of the total count of informative tips under the node;
#	Recall is the fraction of informative tips descending from a given node out of all the informative tips of the entire tree containing the same name

dfPrecision = c()
dfRecall = c()
for(taxon in allTaxa[allTaxa != ""]){
  dfPrecision = cbind(dfPrecision , (DFdescTipName[, taxon] / tipsPerNode ))
  dfRecall = cbind(dfRecall , (DFdescTipName[, taxon] / tipsInTree[taxon] ))
}

#---	STEP 6 ------
#	Calculate F-1
#	The F-measure (F =  2 *􏰁 ((precision *􏰁 recall)/(precision + recall))) (van Rijsbergen, 1979)
#	is then calculated for each internal node for each name at each taxonomic rank in order to determine the optimal internal node for a name.

dfFmeasure = (2* (dfPrecision * dfRecall ))/(dfPrecision + dfRecall)
colnames(dfFmeasure)=allTaxa[allTaxa != ""]

head(dfFmeasure)
head(dfPrecision)
head(dfRecall)


#library(castor)
#install.packages('castor',dependencies=T)
#-----------------------

#nodeDepths: 	Given a rooted phylogenetic tree, calculate the phylogenetic depth of each node (mean distance to its descending tips).

internalNodes = data.frame(cbind('stNode'=names(st),
                                 'nrTips'= castor::count_tips_per_node(tree),
                                 'RED'=castor::get_reds(tree),
                                 'nodeDepths'=castor::get_all_node_depths(tree)),
                           stringsAsFactors=FALSE)


internalNodes$RED = as.numeric(internalNodes$RED)
internalNodes$nodeDepths = as.numeric(internalNodes$nodeDepths)

str(internalNodes)



#----------------
#for each unique name, the internal node with the highest F-measure score for each name is retained.

dfTax2nodes = c()
for(taxon in colnames(dfFmeasure)){

  setectedTaxonF = sort(dfFmeasure[,taxon],decreasing=TRUE)

  #If a tie is encountered, the internal node with the fewest tips is kept.
  if(table(setectedTaxonF == setectedTaxonF[1])['TRUE'] > 1){

    names(setectedTaxonF[setectedTaxonF == setectedTaxonF[1]])
    selNds = internalNodes[internalNodes$stNode %in% names(setectedTaxonF[setectedTaxonF == setectedTaxonF[1]]), ]
    nde = selNds[order(selNds$nrTips),][1,'stNode']
    print(paste(taxon,' ---> tie is encountered, the internal node with the fewest tips is kept'))
  }else{
    nde= names(sort(dfFmeasure[,taxon],decreasing=TRUE)[1])
 }
  dfTax2nodes = rbind(dfTax2nodes,c('taxon'=taxon, 'node'=names(sort(dfFmeasure[,taxon],decreasing=TRUE)[1]),'F'=as.numeric(sort(dfFmeasure[,taxon],decreasing=TRUE)[1])))
}

dfTax2nodes=data.frame(dfTax2nodes,stringsAsFactors=FALSE)
dfTax2nodes$F = as.numeric(dfTax2nodes$F )
head(dfTax2nodes)

#To view a spefic taxon, one can use grep
#	dfTax2nodes[grepl('Pavlovo',dfTax2nodes$taxon),]

#----------------
library(dplyr)
dfPerTaxon = dplyr::right_join(internalNodes, dfTax2nodes, by=c('stNode'= 'node'))
dfPerNode = dplyr::full_join(internalNodes, dfTax2nodes, by=c('stNode'= 'node'))
head(dfPerNode)


#dfTax2nodesUnique = aggregate(as.character(taxon) ~ stNode + nrTips + RED + nodeDepths + F , data = dfPerNode, c)
mergedTaxa = aggregate(taxon ~ stNode + F  , data = dfPerNode, c)
mergedTaxa  = dplyr::right_join(internalNodes, mergedTaxa, by='stNode')

head(mergedTaxa)

#------- node converter to link back to all trees
# using findMCRA function, all descendants from a node

nodeconverter =  c()
for(node in mergedTaxa$stNode){
  nodeconverter = rbind(nodeconverter, c('st'=node,'MRCA'=findMRCA(tree, tips= st[[node]]$tip.label, type=c("node"))))
}
head(nodeconverter)

mergedTaxa$MCRA = as.numeric(as.character(data.frame(nodeconverter)$MRCA))

#------#


tr=tree

#tr$tip.label = paste(tr$tip.label,SINAclassification[tr$tip.label, 'lca_tax_slv'], sep=' - ')

p = ggtree(tr)
tracker=c()
Ftracker = c()
for(i in p$data$node){
  tracker=c(tracker ,as.character(mergedTaxa[mergedTaxa $MCRA==i,][1,'taxon']))
  Ftracker = c(Ftracker ,(mergedTaxa[mergedTaxa $MCRA==i,'F'][1]))
}
p$data$FmeasureTaxon = tracker
p$data$F = Ftracker

q <- p + geom_label_repel(data = p$data[!(p$data$FmeasureTaxon=='NULL'),],aes( label= FmeasureTaxon))+geom_point(data = p$data[!(p$data$FmeasureTaxon=='NULL'),],aes( x,y,color=F))
q

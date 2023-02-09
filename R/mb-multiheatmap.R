
library(tidyr)
library(randomForest)
library(gplots)
library(ggplot2)
library(dplyr)
library(ape)
library(phytools)
library(reshape)
library(phylolm)
#########################################
#	COLORS
#########################################

#---- MANUALLY ADD COLORS AS YOU PLEASE --------------
# assess the grouping factor in the ANVIO_cat file and sort the colors on the basis of
# levels(ANVIO_cat_restricted$group2)

colorssorted = c(
  rgb(57,129,29, maxColorValue=255),#adhaerens
  rgb(166,216,212, maxColorValue=255),#ant
  rgb(118,174,207, maxColorValue=255),#psych
  rgb(242,236,112, maxColorValue=255), #hydrocarb
  rgb(242,169,104, maxColorValue=255),#ex
  rgb(179,77,34, maxColorValue=255),#vinifirmus
  rgb(134,218,129, maxColorValue=255),#alg
  rgb(227,143,221, maxColorValue=255),#lipo
  rgb(131,115,27, maxColorValue=255)#sedim
)

colors = c(
  rgb(57,129,29, maxColorValue=255),#adhaerens
  rgb(134,218,129, maxColorValue=255),#alg
  rgb(166,216,212, maxColorValue=255),#ant
  rgb(242,169,104, maxColorValue=255),#ex
  rgb(242,236,112, maxColorValue=255), #hydrocarb
  rgb(227,143,221, maxColorValue=255),#lipo
  rgb(118,174,207, maxColorValue=255),#psych
  rgb(131,115,27, maxColorValue=255),#sedim
  rgb(179,77,34, maxColorValue=255)#vinifirmus
)

Habitat_colors = c("#E41A1C","#FFD422","#43997A","#658E67","#5D6795","#A35390","#B85F49")

##############

ANVIO <- read.table('~/DATA/MarinobacterGenomics/2018_ProkComp/MARINOBACTER_protein_clusters_summary.txt',header=TRUE,sep="\t") #protein clusters
ANVIO $COG_FUNCTION = sub('\\\t.*', '', ANVIO$COG_FUNCTION)
COG_list = c("D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
sum=0
mat = t(table(droplevels(ANVIO[grep("C", ANVIO $COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
CpcogFull = data.frame(rowSums(mat))
for(i in COG_list){
  mat = t(table(droplevels(ANVIO[grep(i, ANVIO$COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
  CpcogFull = cbind(CpcogFull,rowSums(mat))
}
names(CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X")


COG_list = unique(unlist(strsplit(as.character(unique(ANVIO$COG_FUNCTION_ACC)),'\\|')))
#sum=0
#mat = t(table(droplevels(ANVIO[grep(COG_list[1], ANVIO $COG_FUNCTION_ACC),c('COG_FUNCTION_ACC','genome_name')])))
#CpcogFull = data.frame(rowSums(mat))

mat = t(table(droplevels(ANVIO[grep("C", ANVIO $COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
CpcogFull = data.frame(rowSums(mat))

for(i in COG_list){
  mat = t(table((ANVIO[grep(i, ANVIO$COG_FUNCTION_ACC),c('COG_FUNCTION_ACC','genome_name')])))
  CpcogFull = cbind(CpcogFull,rowSums(mat))
}
CpcogFull = CpcogFull[,2:ncol(CpcogFull)]
colnames(CpcogFull)= COG_list






##############



pruned_tree = read.tree('~/Enrichment_test/pruned_tree')

domains.table <- read.delim(file='~/DATA/MarinobacterGenomics/2018_ProkComp/TCDB/TCDB_melt.tsv',header=T,sep='\t')
domains.wide = spread(domains.table,Genome,Abundance,fill=0)
domains.wide = domains.wide[(-1),]
write.table(domains.wide, file='TCDB_spread_matrix.tsv',sep='\t',row.names=F, quote=F)

head(domains.wide)



#domains are set as rownames
rownames(domains.wide)= domains.wide[,1]
domains.wide = domains.wide[,2:ncol(domains.wide)]

#we remove all rows which sum up to 0
domains.wide = domains.wide[ rowSums(domains.wide)!=0, ]

#we remove all singletons
domains.wide = domains.wide[ rowSums(domains.wide)!=1, ]




metadata  = read.table("~/Enrichment_test/ANVIO_cat3.txt")


heatmap.2(as.matrix(domains.wide), trace='none',margin=c(10,10),col = colorRampPalette(c('white',"darkblue","black",'red','yellow'))(40))


colnames(domains.wide) = gsub('-','_',colnames(domains.wide))
domains.wide = domains.wide[,pruned_tree$tip.label]

#we remove all rows which sum up to 0
domains.wide = domains.wide[ rowSums(domains.wide)!=0, ]

#we remove all singletons
domains.wide = domains.wide[ rowSums(domains.wide)!=1, ]

data = t(domains.wide)
data = data[rownames(metadata),]

combined_data=cbind('total'=rowSums(data),metadata)
ggplot(combined_data,aes(x=group2,y=total))+geom_boxplot()+geom_jitter()+theme_classic()
ggplot(combined_data,aes(x=SS,y=total))+geom_boxplot()+geom_jitter()+theme_classic()

######
obj <- df2genind(data, ncode=1)
obj
tab(obj)
pop(obj) <- metadata$group2
dapc.H3N2 <- dapc(obj, var.contrib = TRUE, scale = FALSE, n.pca = 5, n.da = nPop(obj) - 1)
scatter(dapc.H3N2)
contrib <- loadingplot(dapc.H3N2$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
contrib

dapc.H3N2 <- dapc(obj, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(obj) - 1)
scatter(dapc.H3N2)
contrib <- loadingplot(dapc.H3N2$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
contrib




mod = varpart(data, ~ lifestyle, cophenetic(pruned_tree), data= data.frame(metadata), transfo="hel")



######
pruned_tree2 = drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, rownames(metadata)))

data= CpcogFull
data= PfamData

dataGLM = metadata[,c('group2','SS')]
dataGLM=cbind(dataGLM,data[rownames(dataGLM),])

Results=c()
resTot=c()
for(predictor in as.character(unique(dataGLM$SS))){
  print(predictor)
  dataGLM$selected = ifelse(dataGLM$SS == predictor, 1, 0)
  for(category in colnames(data)){
    selectedDATA = dataGLM[,c('selected', category )]
    colnames(selectedDATA)=c('predictor','trait')
    #m1 = phyloglm(trait~predictor,phy= pruned_tree2, data= selectedDATA, method='logistic_IG10')
    m1 = phylolm(trait~predictor,phy= pruned_tree2, data= selectedDATA, method='BM')

    m1.sum <- summary(m1)
    res <- data.frame(orthogroup.id = category,
                      Estimate = m1.sum$coefficients["predictor",1],
                      SE = m1.sum$coefficients["predictor",2],
                      z.value = m1.sum$coefficients["predictor",3],
                      p.value = m1.sum$coefficients["predictor",4],
                      predictor = predictor)
    resTot = rbind(resTot,res)
  }
}

Pvalues= spread(resTot[,c('predictor','orthogroup.id','p.value')], orthogroup.id,p.value)
rownames(Pvalues) = Pvalues$predictor
Pvalues = Pvalues[,2:ncol(Pvalues)]

adjustedPvalues=sapply(Pvalues, function(x) p.adjust(x, method = "BH", n = nrow(Pvalues)) )
rownames(adjustedPvalues) = rownames(Pvalues)

thresholded = adjustedPvalues<0.05
thresholded = ifelse(thresholded == TRUE, 1, 0)



Estimates = spread(resTot[,c('predictor','orthogroup.id','Estimate')], orthogroup.id,Estimate)
rownames(Estimates) = Estimates $predictor
Estimates = Estimates[,2:ncol(Estimates)]

cleanedEstimates = Estimates* thresholded
cleanedEstimates[is.na(cleanedEstimates)]=0
#heatmap.2(as.matrix(cleanedEstimates),trace='none')
#heatmap.2(as.matrix(cleanedEstimates[,colSums(abs(cleanedEstimates))>0]),trace='none',col= colorRampPalette(c("blue", "white", "red"))(n = 21))
#heatmap.2(as.matrix(cleanedEstimates[,colSums(abs(cleanedEstimates))>2]),trace='none',col= colorRampPalette(c("blue", "white", "red"))(n = 21))

FinalRes = c()
for(ortho in unique(resTot$orthogroup.id)){
  FinalRes = rbind(FinalRes,cbind(resTot[resTot$orthogroup.id == ortho, ],'p.adj_FDB'= p.adjust(resTot[resTot$orthogroup.id == ortho, 'p.value'], method = "BH", n = nrow(Pvalues))))
}



write.table(na.omit(FinalRes), file='~/PhyloLM_CAZY_habitat.txt')
write.table(na.omit(FinalRes), file='~/PhyloLM_TCDB_habitat.txt')
write.table(na.omit(FinalRes), file='~/PhyloLM_Pfam_habitat.txt')
write.table(na.omit(FinalRes), file='~/PhyloLM_OGs_habitat.txt')
write.table(na.omit(FinalRes), file='~/PhyloLM_COGcat_habitat.txt')

write.table(na.omit(FinalRes), file='~/PhyloLM_COGs_habitat.txt')

dataGLM$selected = ifelse(dataGLM$SS == 'Phototroph', 1, 0)


selectedFinal = FinalRes[FinalRes$predictor=='Phototroph' & FinalRes$ p.adj_FDB < 0.05,]
selectedTOP = head(selectedFinal[order(-selectedFinal$Estimate),],20)$orthogroup.id

melt(dataGLM[,c('selected','COG0410','COG1475','COG3284')],id='selected')


selectedTable = melt(dataGLM[,c('selected', as.character(selectedTOP))],id='selected')

ggplot(selectedTable,aes(as.factor(selected), value))+geom_boxplot()+geom_jitter(width = 0.1, height = 0)+theme_classic()+facet_wrap(~ variable,scales="free")


selectedTable = melt(dataGLM[,c('selected','COG0177','COG3141','COG1726','COG1981','COG0605','COG4536')],id='selected')
ggplot(selectedTable,aes(as.factor(selected), value))+geom_boxplot()+geom_jitter(width = 0.5, height = 0)+theme_classic()+facet_wrap(~ variable,scales="free")



######

######################################################################################
#
#	Manhattan - MDS (PcoA)
#
######################################################################################

data.manh <- dist(data,method = "manhattan")

mds.cmdscale <- as.data.frame(cmdscale(as.matrix(data.manh)))
mds.cmdscale$names <- rownames(mds.cmdscale)


ggplot(mds.cmdscale, aes(V1, V2, label=names)) +
  geom_point(aes(colour=factor(as.character(metadata$group2))), size=2.3) +
  geom_text(aes(colour=factor(as.character(metadata$group2))), check_overlap = TRUE, size=2.2,
            hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.005) +
  scale_colour_discrete(name = "groups") +
  labs(x="", y="", title="MDS by Manhattan and cmdscale()") + theme_bw() + scale_color_manual(values= colors)

ggplot(mds.cmdscale, aes(V1, V2, label=names)) +
  geom_point(aes(colour=factor(as.character(metadata$group2))), size=2.3) +
  scale_colour_discrete(name = "groups") +
  labs(x="", y="", title="MDS by Manhattan and cmdscale()") + theme_bw()+ scale_color_manual(values= colors)


ggplot(mds.cmdscale, aes(V1, V2, label=names)) +
  geom_point(aes(colour=factor(as.character(metadata$SS))), size=2.3) +
  scale_colour_discrete(name = "groups") +
  labs(x="", y="", title="MDS by Manhattan and cmdscale()") + theme_bw()+ scale_color_manual(values= Habitat_colors)



######################################################################################
#
#	RANDOM FOREST
#	The Random Forests (RF) [Breiman 2001] algorithm is an increasingly popular machine learning algorithm within statistical genetics.
#	RF is among the moost accurate classification methods to date and is suitable to analyse high dimentional nxp matrices
#	using the randomForest() function in the 'randomForest() package
#
#
#	Reference:
#		doi:  10.2202/1544-6115.1691
#
#######################################################################################


set.seed(42)

#prepare the dataset suitable for randomForest()
RF_data = data.frame(data,'group'=as.character(metadata$group2))
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$group2))

#TRAIN THE RF (look for code in RandomForest_Brieuc_modification.R)
RF_1 <- randomForest(group ~ ., data= RF_data,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 100 ,strata= group)


# Extracts variable importance (Mean Decrease in Gini Index)
# Sorts by variable importance and relevels factors to match ordering
RF_var_importance <- data_frame(variable=setdiff(colnames(RF_data), "group"),
                                importance=as.vector(importance(RF_1,type=2)))

RF_var_importance  = data.frame(RF_var_importance)
str(RF_var_importance)
RF_var_importance <- arrange(RF_var_importance, desc(importance))
RF_var_importance$variable <- factor(RF_var_importance $variable, levels= RF_var_importance $variable)

head(RF_var_importance)
RF_var_importance[1:10,]
p <- ggplot(RF_var_importance, aes(x=variable, weight=importance))
p <- p + geom_bar(color="darkgreen") + ggtitle("Variable Importance from Random Forest Fit")
p <- p + xlab("Protein Clusters") + ylab("Variable Importance (Mean Decrease in Gini Index)")
p <- p +theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p + geom_vline(xintercept = 50,colour="red", linetype="dashed") + scale_y_continuous(expand = c(0, 0))


q <- ggplot(RF_var_importance[1:100,], aes(x=variable, weight=importance))
q <- q + geom_bar(color="darkgreen") + ggtitle("Variable Importance from Random Forest Fit")
q + coord_flip()



RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:50,]$variable)]
RFdata_selected = RFdata_selected[,!colnames(RFdata_selected) %in% c("The.ATP.binding.Cassette..ABC..Superfamily",'The.Tripartite.ATP.independent.Periplasmic.Transporter..TRAP.T..Family','The.H.sup....sup...or.Na.sup....sup..translocating.Bacterial.Flagellar.Motor.ExbBD.Outer.Membrane.Transport.Energizer..Mot.Exb..Superfamily','')]

head(RFdata_selected)
as.character(colnames(data))



pruned_tree$edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(RFdata_selected))
combined_ordered_matrix <- RFdata_selected[new_order,]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
          ,col = colorRampPalette(c('white',"darkblue","black",'red','yellow'))(100),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none",col=colorRampPalette(brewer.pal(9, "Greens"))(100),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)


testMolten = melt(as.matrix((data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)])))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
ggplot(testMolten,aes(x=group2,y=value,color=group2))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values=colors)
ggplot(testMolten,aes(x=SS,y=value,color=SS))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values= Habitat_colors)


heatmap.2(t(dfsorted[,1:20]),trace='none',margin=c(20,20),col = colorRampPalette(c('white',"darkblue","black",'red','yellow'))(40))

#testMolten = melt(data[,grepl(paste('RND','MFS','ATP','ABC','OOP','Amt','Mot','MscS',sep='|'),colnames(data))])
#testMolten = melt(data[,grepl(paste('ABC','DNA','IIISP','IVSP','MTB',sep='|'),colnames(data))])
#testMolten = melt(data[,grepl(paste('Amt','Mot','MscS','TRIC','VIC',sep='|'),colnames(data))])
#testMolten = melt(as.matrix((data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)])))
#testMolten = melt(as.matrix((data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)])))
#dfsorted = data[,order(colSums(-data,na.rm=TRUE))]
#dfsorted[,1:10]
#testMolten = melt(dfsorted[,1:25])

testMolten = melt(as.matrix((data.frame(data)[,c("The.ATP.binding.Cassette..ABC..Superfamily",'The.Tripartite.ATP.independent.Periplasmic.Transporter..TRAP.T..Family','The.H.sup....sup..or.Na.sup....sup..translocating.NADH.Dehydrogenase..NDH..Family ','The.Proton.translocating.Cytochrome.Oxidase..COX..Superfamily', 'The.Fatty.Acid.Transporter..FAT..Family','The.Monovalent.Cation.Proton.Antiporter')])))
RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:28,]$variable)]


testMolten = melt(as.matrix((data.frame(data)[,names(head(sort(-colSums(RFdata_selected)),10))])))

colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')

totalsDF = data.frame(rowSums(data))
totalsDF$name = rownames(totalsDF)
colnames(totalsDF) = c('total','name')

testMolten = right_join(testMolten, totalsDF,by='name')


testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))

br1 = ggplot(testMolten,aes(x=group2,y=value,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)

br1 = ggplot(testMolten,aes(x=group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)


testMolten = melt(as.matrix(data[,names(head(sort(-colSums(data)),10))]))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
totalsDF = data.frame(rowSums(data))
totalsDF$name = rownames(totalsDF)
colnames(totalsDF) = c('total','name')
testMolten = right_join(testMolten, totalsDF,by='name')
testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))

br1 = ggplot(testMolten,aes(x=group2,y=value,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)

br1 = ggplot(testMolten,aes(x=group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)




head(testMolten)
tot1 = ggplot(testMolten,aes(x=group2,y=total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)+ theme(legend.position = "none")
tot2 = ggplot(testMolten,aes(x= SS,y=total,fill= SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)+ theme(legend.position = "none")

multiplot(tot1,tot2,cols=2)



CheckM2_annotated2 = CheckM2_annotated
rownames(CheckM2_annotated2)= CheckM2_annotated2$Bin_Id
totalsDF

totalsDFannotated = cbind(totalsDF ,CheckM2_annotated2[rownames(totalsDF),])
ggplot(totalsDFannotated,aes(x= total_length/1000000,y=total))+geom_point()+geom_smooth(method='lm') +geom_point(aes(color=group2)) +scale_color_manual(values= colors)

####################

RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$SS))
RF_1 <- randomForest(group ~ ., data= RF_data,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 100 ,strata= group)

RF_var_importance <- data_frame(variable=setdiff(colnames(RF_data), "group"),
                                importance=as.vector(importance(RF_1,type=2)))

RF_var_importance  = data.frame(RF_var_importance)
str(RF_var_importance)
RF_var_importance <- arrange(RF_var_importance, desc(importance))
RF_var_importance$variable <- factor(RF_var_importance $variable, levels= RF_var_importance $variable)

RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)]
RFdata_selected = RFdata_selected[,!colnames(RFdata_selected) %in% c("The.ATP.binding.Cassette..ABC..Superfamily")]

combined_ordered_matrix <- RFdata_selected[new_order,]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
          ,col=colorRampPalette(brewer.pal(9, "BuPu"))(25),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)

testMolten = melt(as.matrix((data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)])))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
#ggplot(testMolten,aes(x=group2,y=value,color=group2))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values=colors)
ggplot(testMolten,aes(x=SS,y=value,color=SS))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values= Habitat_colors)


####################

RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$lifestyle))
RF_1 <- randomForest(group ~ ., data= RF_data,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 100 ,strata= group)

RF_var_importance <- data_frame(variable=setdiff(colnames(RF_data), "group"),
                                importance=as.vector(importance(RF_1,type=2)))

RF_var_importance  = data.frame(RF_var_importance)
str(RF_var_importance)
RF_var_importance <- arrange(RF_var_importance, desc(importance))
RF_var_importance$variable <- factor(RF_var_importance $variable, levels= RF_var_importance $variable)

RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)]
RFdata_selected = RFdata_selected[,!colnames(RFdata_selected) %in% c("The.ATP.binding.Cassette..ABC..Superfamily")]

combined_ordered_matrix <- RFdata_selected[new_order,]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
          ,col=colorRampPalette(brewer.pal(9, "BuPu"))(25),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)

testMolten = melt(as.matrix((data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)])))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
ggplot(testMolten,aes(x=group2,y=value,color=group2))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values=colors)
ggplot(testMolten,aes(x=SS,y=value,color=SS))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values= Habitat_colors)




####################
#CAZY
####################
CAZY =read.delim("~/Downloads/dbCAN/dbCAN_parse.txt")

CAZYl = reshape(data= CAZY,idvar="Genome",timevar = "CAZy_domaim",direction="wide")
colnames(CAZYl)= gsub("Abundance.", "", colnames(CAZYl))
CAZYl = CAZYl[,1: 139]
rownames(CAZYl)= CAZYl$Genome
CAZYl = CAZYl[,2:length(colnames(CAZYl))]
head(CAZYl)
CAZYl[is.na(CAZYl)] <- 0

#we remove empty domains
CAZYl = CAZYl[ , colSums(CAZYl)!=0]
rownames(CAZYl) = gsub('-','_', rownames(CAZYl))
#CAZYl = CAZYl[,pruned_tree$tip.label]

data = CAZYl


RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$algicolad))
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$Polards))
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$AlgalAss))



RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$group2))
RF_1 <- randomForest(group ~ ., data= RF_data,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 100 ,strata= group)

RF_var_importance <- data_frame(variable=setdiff(colnames(RF_data), "group"),
                                importance=as.vector(importance(RF_1,type=2)))

RF_var_importance  = data.frame(RF_var_importance)
str(RF_var_importance)
RF_var_importance <- arrange(RF_var_importance, desc(importance))
RF_var_importance$variable <- factor(RF_var_importance $variable, levels= RF_var_importance $variable)

RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:25,]$variable)]
RFdata_selected = RFdata_selected[,!colnames(RFdata_selected) %in% c("GT2","GT4")]

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(RFdata_selected))
combined_ordered_matrix <- RFdata_selected[new_order,]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
          ,col=colorRampPalette(brewer.pal(9, "BuPu"))(28),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)

testMolten = melt(as.matrix((data.frame(data)[,as.character(RF_var_importance[1:50,]$variable)])))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
ggplot(testMolten,aes(x=group2,y=value,color=group2))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values=colors)
ggplot(testMolten,aes(x=SS,y=value,color=SS))+geom_point()+geom_boxplot()+facet_wrap(~domain,scales='free_y')+scale_color_manual(values= Habitat_colors)

rownames(CAZYl)


testMolten = melt(as.matrix((data.frame(data)[,names(head(sort(-colSums(RFdata_selected)),10))])))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
totalsDF = data.frame(rowSums(data))
totalsDF$name = rownames(totalsDF)
colnames(totalsDF) = c('total','name')
testMolten = right_join(testMolten, totalsDF,by='name')
testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))

#br1 = ggplot(testMolten,aes(x=group2,y=value,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
#br2 = ggplot(testMolten,aes(x=SS,y=value,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
#multiplot(br1,br2,cols=1)

br1 = ggplot(testMolten,aes(x=group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)


testMolten = melt(as.matrix(data[,names(head(sort(-colSums(data)),10))]))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
totalsDF = data.frame(rowSums(data))
totalsDF$name = rownames(totalsDF)
colnames(totalsDF) = c('total','name')
testMolten = right_join(testMolten, totalsDF,by='name')
testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))

#br1 = ggplot(testMolten,aes(x=group2,y=value,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
#br2 = ggplot(testMolten,aes(x=SS,y=value,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
#multiplot(br1,br2,cols=1)

br1 = ggplot(testMolten,aes(x=group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)


tCAZYl=t(CAZYl)
cazyGroups = c()
for( i in c('GH','GT','PL','CE','CBM')){
  cazyGroups = cbind(cazyGroups , colSums(tCAZYl[grepl(i,rownames(tCAZYl)),]))
}
colnames(cazyGroups) = c('GH','GT','PL','CE','CBM')
cazyGroups=data.frame(cazyGroups)
cazyGroups$name=rownames(cazyGroups)

testMolten = melt(cazyGroups)



testMolten = right_join(testMolten, metadata,by='name')
testMolten = right_join(testMolten, totalsDF,by='name')
testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))


br1 = ggplot(testMolten,aes(x= group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)+facet_wrap(~ variable,scales='free_y',ncol=5)
br2 = ggplot(testMolten,aes(x= SS,y= value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)+facet_wrap(~ variable,scales='free_y',ncol=5)
multiplot(br1,br2,cols=1)

cazyGroups

###################
# PFAM
###################
#---- PFAM data --------------
PfamData <- read.table('~/DATA/MarinobacterGenomics/2018_ProkComp/PfamScan_output.csv',header=TRUE,sep='\t')
PfamData = table(PfamData$genome, PfamData$domain)
PAPfamData = PfamData
PAPfamData[PAPfamData>1] <- 1

PfamData = as.data.frame.matrix(PfamData)

metadata$algicolad = ifelse(metadata$group2 == "algicola",'yes','no')
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$algicolad))

metadata$Polards = ifelse(metadata$lifestyle == "polar",'yes','no')
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$Polards))

metadata$AlgalAss = ifelse(metadata$lifestyle == "alga",'yes','no')
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$AlgalAss))


selected = c('Cofac','Cytochrom')
RFdata_selected = data[,grepl(paste(selected,collapse='|'),colnames(data))]
combined_ordered_matrix <- RFdata_selected[new_order,]
heatmap.2(as.matrix((combined_ordered_matrix)),trace="none",col=colorRampPalette(brewer.pal(9, "YlOrRd"))(25),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)

RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$group2))
RF_1 <- randomForest(group ~ ., data= RF_data,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 100 ,strata= group)

RF_var_importance <- data_frame(variable=setdiff(colnames(RF_data), "group"),
                                importance=as.vector(importance(RF_1,type=2)))

RF_var_importance  = data.frame(RF_var_importance)
str(RF_var_importance)
RF_var_importance <- arrange(RF_var_importance, desc(importance))
RF_var_importance$variable <- factor(RF_var_importance $variable, levels= RF_var_importance $variable)

RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:50,]$variable)]

#drop out highly ambundant groups
#RFdata_selected = RFdata_selected[,!colnames(RFdata_selected) %in% c("GGDEF","MCPsignal")]



clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(RFdata_selected))
combined_ordered_matrix <- RFdata_selected[new_order,]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
          ,col=colorRampPalette(brewer.pal(9, "YlOrRd"))(15),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)
heatmap.2(as.matrix((combined_ordered_matrix)),trace="none",col=colorRampPalette(brewer.pal(9, "YlOrRd"))(75),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)

testMolten = melt(as.matrix((data.frame(data)[,names(head(sort(-colSums(RFdata_selected)),10))])))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
totalsDF = data.frame(rowSums(data))
totalsDF$name = rownames(totalsDF)
colnames(totalsDF) = c('total','name')
testMolten = right_join(testMolten, totalsDF,by='name')
testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))

#br1 = ggplot(testMolten,aes(x=group2,y=value,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
#br2 = ggplot(testMolten,aes(x=SS,y=value,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
#multiplot(br1,br2,cols=1)

br1 = ggplot(testMolten,aes(x=group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)


testMolten = melt(as.matrix(data[,names(head(sort(-colSums(data)),10))]))
colnames(testMolten)=c('name','domain','value')
testMolten = right_join(testMolten,metadata,by='name')
totalsDF = data.frame(rowSums(data))
totalsDF$name = rownames(totalsDF)
colnames(totalsDF) = c('total','name')
testMolten = right_join(testMolten, totalsDF,by='name')
testMolten $group2  = factor(testMolten $group2,levels=c('adhaerens', 'antarcticus', 'psychro', 'hydrocarbo', 'excellens', 'vinifirmus','algicola', 'lipo','sedim'))

#br1 = ggplot(testMolten,aes(x=group2,y=value,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
#br2 = ggplot(testMolten,aes(x=SS,y=value,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
#multiplot(br1,br2,cols=1)

br1 = ggplot(testMolten,aes(x=group2,y=value/total,fill=group2,color=group2))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= colorssorted)+scale_color_manual(values= colorssorted)
br2 = ggplot(testMolten,aes(x=SS,y=value/total,fill=SS,color= SS))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~domain,scales='free_y',ncol=5)+scale_fill_manual(values= Habitat_colors)+scale_color_manual(values= Habitat_colors)
multiplot(br1,br2,cols=1)


#---- PFAM data --------------
#generate a dataframe that contains only the DUF domains
DUFData = PfamData[,grep("DUF", colnames(PfamData))]
PA_DUFData = DUFData
PA_DUFData[PA_DUFData>1] <- 1
###################
# KEGG
###################




MAPLE= read.delim("~/DATA/MarinobacterGenomics/2018_ProkComp/MAPLE_v1.txt")
rownames(MAPLE)=paste(MAPLE$ID,MAPLE$Name)
head(MAPLE)
MAPLE= MAPLE[,6:ncol(MAPLE)]

MAPLE=t(MAPLE)
dim(MAPLE)

MAPLE  = MAPLE[intersect(rownames(metadata),rownames(MAPLE)),]
data = MAPLE

RF_data = data.frame(data[intersect(rownames(metadata),rownames(MAPLE)),],'group'=as.character(metadata[intersect(rownames(metadata),rownames(MAPLE)),]$group2))
RF_data = data.frame(data[intersect(rownames(metadata),rownames(MAPLE)),],'group'=as.character(metadata[intersect(rownames(metadata),rownames(MAPLE)),]$algicolad))
RF_data = data.frame(data[intersect(rownames(metadata),rownames(MAPLE)),],'group'=as.character(metadata[intersect(rownames(metadata),rownames(MAPLE)),]$Polards))



RF_1 <- randomForest(group ~ ., data= RF_data,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 100 ,strata= group)

RF_var_importance <- data_frame(variable=setdiff(colnames(RF_data), "group"),
                                importance=as.vector(importance(RF_1,type=2)))

RF_var_importance  = data.frame(RF_var_importance)
str(RF_var_importance)
RF_var_importance <- arrange(RF_var_importance, desc(importance))
RF_var_importance$variable <- factor(RF_var_importance $variable, levels= RF_var_importance $variable)

RFdata_selected=data.frame(data)[,as.character(RF_var_importance[1:50,]$variable)]

#drop out highly ambundant groups
#RFdata_selected = RFdata_selected[,!colnames(RFdata_selected) %in% c("GGDEF","MCPsignal")]



clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(RFdata_selected))
combined_ordered_matrix <- RFdata_selected[new_order,]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
          ,col=colorRampPalette(brewer.pal(9, "YlOrRd"))(15),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)
heatmap.2(as.matrix((combined_ordered_matrix)),trace="none",col=colorRampPalette(brewer.pal(9, "YlOrRd"))(75),srtCol=45,cexRow=0.4,cexCol=0.5,main='',margins=c(10,20),Rowv= pruned_tree_d)


metadata$algicolad = ifelse(metadata$group2 == "algicola",'yes','no')
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$algicolad))

metadata$Polards = ifelse(metadata$lifestyle == "polar",'yes','no')
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$Polards))

metadata$AlgalAss = ifelse(metadata$lifestyle == "alga",'yes','no')
RF_data = data.frame(data[rownames(metadata),],'group'=as.character(metadata$AlgalAss))





#################

OGs = read.table('~/Github/ProkComp/data_files/panmatrix.txt')

data = t(OGs)
data = PfamData
data = CAZYl
data = t(domains.wide)
data = MAPLE

PRM_data = data.frame(data[intersect(rownames(metadata),rownames(data)),])
dPRM_data = distManhattan(PRM_data)
#adonis2(dPRM_data ~ group2 * SS, data = metadata[rownames(PRM_data),])
adonis_phylogeny = adonis(dPRM_data ~ group2, data = metadata[rownames(PRM_data),])
adonis_habitat = adonis(dPRM_data ~ SS, data = metadata[rownames(PRM_data),])


anosim_phylogeny = anosim(dPRM_data, metadata[rownames(PRM_data),]$group2)
anosim_phylogeny # take a look at results
summary(anosim_phylogeny)
plot(anosim_phylogeny)
anosim_habitat= anosim(dPRM_data, metadata[rownames(PRM_data),]$SS)
anosim_habitat # take a look at results
summary(anosim_habitat)
plot(anosim_habitat)


dispersion<-betadisper(dPRM_data, group= metadata[rownames(PRM_data),]$group2)
permutest(dispersion)
plot(dispersion, hull=F, ellipse=TRUE) ##sd ellipse
PCoA1=dispersion$vectors[,1]
PCoA2=dispersion$vectors[,2]
PCoA.plot<-cbind(PCoA1, PCoA2, metadata[rownames(PRM_data),])
Habitat.pcoa.p<-ggplot(PCoA.plot, aes(PCoA1, PCoA2))+
  geom_point(position=position_jitter(.1),aes(color=SS))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill=SS))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= Habitat_colors)+scale_fill_manual(values= Habitat_colors)+ theme(legend.position = "none")
Habitat.pcoa.p

average_dispersion_from_centroid = cbind(data.frame(dispersion$distances),metadata)
habitat.disp.p = ggplot(average_dispersion_from_centroid,aes(x=SS, y=dispersion.distances,fill= SS, color=SS))+geom_boxplot(alpha=0.7)+geom_point()+scale_color_manual(values= Habitat_colors)+scale_fill_manual(values= Habitat_colors)+theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "none")



dispersion<-betadisper(dPRM_data, group= metadata[rownames(PRM_data),]$SS)
permutest(dispersion, pairwise=T, permutation=1000)
PCoA1=dispersion$vectors[,1]
PCoA2=dispersion$vectors[,2]
PCoA.plot<-cbind(PCoA1, PCoA2, metadata[rownames(PRM_data),])

phylo.pcoa.p <-ggplot(PCoA.plot, aes(PCoA1, PCoA2))+
  geom_point(position=position_jitter(.1),aes(color=group2))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill= group2))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+ theme(legend.position = "none")
phylo.pcoa.p

average_dispersion_from_centroid = cbind(data.frame(dispersion$distances),metadata)
phylo.disp.p = ggplot(average_dispersion_from_centroid,aes(x= group2, y=dispersion.distances,fill= group2, color= group2))+geom_boxplot(alpha=0.7)+geom_point()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "none")



phyloMDS<-metaMDS(dPRM_data, k=2, trymax=35, autotransform=TRUE) ##k is the number of dimensions
phyloMDS ##metaMDS takes eaither a distance matrix or your community matrix (then requires method for 'distance=')
stressplot(phyloMDS)

NMDS1 <- phyloMDS$points[,1] ##also found using scores(phyloMDS)
NMDS2 <- phyloMDS$points[,2]

mds.plot<-cbind(NMDS1, NMDS2, metadata[rownames(PRM_data),])

#plot ordination
Habitat.p<-ggplot(mds.plot, aes(NMDS1, NMDS2))+
  geom_point(position=position_jitter(.1),aes(color=SS))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill=SS))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= Habitat_colors)+scale_fill_manual(values= Habitat_colors)+ theme(legend.position = "none")
Habitat.p


phylo.p<-ggplot(mds.plot, aes(NMDS1, NMDS2))+
  geom_point(position=position_jitter(.1),aes(color=group2))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill= group2))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+ theme(legend.position = "none")
phylo.p

multiplot(phylo.p, phylo.pcoa.p ,Habitat.p, Habitat.pcoa.p ,cols=2)
adonis_phylogeny
adonis_habitat


multiplot(phylo.disp.p, habitat.disp.p,cols=2)

SIMPER = simper(PRM_data, metadata$group2, permutations=100)




OGs = read.table('~/Github/ProkComp/data_files/panmatrix.txt')

data = t(OGs)


rownames(ANVIO_cat) = ANVIO_cat$name
ANVIO_cat2 = ANVIO_cat[rownames(data),]

PRM_data = data.frame(data)
dPRM_data = distManhattan(PRM_data)


metadata= ANVIO_cat2

adonis_phylogeny = adonis(dPRM_data ~ group2, data = metadata[rownames(PRM_data),])
adonis_habitat = adonis(dPRM_data ~ SS, data = metadata[rownames(PRM_data),])


anosim_phylogeny = anosim(dPRM_data, metadata[rownames(PRM_data),]$group2)
anosim_phylogeny # take a look at results
summary(anosim_phylogeny)
plot(anosim_phylogeny)
anosim_habitat= anosim(dPRM_data, metadata[rownames(PRM_data),]$SS)
anosim_habitat # take a look at results
summary(anosim_habitat)
plot(anosim_habitat)


dispersion<-betadisper(dPRM_data, group= metadata[rownames(PRM_data),]$group2)
permutest(dispersion)
plot(dispersion, hull=F, ellipse=TRUE) ##sd ellipse
PCoA1=dispersion$vectors[,1]
PCoA2=dispersion$vectors[,2]
PCoA.plot<-cbind(PCoA1, PCoA2, metadata[rownames(PRM_data),])
phylo.pcoa.p <-ggplot(PCoA.plot, aes(PCoA1, PCoA2))+
  geom_point(position=position_jitter(.1),aes(color=SS))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill=SS))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+ theme(legend.position = "none")
phylo.pcoa.p

average_dispersion_from_centroid = cbind(data.frame(dispersion$distances),metadata)
habitat.disp.p = ggplot(average_dispersion_from_centroid,aes(x=SS, y=dispersion.distances,fill= SS, color=SS))+geom_boxplot(alpha=0.7)+geom_point()+scale_color_manual(values= Habitat_colors)+scale_fill_manual(values= Habitat_colors)+theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "none")



dispersion<-betadisper(dPRM_data, group= metadata[rownames(PRM_data),]$SS)
permutest(dispersion, pairwise=T, permutation=1000)
PCoA1=dispersion$vectors[,1]
PCoA2=dispersion$vectors[,2]
PCoA.plot<-cbind(PCoA1, PCoA2, metadata[rownames(PRM_data),])

Habitat.pcoa.p <-ggplot(PCoA.plot, aes(PCoA1, PCoA2))+
  geom_point(position=position_jitter(.1),aes(color=group2))+##separates overlapping points
  stat_ellipse(type='t',size =1,alpha=0.2,geom="polygon",aes(fill= group2))+ ##draws 95% confidence interval ellipses
  theme_classic()+scale_color_manual(values= Habitat_colors)+scale_fill_manual(values= Habitat_colors)+ theme(legend.position = "none")
Habitat.pcoa.p

average_dispersion_from_centroid = cbind(data.frame(dispersion$distances),metadata)
phylo.disp.p = ggplot(average_dispersion_from_centroid,aes(x= group2, y=dispersion.distances,fill= group2, color= group2))+geom_boxplot(alpha=0.7)+geom_point()+scale_color_manual(values= colors)+scale_fill_manual(values= colors)+theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "none")

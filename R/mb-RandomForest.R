#-----------------------------------------------#
#	  Random Forest Classifier
#-----------------------------------------------#

#------- Utility Functions! --------#

#plotRF_proxMDS
#===================================================================================#
#' helper function to generate Multi-dimensional Scaling plot of Proximity matrix
#'
#' @param RF
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
plotRF_proxMDS <- function(RF,metadata){
  prox.dist <- dist(1-RF$proximity)
  rownames(metadata) <- metadata$genome
  MDS = cmdscale(prox.dist, eig = TRUE, x.ret=TRUE)
  MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
  MDSdata = data.frame(MDS$points)
  MDSdata$Name = rownames(MDSdata)
  MDSdata$grp = metadata[rownames(MDSdata),'group']
  colnames(MDSdata) = c("x",'y','Name',"group")

  p <- ggplot2::ggplot(MDSdata,aes(x=x,y=y,fill= group,label=Name)) +
    ggplot2::theme_classic() +
    ggplot2::xlab(paste('PCoA1 ', MDS.var.perc[1], '%',sep='')) +
    ggplot2::ylab(paste('PCoA2 ', MDS.var.perc[2], '%',sep='')) +
    ggplot2::labs(fill = "Genus") +
    ggplot2::ggtitle(label='MDS (1 - RF proximity)')+
    ggplot2::geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
    ggplot2::geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
    ggplot2::geom_point(shape = 21, size = 2) +
    fdb_style()
  return(p)
}

#usage
#plotRF_proxMDS(RF,metadata)
#plotRF_proxMDS(RF_perGroup[['cl1']],metadata)
#plotRF_proxMDS(RF_perGroup[['cl2']],metadata)


#plotRF_VarImp
#===================================================================================#
#' helper function to plot variable importance for both MeanDecreaseAcc and MeanDecreaseGini
#'
#' @param RF random forest (RF) object, library(randomForest)
#' @param nvar TOP n of variables contributing to model accuracy (variables are sorted according to importance)
#'
#' @return
#' @export
#'
#' @examples
plotRF_VarImp <- function(RF, nvar=50){
  imp <- data.frame(randomForest::importance(RF))
  imp $features <- rownames( imp )
  imp <- arrange( imp  , desc(MeanDecreaseAccuracy)  )
  g.mdac <-
    ggplot2::ggplot(imp[1:nvar,], aes(x= reorder(features, MeanDecreaseAccuracy), y= MeanDecreaseAccuracy)) +
    ggplot2::geom_segment(color='grey',aes(x= reorder(features, MeanDecreaseAccuracy),xend= reorder(features, MeanDecreaseAccuracy),y=0,yend= MeanDecreaseAccuracy)) +
    ggplot2::geom_point(shape = 21, size = 2,fill='grey') +
    ggplot2::ylab("MeanDecreaseAccuracy") +
    ggplot2::xlab("") +
    ggplot2::coord_flip()+
    ggplot2::theme_classic()+
    ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
      )

  imp <- dplyr::arrange( imp  , desc(MeanDecreaseGini)  )

  g.gini <-
    ggplot2::ggplot(imp[1:nvar,], aes(x=reorder(features, MeanDecreaseGini), y= MeanDecreaseGini)) +
      ggplot2::geom_segment(color='grey',aes(x= reorder(features, MeanDecreaseGini),xend= reorder(features, MeanDecreaseGini),y=0,yend= MeanDecreaseGini)) +
      ggplot2::geom_point(shape = 21, size = 2,fill='grey') +
      ggplot2::ylab("MeanDecreaseGini") +
      ggplot2::xlab("") +
      ggplot2::coord_flip()+
      ggplot2::theme_classic()+
      ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
      )

  p = gridExtra::grid.arrange(g.mdac , g.gini, ncol=2)
  return(p)
}


#usage
#plotRF_VarImp(RF,nvar=50)
#plotRF_VarImp(RF,nvar=10)



#plotRF_VarImpCor
#===================================================================================#
#' helper function to plots the correlation between "MeanDecreaseAcc" and "MeanDecreaseGini"
#'
#' @param RF random forest (RF) object, library(randomForest)
#'
#' @return
#' @export
#'
#' @examples
plotRF_VarImpCor <- function(RF){
  imp <- data.frame(importance(RF))
  imp$features <- rownames( imp )

  p <- ggplot2::ggplot(imp,aes(x= MeanDecreaseAccuracy, MeanDecreaseGini))+
        ggplot2::geom_point(shape = 21, size = 2,fill='grey',alpha=0.5) +
        ggplot2::theme_classic() +
        ggpubr::stat_cor(method = "pearson") +
        ggplot2::theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 11, colour = '#000000'),
            axis.text = element_text(size = 10, colour = '#000000'),
            legend.justification = c(1, 1),
            legend.key.width = unit(0.25, 'cm'),
            legend.key.height = unit(0.55, 'cm'),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 11))
  return(p)
}

#usage
#plotRF_VarImpCor(RF)

#----------------------------------------------------------------------------------
#-- function plotRF_combinedGraph
#   helper function to combine several graph types
#   with:
#     RF: random forest (RF) object, library(randomForest)

plotRF_combined <- function(RF,metadata,nvar=50){
  p.mds = plotRF_proxMDS(RF, metadata=metadata)
  p.varImpCor = plotRF_VarImpCor(RF)
  p.varImp = plotRF_VarImp(RF, nvar)
  p.combined = gridExtra::grid.arrange(
                  gridExtra::arrangeGrob(p.varImp, ncol=1, nrow=1),
                  gridExtra::arrangeGrob(p.varImpCor, p.mds , ncol=2, nrow=1),nrow=2, heights=c(2,1))
  return(p.combined)
}




#==================================================================================
# Main computation
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

Res_pan <- NULL
for(pan.name in c('OGpan')){
#for(pan.name in names(pan.list)){
  message('you are brave... better have time to wait!')
  message(paste0('=====  ',pan.name,'  ====='))

  pan = pan.list[[pan.name]]
  outList_hgm.all = list()
  Res_col <- NULL
  for(sel.col in columns){
    message(paste0('-- ',sel.col))
    metadata <-genome.tbl %>%
      dplyr::left_join(sorted.cliques, by='genome') %>%
      dplyr::left_join(out.colors, by='genome') %>%
      dplyr::filter(!is.na(!!as.name(sel.col))) %>%
      data.frame()
    metadata$grp <- metadata[,sel.col]
    metadata <- metadata %>%
      dplyr::mutate(group=paste0('cl',grp)) %>%
      dplyr::select(-grp)

    if(remove_small_groups == TRUE){
      larger_groups <- metadata %>% group_by(group) %>% tally() %>% filter(n>2) %>% select(group) %>% pull()
      metadata<-metadata %>% filter(group %in% larger_groups)

    }

    pan = pan[metadata$genome,]
    pan = pan[,which(colSums(pan)!=0)]
    pan = pan[,which(colSums(pan)!=1)]
    pan <-pan[,which(apply(pan,2,sd)>0)]


    Res <- NULL
    RF_perGroup = list()

    for(grp in as.character(unique(metadata$group))){
      #grp='cl1'
      print(paste('RF for ',grp,sep=''))
      RFdf <- data.frame(pan,'group'=as.character(metadata$group))

      #binarises the whether genome is part of selected group or not, (2-class, 1 for yes, 0 for no )
      RFdf$group <- as.factor(RFdf$group == grp)
      RF <- randomForest::randomForest(group ~ ., data= RFdf,importance=TRUE ,proximity=TRUE, ntree=5000, strata= group)
      res <-
        c(
          group = grp,
          OOB.error = unname(RF$err.rate[nrow(RF$err.rate),1]),
          ingroup.class.error = RF$confusion[2,'class.error'],
          ingroup.size = sum(RF$confusion[2,c('TRUE','FALSE')]),
          ingroup.success = RF$confusion[2,c('TRUE')],
          ingroup.fail = RF$confusion[2,c('FALSE')],
          outgroup.class.error = RF$confusion[1,'class.error'],
          outgroup.size = sum(RF$confusion[1,c('TRUE','FALSE')]),
          outgroup.success = RF$confusion[1,c('FALSE')],
          outgroup.fail = RF$confusion[1,c('TRUE')]
        )
      Res = rbind(Res, res)
      RF_perGroup[[grp]] = RF
    }
    Res <- data.frame(Res)
    Res$ingroup.size = as.numeric(as.character(Res$ingroup.size))
    Res$ingroup.class.error =  as.numeric(as.character(Res$ingroup.class.error))

    save(RF_perGroup, file=paste0("~/DATA/MarbGenomics/",'multiRF',pan.name,sel.col,remove_small_groups,'.RData'))
    Res_col <- rbind(Res_col, cbind(Res,'level'=sel.col))
    }
  Res_pan <- rbind(Res_pan, cbind(Res_col,'pan'=pan.name))
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

head(Res_pan) %>% data.frame() %>% str()

Res_pan$OOB.error <- as.numeric(as.character(Res_pan$OOB.error))

p.RF <- Res_pan %>% ggplot(aes(level,OOB.error))+
  geom_boxplot(fill='grey',alpha=0.3,outlier.shape = NA)+
  geom_jitter(shape=21,alpha=.7,fill='grey')+
  facet_wrap(~pan,scales='free',nrow=1)+fdb_style() + theme(axis.text.x = element_text(angle = 90,hjust = 1))

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/RF_OOB_error_allgroups_pfam.pdf"),plot=p.RF, width = 13, height = 3,unit='in')

Res_pan <- Res_pan$cutoff <- as.numeric(gsub('sc','',Res_pan$level))
p.RF <- Res_pan %>% ggplot(aes(cutoff,OOB.error,color=pan,fill=pan))+
  geom_smooth(method='lm')+
  geom_point(shape=21,alpha=.7)+facet_wrap(~pan,scales='free',nrow=1)+
  fdb_style() + theme(axis.text.x = element_text(angle = 90,hjust = 1))

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/RF_OOB_error_allgroups_cont_pfam.pdf"),plot=p.RF, width = 13, height = 3,unit='in')




#---------------
# Using stored data
#=====================================================================#
# unpack
pan.name= 'OGpan'
pan.name= 'pfam'

sel.col = "sc0.85"

#Load in the data
load(paste0("~/DATA/MarbGenomics/",'multiRF',pan.name,sel.col,'.RData'))

RF_perGroup
importance.all <- NULL
importance.wide <- rownames(RF_perGroup[[1]]$importance)  %>% data.frame()
for(grp in names(RF_perGroup)){
  set <- RF_perGroup[[grp]]$importance %>% data.frame() %>%
    tibble::rownames_to_column(var='feature.id') %>%
    dplyr::mutate(group=grp)
  importance.all<-rbind(importance.all, set)
  s <- set[,c('MeanDecreaseAccuracy','MeanDecreaseGini')]
  colnames(s) <- paste(grp,c('MDA','MDG'),sep='_')
  importance.wide <- cbind(importance.wide, s)
}


#extract 10 most important variables per group!

importance.highest <- importance.all %>%
  tibble::as_tibble() %>%
  dplyr::group_by(group) %>%
  dplyr::slice_max(order_by = MeanDecreaseAccuracy, n = 20)

p.top10 <- importance.highest %>%
  ggplot2::ggplot(aes(tidytext::reorder_within(feature.id, MeanDecreaseAccuracy, group),MeanDecreaseAccuracy))+
    ggplot2::geom_bar(stat='identity')+
    ggplot2::coord_flip()+
    ggplot2::facet_wrap(~group,scales='free')+
    ggplot2::ylab('feature')+
    fdb_style()

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/RF_top10_",pan.name,"_",sel.col,".pdf"),plot=p.top10, width = 8, height = 16,unit='in')

importance.highest <- importance.all %>%
  tibble::as_tibble() %>%
  dplyr::group_by(group) %>%
  dplyr::slice_max(order_by = MeanDecreaseAccuracy, n = 100)

p.top100 <- importance.highest %>%
  ggplot2::ggplot(aes(tidytext::reorder_within(feature.id, MeanDecreaseAccuracy, group),MeanDecreaseAccuracy))+
  ggplot2::geom_bar(stat='identity')+
  ggplot2::coord_flip()+
  ggplot2::facet_wrap(~group,scales='free')+
  ggplot2::ylab('feature')+
  fdb_style()+
  ggplot2::theme(axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/RF_top100_",pan.name,"_",sel.col,".pdf"),plot=p.top100, width = 8, height = 8,unit='in')




#comparing hypergeometrix test to RF
#====================================================================#
importance.highest <- importance.all %>%
  tibble::as_tibble() %>%
  dplyr::group_by(group) %>%
  dplyr::slice_max(order_by = MeanDecreaseAccuracy, n = 100)

RF.top100.all <- importance.highest %>% select(feature.id) %>%
  pull() %>%
  as.character()

Hypg.enriched <- ANI.hmg.grp.tbl %>%
  filter(level==sel.col) %>%
  filter(fdr.p.value_enriched_binary < 0.01) %>%
  select(feature.id) %>%
  pull() %>%
  as.character()

Hypg.depleted <- ANI.hmg.grp.tbl %>%
  filter(level==sel.col) %>%
  filter(fdr.p.value_depletion_binary < 0.01) %>%
  select(feature.id) %>%
  pull() %>%
  as.character()

Hypg.all <- c(Hypg.depleted,Hypg.enriched) %>% unique()

gplots::venn(list('RF'=RF.top100.all, 'hypg'=Hypg.all))





out.hypg <-  ANI.hmg.grp.tbl %>%
  filter(level==sel.col) %>%
  #filter(fdr.p.value_enriched_binary < 0.01) %>%
  filter(group=='cl8')

out.RF <- importance.all %>%
  tibble::as_tibble() %>% filter(group=='cl8')

out.combined.group <- out.RF %>% left_join(out.hypg,by='feature.id')

out.combined.group %>% ggplot(aes(MeanDecreaseAccuracy,fdr.p.value_enriched_rawcounts))+geom_point()


#----------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------
#
#     2-class RF models (per group),
#       (in group =TRUE, outgroup = FALSE)
#
#---------------------------------------------------------------------------------------------

# Test whether we get better predictions of what is is to belong to cl17 for examlpe whem using a binary class based thing
# Given that our dataset is quite wierdly distributed, (some clades well represented and others not), we try different things conceptually
#   2-class (binary) and multiclass scenarios

#For loop over all the groups, and do a 2-class approach (in group =TRUE, outgroup = FALSE)
#BinaryRF <- function(pan, metadata){
#
#
#}

pan = cazy.pan
#pan = MCRpan




#======== SINGLE

pan = pan[metadata$genome,]
pan = pan[,which(colSums(pan)!=0)]
pan = pan[,which(colSums(pan)!=1)]
pan <-pan[,which(apply(pan,2,sd)>0)]


Res = NULL
RF_perGroup = list()

for(grp in as.character(unique(metadata$group))){
  #grp='cl1'
  print(paste('RF for ',grp,sep=''))
  RFdf <- data.frame(pan,'group'=as.character(metadata$group))

  #binarises the whether genome is part of selected group or not, (2-class, 1 for yes, 0 for no )
  RFdf$group <- as.factor(RFdf$group == grp)
  RF <- randomForest(group ~ ., data= RFdf,importance=TRUE ,proximity=TRUE, ntree=5000, strata= group)
  res <-
    c(
      group = grp,
      OOB.error = unname(RF$err.rate[nrow(RF$err.rate),1]),
      ingroup.class.error = RF$confusion[2,'class.error'],
      ingroup.size = sum(RF$confusion[2,c('TRUE','FALSE')]),
      ingroup.success = RF$confusion[2,c('TRUE')],
      ingroup.fail = RF$confusion[2,c('FALSE')],
      outgroup.class.error = RF$confusion[1,'class.error'],
      outgroup.size = sum(RF$confusion[1,c('TRUE','FALSE')]),
      outgroup.success = RF$confusion[1,c('FALSE')],
      outgroup.fail = RF$confusion[1,c('TRUE')]
    )
  Res = rbind(Res, res)
  RF_perGroup[[grp]] = RF
}
Res <- data.frame(Res)
Res$ingroup.size = as.numeric(as.character(Res$ingroup.size))
Res$ingroup.class.error =  as.numeric(as.character(Res$ingroup.class.error))

save(RF_perGroup, file=paste0("~/DATA/MarbGenomics/",'multiRF','.RData'))










importance.all <- NULL
importance.wide <- rownames(RF_perGroup[[1]]$importance)  %>% data.frame()
for(grp in names(RF_perGroup)){
  set <- RF_perGroup[[grp]]$importance %>% data.frame()
  importance.all<-rbind(importance.all, set %>% mutate(group=grp))
  s <- set[,c('MeanDecreaseAccuracy','MeanDecreaseGini')]
  colnames(s) <- paste(grp,c('MDA','MDG'),sep='_')
  importance.wide <- cbind(importance.wide, s)
}



ggplot(Res, aes(ingroup.size,ingroup.class.error,label=group))+
  geom_point()+geom_label()+fdb_style(aspect.ratio=1)

#look at those groups that have some sort of success rate (say below 25% error rate)
Res %>% filter(ingroup.class.error<0.25)

#well represented groups, say above a genome nr of 5 genomes (approx, 5% of dataset)
wellrep.grps = Res %>% filter(ingroup.size>=5) %>% select(group) %>% pull() %>% as.character() %>% unname()


#------VarImp per group

#Show all RF proximity based MDS plots at once
proxMDS_plotlist = list()
for(grp in names(RF_perGroup)){
  proxMDS_plotlist[[grp]]=plotRF_proxMDS(RF_perGroup[[grp]],metadata)
}
gridExtra::grid.arrange(grobs = proxMDS_plotlist, ncol=4)

#Combined graphics per group

RFcombined_plotlist = list()
for(grp in names(RF_perGroup)){
  RFcombined_plotlist[[grp]] = plotRF_VarImpCor(RF = RF_perGroup[[grp]],metadata=metadata,nvar=25)
}
#Show all RF var importance of key groups
RF_VIP_plotlist = list()
for(grp in c('cl1','cl3','cl4','cl7','cl14')){
  RF_VIP_plotlist[[grp]]=plotRF_VarImp(RF_perGroup[[grp]], nvar=25)
}
gridExtra::grid.arrange(grobs = RF_VIP_plotlist, ncol=3)



#---------------------------------------------------------------------------------------------
# USE SOME OF THE VISUALISATION FUNCTIONS

plotRF_proxMDS(RF = RF_perGroup[['cl14']],metadata=metadata)
plotRF_VarImp(RF= RF_perGroup[['cl14']], nvar=25)
plotRF_VarImpCor(RF = RF_perGroup[['cl14']])
plotRF_combined(RF = RF_perGroup[['cl14']],metadata=metadata, nvar=25)

RF_perGroup[['cl14']]


#----! NEED TO IMPLEMENT COLORATION OF plotRF_VarImp function, enabling the notion that in the in vs outgroup the feature is more or less abundance!




#---------------------------------------------------------------------------------------------

m2n[m2n$module %in% c('M00618', 'M00579', 'M00628', 'M00627', 'M00668'),]
m2n[m2n$module %in% c('M00497', 'M00016', 'M00639', 'M00207', 'M00133', 'M00507', 'M00001'),]



#---- using the multiclass scenatio, I always find bad predictability of certain poorly represented clades


set.seed(1234)
RFdf = data.frame(pan,'group'=as.character(metadata$group))
RFdf %>% count(group)
RF <- randomForest(group ~ ., data= RFdf,importance=TRUE ,proximity=TRUE, ntree=5000, strata= group)

print(RF)
#	out of bag error (OOB.error) is 17.02%
#	OOB.error = 17.02
#	accuracy = 100 - 17.02 or 82.98%....
# "out of bag estimate of error rate" [basically uses the data not seen by the tree to use for calculating the error of that sp. tree (OOB data)]
#look at confusion matrix, and class.error,
#it seems that predictions are very good for many of the classes/groups (approx. 0 class error)
# while for other groups, it is pretty bad, approximating 1



#---------------
# RF optimiser

#optimise mtry

mtry_tests = c(1:20)
val.ntree = 3000
val.ntree = 100

oob.values <- c()
for(i in mtry_tests){
  temp.RF <- randomForest(group ~ ., data= RFdf, mtry = i, ntree= val.ntree, strata= group)
  oob.values = c(oob.values, temp.RF$err.rate[nrow(temp.RF$err.rate),1])
  print(oob.values)
}


df.mtry.res = data.frame(cbind('mtry' = mtry_tests, 'OOB.error' = oob.values))
ggplot(df.mtry.res, aes(mtry, OOB.error)) + geom_line()	+geom_point()

#---------------



attributes(RF)

#look for plateau error rates
plot(RF)
# looks like arround 2000 trees, the model does not improve much...

# to give an idea on what mtry value to use: we can use tuneRF
# you can exrtract it from the t.RF object
t.RF = tuneRF(
  x = RFdf[,!colnames(RFdf)=='group'] ,
  y = RFdf$group,
  stepFactor = 0.05,
  plot=TRUE,
  ntreeTry = 4000,
  trace = TRUE,
  improve = 0.01)


RF2 <- randomForest(
  group ~ .,
  data= RFdf,
  importance=TRUE ,
  proximity=TRUE,
  ntree=4000,
  mtry=41,
  strata= group)

print(RF2)
print(RF)

hist(treesize(RF),
     main='No. of Nodes for the Trees',
     col='grey')

#Variiable importance
# MeanDecreaseAccuracy: 	for example, if we remove COH4307 when making the trees, what is the mean decrease in accuracy. most important variables contributing to the accuracy
#MeanDecreaseGini: "how pure the nodes are at the end of the tree without each variable"
varImpPlot(RF, sort=TRUE, n.var=50, main='top 50 - Variable Importance')


importance(RF) #if appr. 0, then they are "not of much use"
varUsed(RF) # how many times variables were used in the model (should correlate to the var importance estimates)

#-- Partial Dependence Plot
#Partial Dependence plot gives a graphical depiction of the marginal effect of a variable on the class probability (classificaiton) or response (regression)
partialPlot(RF, RFdf, c('COG4307'), 'cl1') #say, if COG4307 is 2, it tend to predict it is cl1
partialPlot(RF, RFdf, 'COG4307', 'cl7') #say, if COG4307 is 2, it tend to predict it is cl1





#-----

oob.error.data <- data.frame(
  Trees=rep(1:nrow(RF$err.rate),times=4),
  Type=rep(c('OOB','cl1','cl7','cl14'),each=nrow(RF$err.rate)),
  Error=c(RF$err.rate[,"OOB"],
          RF$err.rate[,"cl1"],
          RF$err.rate[,"cl7"],
          RF$err.rate[,"cl14"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error))+
  geom_line(aes(color=Type))+theme_classic()+theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
  )


#----------------------------------



#----------------------------------




RF_classify_imp <- as.data.frame( RF $importance )
RF_classify_imp$features <- rownames( RF_classify_imp )
RF_classify_imp_sorted <- arrange( RF_classify_imp  , desc(MeanDecreaseAccuracy)  )

barplot(RF_classify_imp_sorted $MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

barplot(RF_classify_imp_sorted[1:30,"MeanDecreaseAccuracy"], names.arg= RF_classify_imp_sorted[1:30,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=2, main="Classification RF")

#EggNog[RF_classify_imp_sorted[1:30,"features"],"eggNOG_HMM"]


heatmap.2(t(pan[RF_classify_imp_sorted[1:100,"features"],]),trace="none")




scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)


heatmap.2(pan[,RF_classify_imp_sorted[1:500,"features"]], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)







#----- Prediction & confusion matrix - train data
library(caret)

p1 <- predict(RF, RFdf)
caret::confusionMatrix(p1,RFdf$group)

#explore the accuracy of the model, and look at 95%CI


#----



v1 = venn(list('hgm'= allSig_hgm,'RF500'= RF_classify_imp_sorted[1:500,"features"]))
v2 = venn(list('hgm.raw'= all.raw.grp_hgm,'hgm.binary'= all.bin.grp_hgm, 'RF500'= RF_classify_imp_sorted[1:500,"features"]))

intersection_v1 = attr(v1,"intersections")$`hgm:RF500`
intersection_v2 = attr(v2,"intersections")$`hgm.raw:hgm.binary:RF500`


heatmap.2(pan[, intersection_v1], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)

#----



distMan <- distManhattan(pan[, intersection_v1])
distMan = as.matrix(distMan)
distMan = distMan[ order(row.names(distMan)), ]
distMan = distMan[ , order(colnames(distMan))]
#subDistMan = distMan[grepl('Marinobacter|Tamil',rownames(distMan)),grepl('Marinobacter|Tamil',colnames(distMan))]


MDS = cmdscale(distMan, eig = TRUE)
MDSdata = data.frame(MDS$points)
MDSdata$Name = rownames(MDSdata)
MDSdata$grp = metadata[rownames(MDSdata),'group']
colnames(MDSdata) = c("x",'y','Name',"group")


ggplot(MDSdata,aes(x=x,y=y,fill= group,label=Name)) +
  theme_classic() +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  labs(fill = "Genus") +
  ggtitle(label='Manhattan distance')+
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




#-------------------------------------------
#
#
#
#
#
#


#---------------------------------------
#
#   Multi-Class RF models using bootstrapping (tidymodels)
#
#---------------------------------------

RFdf %>%
  group_by(group) %>%
  summarise(across(c(1:20), list(mean = mean, sd = sd)))

# one can see already that things are very uneven
RFdf %>% count(group)

#boostrap sample the dataset
RFfolds = RFdf %>%
  mutate(group = as.factor(group))  %>%
  bootstraps(times=25)


#Set up random forest model
rf_spec <-  rand_forest(trees=1000) %>%
  set_mode('classification') %>%
  set_engine('ranger')

#worflow to use all bootsrtap sample sets
rf_wf <- workflow() %>%
  add_formula(group ~ COG0019 + COG0028 ) %>% add_model(rf_spec)

rf_wf

#run random forest to the bootstrap resamples
doParallel::registerDoParallel()
RF_res <- fit_resamples(
  rf_wf,
  resamples = RFfolds,
  control = control_resamples(save_pred = TRUE)
)

RF_res

# -- evalutation

collect_metrics(RF_res)

# --

RF_res %>%
  collect_predictions() %>%
  group_by(id) %>%
  ppv(group, .pred_class)

RF_res %>%
  collect_predictions() %>%
  group_by(id) %>%
  roc_curve(group, .pred_cl1:.pred_cl9)  %>%
  ggplot(aes(1-specificity, sensitivity, color = id)) +
  geom_abline(lty=2,color='grey80',alpha=-.6) +
  geom_path(show.legend = FALSE, alpha=.6) +
  facet_wrap(~.level, ncol=5)





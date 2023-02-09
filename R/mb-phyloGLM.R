library(phylolm)
#---------------------------------------
#
#   phyloglm
#
#---------------------------------------
# -- since I am using R version 3.6.3 I had to install from source, and download an older version from CRAN
#  install.packages('~/Downloads/phylolm_2.5.tar.gz', repos = NULL, type="source")
#	!!!!!!!!!!!!!
#	frankly, this does not makes much sense when asking for lineage specific enrichment, does it....
#	nonetheless, we can use this to look for habitat specific enrichment!

#fit = phyloglm(gene~group,phy=tree,data= dat,boot=100)


#pan_phyloglm
#=================================================================================#
#' Title
#'
#' @param pan
#' @param tree
#' @param metadata
#' @param grp
#'
#' @return
#' @export
#'
#' @examples
#' pan_phyloglm(pan=cazy.pan, tree=tree, metadata=metadata,grp='cl8')
pan <- m.pan[,1:5]
out.phm <- pan_phyloglm(pan=cazy.pan, tree=tree, metadata=metadata,grp='cl8')

pan_phyloglm <- function(pan, tree, metadata, grp){
  pan = pan[metadata$genome,]
  pan = pan_clean(pan, remove_singletons = TRUE, SD=2)

  tree = drop.tip(tree, tip = setdiff(tree$tip.label,metadata$genome))
  pan = pan[tree$tip.label,]

  dat = metadata[,c('genome','group')]
  dat$grp = as.numeric(dat$group == grp)
  rownames(dat) = dat$genome

  output <- NULL
  for(feature in colnames(pan)){
    dat$gene = pan[metadata$genome, feature]
    dat$gene[dat$gene > 0] = 1
    m1 <- tryCatch(
      phylolm::phyloglm(gene~ grp,phy=tree,data= dat, method = "logistic_IG10")
      , error = function(e) list(coefficients = NA))

    if(is.na(coef(m1)[1])){
      res <- data.frame(
        feature.id = feature,
        Estimate = NA,
        SE = NA,
        z.value = NA,
        p.value = NA,
        group=grp)

    }else{
      m1.sum <- summary(m1)
      res <- data.frame(
        feature.id = feature,
        Estimate = m1.sum$coefficients["grp",1],
        SE = m1.sum$coefficients["grp",2],
        z.value = m1.sum$coefficients["grp",3],
        p.value = m1.sum$coefficients["grp",4],
        group=grp)

      if(m1.sum$coefficients["grp",4]<0.05){
        message(res)
        }
      }
    output=rbind(output, res)
  }
  output$p.fdr <- p.adjust(output$p.value, method='fdr')
  output<-data.frame(output)
  return(output)
}

# usage
# out.phm <- pan_phyloglm(pan=cazy.pan, tree=tree, metadata=metadata,grp='cl8')
#pan_phyloglm.all
#=================================================================================#
#' Title
#'
#' @param pan
#' @param tree
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
#' pan_phyloglm.all(pan=cazy.pan, tree=tree, metadata=metadata)
pan_phyloglm.all <- function(pan, tree, metadata){
  outList_phm = list()
  for(grp in unique(as.character(metadata$group))){
    message(paste('analysing', grp))
    outList_phm[[grp]] = pan_phyloglm(pan[metadata$genome,], tree, metadata, grp)
  }
  phm.grp.tbl <- outList_phm  %>%
    tibble::as_tibble_col() %>%
    tidyr::unnest(cols=c(value))

  return(phm.grp.tbl)
}






#=============
# EXAMPLE

metadata <-genome.tbl %>%
  dplyr::left_join(sorted.cliques, by='genome') %>%
  dplyr::left_join(out.colors, by='genome') %>%
  dplyr::filter(!is.na(!!as.name(selcol)))%>%
  mutate(group=paste0('cl',c0.85)) %>%
  data.frame()
rownames(metadata) <- metadata$genome

#metadata <- out.out.tbl %>%
#  filter(type=='ANIb') %>%
#  filter(variable == 'sc0.8') %>%
#  filter(method=='cluster_edge_betweenness') %>%
#  mutate(group=paste0('cl',value))

sel.tree <- lsTrees.94.rooted[[1]]
metadata_ <- genome.tbl %>% filter(genome %in% sel.tree$tip.label) %>%
  dplyr::left_join(sorted.cliques, by='genome') %>%
  dplyr::left_join(out.colors, by='genome') %>%
  #dplyr::filter(!is.na(!!as.name(selcol)))%>%
  #mutate(group=paste0('cl',c0.85)) %>%
  mutate(group=phylogroup.x) %>%
  data.frame()
rownames(metadata) <- metadata$genome

phglm.out <- pan_phyloglm.all(pan=m.pan[sel.tree$tip.label,], tree=sel.tree, metadata=metadata_)

#h.pan <- t(COGCATpan[,sel.genome])
#h.pan <- h.pan[,2:ncol(h.pan)]
#phglm.out <- pan_phyloglm.all(pan=h.pan, tree=tr.tree, metadata=metadata)
#phglm.out <- pan_phyloglm.all(pan=OG.pan[sel.genome,], tree=tr.tree, metadata=metadata)

write.table(phglm.out,file='~/DATA/MarbGenomics/phglm_MCR.tsv',sep='\t',quote = FALSE)

pan_phyloglm(m.pan[sel.tree$tip.label,], sel.tree, metadata_, "P1")


phglm.out %>% head()
phglm.out %>% filter(p.value<0.05) %>%
  #utate(grp = c(rep('cl1',13),rep('clNA',7),rep('cl2',10),rep('cl3',1),rep('cl4',2),rep('cl5',23),rep('cl6',1),rep('cl7',2),rep('cl8',10),rep('cl9',3))) %>%
  left_join(m2n,by=c('feature.id'='module')) %>%
  #filter(grp == 'cl2') %>%
  ggplot(aes(Estimate, DESCRIPTION)) + geom_bar(stat='identity',aes(fill=p.value)) +
  facet_wrap(~group,ncol=2,scales='free')





#=================================================================================#
# the ultra loop
# be brave, or have time!
#=================================================================================#

#=================================================================================#
# The below runs over all feature abundance tables (OG, KO, etc) as well as over all ANI clique levels
# This takes a while, but generates output files that are stored
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!! WARNING LONG COMPUTATION TIME
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
outdir = "~/DATA/MarbGenomics/"
pan.list = panList
tree <- drop.tip(concat.rps.nuc, tip=outgroup)
columns <- c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')
sorted.cliques = sorted.ANI.cliques
out.colors = out.colors

#!! SET TREE, AS I MAY HAVE MESSED HIM UP
#tree = tr.tree

for(pan.name in names(pan.list)){
  message('you are brave... better have time to wait!')
  message(paste0('=====  ',pan.name,'  ====='))

  pan = pan.list[[pan.name]]
  outList_phm.all = list()
  for(selcol in columns){
    message(paste0('-- ',selcol))
    metadata <-genome.tbl %>%
      dplyr::left_join(sorted.cliques, by='genome') %>%
      dplyr::left_join(out.colors, by='genome') %>%
      dplyr::filter(!is.na(!!as.name(selcol))) %>%
      data.frame()
    metadata$grp <- metadata[,selcol]
    metadata <- metadata %>%
      dplyr::mutate(group=paste0('cl',grp)) %>%
      dplyr::select(-grp)

    #use panHypergeometric.all to run all groups versus non-groups
    phm.grp.tbl <- pan_phyloglm.all(pan, tree, metadata)
    outList_phm.all[[selcol]] = phm.grp.tbl %>% mutate(level=selcol)
  }

  out.tbl <- outList_phm.all  %>%
    tibble::as_tibble_col() %>%
    tidyr::unnest(cols=c(value))

  write.table(out.tbl , file = paste0(outdir,pan.name,'phm_allgroups.tsv'),sep='\t')
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#






combRaw = as.character(output[output$p.value<0.05 & !is.na(output$p.value),]$feature.id)

scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)
heatmap.2(pan[, combRaw], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)

cog_desc[combRaw,]



#=================================================================================#
# ANCESTRAL STATE RECONSTRUCTION OF PAN GENOME
# PREDICTION OF GENE LOSS AND GENE AQUISITION
#=================================================================================#

# Ancestral state reconstruction

tr=tr.tree
pan=cazy.pan
outdir = '~/DATA/MarbGenomics/Graphs/'
pan.name='KO'
metadata <- genome.tbl %>% filter(genome %in% sel.genome)%>% data.frame()

pan.list = list('OG'=OG.pan,
                'KO'=kfm.pan,
                'COG'=COG.pan,
                'CAZY'=cazy.pan,
                'TCDB'=tcdb.pan)


for(pan.name in names(pan.list)){

}



#=====
# Load the output back in and detect which of the features are associated with gene loss or gain events
# assess how these are strunctured in the genomes by proximity cluster detection
# detect the flagellar operon for example, that'd be cool to map here
# we can perhaps pull out som very interesting things,
pan.name='OG'
test <-read.table(paste0('~/DATA/MarbGenomics/',pan.name,'_ancestral_state_per_node','.tsv'))
test2 <-read.table(paste0('~/DATA/MarbGenomics/',pan.name,'_ancestral_state_per_feature','.tsv'))




#==========================================
#next assess character state distribution of isolation sources
metadata %>% select(sourceClass)

#could implement Richardson et al 2018 approach here
grp <- 'SS'
mtdt <- metadata %>% mutate(SS = as.character(SS ))
mtdt[mtdt$SS == '','SS'] <- 'other'

mtdt <- mtdt %>% select(genome,SS ,lon,lat,genus)

x <- mtdt[,'SS']
names(x) <- mtdt[,'genome']


fitER<-ace(x,tr.tree,model="ER",type="discrete")
fitER

plotTree(tr,fsize=0.8,ftype="i")
nodelabels(node=1:tr$Nnode+Ntip(tr),
           pie=fitER$lik.anc,cex=0.5)
tiplabels(pie=to.matrix(x,sort(unique(x))),cex=0.3)
add.simmap.legend(prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tr)),fsize=0.8)



#=== GGTREE APPROACH
habitat_inf <- ace(x, tr, type="discrete", method="ML")
i <- apply(fitER$lik.anc, 1, which.max)
habitat.df <- data.frame(node=1:Nnode(tr, internal.only=F),
                      habitat=c(x, sort(unique(x))[i]))

dp <- ggtree(tr)
dp$data <- dp$data %>% left_join(habitat.df, by='node')
#tr <- habitat.df

dp + aes(color=habitat) + geom_point(aes(color=habitat)) + scale_color_manual(values=cols)#+
  #theme(legend.position=c(.25, .85))


#An alternative approach to the one outline above is to use an MCMC approach to sample character histories from their posterior probability distribution. This is called stochastic character mapping (Huelsenbeck et al. 2003). The model is the same but in this case we get a sample of unambiguous histories for our discrete character's evolution on the tree - rather than a probability distribution for the character at nodes.
#For instance, given the data simulated above - we can generate the stochastic character map as follows:
cols = rev(colors()[seq(2,152,10)])
names(cols) <- c(unique(mtdt$SS))


mtree<-make.simmap(tr,x,model="ER")
plot(mtree,cols,fsize=0.8,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=-max(nodeHeights(tr)),fsize=0.8)


mtrees<-make.simmap(tr,x,model="ER",nsim=100)
par(mfrow=c(10,10))
null<-sapply(mtrees,plot,colors=cols,lwd=1,ftype="off")
pd<-summary(mtrees,plot=FALSE)
plot(pd,fsize=0.6,ftype="i")
plot(pd,fsize=0.6,colors=cols,ftype="i")

## now let's plot a random map, and overlay the posterior probabilities
plot(mtrees[[1]],cols,fsize=0.8,ftype="i")
nodelabels(pie=pd$ace,piecol=cols,cex=0.2)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(tr)),fsize=0.8)

#Finally, since we obtained these inferences under exactly the same model, let's compare the posterior probabilities from stochastic mapping with our marginal ancestral states. In the former case, our probabilities were obtained by sampling from the joint (rather than marginal) probability distribution for the ancestral states.
plot(fitER$lik.anc,pd$ace,xlab="marginal ancestral states",
     ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)

#This tells us that although joint & marginal reconstruction are not the same, the marginal probabilities from joint stochastic sampling and the marginal ancestral states are quite highly correlated - at least in this case study.


pd$ace


i <- apply(pd$ace, 1, which.max)
habitat.df <- data.frame(node=1:Nnode(1:Nnode(tr, internal.only=F), internal.only=F),
                         habitat=c(x, unique(x)[i]))

dp <- ggtree(tr)
dp$data <- dp$data %>% left_join(habitat.df, by='node')
#tr <- habitat.df

dp + aes(color=habitat) + geom_point(aes(color=habitat)) + scale_color_manual(values=cols)#+





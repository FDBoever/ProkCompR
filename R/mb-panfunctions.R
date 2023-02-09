# mb-panfunctions
#
#-----------------------------------------------------------------------------#
#   UTILITY FUNCTIONS   ####
#-----------------------------------------------------------------------------#

#remove_rare
#=================================================================================#
#' removing rare features
#' @param pan abundance/pan table.
#' @param fraction if gene is present below treshold
#'
#' @return
#' @export
#'
#' @examples
#' pan.rare_removed <- remove_rare(pan= cazy.pan, fraction=0.02)

remove_rare <- function( pan ,  fraction=0.05 ) {
  pan <- t(pan)
  row2keep <- c()
  cutoff <- ceiling( fraction * ncol(pan) )
  for ( i in 1:nrow(pan) ) {
    row_nonzero <- length( which( pan[ i , ]  > 0 ) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  output <- t(pan [ row2keep , , drop=F ])
  return(output)
}



# mostAbundantParalogs
#=================================================================================#
#' function to select the most abundant paralogs/variables in the pan-table
#'
#' @param pan specifies pan table/abundance table (genomes x features)
#' @param top_n specifies the number of top hits, most abundant
#'
#' @return
#' @export
#'
#' @examples
#' selected_var = mostAbundantParalogs(pan = COG.pan, top_n = 10)

mostAbundantParalogs <- function(pan, top_n){
  pan = t(pan)
  panRankedParalogs = data.frame(sort(rowSums(pan),decreasing=TRUE))
  colnames(panRankedParalogs) = c('count')
  selected_var = head(rownames(panRankedParalogs), top_n)
  return(selected_var)
}


# pan2ggr
#=================================================================================#
#' Calculates Gene Repertoire relatedness index (GRR)
#' Between two genomes defined as the number of common gene families (the intersection) divided by the number of genes in the smallest genome.
#' It is close to 100% if the gene repertoires are very similar (or one is a subset of the orhter) and lower otherwise
#'
#' @param pan
#'
#' @return
#' @export
#'
#' @examples
#' pan.ggr = pan2ggr(OG.pan)
#' head(pan.ggr)

pan2ggr <- function(pan) {
  g1tr = c()
  pan <- t(pan)
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




# COULD BE DEVELOPED
# Specific to group
#=================================================================================#
# groupSpecific_pan <- function(pan, ingroup){
#   if(!all(ingroup %in% rownames(pan))){
#     stop('ingroup should be subset of rownames(pan)')
#   }else{
#       message('lucky you')
#     }
#   }
#
#
# groupSpecific_pan(pan = OG.pan, ingroup=c('Marinobacter_sp_FDB33','Marinobacter_algicola_DG893'))
#
# selectedgenomes = metadata %>% filter(group=='cl14') %>% select(genome) %>% pull
# groupSpecific_pan(pan = OG.pan, ingroup=selectedgenomes)


###################
# 	pangenome curves
###################

genome.tbl %>% select(genus) %>% pull() %>% unique()

combSP = c()
for (i in genera_selected){
  print(i)
  s_tOF = tOF[grep(paste(i,"_",sep=""),rownames(tOF)),]
  s_tOF = data.frame(s_tOF[,colSums(s_tOF) > 0])
  if (nrow(s_tOF)>1) {
    sp <- specaccum(s_tOF, 'random', permutations=100)
    data_sp <- data.frame(Genomes= sp $sites, Orthologs= sp$richness, SD= sp$sd)
    combSP = rbind(combSP, cbind(data_sp,"Genus"=rep(i,nrow(data_sp))))
    p = ggplot() +
      geom_point(data= data_sp, aes(x= Genomes, y= Orthologs)) +
      geom_line(data= data_sp, aes(x= Genomes, y= Orthologs)) +
      geom_ribbon(data= data_sp ,aes(x= Genomes, ymin=(Orthologs-2*SD),ymax=(Orthologs +2*SD)),alpha=0.2) + 			scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous( expand = c(0, 0)) +
      theme_classic()

    print(p)
  }
}



#=================================================================================#

# SOME RANDOM BRAIN STORM ON MORE FUNCTIONS TO INCLUDE

#   Unique per group
#   Uniquely absent per group
#   pan group summary
#   pan overall summary
#   pan to softpan # look for definition of soft pan and implement
#   paralog expansion #Look for something smart here
#   pan2node


# total number of stuff on a ratio per predicted_genes
# total pan genome identified features

# a function that generates weighted pangenomes?
# weight pangenome on nr of predicted genes per genome?!

#=================================================================================#



# ------------------#



#Rarefraction an power law curves

#Analyses Methods
#' @description
#' Rarefact pangenome or corgenome. Compute the number of genes which belong to
#' the pangenome or to the coregenome, for a number of random permutations of
#' increasingly bigger sample of genomes.
#' @param what One of \code{"pangenome"} or \code{"coregenome"}.
#' @param n.perm The number of permutations to compute (default: 10).
#' @return A \code{matrix}, rows are the number of genomes added, columns are
#' permutations, and the cell number is the number of genes in each category.
rarefact = function(pan, what = 'pangenome', n.perm = 20){
  what <- match.arg(what, c('pangenome', 'coregenome'))
  pm <- pan
  pm[which(pm>1L, arr.ind = TRUE)] <- 1L
  norgs <- nrow(pan)
  rmat <- matrix(0L, nrow = norgs, ncol = n.perm)
  if (what=='pangenome'){
    for (i in seq_len(n.perm)){
      cm <- apply(pm[sample(norgs), ], 2, cumsum)
      rmat[, i] <- rowSums(cm > 0)
    }
  }else{
    sq <- seq_len(norgs)
    for (i in seq_len(n.perm)){
      cm <- apply(pm[sample(norgs), ], 2, cumsum)
      cr <- apply(cm, 2, `==`, sq)
      rmat[, i] <- rowSums(cr)
    }
  }
  rownames(rmat) <- seq_len(norgs)
  colnames(rmat) <- paste0('permut_', seq_len(n.perm))
  rmat
}


rmat <- rarefact(pan, what = 'coregenome', n.perm = 500)

p.cor <- rmat %>% data.frame() %>%
  rownames_to_column(var='n') %>%
  pivot_longer(permut_1:permut_500) %>% mutate(n=as.numeric(as.character(n))) %>%
  ggplot(aes(n,value))+geom_line(aes(group=name),alpha=0.3)+fdb_style(aspect.ratio = 0.5)+xlab('# genomes')+ylab('core genome size')
ggsave('~/DATA/MarbGenomics/Graphs/allMB_core_genome_size.pdf',plot=p.cor, width = 4, height = 3,unit='in')


#ggplot(fitExp$model, aes_string(x = names(fitExp$model)[2], y = names(fitExp$model)[1])) +
#  geom_point() +
#  stat_smooth(method = "lm", col = "red") +
#  labs(title = paste("Adj R2 = ",signif(summary(fitExp)$adj.r.squared, 5),
#                     "Intercept =",signif(fitExp$coef[[1]],5 ),
#                     " Slope =",signif(fitExp$coef[[2]], 5),
#                     " P =",signif(summary(fitExp)$coef[2,4], 5)))


#' @description
#' Compute distance between all pairs of genomes. The default dist method is
#' \code{"bray"} (Bray-Curtis distance). Another used distance method is \code{"jaccard"},
#' but you should set \code{binary = FALSE} (see below) to obtain a meaningful result.
#' See \code{\link[vegan]{vegdist}} for details, this is just a wrapper function.
#' @param method The distance method to use. See \link[vegan]{vegdist}
#' for available methods, and details for each one.
#' @param binary Transform abundance matrix into a presence/absence
#' matrix before computing distance.
#' @param diag Compute diagonals.
#' @param upper Return only the upper diagonal.
#' @param na.rm Pairwise deletion of missing observations when
#' computing dissimilarities.
#' @param ... Other parameters. See \link[vegan]{vegdist} for details.
#' @return A \code{dist} object containing all pairwise dissimilarities between genomes.
dist = function(method = 'bray',
                binary = FALSE,
                diag = FALSE,
                upper = FALSE,
                na.rm = FALSE,
                ...){

  METHODS <- c("manhattan","euclidean","canberra","bray",
               "kulczynski","jaccard","gower","altGower",
               "morisita","horn","mountford","raup",
               "binomial","chao","cao","mahalanobis")
  method <- match.arg(method, METHODS)

  if (method == 'jaccard' & binary == FALSE){
    warning('It is recommended to set binary = TRUE when running dist(method = "jaccard")')
  }

  #vegan::vegdist()
  vegan::vegdist(self$pan_matrix,
                 method = method,
                 diag = diag,
                 upper = upper,
                 na.rm = na.rm,
                 ...)
}


#' @description
#' Fits a power law curve for the pangenome rarefaction simulation.
#' @param raref (Optional) A rarefaction matrix, as returned by \code{rarefact()}.
#' @param ... Further arguments to be passed to \code{rarefact()}. If \code{raref}
#' is missing, it will be computed with default arguments, or with the ones provided here.
#' @return A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
#' with fitted parameters. An attribute \code{"alpha"} is also returned (If
#' \code{alpha>1}, then the pangenome is closed, otherwise is open.
pg_power_law_fit = function(raref, ...){
  # #micropan::heaps()
  # heaps(self$pan_matrix ,n.perm = n.perm)
  if (missing(raref)){
    raref <- self$rarefact(...)
  }
  rm <- melt(raref)
  # Power law linearization:
  # y = K * x ^ delta ==> log(y) = log(K) + delta * log(x)
  fitHeaps <- lm(log(rm$value) ~ log(rm$Var1))
  logK <- summary(fitHeaps)$coef[1]
  delta <- summary(fitHeaps)$coef[2]
  K <- exp(logK)
  alpha <- 1-delta
  ret <- list(formula=NULL,
              params=NULL)
  ret$formula <- function(x) K * x ^ delta
  ret$params <- c(K = K, delta = delta)
  attr(ret, 'alpha') <- alpha
  ret
}

#' @description
#' Fits an exponential decay curve for the coregenome rarefaction simulation.
#' @param raref (Optional) A rarefaction matrix, as returned by \code{rarefact()}.
#' @param pcounts An integer of pseudo-counts. This is used to better fit the function
#' at small numbers, as the linearization method requires to subtract a constant C, which is the
#' coregenome size, from \code{y}. As \code{y} becomes closer to the coregenome size, this operation
#' tends to 0, and its logarithm goes crazy. By default \code{pcounts=10}.
#' @param ... Further arguments to be passed to \code{rarefact()}. If \code{raref}
#' is missing, it will be computed with default arguments, or with the ones provided here.
#' @return A \code{list} of two elements: \code{$formula} with a fitted function, and \code{$params}
#' with fitted intercept and decay parameters.
cg_exp_decay_fit = function(raref, pcounts = 10, ...){
  # Exponential decay linearization:
  # y = A * exp(K * t) + C ==> y - C = K*t + log(A)
  # C == core size
  if (missing(raref)){
    raref <- self$rarefact(what = 'core', ...)
  }
  rm <- melt(raref)
  C <- min(rm$value)
  fitExp <- lm(log(value - C + pcounts) ~ Var1, data= rm)
  A <- exp(summary(fitExp)$coef[1])
  B <- summary(fitExp)$coef[2]
  ret <- list(formula=NULL,
              params=NULL)
  ret$formula <- function(x) A * exp(B * x) + C
  ret$params <- c(A = A, B = B, C = C)
  attr(ret, 'pseudo_counts') <- pcounts
  ret
}

#-----------------------------------------------#

gg_curves = function(what = c('pangenome', 'coregenome'),
                     ...){
  what <- match.arg(what, c('pangenome', 'coregenome'), several.ok = TRUE)
  names(what) <- what
  lrar <- lapply(what, function(x) self$rarefact(what = x))
  lfun <- lapply(what, function(x){
    if (x == 'pangenome'){
      self$pg_power_law_fit(raref = lrar[[x]])#$formula
    }else{
      self$cg_exp_decay_fit(raref = lrar[[x]])#$formula
    }
  })
  lrarm <- lapply(lrar, melt)
  ll <- lapply(what, function(x) {
    olst <- list(data = NULL, formula = NULL)
    lrarm[[x]]$Category <- x
    lrarm[[x]]
  })
  data <- Reduce(rbind, ll)

  #plot
  g <- ggplot(data, aes(x=Var1, y=value, colour=Category)) +
    xlab('Number of genomes') +
    ylab('Number of clusters')
  for (i in seq_along(what)){
    g <- g +
      stat_function(data = data[which(data$Category == what[i]), ],
                    fun = lfun[[what[i]]]$formula)
  }
  g
}













tcdb.pan
cazy.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='feature') %>%
  group_by(genome) %>%
  mutate(total = sum(value)) %>%
  distinct(genome, .keep_all=TRUE)%>%
  left_join(metadata,by=c('genome'='genome')) %>%
  ungroup() %>%
  #filter out the NA ones, perhaps sometimes interesting!
  filter(!(group=="NA") | !(is.na(group))) %>%
  ggplot(aes(x=group,y=total,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  ylab('total feautres detected')+
  fdb_style(aspect.ratio=0.5)



OG.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='feature') %>%
  group_by(genome) %>%
  mutate(total = sum(value)) %>%
  distinct(genome, .keep_all=TRUE)%>%
  left_join(metadata,by=c('genome'='genome')) %>%
  ungroup() %>%
  #filter out the NA ones, perhaps sometimes interesting!
  filter(!(group=="NA") | !(is.na(group))) %>%
  #in that case one could do
  mutate(ratio = total/predicted_genes  ) %>%
  ggplot(aes(x=group,y=ratio,fill=group))+
    geom_boxplot(outlier.shape=NA,alpha=.4)+
    geom_jitter(shape=21,size=2,width = 0.2)+
    ylab('feautres/nr predicted_genes')+
    fdb_style(aspect.ratio=0.5)


OG.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='feature') %>%
  group_by(genome) %>%
  mutate(total = sum(value)) %>%
  distinct(genome, .keep_all=TRUE)%>%
  left_join(metadata,by=c('genome'='genome')) %>%
  ungroup() %>%
  filter(!(group=="NA") | !(is.na(group))) %>%
  mutate(ratio = total/predicted_genes  ) %>% filter(ratio<0.94) %>% select(genome,group,ratio)


# These I would like to inspect a bit more!
annotation.tbl %>%
  filter(genome=='Marinobacter_sp_AC_23') %>%
  select(seqnames, genome, locus_tag, OG, product, cazy_domain, COG, kfm_domain, tcdb_domain) %>%
  filter(is.na(OG)) %>%
  group_by(seqnames) %>% count() %>% data.frame

annotation.tbl %>%
  filter(genome=='Marinobacter_sp_AC_23') %>%
  select(seqnames, genome, locus_tag, OG, product, cazy_domain, COG, kfm_domain, tcdb_domain) %>%
  filter(is.na(OG)) %>% select(seqnames,locus_tag,product) %>% filter(seqnames=='MBPP01000001.1') %>% data.frame()




# SOM INTERESTING PAN MINGLING GOING ON HERE

# SHOW GH SPECIFICALLY
cazyCat.pan = cazyCat.pan[rownames(metadata),]
cazyCat.pan %>% rownames_to_column(var='genome') %>%
  left_join(metadata,by=c('genome'='genome')) %>%
  ggplot(aes(x=group,y=GH,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  fdb_style(aspect.ratio=0.5)

# SHOW all categories
cazyCat.pan %>% rownames_to_column(var='genome') %>%
  left_join(metadata,by=c('genome'='genome')) %>%
  pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT') %>%
  ggplot(aes(x=group,y=value,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  facet_wrap(~cazyCAT,scales='free_y',ncol=1)+
  fdb_style(aspect.ratio=0.3)

# SHOW all bloody domains
cazyCat.pan %>% rownames_to_column(var='genome') %>%
  left_join(metadata,by=c('genome'='genome')) %>%
  pivot_longer(c(GH,GT,PL,CE,CBM),names_to='cazyCAT') %>%
  dplyr::select(genome,cazyCAT, group, value) %>%
  group_by(genome) %>%
  mutate(total = sum(value)) %>% distinct(genome, .keep_all=TRUE)%>%
  ungroup() %>%
  ggplot(aes(x=group,y=total,fill=group))+
  geom_boxplot(outlier.shape=NA,alpha=.4)+
  geom_jitter(shape=21,size=2,width = 0.2)+
  ylab('total CAZY domains detected')+
  fdb_style(aspect.ratio=0.5)




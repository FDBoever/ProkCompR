testMCP <- cbind('x'=c(as.matrix(ANIb.sel)), 'y' =c(GRR.sel)) %>% data.frame()

testMCP %>% data.frame() %>% ggplot(aes(x, y)) + geom_point()


library(mcp)


model = list(
  y~1+x,
  1 + (1|id) ~ 0 + x
)

ex = mcp_example('varying')
fit = mcp(model, ex$data)
plot(fit, facet_by = 'id')


model = list(
  y~1, #intercept
  ~ 0 + x #joint slope
)
# using distance based measures, rather then gene alignments
fit = mcp(model, testMCP)
plot(fit, facet_by = 'id')


model=list(response~1,~0+time,~1+time)
ex=mcp_example('demo')
fit=mcp(model,data=ex$data)


#	loadTree
#=================================================================================#
#' function for loading phylogenetic tree
#'
#' @param inPath specifies the path to the input tree
#' @param clean_genome_names boolean to specify whether to clean the genome names converting unwanted characters to _
#'
#' @return tree in ape phylo object
#' @export
#'
#' @examples
#' tree = loadTree(inPath='~/DATA/MarbGenomics/SCO_genefpair_trimAl/concat.treefile', clean_genome_names=TRUE)
#' tree$tip.label[tree$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"

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


# loadpyANIb
#=================================================================================#
#' reads pyANI output data
#'
#' @param inPath specifies the path to the pyANI output directory (not the file itself)
#' @param clean_genome_names boolean to specify whether to clean the genome names
#'
#' @return dataframe with ANIb results
#' @export
#'
#' @examples
#' ANIb <- loadpyANIb(inPath='~/DATA/MarbGenomics/ANIb_output',clean_genome_names = TRUE)
#' ANI16S = loadpyANIb(inPath='~/DATA/MarbGenomics/ANIb_16Sseq_output',clean_genome_names = TRUE)

loadpyANIb <- function(inPath='', clean_genome_names=TRUE){
  if(inPath==''){
    print("You'd need to specify a directory using inPath")
  }else{
    ANIb = read.delim(paste0(inPath,'/ANIb_percentage_identity.tab'))
    rownames(ANIb)= ANIb[,1]
    ANIb = ANIb[,c(2:(ncol(ANIb)))]

    if(clean_genome_names==TRUE){
      rownames(ANIb) = gsub('-','_',gsub('\\.','_', rownames(ANIb)))
      colnames(ANIb) = gsub('-','_',gsub('\\.','_', colnames(ANIb)))
    }
    return(ANIb)
  }
}


#=================================================================================#


# loadCompareM
#   comparem aai_wf ./Assemblies/ ./CompareM_out/ --cpus 16
#   comparem codon_usage ./Assemblies/ ./CompareM_out/codon_usage --cpus 16
#   comparem aa_usage ./prokka_faa/ ./CompareM_out/aa_usage --cpus 16
#   comparem stop_usage ./Assemblies/ ./CompareM_out/stop_usage --cpus 16
#   comparem kmer_uasge ./Assemblies/ ./CompareM_out/kmer_profile --cpus 16
#=================================================================================#
# ?load the CompareM data (Get rid of the spaces in headers in text-editor!?, I may have fixed this by removing first line!)

#' Loading CompareM resusults
#'
#' @param inPath specifies the path to the CompareM output directory (not the file itself)
#' @param metric to define the CompareM analysis output. can be either 'AAI', 'aa_usage', 'codon_usage' or 'kmer_usage'
#' @param clean_genome_names boolean to specify whether to clean the genome names
#'
#' @return dataframe with the requested CompareM metric data
#' @export
#'
#' @examples
#' AAI = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out/aai', metric='AAI', clean_genome_names=TRUE)
#' aa_usage = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out', metric='aa_usage', clean_genome_names=TRUE)
#' codon_usage = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out', metric='codon_usage', clean_genome_names=TRUE)
#' kmer_usage = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out', metric='kmer_usage', clean_genome_names=TRUE)

loadCompareM <- function(inPath='', metric='AAI', clean_genome_names=TRUE){
  if(inPath==''){
    print("You'd need to specify a directory using inPath")
  }else{
    if(metric=='AAI'){
      CompareM = read.delim(paste0(inPath,'/aai_summary.tsv'),header=FALSE,skip=1)
      colnames(CompareM) = c('Genome_A','Genes_in_A','Genome_B','Genes_in_B','orthologous_genes','Mean_AAI','Std_AAI','Orthologous_fraction')
    }
    if(metric=='aa_usage'){
      CompareM = read.delim(paste0(inPath,'/codon_usage.txt'),header=TRUE,sep='\t')
      rownames(CompareM) = CompareM[,1]
      CompareM = CompareM[,c(2:ncol(CompareM))]
    }
    if(metric=='codon_usage'){
      CompareM = read.delim(paste0(inPath,'/codon_usage.txt'),header=TRUE,sep='\t')
      rownames(CompareM) = CompareM[,1]
      CompareM = CompareM[,c(2:ncol(CompareM))]
    }
    if(metric=='kmer_usage'){
      CompareM = read.delim(paste0(inPath,'/kmer_profile.txt'),header=TRUE,sep='\t')
      rownames(CompareM) = CompareM[,1]
      CompareM = CompareM[,c(2:ncol(CompareM))]
    }

    if(clean_genome_names==TRUE){
      if(metric=='AAI'){
        CompareM$Genome_A = gsub('-','_',gsub('\\.','_', CompareM$Genome_A))
        CompareM$Genome_B = gsub('-','_',gsub('\\.','_', CompareM$Genome_B))
      }
      if(metric %in% c('aa_usage','codon_usage','kmer_usage')){
        rownames(CompareM) = gsub('-','_',gsub('\\.','_', rownames(CompareM)))
        colnames(CompareM) = gsub('-','_',gsub('\\.','_', colnames(CompareM)))
      }
    }
    return(CompareM)
  }
}


#=================================================================================#


##### - not sure
#!!!!!!!!!!!!!!!!!!
# NEED TO SPLIT GENOME_A, GENOME_B, not symmetric!
# Currently works with igraph functions, if this is annoying at a certain point, consider alternatives, perhaps with dplyr
# AAI <- comparem %>%
#   dplyr::group_by(Genome_A) %>% select(Genome_A, Genome_B, !!as.name("Mean_AAI"))%>%
#   tidyr::spread(Genome_B,!!as.name("Mean_AAI"))

#spreadCompareM
#=================================================================================#
#' To convert in a distance formatted table
#'
#' @param tbl table object obtained from loadCompareM() AAI metric
#' @param feature.col to define the column in the original table, for exmaple 'Mean_AAI', or 'Orthologous_fraction'
#'
#' @return a wide formatted pair-wise AAI table
#' @export
#'
#' @examples
#' AAI = spreadCompareM(tbl=comparem, feature.col='Mean_AAI')
#' Ortho.Fraction = spreadCompareM(tbl=comparem, feature.col='Orthologous_fraction')

spreadCompareM <- function(tbl, feature.col='Mean_AAI'){
  dat = tbl[,c("Genome_A","Genome_B",feature.col)]
  g <- igraph::graph.data.frame(dat, directed=FALSE)
  output = igraph::get.adjacency(g, attr=feature.col, sparse=FALSE)

  if(feature.col == 'Mean_AAI'){
    output[output<1] <- 100
  }
  return(output)
}


#interaDistances
#=================================================================================#
#' Function to asses inter and intra-group distances based on pair-wise distance-like tables
#'
#' @param df.dist defines a pairwise, distance-like table with row and column names as genomes
#' @param var_name sets a name, this could be the name of the distance metric
#' @param metadata defines the metadata table from which the group will be extracted
#' @param grp defines the column name of the grouping variable as in the metadata table
#'
#' @return output dataframe with the within and between group distances
#' @export
#'
#' @examples
#' ANIb.intra <- interaDistances(df.dist = ANIb, var_name = "ANIb", metadata = metadata, grp = "group")
#' AAI.intra <- interaDistances(df.dist = AAI, var_name = "AAI",  metadata = metadata, grp = "group")

interaDistances <- function(df.dist, var_name, metadata, grp){
  output <- df.dist %>% data.frame() %>%
    tibble::rownames_to_column(var='genome') %>%
    tidyr::pivot_longer(cols = (-genome), names_to = 'genome_b',values_to= 'value' ) %>%
    dplyr::left_join(metadata %>% select(genome,!!as.name(grp)), by=c('genome'='genome')) %>%
    dplyr::left_join(metadata %>% select(genome,!!as.name(grp)), by=c('genome_b'='genome')) %>%
    dplyr::mutate(grouping =  ifelse(!!as.name(paste0(grp, '.x')) == !!as.name(paste0(grp, '.y')), 'Intra','Inter')) %>%
    dplyr::rename(., group.x = !!as.name(paste0(grp, '.x'))) %>%
    dplyr::rename(., group.y = !!as.name(paste0(grp, '.y'))) %>%
    dplyr::mutate(itself =  ifelse(genome == genome_b, TRUE, FALSE))

    return(output)
}


# Convert to Long format
#concat.sco.nuc = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
#concat.sco.nuc$tip.label[concat.sco.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
#tree = concat.sco.nuc
#tree = nj.AAI
#tree = ape::drop.tip(tree, tip = tree$tip.label[!grepl("Marinobacter|Tamil",tree$tip.label)])


#clusterCophenetic
#=================================================================================#
#' assign a phylogenetic cluster to each tip objectively
#' Inspired by Levy et al 2017
#' It first calculates the cophenetic distances between each pair of tips
#' before clustering, we estimate the optimal number of clusters (k) based on max width of Silhouette
#' next, it clusters tips based on cophenetic distance using PAM clustering
#' NOTE: function depends on ape::cophenetic, factoextra::fviz_nbclust, and fpc::pamk
#' @param tree phylogenetic tree in ape::phylo format
#'
#' @return assigned cluster for each tip
#' @export
#'
#' @examples
#' test.x <- clusterCophenetic(tree)

clusterCophenetic <- function(tree){
  #Root tree using the outgroup?
  #Calculate cophenetic distance
  COPH = ape::cophenetic.phylo(tree)
  silhoutte <- factoextra::fviz_nbclust(COPH, kmeans, method = "silhouette",k.max=25)+
    labs(subtitle = "Silhouette method")

  df_silh <- silhoutte$data
  optimal_k_silh <- as.numeric(as.character(df_silh[which.max(df_silh$y),]$clusters))
  optimal_silh_width <- as.numeric(as.character(df_silh[which.max(df_silh$y),]$y))

  #For testing purposes
  #optimal_k_silh = 20
  pam_out = fpc::pamk(COPH,krange=optimal_k_silh)
  pam_cl = data.frame(pam_out$pamobject$clustering)
  colnames(pam_cl) = c('cluster')
  pam_cl$cluster = paste0('cl',pam_cl$cluster)
  return(pam_cl)
}

#clusterCophenetic_multi
#=================================================================================#
#' function to run clusterCophenetic on multiple trees at once
#'
#' @param treelist a list of trees
#'
#' @return clusters assigned to each tip in each tree
#' @export
#'
#' @examples
#' clusterCophenetic_multi(lsTrees[1:3])

clusterCophenetic_multi <- function(treelist){
  nms <- names(treelist)
  if(length(treelist)>1){
    output <- clusterCophenetic(treelist[[1]])
    s <- sort(rownames(output))
    output <- output[s,]
    q<-2
    for(i in treelist[2:length(treelist)]){
      icl <- clusterCophenetic(i)
      output <- cbind(output,icl[s,])
      colnames(output)[q] <- paste0('tr',q)
      q<-q+1
      }
    rownames(output) = s
    colnames(output)[1] <- paste0('tr',1)
    output = data.frame(output)
    colnames(output) <- nms
    return(output)
    }
}

#multi.drop.tip
#=================================================================================#
#' drop.tip for list of symmetric trees
#'
#' @param treelist list of trees
#' @param tip vector specifying which tip(s) to drop
#'
#' @return a pruned tree (ape::phylo) object
#' @export
#'
#' @examples
#' outgroup <- concat.sco.nuc$tip.label[!grepl("Marinobacter|Tamil",concat.sco.nuc$tip.label)]
#' lsTrees.ingroup <- multi.drop.tip(treelist = lsTrees[1:3], tip = outgroup)
#' clusterCophenetic_multi(lsTrees.ingroup)

multi.drop.tip <- function(treelist, tip){
  nms <- names(treelist)
  lstrees <- list()
  q<-1
  for(i in treelist){
    tre = ape::drop.tip(i, tip = tip)
    lstrees[[q]] <- tre
    q<-q+1
  }
  names(lstrees) <- nms
  class(lstrees) = 'multiPhylo'
  return(lstrees)
}


#multi.root
#=================================================================================#
#' Root multiple trees at ones given an outgroup
#'
#' @param treelist list of trees
#' @param outgroup vector specifying tip names of the outgroup
#'
#' @return treelist object
#' @export
#'
#' @examples
#' lsTrees.rooted <- multi.root(lsTrees, outgroup = outgroup)

multi.root <- function(treelist, outgroup){
  lstrees <- list()
  nms <- names(treelist)
  q=1
  for(i in treelist){
    tre = ape::root(i, outgroup)
    tre = ape::ladderize(tre,right=FALSE)
    lstrees[[q]] <- tre
    q<-q+1
  }
  names(lstrees) <- nms
  class(lstrees) = 'multiPhylo'
  return(lstrees)
}


#check.monoph()
#=================================================================================#
#' Checks if groups are monophyletic given a input tree
#'
#' @param tree phylogenetic tree in ape::phylo format
#' @param taxonomy table with tip.labels and annotation
#'
#' @return result output table
#' @export
#'
#' @examples
#' check.monoph(tree,taxonomy)

check.monoph <- function(tree, taxonomy){
  output <- NULL
  taxa = taxonomy$taxa
  rownames(taxonomy) = taxa
  for(col in colnames(taxonomy)[2:ncol(taxonomy)]){
    message(paste0('--Analysing',col))
    sbdf <- taxonomy[,c('taxa',col)]

    sbdf <- sbdf[!is.na(sbdf[,col]),] #drop NA
    # I COULD IMPLEMENT HERE A ALTERNATIVE METHOD, TO ALSO EXCLUDE NAs (singleton genomes) FROM THE TREE

    sbdf[,col] <- as.factor(sbdf[,col])
    for(i in unique(sbdf[,2])){
      tips <- as.character(sbdf[sbdf[,2]==i,]$taxa)
      #If taxon contains more than one tip check monophyly
      if(length(tips)>1){
        monophyletic <- ape::is.monophyletic(phy=tree, tips= tips)
        output= rbind(output,c('level'=col,'taxon'=i, 'mono'=monophyletic,'nroftips'= length(tips)))
        message(monophyletic)

        #If not monophyletic, this initial screen looks for those tips that upon removal turn the taxon monophyletic
        if(monophyletic==FALSE){
          for(tip in tips){
            MFL= is.monophyletic(phy = drop.tip(tree, tip), tips = tips[tips!=tip])
            if(MFL == TRUE){
              #message(paste('removing', tip, 'turns', i , 'monophyletic'))
            }
          }
        }

      #if taxon contains only one tip, it assigns NA
      }else{
        output= rbind(output,c('level'=col,'taxon'=i, 'mono'=NA,'nroftips'= length(tips)))
      }
    }
  }
  output <- data.frame(output,stringsAsFactors = FALSE)
  output <- output[output$level !='taxa',]
  return(output)
}


#robustMonophyly
#=================================================================================#
#' assess monophyly in multiple trees, by running over all nodes of a master tree (the first in a treelist) to all nodes in all other trees
#' for each node, the offspring tips are selected and checked for monophyly in other trees
#' function returns a handy tibble that can be used as input for ggtree straight away
#' @param treelist
#'
#' @return
#' @export
#'
#' @examples
#' mono.tdtree <- robustMonophyly(lsTrees.rooted)
#' ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=out.tbl%>%dplyr::filter(isTip==FALSE),aes(color=support))


robustMonophyly <- function(treelist){
  A.tree = treelist[[1]]
  A.ggtree <- ggtree(A.tree)
  A.tbl = A.tree %>%  tidytree::as_tibble()
  A.tbl <- A.tbl %>% left_join(A.ggtree$data[,c('node','isTip','x','y','branch','angle')] , by='node')

  #Run over all the nodes, and compare monophyly in trees
  internalNodes <- A.tbl %>% filter(isTip != TRUE) %>% filter(label!='Root') %>% select(node) %>% pull()

  dmono <-NULL
  for(inNode in internalNodes){
    #Detect offspring in tree A
    offspring <- A.tbl %>% tidytree::offspring(inNode) %>%
      dplyr::filter(isTip==TRUE) %>%
      dplyr::select(label) %>% dplyr::pull()

    #Assess whether monophyletic in trees
    A_mono <- ape::is.monophyletic(phy=A.tree,tips=offspring)
    #B_mono <- ape::is.monophyletic(phy=B.tree,tips=offspring)

    v_mono_all <- NULL
    v_nms <- NULL
    for(n.tree in 2:length(treelist)){
      B_mono <- ape::is.monophyletic(phy=treelist[[n.tree]],tips=offspring)
      v_mono_all <- c(v_mono_all, B_mono)
      v_nms <- c(v_nms, paste0('mono_',n.tree))

    }
    names(v_mono_all) <- v_nms

    #Store the output data
    vmono <- c('inNode' = inNode, 'mono_1' = A_mono, v_mono_all , n_offspring = length(offspring))
    dmono <- rbind(dmono,vmono)
  }

  #combine and tidy
  dmono <- dmono %>% data.frame()
  n.trees <- dmono %>% select(starts_with("mono_")) %>% ncol()
  dmono <- dmono %>% mutate(positive = rowSums(select(., starts_with("mono_")))) %>%
    mutate(support = positive/n.trees)
  out.tbl <- A.tbl %>% left_join(dmono,by=c('node'='inNode'))

  return(out.tbl)
}



#=================================================================================#
#### RUN CODE ####
#=================================================================================#

#---- LOAD ANIb FILES
#=================================================================================#
#pyANI results are loaded in using the loadpyANIb function, in case of Marinobacter clean_genome_names is used

ANIb <- loadpyANIb(inPath='~/DATA/MarbGenomics/ANIb_output',
                   clean_genome_names = TRUE)
ANI16S = loadpyANIb(inPath='~/DATA/MarbGenomics/ANIb_16Sseq_output',
                    clean_genome_names = TRUE)

#---- LOAD AAI and other CompareM data
#compareM restuls are loaded in using the loadCompareM function

comparem = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out/aai',
                        metric='AAI',
                        clean_genome_names=TRUE)

aa_usage = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out',
                        metric='aa_usage',
                        clean_genome_names=TRUE)

codon_usage = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out',
                           metric='codon_usage',
                           clean_genome_names=TRUE)

kmer_usage = loadCompareM(inPath='/Users/sa01fd/DATA/MarbGenomics/CompareM_out',
                          metric='kmer_usage',
                          clean_genome_names=TRUE)

#-- Exytacy AAI and OrthoFraction data and convert to pair-wise distance like table
# AAI and Ortho.Fraction are turned into pairwise-genome to genome format using the speadCompareM function
AAI = spreadCompareM(tbl=comparem, feature.col='Mean_AAI')
AAI <- data.frame(AAI/100)

Ortho.Fraction = spreadCompareM(tbl=comparem, feature.col='Orthologous_fraction')
Ortho.Fraction <- data.frame(Ortho.Fraction/100)



# Functional DISTANCES
#=================================================================================#
# generate pair-wise functional distance matrices from feature abundance tables

# Jaccard distance
#--------------------------------------------------#
# distJac <- micropan::distJaccard(OG.pan)
# distJac= as.matrix(distJac)
# distJac = distJac[ order(row.names(distJac)), ]
# distJac = distJac[ , order(colnames(distJac))]

#In case you want to subset it just now
#subDistJac = distJac[grepl('Marinobacter|Tamil',rownames(distJac)),grepl('Marinobacter|Tamil',colnames(distJac))]

# Manhattan distance
#--------------------------------------------------#
distMan <- micropan::distManhattan(OG.pan)
distMan = as.matrix(distMan)
distMan = distMan[ order(row.names(distMan)), ]
distMan = distMan[ , order(colnames(distMan))]

#In case you want to subset it just now
#subDistMan = distMan[grepl('Marinobacter|Tamil',rownames(distMan)),grepl('Marinobacter|Tamil',colnames(distMan))]


# ==== ANALYSIS INCLUDING THE OUTGROUP!

# Phylogenetic trees
#=================================================================================#
# Phylogenetic trees are loaded with the loadTree function

#tree
concat.rps <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins/concat.treefile', clean_genome_names=TRUE)
concat.rps.nuc <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins_pal2nal/concat.treefile', clean_genome_names=TRUE)
concat.sco <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile', clean_genome_names=TRUE)
concat.sco.nuc <- loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
concat.rDNA <- loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees//concat.treefile', clean_genome_names=TRUE)
rDNA16S <- loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees/16S.treefile', clean_genome_names=TRUE)

#Idialy find a solution for this, as this really is odd...
concat.rps$tip.label[concat.rps$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.rps.nuc$tip.label[concat.rps.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.sco$tip.label[concat.sco$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.sco.nuc$tip.label[concat.sco.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.rDNA$tip.label[concat.rDNA$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
rDNA16S$tip.label[rDNA16S$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"

#set master tree
tree <- concat.sco.nuc
SelectedTips = tree$tip.label


#----------------------#
#Filter the genomes that are in the master tree
# in the mean time, they are ordered, handy for later

ANIb.f <- ANIb[SelectedTips,SelectedTips]
AAI.f <- AAI[SelectedTips,SelectedTips]

dMan = distMan[SelectedTips,SelectedTips]
dJac = distJac[SelectedTips,SelectedTips]
dGGR = pan.ggr[SelectedTips,SelectedTips]
dOF <- 1-Ortho.Fraction[SelectedTips,SelectedTips]

# Implicit phylogeny, NJ of ANIb and AAI
#=================================================================================#
# implicit phylogenetic method, usng ANIb and AAI data to reconstruct phylogenetic trees using neighbour-joining algorithm
# depends on phangorn::NJ

#calculate NJ from filtered ANI AAI sets, based on tips in master tree
nj.ANIb = phangorn::NJ(1-as.matrix(ANIb.f))
nj.AAI = phangorn::NJ(1-(as.matrix(AAI.f)))
nj.Man = phangorn::NJ(as.matrix(dMan))
nj.GGR = phangorn::NJ(as.matrix(dGGR))
nj.OF = phangorn::NJ(as.matrix(dOF))

#nj.Jac = NJ(as.matrix(dJac))

#nj.GGR =NJ(as.matrix(dGGR))


#filter the real phylogenies
concat.rps = keep.tip(concat.rps, SelectedTips)
concat.rps.nuc = keep.tip(concat.rps.nuc, SelectedTips)
concat.sco = keep.tip(concat.sco, SelectedTips)
concat.sco.nuc = keep.tip(concat.sco.nuc, SelectedTips)

#-------------------------------------------------------#


#Example of a treelist object
lsTrees =
  list('concat.sco.nuc'=concat.sco.nuc,
       'concat.sco'=concat.sco,
       'concat.rps.nuc'=concat.rps.nuc,
       'concat.rps'=concat.rps,
       'nj.ANIb' = nj.ANIb,
       'nj.AAI' = nj.AAI,
       'nj.Man'= nj.Man,
       'nj.GGR' = nj.GGR
  )

#Define the outgroup
outgroup <- concat.sco.nuc$tip.label[!grepl("Marinobacter|Tamil",concat.sco.nuc$tip.label)]

#Root all trees at outgroup
lsTrees.rooted <- multi.root(lsTrees, outgroup = outgroup)

#midpoint?
lsTrees.rooted <- list()
for(t.name in names(lsTrees)){
  lsTrees.rooted[[t.name]] <- phangorn::midpoint(lsTrees[[t.name]])
  #densi.trees[[t.name]]$edge.length  <- NULL
}
lsTrees.rooted
class(lsTrees.rooted) <- 'multiPhylo'


#drop the outgroup
lsTrees.ingroup <- multi.drop.tip(treelist = lsTrees.rooted, tip = outgroup)

#run clustering approach on each tree
multi.phylo.cluster <- clusterCophenetic_multi(lsTrees.ingroup)
s <- rownames(multi.phylo.cluster)
multi.phylo.cluster[ , colnames(multi.phylo.cluster)] <- data.frame(lapply(multi.phylo.cluster[ , colnames(multi.phylo.cluster)] , factor))


for(col in colnames(multi.phylo.cluster)){
  print(col)
  multi.phylo.cluster[,col]
  n <- length(levels(multi.phylo.cluster[,col]))
  multi.phylo.cluster[,col] = factor(multi.phylo.cluster[,col], levels =paste0('cl',1:n))
}

rownames(multi.phylo.cluster) <- s

#visualise
p.tr <- ggtree(lsTrees.rooted[[1]]) + geom_tiplab(size=2)
gheatmap(p.tr, multi.phylo.cluster, offset = 0, width = 1, colnames = TRUE, colnames_position = "bottom", font.size = 4) + scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(25))


#Check Monophyly
monodat <- multi.phylo.cluster %>% rownames_to_column(var='taxa') %>% select(taxa,concat.sco.nuc)
check.monoph(lsTrees.rooted[[1]] ,monodat) %>% mutate(tree=names(lsTrees.rooted)[1])
check.monoph(lsTrees.rooted[[2]] ,monodat) %>% mutate(tree=names(lsTrees.rooted)[2])
check.monoph(lsTrees.rooted[[3]] ,monodat) %>% mutate(tree=names(lsTrees.rooted)[3])
check.monoph(lsTrees.rooted[[4]] ,monodat) %>% mutate(tree=names(lsTrees.rooted)[4])
check.monoph(lsTrees.rooted[[5]] ,monodat) %>% mutate(tree=names(lsTrees.rooted)[5])
check.monoph(lsTrees.rooted[[6]] ,monodat) %>% mutate(tree=names(lsTrees.rooted)[6])


monodat <- multi.phylo.cluster %>% rownames_to_column(var='taxa')
monodat <- multi.phylo.cluster %>% rownames_to_column(var='taxa')

p.tree <- ggtree(lsTrees.rooted[[1]]) + geom_tiplab(size=2)
p.tree$data <- p.tree$data %>%
  left_join(genome.tbl, by=c('label'='genome')) %>%
  left_join(monodat, by=c('label'='taxa'))


p.tree + geom_point(data=p.tree$data[p.tree$data$isTip==TRUE,], aes(color=SS))
p.tree + geom_point(data=p.tree$data[p.tree$data$isTip==TRUE,], aes(color=nj.AAI))
p.tree + geom_point(data=p.tree$data[p.tree$data$isTip==TRUE,], aes(color=nj.ANIb))
p.tree + geom_point(data=p.tree$data[p.tree$data$isTip==TRUE,], aes(color=concat.sco.nuc))


p.tree + geom_point(data=p.tree$data[p.tree$data$isTip==TRUE,], aes(color=concat.sco.nuc))



#------
p.tree <- ggtree(lsTrees.rooted[[1]],layout='circular',branch.length = 'none')
p.tree <- ggtree(lsTrees.rooted[[1]],branch.length = 'none')
p.tree$data <- p.tree$data %>%
  left_join(genome.tbl, by=c('label'='genome')) %>%
  left_join(monodat, by=c('label'='taxa'))




#=
RFoulds <- as.matrix(phangorn::RF.dist(lsTrees.rooted, normalize = TRUE, check.labels = TRUE))








#--------------#
# TO BE IMPLEMENTED
# VISUALISE FOR EACH TREE, WHICH GROUPS ARE IDENTIFIED,
# CHECK MONOPHYELY OF THAT STUFF, AS THIS WILL BE IMPORANT

p.tr <- ggtree(lsTrees.rooted[[1]]) + geom_tiplab(size=2)
gheatmap(p.tr, multi.phylo.cluster[,colnames(multi.phylo.cluster)[1]], offset = 0, width = 1, colnames = TRUE, colnames_position = "bottom", font.size = 4) + scale_fill_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(25))
#=================================================================================#
# PERHAPHS, IMPLEMENT ANOTHER METHOD
# THAT LOOKS FOR DEEPEST MONOPHYLETIC CLADES SUPPORTED BY ALL TREES
# THIS IS A CHALLENCE SO THINK FIRST IF THIS IS BETTER OR WORSE, MORE SENSIBLE OR NOT...
# MAY YIELD VERY LITTLE INFORMATION BACK
#
# WHAT IS KEY TO TRY AND IMPLEMENT IS THE DEFINITION OF ROBUST WELL RESOLVED MONOPHYLETIC NODES FROM WHICH TO DRAW PHYLOGGROUPS
# REPRESENTING DEEP CLADES, IGNORING PHYLOGENETIC DEPTH,
# OR CUTS AT VARIOUS PHYLOGENETIC DEPTHS,
# SAY CUTS THE TREE AT X phylogenetic depht!
# ROBUST MONOPHYLY IS KEY THO!

# CURRENTLY NOT NEEDED, INTERESTING WHEN WANTING TO CROSS COMPARE, OR DO OTHER RELATED COMPARATIVE THINGS
#prepare for the others as well
#lsTreeTbls <- list()
#for(n.tree in 2:length(treelist)){
#  B.tree = treelist[[n.tree]]
#  B.ggtree <- ggtree(B.tree)
#  B.tbl = B.tree %>%  tidytree::as_tibble()
#  B.tbl <- B.tbl %>% left_join(B.ggtree$data[,c('node','isTip','x','y','branch','angle')] , by='node')
#  lsTreeTbls[[n.tree-1]]<- B.tbl
#  }




mono.tdtree <- robustMonophyly(lsTrees.rooted)

mono.tdtree <- robustMonophyly(lsTrees.rooted[c(2,3,5)])

  ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=out.tbl%>%dplyr::filter(isTip==FALSE),aes(color=support))
  ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=out.tbl%>%dplyr::filter(isTip==FALSE)%>%filter(support>0.4),aes(color=support,size=1.5))

  ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=out.tbl%>%dplyr::filter(isTip==FALSE)%>%filter(support>0.6),aes(color=support,size=1.5))

  ggtree(mono.tdtree) + geom_tiplab(size=2) + geom_point(data=out.tbl%>%dplyr::filter(isTip==FALSE)%>%filter(support==1),aes(color=support,size=1.5))









  ggtree(out.tbl,layout='fan',branch.length = NA) + geom_tiplab(size=2) + geom_point(data=out.tbl%>%dplyr::filter(isTip==FALSE)%>%filter(support>0.6),aes(color=support,size=1.5))


  tidytree::as.phylo(out.tbl)


  A.tbl %>% tidytree::child(120)


  concat.sco.nuc.tbl %>% tidytree::parent(120)
  concat.sco.nuc.tbl %>% tidytree::offspring(119)
  concat.sco.nuc.tbl %>% tidytree::ancestor()
  concat.sco.nuc.tbl %>% tidytree::sibling()
  concat.sco.nuc.tbl %>% tidytree::MCRA()

#=================================================================================#



#Intra phylogroup

intraPG <- interaDistances(df.dist = ANIb.f[metadata$genome,metadata$genome],
                                var_name = "ANIb",
                                metadata = metadata,
                                grp = "phylogroup")

intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(group.x, value)) +
  geom_boxplot(aes(color=group.x),outlier.shape=NA) +
  geom_jitter(aes(color=group.x),shape=21,width=0.2,alpha=0.5) +
  facet_wrap(~grouping)+
  scale_color_manual(values=phylogroup.colors) +
  fdb_style()


aov(value ~ grouping, data=intraPG) %>% summary()
aov(value ~ grouping + group.x, data=intraPG) %>% summary()


intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(group.x, value,color=grouping)) +
  #geom_boxplot(aes(),outlier.shape=NA) +
  geom_jitter(shape=21,alpha=0.5,position=position_jitterdodge()) +
  #facet_wrap(~grouping)+
  #scale_color_manual(values=phylogroup.colors) +
  fdb_style()




intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(grouping, value,color=grouping)) +
  #geom_boxplot(aes(),outlier.shape=NA) +
  geom_jitter(shape=21,alpha=0.5,position=position_jitterdodge()) +
  #facet_wrap(~grouping)+
  #scale_color_manual(values=phylogroup.colors) +
  fdb_style()


intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(value,fill=grouping)) +
  geom_density(alpha=0.5)+
  geom_hline(yintercept=0)+
  geom_point(aes(x=value,y=(as.numeric(as.factor(grouping))*5)-13,color=grouping),
             fill=NA,shape=21,alpha=0.5,
             position=position_jitter(w=0,h=1))+
  fdb_style(aspect.ratio = 0.5) + viridis::scale_fill_viridis(discrete = TRUE) + viridis::scale_color_viridis(discrete = TRUE)



#----------------------------------------------------------------------------------#
### SI FIGURE #####
# Density distribution of pair-wise ANIb values showing non homogeneous distribution of intra phylogroup distances between phylroups

#overal plot, global
p1 <- intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(value,fill=grouping)) +
  geom_density(alpha=0.5)+
  geom_hline(yintercept=0)+
  geom_point(aes(x=value,y=(as.numeric(as.factor(grouping))*5)-13,color=grouping),
             fill=NA,shape=3,alpha=0.5,
             position=position_jitter(w=0,h=1))+
  fdb_style(aspect.ratio = 0.5) +
  scale_fill_manual(values=c('purple','forestgreen')) +
  scale_color_manual(values=c('purple','forestgreen'))+scale_x_continuous(limits=c(0.7,1))+ylab('Density')

#plot with ridges per phylogroup
p2 <- intraGenus%>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x= value , y= group.x)) +
  ggridges::geom_density_ridges(mapping = aes(group=interaction(group.x, grouping), fill=grouping),alpha=0.5) +
  geom_point(aes(color=grouping),
             fill=NA,shape=3,alpha=0.5,
             position=position_jitter(w=0,h=0.1))+
  scale_color_manual(values=c('purple','forestgreen'))+
  scale_fill_manual(values=c('purple','forestgreen'))+
  fdb_style(aspect.ratio = 0.5)+
  scale_y_discrete(limits=rev(levels(intraGenus$group.x)))+
  scale_x_continuous(limits=c(0.7,1))+ylab('Phylogroup')

#combine and save
p.o <- ggpubr::ggarrange(p1,p2,labels=c('A','B'),ncol=1,align='v')#,heights = c(1,1.5))
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/inter_intra_non_homogeneous_1.pdf',plot=p.o, width = 8, height = 5,unit='in')


#plot without ridges
p3 <- intraGenus%>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot2::ggplot(ggplot2::aes(x= value , y= group.x)) +
  ggplot2::geom_point(ggplot2::aes(color=grouping),
             fill=NA,shape=3,alpha=0.5,
             position=position_jitterdodge())+
  ggplot2::scale_color_manual(values=c('purple','forestgreen'))+
  ggplot2::scale_fill_manual(values=c('purple','forestgreen'))+
  fdb_style(aspect.ratio = 0.5)+
  ggplot2::scale_y_discrete(limits=rev(levels(intraGenus$group.x)))+
  ggplot2::scale_x_continuous(limits=c(0.7,1))+ylab('Phylogroup')+
  ggplot2::theme(panel.grid.major.y=element_line(color='grey'))

#combine and save
p.o <- ggpubr::ggarrange(p1,p3,labels=c('A','B'),ncol=1,align='v')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/inter_intra_non_homogeneous_2.pdf',plot=p.o, width = 8, height = 5,unit='in')


#-------------------------------------------------------------------------------#




##### RUN THIS CODE ####
#=================================================================================#
#ANI and AAI comparison, intraDistances!

intraGenus <- interaDistances(df.dist = ANIb.f[rownames(metadata),rownames(metadata)],
                var_name = "ANIb",
                metadata = metadata,
                grp = "phylogroup")

#------ describe distibution within groups
ANIb.intra <- interaDistances(df.dist = ANIb,
                              var_name = "ANIb",
                              metadata = metadata,
                              grp = "phylogroup")

AAI.intra <- interaDistances(df.dist = AAI*100,
                             var_name = "AAI",
                             metadata = metadata,
                             grp = "phylogroup")

AAI.intra %>% dplyr::filter(grepl('Marinobacter|Tamil',genome)) %>%
  dplyr::filter(grepl('Marinobacter|Tamil',genome_b)) %>%
  mutate(tamil=grepl('Tamil',genome)) %>% dplyr::filter(itself==FALSE) %>%ggplot(aes(x = value)) +
  geom_histogram(aes(fill = tamil),
                 binwidth = 0.5, color = 'white', alpha = .8, position="identity") +
  scale_y_continuous(expand=c(0,0)) +
  fdb_style(aspect.ratio=0.5) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+scale_fill_manual(values=c('lightgrey','red'))

AAI.intra %>% dplyr::filter(grepl('Marinobacter|Tamil',genome)) %>%
  dplyr::filter(grepl('Marinobacter|Tamil',genome_b)) %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(value)) +
  geom_density() +
  geom_histogram(bins=50) +
  scale_y_continuous(expand=c(0,0)) +
  fdb_style(aspect.ratio=0.5)

mean.val <- AAI.intra %>% dplyr::filter(itself==FALSE) %>% dplyr::filter(!is.na(grouping)) %>% select(value) %>% pull() %>% mean()

AAI.intra %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x = value)) +
  geom_histogram( aes(y = ..density..), binwidth = 0.5, color = 'white', fill = "grey30", alpha = .3) +
  geom_density(alpha = .3, fill = "antiquewhite3",color ='antiquewhite3')+
  scale_y_continuous(expand=c(0,0)) +
  fdb_style(aspect.ratio=0.5) +
  geom_vline(aes(xintercept=0.90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=0.95),color = "grey30", linetype = "dashed")
  #geom_vline(aes(xintercept=mean.val), color = "red", linetype = "dashed")+
  #geom_text(aes(x = mean.val, y = 0.1, label = mean.val), size = 3, hjust = -.1)
  #geom_text(aes(x=77, label="\nthe strong cars", y=0.05), colour="blue", angle=90)


ANIb.intra %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  dplyr::mutate(value=value*100) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.5, color = "grey30", fill = "white") +
  geom_density(alpha = .2, fill = "antiquewhite3")+
  scale_y_continuous(expand=c(0,0)) +
  fdb_style(aspect.ratio=0.5) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")#+
  #geom_vline(aes(xintercept=mean.val), color = "red", linetype = "dashed")+
  #geom_text(aes(x = mean.val, y = 0.1, label = mean.val), size = 3, hjust = -.1)
  #geom_text(aes(x=77, label="\nthe strong cars", y=0.05), colour="blue", angle=90)

AAI.intra %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(grouping, value)) +
  geom_boxplot(aes(color=grouping),outlier.shape=NULL) +
  geom_jitter(aes(color=grouping, fill=grouping),shape=21,width=0.2,alpha=0.5) +
  fdb_style()

AAI.intra %>% dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(grouping, value)) +
  geom_boxplot(aes(color=grouping),outlier.shape=NULL) +
  geom_jitter(aes(color=grouping, fill=grouping),shape=21,width=0.2,alpha=0.5) +
  fdb_style()+facet_wrap(~group.x,nrow = 2)

AAI.intra %>% dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(group.x, value)) +
  geom_point(aes(color=grouping),alpha=0.5, position=position_jitterdodge()) +
  geom_boxplot(aes(color=grouping),outlier.shape=NULL) +
  fdb_style(aspect.ratio = 0.5)

AAI.intra %>% dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(group.x, value,fill=grouping)) +
  geom_bar(aes(x=group.x, y=value,fill=grouping),stat="summary", fun.y="mean",position="dodge",alpha=0.3) +
  geom_point(aes(color=grouping),alpha=0.5, position=position_jitterdodge()) +
  fdb_style(aspect.ratio = 0.5)



#---------
metrics.merged <- rbind(cbind(ANIb.intra%>%mutate(value=value*100), 'metric'='ANIb'),
      cbind(AAI.intra,'metric'='AAI'))


#---
# visualise histogram and overlayed kernal density plot of both ANI and AAI
p.density.all <- metrics.merged %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.5, color = 'white', fill = "grey30", alpha = .3) +
  geom_density(alpha = .3, fill = "antiquewhite3",color ='antiquewhite3')+
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  facet_wrap(~metric,nrow=2, scales='free_y') +
  fdb_style(aspect.ratio=0.25)

#zoomed
p.density.zoom <- metrics.merged %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  dplyr::filter(value>=85)%>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.25, color = 'white', fill = "grey30",alpha = .3) +
  geom_density(alpha = .3, fill = "antiquewhite3" ,color ='antiquewhite3')+
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30" , linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30" , linetype = "dashed")+
  facet_wrap(~metric,nrow=2, scales='free_y') +
  fdb_style(aspect.ratio=0.25)

p.density.c.c <- grid.arrange(p.density.all, p.density.zoom, nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/ANI_AAI_density_grey.pdf',plot=p.density.c.c, width = 8, height = 4,unit='in')


#---#
#In color now

p.density.all <- metrics.merged %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density.., fill = metric),
                 binwidth = 0.5, alpha=.5,color='white') +
  geom_density(aes(fill=metric,color=metric),alpha = .2)+
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  scale_color_manual(values=c('purple','forestgreen'))+
  scale_fill_manual(values=c('purple','forestgreen'))+
  facet_wrap(~metric,nrow=2, scales='free_y') +
  fdb_style(aspect.ratio=0.25)+ theme(legend.position = "none")

#Zoom in onto the region >85
p.density.zoom <- metrics.merged %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x = value)) +
  xlim(c(85,100))+
  geom_histogram(aes(y = ..density.., fill = metric),
                 binwidth = 0.25, alpha=.5,color='white') +
  geom_density(aes(fill=metric,color=metric),alpha = .2)+
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  scale_color_manual(values=c('purple','forestgreen'))+
  scale_fill_manual(values=c('purple','forestgreen'))+
  facet_wrap(~metric,nrow=2, scales='free_y') +
  fdb_style(aspect.ratio=0.25) + theme(legend.position = "none")

p.density.c.c <- grid.arrange(p.density.all, p.density.zoom, nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/ANI_AAI_density_color.pdf',plot=p.density.c.c, width = 8, height = 4,unit='in')




metrics.merged %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x = value,fill=grouping,color=grouping)) +
  #xlim(c(85,100))+
  geom_histogram(aes(y = ..density..),position='identity',
                 binwidth = 0.25, alpha=.5,color='white') +
  geom_density(alpha = .2)+
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  scale_color_manual(values=c('purple','forestgreen'))+
  scale_fill_manual(values=c('purple','forestgreen'))+
  facet_wrap(~metric,nrow=2, scales='free_y') +
  fdb_style(aspect.ratio=0.25) + theme(legend.position = "none")








#----- Check at Tamilnaduibacter, genome, vs Marinobacter genomes, far away!
p.density.tamil <- metrics.merged %>% dplyr::filter(grepl('Marinobacter|Tamil',genome)) %>%
  dplyr::filter(grepl('Marinobacter|Tamil',genome_b)) %>%
  mutate(tamil=grepl('Tamil',genome)) %>% dplyr::filter(itself==FALSE) %>%ggplot(aes(x = value)) +
  geom_histogram(aes(fill = tamil),
                 binwidth = 0.5, color = 'white', alpha = .8, position="identity") +
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  scale_fill_manual(values=c('lightgrey','red')) + facet_wrap(~metric,nrow=2, scales='free_y') +
  ggtitle('Tamilnaduibacter vs Marinobacter') +
  fdb_style(aspect.ratio=0.25) + theme(legend.position = "none")

#----- Check at Tamilnaduibacter, genome, vs Marinobacter genomes, far away!
p.density.sensu_stricto <- metrics.merged %>% dplyr::filter(grepl('Marinobacter|Tamil',genome)) %>%
  dplyr::filter(grepl('Marinobacter|Tamil',genome_b)) %>%
  mutate(cl19=ifelse(group.x=='cl19',TRUE,FALSE)) %>% dplyr::filter(itself==FALSE) %>%
  filter(!is.na(cl19))%>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = cl19),
                 binwidth = 0.5, color = 'white', alpha = .8, position="identity") +
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  scale_fill_manual(values=c('lightgrey','red')) + facet_wrap(~metric,nrow=2, scales='free_y') +
  ggtitle('Marinobacter_A vs Marinobacter sensu stricto') +
  fdb_style(aspect.ratio=0.25) + theme(legend.position = "none")

p.density.c.p <- grid.arrange(p.density.tamil, p.density.sensu_stricto, nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/ANI_AAI_density_sensu_stricto.pdf',plot=p.density.c.p, width = 8, height = 4,unit='in')






#----------#
metrics.merged %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x = value)) +
  xlim(c(85,100))+
  geom_histogram(aes(y = ..density.., fill = metric),
                 binwidth = 0.25, color = "white",alpha=.5) +
  geom_density(aes(fill=metric,color=metric),alpha = .5, position="identity")+
  scale_y_continuous(expand=c(0,0)) +
  fdb_style(aspect.ratio=0.5) +
  geom_vline(aes(xintercept=90),color = "grey30", linetype = "dashed")+
  geom_vline(aes(xintercept=95),color = "grey30", linetype = "dashed")+
  scale_color_manual(values=c('purple','forestgreen'))+
  scale_fill_manual(values=c('purple','forestgreen'))#+
  #facet_wrap(~metric,nrow=2, scales='free_y')

#-------#

metrics.merged %>% dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(group.x, value)) +
  geom_point(aes(color=grouping),alpha=0.5, position=position_jitterdodge()) +
  geom_boxplot(aes(color=grouping),outlier.shape=NULL) +
  fdb_style(aspect.ratio = 0.5) + facet_wrap(~metric,nrow=2, scales='free_y')


#---- Summarise the table
AAI.intra %>% dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  dplyr::group_by(group.x,grouping) %>%
  dplyr::summarise(mean=mean(value),sd=sd(value), n=n(), se=sd/sqrt(n),var=var(value))


#-------------------------------------------------------#

#Make tibble using tidytree as_tibble function
concat.sco.nuc.tbl <- concat.sco.nuc %>% tidytree::as_tibble()
#Add what you want by joining metadata in
concat.sco.nuc.tbl %>% dplyr::left_join(metadata, by =c('label'='genome'))

#as.treedata(concat.sco.nuc.tbl)
#concat.sco.nuc.tbl
#concat.sco.nuc.tbl %>% tidytree::as.treedata() %>% tidytree::as_tibble()

merged = treeio::merge_tree(concat.sco %>% tidytree::as.treedata(),
                   concat.sco.nuc %>% tidytree::as.treedata())



concat.sco.nuc.tbl %>% tidytree::child(120)
concat.sco.nuc.tbl %>% tidytree::parent(120)
concat.sco.nuc.tbl %>% tidytree::offspring(119)
concat.sco.nuc.tbl %>% tidytree::ancestor()
concat.sco.nuc.tbl %>% tidytree::sibling()
concat.sco.nuc.tbl %>% tidytree::MCRA()

concat.sco.nuc.tbl %>% offspring(110)



#back to tree
tidytree::as.phylo(concat.sco.nuc.tbl)



#---------

concat.rps.nuc = keep.tip(concat.rps.nuc, SelectedTips)
concat.sco = keep.tip(concat.sco, SelectedTips)
concat.sco.nuc = keep.tip(concat.sco.nuc, SelectedTips)
nj.ANIb <- keep.tip(nj.ANIb, SelectedTips)
z <- concat.rps.nuc
x <- concat.sco
y <- concat.sco.nuc
y <- nj.ANIb
#y = nj.Man
p1 <- ggtree(x, layout='rectangular')
p2 <- ggtree(y)

d1 <- p1$data
d2 <- fortify(y)
d3 <- fortify(z)
d2$x <- d2$x + max(d1$x) + 1
d3$x <- d3$x + max(d2$x) + 1

dd = bind_rows(d1, d2, d3) %>%
  filter(!is.na(label))

p1 + geom_tree(data = d2) + geom_tree(data = d3) + geom_tiplab(data=d3) +
  geom_line(aes(x, y, group=label, color=node < 15), data=dd, alpha=.3)



#-------------------#

p1 <- ggtree(x, layout='rectangular') +
  geom_hilight(
    mapping=aes(subset = node %in% c(38, 48, 58, 36),
                node = node,
                fill = as.factor(node)
    )
  ) +
  labs(fill = "clades for tree in left" )

p2 <- ggtree(y)

d1 <- p1$data
d2 <- p2$data

#IF YOU WANT PANGENOME
#d2 <- d2 %>% mutate(branch.length=branch.length/100000)
#d2 <- d2 %>% mutate(x=x/5000)

## reverse x-axis and
## set offset to make the tree in the right hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2) +
  ggnewscale::new_scale_fill() +

  labs(fill = "clades for tree in right" )

dd <- bind_rows(d1, d2) %>%
  filter(!is.na(label)) %>% filter(isTip==TRUE) %>% left_join(metadata, by=c('label'='genome'))

pp + geom_line(aes(x, y, group=label,  color=sourceClass), data=dd) +
  geom_tiplab(size=3) +
  geom_tiplab(data = d2, size=3, hjust=1)

pp + geom_line(aes(x, y, group=label,  color=group), data=dd) +
  geom_point(data = dd,aes(x, y, fill=group), size=2,shape=21)+
  ggtitle('species tree vs ANI.tree')+coord_flip()

pp + geom_line(aes(x, y, group=label,  color=SS), data=dd) +
  geom_point(data = dd,aes(x, y, fill=SS), size=2,shape=21)+
  ggtitle('species tree vs ANI.tree')+coord_flip()

pp + geom_line(aes(x, y, group=label,  color=GC), data=dd) +
  geom_point(data = dd,aes(x, y, fill=GC), size=2,shape=21)+
  ggtitle('species tree vs ANI.tree')+coord_flip()




#-------------------#
#================================================================================
SelectedTips<-sel.genome

#----------------------#
#Filter the genomes that are in the master tree
# in the mean time, they are ordered, handy for later
ANIb.f <- ANIb[SelectedTips,SelectedTips]
AAI.f <- AAI[SelectedTips,SelectedTips]
dMan = distMan[SelectedTips,SelectedTips]
#dJac = distJac[SelectedTips,SelectedTips]
dGGR = pan.ggr[SelectedTips,SelectedTips]
dOF <- Ortho.Fraction[SelectedTips,SelectedTips]


d_kmer_usage <- vegan::vegdist(kmer_usage[SelectedTips,],method='bray', diag=TRUE, upper=TRUE)
d_aa_usage <- vegan::vegdist(aa_usage[SelectedTips,],method='bray', diag=TRUE, upper=TRUE)
d_codon_usage <- vegan::vegdist(codon_usage[SelectedTips,],method='bray', diag=TRUE, upper=TRUE)





# Implicit phylogeny, NJ of ANIb and AAI
#=================================================================================#
# implicit phylogenetic method, usng ANIb and AAI data to reconstruct phylogenetic trees using neighbour-joining algorithm
# depends on phangorn::NJ

#calculate NJ from filtered ANI AAI sets, based on tips in master tree
nj.ANIb = phangorn::NJ(1-as.matrix(ANIb.f))
nj.AAI = phangorn::NJ(1-(as.matrix(AAI.f)))
nj.Man = phangorn::NJ(as.matrix(dMan))
nj.GGR = phangorn::NJ(1-as.matrix(dGGR))
nj.OF = phangorn::NJ(1-as.matrix(dOF))
nj.kmer = phangorn::NJ(as.matrix(d_kmer_usage))
nj.aa_usage = phangorn::NJ(as.matrix(d_aa_usage))
nj.codon_usage = phangorn::NJ(as.matrix(d_codon_usage))

nj.kmer = phangorn::NJ(as.matrix(d_kmer_usage))
nj.aa_usage = phangorn::NJ(as.matrix(d_aa_usage))
nj.codon_usage = phangorn::NJ(as.matrix(d_codon_usage))
nj.KO=phangorn::NJ(as.matrix(kfm.dist[sel.genome,sel.genome]))
nj.CAZY=phangorn::NJ(as.matrix(cazy.dist[sel.genome,sel.genome]))
nj.TCDB=phangorn::NJ(as.matrix(tcdb.dist[sel.genome,sel.genome]))

setdiff(nj.AAI$tip.label,sco.94.concat.tree$tip.label)


# INGROUP
lsTrees =
  list('rib.94.concat'=rib.94.concat.tree.rooted,
       'sco.94.concat'=sco.94.concat.tree.rooted,
       'astral.94'=astral.94.tree.rooted,
       'nj.ANIb' = phytools::midpoint.root(nj.ANIb),
       'nj.AAI' = phytools::midpoint.root(nj.AAI),#,
       'nj.kmer'=phytools::midpoint.root(nj.kmer),
       'nj.Man' = phytools::midpoint.root(nj.Man),#,
       'nj.GGR' = phytools::midpoint.root(nj.GGR),
       'nj.OF'= phytools::midpoint.root(nj.OF),
       'nj.KO' = phytools::midpoint.root(nj.KO),
       'nj.CAZY' = phytools::midpoint.root(nj.CAZY),
       'nj.TCDB' = phytools::midpoint.root(nj.TCDB)

  )

#WITH OUTGROUP
lsTrees =
  list('concat.rps'=concat.rps,
  'concat.rps.nuc'=concat.rps.nuc,
  'concat.sco'=concat.sco,
  'concat.sco.nuc'=concat.sco.nuc,
  'astral.106.tree'=astral.106.tree,
  'nj.ANIb' = nj.ANIb,
  'nj.AAI' = nj.AAI,#,
  'nj.Man' = nj.Man,#,
  #'nj.Jac' = nj.Jac,
  'nj.GGR' = nj.GGR,
  'nj.OF'= nj.OF
  )

#Snipen
#w <- micropan::geneWeights(OG.pan[SelectedTips,],type="shell")
#pan.tree <- micropan::panTree(OG.pan[SelectedTips,], scale=0.1, weights=w)

#dim()


names_of_trees = names(lsTrees)
class(lsTrees) = 'multiPhylo'

#not normalised
#RFoulds <- phytools::multiRF(lsTrees)
#colnames(RFoulds) = names_of_trees
#rownames(RFoulds) = names_of_trees

#normalised
#The normalized Robinson-Foulds distance is derived by dividing d(T_1, T_2) by the maximal possible distance i(T_1) + i(T_2). If both trees are unrooted and binary this value is 2n-6
RFoulds <- as.matrix(phangorn::RF.dist(lsTrees, normalize = TRUE, check.labels = TRUE))
colnames(RFoulds) = names_of_trees
rownames(RFoulds) = names_of_trees

#heatmap.2(RFoulds,trace='none')
p.fouldsHeat <- factoextra::fviz_dist(as.dist(RFoulds), gradient = list(low = "#FFFFFF", high = "#00626a"))+theme(aspect.ratio = 1)
p.fouldsHeat
#gplots::heatmap.2(RFoulds,trace='none')

MDS = cmdscale(RFoulds, eig = TRUE, x.ret=TRUE)
MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
MDSdata = data.frame(MDS$points)
MDSdata$Name = rownames(MDSdata)
MDSdata$type = as.numeric(grepl('concat',MDSdata$Name))
MDSdata$type = ifelse(MDSdata$type == 1, 'alignment','distance')
colnames(MDSdata) = c("x",'y','Name','type')

p <- ggplot(MDSdata,aes(x=x,y=y,fill=type,label=Name)) +
  theme_classic() +
  xlab(paste('PCoA1 (', MDS.var.perc[1], '%',')',sep='')) +
  ylab(paste('PCoA2 (', MDS.var.perc[2], '%',')',sep='')) +
  labs(fill = "Genus") +
  ggtitle(label='Robinsons Foulds')+
  geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_point(shape = 21, size = 2) +
  geom_text_repel()+
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position='none')
p


grDevices::pdf(paste0(outdir,'tree_of_trees','.pdf'), height = 5, width=5)

par(mfrow = c(1, 1))
p.NJfounld<- plot(phangorn::NJ(RFoulds),
     type = "unrooted",
     no.margin = TRUE
)
dev.off()




ggsave('~/DATA/MarbGenomics/Graphs/RobinsonFoulds.pdf',plot=p, width = 3, height = 3,unit='in')
#ggsave('~/DATA/MarbGenomics/Graphs/RobinsonFoulds_NJ.pdf',plot=p.NJfounld, width = 4, height = 4,unit='in')
ggsave('~/DATA/MarbGenomics/Graphs/RobinsonFoulds_heatmap.pdf',plot=p.fouldsHeat, width = 4, height = 4,unit='in')


colors = colorRampPalette(brewer.pal(8, "Set2"))(20)


#====================
# tips are longer in gene content tree
p.tree <- ggtree(sco.94.concat.tree.rooted)+geom_treescale()
p.gene <-ggtree(nj.GGR)+geom_treescale()
p.ov <- grid.arrange(p.tree,p.gene,nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/tipdistances.pdf',plot=p.ov, width = 5, height = 4,unit='in')


tree1 <- sco.94.concat.tree.rooted
tree2 <- nj.GGR


p.list <- list()
for(phgr in phylogroup.names){
  ingroup <- phylogroup.offspring[phylogroup.offspring$phylogroup == phgr,'genome'] %>% as.character()
  sub.p.tree <- ape::keep.tip(tree1,ingroup)
  sub.p.gene <- ape::keep.tip(tree2,ingroup)

  p.tree <- ggtree(sub.p.tree)+geom_treescale()+ggtitle(paste0(phgr,'Sequence similarity'))
  p.gene <-ggtree(sub.p.gene)+geom_treescale()+ggtitle('Gene content')

  p.ov <- grid.arrange(p.tree,p.gene,nrow=1)
  p.list[[phgr]] <- p.ov
}

p.ov <- grid.arrange(grobs=p.list)
ggsave(paste0('~/DATA/MarbGenomics/Graphs/tipdistances_per_P','.pdf'),plot=p.ov, width = 10, height = 10,unit='in')


#====================================================================================#
# recombination driven delay in molecular clock?#
#====================================================================================#

cop.tree <- cophenetic(sco.94.concat.tree)
cop.gene <- cophenetic(nj.Man)
genetreetype<-'nj.man'
dist.df.comp <- cbind('seq'=c(as.matrix(cop.tree)),
'gene'=c(as.matrix(cop.gene[rownames(cop.tree),colnames(cop.tree)])))%>% data.frame()

dist.df.comp %>% ggplot(aes(gene,seq))+geom_point() + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE, color=1)

log.model <-lm(log(seq) ~ gene, dist.df.comp)
exp.model <-lm(seq ~ exp(gene), dist.df.comp)

log.model.df <- data.frame(x = dist.df.comp$gene,
                           y = exp(fitted(log.model)))

ggplot(dist.df.comp, aes(x=gene, y=seq)) +
  geom_point() +
  geom_smooth(method="nls", method.args = list(start=c(a=1,b=1)), se=FALSE, formula= y~a*exp(b*x), colour='orange', linetype = 1) +
  #geom_line(data = log.model.df, aes(x, y, color = "Log Model"), size = 1, linetype = 2) +
  guides(color = guide_legend("Model Type"))+fdb_style()

ggplot(dist.df.comp, aes(x=seq, y=gene)) +
  geom_point() +
  geom_smooth(method="nls", method.args = list(start=c(a=1,b=1)), se=FALSE, formula= y~a*exp(b*x), colour='orange', linetype = 1) +
  #geom_line(data = log.model.df, aes(x, y, color = "Log Model"), size = 1, linetype = 2) +
  guides(color = guide_legend("Model Type"))+fdb_style()+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)

dist.df.comp %>% filter(seq<0.1)%>%
ggplot( aes(x=seq, y=gene))+
  geom_point() +
  geom_smooth(method="nls", method.args = list(start=c(a=1,b=1)), se=FALSE, formula= y~a*exp(b*x), colour='orange', linetype = 1) +
  #geom_line(data = log.model.df, aes(x, y, color = "Log Model"), size = 1, linetype = 2) +
  guides(color = guide_legend("Model Type"))+fdb_style()+
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)




#===== this is interestig
#

cop.tree <- cophenetic(nj.ANIb)
cop.tree <- 1-ANIb[sco.94.concat.tree$tip.label,sco.94.concat.tree$tip.label]
cop.tree <- cophenetic(sco.94.concat.tree)

cop.gene <- cophenetic(nj.GGR)
cop.gene <- GRR.sel

cop.gene <- 1-ANIb
genetreetype <- 'ANI'


cop.tree <- ANIb[sco.94.concat.tree$tip.label,sco.94.concat.tree$tip.label]
cop.gene <- GRR.sel


genetreetype <- 'nj.ggr'
seqtreetype <- 'COPH'

groupdf <- sorted.ANI.cliques %>% select(genome,'sc0.95')
colnames(groupdf) = c('genome','group')
groupdf <- groupdf %>% mutate(group=paste0('cl',group))
groupdf <- groupdf[groupdf$group!='clNA',]

groupdf <- genome.tbl %>% select(genome,phylogroup) %>%
  filter(!is.na(phylogroup)) %>% mutate(group=phylogroup) %>% select(-phylogroup)

df_out <- NULL
for(grp in unique(groupdf$group)){
  if(!is.na(grp)){
    message(grp)
    gnms <- groupdf[groupdf$group==grp,'genome']
    gnms <- gnms[!is.na(gnms)]

    dist.df.comp <- cbind('seq'=c(as.matrix(cop.tree[gnms,gnms])),
                          'gene'=c(as.matrix(cop.gene[gnms,gnms])))%>%
      data.frame() %>%
      mutate(group=grp) %>%
      mutate(level='intra') %>%
      mutate(n=length(gnms))

    dist.df.inter <- cbind('seq'=c(as.matrix(cop.tree[gnms,setdiff(groupdf$genome,gnms)])),
                          'gene'=c(as.matrix(cop.gene[gnms,setdiff(groupdf$genome,gnms)])))%>%
      data.frame() %>%
      mutate(group=grp) %>%
      mutate(level='inter') %>%
      mutate(n=length(setdiff(groupdf$genome,gnms)))

    df_out <- rbind(df_out, dist.df.comp, dist.df.inter)
    #dist.df.comp %>% ggplot(aes(gene,seq))+geom_point() + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE, color=1)

  }
}






formula <- y ~ poly(x, 2, raw = TRUE)
df_out %>% filter(n>2) %>%
  ggplot(aes(gene, seq, color = as.factor(group))) +
  geom_point() +
  stat_smooth(aes(color=as.factor(group)),method = "lm", formula = formula)+fdb_style()+ facet_wrap(~group,nrow=1)

df_out %>% filter(n>2) %>%
  ggplot(aes(seq,gene)) +
  geom_point( aes(color = as.factor(group))) +
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()

df_out %>% filter(n>2) %>% filter(level=='inter') %>%
  ggplot(aes(seq, gene, color = as.factor(group))) +
  geom_point( aes(color = as.factor(group))) +
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()+
  facet_wrap(~group,nrow=1) +
  scale_color_manual(values=phylogroup.colors)+
  geom_smooth(linetype='dashed',alpha=.4,aes(fill=group))+
  scale_fill_manual(values=phylogroup.colors)

p.in1 <- df_out %>% filter(n>2) %>%
  mutate(group=factor(group,levels=c('P1','P3','P4','P5','P7','P11')))%>%
  filter(level=='intra') %>%
  ggplot(aes(seq, gene, color = as.factor(group))) +
  geom_point( aes(color = as.factor(group))) +
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()+
  facet_wrap(~group,nrow=1) +
  scale_color_manual(values=(phylogroup.colors))+
  #scale_color_manual(values=unname(phylogroup.colors))+
  #geom_smooth(linetype='dashed',alpha=.4,aes(fill=group))+
  #scale_fill_manual(values=unname(phylogroup.colors))
  scale_fill_manual(values=(phylogroup.colors))+
  xlab('sequence tree distance')+  ylab('Gene content tree distance')
  ggsave(paste0('~/DATA/MarbGenomics/Graphs/gene_content_vs_sequence_tree_distances2_',seqtreetype,'_',genetreetype,'.pdf'),plot=p.in1, width = 10, height = 3,unit='in')
p.in1 <- df_out %>% filter(n>2) %>%
  mutate(group=factor(group,levels=c('P1','P3','P4','P5','P7','P11')))%>%
  filter(level=='intra') %>%
  ggplot(aes(seq, gene, color = as.factor(group))) +
  geom_point( aes(color = as.factor(group))) +
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()+
  facet_wrap(~group,nrow=2,scales='free') +
  scale_color_manual(values=(phylogroup.colors))+
  #scale_color_manual(values=unname(phylogroup.colors))+
  #geom_smooth(linetype='dashed',alpha=.4,aes(fill=group))+
  #scale_fill_manual(values=unname(phylogroup.colors))
  scale_fill_manual(values=(phylogroup.colors))+
  xlab('sequence tree distance')+  ylab('Gene content tree distance')+theme(legend.position = 'none')

ggsave(paste0('~/DATA/MarbGenomics/Graphs/gene_content_vs_sequence_tree_distances3_',seqtreetype,'_',genetreetype,'.pdf'),plot=p.in1, width = 10, height = 4,unit='in')


p.inter2 <- df_out %>% filter(n>2) %>%
  mutate(group=factor(group,levels=c('P1','P3','P4','P5','P7','P11')))%>%
  filter(level=='inter') %>%
  filter(!is.na(group)) %>%
  ggplot(aes(seq, gene)) +
  geom_point( aes(color = as.factor(group))) +
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F,color='grey30')+fdb_style()+
  geom_smooth(method = "lm", color='grey30')+fdb_style()+
  facet_wrap(~group,nrow=2,scales='free') +
  scale_color_manual(values=(phylogroup.colors))+
  scale_fill_manual(values=(phylogroup.colors))+
  xlab('sequence tree distance')+  ylab('Gene content tree distance')+
  theme(legend.position = 'none')

ggsave(paste0('~/DATA/MarbGenomics/Graphs/gene_content_vs_sequence_INTER_tree_distances2_',seqtreetype,'_',genetreetype,'.pdf'),plot=p.inter2, width = 10, height = 4,unit='in')




p.inter2 <- df_out %>% filter(n>2) %>%
  ggplot(aes(seq, gene, color = as.factor(group),group=level)) +
  geom_point( aes(color = as.factor(group),alpha=level)) +
  geom_smooth(method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()+
  #facet_wrap(~group,nrow=1) +
  #geom_smooth(linetype='dashed',alpha=.4,aes(fill=group))+
  scale_color_manual(values=unname(phylogroup.colors))+
  scale_fill_manual(values=unname(phylogroup.colors))

ggsave(paste0('~/DATA/MarbGenomics/Graphs/gene_content_vs_sequence_both_',seqtreetype,'_',genetreetype,'.pdf'),plot=p.inter2, width = 10, height = 4,unit='in')




#====================≠#

df_out %>% filter(n>2) %>%
  mutate(group=factor(group,levels=c('P1','P3','P4','P5','P7','P11')))%>%
  #filter(level=='intra') %>%
  ggplot(aes(seq, gene)) +
  geom_point(data= . %>% filter(level=='intra'), shape=21,size=0.7, aes(color = as.factor(group))) +
  geom_point(data= . %>% filter(level=='inter'), shape=21,size=0.7, color='grey') +
  geom_smooth(data=. %>% filter(level=='intra') %>% filter(group %in% c('P1','P4','P5','P7')), aes(color = as.factor(group)),method = "nls", formula = y~SSasymp(x,Assym,R0,lrc),se=F)+fdb_style()+
  #facet_wrap(~group,nrow=1) +
  scale_color_manual(values=(phylogroup.colors))+
  #scale_color_manual(values=unname(phylogroup.colors))+
  #geom_smooth(linetype='dashed',alpha=.4,aes(fill=group))+
  #scale_fill_manual(values=unname(phylogroup.colors))
  scale_fill_manual(values=(phylogroup.colors))+
  xlab('sequence tree distance')+  ylab('Gene content tree distance')



df_out %>% filter(n>2) %>%
  mutate(group=factor(group,levels=c('P1','P3','P4','P5','P7','P11')))%>%
  #filter(level=='intra') %>%
  ggplot(aes(seq, gene)) +
  geom_point(data= . %>% filter(level=='intra'), shape=21,size=0.7, aes(color = as.factor(group))) +
  geom_point(data= . %>% filter(level=='inter'), shape=21,size=0.7, color='grey') +
  geom_smooth(data=. %>% filter(level=='intra') %>% filter(group %in% c('P1','P4','P5','P7')), aes(color = as.factor(group)))+fdb_style()+
  #facet_wrap(~group,nrow=1) +
  scale_color_manual(values=(phylogroup.colors))+
  #scale_color_manual(values=unname(phylogroup.colors))+
  #geom_smooth(linetype='dashed',alpha=.4,aes(fill=group))+
  #scale_fill_manual(values=unname(phylogroup.colors))
  scale_fill_manual(values=(phylogroup.colors))+
  xlab('sequence tree distance')+  ylab('Gene content tree distance')










#==================
lsTrees[[1]]

pl.type <- "phylogram"
#h <- 12
h <-7
w <- 12
#h<-10
#w<-10

mtdt <- genome.tbl %>% select(genome, phylogroup) %>% data.frame()
rownames(mtdt) = mtdt$genome

grDevices::pdf(paste0(outdir,'phylogenies_',pl.type,'.pdf'), height = h, width =w)
par(mfrow = c(2, 6))
for(i in names(lsTrees)){
  plot.phylo(lsTrees[[i]],
             type = pl.type,
             show.tip.label = FALSE,
             tip.color = colors[as.character(df.chkm[lsTrees[[i]]$tip.label, 'group'])]
  )
  tiplabels(frame="none", pch=20, col=phylogroup.colors[mtdt[lsTrees[[i]]$tip.label, 'phylogroup']], cex=1.5, fg="transparent")
  title(i)
}
dev.off()


mtdt <- genome.tbl %>% select(genome, phylogroup) %>% data.frame()
rownames(mtdt) = mtdt$genome


p.list=list()
for(i in names(lsTrees)){
  p <- ggtree(lsTrees[[i]])
  p$data <- p$data %>% left_join(genome.tbl, by=c('label'='genome'))
  p<-p + geom_tippoint(aes(color=phylogroup)) +
    scale_color_manual(values=phylogroup.colors)+
    geom_treescale()+
    ggtitle(i)+
    theme(legend.position = 'none')
  p.list[[i]] <- p
}

p.comb <- gridExtra::grid.arrange(grobs=p.list,nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/NJ_midpoint_rooted_similarity.pdf',plot=p.comb, width = 8, height = 8,unit='in')

p.list=list()
for(i in names(lsTrees)){
  p <- ggtree(lsTrees[[i]],layout='circular',branch.length = 'none')
  p$data <- p$data %>% left_join(genome.tbl, by=c('label'='genome'))
  p<-p + geom_tippoint(aes(color=phylogroup)) +
    scale_color_manual(values=phylogroup.colors)+
    geom_treescale()+
    ggtitle(i)+
    theme(legend.position = 'none')
  p.list[[i]] <- p
}

p.comb <- gridExtra::grid.arrange(grobs=p.list,nrow=2)
ggsave('~/DATA/MarbGenomics/Graphs/NJ_midpoint_rooted_similarity_circular2.pdf',plot=p.comb, width = 12, height = 4,unit='in')



#--------------------------------------------#

COPHS <- ape::cophenetic.phylo(concat.sco.nuc)
COPHS

distMan <- distMan[rownames(COPHS),rownames(COPHS)]
#pan.ggr <- pan.ggr[rownames(COPHS),rownames(COPHS)]
ANIb.3 <- ANIb[rownames(COPHS),rownames(COPHS)]
AAI.3 <- AAI[rownames(COPHS),rownames(COPHS)]

dstANI = 1-as.matrix(ANIb.3)
dstAAI = 1-(as.matrix(AAI.3)/100)

dstANI
ingroup <- rownames(dstANI)[grepl('Marinobacter',rownames(dstANI))]

#Quantification of covariation between distance matrices using mantel test
# mantel test statistic (the Pearson correlation between distances)

mantal.1 = vegan::mantel(dstANI, dMan,method='pearson',permutations=999)
mantal.1 = vegan::mantel(COPH, dstANI,method='pearson',permutations=999)

skaawk <- data.frame(cbind('cophenetic'=c(COPHS),'ANIb'=c(as.matrix(dstANI))))
                          # ,'OGman'=c(as.dist(dMan))))
cor.test(skaawk$cophenetic, skaawk$ANIb, method = "pearson")


ade4::mantel.rtest(as.dist(COPH), as.dist(dstANI),nrepet=999)
ade4::mantel.rtest(as.dist(COPH), as.dist(dstAAI),nrepet=999)
ade4::mantel.rtest(as.dist(COPH), as.dist(dMan),nrepet=999)
ade4::mantel.rtest(as.dist(dstANI), as.dist(dMan),nrepet=999)

ade4::mantel.rtest(as.dist(COPH), as.dist(dstANI),nrepet=999)

ade4::mantel.rtest(as.dist(COPH[ingroup,ingroup]), as.dist(dstANI[ingroup,ingroup]),nrepet=999)
ade4::mantel.rtest(as.dist(COPH[ingroup,ingroup]), as.dist(dstAAI[ingroup,ingroup]),nrepet=999)
ade4::mantel.rtest(as.dist(dstANI[ingroup,ingroup]), as.dist(dstAAI[ingroup,ingroup]),nrepet=999)




mantal.1 = vegan::mantel(COPH, distMan,method='pearson',permutations=999)
mantal.1 = vegan::mantel(dstANI, distMan[SelectedTips,SelectedTips],method='pearson',permutations=999)
mantal.1 = vegan::mantel(COPH, dstANI,method='pearson',permutations=999)
mantel.2 = vegan::mantel(ANIb[rownames(geo.dist),rownames(geo.dist)], geo.dist,method='pearson',permutations=999 )
mantel.2 = vegan::mantel(geo.dist,ANIb[rownames(geo.dist),rownames(geo.dist)],method='pearson',permutations=999 )


mantel.correlog(dstANI, dMan)
mantel.correlogram <- vegan::mantel.correlog(geo.dist,ANIb[rownames(geo.dist),rownames(geo.dist)])
plot(mantel.correlogram)


#Do mantel tests!

#Do pearson correlation analysis
df.distances = data.frame(cbind('cophenetic'=c(COPHS),'ANIb'=c(as.matrix(ANIb*100)),'AAI'=c(as.matrix(AAI)),'Manhattan'=c(distMan),'GRR'=c(pan.ggr)))
cop.vs.all = melt(df.distances,id='cophenetic')


p.coph.v.all <- ggplot(cop.vs.all,aes(x= cophenetic,y= value)) +
  geom_hex(bins = 35, colour = NA) +
  viridis::scale_fill_viridis(option="magma",trans = 'sqrt', name = 'Frequency')+
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
  facet_wrap(~variable,nrow=1,scales='free_y',strip.position = "left",
             labeller = as_labeller(c(ANIb='ANIb', AAI='AAI',
                                      Manhattan = "Manhattan distance",
                                      Jaccard = "Jaccard distance",
                                      GRR='GRR')))+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 11))+ylab('')+xlab('Cophenetic distance')


p.coph.v.all
ggsave('~/DATA/MarbGenomics/Graphs/Cophenetic_vs_all.pdf',plot=p.coph.v.all, width = 13, height = 3,unit='in')





#--------------#
-#

df.distances



set.seed(27)
boots <- bootstraps(df.distances, times = 10, apparent = TRUE)
boots

#fit_nls_on_bootstrap <- function(split) {
#  nls(nlsFormula, analysis(split), start = list(k = 1, b = 0))
#}


fit_nls_on_bootstrap <- function(split) {
  nls(ANIb ~ SSasymp(cophenetic, copheneticf, cophenetic0, log_alpha), analysis(split))
}


boot_models <-
  boots %>%
  mutate(model = map(splits, fit_nls_on_bootstrap),
         coef_info = map(model, tidy))


fit_glm_on_bootstrap <- function(split) {
  glm(glmFormula, data= analysis(split))
}

boot_models <-
  boots %>%
  mutate(model = map(splits, fit_glm_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

#Confidence intervals
percentile_intervals <- int_pctl(boot_models, coef_info)
percentile_intervals

#f
p.btstrp.2 <- ggplot(boot_coefs, aes(estimate)) +
  geom_histogram(bins = 30,fill='grey20') +
  facet_wrap( ~ term, scales = "free") +
  geom_vline(aes(xintercept = .lower), data = percentile_intervals, col = "red") +
  geom_vline(aes(xintercept = .upper), data = percentile_intervals, col = "red") +
  scale_y_continuous(expand=c(0,0)) + fdb_style() + theme(panel.spacing=unit(0, 'cm'))

boot_aug <-
  boot_models %>%
  sample_n(200) %>%
  mutate(augmented = map(model, augment)) %>%
  unnest(augmented)

p.btstrp.3 <- ggplot(boot_aug, aes_string(x=x,y=y)) +
  geom_line(aes(y = .fitted, group = id), alpha = .2, col = "red") +
  geom_smooth(data=metadata,aes_string(x=x,y=y),method='lm',se=FALSE,color='darkgrey') +
  geom_point(data=metadata,aes_string(x=x,y=y),size=2, alpha=0.5,color='grey20') +
  fdb_style() + ggpubr::stat_cor()





#--------------

############-------- Not sure if I need this, but it is always nice to compare I guess...
#---- shrink dataset to include the rDNA stuff

concat.rps = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins/concat.treefile', clean_genome_names=TRUE)
concat.rps.nuc = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/RibosomalProteins_pal2nal/concat.treefile', clean_genome_names=TRUE)
concat.sco = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile', clean_genome_names=TRUE)
concat.sco.nuc = loadTree(inPath='~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_pal2nal_trimAl/concat.treefile', clean_genome_names=TRUE)
concat.rDNA = loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees//concat.treefile', clean_genome_names=TRUE)
rDNA16S = loadTree(inPath='~/DATA/MarbGenomics/rRNA.full.oneper/singletrees/16S.treefile', clean_genome_names=TRUE)




#Idialy find a solution for this, as this really is odd...
concat.rps$tip.label[concat.rps$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.rps.nuc$tip.label[concat.rps.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.sco$tip.label[concat.sco$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
concat.sco.nuc$tip.label[concat.sco.nuc$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"

concat.rDNA$tip.label[concat.rDNA$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"
rDNA16S$tip.label[rDNA16S$tip.label =='Marinobacter_pelagius_strain_CGMCC_1'] = "Marinobacter_pelagius_strain_CGMCC_1_6775"


SelectedTips = intersect(concat.rDNA$tip.label, nj.ANIb$tip.label)

concat.rps = keep.tip(concat.rps, SelectedTips)
concat.rps.nuc = keep.tip(concat.rps.nuc, SelectedTips)
concat.sco = keep.tip(concat.sco, SelectedTips)
concat.sco.nuc = keep.tip(concat.sco.nuc, SelectedTips)
nj.ANIb = keep.tip(nj.ANIb, SelectedTips)
nj.Man = keep.tip(nj.Man, SelectedTips)
nj.Jac = keep.tip(nj.Jac, SelectedTips)
nj.GGR = keep.tip(nj.GGR, SelectedTips)




lsTrees =
  list('concat.rps'=concat.rps,
       'concat.rps.nuc'=concat.rps.nuc,
       'concat.sco'=concat.sco,
       'concat.sco.nuc'=concat.sco.nuc,
       'concat.rDNA'=concat.rDNA,
       'rDNA16S'=rDNA16S,
       'nj.ANIb' = nj.ANIb,
       'nj.Man' = nj.Man,
       'nj.Jac' = nj.Jac
#       'nj.GGR' = nj.GGR
      )

names_of_trees = names(lsTrees)
class(lsTrees) = 'multiPhylo'

#not normalised
RFoulds <- phytools::multiRF(lsTrees)
colnames(RFoulds) = names_of_trees
rownames(RFoulds) = names_of_trees

#normalised
#The normalized Robinson-Foulds distance is derived by dividing d(T_1, T_2) by the maximal possible distance i(T_1) + i(T_2). If both trees are unrooted and binary this value is 2n-6
RFoulds <- as.matrix(phangorn::RF.dist(lsTrees, normalize = TRUE, check.labels = TRUE))


#(RFoulds,trace='none')
factoextra::fviz_dist(as.dist(RFoulds), gradient = list(low = "#FFFFFF", high = "#00626a"))+theme(aspect.ratio = 1)
#gplots::heatmap.2(RFoulds,trace='none')



MDS = cmdscale(RFoulds, eig = TRUE, x.ret=TRUE)
MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
MDSdata = data.frame(MDS$points)
MDSdata$Name = rownames(MDSdata)
MDSdata$type = as.numeric(!grepl('nj',MDSdata$Name))
MDSdata$type = ifelse(MDSdata$type == 1, 'alignment','distance')
colnames(MDSdata) = c("x",'y','Name','type')

p <- ggplot(MDSdata,aes(x=x,y=y,fill=type,label=Name)) +
  theme_classic() +
  xlab(paste('PCoA1 (', MDS.var.perc[1], '%',')',sep='')) +
  ylab(paste('PCoA2 (', MDS.var.perc[2], '%',')',sep='')) +
  labs(fill = "Genus") +
  ggtitle(label='Robinsons Foulds')+
  geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_point(shape = 21, size = 2) +
  geom_text_repel()+
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position='none')
p


ggtree(rDNA16S)+geom_tiplab(size=3)
ggtree(nj.ANIb)+geom_tiplab(size=3)
ggtree(concat.rps.nuc)+geom_tiplab(size=3)
ggtree(concat.sco.nuc)+geom_tiplab(size=3)
ggtree(concat.rDNA)+geom_tiplab(size=3)




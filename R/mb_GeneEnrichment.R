#------------------------#
#	mb-EnrichmentTests
#----------#


# CONSIDER GROUPING, AS THIS IS IMPORTANT

# ADD isoelectric Point estimates
# use frequencies of bla
# kmer stuff?


metadata %>% dplyr::count(group)

metadata %>% dplyr::count(group, sort = TRUE)

#can we reduce the number of groups? by including other

metadata_df <- metadata %>% dplyr::transmute(group = case_when(stringr::str_detect(group,"cl14")~"cl14", #detect if substring is ther, ~ call it
                                                 stringr::str_detect(group,"cl7")~"cl7", #
                                                 stringr::str_detect(group,"cl4")~"cl4", #
                                                 stringr::str_detect(group,"cl9")~"cl9", #

                                                 TRUE ~ "Other"), #if not, then call it other
                              BinId, Genome_size, GC, Coding_density, group, sourceClass)

metadata_df <- metadata_df %>%
  dplyr::select(BinId, Genome_size, GC, Coding_density, group,sourceClass) #sourceClass, lon lat, are NA, forms problem in SMOTE



#a mix of categorical and continious data

# Issue with
# smaller group inbalance!
#

#use bootstrap samples, evaluate testing, and training sets
# use this as an estimate of how good it works
# can account for the variablity in group sizes, some are small, too small?
library(tidymodels)
meta_boot <- rsample::bootstraps(metadata_df)


# ---- SMOTE
#made to work with class impalance
#library(themis)
#SMOTE, generates new examples of minority class using nearest neighbours of these cases

#for this to work we need dummy variables, all numeric (so turn the classes into numbers)
# we also need to scale all

metadata_rec <- recipe(group ~ ., data=metadata_df) %>%
  update_role(BinId, new_role = 'id' ) %>% #descide what it is...
  step_dummy(sourceClass) %>% #Convert this predictor to numeric
  step_zv(all_predictors()) %>% #remove all predictors with zero variance
  step_normalize(all_predictors()) %>% #%>% # centered and scaled
  themis::step_smote(group) #SMOTE IT UP -- cheat?!

metadata_prep <- prep(metadata_rec)
juiced = juice(metadata_prep)
juiced
#BLAM 84 COLUMNS!?
#INTERESTING ! gives an intereting
#Factors are now as 1 and 0 in tables!

rf_spec <- rand_forest(trees=1000) %>%
  set_mode("classification") %>%
  set_engine("ranger")

metadata_wf <- workflow() %>%
  add_recipe(metadata_rec) %>%
  add_model(rf_spec)

metadata_wf

metadata_results <- fit_resamples(metadata_wf,
              resamples = meta_boot,
              control = control_resamples(save_pred=TRUE, verbose=TRUE))


metadata_results


#multinomial problems, need different methods, multiclass preformance metrics, multiclass accuracy
# ROC curve, descide how to extent this to multinomial
metadata_results %>%
  collect_metrics()

metadata_results %>%
  collect_predictions()

#all bootstraps
metadata_results %>%
  collect_predictions() %>% conf_mat(group, .pred_class)

# split up per bootstrap
metadata_results %>%
  collect_predictions() %>% group_by(id) %>%
  ppv(group,.pred_class)

#get the estimate

metadata_results %>%
    collect_predictions() %>% group_by(id) %>%
  ppv(group,.pred_class) %>% ggplot(aes(.estimate))+geom_histogram()











#============
#===

#ordered_minimal_chord.df <- ordered_minimal_chord.df[ paste0('cl',1:nrow(ordered_minimal_chord.df)),]




#----------------------------------------------------------#


# Tree object
tree = tree

#pan table
pan = t(COGpan)
pan = OG.pan

#-- different tests run, to see if this works for all of the abundance tables (pan tables)
#pan = MCRpan
#for OFpan
#pan = t(pan)
#pan = data.frame(pan_modules)
#Metadata
metadata = df.chkm


#----------------------------
# ordering step

# based on those rownames in metadata
tree = drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(metadata))])
pan = pan[rownames(metadata),]
pan = pan[,which(colSums(pan)!=0)]

#----------------------------
# Filtering Pan table

#Filter singleton features (uninformative when looking for groups specific enrichments)
pan = pan[,which(colSums(pan)!=1)]
dim(pan)

#Filter out non-varying features (SD=0)
pan <-pan[,which(apply(pan,2,sd)>0)]
dim(pan)


#	pan_clean()
#=================================================================================#
#' Clean up a pan table
#'
#' @param pan
#' @param remove_singletons
#' @param SD
#'
#' @return
#' @export
#'
#' @examples
#' pan_clean(OG.pan, remove_singletons = TRUE, SD=2)

pan_clean <- function(pan, remove_singletons=TRUE, SD=NULL){
  pan = pan[,which(colSums(pan)!=0)]
  if(remove_singletons ==TRUE){
    pan = pan[,which(colSums(pan)!=1)]
  }
  if(!is.null(SD)){
    pan <-pan[,which(apply(pan,SD,sd)>0)]
    }
  return(pan)
}


#	panHypergeometric
#=================================================================================#
#' hypergeometic test functions
#'
#' @param pan
#' @param metadata
#' @param grp
#'
#' @return
#' @export
#'
#' @examples
#' grp='cl1'
#' hmg.tbl <- panHypergeometric(pan, metadata, grp)
#' selected.og <- hmg.tbl %>%
#'   dplyr::mutate_if(is.factor, as.character) %>%
#'   dplyr::filter(group == 'cl8') %>%
#'   dplyr::filter(fdr.p.value_enriched_binary < 0.01) %>% dplyr::select(feature.id) %>% dplyr::pull()

panHypergeometric <- function(pan, metadata, grp){
  pan = pan[metadata$genome,]
  pan = pan_clean(pan, remove_singletons = TRUE, SD=2)

	nfeatures.ingroup <- sum(pan[metadata$group == grp, ]) # total features in the selected group
	nfeatures.outgroup <- sum(pan[metadata$group != grp, ])	# total features in the non-selected genomes

	output <- NULL
	for(feature in 1:ncol(pan)){
		feature.id <- colnames(pan)[feature]
		metadata$feature <- pan[,feature]

    	# Binary version
    pval <- phyper(
    		q = nrow(subset(metadata, group == grp & feature > 0)) - 1,
    		m = sum(metadata$feature > 0),
    		n = nrow(subset(metadata,feature == 0)),
    		k = nrow(subset(metadata, group == grp)),
    		lower.tail=FALSE
    		)

		score <- -log10(pval)

		pval2 <- phyper(
    			q = sum(subset(metadata, group == grp)$feature) - 1,
    			m = nfeatures.ingroup,
    			n = nfeatures.outgroup,
    			k = sum(metadata$feature),
    			lower.tail=FALSE
    			)

		#Test for depletion binary version
		pval_dep <- phyper(
    			q = nrow(subset(metadata, group == grp & feature > 0)),
    			m = sum(metadata$feature > 0),
    			n = nrow(subset(metadata,feature == 0)),
    			k = nrow(subset(metadata, group == grp)),
    			lower.tail=TRUE
    			)

		score_dep <- -log10(pval_dep)

		#Test for depletion Raw counts version
		pval2_dep <- phyper(q = sum(subset(metadata, group == grp)$feature),
                    m = nfeatures.ingroup, n = nfeatures.outgroup, k = sum(metadata$feature),lower.tail=TRUE)
		score2_dep<- -log10(pval2_dep)

		res <- data.frame(
			feature.id = feature.id,
			score_enriched_binary = score,
			p.value_enriched_binary = pval,
			score_depletion_binary=score_dep,
			p.value_depletion_binary=pval_dep,
			score_enriched_rawcounts = -log10(pval2),
			p.value_enriched_rawcounts = pval2,
			score_depletion_rawcounts=score2_dep,
			p.value_depletion_rawcounts=pval2_dep
			)

	metadata$feature <- NULL
	output <- rbind(output,res)
  }

  output<-data.frame(
      feature.id=output$feature.id,
      group = grp,
      score_enriched_binary=output$score_enriched_binary,
      z.score_enriched_binary= (output$score_enriched_binary - mean(output$score_enriched_binary)) / sd(output$score_enriched_binary),
      p.value_enriched_binary=output$p.value_enriched_binary,
      fdr.p.value_enriched_binary = p.adjust(output$p.value_enriched_binary, method = 'fdr'),

      score_depletion_binary=output$score_depletion_binary,
      z.score_depletion_binary=(output$score_depletion_binary - mean(output$score_depletion_binary)) / sd(output$score_depletion_binary),
  		p.value_depletion_binary=output$p.value_depletion_binary,
      fdr.p.value_depletion_binary=p.adjust(output$p.value_depletion_binary, method = 'fdr'),

  		score_enriched_rawcounts=output$score_enriched_rawcounts,
      z.score_enriched_rawcounts=(output$score_enriched_rawcounts - mean(output$score_enriched_rawcounts)) / sd(output$score_enriched_rawcounts),
  		p.value_enriched_rawcounts=output$p.value_enriched_rawcounts,
      fdr.p.value_enriched_rawcounts =p.adjust(output$p.value_enriched_rawcounts, method = 'fdr'),

  		score_depletion_rawcounts=output$score_depletion_rawcounts,
  		z.score_depletion_rawcounts=(output$score_depletion_rawcounts - mean(output$score_depletion_rawcounts)) / sd(output$score_depletion_rawcounts),
  		p.value_depletion_rawcounts=output$p.value_depletion_rawcounts,
      fdr.value_depletion_rawcounts=p.adjust(output$p.value_depletion_rawcounts, method = 'fdr')
		)
  return(output)
}


#panHypergeometric.all
#=================================================================================#
#' for all groups, test group versus non-group
#'
#' @param pan
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
#' hmg.grp.tbl <- panHypergeometric.all(pan, metadata)

panHypergeometric.all <- function(pan, metadata){
  outList_hgm = list()
  for(grp in unique(as.character(metadata$group))){
    message(paste('analysing', grp))
    outList_hgm[[grp]] = panHypergeometric(pan[metadata$genome,], metadata, grp)
  }
  hmg.grp.tbl <- outList_hgm  %>%
    tibble::as_tibble_col() %>%
    tidyr::unnest(cols=c(value))

  return(hmg.grp.tbl)
}



#======
# using the bootrtrap logic


hmg.tbl <- panHypergeometric(bootstraps$pan[[1]], bootstraps$metadata[[1]], 'cl1')

hmg.tbl %>% pivot_longer(cols=contains('fdr.value'),names_to='analysis',values_to='fdr') %>%
  select(feature.id, analysis, fdr) %>% filter(fdr<.05) %>%
  group_by(analysis) %>%
  tally()






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
columns <- c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')


pfam.pan
pan.list <- list('OGpan'=OG.pan)
columns <- c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')

#columns <- c('sc0.85')

sorted.cliques = sorted.ANI.cliques
out.colors = out.colors


#recent development - now removing small groups
remove_small_groups=TRUE

for(pan.name in names(pan.list)){
  message('you are brave... better have time to wait!')
  message(paste0('=====  ',pan.name,'  ====='))

  pan = pan.list[[pan.name]]
  outList_hgm.all = list()
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

    if(remove_small_groups == TRUE){
      larger_groups <- metadata %>% group_by(group) %>% tally() %>% filter(n>2) %>% select(group) %>% pull()
      metadata<-metadata %>% filter(group %in% larger_groups)

    }
    #use panHypergeometric.all to run all groups versus non-groups
    hmg.grp.tbl <- panHypergeometric.all(pan, metadata)
    outList_hgm.all[[selcol]] = hmg.grp.tbl %>% mutate(level=selcol)
  }

  out.tbl <- outList_hgm.all  %>%
    tibble::as_tibble_col() %>%
    tidyr::unnest(cols=c(value))

  write.table(out.tbl , file = paste0(outdir,pan.name,remove_small_groups,'hgm_allgroups.tsv'),sep='\t')
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#



#============================================================================================
# ===== HERE SELECT A SPECIFIC OUTPUT FILE AND ANALYSE
# SAY WE ARE MOST INTERESTED in OGpan variation at ANI>0.8
# we set the pan.name and sel.col accordingly

#Set variables
names(panList)
pan.name = "pfam"
sel.col = 'sc0.85'
pan.name = "OGpan"
sel.col = 'sc0.85'

#Load data
ANI.hmg.grp.tbl <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'hgm_allgroups.tsv'))
ANI.hmg.grp.tbl <- ANI.hmg.grp.tbl %>% as_tibble()
ANI.hmg.grp.tbl$group <- factor(ANI.hmg.grp.tbl$group,levels=paste0('cl',1:length(unique(ANI.hmg.grp.tbl$group))))


#for example extract significatn genes, of specific comp
sig.genes <- ANI.hmg.grp.tbl %>%
  filter(level=='sc0.85') %>%
  filter(fdr.p.value_enriched_rawcounts<0.05) %>%
  filter(group=='cl8') %>% select(feature.id) %>%
  pull %>% as.character()
sig.genes



#Visualise distribution
p.hpg.all <- ANI.hmg.grp.tbl %>%
  dplyr::select(feature.id, group,level, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  ggplot2::ggplot(aes(group, -log10(FDR.p.value))) +
  ggplot2::geom_jitter(aes(fill=-log(FDR.p.value)),shape = 21, size = 2, width = 0.2) +
  facet_wrap(~analysis, ncol=1) +
  geom_hline(yintercept = -log10(0.05), size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  scale_fill_gradientn(colours =brewer.pal(9, 'BuPu'), guide = 'none')  +
  facet_wrap(~level, nrow=1,scales = 'free')+
  fdb_style(aspect.ratio = 0.25)

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",pan.name,'p.val_hgm_allgroups.pdf'),plot=p.hpg.all, width = 15, height = 15,unit='in')


#Visualise up down at specific sel.col level
global.df <- ANI.hmg.grp.tbl %>%
  dplyr::filter(level==sel.col)%>%
  dplyr::select(feature.id, group,level, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  filter(FDR.p.value<0.05)%>% group_by(group) %>% dplyr::count(analysis) %>%
  mutate(direction = ifelse(grepl('depletion',analysis),'-','+')) %>%
  mutate(type = ifelse(grepl('_binary',analysis),'binary','raw'))

p.hpg.bar <- global.df %>%  ggplot(aes(group,n))+
  geom_bar(data=global.df %>% dplyr::filter(direction == '+'),aes(group,n,fill=direction),alpha=0.8, stat="identity", position=position_dodge())+
  geom_bar(data=global.df %>% dplyr::filter(direction == '-'),aes(group,-n,fill=direction),alpha=0.8, stat="identity", position=position_dodge())+
  geom_text(data=global.df %>% dplyr::filter(direction == '+'), aes(label= n), vjust=1.6, color="white", position = position_dodge(0.9), size=2.5)+
  geom_text(data=global.df %>% dplyr::filter(direction == '-'), aes(group,-n,label= n), vjust=-1.6, color="white", position = position_dodge(0.9), size=2.5)+
  facet_wrap(~type, ncol=1)+
  scale_fill_manual(values=c('forestgreen','purple'))+fdb_style()+geom_hline(yintercept = 0)

ggplot2::ggsave(filename=paste0("~/DATA/MarbGenomics/Graphs/",pan.name,'_',sel.col,'hgm_FDR_0.05.pdf'),plot=p.hpg.bar, width = 4, height = 5,unit='in')


#=============

# mingle a bit a annotation table
allannot <- ANI.hmg.grp.tbl %>% filter(level=='sc0.8')%>%
  left_join(annotation.tbl,by=c('feature.id'='OG')) %>%
  select(feature.id,product,kfm_domain,cazy_domain,tcdb_domain) %>%
  unique() %>%
  filter(product!='hypothetical protein') %>%
  dplyr::group_by(feature.id) %>%
  dplyr::mutate(product2 = paste(unique(product), collapse = '; ')) %>%
  dplyr::mutate(kfm_domain2 = paste(unique(kfm_domain), collapse = '; ')) %>%
  dplyr::mutate(cazy_domain2 = paste(unique(cazy_domain), collapse = '; ')) %>%
  dplyr::mutate(tcdb_domain2 = paste(unique(tcdb_domain), collapse = '; ')) %>%
  select(feature.id,product2, kfm_domain2, cazy_domain2, tcdb_domain2) %>% unique()


#annotate the results table
ANI.hmg.grp.tbl %>%
  dplyr::filter(level==sel.col)%>%
  dplyr::select(feature.id, group,level, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  filter(FDR.p.value<0.05) %>%
  left_join(allannot,by=c('feature.id'='feature.id')) %>%
  select(feature.id,group,level,analysis,FDR.p.value,product2,kfm_domain2,cazy_domain2,tcdb_domain2) %>%
  unique() %>%
  filter(group=='cl8') %>%
  #filter(grepl('depletion',analysis))
  filter(grepl('enriched',analysis))


grps <- ANI.hmg.grp.tbl %>%
  dplyr::filter(level==sel.col)%>%
  select(group) %>% pull %>% as.character() %>% unique()

hmg.mcr.out <- NULL
for(grp in grps){
  for(updown in c('enriched_binary','enriched_raw','depletion_binary','depletion_raw')){
    KOs <- ANI.hmg.grp.tbl %>%
      dplyr::filter(level==sel.col)%>%
      dplyr::select(feature.id, group,level, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
      tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
      filter(FDR.p.value<0.05) %>%
      filter(grepl(updown,analysis)) %>%
      filter(group==grp) %>%
      left_join(annotation.tbl,by=c('feature.id'='OG')) %>%
      select(kfm_domain) %>%
      pull() %>%
      unique()

    mcr_kegg <- kegg_MCR(KOs, old_module_table)

    if(nrow(mcr_kegg %>% filter(MCR>50))>0){
      mcr_kegg <- mcr_kegg %>% filter(MCR>50) %>%
        left_join(m2n,by='module') %>%
        mutate(group=grp) %>%
        mutate(analysis=updown)

      hmg.mcr.out <- rbind(hmg.mcr.out, mcr_kegg)
    }
    }
}


hmg.mcr.out %>% ggplot(aes(analysis,DESCRIPTION))+geom_point() + facet_wrap(~group,nrow=1)

hmg.mcr.out %>%
  ggplot2::ggplot(aes(analysis,DESCRIPTION))+
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~group,nrow=1) +
  fdb_style(aspect.ratio=10) +
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))







ANI.hmg.grp.tbl %>% filter(level=='sc0.8')%>%
  #filter(group=='cl8') %>%
  filter(fdr.p.value_depletion_binary<.05)





ANI.hmg.grp.tbl %>% filter(level=='sc0.8')%>%
  filter(group=='cl1') %>%
  filter(fdr.p.value_enriched_rawcounts<.05) %>%
  left_join(annotation.tbl,by=c('feature.id'='OG')) %>%
  select(feature.id,product,kfm_domain,cazy_domain,tcdb_domain) %>%
  unique() %>%
  filter(product!='hypothetical protein') %>%
  dplyr::group_by(feature.id) %>%
  dplyr::summarise(product2 = paste(unique(product), collapse = '; ')) %>%
  filter(!grepl('transpor',product2))


ANI.hmg.grp.tbl %>% filter(level=='sc0.8')%>%
  filter(group=='cl2') %>%
  filter(fdr.p.value_enriched_rawcounts<.05) %>%
  left_join(annotation.tbl,by=c('feature.id'='OG')) %>%
  filter(cazy_domain !='NA') %>%
  select(feature.id,product,kfm_domain,cazy_domain,tcdb_domain) %>%
  dplyr::group_by(feature.id) %>%
  dplyr::summarise(product2 = paste(unique(product), collapse = '; '))


         slice(20:100)

data.frame()


feature.id



#===================================================
#from Random forest
rf.all <- importance.all %>%
  mutate(comb = paste0(feature.id,'_',group))

HGM_RF_COMB <- ANI.hmg.grp.tbl %>%
  filter(level=='sc0.85') %>%
  mutate(comb = paste0(feature.id,'_',group)) %>%
  left_join( rf.all,by='comb')



HGM_RF_COMB %>% ggplot(aes(fdr.p.value_enriched_binary,MeanDecreaseAccuracy)) +geom_point(shape=21)+fdb_style()

HGM_RF_COMB %>%
  ggplot(aes(z.score_depletion_binary,z.score_depletion_rawcounts)) +
  geom_point(aes(color=fdr.p.value_depletion_binary),shape=21)+
  fdb_style()


HGM_RF_COMB %>% filter(group.x=='cl8')




#=============
# Check at the sets of signficicantly identified Orthologs
# seems like, there is an overlap, some are depleted in one group, and enriched in another, so one picks them up twice

#Overall Diagram
p.cutoff <- 0.05
selcol <- 'sc0.85'


sig_enriched_binary <- ANI.hmg.grp.tbl %>% filter(fdr.p.value_enriched_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_enriched_raw <- ANI.hmg.grp.tbl %>% filter(fdr.p.value_enriched_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_depleted_binary <- ANI.hmg.grp.tbl %>% filter(fdr.p.value_depletion_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
sig_depleted_raw <- ANI.hmg.grp.tbl %>% filter(fdr.value_depletion_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()

cOutSig <- c(sig_enriched_binary,sig_enriched_raw,sig_depleted_binary,sig_depleted_raw) %>% unique()
length(cOutSig)


lsOutSig <- list('bin_enriched'=sig_enriched_binary,
                 'raw_enriched'=sig_enriched_raw,
                 'bin_depleted'=sig_depleted_binary,
                 'raw_depleted'=sig_depleted_raw)

gplots::venn(lsOutSig[1:2])
gplots::venn(lsOutSig)


grps <-ANI.hmg.grp.tbl  %>% filter(level==selcol) %>% select(group) %>% unique() %>% pull %>% as.character() %>% sort()


heatmap.2(as.matrix(pfam.pan[,cOutSig]),trace='none')



# PLOT ALL EULER DIAGRAMS
library(eulerr)



p.list <- list()
for(grp in grps){
  message(grp)
  #grp='cl4'
  #Group Specific venn diagrams
  sig_enriched_binary <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.p.value_enriched_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
  sig_enriched_raw <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.p.value_enriched_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
  sig_depleted_binary <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.p.value_depletion_binary <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()
  sig_depleted_raw <- ANI.hmg.grp.tbl %>% filter(group==grp) %>% filter(fdr.value_depletion_rawcounts <p.cutoff) %>% filter(level==selcol) %>% select(feature.id) %>% pull() %>% as.character() %>% unique()

  cOutSig <- c(sig_enriched_binary,sig_enriched_raw,sig_depleted_binary,sig_depleted_raw) %>% unique()
  length(cOutSig)

  if(length(cOutSig)>0){
    lsOutSig <- list('bin_enriched'=sig_enriched_binary,
                     'raw_enriched'=sig_enriched_raw,
                     'bin_depleted'=sig_depleted_binary,
                     'raw_depleted'=sig_depleted_raw)
   ##gplots::venn(lsOutSig[1:2],show.plot=FALSE)
    fit2 = eulerr::euler(lsOutSig)
    #dev.off()
    p.list[[grp]] <- plot(fit2,
         quantities = TRUE,
         labels = list(font = 4),
         fills = list(fill = c("cadetblue2", "cadetblue",'coral', "coral3")),
         main=grp)

  }
}

p.comb <- grid.arrange(grobs=p.list)
ggsave(paste0('~/DATA/MarbGenomics/Graphs/,',pan.name,',hgm_0.8_pergroup_overlap.pdf'),plot=p.comb, width = 8, height = 8,unit='in')






p.list <- plot.range_pcoa(pan = OG.pan[,cOutSig],
                          dist.method='Manhattan',
                          sorted.cliques=sorted.ANI.cliques,
                          columns=c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98'),
                          shrink=TRUE)

p.ani.pcoa <-grid.arrange(grobs = p.list, ncol = 6)

p.ani.pcoa







#==============================================================#

# # or manual
# clique2representative <- c('cl1'='Marinobacter_adhaerens_HP15',
#                           'cl2'='Marinobacter_psychrophilus_20041',
#                           'cl3'='Marinobacter_antarcticus_CGMCC1_10835',
#                           'cl4'='Marinobacter_piscensis_Abdou3 ',
#                           'cl5'='Marinobacter_hydrocarbonoclasticus_ATCC_49840',
#                           'cl6'='Marinobacter_pelagius_strain_CGMCC_1_6775',
#                           'cl7'='Marinobacter_lipolyticus_SM19',
#                           'cl8'='Marinobacter_algicola_DG893',
#                           'cl7'='Marinobacter_lipolyticus_SM19',
#                           'cl9'='Marinobacter_algicola_DG893',
#                           'cl8'='Marinobacter_algicola_DG893',
#                           'cl9'='Marinobacter_sp_YJ_S3_2')


#-------------------------------------------------------------------------------#
# THE ABOVE LONG COMPUTATION EXPORTS OUTPUT FILES!
# WE LOAD THEM AT THIS STAGE BACK IN, AND CROSS ANALYSE THEM!
#-------------------------------------------------------------------------------#

##---- ANALYSIS THAT RUNS OVER ALL PAN GENOME TABLES AND ALL ANI CLIQUE LEVELS
#   MIND THAT hmg.grp.tbl will hold a lot of information (loaded output from previous runs)
#   It is key to analyse them in a holistic manner, need to merge it with the RF output DATA

# THINK OF A WAY TO GENERATE NEAT OUTPUT TABLES
# CAN WE WRANGLE THEM SO THAT WE GET SIGNIFICANCE PER FEATURE, AND OUTPUT THAT AS SUPPLEMENTARY TABLES!

#### THE FIRSY LOOPS BELOW, ARE BUILT TO LOOK FOR GENE PROXIMITY IN ENRICHED GENES IN GENOME REPRESENTATIVE

panList = list('OGpan'=OG.pan, 'COGpan'=COG.pan, 'KOpan'=kfm.pan, 'CAZYpan'=cazy.pan, 'TCDBpan' = tcdb.pan)
outdir = '~/DATA/MarbGenomics/Graphs/'
pan.list=panList

output.tbl <- NULL
for(sel.col in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
    #Select relevant metadata
    mtd <- genome.tbl %>% left_join(sorted.ANI.cliques, by='genome') %>%
      left_join(out.colors, by='genome') %>%
      data.frame()
    mtd$grp <- mtd[,sel.col]
    mtd <- mtd %>% mutate(group=paste0('cl',grp))

    #set color scheme
    colors <- mtd %>% select(paste0('col.',sel.col), group) %>% unique()
    colors <-colors[order(colors$group),]
    grid.colors <- as.character(colors[,1])
    names(grid.colors) <- as.character(colors[,2])
    clique.colors <- grid.colors
    for(pan.name in names(pan.list)){
      # --Load data
      hmg.grp.tbl <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'hgm_allgroups.tsv'))
      hmg.grp.tbl <- hmg.grp.tbl %>% as_tibble()
      hmg.grp.tbl$group <- factor(hmg.grp.tbl$group,levels=paste0('cl',1:length(unique(hmg.grp.tbl$group))))

      # select representative genomes
      # GENERATE A LITTLE TABLE WITH CLIQUE REPRESENTATIVES
      representatives <- mtd %>% select(genome,group,TypeStrain,Completeness) %>%
        dplyr::arrange(genome,desc(Completeness),desc(TypeStrain),group_by=group) %>%
        group_by(group) %>% slice(1) %>% select(genome, group)

      clique2representative <- representatives$genome
      names(clique2representative)  <- representatives$group

      for(sel.group in names(clique2representative)[names(clique2representative)!='clNA']){
        sel.genome <- unname(clique2representative[sel.group])

        sel.features <- hmg.grp.tbl %>%
          filter(level==sel.col) %>%
          filter(group==sel.group) %>%
          filter(fdr.p.value_enriched_binary<0.05) %>%
          select(feature.id) %>%
          pull() %>%
          as.character()

        pan.name2col <- c("OG", "COG", "kfm_domain", "cazy_domain", "tcdb_domain")
        names(pan.name2col) <- c("OGpan", "COGpan", "KOpan", "CAZYpan", "TCDBpan")

        if(length(sel.features)!=0){
          hgm.bin.enriched <- annotation.tbl %>%
            dplyr::filter(genome == sel.genome) %>%
            dplyr::filter(!!as.symbol(unname(pan.name2col[pan.name])) %in% sel.features) %>%
            dplyr::select(locus_tag) %>%
            dplyr::pull()

          out.prox <- proximity_clusters(annotation.tbl, hgm.bin.enriched,distance=3000 )

          if(!is.null(out.prox)){
            #Select all >3 large proximity clusters
            prox.large <- out.prox %>% dplyr::group_by(prox.cluster) %>%
              select(prox.cluster) %>%
              dplyr::count() %>%
              dplyr::filter(n>3) %>%
              dplyr::select(prox.cluster) %>%
              dplyr::pull() %>% as.character()

            if(length(prox.large)>0){
              out.tbl <- NULL
              for(prox.clust in prox.large){
                message(prox.clust)
                cluster.loci <- out.prox %>%
                  filter(prox.cluster == prox.clust) %>%
                  select(locus_tag) %>%
                  pull()

                sel.annotations <- annotation.tbl %>%
                  filter(locus_tag %in% cluster.loci) %>%
                  mutate(prox.cluster = prox.clust)

                selected.summary <- sel.annotations %>%
                  select(locus_tag, product, OG, COG, gene, kfm_domain, cazy_domain,tcdb_domain)

                print(selected.summary)

                gaps.filled <- fillgaps(in.cluster=sel.annotations) %>%
                  mutate(hit = ifelse(locus_tag %in% sel.annotations$locus_tag, TRUE, FALSE)) %>%
                  mutate(prox.cluster = prox.clust) %>%
                  dplyr::mutate(direction=ifelse(strand=='+',1,-1))

                out.tbl <- rbind(out.tbl, gaps.filled)
              }
              out.tbl <- out.tbl %>% mutate(comb = paste0(prox.cluster," " ,seqnames))
              p.plt <- ggplot2::ggplot(out.tbl, aes(xmin = start, xmax = end, y = comb,forward=direction,label=Name,fill=prox.cluster,alpha=as.numeric(hit))) +
                gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
                gggenes::geom_gene_label(align = "left")+
                facet_wrap(~ prox.cluster, scales = "free", ncol = 1) +
                theme_genes()+
                theme(legend.position = "none")+scale_alpha(range=c(0.5,1))+ggtitle(paste0(sel.genome,' - ',sel.col, ' - ', sel.group ))

              ggplot2::ggsave(filename=paste0(outdir,'proximity_cluster_',pan.name,'_',sel.genome,'_',sel.col,'.pdf'),plot=p.plt, width = 10, height = 25,unit='in')

              out.tbl <- out.tbl %>% mutate(pan = pan.name) %>% mutate(sel.col=sel.col) %>% mutate(group=sel.group)
              output.tbl <- rbind(output.tbl, out.tbl)
            }
          }
       }
    }
  }
}

write.table(output.tbl, paste0(outdir,'proximity_cluster_ALL_COMBINED','.tsv'),sep='\t')
output.tbl <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED','.tsv'),sep='\t')

output.tbl %>% filter(sel.col=='sc0.85') %>%
  filter(group=='cl11') %>%
  filter(pan=='OGpan') %>%
  filter(product != 'hypothetical protein') %>%
  filter(prox.cluster=='cluster20')




#==================================================
#once this is run, we can explore the output.tbl

#-------------------------------------------------
# COMBINE SIGNIFICANT OG, COG, KOFAM, CAZY and TCDB to identify clusters, as oposed to seperate
for(sel.col2 in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  groups <- output.tbl %>% filter(hit==TRUE) %>%
    filter(sel.col==sel.col2) %>% select(group) %>% unique() %>% pull() %>% as.character()

  for(sel.group in groups){
  sel.loci <- output.tbl %>% filter(hit==TRUE) %>%
    filter(sel.col==sel.col2) %>%
    filter(group==sel.group) %>% select(locus_tag)%>%
    pull() %>% unique() %>%
    as.character()

  sel.genome <- output.tbl %>% filter(hit==TRUE) %>%
    filter(sel.col==sel.col2) %>%
    filter(group==sel.group) %>% select(genome)%>%
    unique() %>%
    pull() %>%
    as.character()

  if(length(sel.loci)!=0){
    out.prox <- proximity_clusters(annotation.tbl, sel.loci,distance=3000 )

    if(!is.null(out.prox)){
      #Select all >3 large proximity clusters
      prox.large <- out.prox %>% dplyr::group_by(prox.cluster) %>%
        select(prox.cluster) %>%
        dplyr::count() %>%
        dplyr::filter(n>3) %>%
        dplyr::select(prox.cluster) %>%
        dplyr::pull() %>% as.character()

      if(length(prox.large)>0){
        out.tbl <- NULL
        for(prox.clust in prox.large){
          message(prox.clust)
          cluster.loci <- out.prox %>%
            filter(prox.cluster == prox.clust) %>%
            select(locus_tag) %>%
            pull()

          sel.annotations <- annotation.tbl %>%
            filter(locus_tag %in% cluster.loci) %>%
            mutate(prox.cluster = prox.clust)

          selected.summary <- sel.annotations %>%
            select(locus_tag, product, OG, COG, gene, kfm_domain, cazy_domain,tcdb_domain)

          #print(selected.summary)

          gaps.filled <- fillgaps(in.cluster=sel.annotations) %>%
            mutate(hit = ifelse(locus_tag %in% sel.annotations$locus_tag, TRUE, FALSE)) %>%
            mutate(prox.cluster = prox.clust) %>%
            dplyr::mutate(direction=ifelse(strand=='+',1,-1))

          out.tbl <- rbind(out.tbl, gaps.filled)
        }
        out.tbl <- out.tbl %>% mutate(comb = paste0(prox.cluster," " ,seqnames))
        p.plt <- ggplot2::ggplot(out.tbl, aes(xmin = start, xmax = end, y = comb,forward=direction,label=Name,fill=prox.cluster,alpha=as.numeric(hit))) +
          gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
          gggenes::geom_gene_label(align = "left")+
          facet_wrap(~ prox.cluster, scales = "free", ncol = 1) +
          theme_genes()+
          theme(legend.position = "none")+scale_alpha(range=c(0.5,1))+ggtitle(paste0(sel.genome,' - ',sel.col2, ' - ', sel.group ))

        ggplot2::ggsave(filename=paste0(outdir,'proximity_cluster_ALL_COMBINED_',sel.genome,'_',sel.group,'_',sel.col2,'.pdf'),plot=p.plt, width = 10, height = 25,unit='in')
        write.table(out.tbl, paste0(outdir,'proximity_cluster_ALL_COMBINED_',sel.genome,'_',sel.group,'_',sel.col2,'.tsv'),sep='\t')

        #out.tbl <- out.tbl %>% mutate(pan = pan.name) %>% mutate(sel.col=sel.col) %>% mutate(group=sel.group)
        #output.tbl <- rbind(output.tbl, out.tbl)
        }
      }
    }
  }
}


  output.tbl2 <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED_Marinobacter_algicola_DG893_cl8_sc0.8','.tsv'),sep='\t')
  output.tbl2 %>% filter(prox.cluster=='cluster12')

  output.tbl2 <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED_Marinobacter_adhaerens_HP15_cl1_sc0.8','.tsv'),sep='\t')
  output.tbl2 %>% filter(prox.cluster=='cluster20')

  output.tbl2 <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED_Marinobacter_gelidimuriae_BF04_CF_4_cl2_sc0.8','.tsv'),sep='\t')
  output.tbl2 %>% filter(prox.cluster=='cluster4')

  output.tbl2 <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED_Marinobacter_daepoensis_DSM_16072_cl5_sc0.8','.tsv'),sep='\t')
  output.tbl2 %>% filter(prox.cluster=='cluster21')

  output.tbl2 <- read.delim(paste0(outdir,'proximity_cluster_ALL_COMBINED_Marinobacter_antarcticus_CGMCC1_10835_cl3_sc0.8','.tsv'),sep='\t')
  output.tbl2 %>% filter(prox.cluster=='cluster1')

  OG
  OG0003315
  annotation.tbl %>% filter(OG == 'OG0003315') %>% select(product)

  annotation.tbl %>% filter(OG == 'OG0003315') %>% select(genome) %>% pull


#===================================================================================#
#Sadly clique names are not comparable between clique levels
# alternatively we compare PAN HITS with one another

#run over all pans, and extract a neat output table for plotting and reporting

names(pan.list)
p.names <- names(pan.list)
p.names <- c('OGpan',"KOpan","COGpan","pfam","CAZY.pan","TCDBpan")
#

global.all.df <- NULL
for(pan.name in p.names){
  message('loading ',pan.name,',...')

  hmg.grp.tbl <- read.delim( paste0("~/DATA/MarbGenomics/",pan.name,'hgm_allgroups.tsv'))
  hmg.grp.tbl <- hmg.grp.tbl %>% as_tibble()
  hmg.grp.tbl$group <- factor(hmg.grp.tbl$group,levels=paste0('cl',1:length(unique(hmg.grp.tbl$group))))

  global.df <- hmg.grp.tbl %>%
    dplyr::select(feature.id, group,level, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
    tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
    filter(FDR.p.value<0.05)%>% group_by(group,level) %>% dplyr::count(analysis) %>%
    mutate(direction = ifelse(grepl('depletion',analysis),'-','+')) %>%
    mutate(type = ifelse(grepl('_binary',analysis),'binary','raw')) %>%
    mutate(pan=pan.name) %>% ungroup()

  global.all.df <- rbind(global.all.df,global.df)
}


#THIS IS A COOL TABLE
global.all.df


global.all.df %>% filter(level=='sc0.85') %>% ggplot(aes(group,n)) +
  geom_bar(aes(group,n,fill=direction),alpha=0.8, stat="identity", width=0.7,position=position_dodge(preserve = 'single')) + facet_grid(pan~type,scales='free_y') + fdb_style(aspect.ratio=0.5)



#NOW RUN OVER ALL LEVELS AND COMPARE
for(sel.col in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  message('visualising ',sel.col,',...')
  plot.df <- global.all.df %>% filter(level==sel.col)
  plot.df$pan <- factor(plot.df$pan, levels=names(pan.list))

  p.hpg.bar <- plot.df %>% ggplot(aes(group,n))+
    geom_bar(data=plot.df %>% dplyr::filter(direction == '+'),aes(group,n,fill=direction),alpha=0.8, stat="identity", position=position_dodge())+
    geom_bar(data=plot.df %>% dplyr::filter(direction == '-'),aes(group,-n,fill=direction),alpha=0.8, stat="identity", position=position_dodge())+
    geom_text(data=plot.df %>% dplyr::filter(direction == '+'), aes(label= n), vjust=1.6, color="white", position = position_dodge(0.9), size=2.5)+
    geom_text(data=plot.df %>% dplyr::filter(direction == '-'), aes(group,-n,label= n), vjust=-1.6, color="white", position = position_dodge(0.9), size=2.5)+
    facet_grid(pan~type,scale='free_y')+
    scale_fill_manual(values=c('forestgreen','purple'))+fdb_style(aspect.ratio=0.5)+
    geom_hline(yintercept = 0,size=.5)+
    ggtitle(paste0(sel.col))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90))

  ggplot2::ggsave(filename=paste0(outdir,'HPG_FDR0.05_',sel.col,'.pdf'),plot=p.hpg.bar, width = length(unique(as.character(plot.df$group)))*0.8, height = 8,unit='in')
}

#Using color scheme from the sc0.8 level
for(sel.col in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  message('visualising ',sel.col,',...')
  plot.df <- global.all.df %>% filter(level==sel.col)
  plot.df$pan <- factor(plot.df$pan, levels=names(pan.list))

  #------- select colors -------
  grid.colors <- cliquelevel2colors(sorted.ANI.cliques,out.colors, sel.col)

  p.hpg.bar <- plot.df %>% mutate(direction=ifelse(direction=='+',0,1)) %>% ggplot(aes(group,n))+
    geom_bar(data=plot.df %>% dplyr::filter(direction == '+'),aes(group,n,fill=group), stat="identity", position=position_dodge())+
    geom_bar(data=plot.df %>% dplyr::filter(direction == '-'),aes(group,-n,fill=group),stat="identity", position=position_dodge())+
    geom_text(data=plot.df %>% dplyr::filter(direction == '+'), aes(label= n), vjust=1.6, color="white", position = position_dodge(0.9), size=2.5)+
    geom_text(data=plot.df %>% dplyr::filter(direction == '-'), aes(group,-n,label= n), vjust=-1.6, color="white", position = position_dodge(0.9), size=2.5)+
    facet_grid(pan~type,scale='free_y')+
    scale_fill_manual(values=grid.colors)+
    fdb_style(aspect.ratio=0.5)+
    geom_hline(yintercept = 0,size=.5)+
    ggtitle(paste0(sel.col))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90))#+
    #scale_alpha(range=c(0,1))

  ggplot2::ggsave(filename=paste0(outdir,'HPG_FDR0.05_color_',sel.col,'.pdf'),plot=p.hpg.bar, width = length(unique(as.character(plot.df$group)))*0.8, height = 8,unit='in')
}


#============
# plot some score correlations

#Using color scheme from the sc0.8 level
for(sel.col in c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  message('visualising ',sel.col,',...')
  #------- select colors -------
  grid.colors <- cliquelevel2colors(sorted.ANI.cliques,out.colors, sel.col)

  p.enriched <- hmg.grp.tbl %>% filter(level==sel.col) %>% ggplot(aes(score_enriched_binary,score_enriched_rawcounts))+geom_point(shape=21,fill='black',alpha=.5)+facet_grid(~group,scales='free')+fdb_style()
  p.depleted <- hmg.grp.tbl %>% filter(level==sel.col) %>% ggplot(aes(score_depletion_binary,score_depletion_rawcounts))+geom_point(shape=21,fill='black',alpha=.5)+facet_grid(~group,scales='free')+fdb_style()
  p.combined <- gridExtra::grid.arrange(p.enriched, p.depleted,ncol=1)
  ggplot2::ggsave(filename=paste0(outdir,'HPG_cor_color_',sel.col,'.pdf'),plot=p.combined, width = length(unique(as.character(plot.df$group)))*2, height = 5,unit='in')
}


#===================================================================================#

# VENN Diagrams

#Overall plot (lumped, all linaege associated hits together)
enriched.binary <- hmg.grp.tbl %>% filter(level=='sc0.8') %>% filter(fdr.p.value_enriched_binary<0.05) %>% select(feature.id) %>% pull() %>% as.character()
enriched.raw <- hmg.grp.tbl %>% filter(level=='sc0.8') %>% filter(fdr.p.value_enriched_rawcounts<0.05) %>% select(feature.id) %>% pull() %>% as.character()
diff.list <- list('enriched.binary'=enriched.binary,'enriched.raw'=enriched.raw)
gplots::venn(diff.list)

#plots for each clique
enriched.binary <- hmg.grp.tbl %>%
  dplyr::filter(level=='sc0.8') %>%
  dplyr::filter(group=='cl8') %>%
  dplyr::filter(fdr.p.value_enriched_binary<0.05) %>%
  dplyr::select(feature.id) %>%
  dplyr::pull() %>% as.character()

enriched.raw <- hmg.grp.tbl %>%
  dplyr::filter(level=='sc0.8') %>%
  dplyr::filter(group=='cl8') %>%
  dplyr::filter(fdr.p.value_enriched_rawcounts<0.05) %>%
  dplyr::select(feature.id) %>%
  dplyr::pull() %>%
  as.character()

diff.list <- list('enriched.binary'=enriched.binary,'enriched.raw'=enriched.raw)
gplots::venn(diff.list)








p.hpg.bar <- global.df %>%  ggplot(aes(group,n))+
  geom_bar(data=global.df %>% dplyr::filter(direction == '+'),aes(group,n,fill=direction),alpha=0.8, stat="identity", position=position_dodge())+
  geom_bar(data=global.df %>% dplyr::filter(direction == '-'),aes(group,-n,fill=direction),alpha=0.8, stat="identity", position=position_dodge())+
  geom_text(data=global.df %>% dplyr::filter(direction == '+'), aes(label= n), vjust=1.6, color="white", position = position_dodge(0.9), size=2.5)+
  geom_text(data=global.df %>% dplyr::filter(direction == '-'), aes(group,-n,label= n), vjust=-1.6, color="white", position = position_dodge(0.9), size=2.5)+
  facet_grid(level~type)+
  scale_fill_manual(values=c('forestgreen','purple'))+fdb_style(aspect.ratio=0.5)+geom_hline(yintercept = 0)

ggplot2::ggsave(filename=paste0(outdir,'HPG_FDR0.05_',pan.name,'.pdf'),plot=p.hpg.bar, width = 10, height = 10,unit='in')




#==================





#===================================================================================#

#Select the one you want
# store as extracted goodies

#Selected.features
sel.group = 'cl8'

sel.features <- ANI.hmg.grp.tbl %>%
  filter(level==sel.col) %>%
  filter(group==sel.group) %>%
  filter(fdr.p.value_enriched_binary<0.05) %>%
  select(feature.id) %>%
  pull() %>%
  as.character()

#===== CURRENTLY LINK TO DETECT CLUSTERS















p.hpg.all <- ANI.hmg.grp.tbl %>%
  dplyr::select(feature.id, group,level, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  ggplot2::ggplot(aes(group, -log10(FDR.p.value))) +
  ggplot2::geom_jitter(aes(fill=group,alpha=-log(FDR.p.value)),shape = 21, size = 2, width = 0.2) +
  facet_wrap(~analysis, ncol=1) +
  geom_hline(yintercept = -log10(0.05), size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  #coord_flip() +
  #scale_fill_gradientn(colours =brewer.pal(9, 'BuPu'), guide = F)  +
  facet_wrap(~level, nrow=1,scales = 'free')+
  fdb_style(aspect.ratio = 0.25)


ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/hpg_allANI2.pdf',plot=p.hpg.all, width = 15, height = 15,unit='in')


#===================================================================================#
# ran over multi, ANI levels
#===================================================================================#




#### RUN CODE ####
#-----------------------------------------------------
#    Visualise multi-level results
#-----------------------------------------------------

MCR.pan <- read.delim(file = "~/DATA/MarbGenomics/MCR.pan.tsv",sep='\t')
#
rnmw <- rownames(MCR.pan)
MCR.pan <- sapply( MCR.pan, as.numeric )
rownames(MCR.pan) <- rnmw

#SET METADATA
metadata <-genome.tbl %>% left_join(sorted.ANI.cliques, by='genome') %>%
  left_join(out.colors, by='genome') %>% filter(!is.na(sc0.8)) %>%
  data.frame()
metadata$grp <- metadata[,'sc0.8']
metadata <- metadata %>% mutate(group=paste0('cl',grp)) %>% select(-grp)

#SET PAN
pan = MCR.pan[metadata$genome,]
#IF MCR, Remove those where min value < 50
pan = pan %>% data.frame() %>% select_if(~min(., na.rm = TRUE) <= 50)

#use panHypergeometric.all to run all groups versus non-groups
hmg.grp.tbl <- panHypergeometric.all(pan, metadata)

#set the order of groups
hmg.grp.tbl$group <- factor(hmg.grp.tbl$group, levels=c(paste0('cl',1:length(unique(hmg.grp.tbl$group)))))

#Use this tibble to visualise the FDR-p-value distribution for each group and analysis
hmg.grp.tbl %>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  ggplot2::ggplot(aes(group, -log10(FDR.p.value))) +
      ggplot2::geom_jitter(aes(fill=-log(FDR.p.value)),shape = 21, size = 2, width = 0.2) +
      facet_wrap(~analysis, ncol=1) +
      geom_hline(yintercept = -log10(0.05), size = 0.25, colour = '#bdbdbd') +
      geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
     #coord_flip() +
      scale_fill_gradientn(colours =brewer.pal(9, 'BuPu'), guide = F)  +
      fdb_style(aspect.ratio = 0.25)

#select only the FDR-pvalues, could be used as export table at some point
hmg.grp.tbl%>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts)

hmg.grp.tbl %>% group_by(group) %>% filter(p.value_enriched_binary < 0.01)


hmg.grp.tbl %>% ggplot2::ggplot(aes(-log10(fdr.p.value_enriched_binary), -log10(fdr.p.value_enriched_rawcounts), fill=group)) +
  ggplot2::geom_point(shape=21,size=2,alpha=0.5)+
  ggplot2::geom_hline(yintercept = -log10(0.05), size = 0.25, colour = '#bdbdbd') +
  ggplot2::geom_vline(xintercept = -log10(0.05), size = 0.25, colour = '#bdbdbd') +
  ggplot2::facet_wrap(~group) +
  fdb_style()


#lets zoom in
#which features are most strongly associated with groups
hmg.grp.tbl %>%
  mutate_if(is.factor, as.character) %>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  filter(analysis == 'fdr.p.value_enriched_binary') %>%
  group_by(group) %>%
  top_n(8, -log10(FDR.p.value)) %>%
  ungroup() %>%
  ggplot(aes(-log10(FDR.p.value), feature.id )) +
  geom_col()+
  facet_wrap(~group, scales="free_y")


hmg.grp.tbl %>%
  mutate_if(is.factor, as.character) %>%
  filter(group == 'cl8') %>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  filter(analysis == 'fdr.p.value_enriched_binary') %>%
  group_by(group) %>%
  top_n(15, -log10(FDR.p.value)) %>%
  ungroup() %>%
  #mutate_at(feature.id, funs(as.character)) %>%
  ggplot(aes(-log10(FDR.p.value), reorder(feature.id, -log10(FDR.p.value) ) )) +
    geom_col()+
    ylab('')+
    facet_wrap(~group, scales="free_y") + fdb_style(aspect.ratio=2.5)


sel.group = 'cl5'
selected_features <- hmg.grp.tbl %>%
  mutate_if(is.factor, as.character) %>%
  filter(group == sel.group) %>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  filter(analysis == 'fdr.p.value_enriched_binary') %>%
  group_by(group) %>%
  top_n(100, -log10(FDR.p.value)) %>% ungroup() %>% select(feature.id) %>% pull()

selected_features <- hmg.grp.tbl %>% filter(group==sel.group) %>%
  #filter(fdr.p.value_enriched_rawcounts < 0.01) %>%
  #fdr.value_depletion_rawcounts
  filter(fdr.p.value_enriched_rawcounts < 0.00000001) %>%
  top_n(200, -log10(fdr.p.value_enriched_rawcounts)) %>% select(feature.id) %>% pull() %>% as.character()

selected_features

##==========
m2n[m2n$module %in% selected_features,]



annotation.tbl %>% filter(OG %in% selected_features) %>% filter(genome =='Marinobacter_algicola_DG893') %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain) %>% data.frame()
annotation.tbl %>% filter(OG %in% selected_features) %>% filter(genome =='Marinobacter_adhaerens_HP15') %>% select(locus_tag, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain) %>% data.frame()


annotation.tbl %>% filter(OG %in% selected_features) %>% filter(genome =='Marinobacter_hydrocarbonoclasticus_VT8') %>% select(locus_tag, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain) %>% data.frame()


annotation.tbl %>% filter(cazy_domain %in% selected_features[2]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[3]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[4]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[5]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[6]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[7]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[8]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[9]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)

annotation.tbl %>% filter(OG %in% selected_features[10]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[11]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[12]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[13]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[14]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[15]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[16]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[17]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[18]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)

annotation.tbl %>% filter(OG %in% selected_features[19]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)
annotation.tbl %>% filter(OG %in% selected_features[20]) %>% select(genome, product, OG, COG, kfm_domain,cazy_domain, tcdb_domain)







#----- in case of MCR data, we can link them to m2n module conversion file
# Idially we have a feature.id to description for each of the pan tables, that would be handy, to track it as object in the functions
# so we can automatically annotate the below plot with different input files

#Below is hard coded example for Kegg module completion...
hmg.grp.tbl %>%
  mutate_if(is.factor, as.character) %>%
  filter(group == 'cl3') %>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  filter(analysis == 'fdr.p.value_enriched_binary') %>%
  group_by(group) %>%
  top_n(20, -log10(FDR.p.value)) %>% ungroup() %>% head(n=20) %>%
  left_join(m2n, by=c('feature.id'='module')) %>%
  mutate(fulldesc = paste(DESCRIPTION,feature.id,sep=' - ')) %>%
  ggplot(aes(-log10(FDR.p.value), reorder(fulldesc, -log10(FDR.p.value) ) )) +
    geom_col(aes(fill=TYPE))+
    ylab('')+
    facet_wrap(~group, scales="free_y") + fdb_style(aspect.ratio=2.5)




hmg.grp.tbl %>%
  mutate_if(is.factor, as.character) %>%
  filter(group == 'cl3') %>%
  dplyr::select(feature.id, group, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, fdr.p.value_depletion_binary, fdr.value_depletion_rawcounts) %>%
  tidyr::pivot_longer(cols=fdr.p.value_enriched_binary:fdr.value_depletion_rawcounts, names_to='analysis',values_to='FDR.p.value') %>%
  #filter(analysis == 'fdr.p.value_enriched_binary') %>%
  group_by(analysis) %>%
  top_n(20, -log10(FDR.p.value)) %>% ungroup() %>% #head(n=20) %>%
  left_join(m2n, by=c('feature.id'='module')) %>%
  mutate(fulldesc = paste(DESCRIPTION,feature.id,sep=' - ')) %>%
  ggplot(aes(-log10(FDR.p.value), reorder(fulldesc, -log10(FDR.p.value) ) )) +
  geom_col(aes(fill=TYPE))+
  ylab('')+
  facet_wrap(~analysis, scales="free_y",ncol=1) + fdb_style(aspect.ratio=1)



#--------------------------------------------




rawEnrList_hgm = list()
rawDeplList_hgm = list()
binEnrList_hgm = list()
binDeplList_hgm = list()

for(grp in names(outList_hgm)){
	print(paste('filtering results from', grp))
	seldf = outList_hgm[[grp]]
	rawEnrList_hgm[[grp]] = seldf[seldf $fdr.p.value_enriched_rawcounts < 0.05,]
	rawDeplList_hgm[[grp]] = seldf[seldf $fdr.p.value_depletion_rawcounts < 0.05,]
	binEnrList_hgm[[grp]] = seldf[seldf $fdr.p.value_enriched_binary < 0.05,]
	binDeplList_hgm[[grp]] = seldf[seldf $fdr.p.value_depletion_binary < 0.05,]
}


as.character(rawEnrList_hgm[['cl14']]$feature.id)
as.character(rawDeplList_hgm[['cl14']]$feature.id)
as.character(binEnrList_hgm[['cl14']]$feature.id)
as.character(binDeplList_hgm[['cl14']]$feature.id)

cog_desc[as.character(rawEnrList_hgm[['cl14']]$feature.id),]
cog_desc[as.character(rawDeplList_hgm[['cl14']]$feature.id),]
cog_desc[as.character(binEnrList_hgm[['cl14']]$feature.id),]
cog_desc[as.character(binDeplList_hgm[['cl14']]$feature.id),]


#look at algicola clade, inspect binary enriched
binEnrList_hgm[['cl14']]

clade_specific_enrichment <- binEnrList_hgm[['cl14']] %>%
    select(feature.id, score_enriched_binary, fdr.p.value_enriched_binary)

clade_specific_enrichment <- binEnrList_hgm[['cl14']] %>%
  select(feature.id, score_enriched_binary, fdr.p.value_enriched_binary)

binEnrList_hgm

lsOrthogroup["OG0004320", ]


scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)
heatmap.2(pan[, as.character(clade_specific_enrichment$feature.id)], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)



cnt.raw.EnrList_hgm = unlist(lapply(rawEnrList_hgm, function(x){nrow(x)}))
cnt.raw.DeplList_hgm = unlist(lapply(rawDeplList_hgm, function(x){nrow(x)}))
cnt.bin.EnrList_hgm = unlist(lapply(binEnrList_hgm, function(x){nrow(x)}))
cnt.bin.DeplList_hgm = unlist(lapply(binDeplList_hgm, function(x){nrow(x)}))



features.raw.Enr_hmg = lapply(rawEnrList_hgm, function(x){as.character(x$feature.id)})
features.raw.Depl_hmg = lapply(rawDeplList_hgm, function(x){as.character(x$feature.id)})
features.bin.Enr_hmg = lapply(binEnrList_hgm, function(x){as.character(x$feature.id)})
features.bin.Depl_hmg = lapply(binDeplList_hgm, function(x){as.character(x$feature.id)})



allgrp_raw.Enr = unname(unlist(features.raw.Enr_hmg))
allgrp_raw.Depl = unname(unlist(features.raw.Depl_hmg))

allgrp_bin.Enr = unname(unlist(features.bin.Enr_hmg))
allgrp_bin.Depl = unname(unlist(features.bin.Depl_hmg))


all.raw.grp_hgm = unique(c(allgrp_raw.Enr, allgrp_raw.Depl))
all.bin.grp_hgm = unique(c(allgrp_bin.Enr, allgrp_bin.Depl))

venn(list('raw'= all.raw.grp_hgm,'binary'= all.bin.grp_hgm))


allSig_hgm = unique(c(all.raw.grp_hgm, all.bin.grp_hgm))





scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)
heatmap.2(pan[, all.raw.grp_hgm], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)

scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)
heatmap.2(pan[, all.bin.grp_hgm], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)


cog_desc[]


lapply(features.raw.Enr_hmg, function(x){cog_desc[x,]})
lapply(features.raw.Depl_hmg, function(x){cog_desc[x,]})
lapply(features.bin.Enr_hmg, function(x){cog_desc[x,]})
lapply(features.bin.Depl_hmg, function(x){cog_desc[x,]})






output = outList_hgm[['cl8']]

ggplot(output,aes(x= -log10(fdr.p.value_enriched_binary), -log10(fdr.p.value_enriched_rawcounts)))+
	geom_point(shape = 21, size = 2,fill='grey',alpha=0.5) +
	theme_classic() +
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




rawEnriched = as.character(output[output$fdr.p.value_enriched_rawcounts < 0.05,]$feature.id)
rawDepleted = as.character(output[output$fdr.p.value_depletion_rawcounts < 0.05,]$feature.id)

#rawEnriched = as.character(output[output$p.value_enriched_rawcounts < 0.001,]$feature.id)
#rawDepleted = as.character(output[output$p.value_depletion_rawcounts < 0.001,]$feature.id)

combRaw = unique(c(rawEnriched, rawDepleted))

scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)
heatmap.2(pan[, combRaw], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2)


#annotation file

sigoutput = output[output$feature.id %in% combRaw,]

sigoutput$sigRawEnriched = as.character(sigoutput$feature.id) %in% rawEnriched
sigoutput$sigRawDepleted = as.character(sigoutput$feature.id) %in% rawDepleted

sigoutput  = cbind(sigoutput, cog_desc[as.character(sigoutput$feature.id),])

summarySigoutput = sigoutput[,c('COG', 'name', 'sigRawEnriched', 'sigRawDepleted')]

summarySigoutput[sigoutput$sigRawEnriched == TRUE,]
summarySigoutput[sigoutput$sigRawDepleted == TRUE,]












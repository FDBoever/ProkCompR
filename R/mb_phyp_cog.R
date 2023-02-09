
# TEST CALCULATE COG CATEGORY ENRICHMENT IN EACH OF THE ENRICHED GENE SETS

groups <-c("cl1","cl10","cl12","cl3","cl4","cl5","cl8")


sel.grp = 'cl12'
sel.COG.all <- cbind('func'=c('J','K','L','D','V','T','M','N','U','O','C','G','E','F','H','I','P','Q','R','S')) %>%
  data.frame()

ls.selected.f <- list()
for(sel.grp in groups){
  #extract top 250 OGs ranked according to variable importance of RF
  RF.250 <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
    dplyr::slice(1:250) %>%
    dplyr::select(feature.id.x) %>%
    dplyr::pull() %>%
    as.character()

  #combine HGM and RF data
  out.mat <- HGM_RF_COMB %>%
    dplyr::filter(group.x==sel.grp) %>%
    dplyr::select(feature.id.x, group.x, score_enriched_binary, fdr.p.value_enriched_binary, fdr.p.value_enriched_rawcounts, score_enriched_rawcounts, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
    dplyr::mutate(Enrich_bin = ifelse(fdr.p.value_enriched_binary<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(Enrich_raw = ifelse(fdr.p.value_enriched_rawcounts<0.05,TRUE,FALSE)) %>%
    dplyr::mutate(RFtop250 = ifelse(feature.id.x %in% RF.250,TRUE,FALSE))


  # select features that are enriched (HGM_FDR<.05) in either raw and binary counts
  selected.features <-
    out.mat %>%
    dplyr::filter(Enrich_raw ==TRUE | Enrich_bin ==TRUE)  %>%
    dplyr::select(feature.id.x) %>%
    dplyr::pull() %>%
    as.character()

  #incluging 250 RF
  selected.features <- unique(RF.250,selected.features)

  ls.selected.f[[sel.grp]] <- selected.features


  sel.COG <- annotation.tbl %>%
    dplyr::filter(type=='CDS') %>%
    dplyr::filter((OG %in% selected.features)) %>%
    select(OG, COG) %>%
    filter(!is.na(COG)) %>% unique() %>% arrange(desc(COG)) %>% group_by(OG) %>% slice(1) %>% ungroup() %>%
    dplyr::left_join(cog_func, by='COG') %>% group_by(func) %>% dplyr::count()%>%
    dplyr::mutate(n_acc=n) %>%
    dplyr::select(-n)

  #adding another
  sel.COG.all <- sel.COG.all %>% dplyr::left_join(sel.COG, by='func' ) %>% data.frame()

}



colnames(sel.COG.all) = c('func',groups)




COG.pan.df <- sel.COG.all



COG.pan.df$func <- as.character(COG.pan.df$func)

COG.pan.df[COG.pan.df$func=='',]
COG.pan.df[is.na(COG.pan.df$func),]

COG.pan.df[is.na(COG.pan.df$func),'func'] <- ''
COG.pan.df[COG.pan.df$func=='','func'] <- c('unassigned','unassigned')

COG.pan.df <- COG.pan.df %>% group_by(func)%>%
  summarise_each(funs(sum)) %>% data.frame()

#tidy up , this is marinobacter specific, needs a function
#COG.pan.df[is.na(COG.pan.df$n_singc),'n_singc'] <- c(0,0)
#COG.pan.df


#ok now we should be able to build contigency tables!

pan.class.df <- COG.pan.df
rownames(pan.class.df) <- pan.class.df$func
pan.class.df<- pan.class.df[,2:ncol(pan.class.df)]
pan.class.df <- t(pan.class.df)

rw.nm <- 'n_acc'

df <- pan.class.df

#filter the categories, to those that make sence (get rid of unassigned ones?)
#or low characterized
df <- df[,c('J','K','L','D','V','T','M','N','U','O','C','G','E','F','H','I','P','Q','R','S')]
df[is.na(df)] <- 0

outout <- NULL
for(rw.nm in rownames(df)){
  nfeatures.ingroup <- sum(df[rownames(df) == rw.nm, ]) # total features in the selected group
  nfeatures.outgroup <- sum(df[rownames(df) != rw.nm, ])	# total features in the non-selected genomes

  df.met <- data.frame(cbind('type'=rownames(df)))

  output <- NULL
  for(feature in 1:ncol(df)){
    feature.id <- colnames(df)[feature]
    df.met$feature <- df[,feature]

    # Binary version
    pval <- phyper(
      q = nrow(subset(df.met, type == rw.nm & feature > 0)) - 1,
      m = sum(df.met$feature > 0),
      n = nrow(subset(df.met,feature == 0)),
      k = nrow(subset(df.met, type == rw.nm)),
      lower.tail=FALSE
    )

    score <- -log10(pval)

    pval2 <- phyper(
      q = sum(subset(df.met, type == rw.nm)$feature) - 1,
      m = nfeatures.ingroup,
      n = nfeatures.outgroup,
      k = sum(df.met$feature),
      lower.tail=FALSE
    )

    #Test for depletion binary version
    pval_dep <- phyper(
      q = nrow(subset(df.met, type == rw.nm & feature > 0)),
      m = sum(metadata$feature > 0),
      n = nrow(subset(df.met,feature == 0)),
      k = nrow(subset(df.met, type == rw.nm)),
      lower.tail=TRUE
    )

    score_dep <- -log10(pval_dep)

    #Test for depletion Raw counts version
    pval2_dep <- phyper(q = sum(subset(df.met, type == rw.nm)$feature),
                        m = nfeatures.ingroup, n = nfeatures.outgroup, k = sum(df.met$feature),lower.tail=TRUE)
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
    group = rw.nm,
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

  outout <- rbind(outout,output)

}

outout %>% filter(fdr.p.value_enriched_rawcounts<0.05)

df_plot <- df %>% data.frame() %>% tibble::rownames_to_column(var='type') %>% melt(.,id="type") %>% mutate(comb = paste0(type,variable))
df_sign <- outout %>% mutate(comb=paste0(group,feature.id))
df_plot <- df_plot %>% left_join(df_sign,by='comb') %>% mutate(p.value = fdr.p.value_enriched_rawcounts) %>% select(type,variable,value,comb,p.value)

df_plot$stars <- cut(df_plot$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

#calculate proportions per group
#df1 <- df_plot %>%  filter(type=='cl1') %>%
#  mutate(per=value/sum(value))
#df2 <- df_plot %>%  filter(type=='cl3') %>%
#  mutate(per=value/sum(value))
#df3 <- df_plot %>%  filter(type=='cl8') %>%
#  mutate(per=value/sum(value))

df_plot_out <- c()
for(sel.grp in groups){
  df1 <- df_plot %>%  filter(type==sel.grp) %>%
    mutate(per=value/sum(value))
  df_plot_out <- rbind(df_plot_out,df1)
}

#merge together
#df_plot <- rbind(df1,df2,df3)
#df_plot_out$type <- factor(df_plot$type,levels=rev(groups))

#plot
p.hyperg <- df_plot_out %>%
  ggplot2::ggplot(aes(x=variable, y=type, fill=per))+
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "#2C7BB6", na.value = NA)+
  ggplot2::geom_text(aes(label=stars), color="black", size=5) +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::labs(y=NULL, x=NULL) +
  fdb_style(aspect.ratio =0.15)


p.hyperg

#ggsave('~/DATA/MarbGenomics/Graphs/allMB_COGCAT_hyperg_pan.pdf',plot=p.hyperg, width = 8, height = 2,unit='in')









annotation.quality.tbl %>%
  dplyr::filter((OG %in% ls.selected.f[['cl1']])) %>%
  select(OG, COG,product,kfm_domain,ko_annotation) %>%
  dplyr::left_join(cog_func, by='COG') %>%
  filter(func=='C') %>% unique() %>% group_by(OG) %>% slice(1) %>% arrange(OG)


annotation.quality.tbl %>%
  dplyr::filter((OG %in% ls.selected.f[['cl3']])) %>%
  select(OG, COG,product,kfm_domain,ko_annotation) %>%
  dplyr::left_join(cog_func, by='COG') %>%
  filter(func=='C') %>% unique() %>% group_by(OG) %>% slice(1) %>% arrange(OG)




annotation.tbl %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter((OG %in% ls.selected.f[['cl8']])) %>%
  select(OG, COG,product) %>% dplyr::left_join(cog_func, by='COG') %>% filter(func=='T') %>% unique()


annotation.quality.tbl %>%
#  dplyr::filter(type=='CDS') %>%
  dplyr::filter((OG %in% ls.selected.f[['cl8']])) %>%
  select(OG, COG,product,kfm_domain,ko_annotation,pfam_name) %>%
  dplyr::left_join(cog_func, by='COG') %>%
  filter(func=='C') %>% unique()


annotation.tbl %>%
  dplyr::filter(type=='CDS') %>%
  dplyr::filter((OG %in% ls.selected.f[['cl5']])) %>%
  select(OG, COG,product) %>%
  dplyr::left_join(cog_func, by='COG') %>%
  filter(func=='C') %>% unique()








#==

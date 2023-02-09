# DATA FROM MAGICLAMP HMM

magiclamp_output <- read.delim('~/DATA/MarbGenomics/summary.csv',sep=',')


magiclamp_output %>%
  mutate(ID =gsub('-','_',gsub('\\.','',gsub("\\.faa","",cell)))) %>%
  filter(ID=='Marinobacter_algicola_DG893') %>% filter(substrate=='nitrogen') %>% left_join(annotation.tbl,by=c('ORF'='locus_tag'))

magiclamp_output%>%
  mutate(ID =gsub('-','_',gsub('\\.','',gsub("\\.faa","",cell)))) %>%
  filter(ID=='Marinobacter_algicola_DG893') %>%
  left_join(annotation.tbl,by=c('ORF'='locus_tag')) %>% select(OG,ORF,HMM, reaction, Name, product)

annotation.tbl %>% filter(OG=='OG0002463')



magiclamp_output %>%  left_join(annotation.tbl,by=c('ORF'='locus_tag')) %>%
  select(OG,HMM,substrate) %>%
  filter(!is.na(OG)) %>%
  filter(substrate=='nitrogen')%>%
  unique()


magiclamp_output %>%   select(HMM,reaction,substrate) %>%
  unique()

#====
#Module extraction, of manually curated sets
og <- magiclamp_output %>%  left_join(annotation.tbl,by=c('ORF'='locus_tag')) %>%
  select(OG,HMM,substrate) %>%
  filter(!is.na(OG)) %>%
  filter(substrate=='nitrogen')%>%
  unique() %>% select(OG) %>% pull()

og <-  magiclamp_output %>%  left_join(annotation.tbl,by=c('ORF'='locus_tag')) %>%
  select(OG,HMM,reaction,substrate) %>%
  filter(!is.na(OG)) %>%
  filter(substrate=='sulfur')%>%
  unique() %>% select(OG) %>% pull()

og <-  magiclamp_output %>%  left_join(annotation.tbl,by=c('ORF'='locus_tag')) %>%
  select(OG,HMM,reaction,substrate) %>%
  filter(!is.na(OG)) %>%
  #filter(substrate=='sulfur')%>%
  unique() %>% select(OG) %>% pull()

og <- annotation.tbl %>%
  filter(grepl('flag',product,ignore.case=T)) %>%
  filter(genome=='Marinobacter_adhaerens_HP15') %>% data.frame() %>% select(OG) %>% unique()%>%
  #filter(substrate=='sulfur')%>%
  filter(!is.na(OG)) %>%
  unique() %>% select(OG) %>% pull()

og <-annotation.tbl %>%
  filter(grepl('transporter permease',product,ignore.case=T)) %>%
  data.frame() %>% select(OG) %>% unique()%>%
  filter(!is.na(OG)) %>%
    unique() %>% select(OG) %>% pull()

og <-annotation.tbl %>%
  filter(grepl('NADH',product,ignore.case=T)) %>%
  data.frame() %>% select(OG) %>% unique()%>%
  filter(!is.na(OG)) %>%
  unique() %>% select(OG) %>% pull()

og <-annotation.tbl %>%
  filter(grepl('pho|psiE',Name,ignore.case=T)) %>%
  data.frame() %>% select(OG) %>% unique()%>%
  filter(!is.na(OG)) %>%
  unique() %>% select(OG) %>% pull()

og <-annotation.tbl %>%
  filter(grepl('phyt',product,ignore.case=T)) %>%
  data.frame() %>% select(OG) %>% unique()%>%
  filter(!is.na(OG)) %>%
  unique() %>% select(OG) %>% pull()


#ATP_Synthase
m1 <- annotation.tbl %>%
  filter(grepl('ATP synthase',product,ignore.case=T)) %>%
  filter(genome=='Marinobacter_adhaerens_HP15') %>% data.frame() %>%
  select(OG,Name,product) %>%
  unique() %>%
  filter(OG != 'OG0002022') %>% mutate(module='ATP synthase')

m1 <- annotation.tbl %>% filter(OG %in% og) %>%
  select(OG,Name,product) %>%
  unique()%>% mutate(module='nitrogen')


#order tips
is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]
ordered_tips <- tree$tip.label[ordered_tips]

mglmp_df <- annotation.tbl %>% filter(OG %in% m1$OG) %>%
  select(OG, genome) %>% left_join(m1,by='OG') %>%
  group_by(OG, genome, Name, product, module) %>%
  tally() %>% ungroup() %>%
  filter(genome %in% ordered_tips) %>%
  mutate(genome=factor(genome,level=rev(ordered_tips)))

p.magiclamp <- mglmp_df %>% ggplot(aes(Name,genome))+geom_point(aes(color=n))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~module,scales='free_x',space = "free_x")
p.magiclamp

mglmp_df %>% ggplot(aes(Name,genome))+geom_point(aes(color=n,fill=n),shape=21)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~module,scales='free_x',space = "free_x")


p.magiclamp <- mglmp_df %>% ggplot(aes(OG,genome))+geom_point(aes(color=n))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~module,scales='free_x',space = "free_x")
p.magiclamp

ggsave('~/DATA/MarbGenomics/Graphs/magiclamp_fun_for_joy.pdf',p.magiclamp ,width=30,height = 15,units='in')



mglmp_df <- mglmp_df%>%mutate(id=genome) %>% select(id,Name,product,module,n)

ggtree(tree) %>% facet_plot( panel = 'Stacked Barplot', data = mglmp_df,
           geom = geom_point,
           mapping = aes(x = Name, fill = n,color=n),
           stat='identity' )




annotation.tbl %>% filter(OG=='OG0003067')%>% select(product,Name)
#2-haloacid dehalogenase hadL
annotation.tbl %>% filter(OG=='OG0003517')%>% select(product,Name)
annotation.tbl %>% filter(OG=='OG0002197')%>% select(product,Name)

OG0003517

magiclamp_output %>% filter(cell=='Marinobacter_hydrocarbonoclasticus_114E_o.faa')
magiclamp_output %>%
  filter(cell=='Marinobacter_algicola_DG893.faa') %>%
  filter(substrate=='carbon-monoxide')

magiclamp_output %>%
  mutate(ID =gsub('-','_',gsub('\\.','',gsub("\\.faa","",cell)))) %>%
  filter(ID=='Marinobacter_algicola_DG893') %>%
  left_join(annotation.tbl, by=c('ORF'='locus_tag'))%>%
  select(Name, product, ORF, HMM, reaction, substrate, evalue,bitscore) %>%
  as_tibble() %>%
  filter(substrate=='methane')


C1compounds
magiclamp_output %>%
  mutate(ID =gsub('-','_',gsub('\\.','',gsub("\\.faa","",cell)))) %>%
  filter(ID=='Marinobacter_subterrani_JG233') %>%
  left_join(annotation.tbl, by=c('ORF'='locus_tag'))%>%
  select(Name, product, ORF, HMM, reaction, substrate, evalue,bitscore) %>%
  as_tibble() %>%
  filter(substrate=='iron')

magiclamp_output %>%
  mutate(ID =gsub('-','_',gsub('\\.','',gsub("\\.faa","",cell)))) %>%
  filter(ID=='Marinobacter_subterrani_JG233') %>%
  left_join(annotation.tbl, by=c('ORF'='locus_tag'))%>%
  select(Name, product, ORF, HMM)



magiclamp_output %>% group_by(cell,HMM,reaction) %>% tally() %>% ggplot(aes(cell,HMM))+geom_point(aes(color=n))


mglmp_df <- magiclamp_output %>%
  group_by(cell,HMM,reaction,substrate) %>%
  tally() %>% mutate(ID=gsub("\\.faa","",cell)) %>% ungroup() %>% select(ID,HMM, reaction,substrate,n) %>%
  mutate(ID =gsub('-','_',gsub('\\.','',ID))) %>%
  filter(ID %in% tree$tip.label) %>% mutate(ID=factor(ID,level=rev(tree$tip.label)))

p.magiclamp <- mglmp_df %>% ggplot(aes(HMM,ID))+geom_point(aes(color=n))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~substrate,scales='free_x',space = "free_x")


ggsave('~/DATA/MarbGenomics/Graphs/biochemcycles_MB.pdf',p.magiclamp ,width=15,height = 15,units='in')


p.magiclamp2 <- mglmp_df %>% ggplot(aes(HMM,ID))+geom_point(aes(color=n))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(angle = 90, vjust = 0.5, hjust=0))+
  facet_grid(~substrate+reaction,scales='free_x',space = "free_x")

ggsave('~/DATA/MarbGenomics/Graphs/biochemcycles_MB_per_reaction.pdf',p.magiclamp2 ,width=15,height = 20,units='in')





tree

p3 <- ggtree(tree)
facet_plot(p3, 'Trait', data = mglmp_df, geom=geom_point, mapping=aes(x=n))

facet_plot(p2, panel = "Vehicles", data = mglmp_df, geom=geom_tile, mapping = aes(x=HMM, fill = n))




ggtree(sco.94.concat.tree.rooted, branch.length='none') +
      geom_tiplab() +
      ggtree::geom_facet(
        mapping = aes(x=as.numeric(HMM), fill = n),
        data = mglmp_df,
        geom = geom_tile,
        panel = 'Alignment')



ggtree(mergedGP) +
  geom_tippoint(aes(color=Phylum), size=1.5) +
  geom_facet(mapping = aes(x=val,group=label,
                           fill=Phylum),
             data = melt_simple,
             geom = geom_density_ridges,
             panel="Abundance",
             color='grey80', lwd=.3)





facet_plot(p3, 'heatmap', mglmp_df, geom_tile, aes(x=n, fill=n)) + theme_tree2(legend.position='right')

tr <- rtree(10)
dd = data.frame(id=tr$tip.label, value=abs(rnorm(10)))
p <- ggtree(tr)
facet_plot(p, 'Trait', data = dd, geom=geom_point, mapping=aes(x=value))



#================================
# FeGenie from hmm.tblout

fegenie_hmm_cutoffs <- read.delim('/Users/sa01fd/DATA/MarbGenomics/AnnotateMachine/fegenie_HMM-bitcutoffs.txt',header=TRUE,sep='\t')

out.all <- NULL
dirlist <- list.dirs("/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results")
dirlist <- dirlist[2:length(dirlist)]
for(i in dirlist){
  files<-list.files(i, include.dirs = TRUE)
  genome <- gsub('\\.','_',gsub('-','_',gsub(".faa-HMM","",basename(i))))
  message(genome)
  for(file in files){
    inPath = paste(i,file,sep='/')
    #inPath='/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Hahella_chejuensis_KCTC_2396.faa-HMM/FbpC-iron_uptake_ATPase-rep.hmm.tblout'
    hmmName <- gsub(".hmm.tblout","",basename(inPath))
    hmm_cutoff <- fegenie_hmm_cutoffs[fegenie_hmm_cutoffs$HMM==hmmName,'soft_bitscore_cutoff']
    hmm_cutoff <- as.numeric(paste0("1e-",hmm_cutoff))

    hmmloaded <- read_hmmsearch_tblout(inPath = inPath)

    filtered <- hmmloaded %>%
      mutate(sequence_evalue = as.numeric(sequence_evalue)) %>%
      mutate(domain_evalue = as.numeric(domain_evalue)) %>%
      filter(sequence_evalue<hmm_cutoff) %>%
      filter(domain_evalue<hmm_cutoff)


    if(nrow(filtered)>0){
      out.all <- rbind(out.all, cbind(filtered,genome))

    }
  }
}

dim(out.all)

fegenie_pan <- out.all %>% group_by(query_name, genome) %>% dplyr::count() %>%
  spread(query_name,n, fill=0) %>%
  data.frame()

rownames(fegenie_pan)= fegenie_pan[,1]
fegenie_pan = fegenie_pan[,c(2:ncol(fegenie_pan))]

head(fegenie_pan)

fegenie_hmm_2_desc <- out.all %>% select(query_name, description) %>% unique()




sel.hmm <- out.all %>% filter(grepl('vibrioferrin',query_name)) %>% select(query_name) %>% unique() %>% pull()

sel.pan <- fegenie_pan[,grepl('vibrioferrin_biosynthesis',colnames(fegenie_pan))]

fegenie_pan_pres_abs <- fegenie_pan
fegenie_pan_pres_abs[fegenie_pan_pres_abs>1] <- 1


fegenie_long <- fegenie_pan %>% rownames_to_column(var='genome') %>% pivot_longer(cols=c(-genome),names_to='hmm',values_to='value')
fegenie_long <- sel.pan %>% rownames_to_column(var='genome') %>% pivot_longer(cols=c(-genome),names_to='hmm',values_to='value')
fegenie_long <- fegenie_pan_pres_abs %>% rownames_to_column(var='genome') %>% pivot_longer(cols=c(-genome),names_to='hmm',values_to='value')


p.fegenie <- fegenie_long %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  #mutate(value=na_if(value,0)) %>%
  ggplot(aes(hmm,genome)) +
  geom_point(aes(fill=value),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="white", high="black",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#+

p.fegenie


ggsave('~/DATA/MarbGenomics/Graphs/FeGenie_.pdf',p.fegenie ,width=10,height = 15,units='in')




#=======

catdir <- list.dirs('/Users/sa01fd/Downloads/FeGenie-master/hmms/iron')
catdir <- catdir[!grepl('alignments',catdir)]
catdir <- catdir[2:length(catdir)]

fegenie_cat <- NULL
for(cdir in catdir){
  files <- list.files(cdir,pattern='.hmm')
  hmms <- gsub('.hmm','',files)
  if(length(hmms > 0 )){
    fegenie_cat <- rbind(fegenie_cat, cbind('category'=basename(cdir),
                                            'hmm'=hmms))

  }
}
fegenie_cat <- fegenie_cat %>% as_tibble()

fegenie_long <- fegenie_long %>% left_join(fegenie_cat, by='hmm')

fegenie_long %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  #mutate(value=na_if(value,0)) %>%
  ggplot(aes(hmm,genome)) +
  geom_point(aes(fill=value),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="white", high="black",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~category,scales='free_x',space = "free_x")

#======





fegenie_pan




#hmmloaded <- read_hmmsearch_tblout('/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Hahella_chejuensis_KCTC_2396.faa-HMM/FbpC-iron_uptake_ATPase-rep.hmm.tblout')

inPath='/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Hahella_chejuensis_KCTC_2396.faa-HMM/FbpC-iron_uptake_ATPase-rep.hmm.tblout'
hmmName <- gsub(".hmm.tblout","",basename(inPath))
hmm_cutoff <- fegenie_hmm_cutoffs[fegenie_hmm_cutoffs$HMM==hmmName,'soft_bitscore_cutoff']
hmm_cutoff <- as.numeric(paste0("1e-",hmm_cutoff))

hmmloaded <- read_hmmsearch_tblout(inPath = inPath)

hmmloaded %>%
  mutate(sequence_evalue = as.numeric(sequence_evalue)) %>%
  mutate(domain_evalue = as.numeric(domain_evalue)) %>%
  filter(sequence_evalue<hmm_cutoff) %>%
  filter(domain_evalue<hmm_cutoff)




hmmloaded <- read_hmmsearch_tblout('/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Marinobacter_algicola_DG893.faa-HMM/PvsC_vibrioferrin_biosynthesis-rep.hmm.tblout')


read_hmmsearch_tblout



readr::read_tsv('/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Hahella_chejuensis_KCTC_2396.faa-HMM/FbpC-iron_uptake_ATPase-rep.hmm.tblout',
                col_names=c("target","accession","query_name","accession","sEvalue","sscore","sbias","dEvalue","dscore","dbias","exp","reg","clu","ov","env","dom","rep","inc","description"), comment='#', na='-')




readr::read_tsv('/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Hahella_chejuensis_KCTC_2396.faa-HMM/FbpC-iron_uptake_ATPase-rep.hmm.tblout', comment='#', na='-')

readr::read_tsv('/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Hahella_chejuensis_KCTC_2396.faa-HMM/FbpC-iron_uptake_ATPase-rep.hmm.tblout', comment='#', na='-')

inPath = '/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Marinobacter_algicola_DG893.faa-HMM/PvsA_vibrioferrin_biosynthesis-rep.hmm.tblout'

"target","accession","query_name","accession2","sEvalue","sscore","sbias","dEvalue","dscore","dbias",
"exp","reg","clu","ov","env","dom","rep","inc","description"




#hmm_tblout2df
#===========================================================================================
#' Title
#'
#' @param inPath
#'
#' @return
#' @export
#'
#' @examples
#' hmmloaded <- read_hmmsearch_tblout('/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/Marinobacter_algicola_DG893.faa-HMM/PvsC_vibrioferrin_biosynthesis-rep.hmm.tblout')

read_hmmsearch_tblout <- function(inPath){
  col_types <-
    readr::cols(
      target              = readr::col_character(),
      accession           = readr::col_character(),
      query_name          = readr::col_character(),
      accession2          = readr::col_character(),
      sequence_evalue     = readr::col_double(),
      sequence_score      = readr::col_double(),
      sequence_bias       = readr::col_double(),
      domain_evalue       = readr::col_double(),
      domain_score        = readr::col_double(),
      domain_bias         = readr::col_double(),
      exp                 = readr::col_integer(),
      reg                 = readr::col_integer(),
      clu                 = readr::col_integer(),
      ov                  = readr::col_integer(),
      env                 = readr::col_integer(),
      dom                 = readr::col_integer(),
      rep                 = readr::col_integer(),
      inc                 = readr::col_integer(),
      description         = readr::col_character()
    )

  N <- length(col_types$cols)

  output = readr::read_lines(inPath) %>%
    sub(
      pattern = sprintf("(%s) *(.*)", paste0(rep('\\S+', N-1), collapse=" +")),
      replacement = '\\1\t\\2',
      perl = TRUE
    ) %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(col_names=c('X', 'description'), comment='#', na='-') %>%
    tidyr::separate(.data$X, head(names(col_types$cols), -1), sep=' +') #%>%
    #readr::type_convert(col_types=col_types)
  return(output)
}



#hmm_tblout2df
#===========================================================================================
#' Title
#'
#' @param inPath
#' @param evalue.cutoff
#' @param clean_genome_names
#'
#' @return
#' @export
#'
#' @examples
#' fegenie_out <- hmm_tblout2df(inPath='/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/', evalue.cutoff = 1.e10, clean_genome_names=TRUE)

hmm_tblout2df <-function(inPath='', evalue.cutoff = 1.e10, clean_genome_names=TRUE){
  hmm.df <- NULL
  genomDirs <- dir(inPath, full.names = T)

  for(x in genomDirs){
    if(clean_genome_names==TRUE){
      genome = gsub('_faa_HMM','',gsub('-','_',gsub('\\.','_', gsub('.out','',basename(x)))))
    }else{
      genome = gsub('_faa_HMM','',basename(x))
    }
    message(paste0('analysing: ',genome))

    hmmFiles <- list.files(path=x,full.names=TRUE)
    hmmFiles <- hmmFiles[1:(length(hmmFiles)-1)]

    for(y in hmmFiles){
      #Load the file
      hmm <- gsub("\\.hmm\\.tblout","",basename(y))
      hmm.out <- read_hmmsearch_tblout(inPath = y)
      if(nrow(hmm.out)>0){
        hmm.out <- hmm.out %>% dplyr::mutate(genome=genome) %>% dplyr::mutate(HMM=hmm)
        #Filter out low confidence hits based on evalue
        hmm.out <- hmm.out %>% dplyr::filter(sequence_evalue <= evalue.cutoff)
        hmm.df <- rbind(hmm.df, hmm.out)
      }

      }
  }
  return(hmm.df)
}






fegenie_out <- hmm_tblout2df(inPath='/Users/sa01fd/DATA/MarbGenomics/fegenie_out/HMM_results/',
                             evalue.cutoff = 1.e30,
                             clean_genome_names=TRUE)


p.magiclamp3 <- fegenie_out %>% group_by(genome,HMM) %>% tally() %>%
  filter(genome %in% tree$tip.label) %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  left_join(fegenie_cat, by=c('HMM'='hmm'))%>%
  mutate(genome=factor(genome,level=rev(tree$tip.label))) %>%
  ggplot(aes(HMM,genome))+geom_point(aes(color=n)) +
  #theme_minimal() +
  facet_grid(~category,scales='free_x',space = "free_x")+
  scale_fill_continuous(low="white", high="black",
                        guide="colorbar",na.value="white")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  ggsave('~/DATA/MarbGenomics/Graphs/iron_MB_fegenie.pdf',p.magiclamp3 ,width=15,height = 20,units='in')



fegenie_out %>% group_by(genome,HMM) %>% tally()   %>%
  filter(genome %in% tree$tip.label) %>%
  filter(grepl('Pvs',HMM)) %>%
  mutate(genome=factor(genome,level=rev(tree$tip.label))) %>%
  ggplot(aes(HMM,genome))+geom_point(aes(color=n))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
#facet_grid(~substrate,scales='free_x',space = "free_x")



#------------------------------------------
#
#	PREPARE
#
#------------------------------------------

#Module completion prediction uses module format as in KEGG description files
#   http://rest.kegg.jp/get/M00009


#some examples were tested
#modules = c('M00001','M00002','M00003','M00009','M00010','M00144','M00530')
#def = c('(K00844,K12407,K00845,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)',
#        'K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)',
#        '(K01596,K01610) K01689 (K01834,K15633,K15634,K15635) K00927 (K00134,K00150) K01803 ((K01623,K01624,K11645) (K03841,K02446,K11532,K01086,K04041),K01622)',
#        '(K01647,K05942) (K01681,K01682) (K00031,K00030) (K00164+K00658+K00382,K00174+K00175-K00177-K00176) (K01902+K01903,K01899+K01900,K18118) (K00234+K00235+K00236+K00237,K00239+K00240+K00241-(K00242,K18859,K18860),K00244+K00245+K00246-K00247) (K01676,K01679,K01677+K01678) (K00026,K00025,K00024,K00116)',
#        '(K01647,K05942) (K01681,K01682) (K00031,K00030)',
#        'K00330+(K00331+K00332+K00333,K00331+K13378,K13380)+K00334+K00335+K00336+K00337+K00338+K00339+K00340+(K00341+K00342)+K00343',
#        '(K00370+K00371+K00374,K02567+K02568) (K00362+K00363,K03385+K15876)'
#)

#module_table = data.frame(cbind('module_id'= modules,'DEFINITION'= def))
#module_table$DEFINITION = as.character(module_table$DEFINITION)

#show
#module_table$DEFINITION <- gsub("[;]"," ",module_table$DEFINITION)

#------------------------------------------
#
#
#------------------------------------------

#	http://rest.kegg.jp/link/module/ko
ko2mod = read.delim('~/DATA/MarbGenomics/AnnotateMachine/ko2module.txt',header=FALSE)
colnames(ko2mod) = c('ko','module')
ko2mod$ko = gsub('ko\\:','',ko2mod$ko)
ko2mod$module = gsub('md\\:','',ko2mod$module)
allmods = unique(ko2mod[,'module'])

#download from http://rest.kegg.jp/list/ko
ko2ann = read.delim('~/DATA/MarbGenomics/AnnotateMachine/KOlist.txt',header=FALSE,quote="")
colnames(ko2ann) = c('ko','ko_annotation')
ko2ann$ko = gsub('ko\\:','',ko2ann$ko)
head(ko2ann)

#--------------------------------------------------
#	Load the old KEGG modules (unfortionately not continued, so not available on KEGG API)

old_module_table = read.delim('~/DATA/MarbGenomics/AnnotateMachine/old_modules.txt',header=FALSE)
colnames(old_module_table) = c('module','type','DEFINITION')

m2n = read.delim('~/DATA/MarbGenomics/AnnotateMachine/old_modules2names.txt',header=FALSE)
colnames(m2n) = c("MOD_F", 'TYPE','DESCRIPTION')

m2n$module = gsub("\\_.*","", m2n$MOD_F)

m2n2 = read.delim('~/DATA/MarbGenomics/AnnotateMachine/recent_kegg_mods.txt',header=TRUE)
m2n2 <- m2n %>% left_join(m2n2, by=c("module"='Module.ID'))


#-------------------------------------------------
# Simple Function to obtain module information from KEGG REST API, powered by library(KEGGREST)
#	You need to have KEGGREST installed!
#	library(KEGGREST)


#keggREST_module <- function(modules){
#  output = c()
#  print('using library(KEGRGEST) to access KEGG REST API')
#  for(mod in modules){
#    print(paste(c('querying KEGG: ', mod),collapse=''))
#    query <- KEGGREST::keggGet(mod)
#    output = rbind(output, c('module'=mod, unlist(c(query[[1]][1],query[[1]][2],query[[1]][3],query[[1]][5]))))
#  }
#  output = data.frame(output)
#  output$module=as.character(output$module)
#  output$DEFINITION =as.character(output$DEFINITION)
#
#  return(output)
#}

#Usage
#module_table = keggREST_module(allmods)


#---------------------------------------
#
#	Helper functions for module completion calculation
#
#---------------------------------------

misc_module_subgroup_indexing <- function(DEFINITION,...){
  all_chars   <- strsplit(DEFINITION,split="+")[[1]] # split all characters (split by regular expression)
  Group_index <- numeric(length(all_chars)+1) # LEAVE ONE POSITION BEFORE
  # zeros
  for(C in 1:length(all_chars)){
    if(all_chars[C]=="("){
      Group_index[C+1] <- Group_index[C] + 1
    } else if(C>1&&all_chars[C-1]==")"){
      Group_index[C+1] <- Group_index[C] - 1  # If the previous one was a ')', change group, unless its a '('
    } else {
      Group_index[C+1] <- Group_index[C]
    }
  }
  # Remove first zero
  Group_index <- Group_index[-1]
  return(Group_index)
}


misc_module_definition_check <- function(DEFINITION,...){

  # DEFINE SUBGROUPING STRUCTURE
  DEFINITION <- gsub(" ,",",",DEFINITION) # CATCH LOOSE " ," (typo in definition catched 20180228)
  DEFINITION <- gsub(", ",",",DEFINITION)
  DEFINITION <- gsub(" +"," ",DEFINITION)
  DEFINITION <- gsub(" -- "," ",DEFINITION)
  DEFINITION <- gsub("^-- ","",DEFINITION)
  DEFINITION <- gsub(" --$","",DEFINITION)

  group_index <- misc_module_subgroup_indexing(DEFINITION)
  all_chars   <- strsplit(DEFINITION,split="+")[[1]] # split all characters (with regular expression)

  #### SPLIT WHERE ZEROS MATCH SPACES: BLOCKS!
  # IF NO SPACES MATCH A ZERO, SINGLE BLOCK!
  if(length(intersect(which(group_index==0),which(all_chars==" ")))==0){
    BLOCKS <- gsub(" ","+",DEFINITION)
  }else{
    split_here  <- intersect(which(group_index==0),which(all_chars==" "))
    nBlocks     <- length(split_here)+1
    BLOCKS      <- character(nBlocks)

    for(SG in 1:nBlocks){
      if(SG==1){
        thisSubBlock <- gsub(" ","+",substr(DEFINITION,1,split_here[SG]-1))
      } else if(SG>1&&SG<nBlocks){
        thisSubBlock <- gsub(" ","+",substr(DEFINITION,split_here[SG-1]+1,split_here[SG]-1))
      } else {
        thisSubBlock <- gsub(" ","+",substr(DEFINITION,split_here[SG-1]+1,nchar(DEFINITION)))
      }
      # REMOVE BRACKETS WHEN THE LAST CLOSES THE FIRST
      chars   <- strsplit(thisSubBlock,split="+")[[1]] # split all characters
      gp_ind  <- misc_module_subgroup_indexing(chars)

      if(sum(gp_ind==0)==0&&chars[1]=="("&&chars[length(chars)]==")"){
        # YES! - BALANCED
        # --> previous assumption that if they matched, they must be removed because they are flanking block! BUt only true if there is one at beginning and one at end!
        # 20180228: ADDED CHECK TO SEE THAT THERE ARE BRACKETS AT BEGINNING AND END!
        thisSubBlock <- gsub("^[()]","",thisSubBlock)
        BLOCKS[SG]    <- gsub("[)]$","",thisSubBlock)
      }else{
        BLOCKS[SG]    <- thisSubBlock
      }
    }
  }

  DEFINITION_OUT <- paste(BLOCKS,collapse = " ")

  # TRANSFORM DEFINITION TO LOGICAL EXPRESSION
  DEFINITION_OUT <- gsub("+","&",DEFINITION_OUT,fixed = T)
  DEFINITION_OUT <- gsub(",","|",DEFINITION_OUT,fixed = T)

  return(DEFINITION_OUT)
}


#------------------------------------------
#
#	MODULE COMPLETION RATIO PREDICTION
#
#------------------------------------------

#Function to calculate Module compeltion ratio given a list of KOs and a module_table

kegg_MCR <- function(KOs, module_table){
  output =  c()
  for(def in module_table$DEFINITION){
    #print(def)

    #Clean DEFINITION, by removing a potential whitespace at the end or start
    if(endsWith(def,' ')){
      def = substr(def,0,(nchar(def)-1))
    }

    if(startsWith(def," ")){
      def = substr(def,2,nchar(def))
    }

    #further clean DEFINITION
    cdef <- misc_module_definition_check(def)
    #print(cdef)
    blocks    <- strsplit(cdef,split = " ")[[1]]
    nr_of_blocks =  length(blocks)
    count_passed = 0

    for(blck in blocks){
      # blocks are either
      #	- horizontal only, equally suitable KOs, can be any of those --- K01647|K05942
      #	- vertical and horizontal (AND and OR) --- K01902&K01903|K01899&K01900|K18118
      #	- nested vertical and horizontal (AND and OR) --- K00239&K00240&K00241-(K00242|K18859|K18860)|K00244&K00245&K00246-K00247

      # We first asses wheter the current block is of type
      if(	(length(grep('&',blck))>0 |
           length(grep('\\(',blck))>0 |
           length(grep('\\)',blck))>0 |
           length(grep('\\+',blck))>0 |
           length(grep('\\-',blck))>0)
      ){
        #Check whether the block is nested

        if(length(grep('\\(', blck))>0){
          #------------------
          # vertical, nested
          #------------------

          blck = sub('&&','&', blck)
          blck = gsub('\\-','&', blck)
          blck = gsub('\\+','&', blck)

          #print('nested structure detected')
          #we first change all the "|" within paranthesis to "." allowing a alternatives split
          blck = gsub("\\|(?=[^()]*\\))", ".", blck, perl=TRUE)

          alternatives = unlist(strsplit(blck,split='\\|'))

          #--- test the alternatives
          check_alternatives = c()

          #The non nested ones

          for(alt in alternatives[!grepl('\\(', alternatives)]){
            #print(alt)
            mandatoryKOs = unlist(strsplit(alt,split='&'))

            #Check If All are TRUE
            allhit = all(ifelse(mandatoryKOs %in% KOs,TRUE,FALSE))
            check_alternatives = c(check_alternatives, allhit)
          }

          #The nested ones
          for(alt in alternatives[grepl('\\(', alternatives)]){
            nestedhits = c()
            mandatoryKOs = unlist(strsplit(alt,split='&'))

            #select nested
            nested = mandatoryKOs[grepl('\\(', mandatoryKOs)]
            for(nst in nested){
              #print(nst)
              nst = nst %>% gsub("\\(",'',.) %>% #remove brackets
                gsub("\\)",'',.) %>%
                gsub("\\.",'|',.) #return . to | to allow regular expression

              #check if any of these are present
              nestedhits = c(nestedhits, ifelse(length(grep(nst,KOs))>0, TRUE, FALSE))
            }

            #check the non nested ones
            nonnestedhit = mandatoryKOs[!grepl('\\(', mandatoryKOs)] %in% KOs
            combinedhits = c(nestedhits ,nonnestedhit)

            #Check If All are TRUE
            check_alternatives = c(check_alternatives, all(combinedhits))
          }

          #now check if check_alternatives holds at least one TRUE (one route through the reaction)
          if(any(check_alternatives)){
            count_passed = count_passed +1
          }else{
            #print('not passed')
          }
        }else{
          #------------------
          # vertical, non-nested
          #------------------
          #Cleanup, change +s odd "-"s and doubled &&s to single & reprententing "AND"
          blck = sub('&&','&', blck)
          blck = gsub('\\-','&', blck)
          blck = gsub('\\+','&', blck)

          #Here, the "|" splits the alternative routes
          alternatives = unlist(strsplit(blck,split='\\|'))

          #Run over each alternative, and check whether we have the KOs listed
          check_alternatives = c()
          for(alt in alternatives){
            mandatoryKOs = unlist(strsplit(alt,split='&'))

            #Check If All are TRUE
            allhit = all(ifelse(mandatoryKOs %in% KOs,TRUE,FALSE))
            check_alternatives = c(check_alternatives, allhit)
          }

          #now check if check_alternatives holds at least one TRUE (one route through the reaction)
          if(any(check_alternatives)){
            count_passed = count_passed +1
          }else{
            #print('not passed')
          }

        }

        #------------------
        # Horizontal Only
        #------------------

      }else{
        #print('clean - we can just simply look this use grep using this block')
        #print('any of these KOs can be present for the reaction/block to be accepted')
        if(length(grep(blck,KOs))>0 ){
          #print('PASSED')
          count_passed = count_passed +1
        } else {
          #print('not passed')
        }
      }

    } #----- End of BLOCKS for loop


    output = rbind(output, c('passed' = count_passed, 'total' = nr_of_blocks, 'MCR' = (count_passed/nr_of_blocks)*100 ))
    #print(paste(c(count_passed, '/', nr_of_blocks),collapse=''))
  }
  output = data.frame(cbind('module' = as.character(module_table$module) ,output),stringsAsFactors=FALSE)
  output$passed = as.numeric(as.character(output$passed))
  output$total = as.numeric(as.character(output$total))
  output$MCR = as.numeric(as.character(output$MCR))

  return(output)
}




#Example Usage
#inputHMM = read.delim('~/DATA/MarbGenomics/dbCAN_KOfam/Hahella_chejuensis_KCTC_2396_dbCAN_unique_for_cds_extraction.txt')
#KOs = as.character(inputHMM[,2])
#MCR2 = kegg_MCR(KOs, old_module_table)
#completeModules = MCR2[MCR2$MCR == 100,]
#dim(completeModules)

#--------------------------------------------------------------------------------------------
# kofam2MCR
#   function that goes and do it for you, specify the folder where you stored the KOfam output
#   this function will go and read files from dbCAN analysis based on D. Greens pipeline
#   this will take the KOs from there and go and predict MCR with the kegg_MRC function
#   kofan2MCR converts t
#-------------------------------------------------------------------------------------------

kofam2MCR <- function(inPath=NULL, module_table = old_module_table, clean_genome_names=TRUE){
  hmmFiles <- dir(inPath, pattern = '_unique_for_cds_extraction.txt', full.names = T)
  output = c()
  genomeList = c()
  for(x in hmmFiles){
    inputHMM = read.delim(x,header=FALSE)
    if(clean_genome_names==TRUE){
      genome = gsub('-','_',gsub('\\.','_', gsub('_dbCAN_unique_for_cds_extraction.txt','',basename(x))))
    }else{
      genome = gsub('_dbCAN_unique_for_cds_extraction.txt','',basename(x))
    }
    print(paste0('loading: ',genome))
    KOs = as.character(inputHMM[,2])
    MCR2 = kegg_MCR(KOs, module_table)
    output = cbind(output, MCR2$MCR)
    genomeList = c(genomeList, genome)
  }
  colnames(output) = genomeList
  rownames(output) = MCR2$module

  #format to have clean data.frame in "pan" format
  output = data.frame(t(output))
  return(output)
}


#usage
#output = kofam2MCR(inPath='~/DATA/MarbGenomics/dbCAN_KOfam/', clean_genome_names=TRUE)

#==================================================================================================
# kofam output from METABOLIC
# function that load in the metabolic KO data form KO like formatted files


KO_from_metabolic <- function(path){
  files<-list.files(path=path,
                    pattern='.hits.txt',
                    full.names=TRUE)
  KO_all <- NULL
  for(f in files){
    gnm = base::gsub('_hits_txt','',base::gsub('\\.','_',base::gsub('-','_',base::basename(f))))
    message('analysing ',gnm)
    KO <- utils::read.delim(f,header=FALSE)
    colnames(KO) <- c('KO','locus_tag')
    KO <- KO %>%
      dplyr::mutate(locus_tag = strsplit(as.character(locus_tag), ",")) %>%
      tidyr::unnest(locus_tag)
    KO_all = base::rbind(KO_all, KO %>% dplyr::mutate(genome=gnm))
  }
  return(KO_all)
}

# use function to read them all
KO_all <- KO_from_metabolic(path='/Users/sa01fd/DATA/MarbGenomics/metabolics_all2/KEGG_identifier_result')
#KO_all <- KO_all %>% filter(locus_tag!='')
head(KO_all)

#lc_2_gnm <- annotation.tbl %>% select(locus_tag, genome)
#KO_all %>% select(-genome) %>%




KO_annotated <- KO_all %>%
  left_join(ko2mod, by=c("KO"="ko")) %>%
  left_join(ko2ann, by=c("KO"="ko")) %>%
  left_join(m2n, by="module")



KO_annotated %>%
  filter(grepl('oxidoreductase',DESCRIPTION,ignore.case=T)) %>%
  ggplot(aes(ko_annotation,genome))+
  geom_tile()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p.l <- KO_annotated %>%
  mutate(genome = factor(genome,levels=rev(tip_order(lsTrees.106.rooted[[1]])))) %>%
  #filter(grepl('Multidrug resistance',ko_annotation ,ignore.case=T)) %>%
  filter(grepl('oxidoreductase',DESCRIPTION,ignore.case=T)) %>%
   ggplot(aes(ko_annotation,genome))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10)) + facet_wrap(~DESCRIPTION, scales='free_x')

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','tetstuff','.pdf'),plot=p.l, width = 10, height = 20,unit='in')



# Function to convert this into a pan formatted table
metabolic2MCRpan <- function(KO_all,module_table){
  output = c()
  genomeList = c()
  for(gnm in unique(KO_all$genome)){
    message(paste0('loading: ',gnm))
    KOs = KO_all  %>% filter(genome == gnm) %>% select(KO) %>% pull()
    MCR2 = kegg_MCR(KOs, module_table)
    output = cbind(output, MCR2$MCR)
    genomeList = c(genomeList, gnm)
  }
  colnames(output) = genomeList
  rownames(output) = MCR2$module

  #format to have clean data.frame in "pan" format
  output = data.frame(t(output))
  return(output)
}


MCR.all <- metabolic2MCRpan(KO_all, old_module_table)
rownames(MCR.all)[rownames(MCR.all)=='Marinobacter_pelagius_strain_CGMCC_16775'] = 'Marinobacter_pelagius_strain_CGMCC_1_6775'
write.table(MCR.all , file = "~/DATA/MarbGenomics/MCR_metabolic.pan.tsv",sep='\t')

MCR.pan <- MCR.all


KOs = as.character(KO[,1])
MCR2 = kegg_MCR(KOs, old_module_table)
completeModules = MCR2[MCR2$MCR == 100,]
dim(completeModules)




#build KO pan from metabolic
kfm.pan <- KO_annotated %>%
  dplyr::select(genome, KO) %>%
  dplyr::group_by(genome) %>%
  dplyr::count(KO) %>%
  spread(KO,n, fill=0) %>%
  data.frame()

rownames(kfm.pan) <- kfm.pan$genome
kfm.pan <- kfm.pan[,2:ncol(kfm.pan)]
kfm.pan[sel.genome,1:10]







#==================================================================================================
#  DONWSTREAM ANALYSIS
#==================================================================================================
#==================================================================================================
#!!!!!!!!! IF YOU HAVE TIME ELSE LOAD IT IN!
#MCRpan = kofam2MCR(inPath='~/DATA/MarbGenomics/dbCAN_KOfam/', module_table = old_module_table, clean_genome_names=TRUE)
#write.table(MCRpan , file = "~/DATA/MarbGenomics/MCR.pan.tsv",sep='\t')

#==================================================================================================


#load in kegg module metadata
m2n_filtered = m2n[,c('module','TYPE','DESCRIPTION')] %>% unique()

# LOAD PREVIOUS RUN
MCR.pan <- read.table("~/DATA/MarbGenomics/MCR.pan.tsv",sep='\t')
MCR.pan <- read.table("~/DATA/MarbGenomics/MCR_metabolic.pan.tsv",sep='\t')

#m.pan <- MCR.pan[sel.genome,]
m.pan <- MCR.pan




#m.pan <- MCR.pan[grepl('ydroc',rownames(MCR.pan)),]

nr_of_modules = ncol(m.pan)
core_modules_matrix = m.pan %>% as.data.frame() %>% filter(rowSums(.) == (100* nr_of_modules))

m.pan = m.pan %>% as.data.frame() %>% filter(rowSums(.) != 0)


# filter it some more
# WARNING<  NOT SURE ABOUT THIS STEP OF 75%
m.pan[m.pan<50]=0

m.pan = m.pan %>% as.data.frame() %>% filter(rowSums(.) != 0)
m.pan <- m.pan[,colSums(m.pan)!=0]

core.m.pan <- m.pan[,colSums(m.pan)==100*nrow(m.pan)]

Marinobacter_core.pam <- m.pan[,colSums(m.pan[grepl('Marinobacter',rownames(m.pan)),])==100*nrow(m.pan[grepl('Marinobacter',rownames(m.pan)),])]
Marinobacter_absent.pam <- m.pan[,colSums(m.pan[grepl('Marinobacter',rownames(m.pan)),])==0*nrow(m.pan[grepl('Marinobacter',rownames(m.pan)),])]


m.pan <- m.pan[,colSums(m.pan)!=100*nrow(m.pan)]

#core_modules = rownames(core_modules_matrix)
#core_modules_annotated = m2n_filtered %>% filter(module %in% core_modules)
#m2n_filtered %>% filter(module %in% 'M00005')


p.metabolic <- metabolic.tbl.long %>%
  filter(Gene.abbreviation %in% nonzero.genes) %>%
  filter(!grepl('acyl',Gene.abbreviation)) %>%
  filter(genome %in% tip.order) %>%
  mutate(value=na_if(value,0)) %>%


#======
tip.order <- tip_order(ladderize(lsTrees.106.rooted[[1]]))

#mcr.annotated <- mp.pan
m.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='module',values_to='MCR') %>%
  left_join(m2n_filtered, by='module') %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
  filter(!is.na(TYPE)) %>%
  ggplot(aes(module,genome)) +
  geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~TYPE,scales='free_x',space = "free_x")


p.mcr.pan <- m.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='module',values_to='MCR') %>%
  left_join(m2n_filtered, by='module') %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
  filter(!is.na(TYPE)) %>%
  ggplot(aes(DESCRIPTION,genome)) +
  geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="black",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#coord_flip()+
  facet_grid(~TYPE,scales='free_x',space = "free_x")

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_update','.pdf'),plot=p.mcr.pan, width = 30, height = 20,unit='in')




#===
selected.type = 'Complex'
selected.type = 'Pathway'
selected.type = 'FuncSet'


selected_modules <- m.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='module',values_to='MCR') %>%
  left_join(m2n_filtered, by='module') %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
    filter(!is.na(TYPE)) %>%
  filter(TYPE==selected.type) %>%
  filter(grepl("Transport",DESCRIPTION,ignore.case=TRUE)) %>%
  select(module) %>%
  pull() %>% unique

selected_modules <- m2n2 %>% filter(grepl("ATP synthesis",Module.Category, ignore.case=TRUE)) %>% select(module) %>% pull()
selected_modules <- intersect(colnames(m.pan), selected_modules)

module.clust <- hclust(dist(t(m.pan[,selected_modules])))
module.order <- module.clust$labels[module.clust$order]

namesConv <- m2n_filtered$DESCRIPTION
names(namesConv) = m2n_filtered$module
description.order <- unique(as.character(unname(namesConv[module.order])))

p.mcr.pan <- m.pan %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='module',values_to='MCR') %>%
  left_join(m2n_filtered, by='module') %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
  filter(genome %in% tip.order) %>%
  filter(module %in% selected_modules) %>%
  mutate(DESCRIPTION = factor(DESCRIPTION,level=description.order))%>%
  ggplot(aes(DESCRIPTION,genome)) +
  geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="black",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#coord_flip()+
  facet_grid(~TYPE,scales='free_x',space = "free_x")


if(selected.type == 'Pathway'){
  ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_',selected.type,'.pdf'),plot=p.mcr.pan, width = 22, height = 20,unit='in')
}else{
  ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_',selected.type,'.pdf'),plot=p.mcr.pan, width = 10, height = 15,unit='in')
}







#=======

all_categories <- m2n2 %>% select(Module.Category) %>% filter(!is.na(Module.Category)) %>% filter(Module.Category!="") %>% pull() %>% as.character() %>% unique()


for(cat in all_categories){
  message("plotting  ", cat)
  selected_modules <- m2n2 %>% filter(grepl(cat,Module.Category, ignore.case=TRUE)) %>% select(module) %>% pull()
  selected_modules <- intersect(colnames(m.pan), selected_modules)

  if(length(selected_modules)>0){
    if(length(selected_modules)>1){
      module.clust <- hclust(dist(t(m.pan[,selected_modules])))
      module.order <- module.clust$labels[module.clust$order]

      namesConv <- m2n_filtered$DESCRIPTION
      names(namesConv) = m2n_filtered$module
      description.order <- unique(as.character(unname(namesConv[module.order])))

      p.mcr.pan <- m.pan %>% rownames_to_column(var='genome') %>%
        pivot_longer(!genome,names_to='module',values_to='MCR') %>%
        left_join(m2n_filtered, by='module') %>%
        mutate(genome = factor(genome,level=rev(tip.order)))%>%
        mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
        filter(genome %in% tip.order) %>%
        filter(module %in% selected_modules) %>%
        mutate(DESCRIPTION = factor(DESCRIPTION,level=description.order))%>%
        ggplot(aes(DESCRIPTION,genome)) +
        geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
        theme_minimal() +
        scale_fill_continuous(low="darkgrey", high="black",
                              guide="colorbar",na.value="white")+
        #scale_color_manual(na.value='white')
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        facet_grid(~TYPE,scales='free_x',space = "free_x")
    }else{
      p.mcr.pan <- m.pan %>% rownames_to_column(var='genome') %>%
        pivot_longer(!genome,names_to='module',values_to='MCR') %>%
        left_join(m2n_filtered, by='module') %>%
        mutate(genome = factor(genome,level=rev(tip.order)))%>%
        mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
        filter(genome %in% tip.order) %>%
        filter(module %in% selected_modules) %>%
        ggplot(aes(DESCRIPTION,genome)) +
        geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
        theme_minimal() +
        scale_fill_continuous(low="darkgrey", high="black",
                              guide="colorbar",na.value="white")+
        #scale_color_manual(na.value='white')
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#coord_flip()+
        facet_grid(~TYPE,scales='free_x',space = "free_x")
    }


    cat <- gsub('\\/','_',cat)
    ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_',cat,'.pdf'),plot=p.mcr.pan, width = 10, height = 15,unit='in')

  }else{
    message('-- ',cat,': nothing detected ')
  }
}

#====

p.core.mcr.pan <- Marinobacter_core.pam %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='module',values_to='MCR') %>%
  left_join(m2n_filtered, by='module') %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
  ggplot(aes(DESCRIPTION,genome)) +
  geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#coord_flip()+
  facet_grid(~TYPE,scales='free_x',space = "free_x")

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_pan_update_core','.pdf'),plot=p.core.mcr.pan, width = 10, height = 20,unit='in')

p.abs.mcr.pan <- Marinobacter_absent.pam %>% rownames_to_column(var='genome') %>%
  pivot_longer(!genome,names_to='module',values_to='MCR') %>%
  left_join(m2n_filtered, by='module') %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(MCR=ifelse(MCR==0,NA,MCR))%>%
  ggplot(aes(DESCRIPTION,genome)) +
  geom_point(aes(fill=MCR),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
  #scale_color_manual(na.value='white')
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+#coord_flip()+
  facet_grid(~TYPE,scales='free_x',space = "free_x")

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','MCR_absent_in_MB','.pdf'),plot=p.abs.mcr.pan, width = 10, height = 20,unit='in')



#----------------------------------------------

#library(ggfortify)

df <- m.pan
pca_res <- prcomp(df)
#	pca_res <- prcomp(df, scale. = TRUE) --> in this case we don't need scaling, because it is all 0 to 100
library(ggfortify)
autoplot(pca_res, label = TRUE, label.size = 3, loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)

df_out <- as.data.frame(pca_res$x)
df_out <- df_out %>% rownames_to_column(var='genome') %>% left_join(genome.tbl, by='genome')

df_out %>% ggplot(aes(PC1,PC2,color=phylogroup))+geom_point()+scale_color_manual(values=phylogroup.colors)

#
#	Clade specific complete modules
#
#----------------------------------------------


# This is a method I invent, to run over nodes in a tree, identify the core modules (100%), identify the unique modules (only in 100% shared in this clade), and those that are absent in this clade but present in all others.
#this method is a modification of what I intend to implement for pangenomes in general, so perhaps, lets make it an presence absence table
#hard cuts, say all below 100, as absent, 100 as present...

#come up with a better name, now
#node_pangenomics


#FROM MODULES TO PRESENCE ABSENCE TABLE

pan_modules =  t(non_zero_modulues)

pan_modules.presabs = pan_modules
pan_modules.presabs[pan_modules.presabs < 100] <- 0
pan_modules.presabs[pan_modules.presabs == 100] <- 1
pan_modules.presabs = pan_modules.presabs[,colSums(pan_modules.presabs)>0]


#MAIN INPUTS
pan.presabs = pan_modules.presabs
tree = read.tree("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile")

#library(phytools)

outgroup = tree$tip.label[!grepl("Marinobacter", tree$tip.label)]
outgroup = outgroup[! outgroup %in% 'Tamilnaduibacter_salinus_Mi_7']

tree.rooted = root(tree, outgroup = outgroup, resolve.root = TRUE)
p <- ggtree(tree.rooted) + geom_tiplab()




#--------------------

treeTibble = p$data %>% select(c(parent,node, branch.length,label,isTip))
nodesNumbers = treeTibble %>% filter(isTip == FALSE) %>% select(node) %>% pull
tipNumbers = treeTibble %>% filter(isTip == TRUE) %>% select(node) %>% pull
tips = treeTibble %>% filter(isTip == TRUE) %>% select(label) %>% pull
names(tips) = tipNumbers


#need to clean up the dots!
rownames(pan.presabs) = gsub('\\.','_',rownames(pan.presabs))
#Ok this is silly... and needs reconsideration! but apparently, the tree tips dropped after CGMCC_1

rownames(pan.presabs)[rownames(pan.presabs)=='Marinobacter_pelagius_strain_CGMCC_1_6775'] = "Marinobacter_pelagius_strain_CGMCC_1"


#pan.presabs = pan
#pan.presabs[pan.presabs > 0] <- 1
#pan.presabs = t(pan.presabs)

for(node in nodesNumbers){
  descendants = getDescendants(tree, node, curr=NULL)
  if(length(descendants)!=0){

    desc_tips = descendants[descendants %in%  tipNumbers]
    desc_tipNames = tips[desc_tips]
    nr_of_desc_tips = length(desc_tips)


    pan.ingroup = pan.presabs[unname(desc_tipNames), ]
    pan.outgroup = pan.presabs[!(rownames(pan.presabs) %in% desc_tipNames),]

    #Those functional units that are consistently present in the ingroup, CORE functional units
    ingroup.core = pan.ingroup[,colSums(pan.ingroup)==nrow(pan.ingroup)]
    ingroup.core.functions = colnames(ingroup.core)

    #Those functional units which are consistently absent in the ingroup
    ingroup.absent = pan.ingroup[,colSums(pan.ingroup)==0]
    ingroup.absent.functions = colnames(ingroup.absent)

    #CORE in outrgroup
    outgroup.core = pan.outgroup[,colSums(pan.outgroup)==nrow(pan.outgroup)]
    outgroup.core.functions = colnames(outgroup.core)

    #ABSENT in outgroup
    outgroup.absent = pan.outgroup[,colSums(pan.outgroup)==0]
    outgroup.absent.functions = colnames(outgroup.absent)

    #------
    # Uniquely present in the ingroup as compared to outgroup
    uniquely.present.ingroup.functions = ingroup.core.functions[ingroup.core.functions %in% outgroup.absent.functions]

    #Uniquely present in outgroup as compared to ingroup
    uniquely.present.outgroup.functions = outgroup.core.functions[outgroup.core.functions %in% ingroup.absent.functions]

    if(length(uniquely.present.ingroup.functions) > 0){
      print('-----UNIQUE INGROUP-----------------------------------------------------')
      print(paste(c('UNIQUE_INGROUP',node, uniquely.present.ingroup.functions),collapse='; '))
      print('INGROUP:')
      print(unname(desc_tipNames))
      print('OUTGROUP:')
      print(rownames(pan.outgroup))

    }
    if(length(uniquely.present.outgroup.functions) > 0){
      print('-----UNIQUE OUTGROUP-----------------------------------------------------')
      print(paste(c('UNIQUE_OUTGROUP',node, uniquely.present.outgroup.functions),collapse='; '))
      print('INGROUP:')
      print(unname(desc_tipNames))
      print('OUTGROUP:')
      print(rownames(pan.outgroup))
    }
  }

  # 	Get Sister Nodes! and see if the sister node has dropped it or not
  #	getSisters(tree,200,mode="number")

  sister_node = getSisters(tree, node, mode="number")
  if(length(sister_node)!=0){
    sister_descendants = getDescendants(tree, sister_node, curr=NULL)

    sister_desc_tips = sister_descendants[sister_descendants %in%  tipNumbers]
    sister_desc_tipNames = tips[sister_desc_tips]
    sister_nr_of_desc_tips = length(sister_desc_tips)

    pan.sister.ingroup = pan.presabs[unname(sister_desc_tipNames), ]

    if(!is.null(dim(pan.sister.ingroup))){
      #Those functional units that are consistently present in the ingroup, CORE functional units
      sister.core = pan.ingroup[,colSums(pan.sister.ingroup)==nrow(pan.sister.ingroup)]
      sister.core.functions = colnames(sister.core)

      #Those functional units which are consistently absent in the ingroup
      sister.absent = pan.sister.ingroup[,colSums(pan.sister.ingroup)==0]
      sister.absent.functions = colnames(sister.absent)

      #------
      # Uniquely present in the ingroup as compared to outgroup
      present.ingroup.absent.sister.functions = ingroup.core.functions[ingroup.core.functions %in% sister.absent.functions]

      #Uniquely present in outgroup as compared to ingroup
      present.sister.absent.ingroup.functions = sister.core.functions[sister.core.functions %in% ingroup.absent.functions]

      if(length(uniquely.present.ingroup.functions) > 0){
        print('----- CORE ingroup / absent Sister -----------------------------------------------------')
        print(paste(c('CORE ingroup/absent Sister',node, uniquely.present.ingroup.functions),collapse='; '))
        print('INGROUP:')
        print(unname(desc_tipNames))
        print('SISTER:')
        print(rownames(pan.sister.ingroup))
      }
      if(length(uniquely.present.outgroup.functions) > 0){
        print('----- CORE sister / absent ingroup -----------------------------------------------------')
        print(paste(c('CORE sister / absent ingroup',node, uniquely.present.outgroup.functions),collapse='; '))
        print('INGROUP:')
        print(unname(desc_tipNames))
        print('SISTER:')
        print(rownames(pan.sister.ingroup))
      }
    }else{
      print('sister lineage is a singleton')
    }

  }else{print('no sister')}
}

m2n_filtered %>% filter(module %in% c('M00497'))


getDescendants(tree, 128, curr=NULL)


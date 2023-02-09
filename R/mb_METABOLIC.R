##-======

metabolic.tbl <- read.delim('/Users/sa01fd/DATA/MarbGenomics/metabolics_all2/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv') %>%
  as_tibble()


metabolic.tbl %>% select(Marinobacter_adhaerens_HP15.Hit.numbers, Gene.name, Hahella_ganghwensis_DSM_17046.Hits)


metabolic.tbl.long <- metabolic.tbl %>%
  pivot_longer(cols=contains('Hit.numbers'),names_to='genome') %>%
  select(Category:Hmm.detecting.threshold,genome,value) %>%
  mutate(genome=gsub('\\.Hit\\.numbers',"",genome)) %>%
  mutate(genome=gsub('\\.',"_",gsub('-',"_",genome)))

metabolic.tot.counts <- metabolic.tbl.long %>%
  group_by(Gene.abbreviation,Corresponding.KO) %>%
  dplyr::summarise(n=sum(value))

nonzero.genes <- metabolic.tot.counts %>% dplyr::filter(n>0) %>% select(Gene.abbreviation) %>% pull() %>% as.character()



metabolic.tot.counts %>% dplyr::filter(n>0) %>% select(Gene.abbreviation) %>% pull() %>% as.character()



#function to collect the tip order of a tree!
tip_order <- function(tree){
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  ordered_tips <- tree$tip.label[ordered_tips]
  return(ordered_tips)
}

tip.order <- tip_order(sco.106.concat.nuc.tree)

p.metabolic <- metabolic.tbl.long %>%
  filter(Gene.abbreviation %in% nonzero.genes) %>%
  filter(!grepl('acyl',Gene.abbreviation)) %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  mutate(value=na_if(value,0)) %>%
  ggplot(aes(Gene.abbreviation,genome)) +
    geom_point(aes(fill=value),size=2,shape=21,color='grey')+
    theme_minimal() +
    scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
    #scale_color_manual(na.value='white')
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_grid(~Category,scales='free_x',space = "free_x")

ggsave('~/DATA/MarbGenomics/Graphs/MATABOLIC_fun_for_joy.pdf',p.metabolic ,width=20,height = 15,units='in')


#Caultulate percentages per linaege based on presence absence of a gene
p.metabolic2 <- metabolic.tbl.long %>%
  filter(Gene.abbreviation %in% nonzero.genes) %>%
  filter(!grepl('acyl',Gene.abbreviation)) %>%
  filter(genome %in% tip.order) %>% left_join(genome.tbl %>%
  select(genome,phylogroup),by='genome') %>%
  mutate(value=ifelse(value>0,1,0))%>%
  mutate(gn.count =1)%>%
  group_by(phylogroup,Gene.abbreviation,Category) %>%
  dplyr::summarise(freq=sum(value)/sum(gn.count)) %>% ungroup %>%
  ggplot(aes(Gene.abbreviation,phylogroup)) +
  geom_point(aes(fill=phylogroup,alpha=freq),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_manual(values=phylogroup.colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(~Category,scales='free_x',space = "free_x")

ggsave('~/DATA/MarbGenomics/Graphs/MATABOLIC_per_phylogroup.pdf',p.metabolic2 ,width=20,height = 6.5,units='in')



mglmp_df <- annotation.tbl %>% filter(OG %in% m1$OG) %>%
  select(OG, genome) %>% left_join(m1,by='OG') %>%
  group_by(OG, genome, Name, product, module) %>%
  tally() %>% ungroup() %>%
  filter(genome %in% ordered_tips) %>%
  mutate(genome=factor(genome,level=rev(ordered_tips)))

#=====

#MANUALLY BUILDING HYDROCARBON METABOLISM
# notes on alkB,
# removal of Alpha-ketoglutarate-dependent dioxygenase AlkB

annotation.tbl %>% filter(grepl("alkb",Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG) %>%
  filter(!grepl("ketoglutarate",product,ignore.case = TRUE)) %>%
  mutate(Name=gsub('_.*', '', Name)) %>% group_by(genome, product, Name) %>% tally() %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  ggplot(aes(Name,genome)) +
  geom_point(aes(fill=n),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
 # scale_fill_manual(values=phylogroup.colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
  #facet_grid(~Category,scales='free_x',space = "free_x")

annotation.tbl %>% filter(grepl("adh1",Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG)



#========== CHECK THIS AS I ONLY PICK UP LITTLE
#adh1 long-chain alcohol dehydrogenase

sel.ogs <- annotation.tbl %>% filter(grepl("adh1",Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG) %>% select(OG) %>% filter(!is.na(OG)) %>% unique() %>% pull


annotation.tbl %>% filter(grepl("adh1",Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG)%>%
  mutate(Name=gsub('_.*', '', Name)) %>% group_by(genome, product, Name) %>% tally() %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  ggplot(aes(Name,genome)) +
  geom_point(aes(fill=n),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
  # scale_fill_manual(values=phylogroup.colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
#facet_grid(~Category,scales='free_x',space = "free_x")



#
# Build strategy, say check all the ogs.
# see if there is consistence in annotation

search <- 'FadD' #Checked and ready to go
search <- 'FadE'
search <- 'PaaF'
search <- 'FadJ'
search <- 'FadA'

#Naphthalene
search <- 'NahA' # no hits, need other method (expected in MB)
search <- 'NahB'
search <- 'NahC'
search <- 'NahE'
search <- 'NahF'
#sel.ogs <- annotation.tbl %>% filter(grepl('salic',product,ignore.case = TRUE)) %>%
#  select(genome, product, Name, locus_tag, OG)
#sel.ogs

# Catechol
# checked and ready
search <- 'catA|catB|catC'

#Phenathrene ---- UNCLEAR
search <- 'NidA' #nope
search <- 'PhnB|PhdE' #nope
search <- 'PhdF'#nope
search <- 'PhdG'#nope
search <- 'NidD'#nope
search <- 'PhdJ'#nope
search <- 'PhdK' #--------? [perhaps]


# Phthalic acid
search <- 'pht'

# Protocatechuate
search <- 'pca'

# so far for the nature microbio paper, how shit
#  lets continue with my own inferences

#paa cluster
search <- 'paa' # splendid


search <- 'nuo' # STRANGE CHECK BETTER
search <- 'atp'


#T6SS UTTER FAILURE
search <- 'IcmF|ImpA|ImpB|ImpC|ImpG|ImpH|ImpI|ImpJ|ImpK|VasD|VasH|VasI|Vrg'
search.prod <- 'Type VI'
search.prod <- 'secretion'

#
search.prod <-"type II secretion"
"type I secretion"

search.prod <- "Flagellar secretion"

search.prod <- "NADH"
search.prod <- "superoxide"
search.prod <- "superoxide"


search.prod <- "c551|peroxidas|Catalase|superoxide dismutase"
search.prod <- "phosphate|phosphatase|Apolipoprotein|transhydrogenase"


#TCA
search <- 'acnB|CLb|fumA|FumHan|gltA|IDHI|lpdA|lpdA2|mdh|sdhA|sdhB|sucA|sucB|sucC|sucD'


search <- 'waa'
search <- 'cys|sir'

search.prod <- "glycerol"
search <- 'glp'


#=========
#SPECIAL CASE FOR TRANSPOSASE, NO NAMES, generate them absed on family
search.prod <- "transposase"

search <- annotation.tbl %>% filter(grepl(search.prod,product,ignore.case = TRUE)) %>%
  select(product, Name, locus_tag, OG) %>% unique %>% select(Name) %>% filter(!is.na(Name))%>% pull() %>% unique %>% paste(collapse="|")

annotation.tbl %>% filter(grepl(search.prod,product,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG) %>% mutate(Name = gsub('family transposase',"",product))%>%
  mutate(Name = gsub(' .*',"",Name)) %>% mutate(Name=gsub('_.*', '', Name)) %>%
  group_by(genome, product, Name) %>% tally() %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  ggplot(aes(Name,genome)) +
  geom_point(aes(fill=n),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
  # scale_fill_manual(values=phylogroup.colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+


#=========
#~Extract from product
search <- annotation.tbl %>% filter(grepl(search.prod,product,ignore.case = TRUE)) %>%
  select(product, Name, locus_tag, OG) %>% unique %>% select(Name) %>% filter(!is.na(Name))%>% pull() %>% unique %>% paste(collapse="|")


annotation.tbl %>% filter(grepl(search.prod,product,ignore.case = TRUE)) %>%
  select(product, Name, locus_tag, OG) %>% unique


annotation.tbl %>% filter(grepl(search.prod,product,ignore.case = TRUE))  %>%
  select(genome, product, Name, locus_tag, OG) %>% select(OG) %>% filter(!is.na(OG)) %>% unique() %>% pull
annotation.tbl %>% filter(OG %in% sel.ogs) %>% select(genome, product, Name, locus_tag, OG) %>% data.frame %>% arrange(OG)
#check consistency


sel.ogs <- annotation.tbl %>% filter(grepl(search,Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG) %>% select(OG) %>% filter(!is.na(OG)) %>% unique() %>% pull
annotation.tbl %>% filter(OG %in% sel.ogs) %>% select(genome, product, Name, locus_tag, OG) %>% data.frame %>% arrange(OG)
#check consistency

annotation.tbl %>% filter(OG %in% sel.ogs) %>%
  select(genome, product, Name, locus_tag, OG) %>% data.frame %>%
  arrange(OG) %>% group_by(OG,product,Name) %>%
  tally()

#  show based on product annotation
annotation.tbl %>% filter(grepl(search,Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG)

# show based on OrthoFinder clusters
annotation.tbl %>% filter(OG %in% sel.ogs) %>%
  select(genome, product, Name, locus_tag, OG)


annotation.tbl %>% filter(grepl(search,Name,ignore.case = TRUE)) %>%
  select(genome, product, Name, locus_tag, OG) %>%
  mutate(Name=gsub('_.*', '', Name)) %>% group_by(genome, product, Name) %>% tally() %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=rev(tip.order)))%>%
  ggplot(aes(Name,genome)) +
  geom_point(aes(fill=n),size=2,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="red",
                        guide="colorbar",na.value="white")+
  # scale_fill_manual(values=phylogroup.colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
#facet_grid(~Category,scales='free_x',space = "free_x")



#=============





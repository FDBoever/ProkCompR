#

library(dplyr)

bacdiv = read.delim('~/DATA/MarbGenomics/BacDiv_TaxOrgigin.txt')
habitat = read.delim('~/DATA/MarbGenomics/habitat_undeveloped.txt')


bacdiv = as_tibble(bacdiv)
habitat = as_tibble(habitat)

bacdiv = bacdiv %>% filter(genome_name_dotted!='')
habitat  = habitat %>% filter(genome_name_dotted!='')



mergedEnv = habitat %>% left_join( bacdiv ,by = 'genome_name_dotted')

write.table(mergedEnv,file='~/DATA/MarbGenomics/mergedEnv.txt',quote=FALSE,sep='\t')




mergedEnv = read.delim('~/DATA/MarbGenomics/mergedEnv.txt')

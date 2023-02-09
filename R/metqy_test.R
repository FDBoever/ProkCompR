library(MetQy)

data("data_example_KOnumbers_vector") # Load data
OUT <- query_missingGenes_from_module(data_example_KOnumbers_vector,"M00010")
OUT


OUT <- query_missingGenes_from_module(data_example_KOnumbers_vector,"M00046")
OUT


selKO <- KO_annotated %>% filter(genome=='Marinobacter_algicola_DG893') %>% select(KO) %>% pull %>% unique()

seldf <- data.frame(cbind('name'='DG',genes=paste(selKO,collapse=';')))
seldf$gene = as.character(seldf$gene )

query_genomes_to_modules(GENOME_INFO = seldf
                        )


testdf <-  KO_annotated %>% select(genome,KO) %>%
  group_by(genome) %>%
  summarise(test = toString(KO)) %>%
  ungroup()

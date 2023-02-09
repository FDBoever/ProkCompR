
#-------------------------------------------------------------------------------------#

ncbimetadata <-read.delim('/Users/sa01fd/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/NCBI_assembly_metadata.txt',header=FALSE,quote="")

genome2ncbi <- ncbimetadata %>%
  dplyr::left_join(genome.tbl, by=c('V19'='Accession')) %>%
  dplyr::select(genome, V3, V4, V7, V37) %>%
  dplyr::filter(!is.na(genome)) %>% mutate(acc_full = paste0(V4,'_',V7)) %>%
  dplyr::mutate(url_faa = paste0(V37,'/',acc_full,'_protein.faa.gz')) %>%
  dplyr::mutate(url_feature = paste0(V37,'/',acc_full,'_feature_table.txt.gz')) %>%
  dplyr::mutate(url_gbff = paste0(V37,'/',acc_full,'_genomic.gbff.gz'))


#-------------------------------------------------------------------------------------#

out.gnm.tbl <- genome.tbl %>%
  dplyr::left_join(ncbimetadata, by=c('Accession'='V19')) %>%
  dplyr::filter(manualy_removed == FALSE) %>%
  dplyr::select(-Translation_table,-manualy_removed,-c0,c1,c2,c3,c4,c5plus) %>%
  mutate(acc_full = paste0(V4,'_',V7)) %>%
  dplyr::mutate(url_faa = paste0(V37,'/',acc_full,'_protein.faa.gz')) %>%
  dplyr::mutate(url_feature = paste0(V37,'/',acc_full,'_feature_table.txt.gz')) %>%
  dplyr::mutate(url_gbff = paste0(V37,'/',acc_full,'_genomic.gbff.gz'))

write.table(x=out.gnm.tbl,
            file="/Users/sa01fd/CHAPTERS/Marinobacter_manuscript_2022/Figures/SupplementaryTables/SI_table_1.txt",
            quote=FALSE,
            row.names = FALSE,
            sep='\t')

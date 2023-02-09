### BLAST BACK ####
# not some genome names in fasta have '-' or '.' instead of '_' - this was dealth with by manual checking

download_until_success <- function(url, destfile, maxcount = 5) {
  count <- 0
  repeat{
    Sys.sleep(0.5)
    try(download.file(url, destfile))
    count <- count + 1
    if (file.exists(destfile) || count >= maxcount)
      break
  }
}


# OBTAIN ASSEMBLY METADATA FROM EDIRECT SCRIPTS

ncbimetadata <-read.delim('/Users/sa01fd/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/NCBI_assembly_metadata.txt',header=FALSE,quote="")

genome2ncbi <- ncbimetadata %>%
  dplyr::left_join(genome.tbl, by=c('V19'='Accession')) %>%
  dplyr::select(genome, V3, V4, V7, V37) %>%
  dplyr::filter(!is.na(genome)) %>% mutate(acc_full = paste0(V4,'_',V7)) %>%
  dplyr::mutate(url_faa = paste0(V37,'/',acc_full,'_protein.faa.gz')) %>%
  dplyr::mutate(url_feature = paste0(V37,'/',acc_full,'_feature_table.txt.gz')) %>%
  dplyr::mutate(url_gbff = paste0(V37,'/',acc_full,'_genomic.gbff.gz'))


#### DOWNLOAD PROTEIN DATA FROM NCBI ####
# generate directories for each genome in blast_back directory
# copy the prokka genome to be blasted inside the directory
# download protein.faa, feature table and gbff files from NCBI
# R.utils::gunzip used to unzip the downloaded files

#set prokka_dir from which to copy files
prokka_dir = '~/DATA/MarbGenomics/prokka_faa'


# Generate output directory 'blast_back'
outdir <- '~/DATA/MarbGenomics/blast_back'
dir.create(outdir)

for(i in 1:nrow(genome2ncbi)){
  selected.g <- genome2ncbi[i,'genome']
  dir.create(paste0(outdir,'/',selected.g))

  #copy from prokka dir
  tryCatch(
    file.copy(from=paste0(prokka_dir,'/',selected.g,'.faa'),to=paste0(outdir,'/',selected.g)),
    error = function(c) "error")

  tryCatch({
          download_until_success(
                        url=gsub(' ','_',genome2ncbi[i,'url_faa']),
                        destfile = paste0(outdir,'/', selected.g,'/',gsub(' ','_',genome2ncbi[i,'acc_full']),'_protein.faa.gz'))
           R.utils::gunzip(paste0(outdir,'/',selected.g,'/',gsub(' ','_',genome2ncbi[i,'acc_full']),'_protein.faa.gz'))
           },
           error = function(c) "error")
  tryCatch({
          download_until_success(
                        url=gsub(' ','_',genome2ncbi[i,'url_feature']),
                        destfile = paste0(outdir,'/',selected.g,'/',gsub(' ','_',genome2ncbi[i,'acc_full']),'_feature_table.txt.gz'))
            R.utils::gunzip(paste0(outdir,'/',selected.g,'/',gsub(' ','_',genome2ncbi[i,'acc_full']),'_feature_table.txt.gz'))
            },
           error = function(c) "error")
  tryCatch({
          download_until_success(
                        url=gsub(' ','_',genome2ncbi[i,'url_gbff']),
                         destfile = paste0(outdir,'/',selected.g,'/',gsub(' ','_',genome2ncbi[i,'acc_full']),'_genomic.gbff.gz'))
            R.utils::gunzip(paste0(outdir,'/',selected.g,'/',gsub(' ','_',genome2ncbi[i,'acc_full']),'_genomic.gbff.gz'))
          },
           error = function(c) "error")
}






#==============================================================================#
#### read blast output ####


read.blastout <- read.delim('~/DATA/MarbGenomics/blastblack/GCF_000170835.1_ASM17083v1.results.out',header=FALSE,quote="")

read.blastout %>% group_by(V1) %>% arrange(V3) %>% slice(1)


list.out.files <- list.files(path = outdir,
                             pattern='results.out',
                             full.names=TRUE,
                             include.dirs = TRUE,
                             recursive=TRUE)



#==============================================================================#
### READING AND PARSING THE DATA #####
#library('genbankr')

list.out.dirs <- list.dirs(path=outdir, full.names=TRUE,recursive=TRUE)
list.out.dirs <- list.out.dirs[2:length(list.out.dirs)]


prokka2ncbi.list <- list()
ncbifeautures.list <- list()
ncbigbff.list <- list()

for(dir in list.out.dirs){
  s.gnm <- basename(dir)

  list.out.files <- list.files(path = dir,
                               pattern='results.out',
                               full.names=TRUE,
                               include.dirs = TRUE,
                               recursive=TRUE)
  if(length(list.out.files)>0){
    for(res.path in list.out.files){
      message(res.path)
      read.blastout <- read.delim(res.path,header=FALSE,quote="")
      read.blastout <- read.blastout %>% group_by(V1) %>% arrange(V3) %>% slice(1)
      prokka2ncbi.link <- read.blastout %>% select(V1,V2,V3)
      colnames(prokka2ncbi.link) =c('prokka_locus_tag','protein_id','evalue')

      prokka2ncbi.list[[s.gnm]] <- prokka2ncbi.link

      full_accession <- sub('.results.out','',basename(res.path))
      }
    }

    feauture.path <- list.files(path = dir,
                                 pattern='_feature_table.txt',
                                 full.names=TRUE,
                                 include.dirs = TRUE,
                                 recursive=TRUE)
    feauture.path <- feauture.path[!grepl('.gz',feauture.path)]
    if(length(feauture.path)>0){
      if(file.exists(feauture.path)){
        message('file does exist!')
        feature.out <- read.delim(feauture.path,quote='')

        #extract old locus tag from attributes column, and compile a translation table with columns: locus_tag, old_locus_tag
        lc2oldlc <- feature.out %>%
          select(locus_tag, attributes) %>%
          filter(grepl('old_locus_tag',attributes)) %>%
          mutate(old_locus_tag=gsub('old_locus_tag=','',attributes)) %>%
          select(locus_tag, old_locus_tag)


        feature.out <- feature.out %>%
          filter(X..feature=='CDS') %>%
          left_join(lc2oldlc, by='locus_tag')

        ncbifeautures.list[[s.gnm]] <- feature.out
      }
    }


    # gbff.path <- list.files(path = dir,
    #                             pattern='_genomic.gbff',
    #                             full.names=TRUE,
    #                             include.dirs = TRUE,
    #                             recursive=TRUE)
    # gbff.path <- gbff.path[!grepl('.gz',gbff.path)]
    # if(length(gbff.path>0)){
    #   if(file.exists(gbff.path)){
    #     message('reading genbank file')
    #     tryCatch({
    #       ncbi_gff = genbankr::readGenBank(gbff.path)
    #       ncbi_gff <- ncbi_gff %>%
    #         genbankr::cds() %>%
    #         data.frame() %>%
    #         tibble::tibble() %>%
    #         dplyr::select(-translation)
    #       colnames(ncbi_gff) <- paste0('ncbi_',colnames(ncbi_gff))
    #
    #       ncbi_gff <- ncbi_gff %>%
    #         filter(!is.na(ncbi_protein_id)) %>%
    #         filter(ncbi_protein_id!='') %>%
    #         group_by(ncbi_protein_id) %>%
    #         summarise(ncbi_locus_tag=paste(unique(ncbi_locus_tag),collapse=','),
    #                   ncbi_old_locus_tag=paste(unique(ncbi_old_locus_tag),collapse=','),
    #                   ncbi_product=paste(unique(ncbi_product),collapse=','),
    #                   ncbi_gene=paste(unique(ncbi_gene),collapse=','),
    #                   ncbi_EC_number=paste(unique(ncbi_EC_number),collapse=','),
    #                   ncbi_gene_synonym=paste(unique(ncbi_gene_synonym),collapse=',')
    #         ) %>%
    #         arrange(ncbi_locus_tag)
    #
    #
    #       ncbigbff.list[[s.gnm]] <- ncbi_gff
    #     },
    #     error = function(c) "error")
    #   }
    # }
}





#annotate prokka locus_tags with NCBI loci and annotations
prokka2ncbi.annotated.list <- list()
for(s.gnm in names(ncbifeautures.list)){
  if(!is.null(ncbifeautures.list[[s.gnm]])){
    if(!is.null(prokka2ncbi.list[[s.gnm]])){
      message('annotating ', s.gnm)
      prot.id.tbl <- ncbifeautures.list[[s.gnm]] %>%
        dplyr::filter(product_accession != '') %>%
        dplyr::select(product_accession,non.redundant_refseq,name, symbol, GeneID, locus_tag,old_locus_tag) %>%
        dplyr::group_by(product_accession) %>% dplyr::summarise(non.redundant_refseq=paste0(unique(non.redundant_refseq),collapse=';'),
                                                                ncbi_product=paste0(unique(name),collapse=';'),
                                                                ncbi_symbol=paste0(unique(symbol),collapse=';'),
                                                                ncbi_GeneID=paste0(unique(GeneID),collapse=';'),
                                                                ncbi_locus_tag=paste0(locus_tag,collapse=';'),
                                                                ncbi_old_locus_tag=paste0(old_locus_tag,collapse=';')) %>%
        unique() %>%
        dplyr::arrange(ncbi_locus_tag)

      prokka2ncbi.annotated <- prokka2ncbi.list[[s.gnm]] %>%
        dplyr::left_join(prot.id.tbl,by=c('protein_id'='product_accession'))

      prokka2ncbi.annotated.list[[s.gnm]] <- prokka2ncbi.annotated
    }
  }

}





#=============





#----------------------------------------------#
annotation.tbl %>%
  filter(genome=='Marinobacter_algicola_DG893') %>%
  left_join(prokka2ncbi.annotated.list[['Marinobacter_algicola_DG893']], by=c('locus_tag'='prokka_locus_tag')) %>%
  select(product,ncbi_product,locus_tag,ncbi_locus_tag,gene,ncbi_symbol)


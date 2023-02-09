#   CREATE ANOTATION TABLE
#
#-----------------------------------------------------------------------------#
#   UTILITY FUNCTIONS   ####
#-----------------------------------------------------------------------------#

# loadPorkkaGff
#=================================================================================#
#' Loading Prokka derived .gff files into a list object
#' specify the path to a prokka output folder, will collect .gffs recursively in all subfolders
#' requires library(rtracklayer)
#'
#' @param inPath path to prokka output folder
#' @param clean_genome_names boolean, when true replaces "-" and "." to "_"
#'
#' @return
#' @export
#'
#' @examples
#' prokka_gff_perGenome = loadPorkkaGff(inPath="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output/", clean_genome_names=TRUE)

loadPorkkaGff <- function(inPath="", clean_genome_names=TRUE){
  if(inPath==""){
    stop("You'd need to specify a directory using inPath")
  }else{
    prokka_gff_files <- list.files(
      path=inPath,
      pattern = ".gff$",
      recursive = TRUE,
      include.dirs=TRUE,
      full.names=TRUE)

    prokka_gff_perGenome = list()
    message("Loading files, this may take a minute")
    for(x in prokka_gff_files){
      prokka_gff = rtracklayer::import.gff(x) %>% as.data.frame()

      if(clean_genome_names==TRUE){
        genome = gsub("-","_",gsub("\\.","_", gsub(".gff","",basename(x))))
      }else{
        genome = gsub(".gff","",basename(x))
      }

      prokka_gff_perGenome[[genome]]=cbind(prokka_gff,"genome"=genome)
    }
    return(prokka_gff_perGenome)
  }
}


# load_fasta_sequences
#=================================================================================#
#' Function to load sequence data from fasta files
#' point to folder and read in any fasta file in a specified directory (recursive)
#' Function outputs a tibble formatted file with headers, sequences, and locus_tags
#' @param inPath specifies the path to the directory were the input files are located
#' @param suffix specifies the file extension of the input file (faa, fna, ffn, fasta, etc)
#' @param recursive boolian to define whether you want to look recursively in the input directory
#'
#' @return
#' @export
#'
#' @examples
#' prokka.faa = load_fasta_sequences(inPath="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output/",
#'               suffix="faa",
#'               recursive=TRUE)
#'prokka.ffn = load_fasta_sequences(inPath="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output/",
#'               suffix="ffn",
#'               recursive=TRUE)

load_fasta_sequences <- function(inPath="", suffix= "faa",recursive=TRUE){
  if(inPath==""){
    stop("You'd need to specify a directory using inPath")
  }else{
    fastaFiles <- list.files(
      path=inPath,
      pattern = paste0(".",suffix,"$"),
      recursive = recursive,
      include.dirs=TRUE,
      full.names=TRUE)
    message("Loading files, this may take a minute")

    lapply(fastaFiles, function(fastaFile){
      dfGenomeSeq <- data.frame(readFasta(fastaFile))
      #row.names(dfGenomeSeq) <- sub("^([^ ]+).*", "\\1", dfGenomeSeq$Header, perl=T)
      return(dfGenomeSeq)
    }) -> allSequences

    #Convert tidy tibble, and remove text after white space in genome header
    #could improve this by making it optional, say " " or "."
    allSequences <- allSequences %>%
            as_tibble_col() %>%
            unnest(cols=c(value)) %>%
            mutate(locus_tag=sub("^([^ ]+).*", "\\1", Header, perl=T))

    return(allSequences)
  }
}



# loadOGs
#=================================================================================#
#' Function to load OF orthogroups (OGs) outputs a (nx2) tibble with OG per locus_tags
#'
#' @param inPath specifies the path to the directory were the OrthoFinder output files are located
#'
#' @return
#' @export
#'
#' @examples
#' #OG.tbl = loadOGs(inPath = '~/DATA/MarbGenomics/Ortho/OrthoFinder/Results_Oct23_3/Orthogroups/')

loadOGs <- function(inPath) {
  OG.path <- paste0(inPath, "/Orthogroups.tsv")
  #unassigned.path <- paste0(inPath, "/Orthogroups_UnassignedGenes.tsv")

  if (!(file.exists(OG.path) )) { #& file.exists(unassigned.path)
    stop("No file(s), check that you added the right path")
  }

  dfOrthogroup <- read.table(OG.path , header=T, sep = "\t", row.names = 1, stringsAsFactors=F)
  dfsingle <-dfOrthogroup %>% unite(col,sep =', ')
  lsOrthogroup <- strsplit( as.matrix(dfsingle),", ")
  names(lsOrthogroup) <- rownames(dfOrthogroup)
  lsOrthogroup <- lapply(lsOrthogroup, function(z){ z[!is.na(z) & z != ""]})

  OG.tbl <- enframe(lsOrthogroup) %>% filter(!is.na(name))  %>% # creates the 'value' as a `list` column
    mutate(value = map(value, as.character)) %>% group_by(name) %>%
    unnest_longer(col=c(value)) %>%
    ungroup() %>%
    filter(!is.na(value)) %>%
    dplyr::rename("OG" = name) %>%
    dplyr::rename("locus_tag" = value)
  return(OG.tbl)
}


#-----------------------------------------------------------------------------#
#   RUN CODE   ####
#-----------------------------------------------------------------------------#

# BUILD ANNOTATION TABLE
#=================================================================================#


# STEP 1 - START THE ANNOTATION TABLE
#-------------------------------------#
# In this case, we start from prokka gffs, giving the location on the genome/contigs!

#Load Prokka annotations as gff file
prokka_gff_perGenome = loadPorkkaGff(inPath="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output/", clean_genome_names=TRUE)

#Start building tha annotation.tbl with Prokka ggfs as base
annotation.tbl <-prokka_gff_perGenome %>%
    tibble::as_tibble_col() %>%
    tibble::unnest(c(value))

#cleanup
rm(prokka_gff_perGenome)


#- STEP 2 - ADD OF orthogroup data!
# -------------------------------------------------------------#
# add the orthogroup as identified by OrthoFinder

#use function to load the OF from tsv
OG.tbl <- loadOGs(inPath = '~/DATA/MarbGenomics/Ortho/OrthoFinder/Results_Oct23_3/Orthogroups/')

#merge with the annotation table
annotation.tbl <- annotation.tbl  %>%
  dplyr::left_join(OG.tbl, by = "locus_tag")

#Have a look
annotation.tbl %>%
  dplyr::select(locus_tag, product, genome, OG) %>%
  dplyr::filter(OG =='OG0000464')

#cleanup
rm(OG.tbl)


#----------------------------------------------#
#---- prokka specific, COG extraction....
# COG data is stored in db_xref column

annotation.tbl <- annotation.tbl %>%
  tidyr::separate(db_xref, c('tmpcog','COG'), sep=':') %>%
  dplyr::select(!tmpcog)



# STEP 3 - Once we obtained hmm.tbls from our analsys, we can merge them to the annotation.tbl
#-------------------------------------#
#kfm.tbl
annotation.tbl <- annotation.tbl  %>%
  dplyr::left_join(kfm.tbl %>%
  dplyr::select(domain_name,query_name,sequence_evalue), by=c("locus_tag"="query_name"))%>%
  dplyr::group_by(locus_tag) %>%
  dplyr::mutate("kfm_domain" = paste0(domain_name, collapse = "; ")) %>%
  dplyr::slice(1) %>%
  dplyr::select(-domain_name) %>%
  dplyr::rename("kfm_evalue" = sequence_evalue) %>%
  dplyr::ungroup()

#cazy.tbl
annotation.tbl <- annotation.tbl  %>%
  dplyr::left_join(cazy.tbl %>% dplyr::select(domain_name,query_name,sequence_evalue), by=c("locus_tag"="query_name"))%>%
  dplyr::group_by(locus_tag) %>%
  dplyr::mutate("cazy_domain" = paste0(domain_name, collapse = "; ")) %>%
  dplyr::slice(1) %>%
  dplyr::select(-domain_name) %>%
  dplyr::rename("cazy_evalue" = sequence_evalue) %>%
  dplyr::ungroup()

#tcdb.tbl
annotation.tbl <- annotation.tbl  %>%
  dplyr::left_join(tcdb.tbl %>%
  dplyr::select(domain_name,query_name,sequence_evalue), by=c("locus_tag"="query_name"))%>%
  dplyr::group_by(locus_tag) %>%
  dplyr::mutate("tcdb_domain" = paste0(domain_name, collapse = "; ")) %>%
  dplyr::slice(1) %>%
  dplyr::select(-domain_name) %>%
  dplyr::rename("tcdb_evalue" = sequence_evalue) %>%
  dplyr::ungroup()


#===================================================================≠≠≠#
# Some usage examples
#Have a look at what we just added look right

annotation.tbl %>% select(locus_tag,
                          genome,
                          product,
                          Name,
                          gene,
                          eC_number,
                          kfm_domain,
                          cazy_domain,
                          tcdb_domain) %>%
                    filter(cazy_domain!="NA") %>%
  filter(genome=="Marinobacter_sp_FDB33")


#Filter using GREP
annotation.tbl %>% select(locus_tag, genome, product)  %>%
  filter(grepl("16S ribosomal",product))

annotation.tbl %>% select(locus_tag, genome, product,cazy_domain)  %>%
  filter(grepl("PL",cazy_domain)) %>% filter(genome=="Marinobacter_sp_FDB33")






#- STEP 4 = ADD SEQUENCES - load_fasta_sequences()
# -------------------------------------------------------------#

#Load the sequences from Prokka files
prokka.faa = load_fasta_sequences(inPath="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output/",
                                  suffix="faa")
prokka.ffn = load_fasta_sequences(inPath="~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/prokka_output/",
                                  suffix="ffn")

# add the amino-acid sequences to the annotation table
annotation.tbl <- annotation.tbl  %>%
  dplyr::left_join(prokka.faa %>% select(locus_tag,Sequence), by=c("locus_tag"="locus_tag"))%>%
  dplyr::rename(aa_sequence = Sequence)

# add the nucleotide sequences sequences to the annotation table
annotation.tbl <- annotation.tbl  %>%
  dplyr::left_join(prokka.ffn %>% dplyr::select(locus_tag,Sequence), by=c("locus_tag"="locus_tag"))%>%
  dplyr::rename(nuc_sequence = Sequence)

#Have a look at the sequences!
annotation.tbl %>% select(locus_tag,aa_sequence,nuc_sequence)


#cleanup
rm(prokka.faa)
rm(prokka.ffn)



#SQL database####

#-------------------------------------------------#
# Add to database!
# GIven the size of this tibble, it is sensible to use SQL
# writing and wrangling this around is much faster when using the tblAnnotation instead of annotation.tbl
# especially because we loaded in the sequences as well
#--- In the marinobacter example the database is >800MB big....

#Think about inference, as this is stored as a list in the tibble, dbWriteTable hates this
#annotation.tbl %>% select(-inference)

#make connection
con_sqlite <- DBI::dbConnect(
  RSQLite::SQLite(),
  "MarinobacterGenome.db"
)

#write tibble to tblAnnotation
DBI::dbWriteTable(con_sqlite, "tblAnnotation", annotation.tbl %>% select(-inference) , overwrite = TRUE, row.names=F)

#------------#
# load it

pool <- pool::dbPool(RSQLite::SQLite(), dbname = "MarinobacterGenome.db")

#assign variables to pool connections
#tblAssembly <- tbl(pool, "tblAssembly")#%>% unique()
#tblBioSample <- tbl(pool, "tblBioSample")#%>% unique()
#tblBioProject  <- tbl(pool, "tblBioProject")#%>% unique()

tblAnnotation  <- tbl(pool, "tblAnnotation")#%>% unique()
tblAnnotation

annotation.tbl <- tblAnnotation %>% dplyr::select(-aa_sequence,-nuc_sequence) %>% collect()

#-----------------------------------------------------------------------------#
#   TIME CONNECTION    ####
#-----------------------------------------------------------------------------#

# Speed testing for reading from the +800MB object annotation.tbl vs. SQL
#=================================================================================#

starttime = Sys.time()
tblAnnotation %>% select(locus_tag, genome, product,cazy_domain) %>%
  filter(genome=="Marinobacter_sp_FDB33") %>%
  collect() %>%
  filter(grepl("PL",cazy_domain))
runtime1 = Sys.time() - starttime
runtime1

#----------------------------------------#

starttime = Sys.time()
annotation.tbl %>% select(locus_tag, genome, product,cazy_domain) %>%
  filter(genome=="Marinobacter_sp_FDB33") %>%
  collect() %>%
  filter(grepl("PL",cazy_domain))
runtime2 = Sys.time() - starttime
runtime2


#-------------------------------------------------------------#
#
# Bravo, now indeed we can create the perfect sequence extractor!
#
# -------------------------------------------------------------#

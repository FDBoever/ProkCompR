
# LOAD GENOMIC ONES, AND SPLIT INTO BLOCKS OF 10 FOR IMNGS
inPath.contigs <- '/Users/sa01fd/DATA/MarbGenomics/16Sseq/all16S.fasta'
contigs <- Biostrings::readDNAStringSet(inPath.contigs)
n_contigs <- length(contigs)
contigs.tbl <- tibble::tibble(name=names(contigs),
                             sequence=paste(contigs),
                             GC = c(Biostrings::letterFrequency(contigs, letters = "GC", as.prob = TRUE)),
                             length = nchar(sequence),
                             length_threshold = ifelse(length>length_threshold-1, 'long','short'))


contigs.tbl %>%  dplyr::slice(1:10) %>% dplyr::select(name)

writeFasta(data= contigs.tbl %>% slice(1:10) %>% select(name, sequence), filename = 'mb_imngs_1.fasta')
writeFasta(data= contigs.tbl %>% slice(11:20) %>% select(name, sequence), filename = 'mb_imngs_2.fasta')
writeFasta(data= contigs.tbl %>% slice(21:30) %>% select(name, sequence), filename = 'mb_imngs_3.fasta')
writeFasta(data= contigs.tbl %>% slice(31:40) %>% select(name, sequence), filename = 'mb_imngs_4.fasta')
writeFasta(data= contigs.tbl %>% slice(41:50) %>% select(name, sequence), filename = 'mb_imngs_5.fasta')
writeFasta(data= contigs.tbl %>% slice(51:60) %>% select(name, sequence), filename = 'mb_imngs_6.fasta')
writeFasta(data= contigs.tbl %>% slice(61:62) %>% select(name, sequence), filename = 'mb_imngs_6.fasta')




cat accession_numbers.txt | xargs -n 1 sh -c 'efetch -db nuccore -id $0 -format gbc | xtract -pattern INSDSeq -tab "\t" -sep "; " -def "NA" -element INSDSeq_locus INSDSeq_length INSDSeq_strandedness INSDSeq_moltype INSDSeq_topology INSDSeq_division INSDSeq_update-date INSDSeq_create-date INSDSeq_definition INSDSeq_primary-accession INSDSeq_accession-version INSDSeq_source INSDSeq_organism INSDSeq_taxonomy INSDSeq_sequence INSDReference_title INSDXref_dbname INSDXref_id INSDReference_pubmed INSDReference_journal' >  m16.insd.metadata.txt
cat accession_numbers.txt | xargs -n 1 sh -c 'efetch -db nuccore -id $0 -format gbc | xtract -insd source organism mol_type isolation_source db_xref isolate host culture_collection strain clone environmental_sample country lat_lon collection_date collected_by type_material PCR_primers' >  m16.source.metadata.txt




#-------------------------------------------------------
# PREPARATION STEP TO FILTER OUT SILVA SEQEUCNES
inPath.contigs <- '/Users/sa01fd/DATA/MarbGenomics/16Sseq/including_silva.fasta'
length_threshold <- 1300

message(':: loading fasta input')
contigs <- Biostrings::readDNAStringSet(inPath.contigs)
n_contigs <- length(contigs)
message('  - nr of contigs: ',n_contigs)

message(':: compiling tibble')
contigs.tbl <- tibble::tibble(name=names(contigs),
                              sequence=paste(contigs),
                              GC = c(Biostrings::letterFrequency(contigs, letters = "GC", as.prob = TRUE)),
                              length = nchar(sequence),
                              length_threshold = ifelse(length>length_threshold-1, 'long','short'))


contigs.filtered <- contigs.tbl %>%
  dplyr::filter(length>1300) %>% # filter out any sequence below 1200 nt
  dplyr::filter(!grepl('N',sequence)) %>% # filte any sequence with an N
  dplyr::filter(!grepl('Marinobacterium',name)) %>% # filter out marinobacterium sequences
  dplyr::arrange(desc(length)) #sorted on length needed/prefered by USEARCH




writeFasta(data= contigs.filtered %>% select(name,sequence), filename = 'filtered_16s.fasta')


# DECIPHER STEP
inPath.decipher <- '/Users/sa01fd/DATA/MarbGenomics/16Sseq/DECIPHER_out.txt'

decipher <- read.delim(inPath.decipher)

head(decipher)
decipher %>%
  dplyr::group_by(Result) %>%
  dplyr::tally()

decipher %>%
  dplyr::filter(Result=='Not deciphered to be a chimera') %>%
  dplyr::group_by(Group)%>%
  dplyr::tally() %>%

filtered.decipher <- decipher %>%
  dplyr::filter(Result=='Not deciphered to be a chimera') %>%
  dplyr::filter(Group %in% c('Marinobacter','Marinobacter_hydrocarbonoclasticus','Marinobacter_salarius','Oceanospirillaceae','Oceanospirillales','Hahellaceae')) %>%
  dplyr::select(Name) %>%
  #dplyr::mutate(Name=gsub(" .*", "",Name)) %>%
  dplyr::pull()

contigs.deciphered <-  contigs.filtered %>% dplyr::filter(name %in% filtered.decipher)

writeFasta(data= contigs.deciphered %>%
             dplyr::select(name,sequence), filename = 'filtered_decipher_16s.fasta')



#
# PLAYING TO LOAD THE ALIGNMENT, AND SEE IF I CAN EXTYACT ACCESSION NUMBERS TO BE QUERIED AGAINST NCBI
#===========================


inPath.contigs <- '/Users/sa01fd/DATA/MarbGenomics/16Sseq/16S_qual_SINA_TrimAl.fasta'

message(':: loading fasta input')
contigs <- Biostrings::readDNAStringSet(inPath.contigs)
n_contigs <- length(contigs)
message('  - nr of contigs: ',n_contigs)

message(':: compiling tibble')
contigs.tbl <- tibble::tibble(name=names(contigs),
                              sequence=paste(contigs),
                              GC = c(Biostrings::letterFrequency(contigs, letters = "GC", as.prob = TRUE)),
                              length = nchar(sequence),
                              length_threshold = ifelse(length>length_threshold-1, 'long','short'))


contigs.filtered <- contigs.tbl %>% filter(length>1300) %>% # filter out any sequence below 1200 nt
  filter(!grepl('N',sequence)) %>% # filte any sequence with an N
  filter(!grepl('Marinobacterium',name)) %>% # filter out marinobacterium sequences
  arrange(desc(length)) #sorted on length needed/prefered by USEARCH



accns <- contigs.filtered %>% mutate(accession=sub("\\..*", "", name)) %>%
  filter(!grepl("Marinobacter|Halospina|Oleiphilus|Hahella|Oleispira|Thalassolituus|Mangrovitalea",accession)) %>%
  select(accession) %>% pull()





write(accns, "accession_numbers.txt")


in.tree <- '/Users/sa01fd/DATA/MarbGenomics/16Sseq/16S_qual_SINA_TrimAl_deciphered.fasta.treefile'

tree <- ape::read.tree(in.tree)
tree.rooted <- phytools::midpoint.root(tree)

drops <- tree.rooted$tip.label[grepl('HG004177|EU907917|KR088715|HG004182|KR088715|GQ250093|FJ973518|KF201610|KC573502',tree.rooted$tip.label)]
tree.rooted <- ape::drop.tip(tree.rooted,tip=drops)


outgroup <- tree.rooted$tip.label[grepl('Oleispir|Thalassolituu|Oleiphilu|Hahell|Mangrovital',tree.rooted$tip.label)]


phytools::findMRCA(tree.rooted, tips =outgroup)
tree.rooted <- phytools::reroot(tree.rooted, 713)

p.1 <- ggtree::ggtree(tree.rooted,layout='circular') + ggtree::geom_tiplab(align = TRUE,size=3,aes(color=ifelse(grepl('Marinobacter',label),'m','p')))
ggtree::open_tree(p.1,angle=180)

p.2 <- ggtree::ggtree(tree.rooted) + ggtree::geom_tiplab(align = TRUE,size=2,aes(color=ifelse(grepl('Marinobacter',label),'m','p')))
ggsave('tree.16S.initial.pdf',p.2,width = 15,height = 30)
ncbi.meta <- read.delim('/Users/sa01fd/DATA/MarbGenomics/16Sseq/m16.source.metadata.txt',header=F) %>% tibble::tibble()


p.1 <- ggtree::ggtree(tree.rooted,layout='circular') + ggtree::geom_tiplab(align = TRUE,size=3,aes(color=ifelse(grepl('Marinobacter',label),'m','p')))
ggtree::open_tree(p.1,angle=180)


tree.dropped.outgroup <- ape::drop.tip(tree.rooted,tip=outgroup)
p.1 <- ggtree::ggtree(ape::ladderize(phytools::midpoint.root(tree.dropped.outgroup)),layout='circular') + ggtree::geom_tiplab(align = TRUE,size=3,aes(color=ifelse(grepl('Marinobacter',label),'m','p')))
ggtree::open_tree(p.1,angle=180)



ape::ladderize(tree.dropped.outgroup)

#reading datafiles
NCBImetadata_inds = read.delim('/Users/sa01fd/DATA/MarbGenomics/16Sseq/m16.insd.metadata.txt',sep='\t',header = FALSE,stringsAsFactors = FALSE)
NCBImetadata_source = read.delim('/Users/sa01fd/DATA/MarbGenomics/16Sseq/m16.source.metadata.txt',sep='\t',header = FALSE ,stringsAsFactors = FALSE)

#-------------------------------
#	Preparing before merging
#--------------------------------
#NCBImetadata_inds
colnames(NCBImetadata_inds)=c("locus","length","strandedness","moltype","topology","division","update_date","create_date","definition","primary_accession","accession_version","source","organism","taxonomy","sequence","Reference_title","Xref_dbname","Xref_id","Reference_pubmed","Reference_journal")
head(NCBImetadata_inds)

NCBImetadata_inds = NCBImetadata_inds %>% unique()
#rownames(NCBImetadata_inds) = NCBImetadata_inds[,1]
NCBImetadata_inds = NCBImetadata_inds %>% tibble::as_tibble()


#--------------------------------
#NCBImetadata_source
#organism mol_type isolation_source db_xref isolate culture_collection strain clone environmental_sample country lat_lon collection_date collected_by type_material PCR_primers
colnames(NCBImetadata_source)=c("accession_version","organism","mol_type","isolation_source","db_xref", 'isolate','host','culture_collection','strain',"clone","environmental_sample","country","lat_lon","collection_date","collected_by","type_material",'PCR_primers')
head(NCBImetadata_source)

NCBImetadata_source = NCBImetadata_source %>% unique()
#rownames(NCBImetadata_source) = NCBImetadata_source[,1]
NCBImetadata_source = NCBImetadata_source %>% tibble::as_tibble()


#-------------------------------
#	Merging tibbles
#--------------------------------
mergedNCBImetadata = NCBImetadata_inds %>% dplyr::left_join(NCBImetadata_source %>%
                                                              dplyr::select(-organism), by = c("accession_version" = "accession_version"))


mergedNCBImetadata %>%
  dplyr::filter(accession_version %in% sub('^([^.]+.[^.]+).*', '\\1', tree.rooted$tip.label))  %>%
  dplyr::select(accession_version, organism, type_material)


p.1 <- ggtree::ggtree(ape::ladderize((tree.dropped.outgroup)),layout='circular')

p.1$ data <- p.1$data %>%
  dplyr::mutate(name = label) %>%
  tidyr::separate(name, c("A","B",'C'), sep = "([.?:])") %>%
  dplyr::mutate(B=as.numeric(B)) %>%
  dplyr::mutate(g = ifelse(B>1,'y','n')) %>% #dplyr::filter(B>1) %>%
  dplyr::left_join(mergedNCBImetadata, by=c('A'='locus')) %>%
  #dplyr::select(label, organism, isolation_source, strain) %>% data.frame()
  dplyr::mutate(label=paste0(label,'-',organism))


p.4 <- p.1 + geom_tiplab(size=3, aes(color=g),align=T)
ggsave('tree.16S.initial_t.pdf',p.4,width = 30,height = 30)

mergedNCBImetadata




mergedNCBImetadata %>%
  select(accession_version, organism, type_material, lat_lon, isolation_source) %>%
  filter(isolation_source !='-') %>%
  group_by(isolation_source) %>%
  filter(n>1) %>%
  tally() %>% ggplot(aes(reorder(isolation_source,n),n)) + geom_col() + coord_flip()

mergedNCBImetadata %>%
  select(accession_version, organism, type_material, lat_lon, isolation_source) %>%
  filter(lat_lon !='-')



##==========DOWNLOADS FROM ARB
inPath.contigs <- '/Users/sa01fd/Downloads/arb-silva.de_2021-08-26_id1046107_tax_silva.fasta'
#
length_threshold <- 1299

#loading fasta into DNAStringSet
message(':: loading fasta input')
contigs <- Biostrings::readDNAStringSet(inPath.contigs)
n_contigs <- length(contigs)
message('  - nr of contigs: ',n_contigs)

message(':: compiling tibble')
contigs.tbl <- tibble::tibble(name=names(contigs),
                              sequence=paste(contigs),
                              GC = c(Biostrings::letterFrequency(contigs, letters = "GC", as.prob = TRUE)),
                              length = nchar(sequence),
                              length_threshold = ifelse(length>length_threshold-1, 'long','short'))


contigs.filtered <- contigs.tbl %>% filter(length>1200) %>% # filter out any sequence below 1200 nt
  filter(!grepl('N',sequence)) %>% # filte any sequence with an N
  filter(!grepl('Marinobacterium',name)) %>% # filter out marinobacterium sequences
  arrange(desc(length)) #sorted on length needed/prefered by USEARCH

writeFasta(data= contigs.filtered %>% select(name,sequence), filename = 'test_marinobacter_silva.fasta')

contigs.filtered %>% ggplot(aes(length, GC)) + geom_point()
contigs.filtered %>% ggplot(aes(length, GC)) + geom_point()


contigs.filtered %>% arrange(desc(GC))

contigs.tbl %>% ggplot(aes(length)) + geom_histogram()





#============================================





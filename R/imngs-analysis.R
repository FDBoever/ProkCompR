#install.packages('parzer')

library(dplyr)

#IMNGS
#=====================================================================
inPath.contigs <- '/Users/sa01fd/DATA/MarbGenomics/16Sseq/all16S.fasta'
#inPath.contigs <- '/Users/frederik/DATA/IMNGS/all16S.fasta'


length_threshold <- 1200

contigs <- Biostrings::readDNAStringSet(inPath.contigs)
n_contigs <- length(contigs)
contigs.tbl <- tibble::tibble(name=names(contigs),
                              sequence=paste(contigs),
                              GC = c(Biostrings::letterFrequency(contigs, letters = "GC", as.prob = TRUE)),
                              length = nchar(sequence),
                              length_threshold = ifelse(length>length_threshold-1, 'long','short'))


mb_imngs_1.names <- contigs.tbl %>%
  dplyr::slice(1:10) %>%
  dplyr::select(name) %>%
  dplyr::pull()

mb_imngs_2.names <- contigs.tbl %>%
  dplyr::slice(11:20) %>%
  dplyr::select(name) %>%
  dplyr::pull()

mb_imngs_3.names <- contigs.tbl %>%
  dplyr::slice(21:30) %>%
  dplyr::select(name) %>%
  dplyr::pull()

mb_imngs_4.names <- contigs.tbl %>%
  dplyr::slice(31:40) %>%
  dplyr::select(name) %>%
  dplyr::pull()

mb_imngs_5.names <- contigs.tbl %>%
  dplyr::slice(41:50) %>%
  dplyr::select(name) %>%
  dplyr::pull()

mb_imngs_6.names <- contigs.tbl %>%
  dplyr::slice(51:60) %>%
  dplyr::select(name) %>%
  dplyr::pull()

mb_imngs_7.names <- contigs.tbl %>%
  dplyr::slice(61:62) %>%
  dplyr::select(name) %>%
  dplyr::pull()


#=====================================================================
# value column in percent, for example, 0.1 is 0.1% of total number of sequences in that sample

#SRA-derroved sample is considered positive if the query-like sequences sum up to more than x% of the total number of sequences in that sample (abundance >0.01%)
colnames(imngs.sample_table[,c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10')])

inPath.imngs <- '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3264/'
inPath.imngs <- '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3265/'
inPath.imngs <- '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3266/'

inPath.imngs <- '/Users/frederik/DATA/IMNGS/paraller_sim_query_3267/'


inPath.imngs.m <- c(
          mb_imngs_1 = '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3264/',
					mb_imngs_2 = '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3265/',
					mb_imngs_3 = '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3266/',
					mb_imngs_4 = '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3267/',
					mb_imngs_5 = '/Users/sa01fd/DATA/MarbGenomics/IMNGS/paraller_sim_query_3281/')

names.imngs.m <- list(mb_imngs_1 = mb_imngs_1.names,
					mb_imngs_2 = mb_imngs_2.names,
					mb_imngs_3 = mb_imngs_3.names,
					mb_imngs_4 = mb_imngs_4.names,
					mb_imngs_5 = mb_imngs_5.names)

imngs.RA.l <- c()

for(i in names(inPath.imngs.m)){
	message(i)
	imngs.sample_table <- read.delim(paste0(unname(inPath.imngs.m[i]),'IMNGS_99/IMNGS_99_counts_matrix.tab'))
	colnames(imngs.sample_table)[4:13] <- names.imngs.m[[i]]

	imngs.RA <- read.delim(paste0(unname(inPath.imngs.m[i]),'IMNGS_99/IMNGS_99_abun_0_matrix.tab'))
	colnames(imngs.RA)[2:11] <- names.imngs.m[[i]]

	imngs.RA.tmp <- imngs.RA %>% tidyr::pivot_longer(-Sample_Name)%>%
	dplyr::left_join(imngs.sample_table %>%
                     dplyr::select(Sample_Name, Description), by='Sample_Name')

    imngs.RA.l <- rbind(imngs.RA.l, imngs.RA.tmp)

}




head(imngs.RA.l)


imngs.RA.l %>%
  dplyr::filter(value>0.01) %>%
  ggplot2::ggplot(ggplot2::aes(x=name,y=value, color = Description)) +
  ggplot2::geom_jitter() +
  ggplot2::coord_flip() +
  fdb_style(aspect.ratio=2)

imngs.RA.l %>%
  dplyr::filter(value>0.01) %>%
  dplyr::group_by(name,Description) %>%
  dplyr::summarise(sum=sum(value))

imngs.RA.l %>%
  dplyr::filter(value>0.01) %>%
  dplyr::group_by(name,Description) %>%
  dplyr::summarise(sum=sum(value)) %>%
  ggplot2::ggplot(ggplot2::aes(x=name,sum,fill=Description)) +
  ggplot2::geom_bar(position="stack", stat="identity")+
  ggplot2::coord_flip() +
  fdb_style(aspect.ratio=2)


imngs.RA.l %>%
  dplyr::filter(value>0.01) %>%
  dplyr::group_by(name,Description) %>%
  dplyr::summarise(sum=sum(value)) %>%
  ggplot2::ggplot(ggplot2::aes(x=name,sum,fill=Description)) +
  ggplot2::geom_bar(position="fill", stat="identity")+
  ggplot2::coord_flip() +
  fdb_style(aspect.ratio=2)

imngs.RA.l %>%
  dplyr::filter(value>0.01) %>%
  dplyr::group_by(name,Description) %>%
  dplyr::summarise(sum=sum(value)) %>%
  ggplot2::ggplot(ggplot2::aes(x=Description,sum,fill=Description)) +
  ggplot2::geom_bar(position="fill", stat="identity") +
  #fdb_style(aspect.ratio=.5)+
  ggplot2::facet_wrap(~name,ncol=1,scales='free')+
  ggplot2::coord_flip()

imngs.RA.l %>%
  dplyr::filter(value>0.01) %>%
  dplyr::filter(name %in% c('Marinobacter_algicola_DG893','Marinobacter_adhaerens_HP15','Marinobacter_hydrocarbonoclasticus_ATCC_49840',"Marinobacter_psychrophilus_20041","Marinobacter_sp_BSs20148","Marinobacter_sp_ELB17" )) %>%
  dplyr::group_by(name,Description) %>%
  dplyr::summarise(sum=sum(value)) %>%
  ggplot2::ggplot(ggplot2::aes(x=Description,sum,fill=Description)) +
  ggplot2::geom_bar(stat="identity") +
  #fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~name,ncol=2,scales='free')+
  ggplot2::coord_flip()


imngs.RA.l %>%
  dplyr::left_join(imngs.sample_table %>%
                     dplyr::select(Sample_Name, Description), by='Sample_Name') %>%
  dplyr::filter(value>0.01) %>%
  dplyr::filter(name %in% c('Marinobacter_algicola_DG893','Marinobacter_adhaerens_HP15')) %>%
  dplyr::group_by(name,Description) %>%
  #dplyr::summarise(sum=sum(value)) %>%
  ggplot(aes(x=Description,value,color=Description)) +
  geom_point() +
  fdb_style(aspect.ratio=2)+facet_wrap(~name,ncol=2,scales='free')+coord_flip()









imngs.RA.l %>%
  dplyr::left_join(imngs.sample_table %>%
                     dplyr::select(Sample_Name, Description), by='Sample_Name') %>% dplyr::arrange(desc(value)) %>%   dplyr::filter(name %in% c('Marinobacter_algicola_DG893')) %>%
  dplyr::filter(value>0.01)




#sa01fd@salmonella:/media/data/SANDBOXES/sa01fd/mb16S/imngs$ nohup bash run_screen_sra.txt > nohup.sra.lonlat &
#
# sa01fd@salmonella:/media/data/SANDBOXES/sa01fd/mb16S/imngs$ nohup bash run_screen_sra.txt > nohup.sra.lonlat &

inPath.lonlat <- '/Users/sa01fd/DATA/MarbGenomics/IMNGS/sra_lon_lat_dave_out.tsv'
inPath.lonlat <- '/Users/frederik/DATA/IMNGS/sra_lon_lat_dave_out.tsv'

sra.lonlat <- read.delim(inPath.lonlat,header=TRUE,sep='\t')
dim(sra.lonlat)
head(sra.lonlat)

# For those that are separated by a _
latlon1 <- sra.lonlat %>%
  dplyr::filter(LAT_LONG !='') %>%
  dplyr::filter(grepl('_',LAT_LONG)) %>%
  tidyr::separate(col=LAT_LONG,into=c('lat','lon'),sep='_') %>%
  dplyr::mutate(lat=parzer::parse_lat(lat),
                lon=parzer::parse_lon(lon))%>%
  dplyr::select(SRA, lat, lon,Title, Study_Title,Bioproject)
head(latlon1)

#for those that are annotated as N/S E/W
latlon2 <- sra.lonlat %>% tibble::tibble() %>%
  dplyr::filter(LAT_LONG !='') %>%
  dplyr::filter(!grepl('_',LAT_LONG)) %>%
  tidyr::separate(col=LAT_LONG,into=c('lat','lon'),sep='N|S',remove=FALSE) %>%
  tidyr::separate(col=lon,into=c('lon','londir'),sep='E|W',remove=FALSE) %>%
  dplyr::mutate(lat=as.numeric(lat),
                lon=as.numeric(lon)) %>%
    dplyr::mutate(lat = ifelse(grepl('N',LAT_LONG),lat,-lat)) %>%
  dplyr::mutate(lon = ifelse(grepl('E',LAT_LONG),lon,-lon)) %>%
  dplyr::select(SRA, lat, lon,Title, Study_Title,Bioproject)
head(latlon2)





sra.coords <- rbind(latlon1, latlon2)






#================================================#

imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  	dplyr::filter(value>0.01) %>%
	dplyr::filter(!is.na(Bioproject)) %>%
	dplyr::group_by(Bioproject) %>% dplyr::tally() %>%
	dplyr::arrange(desc(n)) %>%
	ggplot2::ggplot(ggplot2::aes(x=reorder(Bioproject,n ),y=n)) +
	ggplot2::geom_bar(stat="identity") +
	ggplot2::coord_flip()

imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  	dplyr::filter(value>0.01) %>%
	dplyr::filter(!is.na(Sample_Name)) %>%
	dplyr::group_by(Sample_Name) %>% dplyr::tally() %>%
	dplyr::arrange(desc(n)) %>%
	ggplot2::ggplot(ggplot2::aes(x=reorder(Sample_Name,n ),y=n)) +
	ggplot2::geom_bar(stat="identity") +
	ggplot2::coord_flip()


imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  	dplyr::filter(value>0.01) %>%
	dplyr::filter(!is.na(Bioproject)) %>%
	dplyr::filter(!is.na(Sample_Name)) %>%
	dplyr::group_by(Sample_Name) %>% dplyr::tally() %>%
	dplyr::arrange(desc(n))

imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
    	dplyr::filter(value>0.01) %>%
  	dplyr::filter(Sample_Name == 'SRR2041075')

imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
    	dplyr::filter(value>0.01) %>%
  	dplyr::filter(grepl('arctic',Study_Title))


imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  	dplyr::filter(value>0.01) %>%
	dplyr::filter(!is.na(Title)) %>%
	dplyr::filter(!is.na(Title)) %>%
	dplyr::group_by(Title) %>% dplyr::tally() %>%
	dplyr::arrange(desc(n))%>%
  	dplyr::filter(grepl('arctic',Study_Title))

#================================================#

  imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  	dplyr::filter(grepl('amplicon',Title, ignore.case=TRUE))



#library(ggmap)
#library(ggplot2)

mapWorld <- ggplot2::borders("world", colour="black", fill="grey95", size = 0.1) # create a layer of borders

p = ggplot2::ggplot() +
  mapWorld +
  ggplot2::geom_point(data=sra.coords , ggplot2::aes(x=lon, y=lat), color='cadetblue3',size=2,alpha=0.3) +
  ggplot2::geom_point(data=sra.coords , ggplot2::aes(x=lon, y=lat), color='dodgerblue4',size=2,shape=21) +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::theme_bw() +
  ggplot2::coord_quickmap() +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.05, linetype = "dashed", colour = "grey50"))
p




imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  dplyr::filter(name %in% c('Marinobacter_algicola_DG893','Marinobacter_adhaerens_HP15','Marinobacter_hydrocarbonoclasticus_ATCC_49840',"Marinobacter_psychrophilus_20041","Marinobacter_sp_BSs20148","Marinobacter_sp_ELB17","Marinobacter_sp_Arc7_DN_1","Marinobacter_sp_AC-23","Marinobacter_lutaoensis_T5054")) %>%
  dplyr::filter(value>0.001) %>%
  ggplot2::ggplot() +
  mapWorld +
  ggplot2::geom_point(ggplot2::aes(x=lon, y=lat, color=name,size=value),alpha=0.3) +
  ggplot2::geom_point(ggplot2::aes(x=lon, y=lat, color=name,size=value),shape=21) +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::coord_quickmap() +
  ggplot2::facet_wrap(~name)




imngs.RA.l %>%
  dplyr::left_join(sra.coords, by=c('Sample_Name'='SRA')) %>%
  #dplyr::filter(name %in% c('Marinobacter_adhaerens_HP15','Marinobacter_algicola_DG893')) %>%
  dplyr::filter(value>0.001) %>%
  ggplot2::ggplot() +
  mapWorld +
  ggplot2::geom_point(ggplot2::aes(x=lon, y=lat, color=name,size=value),alpha=0.3) +
  ggplot2::geom_point(ggplot2::aes(x=lon, y=lat, color=name,size=value),shape=21) +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::coord_quickmap() #+ ggplot2::facet_wrap(~name)









positive_sra <- imngs.RA.l %>% dplyr::filter(value!=0) %>% dplyr::select(Sample_Name) %>% unique() %>% dplyr::pull() %>% as.character()
write(positive_sra,file='sra_ids.txt')




imngs <- read.delim('/Users/sa01fd/Downloads/paraller_sim_query_3262/IMNGS_99/IMNGS_99_counts_matrix.tab')
imngs %>% ggplot() + geom_point(aes(Size, X1, color=Description))
imngs %>% ggplot() + geom_point(aes(1, X1, color=Description))
imngs %>% filter(X1>1000)%>% ggplot() + geom_histogram(aes(X1))
imngs %>% filter(X1>3)%>% ggplot()  + geom_point(aes(Description, X1)) + coord_flip()


#x1 quiet abundance in algae metagenome
imngs.summed <- imngs %>% dplyr::group_by(Description) %>% dplyr::summarise(X1 = sum(X1),
                                                              X2 = sum(X2),
                                                              X3 = sum(X3),
                                                              X4 = sum(X4),
                                                              X5 = sum(X5),
                                                              X6 = sum(X6),
                                                              X7 = sum(X7),
                                                              X8 = sum(X8),
                                                              X9 = sum(X9),
                                                              X10 = sum(X10)) %>%

  arrange(desc(X4))

imngs.filtered <- imngs.summed %>%
  filter(!Description %in% c('Bacteria','metagenome','milk metagenome')) %>%
  mutate(summed = rowSums(across(X1:X10))) %>%
  arrange(desc(summed)) %>%
  head(20)

imngs.filtered

imngs.matrix <- as.matrix(imngs.filtered[,2:11])
rownames(imngs.matrix) <- as.character(imngs.filtered$Description)
heatmap.2(imngs.matrix,trace='none',scale='row',margins=c(10,10))
heatmap.2(imngs.matrix,trace='none',margins=c(10,10))


'HP15','DG893','ATCC_49840','SM19','D15_8W','ps20041','SMR5','BSs20148','FB1','T5054'


imngs %>% group_by(Description) %>% summarise(X2 = sum(X2)) %>% arrange(desc(X2))
imngs %>% group_by(Description) %>% summarise(X3 = sum(X3)) %>% arrange(desc(X3))
imngs %>% group_by(Description) %>% summarise(X3 = sum(X3)) %>% arrange(desc(X3))






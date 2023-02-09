if(!require(tibble)) install.packages("tibble", repos = "http://cran.us.r-project.org",dependencies=T)
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org",dependencies=T)
if(!require(RMySQL)) install.packages("RMySQL", repos = "http://cran.us.r-project.org",dependencies=T)
if(!require(RSQLite)) install.packages("RSQLite", repos = "http://cran.us.r-project.org",dependencies=T)

#setwd('~/DATA/MarbGenomics')


con_sqlite <- dbConnect(
		RSQLite::SQLite(),
		"MarinobacterGenome.db"
	)


#---------------------------

df = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/NCBI_assembly_metadata.txt',sep='\t',header = FALSE,stringsAsFactors = FALSE)
colnames(df) = c('RsUid','GbUid','AssemblyAccession','LastMajorReleaseAccession','LatestAccession','ChainId','AssemblyName','UCSCName','EnsemblName','Taxid','Organism','SpeciesTaxid','AssemblyType','AssemblyStatus','AssemblyStatusSort','WGS','BioprojectAccn','BioprojectId','BioSampleAccn','Sub_type','Sub_value','BioSampleId','InfraspeciesList','Isolate','Coverage','AsmReleaseDate_GenBank','AsmReleaseDate_RefSeq','SubmissionDate','SubmitterOrganization','RefSeq_category','FromType','Genbank','RefSeq','ContigN50','ScaffoldN50','FtpPath_GenBank','FtpPath_RefSeq','FtpPath_Assembly_rpt','FtpPath_Stats_rpt','SortOrder')
df = df[,c('AssemblyAccession',colnames(df)[colnames(df) != 'AssemblyAccession'])]
dbWriteTable(con_sqlite, "tblAssembly", df, overwrite = TRUE, row.names=F)

df = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/NCBI_bioproject_metadata.txt',sep='\t',header = FALSE,stringsAsFactors = FALSE)
colnames(df) = c('Id','TaxId','Project_Acc','Project_Type','Project_Data_Type','Sort_By_ProjectType','Sort_By_DataType','Sort_By_Organism','Project_Subtype','Project_Target_Scope','Project_Target_Material','Project_ObjectivesType','Registration_Date','Project_Name','Project_Title','Project_Description','Keyword','Relevance_Agricultural','Relevance_Medical','Relevance_Environmental','Relevance_Evolution','Relevance_Model','Relevance_Other','Organism_Name','Organism_Label','Sequencing_Status','string','Supergroup')
df = df[,c('Project_Acc',colnames(df)[colnames(df) != 'Project_Acc'])]
dbWriteTable(con_sqlite, "tblBioProject", df, overwrite = TRUE, row.names=F)

df1 = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/NCBI_biosample_metadata.txt',sep='\t',header = FALSE,stringsAsFactors = FALSE)
colnames(df1) = c('Id','Title','Accession','Date','PublicationDate','ModificationDate','Organization','Taxonomy','Organism','SourceSample','SampleData','Identifiers','Infraspecies','Package','SortKey')
df2 = read.delim('~/DATA/OrthoFinderMarinobacter/RecentAnalysis/Annotations2/NCBI_biosample_metadata_attributes.txt',sep='\t',header = FALSE,stringsAsFactors = FALSE)
colnames(df2) = c('Id','Sample_name','CFSAN','SRA','strain','isolate','isolate_name_alias','serovar','serotype','collection_date','lat_lon','geo_loc_name','host','isolation_source','attribute_package','IFSAC_Category','FoodOn_Ontology_Term','collected_by','bioprojectId')
df = df1 %>% left_join(df2,by=c('Accession'='CFSAN'))
df = df[,c('Accession',colnames(df)[colnames(df) != 'Accession'])]
dbWriteTable(con_sqlite, "tblBioSample", df, overwrite = TRUE, row.names=F)

#---------------------------

annotation.tbl


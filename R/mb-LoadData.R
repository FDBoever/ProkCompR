pool <- dbPool(RSQLite::SQLite(), dbname = "MarinobacterGenome.db")
pool <- dbPool(RSQLite::SQLite(), dbname = "~/DATA/MarbGenomics/MarinobacterGenome.db")

#assign variables to pool connections
tblAssembly <- tbl(pool, "tblAssembly")#%>% unique()
tblBioSample <- tbl(pool, "tblBioSample")#%>% unique()
tblBioProject  <- tbl(pool, "tblBioProject")#%>% unique()

tblAssembly <- tblAssembly %>% collect()

#----
 write.table(tblBioSample %>% select(Organism, strain,lat_lon,geo_loc_name, host, isolation_source) %>% collect  , file='~/DATA/MarbGenomics/habitat.txt',sep='\t',quote=FALSE,row.names=FALSE)





 tblAssembly


tblBioSample %>% select(Organism, strain,lat_lon,geo_loc_name, host, isolation_source)
tblBioSample %>% select(Organism, strain,lat_lon,geo_loc_name, host, isolation_source) %>%
tblBioSample %>% select(host) %>% pull() %>% unique()

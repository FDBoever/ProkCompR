#


 loadPATRICtxt <- function(path){
  files<-list.files(path=path,
                    pattern='.txt',
                    full.names=TRUE)

  ff <- NULL
  for(f in files){
      ff <- base::rbind(ff, utils::read.delim(f,header=TRUE))
  }
  ff <- as_tibble(ff)
  return(ff)
 }


 patric_txt <- loadPATRICtxt(path = "~/DATA/bac_adapt_alg/PATRIC_txt")
 patric_txt %>% unique() %>% select(Genome.Name) %>% data.frame()

 patric_txt %>% group_by(Host.Name) %>% tally() %>% data.frame()

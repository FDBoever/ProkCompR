# ANGST maker
# Retrieve gene families with at least 3 proteins
angst.OG <- annotation.tbl %>%
  filter(genome %in% sel.genome) %>%
  filter(type=='CDS') %>%
  group_by(OG) %>% tally() %>%
  filter(n>3) %>%
  select(OG) %>%
  pull()

#How many do we have
length(angst.OG)

og.sequences <- tblAnnotation %>%
  filter(genome %in% sel.genome) %>%
  filter(OG %in% angst.OG) %>%
  select(genome,OG,locus_tag,aa_sequence,nuc_sequence) %>% collect()

#Angst requires removal of _ and .
og.sequences <- og.sequences %>%
  mutate(gnm = gsub('_','',gsub('\\.','',genome))) %>%
  mutate(loctg = gsub('_','',gsub('\\.','',locus_tag))) %>%
  mutate(hdr = paste(gnm,loctg,sep='_'))


#write each SCO down
for(i in 1:length(angst.OG)){
  message('writing ',angst.OG[i])
  #save faa
  sel.og.seq <- og.sequences %>% filter(OG == angst.OG[i]) %>% select(hdr,aa_sequence) %>% data.frame()
  f.name = paste0("~/DATA/MarbGenomics/angst_loci/aa/",angst.OG[i],'.faa')
  writeFasta(data=sel.og.seq, filename=f.name)

  #save fna
  sel.og.seq <- og.sequences %>% filter(OG == angst.OG[i]) %>% select(hdr,nuc_sequence) %>% data.frame()
  f.name = paste0("~/DATA/MarbGenomics/angst_loci/nuc/",angst.OG[i],'.fna')
  writeFasta(data=sel.og.seq, filename=f.name)
}


blastp

blastp <- read.delim('~/DATA/MarbGenomics/testblastp.txt',header=FALSE)

hits <- blastp %>% select(V2) %>% pull() %>% as.character()

annotation.tbl %>% filter(locus_tag %in% hits) %>% select(locus_tag, genome, product, OG)


hittedOGs <- annotation.tbl %>% filter(locus_tag %in% hits) %>%
  select(OG) %>% pull() %>%
  as.character()

annotation.tbl %>% filter(OG %in% hittedOGs) %>%
  select(locus_tag, genome, product, OG) %>%
  filter(genome=='Marinobacter_algicola_DG893')


#=================
title_suffix="NAD|NADH"#!!!!!!
title_suffix="transport|import|permease"
title_suffix="chemotaxis|flag"
title_suffix="ATP synthase|ATPase"
title_suffix="catalase|superoxide|dismutase"

annotation.tbl %>%
  #filter(grepl('transport|import|permease', product, ignore.case=TRUE)) %>%
  filter(grepl(title_suffix, product, ignore.case=TRUE)) %>%
  select(locus_tag, genome, product, OG) %>%
  filter(genome=='Marinobacter_algicola_DG893') %>%
  data.frame()



transportome <- annotation.tbl %>%
  #filter(grepl('transport|import|permease', product, ignore.case=TRUE)) %>%
  filter(grepl(title_suffix, product, ignore.case=TRUE)) %>%
  select(locus_tag, genome, product, OG) %>%
  select(OG) %>%
  filter(!is.na(OG)) %>%
  pull() %>%
  as.character() %>%
  unique()

#OG.pan[sel.genome,transportome]
dist.pan <- pan2dist(pan = OG.pan[sel.genome,transportome],
         dist.method="Manhattan")
COGpan.p.pcoa <- pan2pcoa(dist.pan = dist.pan,
                          metadata=metdad,
                          color="group",
                          title_suffix="COG")#,
                          #color_scale=clique.colors)

metadata=metdad
color="group"
MDS = cmdscale(dist.pan, eig = TRUE)
MDSdata = data.frame(MDS$points)
MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
MDSdata$Name = rownames(MDSdata)
#MDSdata$Group = colorGroups[MDSdata$Name]
colnames(MDSdata) = c("x","y","Name")

MDSdata <- MDSdata %>% left_join(metadata,by=c('Name'='genome'))

p <- ggplot(MDSdata,aes(x=x,y=y,fill=group,label=genome))
p <- p +
  theme_classic() +
  xlab(paste("PCoA1 (", MDS.var.perc[1], "%)",sep="")) +
  ylab(paste("PCoA2 (", MDS.var.perc[2], "%)",sep="")) +
  ggtitle(label=title_suffix)+
  geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
  geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
  geom_point(shape = 21, size = 2) +
  fdb_style() + theme(legend.position = "none")

p

og.p <- OG.pan[sel.genome,transportome]
#og.p[og.p==0] <- NaN

gplots::heatmap.2(as.matrix(og.p), trace ="none", margins = c(2,10), cexRow = 0.5)

#==============================









transportome <
  filter(genome=='Marinobacter_algicola_DG893')

#Multiblast vis
# ~/Genomics/ncbi-blast-2.7.1+/bin/blastp -query ./pimped2.fasta -db protDB_marinobacter -out ./pimped2.results.out -outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps evalue bitscore qcovs qcovhsp sseq" -evalue 1e-30
#~/Genomics/ncbi-blast-2.7.1+/bin/blastp -query ./vibrioferrin_algicola.fasta -db protDB_marinobacter -out ./vib_alg.results.out -outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps evalue bitscore qcovs qcovhsp sseq" -evalue 1e-30


##
#
#
#Possibly need to remove T6SS, as outgroup does have them ?
# T6SS
# beta lactamase... whats the point, it seems from everywhere


#"hp15 2692" OR "hp15 2693" OR "hp15 2694" OR "hp15 2695" OR "hp15 2696" OR "hp15 2697" OR "hp15 2698" OR "hp15 2699" OR "hp15 2700" OR "HP15_2708"


balst.res <- read.delim('~/DATA/MarbGenomics/pimped2.results.out',header=FALSE)
balst.res <- read.delim('~/DATA/MarbGenomics/vib_alg.results.out',header=FALSE)


colnames(balst.res) =c('qseqid','sseqid','slen', 'qstart', 'qend', 'length', 'mismatch', 'gapopen', 'gaps', 'evalue', 'bitscore',  'qcovs', 'qcovhsp', 'sseq')
balst.res <- balst.res %>% as_tibble()

loci.ann <- read.delim('~/DATA/MarbGenomics/grouped_loci.txt',header=TRUE)

balst.res.ann <- balst.res %>%
  tidyr::separate(qseqid, into=c('tr','qentry','qentryname'),sep='\\|') %>%
  select(-tr) %>% left_join(annotation.tbl, by=c('sseqid'='locus_tag')) %>%
  left_join(loci.ann, by=c('qentry'='Entry'))

balst.res.ann

module.sizes <- loci.ann %>% group_by(Group) %>% tally


plot.table <- balst.res.ann %>% group_by(genome, qentry, Group) %>%
  tally() %>%
  mutate(n=ifelse(n>=1,1,0)) %>% #Sets multiple hits to one
  ungroup() %>%
  group_by(Group, genome) %>%
  dplyr::summarise(sum=sum(n)) %>% #counts per module
  left_join(module.sizes,by='Group')%>%
  mutate(mcr = sum/n) %>% select(Group, genome, mcr) %>%
  pivot_wider(id_cols=Group, names_from=genome,values_from=mcr) %>% #pivot_wider and longer sequenially to add NAs
  pivot_longer(-Group,names_to='genome',values_to='mcr') %>%
  mutate(Group=factor(Group,levels =(c("RnfABCDGE","rnfABDGE","nqrABCDEF","nuoABCEFGHIJKLMN" , "ccmABCD","F_ATPase_I","F_ATPase_II", "arcABCD","H2ase","mnh_antiport","pcrABCD","Xanthorohopsin" ,"amt","Nrt","NrtABC","ethanolamine","NirBD","NarGHI","NorBD","N2OR","NirS","Nitroreductase","nifH","urea_amidolyase","urease_complex" ,"phosphonate","dsrEFHC","SorAB","SQR","tauACD","alkL","aupAB","alkB","almA","cP450", "2-nitropropane dioxygenase ", "anaerobic_hydrocarb","Benzene_degradation", "Benzoate_degradation", "Catechol_MetaCleavage","Catechol_OrthoCleavage","cyclohexanone monooxygenase","atuABCDEFGH","NPD","alginate biosynthesis","capD","pelABCDEFG","polysaccharide" ,"Arsenic_resist","Copper_resistance","CzcCBA","Mercuric_resistance","Terullium_resistance" ,"B1_transport","B12_1","B12_2","B12_transp","PQQABCDE","beta_lactamase","Desferox", "hemin_transp","iron_hydroxamate_tp" ,"petrobactin" ,"petro_uptake","pvsABCDE","pvuABCD","glcDEF", "carbohydrate_tp", "maltose_tp","PTS_fructose" ,"sugar_tp","glycogen","phaABC" ,"Fimb_comp","Flp","PilMNOPQ", "T6SS"))))



p.bast.out <- plot.table  %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=tip.order))%>%
  mutate(mcr=ifelse(mcr<0.5,NA,mcr)) %>%
  ggplot(aes(Group,genome))+
  geom_point(aes(fill=mcr),size=5,shape=21,color='grey')+
  theme_minimal() +
  scale_fill_continuous(low="darkgrey", high="black",
                        guide="colorbar",na.value="white")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+#coord_flip()+
  #facet_grid(~TYPE,scales='free_x',space = "free_x")
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','p.bast.out','.pdf'),plot=p.bast.out, width = 20, height = 20,unit='in')



module.clust <- hclust(dist(t(m.pan[,selected_modules])))
module.order <- module.clust$labels[module.clust$order]

namesConv <- m2n_filtered$DESCRIPTION
names(namesConv) = m2n_filtered$module
description.order <- unique(as.character(unname(namesConv[module.order])))






#====

dotplot <- plot.table  %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=tip.order))%>%
  mutate(mcr=ifelse(mcr<0.5,NA,mcr)) %>%
  ggplot(aes(Group,genome))+
  geom_point(aes(fill=mcr, size=mcr),shape=21,color='grey')+
  theme_minimal() +
  scale_fill_viridis_c(name = '') +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')+
  scale_y_discrete(position = "right")


ggtree_plot <- ggtree::ggtree(lsTrees.94.rooted[[1]])+geom_tiplab(align=TRUE, size=0)
ggtree_plot



cowplot::plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')
cowplot::plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')




#======================

mat <- plot.table %>%
  filter(genome %in% tip.order) %>%
  select(-Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = genome, values_from = mcr) %>%
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Group  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
mat[is.na(mat)]=0

v_clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
############ NOTICE THE t() above)

ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()


dotplot <- plot.table  %>%
  filter(genome %in% tip.order) %>%
  mutate(genome = factor(genome,level=tip.order))%>%
  mutate(Group = factor(Group, levels = v_clust$labels[v_clust$order])) %>%
  mutate(mcr=ifelse(mcr<0.5,NA,mcr)) %>%
  ggplot(aes(Group,genome))+
  geom_point(aes(fill=mcr, size=mcr),shape=21,color='grey')+
  theme_minimal() +
  scale_fill_viridis_c(name = '') +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')+
  scale_y_discrete(position = "right")

ggtree_plot_col <- ggtree_plot_col + ggtree::xlim2(dotplot)

labels <- ggplot(plot.table %>%
                   mutate(`Cell Type` = Group,
                          Group = factor(Group, levels = v_clust$labels[v_clust$order])),
                 aes(x = Group, y = 1, fill = `Cell Type`)) +
  geom_tile() +
  scale_fill_brewer(palette = 'Set1') +
  theme_void() + theme(legend.position = "none")

legend <- cowplot::plot_grid(cowplot::get_legend(labels + theme(legend.position="bottom")))

p.bast.out<- patchwork::plot_spacer() + patchwork::plot_spacer() + ggtree_plot_col +
patchwork::plot_spacer() + patchwork::plot_spacer() + labels +
patchwork::plot_spacer() + patchwork::plot_spacer() + patchwork::plot_spacer() +
ggtree_plot + patchwork::plot_spacer() + dotplot +
patchwork::plot_spacer() + patchwork::plot_spacer() + legend +
patchwork::plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 0.1, -0.1, 4, 1)) + theme(legend.position = 'none')

ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/','p.bast.out_testing','.pdf'),plot=p.bast.out, width = 30, height = 20,unit='in')






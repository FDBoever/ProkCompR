
library(ggtree)
library(ape)
library(data.table)
library(RColorBrewer)
library(phytools)
library(ggplot2)
###########################
tree = read.tree(file='~/DATA/phyloCCAP/trees/18S_SINA_FastTree.tree')

#rownames(alkB_cytoscape) = alkB_cytoscape$name
#mtdat = read.delim('~/DATA/phyloCCAP/EukProk/EukProk_rRNA_aligned.csv',sep=',')
#mtdat = read.delim('~/DATA/phyloCCAP/PR2_18S_annotated_Decipher.txt',sep='\t')

#sel.col = 'sc0.8'

mtdat <-genome.tbl %>%
  dplyr::left_join(sorted.cliques, by='genome') %>%
  dplyr::left_join(out.colors, by='genome') %>%
  #dplyr::filter(!is.na(!!as.name(sel.col))) %>%
  data.frame()
mtdat$grp <- mtdat[,sel.col]
mtdat <- mtdat %>%
  dplyr::mutate(group=paste0('cl',grp)) %>%
  dplyr::select(-grp)
mtdat[mtdat$group=='clNA','group'] <- NA


rownames(mtdat)<- mtdat$genome

#tree = drop.tip(mdTree,setdiff(mdTree$tip.label,as.character(mtdat[!is.na(mtdat $order),'name'])))
tree = sco.106.concat.tree
mtdat = mtdat[tree $tip.label,]

tree = ape::ladderize(ape::root(tree,outgroup),right=TRUE)

rankSelected_1 = 'group'
rankSelected_2 = 'group'
selectedBarVariable = 'GC'


my_info <- data.table(tip_lbs = mtdat$genome,
                      groupA = mtdat[,c(rankSelected_1)],
                      groupB = mtdat[,c(rankSelected_2)],
                      val = mtdat[,c(selectedBarVariable)])
grA <- split(my_info$tip_lbs, my_info$groupA)
tree_grA <- ggtree::groupOTU(tree, grA)
str(tree_grA)


levels(attributes(tree_grA)$group)[1] <- levels(attributes(tree_grA)$group)[2]

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#fcolors = c("#999999","#629363","#3A85A8","#C4625D","#E41A1C","#FFC81D","#BF862B","#EB7AA9")

fcolors = getPalette(length(unique(c(as.character(my_info$groupA), as.character(my_info$groupB))))+1)
fcolors2 = getPalette(length(unique(my_info$groupB)))

#fcolors = getPalette(20)
#A_fcolors = getPalette(unique(my_info$groupA))


tree_plot <-
  ggtree(tr = tree_grA,
         branch.length='none',
         # color by group attribute, check str(tree_grA)
         mapping = aes(color = group),
         layout  = 'circular',
         # set line thikness
         size = 0.4) + scale_color_manual(name = 'Group A',
                                          values = fcolors) #+ theme(legend.position='none')
tree_plot

tree_plot <- tree_plot %<+% my_info

tree_dt <- data.table(tree_plot$data)
head(tree_dt)

# select only the tip labels and order by coord y
tree_dt <- tree_dt[isTip == TRUE][order(y)]

# Make table with y cords for each groupB in consecutive order;
# this helps for drawing & labeling segments.
# Note the usage of "rleid" function, which is a grouping ID generator,
# needed because consecutive rows of an identical reoccurring group must form a unique group.
coord_groups <- tree_dt[, .(y1 = y[1],
                            y2 = y[.N],
                            angle = mean(angle),
                            n = .N), # optional - helps with counting
                        by = .(groupB,
                               id_gr = rleid(groupB,
                                             prefix = "grp"))]
coord_groups

# Compute the middle y - will be used for placing the group label;
# similarly the mean angle was computed above already.
coord_groups[, y_mid := rowMeans(.SD), .SDcols = c("y1", "y2")]

# For one-record groups where y1=y2, adjust their y coordinates so that a segment gets drawn;
# If not, no segment get drawn for such cases.
# To force a segment add and subtract a small amount (try & error until seems ok).
# Prefer a smaller values since with bigger ones you risk to exaggerate the segments.
coord_groups[, y1_adj := ifelse(y1 == y2, y1 - 0.1, y1)]
coord_groups[, y2_adj := ifelse(y1 == y2, y2 + 0.1, y2)]

# Labels need angle adjustment for cases between 90 and 270 dg
coord_groups[, angle_adj := ifelse(angle %between% c(90, 180),
                                   yes = angle + 180,
                                   no = ifelse(angle > 180 & angle <= 270,
                                               yes = angle - 180,
                                               no = angle))]

# Labels with angles between 90 and 270 dg
# need change of horizontal adjustment argument from 0 to 1.
coord_groups[, hjust_adj := ifelse(angle %between% c(90, 270), yes = 1L, no = 0L)]

# if needed, coloring could be binary
# coord_groups[, col := ifelse(.I%%2, 0.5, 1)]
coord_groups
colourCount = length(unique(coord_groups$groupB))


# Define variable to control x coordinate of segments & labels
my_x <- max(tree_dt$x) + 1.5 #+ 0.05

tree_labeled <-
  tree_plot +
  geom_segment(data = coord_groups,aes(x = my_x, y = y1_adj, xend = my_x, yend = y2_adj, color = groupB),lineend = "butt",size = 3) +
  geom_text(data = coord_groups,aes(x = my_x,y = y_mid,angle = angle_adj,hjust = hjust_adj,label = groupB),vjust = 0.5, size  = 2.5,nudge_x = 0.05, color = "black") +
  theme(
    text = element_text(size = 10, family = "sans"),
    legend.justification = c(1,0),
    #   legend.position = 'c(1.1, 0.05)',
    #    legend.position = 'c(1.1, 0.05)',

    #legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.key.height = unit(4, "mm"),
    plot.margin = unit(c(t = 0.2, r = 1.3, b = -0.35, l = -0.2), "cm")
  )
tree_labeled


#tree_labeled$data <- tree_labeled$data %<+% mtdat
mtdat['Marinobacter_salsuginis_SD_14B','Genome_size'] <- NA
mtdat[outgroup,'Genome_size'] <- rep(NA, length(mtdat[outgroup,'Genome_size']))
mtdat[outgroup,'GC'] <- rep(NA, length(mtdat[outgroup,'GC']))
mtdat[outgroup,'Coding_density'] <- rep(NA, length(mtdat[outgroup,'Coding_density']))
mtdat[outgroup,'b'] <- rep(NA, length(mtdat[outgroup,'b']))


tree_labeled <- tree_labeled %<+% mtdat
tree_labeled$data<-tree_labeled$data %>% mutate(group = group.y)




#================================================================================#


tipsize = 2
cols <- c("100" = "black", "95-99" = "grey40", "85-95" = "grey", "75-85" = "lightgrey", "<75" = "white")

circ <- ggtree::ggtree(tr.tree, layout="circular") +
  geom_tiplab(align=TRUE,size=0,linetype='dashed',linesize=0.3)

circ$data$bootstrap <- as.numeric(circ$data$label)
circ <- circ + ggtree::geom_point2(aes(subset=!is.na(bootstrap) & as.numeric(bootstrap) == 100,label = "100", fill="100") ,pch=21) +
  ggtree::geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) >95 & bootstrap < 100, label="95-99", fill="95-99"), pch=21) +
  ggtree::geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) > 85 & bootstrap < 95, label="85-95", fill="85-95"), pch=21) +
  ggtree::geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) > 75 & bootstrap < 85, label="75-85", fill="75-85") ,pch=21) +
  ggtree::geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) < 75 , label="<75", fill="<75") ,pch=21) +
  scale_fill_manual(name = "bootstrap", breaks=c("100", "95-99", "85-95", "75-85", "<75"), values = cols, labels = c("100", "95-99", "85-95", "75-85", "<75"))


circ$data <- circ$data %>% left_join(mtdat, by=c('label'='genome'))

tree_labeled <- circ


library(ggstar)
library(ggnewscale)
#=======================

c.offset=0.06
c.size=0
c.width=0.05

p3 <- tree_labeled +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=ggstar::geom_star,
    mapping=aes(y=label, fill=TypeStrain, starshape=TypeStrain),
    size=2.5,
    starstroke=0,
    pwidth=0.1,
    inherit.aes = FALSE,
    offset=0.03,
    alpha=0.5
    )+
  scale_fill_manual(
    name="Type strain",
    values=c("white", "black", "#FF0100")) +
  scale_starshape_manual(
    values=c(1,28),
    guide="none"
  ) +
  new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
             mapping=aes(fill=Genome_size),
             color = "white", offset = 0.05,size = c.size,width=c.width)+
  scale_fill_distiller( palette = 'Oranges',direction = 1)+

  new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
             mapping=aes(fill= Coding_density),
             color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_distiller( palette = 'Blues',direction = 1)+
  new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
             mapping=aes(fill= GC),
             color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_distiller( palette = 'Greens',direction = 1)+

  new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
             mapping=aes(fill= b),
             color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_distiller( palette = 'Purples',direction = 1)


#================================================================================#
clique.colors <- cliquelevel2colors(sorted.ANI.cliques,out.colors,'sc0.8')

c.offset=0.04
c.size=0 #turn off line width
c.width=0.03

p4 <- p3 +  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=group),
                          color = "white", offset=0.08,size=c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  #scale_fill_manual(values=c('red','grey','grey20','grey30','grey40','grey50','grey60','grey70','grey80','black'))+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=group),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=clique.colors)+
  ggnewscale::new_scale_fill() +

  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=group,alpha=-as.numeric(as.character(cc0.85))),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=clique.colors)+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=group,alpha=-as.numeric(as.character(cc0.9))),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=clique.colors)+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=group,alpha=-as.numeric(as.character(cc0.95))),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+scale_alpha(na.value=0)+
  scale_fill_manual(values=clique.colors)+

  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(geom=geom_tile,
                          mapping=aes(fill=group,alpha=-as.numeric(as.character(cc0.98))),
                          color = "white", offset = c.offset,size = c.size,width=c.width)+
  scale_fill_discrete(na.value = 'white')+
  scale_fill_manual(values=clique.colors)+
  scale_alpha(na.value=0,range=c(0.5,1))


p7 <- p4 + geom_treescale(fontsize=2, linesize=0.3) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"))

#generate different layouts
p.out.circular <- ggtree::open_tree(p7, 15)
p.out.opened <- ggtree::open_tree(p7, 180)
p.out.rotated <- ggtree::rotate_tree(p.out.opened, -90)
p.out.rectangular <- p7 + layout_rectangular()

# Save down what you want
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ML_SCO_annotated.circular.pdf',plot=p.out.circular, width = 15, height = 15,unit='in')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ML_SCO_annotated.opened.pdf',plot=p.out.opened, width = 15, height = 15,unit='in')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ML_SCO_annotated.rotated.pdf',plot=p.out.rotated, width = 15, height = 15,unit='in')
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ML_SCO_annotated.rectangular.pdf',plot=p.out.rectangular, width = 7, height = 15,unit='in')

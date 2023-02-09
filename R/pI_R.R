



#play arround howmany of the plots you want
multiplot(plotlist = myplots[1:12],cols=6)
multiplot(plotlist = myplots,cols=9)
multiplot(plotlist = myplots2,cols=9)
multiplot(plotlist = myplots,cols=6)
multiplot(plotlist = myplots,cols=9)


##################################################
# PI BIASS per genome

q = ggplot(bias_data, aes(x=reorder(genome,bias), y=bias)) +
    geom_point() +      # Thinner lines
    xlab("Genome") +
    ylab("pI bias") +
    ggtitle("pI bias") +
    theme_bw()
q + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) + coord_cartesian(ylim = c(45, 70)) + coord_flip()


bias_data$genome <- factor(bias_data $genome, levels = levels(pI$genome))


q = ggplot(bias_data, aes(x=genome, y=bias)) +
    geom_point() +      # Thinner lines
    xlab("Genome") +
    ylab("pI bias") +
    ggtitle("pI bias") +
    theme_bw()
q + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) + coord_cartesian(ylim = c(45, 70)) + coord_flip()








#show violin graph
 p <- ggplot(pI, aes(factor(genome), pI))
p + geom_violin()+theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) + coord_cartesian(ylim = c(45, 70)) + coord_flip()


##################################################

# Select genomes of interest from the set
pIsub1 =subset(pI, genome == "Marinobacter_lipolyticus_SM19.faa")
pIsub2 =subset(pI, genome == "Marinobacter_fuscus_NH169_3.faa")
rbind(pIsub1,pIsub2)


g2 = ggplot(rbind(pIsub1,pIsub2), aes(x= pI, color= genome, fill= genome)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.2,binwidth=0.1)+
geom_density(alpha=0.2)+scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
labs(title="",x="pI", y = "Density")+
theme_classic()






bias_annotation = data.frame(bias_data[,c('b')])
rownames(bias_annotation)= bias_data $genome
colnames(bias_annotation)='b'
bias_annotation$b = as.numeric(as.character(bias_annotation$b))
p <- ggtree(tree, layout = "circular", open.angle = 50,branch.length="none", size=0.5)

p2 = open_tree(p, 180)

gheatmap(p2, bias_annotation,offset=-1, width=0.2)+
  scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdBu')),name='pI bias')+theme(
		aspect.ratio = 1,
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
        #legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.key.width = unit(0.25, 'cm'),
        legend.key.height = unit(0.55, 'cm'),
        legend.text = element_text(size = 10),
        #legend.text.align = 1,
        legend.title = element_text(size = 11))






p <- ggtree(tree, layout = "circular", open.angle = 50, size=0.5)

p2 = open_tree(p, 180)

gheatmap(p2, bias_annotation,offset=-1, width=0.2)+
  scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdBu')),name='pI bias')+theme(
		aspect.ratio = 1,
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
        #legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.key.width = unit(0.25, 'cm'),
        legend.key.height = unit(0.55, 'cm'),
        legend.text = element_text(size = 10),
        #legend.text.align = 1,
        legend.title = element_text(size = 11))







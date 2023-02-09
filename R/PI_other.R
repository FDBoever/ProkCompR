library(ggplot2)
library(plyr)
library(dplyr)

################################
# Set colors for Marinobacter story
################################

colors = c(
rgb(57,129,29, maxColorValue=255),#adhaerens
rgb(134,218,129, maxColorValue=255),#alg
rgb(166,216,212, maxColorValue=255),#ant
rgb(242,169,104, maxColorValue=255),#ex

rgb(242,236,112, maxColorValue=255), #hydrocarb
rgb(227,143,221, maxColorValue=255),#lipo

rgb(137,137,137, maxColorValue=255),#other
rgb(118,174,207, maxColorValue=255),#psych
rgb(131,115,27, maxColorValue=255),#sedim
rgb(179,77,34, maxColorValue=255)#vinifirmus
)


#Read file
#pI = read.table("~/Genomics/ProkComp/pI.txt",header=TRUE)
#pI = read.table("~/Genomics/ProkComp/pIbacterioplankton.txt",header=TRUE)

ProtAna=  read.delim('~/DATA/pI_final.csv')
pI = ProtAna
pI = pI[!pI$genome %in% c('Hahella_ganghwensis.faa', 'Hahella_sp_CCB-MM4.faa','Oleiphilus_messinensis.faa'),]


ProtAnaCORE= read.delim('~/DATA/pI_final__core.csv')
head(ProtAnaCORE)
pI = ProtAnaCORE

ANVIO_cat = read.table('~/DATA/MarinobacterGenomics/2018_ProkComp/ANVIO_CAT.txt',header=TRUE,sep="\t")
ANVIO_cat = ANVIO_cat[ANVIO_cat$name %in% unique(pI$genome),]
colnames(ANVIO_cat) = c('genome','group','group2','lifestyle','SS','sourceClass','lat','lon')

pI = right_join(pI, ANVIO_cat, by = 'genome')
pI$acidic <- ifelse(pI$pI > 7.5 ,"basic", "acidic")


ggplot(pI, aes(pI, instability))+geom_point()
ggplot(pI, aes(instability, Aromaticity))+geom_point()
ggplot(pI, aes(instability, Aromaticity))+geom_point(aes(color=group2))+scale_color_manual(values=colors)


ggplot(pI, aes(instability, group2))+geom_point(aes(color=group2))+scale_color_manual(values=colors)
ggplot(pI, aes(Aromaticity, group2))+geom_point(aes(color=group2))+scale_color_manual(values=colors)
ggplot(pI[pI$pI>7.5,], aes(pI, group2))+geom_point(aes(color=group2))+scale_color_manual(values=colors)


average_pI = aggregate(pI, list(pI$group2,pI$PC), mean)

#average_pI <- aggregate(pI, by = list(pI $id1, pI $id2), mean)
head(average_pI)
ggplot(average_pI, aes(instability, Group.1))+geom_point(aes(color= Group.1))+scale_color_manual(values=colors)
ggplot(average_pI, aes(pI, Group.1))+geom_point(aes(color= Group.1))+scale_color_manual(values=colors)
ggplot(average_pI, aes(x=pI, color= Group.1))+geom_histogram(fill='white',alpha=0.3,position='identity')+scale_color_manual(values=colors)

average_pI $acidic <- ifelse(average_pI $pI > 7.5 ,"basic", "acidic")
ggplot(average_pI, aes(x=pI, y= Group.1,color=acidic))+geom_point()


ggplot(average_pI, aes(Group.1))+geom_bar(aes(fill=acidic), stat="identity")


ggplot(average_pI, aes(Group.1,fill = acidic))+geom_bar(stat="count", position = "dodge")+coord_flip()+scale_fill_manual(values=c('black','grey'))+theme_classic()

ggplot(average_pI[average_pI$Group.1 != 'other',], aes(Group.1,fill = acidic))+geom_bar(stat="count")+scale_fill_manual(values=c('black','grey'))+ scale_y_continuous(expand = c(0, 0))+xlab('')+ylab('nr of PCs')+theme_classic()+coord_flip()+ labs(title = "CORE genes", fill = "")



pI$genome <- factor(pI$genome, levels = tree3$tip.label[ordered_tips][4:61])

pI = pI[!is.na(pI$genome),]
ggplot(pI, aes(genome,fill = acidic))+geom_bar(stat="count")+scale_fill_manual(values=c('black','grey'))+ scale_y_continuous(expand = c(0, 0))+xlab('')+ylab('nr of PCs')+theme_classic()+coord_flip()+ labs(title = "CORE genes", fill = "")





##################################################
######## VALUABLE FUNCTIONS #########
# multiplot() will allow you to visualize multiple plots in a grid arrangment

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##################################################



#add additional column to classify protein as basic vs acidic
#play with the threshold number
pI$acidic <- ifelse(pI$pI > 7.5 ,"basic", "acidic")


#for up to 3 proteomes
#g2 = ggplot(pI, aes(x= pI, color= genome, fill= genome)) +
#geom_histogram(aes(y=..density..), position="identity", alpha=0.2)+
#geom_density(alpha=0.2)+
#scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#labs(title="",x="pI", y = "Density")+
#theme_classic()

#global view, shows all the proteomes on one graph
g1 = ggplot(pI, aes(x= pI, color= "black", fill= "black")) +
geom_histogram(aes(y=..density..), position="identity", alpha=1,binwidth = 0.1)+ scale_y_continuous(expand = c(0, 0))+scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black"))+ 
labs(title= strBias ,x="class of pI", y = "Density")+
theme_classic()+theme(legend.position="none")


g2 = ggplot(pI, aes(x=pI,y= log(length),color = acidic)) + geom_point()+
scale_color_manual(values=c("grey", "black")) + scale_fill_manual(values=c("grey", "black"))+ labs(title=" ",x="pI", y = "log Length")+ theme_classic()+theme(legend.position="none")


multiplot(g1,g2,cols=2)



##################################################
# plot visialisation per proteome


myplots <- list()
myplots2 <- list()
geno= c()
bias = c()

i = 1
for (organism in unique(pI$genome)){
	pIsub =subset(pI, genome ==organism)
	Nbasic = length(which(pIsub$acidic == 'basic'))
	Nacidic = length(which(pIsub$acidic == 'acidic'))
	pIBias = (Nbasic-Nacidic)/(Nbasic+Nacidic)*100
	strBias = paste(organism,":",round(pIBias,2),"%",sep=" ")
	
	bias = c(bias, pIBias)
	geno = c(geno, organism)
	p1 <- ggplot(pIsub, aes(x= pI, color= "black", fill= "black")) + geom_histogram(aes(y=..density..), position="identity", alpha=1,binwidth = 0.1) + scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black"))+ labs(title= strBias, x="class of pI", y = "Density")+ theme_classic()+theme(legend.position="none",plot.title = element_text(size=6))
	g2 = ggplot(pIsub, aes(x=pI,y= log(length),color = acidic)) + geom_point()+scale_color_manual(values=c("grey", "black")) + scale_fill_manual(values=c("grey", "black"))+ labs(title= strBias,x="pI", y = "log Length")+ theme_classic()+theme(legend.position="none",plot.title = element_text(size=6))
	myplots[[i]] <- p1 
	myplots2[[i]] <- g2 
	i = i + 1
}
#stores all the calculated b-values in a new data.frame
bias_data = data.frame(cbind(genome=geno,b=bias))


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


bias_data $genome <- factor(bias_data $genome, levels = levels(pI$genome))


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
pIsub2 =subset(pI, genome == "Marinobacter_lipoliticus_BF04_CF-4.faa")
rbind(pIsub1,pIsub2)


g2 = ggplot(rbind(pIsub1,pIsub2), aes(x= pI, color= genome, fill= genome)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.2,binwidth=0.1)+
geom_density(alpha=0.2)+scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
labs(title="",x="pI", y = "Density")+
theme_classic()

g3 = ggplot(rbind(pIsub1,pIsub2), aes(x= pI, color= genome, fill= genome)) +
geom_density(alpha=0.2)+scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
labs(title="",x="pI", y = "Density")+
theme_classic()

multiplot(g2,g3,cols=2)


AAcomp = read.table("/Users/sa01fd/Genomics/ProkComp/pI_AA.csv",header=TRUE)





data1=data.frame(cbind(AAcomp[1:2], AAcomp[6:25]))
head(data1)

longdf=melt(data1,id.vars=c("PC","genome"))
head(longdf)

df=ddply(longdf,.(genome,variable),summarize,mean(value))
colnames(df)[length(df)]="value"
df

df = right_join(df, ANVIO_cat, by = 'genome')

head(df)


ggplot(data=df,aes(x=variable,y=value,group= group2,colour= group2))+
        geom_point()+scale_color_manual(values=colors)


average_AA = aggregate(df, list(df $group2, df $variable), mean)
#average_AA = aggregate(df, list(df $genome, df $variable), mean)

ggplot(data= average_AA,aes(x=  Group.1,y=value,group= value,colour= Group.1))+
        geom_point()+scale_color_manual(values=colors)+facet_wrap(~ Group.2,sacles='free_y')

ggplot(data= average_AA,aes(x= Group.2,y=value,group= value,colour= Group.1))+
        geom_point()+scale_color_manual(values=colors)+coord_polar()+geom_polygon(alpha=0.001)



ggplot(data=df,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.001)+coord_polar()

ggplot(data=df,aes(x=variable,y=value,group=genome))+
        geom_point()+geom_polygon(alpha=0.001)+coord_polar()



ggplot(data=df,aes(x=variable,y=value,group=genome,colour=genome))+
        geom_point()+geom_polygon(alpha=0.1)+coord_polar()


rescale_df=function(data,groupvar=NULL){
        if(is.null(groupvar)) df=data
        else df=data[,-which(names(data) %in% groupvar)]
        
        select=sapply(df,is.numeric)
        df[select]=lapply(df[select], scales::rescale)
        if(!is.null(groupvar)) {
                df=cbind(df,data[[groupvar]])
                colnames(df)[length(df)]=groupvar
        }        
        df
}

rescaled=rescale_df(data1)
head(rescaled)


longdf2=melt(rescaled,id.vars="genome")
head(longdf2)


df2=ddply(longdf2,.(genome,variable),summarize,mean(value))
colnames(df2)[length(df2)]="value"
df2

unique(df2$genome)
ggplot(data=df2[df2$genome %in%c('Marinobacter_hydrocarbonoclasticus_ASM1536.faa','Marinobacter_psychrophilus_20041.faa','Marinobacter_vinifirmus_FB1.faa'),],aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.4)+coord_polar()


coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
        theta <- match.arg(theta, c("x", "y"))
        r <- if (theta == "x") 
                "y"
        else "x"
        ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
                direction = sign(direction),
                is_linear = function(coord) TRUE)
}



ggplot(data=df2,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.4)+coord_radar()


ggplot(data=df2,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.4)+coord_radar()+ylim(0,1)+
        theme(legend.position="bottom")+xlab("")+ylab("")

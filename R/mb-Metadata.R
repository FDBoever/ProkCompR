#------------------------------------------------------------
#	Genome Metadata
#	Marinobacter comparative genomics
#	Frederik De Boever, 26-nov-2020
#------------------------------------------------------------
#

#==================================================#
# Phylogenetic signal

#---------------
phylosignal.df <- genome.tbl %>% mutate(hemisphere=ifelse(lat>=0,'North','South')) %>%
  filter(grepl('Marinobacter|Tamil',genome))%>%
  filter(!is.na(SS)) %>%
  filter(SS!='') %>%
  select(genome,SS)


#collect tree
tr.ph <- lsTrees.106.rooted[['sco.106.concat']]

#drop if tip if not in phylosignal
tr.ph.s <- ape::drop.tip(tr.ph,
              tip = tr.ph$tip.label[!grepl(paste0(phylosignal.df$genome,collapse = '|'), tr.ph$tip.label)])


#we rescale the tree with lambda=0 to obtain a starlike tree
tree <- tr.ph.s
tree0 <- geiger::rescale(tree, model = "lambda", 0)


# load trait data
traits <- as.character(phylosignal.df$SS)
names(traits) <- phylosignal.df$genome


geiger::name.check(tree, traits)



#Phylogenetic signal with Pagel's lambda
#

# Use fitDiscrete (geiger) using an equal rates (ER) model for isolation source
# Geiger0fitter comparative model of discret data

# lambda 1 and lambda 0 trees combined with trait data
trait_geiger <- geiger::treedata(tree,traits)
trait_geiger_0 <- geiger::treedata(tree0,traits)


# Comparison of different models for trait evolution
trait_lambda_ER <- fitDiscrete(trait_geiger$phy, traits, type="discrete", model = "ER", transform='lambda', niter = 1000)
trait_lambda0_ER <- fitDiscrete(trait_geiger_0$phy, traits, type="discrete", model = "ER", transform='lambda', niter = 1000)


#likelihood ratio test  lambda vs lambda 0
# p > 0.05 -- Not rejects ER-model
trait_d_lambda_vs_lambda0 <- abs(2*(trait_lambda0_ER$opt$lnL-trait_lambda_ER$opt$lnL)) # 1.543
trait_p_lambda_vs_lambda0 <- pchisq(trait_d_lambda_vs_lambda0, 3-1, lower.tail=FALSE)










# BIOGEOGRAPHY
dat2 <- genome.tbl %>% mutate(hemisphere=ifelse(lat>=0,'North','South')) %>%
  filter(grepl('Marinobacter|Tamil',genome))%>%
  filter(!is.na(lat)) %>%
  filter(quality_class !='low')
#==================================================
# ABSOLUTE LATITUDE

#---
#polynomial function to capture hump shaped data
mod <- lm(GC~I(abs(lat))+I(abs(lat)^2), data=dat2)
m=mod

# function to write equation
# now only liniar model, need to modify this

lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));

  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  }
  as.character(as.expression(eq));
}

lm_eqn(mod)

p.gc.lat <- dat2 %>%
  ggplot(aes(x=abs(lat),y=GC))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2),color='grey20',size=0.5)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)))+
  scale_shape_manual(values = c(21,22),guide=FALSE)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdBu'), guide = F) +
  xlab("Absolute latitude") +
  ylab('%GC content')+fdb_style()

#mod <- lm(GC~I(abs(lat))+I(abs(lat)^2), data=dat2)
#p.gc.lat <- p.gc.lat + ggplot2::geom_text(x = 10, y = 61, label = paste(lm_eqn(mod)), parse = TRUE)

ggsave('~/DATA/MarbGenomics/Graphs/GC_vs_LAT.pdf',plot=p.gc.lat, width = 4, height = 3,unit='in')

#since plotting a polynomial the LM looks like this
summary(lm(GC~I(abs(lat))+I(abs(lat)^2), data=dat2))

#-----------------
p.cd.lat <- dat2 %>%
  ggplot(aes(x=abs(lat),y=Coding_density))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2),color='grey20',size=0.5)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)))+
  scale_shape_manual(values = c(21,22), guide = F)+
  scale_fill_gradientn(colours =  RColorBrewer::brewer.pal(11, 'RdBu'), guide = F) +
  xlab("Absolute latitude") +
  ylab('Coding density')+fdb_style()

ggsave('~/DATA/MarbGenomics/Graphs/Coding_density_vs_LAT.pdf',plot=p.cd.lat, width = 4, height = 3,unit='in')

p.cd.lat2 <- dat2 %>%
  ggplot(aes(x=abs(lat),y=Coding_density))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2),color='grey20',size=0.5)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)), guide = F)+
  scale_shape_manual(values = c(21,22), guide = F)+
  scale_fill_gradientn(colours =  RColorBrewer::brewer.pal(11, 'RdBu')) +
  xlab("Absolute latitude") +
  ylab('Coding density')+fdb_style()

ggsave('~/DATA/MarbGenomics/Graphs/Coding_density_with_legends.pdf',plot=p.cd.lat2, width = 4, height = 3,unit='in')


#since plotting a polynomial the LM looks like this
summary(lm(Coding_density~I(abs(lat))+I(abs(lat)^2), data=dat2))

#-----------#
p.b.lat <- dat2 %>% mutate(pIbias=as.numeric(b)) %>%
  ggplot(aes(x=abs(lat),y=pIbias))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2),color='grey20',size=0.5, guide = F)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)), guide = F)+
  scale_shape_manual(values = c(21,22), guide = F)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdBu'), guide = F) +
  xlab("Absolute latitude") +
  ylab('pI bias')+fdb_style()

ggsave('~/DATA/MarbGenomics/Graphs/pIbias_vs_LAT.pdf',plot=p.b.lat, width = 4, height = 3,unit='in')

#since plotting a polynomial the LM looks like this
summary(lm(as.numeric(b)~I(abs(lat))+I(abs(lat)^2), data=dat2))

p.c <- gridExtra::grid.arrange(p.gc.lat, p.cd.lat, p.b.lat,ncol=3)
ggsave('~/DATA/MarbGenomics/Graphs/abs_lat_vs_characteristics.pdf',plot=p.c, width = 7, height = 3,unit='in')


#=====================================#
# Correlate them among one another

p.cor.lat.1 <- dat2 %>%
  ggplot(aes(x=Coding_density,y=GC))+
  geom_smooth(method = "lm",color='grey20',size=0.5)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)))+
  scale_shape_manual(values = c(21,22), guide = F)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdBu'), guide = F) +
  xlab("Coding density") +
  ylab('%GC content')+ggpubr::stat_cor()+
  fdb_style()

p.cor.lat.2 <- dat2 %>% mutate(pIbias=as.numeric(b)) %>%
  ggplot(aes(x=Coding_density,y=pIbias))+
  geom_smooth(method = "lm",color='grey20',size=0.5)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)))+
  scale_shape_manual(values = c(21,22), guide = F)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdBu'), guide = F) +
  xlab("Coding density") +
  ylab('pI bias')+ggpubr::stat_cor()+
  fdb_style()

p.cor.lat.3 <- dat2 %>% mutate(pIbias=as.numeric(b)) %>%
  ggplot(aes(x=pIbias,y=GC))+
  geom_smooth(method = "lm",color='grey20',size=0.5)+
  geom_point(size=2,alpha=0.8,aes(shape=hemisphere,fill=abs(lat)))+
  scale_shape_manual(values = c(21,22), guide = F)+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdBu'), guide = F) +
  xlab("pI bias") +
  ylab('%GC content')+ggpubr::stat_cor()+
  fdb_style()

p.c <- gridExtra::grid.arrange(p.cor.lat.1,p.cor.lat.2, p.cor.lat.3,ncol=3)
ggsave('~/DATA/MarbGenomics/Graphs/interCOR_GC_pI_CodingDensity.pdf',plot=p.c, width = 7, height = 3,unit='in')


#===========#



# metadata has been collected from NCBI BioSample, IMG database, BacDive, LPSN, and relevant scientific literature
mergedEnv = read.delim('~/DATA/MarbGenomics/mergedEnv.txt')

#consistent formattting of genome names (and set as rownames)
mergedEnv$genome_name = gsub('-','_',gsub('\\.','_', mergedEnv$genome_name_dotted))
rownames(mergedEnv) = mergedEnv$genome_name

#HQ_genomes object retrieved from the CheckM quality filtering
HQ_genomes


mergedEnv[HQ_genomes,]
mergedEnv[HQ_genomes,]

envVariables = c("TypeStrain","lifestyle","SS","sourceClass","lat","lon","geo_loc_name","host","isolation_source","origin.Continent")

mergedEnv[HQ_genomes, envVariables]


hostAssociatedMB = mergedEnv[mergedEnv$host!='','genome_name']
length(hostAssociatedMB)


metadata = cbind(metadata, mergedEnv[rownames(metadata), envVariables])
metadata$hostAssociated = metadata$host != ''




#---------------------------------
#Explore
metadata %>% dplyr::group_by(sourceClass) %>%
  dplyr::count() %>% data.frame()




#-----------------------------------
LatitudonalRange <- genome.tbl %>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::group_by(phylogroup) %>%
  dplyr::summarise(minLat = min(lat), minLatName = genome[which.min(lat)],
            maxLat = max(lat), maxLatName = genome[which.max(lat)]) %>%
  dplyr::mutate(rangeLat = maxLat-minLat)

LatitudonalRange %>%
  ggplot(aes(phylogroup,rangeLat)) +
    geom_col()


LatitudonalRange %>%
  dplyr::filter(!is.na(phylogroup)) %>%
  ggplot() +
  geom_point(aes(minLat,phylogroup,color=phylogroup))+
  geom_point(aes(maxLat,phylogroup,color=phylogroup))+
  geom_segment(aes(x = minLat,
                   y = phylogroup,
                   xend = maxLat,
                   yend = phylogroup,
                   color=phylogroup),size = 1)+
  scale_color_manual(values=phylogroup.colors) + fdb_style()



plat <- genome.tbl %>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::filter(!is.na(phylogroup)) %>%
  ggplot(aes(phylogroup,lat)) +
    geom_violin(aes(fill=phylogroup),color=NA,alpha=0.5,width=2)+
    geom_jitter(aes(color=phylogroup))+
    scale_color_manual(values=phylogroup.colors) +
    scale_fill_manual(values=phylogroup.colors) +
  geom_hline(yintercept=0,linetype='dashed',color='grey70') +
  geom_hline(yintercept=60,linetype='dashed',color='grey70') +
  geom_hline(yintercept=-60,linetype='dashed',color='grey70') +
  scale_y_continuous(breaks=c(-60,-30,0,30,60))+
  fdb_style(aspect.ratio=0.75)+theme(legend.position = 'none')
plat
ggsave('~/DATA/MarbGenomics/Graphs/latidutonal_distribution_phylogroups.pdf',plot=plat, width = 4, height = 4,unit='in')


genome.tbl %>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::filter(!is.na(phylogroup)) %>%
  mutate(abslat = abs(lat)) %>%
  ggplot(aes(phylogroup,abslat)) +
  geom_violin(aes(fill=phylogroup),color=NA,alpha=0.5,width=20)+
  geom_jitter(aes(color=phylogroup))+
  scale_color_manual(values=phylogroup.colors) +
  scale_fill_manual(values=phylogroup.colors) +
  fdb_style(aspect.ratio=1.5)




#------------------------------------#

CoordinateRange <- metadata %>% tibble::as.tibble() %>%
  filter(!is.na(lon)) %>%
  group_by(phylogroup) %>%
  summarize(minLat = min(lat), minLatName = genome[which.min(lat)],
            maxLat = max(lat), maxLatName = genome[which.max(lat)],
            minLon = min(lon), minLatName = genome[which.min(lon)],
            maxLon = max(lon), maxLatName = genome[which.max(lon)]) %>%
  mutate(rangeLat = maxLat-minLat, rangeLon = maxLon-minLon)

CoordinateRange %>%
  ggplot(aes(phylogroup,rangeLon)) +
  geom_col()



metadata %>% as.tibble() %>% group_by(SS) %>% count()
metadata %>% as.tibble() %>% group_by(SS) %>% count()

#Filter out the low frequency habitat classes
filtered_on_habitat <- metadata %>%
  tibble::as.tibble() %>%
  dplyr::filter(SS !="") %>%
  dplyr::group_by(SS) %>%
  dplyr::filter(n() > 5L) #%>%
  #dplyr::count()

#Filter out the low frequency habitat classes
filtered_on_SC <- metadata %>%
  tibble::as.tibble() %>%
  dplyr::filter(sourceClass !="") %>%
  dplyr::group_by(sourceClass) %>%
  dplyr::filter(n() > 5L) #%>%
#dplyr::count()


filtered_on_group <- metadata %>%
  tibble::as.tibble() %>%
  dplyr::filter(is.na(phylogroup)) %>%
  dplyr::group_by(phylogroup) %>%
  dplyr::filter(n() > 5L)#%>%
  #dplyr::count()


Filtered_on_both <-filtered_on_habitat %>% filter(genome %in% filtered_on_group$genome)




ggplot(Filtered_on_both, aes(x= SS, y = selgene, fill=SS)) +
  #xlim(c(round(min(df.chkm$GC)) - 2, round(max(df.chkm$GC)) + 2)) +
 # geom_hline(yintercept = mean(df.chkm$GC), size = 0.25, colour = '#bdbdbd') +
  geom_boxplot(alpha=0.3,outlier.size=-1)+
  geom_jitter(height=0,width=0.2,size=2,shape=21,alpha=0.7) +
  ylab('GC-content (%)') +
  xlab('clade') +
  theme_classic() +
  theme(
    aspect.ratio = 1/2,
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    #legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 10),
    #legend.text.align = 1,
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))


ggplot(filtered_on_habitat,aes(Genome_size, Coding_density,fill= SS)) +
  xlab('Genome size (bp)') +
  ylab('Coding density (%)') +
  geom_point(shape = 21, size = 2) +
  stat_ellipse(aes(color=SS),level=0.75) +
  theme_classic() +
  theme(
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


#-----------------------------------


#=========================================================================
# GEOGRAPHICAL DISTANCE
#=========================================================================


geo.dist <- geosphere::distm(dat2[,c('lon','lat')])
rownames(geo.dist) <- dat2$genome
colnames(geo.dist) <- dat2$genome

dfc <- cbind('geo'=c(geo.dist), 'ani'=c(as.matrix(ANIb.f[rownames(geo.dist),rownames(geo.dist)])))
dfc %>% as_tibble() %>% ggplot(aes(x=geo,y=ani))+geom_point()+fdb_style()


dfc %>% as_tibble() %>% ggplot(aes(x=geo,y=ani))+geom_hex(bins = 35, colour = NA) +fdb_style()

dfc %>% as_tibble() %>% ggplot(aes(x=geo,y=ani))+geom_hex(bins = 100, colour = NA) +fdb_style()


dfc %>% as_tibble() %>% ggplot(aes(x=geo/1000,y=ani))+ geom_point(aes(color=ani)) +fdb_style()


dfc <- dfc %>% data.frame() %>% mutate(ani_group = ifelse(ani>0.98,'0.98',
                                  ifelse(ani>0.95,'0.95',
                                         ifelse(ani>0.90,'0.90',
                                                ifelse(ani>0.85,'0.85',
                                                       ifelse(ani>0.80,'0.80',
                                                              ifelse(ani>0.75,'0.75','0.7'))
                                                       )))))



p.ani_geo <- dfc %>%
  tibble::as_tibble() %>%
  dplyr::filter(ani<1) %>%
  ggplot2::ggplot(ggplot2::aes(x=geo/1000,y=ani))+
  ggplot2::geom_point(shape=21,ggplot2::aes(color=ani)) +
  viridis::scale_color_viridis(limits=c(0.7,1))+
  ggplot2::ylab('ANIb') +
  ggplot2::xlab('Geographical distance (km)') +
  fdb_style()

ggplot2::ggsave('~/DATA/MarbGenomics/Graphs/ani_vs_geographical_dist.pdf',plot=p.ani_geo, width = 4, height = 4,unit='in')





df2 %>% filter(ANIb>0.95) %>%
  ggplot2::ggplot(aes(x=COPH, y=GRR)) +
  ggplot2::geom_point(aes(x=COPH, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',shape=21) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=(length.bias - b.min) / b.diff * a.diff + a.min),color='grey',se=F) +
  ggplot2::geom_point(shape=21,aes(color=ANIb)) +
  ggplot2::geom_smooth(data=modeldat,aes(x=COPH, y=GRR),color='black',se=F,size=1.5) +
  scale_y_continuous(sec.axis = sec_axis( ~((. -a.min) * b.diff / a.diff) + b.min,
                                          name = "Length Bias"))+
  fdb_style(aspect.ratio=1)+ ylab('GRR') + viridis::scale_color_viridis(limits=c(0.7,1))



commsimi <- ANIb[rownames(geo.dist),colnames(geo.dist)]
geodis <- geo.dist
geoscale=c(0,Inf)
nperm=20

###--- distance decay
disdecay<-function(commsimi, geodis, geoscale=c(0,Inf),nperm=999){
  Y<-commsimi[row(commsimi)>col(commsimi)]
  X<-geodis[row(geodis)>col(geodis)]
  minY<-min(Y[Y>0])
  Y[Y<=0]<-minY
  logY<-log10(Y)
  logX<-log10(X+1)
  logY1<-logY[X>=geoscale[1] & X<geoscale[2]]
  logX1<-logX[X>=geoscale[1] & X<geoscale[2]]

  logX1<-cbind(rep(1,length(logX1)),logX1)

  XX <- crossprod(logX1)
  XX <- solve(XX)

  # will need to calculate Xy for each permutation
  XY <- crossprod(logX1, logY1)

  # regression coefficients
  b <- XX %*% XY
  #slop
  slop.obs<-b[2]

  perm.slope<-sapply(1:nperm,function(i){
    newY <- ecodist::full(logY)
    newSample <- sample(nrow(newY))
    newY <- newY[newSample, newSample]
    logYnew <- ecodist::lower(newY)
    logYnew<-logYnew[X>=geoscale[1] & X<geoscale[2]]

    # will need to calculate Xy for each permutation
    XY <- crossprod(logX1, logYnew)

    # regression coefficients
    b <- XX %*% XY
    #slop
    b[2]
  })

  permSlope.sd<-sd(perm.slope)
  slope.rank<-rank(c(slop.obs,perm.slope),ties.method="first")[1]
  ttest.pvalue<-t.test(perm.slope,mu=slop.obs)$p.value

  c(slope.obs=slop.obs,permSlp.mean=mean(perm.slope),permSlope.sd=permSlope.sd,slope.rank=slope.rank,nperm=nperm,ttest.pvalue=ttest.pvalue)
}








MDS = cmdscale(geo.dist, eig = TRUE)
MDSdata = data.frame(MDS$points)
MDS.var.perc <- round(MDS$eig/sum(MDS$eig)*100,1)
MDSdata <- MDSdata %>%
  rownames_to_column(var='genome') %>%
  left_join(metadata, by='genome')


p = ggplot(MDSdata,aes(x=X1,y=X2,fill=abs(lat),label=genome)) +
  theme_classic() +
  xlab(paste("PCoA1 (", MDS.var.perc[1], "%)",sep="")) +
  ylab(paste("PCoA2 (", MDS.var.perc[2], "%)",sep="")) +
  geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
  geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
  geom_point(shape = 21, size = 2) +
  fdb_style()

p


#------------------------------------------------------#
# PLOT MAPS
#-------------------------------------------------------#
# FLAT MAP, GLOBAL, EASY

library(ggplot2)
library(dplyr)
require(maps)
require(viridis)


#' Title
#'
#' @param metadata
#' @param color.value
#' @param discrete
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
plot.world <- function(metadata,color.value='abs(lat)',discrete=TRUE,colors=NULL){
  world <- ggplot2::map_data("world")

  p.w <-  ggplot() +
      geom_map(
        data = world, map = world,
        aes(long, lat, map_id = region),
        color = "grey70",
        fill = "gray90", size = 0.2
      ) +
      ylim(c(-90,90)) +
      xlim(c(-180,180)) +
    geom_hline(yintercept=0,size = 0.25,linetype='dashed',color='grey70') +
    geom_hline(yintercept=60,size = 0.25,linetype='dashed',color='grey70') +
    geom_hline(yintercept=-60,size = 0.25,linetype='dashed',color='grey70') +
    geom_point(
        data = metadata,
        aes_string('lon','lat',fill= color.value),
        #alpha = 0.7,
        shape = 21, size = 2.5
      )+ theme_void()+
      #geom_hline(yintercept = 0, size = 0.25, colour = "grey70") +
      #geom_hline(yintercept=66.5,size = 0.25,linetype='dashed',color='deepskyblue3') +

    scale_y_continuous(breaks=c(-60,-30,0,30,60))+
      theme(
        aspect.ratio = 0.5,
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 10, colour = '#000000'),
        legend.justification = c(1, 1),
        legend.key.width = unit(0.25, 'cm'),
        legend.key.height = unit(0.55, 'cm'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))+

    if(discrete==FALSE){
      p.w <- p.w + scale_fill_gradientn(colours = brewer.pal(11, 'RdBu'), guide = F)
    }
    if(!is.null(colors)){
      p.w <- p.w + scale_fill_manual( values = colors)
    }

  return(p.w)
}

#p.wm.lat <- plot.world(metadata,color.value = 'abs(lat)',discrete = FALSE)
#ggsave('~/DATA/MarbGenomics/Graphs/worlmap_absolute_latitude.pdf',plot=p.wm.lat , width = 7, height = 3,unit='in')

p.wm.ani <- plot.world(metadata,color.value = 'phylogroup',discrete = TRUE,phylogroup.colors)
ggsave('~/DATA/MarbGenomics/Graphs/worlmap_phylogroup_clique.pdf',plot=p.wm.ani , width = 7, height = 3,unit='in')

p.wm.habitat <- plot.world(metadata,color.value = 'SS',discrete = TRUE)
ggsave('~/DATA/MarbGenomics/Graphs/worlmap_habitat.pdf',plot=p.wm.habitat , width = 7, height = 3,unit='in')


p.wm.ani <- p.wm.ani + coord_fixed(1.3) + theme_void() + theme(legend.position = 'none')
p.comb <- grid.arrange(arrangeGrob(p.wm.ani, plat,nrow=1, widths=c(1.5,1)))
ggsave('~/DATA/MarbGenomics/Graphs/worlmap_habitat_with_latplot.pdf',plot=p.comb , width = 8, height = 3,unit='in')

#----
# WITH COORD MAP

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0.1
  ) +
  ylim(c(-100,100)) +
  xlim(c(-180,180)) +
  geom_point(
    data = metadata,
    aes(lon,lat,fill= SS),
    alpha = 0.7,
    shape = 21, size = 2
  )+ theme_void()+
  theme(
    aspect.ratio = 0.5,
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    #legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 10),
    #legend.text.align = 1,
    legend.title = element_text(size = 11))+
  coord_map()


#

#------------------------------------------------------#
# ON A SPHERE
#-------------------------------------------------------#
# FLAT MAP, GLOBAL, EASY

#----------------------------
### Getting basemap shapefile
library(ggplot2)
library(rgdal)
library(ggmap)
library(sp)
library(dplyr)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#### Breaking apart all the values
x = c(-21.5,19.0,-161.5,-147.5)
y = c(70.5,74.5,58.5,60.5)
Type =c("A","B","A","B")

### Creating spatial LAT LONG coordinates, which will be converted to Lambert Conformal Conic Projection below
dat <- data.frame(lon = x, lat = y)
dat <- metadata %>% select(lon,lat) %>% data.frame() %>% filter(!is.na(lon)) %>% filter(lat>50)

#### Creating LAT LONG SpatialPoints
coordinates(dat) = c("lon", "lat")
proj4string(dat) <- CRS("+init=epsg:4326")

#### The coordinate reference system, that is used in your shapefile. Will use this when converting the spatial points
polar = "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
b <- bbox(dat)

### Converting the LAT LONG to the polar CRS shown above
polar_dat = spTransform(dat, polar)
polar_dat = as.data.frame(polar_dat)

#### Adding the Type column back to the data frame, with the new polar coordinates

polar_dat = data.frame(polar_dat, 'group' = metadata[rownames(polar_dat),'phylogroup'])





spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)
ggplot(wmap_df_robin, aes(long,lat, group=group)) +
  geom_polygon() +
  labs(title="World map (robinson)") +
  coord_equal() +
  theme_opts

#----------
#install.packages(c('mapview','lwgeom'),dependencies=TRUE)
library(sf)
library(ggplot2)
library(mapview)
library(lwgeom)
library(rnaturalearth)

# world data
world <- rnaturalearth::ne_countries(scale = 'small', returnclass = 'sf')

# Fix polygons so they don't get cut in ortho projection
world  <- st_cast(world, 'MULTILINESTRING') %>%
  st_cast('LINESTRING', do_split=TRUE) %>%
  mutate(npts = mapview::npts(geometry, by_feature = TRUE)) %>%
  st_cast('POLYGON')


#===
spTransform(world, CRS("+proj=robin"))



#---------------------------------------------#
# NORTH POLE

dat <- metadata %>% dplyr::select(lon,lat,phylogroup) %>% data.frame() %>% filter(!is.na(lon)) %>% filter(lat>0)

events_sf <- dat %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
events_sf

crs_val = "+proj=ortho +lat_0=90 +lon_0=0"

events_transformed <- events_sf %>%
  st_transform(crs = crs_val)

# The coordinates are different now!
events_transformed
st_coordinates(events_transformed)

# Bind the coordinates to the full data frame
events_transformed_with_lat_lon <- cbind(events_transformed, st_coordinates(events_transformed))

# The X and Y points are now extracted!
events_transformed_with_lat_lon

north.p <- ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  # geom_sf(data = whatever_shapefile_youre_using) +
  geom_point(data = events_transformed_with_lat_lon, aes(x = X, y = Y, fill=phylogroup), size = 3, shape=21) +
  coord_sf(crs = crs_val)+
  theme_minimal() +
  scale_fill_manual(values=phylogroup.colors)


#--------------------------------------------#
# SOUTH POLE


dat <- metadata %>% dplyr::select(lon,lat,phylogroup) %>% data.frame() %>% filter(!is.na(lon)) #%>% filter(lat<0)

events_sf <- dat %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
events_sf

crs_val = "+proj=ortho +lat_0=-90 +lon_0=0"

events_transformed <- events_sf %>%
  st_transform(crs = crs_val)

# The coordinates are different now!
events_transformed
st_coordinates(events_transformed)

# Bind the coordinates to the full data frame
events_transformed_with_lat_lon <- cbind(events_transformed, st_coordinates(events_transformed))

# The X and Y points are now extracted!
events_transformed_with_lat_lon

#drop empty points
events_transformed_with_lat_lon <- events_transformed_with_lat_lon[!st_is_empty(events_transformed_with_lat_lon),,drop=FALSE]

south.p <- ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  geom_point(data = events_transformed_with_lat_lon, aes(x = X, y = Y, fill=phylogroup), size = 3, shape=21) +
  coord_sf(crs = crs_val)+theme_minimal()+
  scale_fill_manual(values=phylogroup.colors)

gridExtra::grid.arrange(north.p,south.p,nrow=1)

#-----


##### FINAL
# THE REMOVE POINT THING ALLOWS NOT TO USE LAT FILTER <0
dat <- metadata %>% dplyr::select(lon,lat,phylogroup) %>% data.frame() %>% filter(!is.na(lon)) #%>% filter(lat>0)

#dat <- metadata %>% dplyr::select(lon,lat,group) %>% data.frame() %>% filter(!is.na(lon)) %>% filter(lat>40)

events_sf <- dat %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
events_sf

crs_val = "+proj=ortho +lat_0=90 +lon_0=0"
crs_val = "+proj=ortho +lat_0=-90 +lon_0=0"

crs_val = "+proj=ortho +lat_0=0 +lon_0=0"
#for US
crs_val = "+proj=ortho +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0"

#For China
crs_val = "+proj=ortho +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=96 +x_0=0 +y_0=0"
crs_val = "+proj=ortho +lat_1=-50 +lat_2=-40 +lat_0=-70 +lon_0=0 +x_0=0 +y_0=0"

#crs_val = "+proj=ortho +lat_1=-50 +lat_2=-70 +lon_0=0 +x_0=0 +y_0=0"


events_transformed <- events_sf %>%
  st_transform(crs = crs_val)

# The coordinates are different now!
#events_transformed
#st_coordinates(events_transformed)

# Bind the coordinates to the full data frame
events_transformed_with_lat_lon <- cbind(events_transformed, st_coordinates(events_transformed))

# The X and Y points are now extracted!
#events_transformed_with_lat_lon

#drop empty points
events_transformed_with_lat_lon <- events_transformed_with_lat_lon[!st_is_empty(events_transformed_with_lat_lon),,drop=FALSE]

anlge.p.1 <- ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  geom_point(data = events_transformed_with_lat_lon, aes(x = X, y = Y, fill=phylogroup), size = 3, shape=21) +
  coord_sf(crs = crs_val)+theme_minimal()+  scale_fill_manual(values=phylogroup.colors)






#------------------------------------------#

#Angle 1 and angle 2
#crs_val = "+proj=ortho +lat_1=60 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0"
crs_val = "+proj=ortho +lat_1=80 +lat_2=90 +lat_0=70 +lon_0=-45 +x_0=0 +y_0=0"

events_transformed <- events_sf %>%
  st_transform(crs = crs_val)

events_transformed_with_lat_lon <- cbind(events_transformed, st_coordinates(events_transformed))
events_transformed_with_lat_lon <- events_transformed_with_lat_lon[!st_is_empty(events_transformed_with_lat_lon),,drop=FALSE]

anlge.p.1 <- ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  geom_point(data = events_transformed_with_lat_lon, aes(x = X, y = Y, fill=phylogroup), size = 3, shape=21) +
  coord_sf(crs = crs_val)+theme_minimal()  + scale_fill_manual(values=phylogroup.colors)+theme(legend.position='none')




#
crs_val = "+proj=ortho +lat_1=-50 +lat_2=-40 +lat_0=-70 +lon_0=96 +x_0=0 +y_0=0"

events_transformed <- events_sf %>%
  st_transform(crs = crs_val)

events_transformed_with_lat_lon <- cbind(events_transformed, st_coordinates(events_transformed))
events_transformed_with_lat_lon <- events_transformed_with_lat_lon[!st_is_empty(events_transformed_with_lat_lon),,drop=FALSE]

anlge.p.2 <- ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  geom_point(data = events_transformed_with_lat_lon, aes(x = X, y = Y, fill=phylogroup), size = 3, shape=21) +
  coord_sf(crs = crs_val)+theme_minimal() + scale_fill_manual(values=phylogroup.colors)+theme(legend.position='none')

p.comp <- gridExtra::grid.arrange(anlge.p.1,anlge.p.2,nrow=1)
ggsave('~/DATA/MarbGenomics/Graphs/worlmap_rotated_hemispheres_habitat.pdf',plot=p.comp , width = 7, height = 3,unit='in')





#-----------------
#Show specific region of the earth, say China sea

ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  #geom_text(data = world, aes(X, Y, label = name), size = 5) +
  geom_point(data = dat, aes(x = lon, y = lat, fill=phylogroup), size = 3, shape=21) +
  coord_sf(xlim = c(80, 160), ylim = c(-15, 50), expand = FALSE)+theme_minimal()+ scale_fill_manual(values=phylogroup.colors)+theme(legend.position='none')


ggplot() +
  geom_sf(data=world, color="gray80",size=0.5) +
  #geom_text(data = world, aes(X, Y, label = name), size = 5) +
  geom_point(data = dat, aes(x = lon, y = lat, fill=phylogroup), size = 3, shape=21) +
  coord_sf(xlim = c(115, 135), ylim = c(28, 41), expand = FALSE)+theme_minimal()+ scale_fill_manual(values=phylogroup.colors)+theme(legend.position='none')



#============================

#######################
#http://rpsychologist.com/working-with-shapefiles-projections-and-world-maps-in-ggplot
library(rgdal)
library(ggplot2)

setwd('~')

# read shapefile
wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")
# convert to dataframe
wmap_df <- fortify(wmap)

# create a blank ggplot theme
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill="#e6e8ed"),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=22)))

# plot map
ggplot(wmap_df, aes(long,lat, group=group)) +
  geom_polygon() +
  labs(title="World map (longlat)") +
  coord_equal() +
  theme_opts

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)
ggplot(wmap_df_robin, aes(long,lat, group=group)) +
  geom_polygon() +
  labs(title="World map (robinson)") +
  coord_equal() +
  theme_opts

ggplot(wmap_df_robin, aes(long,lat, group=group, fill=hole)) +
  geom_polygon() +
  labs(title="World map (Robinson)") +
  coord_equal() +
  theme_opts +
  scale_fill_manual(values=c("#262626", "#e6e8ed"), guide="none") # change colors & remove legend

# add graticule and bounding box (longlat)
grat <- rgdal::readOGR("ne_110m_graticules_all", layer="ne_110m_graticules_15")
grat_df <- fortify(grat)

bbox <- readOGR("ne_110m_graticules_all", layer="ne_110m_wgs84_bounding_box")
bbox_df<- fortify(bbox)

ggplot(bbox_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_df, aes(long,lat, group=group, fill=hole)) +
  geom_path(data=grat_df, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map + graticule (longlat)") +
  coord_equal() +
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none") # change colors & remove legend


grat_robin <- spTransform(grat, CRS("+proj=robin"))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <- spTransform(bbox, CRS("+proj=robin"))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_df_robin, aes(long,lat, group=group, fill=hole)) +
  #geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map (Robinson)") +
  coord_equal() +
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none") # change colors & remove legend

countries <- readOGR("ne_110m_admin_0_countries", layer="ne_110m_admin_0_countries")
countries_robin <- spTransform(countries, CRS("+init=ESRI:54030"))
countries_robin_df <- fortify(countries_robin)

places_robin_df <- project(cbind(metadata$lon, metadata$lat), proj="+init=ESRI:54030")
places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("lon", "lat")
places_robin_df$SS <- metadata$SS
places_robin_df$phylogroup <- metadata$phylogroup

ggmap = ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole)) +
  geom_path(data=countries_robin_df, aes(long,lat, group=group, fill=hole), color="white", size=0.3) +
  #  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  coord_equal() +
  theme_opts   +
  scale_fill_manual(values=c("lightgray", "white"), guide="none") # change colors & remove legend

ggmap +    geom_point(data=places_robin_df, aes(LONGITUDE, LATITUDE, group=SS, colour=SS,size=2)) +scale_colour_manual(values=Habitat_colors)





#====================================#

SRA_merged <-
  sra.coords  %>%
  dplyr::left_join(ot, by="SRA") %>%
  dplyr::left_join(imngs.RA.l, by=c('SRA'='Sample_Name')) %>% filter(value>0)

SRA_projected <- SRA_merged
SRA_projected[,c('lon','lat')] <- project(cbind(SRA_projected$lon, SRA_projected$lat), proj="+init=ESRI:54030")
SRA_projected$abslat = abs(SRA_merged$lat)

imngs.RA.l %>% filter(value!=0) %>% select(Sample_Name) %>% unique()
imngs.RA.l %>% filter(value>0.0001) %>% select(Sample_Name) %>% unique()


gplt <- ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole),alpha=0.7) +
  geom_path(data=countries_robin_df, aes(long,lat, group=group, fill=hole), color="white", size=0.3) +
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="lightgrey",size=0.3) +
  scale_fill_manual(values=c("grey", "white"), guide="none")+
  geom_point(data=SRA_projected %>% filter(value>0.00001) %>%
               select(lon,lat,abslat) %>%
               unique(),
             aes(lon,lat, color=abslat),size=1,shape=4,alpha=0.7) +
  geom_hline(yintercept=0,color='darkgrey',linetype='dashed')+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdYlBu')) +
  ggnewscale::new_scale_color()+
  ggnewscale::new_scale_fill()+
  geom_point(data=places_robin_df, aes(lon,lat, fill=phylogroup),shape=21,size=3) +
  coord_equal() +
  theme_opts   +
  scale_fill_manual(values=phylogroup.colors)# change colors & remove legend

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/globe_test.pdf',plot=gplt, width = 8, height = 4,unit='in')






gplt <- ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole),alpha=0.7) +
  geom_polygon(data=bbox_robin_df, aes(long,lat),color='darkgrey',fill='transparent',size=0.3) +
  #geom_path(data=countries_robin_df, aes(long,lat, group=group, fill=hole), color="white", size=0.5) +
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed",color="lightgrey",size=0.3) +
  scale_fill_manual(values=c("lightgrey", "grey20"), guide="none")+
  geom_point(data=SRA_projected %>% filter(value>0.00001) %>%
               select(lon,lat,abslat) %>%
               unique(),
             aes(lon,lat),size=0.5,shape=4,alpha=0.7) +
  #geom_hline(yintercept=0,color='grey',linetype='dashed')+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdYlBu')) +
  ggnewscale::new_scale_color()+
  ggnewscale::new_scale_fill()+
  geom_point(data=places_robin_df, aes(lon,lat, fill=phylogroup),color='black',shape=21,size=3) +
  coord_equal() +
  theme_opts   +
  theme(plot.background = element_rect(fill="white"))+
  scale_fill_manual(values=phylogroup.colors)+theme(legend.position = 'none')# change colors & remove legend

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/globe_test2.pdf',plot=gplt, width = 8, height = 4,unit='in')



#---- grid line labels
lbl.Y <- data.frame(lon=rep(c(-180,180),each=13),
                    lat=rep(seq(from=90, to=-90, by=-15), times=2))

lbl.Y$dir <- ifelse(lbl.Y$lat ==0,'',ifelse(lbl.Y$lat>0,'N','S'))
lbl.Y$lbl <- paste0(abs(lbl.Y$lat), lbl.Y$dir)

lbl.X <- data.frame(lon=rep(seq(160,-160,by=-20),times=2),
                    lat=rep(c(from=-85, to=85), each=17))

lbl.X$dir <- ifelse(lbl.X$lon ==0,'',ifelse(lbl.X$lon>0,'E','W'))
lbl.X$lbl <- paste0(abs(lbl.X$lon), lbl.Y$dir)


prj.coord <- project(cbind(lbl.Y$lon, lbl.Y$lat), proj="+init=ESRI:54030")
lbl.Y.prj <- cbind(prj.coord, lbl.Y)
names(lbl.Y.prj)[1:2] <- c('X.prj','Y.prj')

prj.coord <- project(cbind(lbl.X$lon, lbl.X$lat), proj="+init=ESRI:54030")
lbl.X.prj <- cbind(prj.coord, lbl.X)
names(lbl.X.prj)[1:2] <- c('X.prj','Y.prj')


gplt <- ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole),alpha=0.7) +
  geom_polygon(data=bbox_robin_df, aes(long,lat),color='darkgrey',fill='transparent',size=0.3) +
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed",color="lightgrey",size=0.3) +
  scale_fill_manual(values=c("lightgrey", "grey20"), guide="none")+
  geom_point(data=SRA_projected %>% filter(value>0.00001) %>%
               select(lon,lat,abslat) %>%
               unique(),
             aes(lon,lat),size=0.5,shape=4,alpha=0.7) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdYlBu')) +
  ggnewscale::new_scale_color()+
  ggnewscale::new_scale_fill()+
  geom_point(data=places_robin_df, aes(lon,lat, fill=phylogroup),color='black',shape=21,size=3) +
  #geom_text(data=lbl.X.prj, aes(x=X.prj, y=Y.prj,label=lbl),color='grey50',size=2) +
  geom_text(data=lbl.Y.prj, aes(x=X.prj, y=Y.prj,label=lbl),color='grey50',size=2) +

  coord_equal() +
  theme_opts   +
  theme(plot.background = element_rect(fill="white"))+
  scale_fill_manual(values=phylogroup.colors)+theme(legend.position = 'none')# change colors & remove legend

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/globe_test3.pdf',plot=gplt, width = 8, height = 4,unit='in')







gplt <- ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole),alpha=0.7) +
  geom_polygon(data=bbox_robin_df, aes(long,lat),color='darkgrey',fill='transparent',size=0.3) +
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed",color="lightgrey",size=0.3) +
  scale_fill_manual(values=c("lightgrey", "grey20"), guide="none")+
  geom_point(data=SRA_projected %>% filter(value>0.0001) %>%
               select(lon,lat,abslat) %>%
               unique(),
             aes(lon,lat),size=0.5,shape=4,alpha=0.7) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdYlBu')) +
  ggnewscale::new_scale_color()+
  ggnewscale::new_scale_fill()+
  geom_point(data=places_robin_df, aes(lon,lat, fill=phylogroup),color='black',shape=21,size=3) +
  #geom_text(data=lbl.X.prj, aes(x=X.prj, y=Y.prj,label=lbl),color='grey50',size=2) +
  geom_text(data=lbl.Y.prj, aes(x=X.prj, y=Y.prj,label=lbl),color='grey50',size=2) +

  coord_equal() +
  theme_opts   +
  theme(plot.background = element_rect(fill="white"))+
  scale_fill_manual(values=phylogroup.colors)+theme(legend.position = 'none')# change colors & remove legend

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/globe_test_0.0001.pdf',plot=gplt, width = 8, height = 4,unit='in')











top_titles <- SRA_merged %>% filter(value>0.0001) %>%
  group_by(Study_Title,Bioproject) %>%
  summarise(mean=mean(value),sum=sum(value)) %>%
  arrange(desc(sum)) %>%
  data.frame() %>% head(75) %>% select(Study_Title) %>% pull() %>% as.character()


#SRA_merged %>%
#  group_by(Study_Title,Bioproject) %>% filter(value>0.0001) %>%
#  data.frame() %>% ggplot(aes(value,Study_Title))+geom_col()



SRA_merged %>%
  group_by(Study_Title,Bioproject) %>% filter(value>0.0001) %>%
  filter(Study_Title %in% top_titles) %>%
  data.frame() %>% ggplot(aes(value,reorder(Study_Title,value)))+geom_jitter(aes(color=Description))


SRA_merged %>% filter(value>0.0001) %>% group_by(Description) %>% tally() %>% arrange(desc(n))



p.host_as <- SRA_merged %>%
  dplyr::mutate(Study_Title = gsub('>','',gsub('</','',Study_Title))) %>%
  dplyr::mutate(Study_Title = paste0(Study_Title,'- ',Bioproject)) %>%
  dplyr::filter(grepl(paste(
    c('dinoflagellate','alga','plant','rhizosphere','sea stars','coral','fish','sea anemone','epibiont','endophyte','sponge','oyster','zebrafish','shrimp'),collapse='|'),
    Description,ignore.case=TRUE)) %>%
  dplyr::filter(value>0.0001) %>%
  data.frame() %>%
  ggplot2::ggplot(ggplot2::aes(value,reorder(Study_Title,value)))+
  ggplot2::geom_jitter(ggplot2::aes(color=Description),width = 0.1,shape=21)+
  ggplot2::facet_grid(Description~., space='free', scales='free_y') +
  ggplot2::theme(strip.text.y=element_text(angle=0))+
  ggplot2::scale_x_log10() +
  ggplot2::xlab('Relative Abundance (%)') +
  ggplot2::ylab('') +
  scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,'Blues'))(30)[15:30])+
  ggplot2::theme(legend.position = 'none')

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/SRA_host_associated2.pdf',
                plot=p.host_as,
                width = 360, height = 360, unit='mm')




SRA_merged %>%
  filter(grepl(paste(
    c('Tibet','Lake','altitu','mountain'),collapse='|'),
    Description,ignore.case=TRUE)) %>% filter(value>0.0001) %>%
  #  filter(Study_Title %in% top_titles) %>%
  data.frame() %>% ggplot(aes(value,reorder(Study_Title,value)))+geom_jitter(aes(color=Description))+
  facet_grid(Description~., space='free', scales='free_y') + theme(strip.text.y=element_text(angle=0))+
  #coord_trans(x='log10')
  scale_x_log10()


#---------- random stuff?
#1  PRJNA385736  5491 Marine metagenomes Metagenome</
#2  PRJNA429259  3587 EEPEND: Microbiome and bacterioplankton rRNA gene sequence data collected from Gulf of Mexico seawater samples. Cruises DP03 and DP04 from Jan 2016 - December 2016
#3  PRJNA414409  1707 ADDOME
#4  PRJNA284506  1527 Aquarium Targeted loci environmental</
#5  PRJNA419719  1130 ADDOMEx Tier 3 Experiments: Mesocosm Si with Gulf of Mexico coastal waters</
#6  PRJNA320765  1069 ADDOMEx Tier 3 Experiments: Mesocosm #1 and #2
#7  PRJNA347350   955 >Warming effect on Zostera Marina metagenome</
#8  PRJNA522679   904 hycosphere Bacteria Diversity of Pyropia with yellow sport disease Raw sequence reads</
#9    PRJDB3419   802 >Bacterial Biogeography of Coastal Water in the Northern East China Sea</
#10 PRJNA325151   598 >seawater metagenome Targeted loci environmental</


SRA_merged %>%
  group_by(Study_Title,Bioproject) %>% filter(value>0.0001) %>%
  filter(Study_Title %in% top_titles) %>%
  data.frame() %>% ggplot(aes(value,reorder(Study_Title,value)))+geom_jitter(aes(color=Description))+facet_wrap(~Description,ncol=1, scales='free_y')



SRA_merged %>%
  group_by(Study_Title,Bioproject)%>% filter(value>0.0001) %>%
  filter(Study_Title %in% top_titles) %>%
  #filter(value>0.0001) %>%
  ggplot(aes(x=value))+geom_density()





SRA_merged %>% select(-assembl,-ph, -name) %>%
  arrange(desc(value)) %>%
  group_by(SRA,lat,lon,Title,Study_Title,Bioproject,value,Description) %>%
  summarise(max=max(value))





SRA_merged %>% select(lon,lat) %>% unique() %>% ggplot(aes(x=lat))+geom_density()

SRA_merged %>%
  group_by(Study_Title,Bioproject) %>% arrange(desc(value)) %>% data.frame()




outdir

SRA_merged
write.table(SRA_merged, file = paste0(outdir,'/SRA_merged_all.txt'), sep='\t', row.names=FALSE, quote = FALSE)





sra.coords  %>%
  dplyr::left_join(ot, by="SRA") %>%
  dplyr::left_join(imngs.RA.l, by=c('SRA'='Sample_Name')) %>% filter(value>0) %>%
  #filter(assembl %in% single_ot) %>%
  #dplyr::filter(SRA %in%  selected_SRA) %>%
  ggplot2::ggplot() +
  mapWorld +
  ggplot2::geom_point(ggplot2::aes(x=lon, y=lat,color=ph,size=value),alpha=0.3) +
  ggplot2::geom_point(ggplot2::aes(x=lon, y=lat,color=ph,size=value),shape=21) +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::coord_quickmap()+facet_wrap(~ph)







ggplot(bbox_robin_df, aes(long,lat, group=group)) +
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole),alpha=0.7) +
  geom_path(data=countries_robin_df, aes(long,lat, group=group, fill=hole), color="white", size=0.3) +
  #  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  scale_fill_manual(values=c("grey", "white"), guide="none")+
  geom_point(data=SRA_projected %>% filter(value>0.00001) %>%
               select(lon,lat,abslat) %>%
               unique(),
             aes(lon,lat),size=1,shape=4,alpha=0.7) +
  geom_hline(yintercept=0,color='darkgrey',linetype='dashed')+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, 'RdYlBu')) +
  ggnewscale::new_scale_color()+
  ggnewscale::new_scale_fill()+
  geom_point(data=places_robin_df, aes(lon,lat, fill=phylogroup),shape=21,size=3) +
  coord_equal() +
  theme_opts   +
  scale_fill_manual(values=phylogroup.colors)# change colors & remove legend



#just some analysis


#dfc <- cbind('geo'=c(geo.dist), 'ani'=c(dstANI[rownames(geo.dist),rownames(geo.dist)]))
#gplot(dfc,aes(x=geo,y=ani))+geom_point()+fdb_style()
#ggplot(dfc,aes(x=geo,y=ani))+geom_hex(bins = 35, colour = NA) +fdb_style()


#Let's do some modeling!?

#augmented bootstrap shitness
# of GML, NLS you name it


data <- genome.tbl %>% filter(genus=='Marinobacter') %>% filter(quality_class=='high') %>% mutate(distequator = abs(lat)) %>% left_join(bias_data,by='genome')
data$b <- as.numeric(as.character(data$b))
y <- 'Coding_density'
x <- 'Genome_size'


y <- 'GC'
x <- 'distequator'

y <- 'Coding_density'
x <- 'distequator'

y <- 'Coding_density'
x <- 'GC'

y <- 'b'
x <- 'GC'

y <- 'b'
x <- 'distequator'

df.distances$


#=========
#exponential decay
#nlsfit <- nls(ANIb ~ SSasymp(COPH, yf, y0, log_alpha), data = data)


# having loads of data points (comparisons, we could use a smarter way, bootstrap it)

x <- 'COPH'
y <- 'ANIb'
df.distances$COPH <- 1-df.distances$COPH
data = df.distances
data = df.distances[df.distances$ANIb>0.95,]

#mytarget <- "group"
glmFormula <- as.formula(paste0(y,'~',x))
#nlsFormula <- as.formula(paste0(y, "~ k / ",x,'+ b'))
nlsFormula <- as.formula(paste0(y,' ~ SSasymp(',x,', yf, y0, log_alpha)'))


#Inspect
ggplot(data, aes_string(x=x,y=y)) +
  geom_point()

#Run models
#nlsfit <- nls(nlsFormula, data, start = list(k = 1, b = 0))
nlsfit <- nls(nlsFormula, data = data)
glmfit <- glm(glmFormula , data= data)

summary(nlsfit)
summary(glmfit)

summ.table <- do.call(rbind, lapply(list(glmfit, glmfit), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")
reported.table
summ.table <- do.call(rbind, lapply(list(nlsfit), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")
reported.table

ggplot(data, aes_string(x=x,y=y)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_hex(bins=70,size=0) +
  geom_line(aes(y = predict(nlsfit)),color='orange')+
  #geom_line(aes(y = predict(glmfit)))+
  geom_line(aes(y = predict(glmfit)),color='cyan3')+
  scale_fill_continuous(type = "viridis") +
  fdb_style()


set.seed(27)
boots <- bootstraps(data, times = 500, apparent = TRUE)
boots

fit_nls_on_bootstrap <- function(split) {
  nls(nlsFormula, analysis(split))
}

boot_models <-
  boots %>%
  mutate(model = map(splits, fit_nls_on_bootstrap),
         coef_info = map(model, tidy))


#====


fit_glm_on_bootstrap <- function(split) {
  glm(glmFormula, data= analysis(split))
}

boot_models <-
  boots %>%
  mutate(model = map(splits, fit_glm_on_bootstrap),
         coef_info = map(model, tidy))

###

boot_coefs <-
  boot_models %>%
  unnest(coef_info)

#Confidence intervals
percentile_intervals <- int_pctl(boot_models, coef_info)
percentile_intervals

#
p.btstrp.2 <- ggplot(boot_coefs, aes(estimate)) +
  geom_histogram(bins = 30,fill='grey20') +
  facet_wrap( ~ term, scales = "free") +
  geom_vline(aes(xintercept = .lower), data = percentile_intervals, col = "red") +
  geom_vline(aes(xintercept = .upper), data = percentile_intervals, col = "red") +
  scale_y_continuous(expand=c(0,0)) + fdb_style() + theme(panel.spacing=unit(0, 'cm'))

boot_aug <-
  boot_models %>%
  sample_n(200) %>%
  mutate(augmented = map(model, augment)) %>%
  unnest(augmented)

p.btstrp.3 <- ggplot(boot_aug, aes_string(x=x,y=y)) +
  #geom_smooth(data=data,aes_string(x=x,y=y),method='lm',se=FALSE,color='darkgrey') +
  geom_point(data=data,aes_string(x=x,y=y),size=2, alpha=0.5,color='grey20') +
  geom_line(aes(y = .fitted, group = id), alpha = .2, col = "red") +
  fdb_style() + ggpubr::stat_cor()

nls_out <- int_pctl(boot_models,coef_info)

nls_out





##====================================
#compare models AIC score wise











p.btstrp<- grid.arrange(p.btstrp.2, p.btstrp.3,nrow=1)
filetype <- 'pdf'
filename <- paste0('bootstrap_',y,'_',x,'.',filetype)
outdir <- "~/DATA/MarbGenomics/Graphs/"
outfile <- paste0(outdir,filename)
ggsave(outfile,plot=p.btstrp, width = 6, height = 2.5,unit='in')

filename <- paste0('bootstrap_minimal_',y,'_',x,'.',filetype)
outfile <- paste0(outdir,filename)
ggsave(outfile,plot=p.btstrp.3, width = 2.5, height = 2.5,unit='in')















fit_spline_on_bootstrap <- function(split) {
  data <- analysis(split)
  smooth.spline(data$Coding_density, data$Genome_size, df = 4)
}

boot_splines <-
  boots %>%
  sample_n(200) %>%
  mutate(spline = map(splits, fit_spline_on_bootstrap),
         aug_train = map(spline, augment))

splines_aug <-
  boot_splines %>%
  unnest(aug_train)

ggplot(splines_aug, aes(x, y)) +
  geom_line(aes(y = .fitted, group = id), alpha = 0.2, col = "blue") +
  geom_point()+fdb_style()









#------- PCA STUFF AND SHEIT

#-------------
#lets play with Dimensionality reduction in tidyverse, tidymodels etc
#visualisation, understanding the dataset

selectedVar.chkm = df.chkm %>% dplyr::select(BinId,Genome_size,group)
df.pan = cbind(selectedVar.chkm,pan[selectedVar.chkm$BinId,])

#tidymodels

#library(tidymodels)
# Unsupervised machine learning, (there is no already known outcome)
df.pca.chkm <- df.chkm %>% dplyr::select(BinId,group,Completeness:Genome_size,GC,Coding_density)

pca_rec <- recipe(~.,data=df.pca.chkm) %>%
  update_role(BinId, group, new_role = 'id' ) %>% #descide what it is...
  step_normalize(all_predictors()) %>% #  center and scale
  step_pca(all_predictors())
pca_rec


#Define the PCA recipe, update_role for desciding what are predictors and what are not
pca_rec <- recipe(~.,data=df.pan) %>%
  update_role(BinId, Genome_size, group, new_role = 'id' ) %>% #descide what it is...
  step_normalize(all_predictors()) %>% #  center and scale
  step_pca(all_predictors())
pca_rec

#prepare recipe
pca_prep <- prep(pca_rec)
pca_prep

#gives variable values per principal component
tidied_pca <- tidy(pca_prep, 2 ) #select two for the principle comonent analysis step, not step_normalise

#we can visualise some
#we have way to many to vi niceely?

tidied_pca %>% filter(component %in% paste0("PC",1:5)) %>%
  mutate(component =  fct_inorder(component)) %>%
  ggplot(aes(value, terms, fill=terms))+
    geom_col(show.legend=FALSE)+
    facet_wrap(~component, nrow=1) +
    labs(y=NULL)

#sick plot, because we can see, which variables go together, which do not, can we extract more info from this?
# see further

#lets zoom in
#What contributes the most to each component?
tidied_pca %>%
  filter(component %in% paste0("PC",1:4)) %>%
  group_by(component) %>%
  top_n(8, abs(value)) %>%
  ungroup() %>%
  ggplot(aes(abs(value), terms, fill=value >0)) +
  geom_col()+
  facet_wrap(~component, scales="free_y")

#Retrieve the analysis
juice(pca_prep)

#(why only 5 dimensions?)
juice(pca_prep) %>%
  ggplot(aes(PC1, PC2))+
  geom_point(aes(fill=group),shape=21, size=2)+ theme_classic() +
  xlab('PC 1') +
  ylab('PC 2') +
  geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  fdb_style()


# LOADINGS FUN!!!!!!
# say 5 first principal components
df.pcloadings <- tidied_pca %>%
  filter(component %in% paste0("PC",1:5)) %>%
  pivot_wider(names_from = component, values_from = value, values_fill=0)

df.pcloadings %>%  ggplot(aes(PC1, PC2))+
  geom_point(alpha=0.7, size=2)+
  theme_classic() +
  xlab('PC 1') +
  ylab('PC 2')
#umap_rec <- recipe(~.,data=df.pan) %>%
#library(embed)

#Define the PCA recipe, update_role for desciding what are predictors and what are not
umap_rec <- recipe(~.,data=df.pcloadings) %>%
  update_role(terms, id, new_role = 'id' ) %>% #descide what it is...
  step_normalize(all_predictors()) %>% #  center and scale
  step_umap(all_predictors())
umap_rec

#prepare recipe
umap_prep <- prep(umap_rec)
umap_prep

#juice(umap_prep) %>%
#  ggplot(aes(PC2, PC3))+
#  geom_point(aes(color=group),alpha=0.7, size=2)

#Retrieve the analysis
juice(umap_prep)


juice(umap_prep) %>%
  ggplot(aes(umap_1,umap_2))+
  geom_point(alpha=0.7, size=2)+
  theme_classic() +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_point(shape = 21, size = 2) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 11, colour = '#000000'),
    axis.text = element_text(size = 10, colour = '#000000'),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11))

clustered_feautures <- juice(umap_prep) %>%
  filter(umap_1 < -20) %>%
  pull(terms) %>%
  as.character

cog_desc[clustered_feautures,]


scaleRYG <- colorRampPalette(c("white","blue",'red'), space = "rgb")(30)
pan.hm <- pan
pan.hm[pan.hm == 0] <- NA

heatmap.2(pan.hm[,clustered_feautures], col=scaleRYG, margins = c(7,10), density.info = "none", trace = "none", lhei = c(2,6), xlab = "Identifier", ylab = "Rows",cexRow=0.2,na.color = "grey")
#looks a bit useless? very little consistency to me, maybe some but little
# It seems to fucking pick up the odd one again, many diupplications, look a lot like a error our thing, or is this real and does it represent a dupplication of an entire section?
# study this through prokka.tsv or alternative, to understand if it is all a adjacent region, duplication, or if it is scattered arount the genome, and why is this, look at trying to understand more

COG1995
pan.hm %>% data.frame() %>% filter(COG1995>1) %>% dplyr::select(clustered_feautures)

#---------------
#---------------

# USE UMAP INSTEAD OF PCA, genomes in pan genome space
#-----



#-----
# DIfferent unsupervised dimensionality reduction

## UMAP
#based on topology, manifold approximations

#library(embed)

#Define the PCA recipe, update_role for desciding what are predictors and what are not
umap_rec <- recipe(~.,data=df.pan) %>%
  update_role(BinId, Genome_size, group, new_role = 'id' ) %>% #descide what it is...
  step_normalize(all_predictors()) %>% #  center and scale
  step_umap(all_predictors())
umap_rec

#prepare recipe
umap_prep <- prep(umap_rec)
umap_prep

#juice(umap_prep) %>%
#  ggplot(aes(PC2, PC3))+
#  geom_point(aes(color=group),alpha=0.7, size=2)

#Retrieve the analysis
juice(umap_prep)

#(why only 5 dimensions?)
juice(umap_prep) %>%
  ggplot(aes(umap_1,umap_2))+
    geom_point(aes(color=group),alpha=0.7, size=2)+
    theme_classic() +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
    geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
    geom_point(shape = 21, size = 2) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 11, colour = '#000000'),
      axis.text = element_text(size = 10, colour = '#000000'),
      legend.justification = c(1, 1),
      legend.key.width = unit(0.25, 'cm'),
      legend.key.height = unit(0.55, 'cm'),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11))


#instead of smooth distributions, we have clumpy distibutions

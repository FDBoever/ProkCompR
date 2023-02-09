###
#COVARIANCE IN PANGENOME


# ONE CAN CLUSTER, iGRAPH OR DIMENSIONALITY REDUCE pangenomes gene wise, as oposed to per genome
# this can make us understand which gene families co-occur, co-vary and is an alternative to hypothesis testing or feature detection when using factor predictors
# This may be a nice addition to pangenome analysis, although not yet commonly used
# perhaps too far fetched for the Marinobacter genomics paper, depending what we want to say I guess..


#tsne_model_1 = Rtsne::Rtsne(as.matrix(t(OG.pan)), check_duplicates=FALSE, pca=TRUE, perplexity=25, theta=0.5, dims=2)
#d_tsne_1 = as.data.frame(tsne_model_1$Y)


# COLOR A GENE-WISE TSNE ACCORDING TO RF IMPORTANCE

pan = cazy.pan
lookup_table = RF$importance

pan = pan[,rownames(RF$importance)]

tsne_model_1 = Rtsne::Rtsne(as.matrix(t(pan)), check_duplicates=FALSE, pca=TRUE, perplexity=25, theta=0.5, dims=2)
d_tsne_1 = as.data.frame(tsne_model_1$Y)

d_tsne_1$Name = colnames(pan)
d_tsne_1$Group = colorGroups[d_tsne_1$Name]

d_tsne_1 = cbind(d_tsne_1,lookup_table[d_tsne_1$Name,])


ggplot(d_tsne_1,aes(x=V1,y=V2)) +
  theme_classic() +
  xlab(paste("tSNE1",sep="")) +
  ylab(paste("tSNE2",sep="")) +geom_point(shape=21,size=2, aes(fill=(MeanDecreaseAccuracy)))+fdb_style() + scale_fill_gradientn(colours = brewer.pal(11, 'RdBu'))


#----------------------------#
#
#   CompareM methodology
#
#-----------------------------#



df_pca <- prcomp(codon_usage)

df_pca <- prcomp(kmer_profile[intersect(rownames(kmer_profile),rownames(metadata)),])

df_out <- as.data.frame(df_pca$x)
#df_out$group <- sapply( strsplit(as.character(row.names(df)), "_"), "[[", 1 )
p<-ggplot(df_out,aes(x=PC1,y=PC2 ))
p<-p+geom_point()
p


combined <- cbind(aa_usage[intersect(rownames(aa_usage),rownames(metadata)),],
              codon_usage[intersect(rownames(aa_usage),rownames(metadata)),],
              kmer_profile[intersect(rownames(aa_usage),rownames(metadata)),])
df_pca <- prcomp(combined)
df_out <- as.data.frame(df_pca$x)
p<-ggplot(df_out,aes(x=PC1,y=PC2 ))
p<-p+geom_point()
p



#CompareM =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/AAI_oceanospirilalles.tsv",header=TRUE)

#show graphically
ggplot(data=CompareM,aes(x= Mean_AAI,y=Orthologous_fraction))+geom_point()+geom_smooth(col="darkgreen",method="lm",se = TRUE)+theme_classic()
ggplot(data=CompareM,aes(x= Mean_AAI,y= orthologous_genes))+geom_point()+geom_point()+geom_smooth(col="darkgreen",method="lm",se = TRUE)+theme_classic()

# correlations
cor(CompareM$Mean_AAI, CompareM$orthologous_genes)
cor(CompareM$Mean_AAI, CompareM$Orthologous_fraction)

# statistics of linear regression
summary(lm(CompareM$orthologous_genes ~ CompareM$Mean_AAI))
summary(lm(CompareM$Orthologous_fraction ~ CompareM$Mean_AAI))

#transform the data in a matrix format
dat = CompareM[,c("Genome_A","Genome_B","Mean_AAI")]
g <- graph.data.frame(dat, directed=FALSE)
AAI = get.adjacency(g, attr="Mean_AAI", sparse=FALSE)

#genomes compared to themselves have 100 as a value
AAI[AAI<1] <- 100

#heatmap
heatmap.2(AAI,trace="none",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

#PCA
CompM.pca <- PCA(AAI)
plot(CompM.pca,cex=0.5)

#same, but for the Orhologous fraction
dat = CompareM[,c("Genome_A","Genome_B","Orthologous_fraction")]
g <- graph.data.frame(dat, directed=FALSE)
OF = get.adjacency(g, attr="Orthologous_fraction", sparse=FALSE)
OF[OF <1] <- 100
heatmap.2(OF,trace="none",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

CompM.pca <- PCA(OF)
plot(CompM.pca,cex=0.5)

#same, but for the Orhologous fraction
dat = CompareM[,c("Genome_A","Genome_B","orthologous_genes")]
g <- graph.data.frame(dat, directed=FALSE)
Og = get.adjacency(g, attr="orthologous_genes", sparse=FALSE)
Og[Og<1] <- NA
heatmap.2(Og,trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

CompM.pca <- PCA(Og)
plot(CompM.pca,cex=0.5)


######### - statistics - ##########

# Mantel test
#H0 - matrices are different, Ha - matrices correlate
mantel.rtest(as.dist(AAI), as.dist(OF), nrepet = 999)

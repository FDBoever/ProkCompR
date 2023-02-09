library(ggtree)
library(phyloseq)
library(ggjoy)
library(dplyr)

data("GlobalPatterns")
GP <- GlobalPatterns
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% 
  c("Feces", "Skin") 

mergedGP <- merge_samples(GP, "SampleType")
mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Order") 

melt_simple <- psmelt(mergedGP) %>% 
  filter(Abundance < 120) %>% 
  select(OTU, val=Abundance)

p <- ggtree(mergedGP) + 
  geom_tippoint(aes(color=Phylum), size=1.5)

facet_plot(p, panel="Abundance", data=melt_simple, 
           geom_joy, mapping = aes(x=val,group=label, 
                                        fill=Phylum), 
           color='grey80', lwd=.3)
           

#---------------------------           


 
           
tree <- rtree(30)

# Make the original plot
p <- ggtree(tree)


# generate some random values for each tip label in the data
d1 <- data.frame(id=tree$tip.label, val=rnorm(30, sd=3))

# Make a second plot with the original, naming the new plot "dot", 
# using the data you just created, with a point geom.
p2 <- facet_plot(p, panel="dot", data=d1, geom=geom_point, aes(x=val), color='red3')

# Make some more data with another random value.
d2 <- data.frame(id=tree$tip.label, value = abs(rnorm(30, mean=100, sd=50)))

# Now add to that second plot, this time using the new d2 data above, 
# This time showing a bar segment, size 3, colored blue.
p3 <- facet_plot(p2, panel='bar', data=d2, geom=geom_segment, 
           aes(x=0, xend=value, y=y, yend=y), size=3, color='blue4') 

# Show all three plots with a scale
p3 + theme_tree2()





concat_tree = read.tree("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile")
COPH = cophenetic(concat_tree)

p <- ggtree(concat_tree) + geom_tiplab(size=2)

p = ggtree(concat_tree) + geom_text2(size=2,hjust=0.5,vjust=0.2 aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))+	geom_tiplab(size=2)





head(bias_data)
bias_data$genome = gsub('.faa','',bias_data$genome)
rownames(bias_data) = bias_data$genome
dfSorted = bias_data[concat_tree $tip.label,]
colnames(dfSorted)=c('id','b')
# using the data you just created, with a point geom.
dfSorted$b=as.numeric(as.character(dfSorted$b))

p2 <- facet_plot(p, panel="dot", data= dfSorted, geom=geom_point, aes(x=b), color='red3')

p2



head(bias_data)

          
          
          


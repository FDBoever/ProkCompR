#-----------------------------------------------#
# Mantel tests
#-----------------------------------------------#

#NEED TO TIDY UP, CALCULATE ALL DISTANCE BASED MEASURES, AND RUN MANY DIFFERENT TESTS
# MANTEL ~ gene content , vs COPH, vs stuff like .....

dstANI = 1-as.matrix(ANIb)
dMan = distMan[SelectedTips,SelectedTips]
COPH = COPH[SelectedTips,SelectedTips]
distMan

mantal.1 = vegan::mantel(dstANI, dMan,method='pearson',permutations=999)
mantal.1 = vegan::mantel(COPH, dstANI,method='pearson',permutations=999)
mantal.1 = vegan::mantel(COPH, distMan,method='pearson',permutations=999)
mantal.1 = vegan::mantel(dstANI, distMan[SelectedTips,SelectedTips],method='pearson',permutations=999)
mantal.1 = vegan::mantel(COPH, dstANI,method='pearson',permutations=999)
mantel.2 = vegan::mantel(ANIb[rownames(geo.dist),rownames(geo.dist)], geo.dist,method='pearson',permutations=999 )
mantel.2 = vegan::mantel(geo.dist,ANIb[rownames(geo.dist),rownames(geo.dist)],method='pearson',permutations=999 )


mantel.correlog(dstANI, dMan)
mantel.correlogram <- vegan::mantel.correlog(geo.dist,ANIb[rownames(geo.dist),rownames(geo.dist)])
plot(mantel.correlogram)

#=========================#

vegan::varpart(geo.dist,ANIb[rownames(geo.dist),rownames(geo.dist)], ~ SS, COPH[rownames(geo.dist),rownames(geo.dist)],data=metadata[rownames(geo.dist),])

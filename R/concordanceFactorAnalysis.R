#http://www.robertlanfear.com/blog/files/concordance_factors.html

library(viridis)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(GGally)
library(entropy)

# read the data
#d = read.delim("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concord.cf.stat", header = T, comment.char='#')
#d = read.delim("~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair/concord.cf.stat", header = T, comment.char='#')


#fort he  96 tree
d = read.delim("~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair_trimAl/concord.cf.stat", header = T, comment.char='#')
t <- read.tree("~/DATA/MarbGenomics/MB_SCO_Protein_alignments_mafft_genafpair_trimAl/concat.treefile")

#for the 106 tree
d = read.delim("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concord.cf.stat", header = T, comment.char='#')
t <- read.tree("~/DATA/MarbGenomics/iqtree_out/SCO_filtered_mafft_genafpair_trimAl/concat.treefile")

ggtree(t)+geom_tiplab()

names(d)[18] = "bootstrap"
names(d)[19] = "branchlength"


# find a branch of interest
d[d$ID==115,]

#gCF= 26.19
#gN (number of decisive trees)
#gCF_N = 33 (the number of trees that showed it)

#gCF_N = 33 (33 trees supported this branching)
#gDF1_N = 3 (3 trees resolve the 2nd resolution of the clades around that branch)
#gDF2_N = 0 (0 trees resolve the 3rth)
#gDFP_N = 90 --> the other trees had a topoplogy that wasnt any one of the three possible arrangements of the four clades around this branch, which is a typical signal of noisy single locus trees

#--> Most common resulution of this branch in the gene trees is NOY the one we see in the ML concatenated trees
# at least two reasons why that might be the case -  we might be in the "anomaly zone", or it might be that the informatice sietes for this tree are scattered around in all the noisy trees, such that individial gene trees really arent giging us much useful information here
#The sCF and sDF values suggest that the latter explanation is true.
#highlights a limitation of ethods that look to build a species tree from pre-estimated gene trees, when the trees are noisy these methods will struggle
#(which indicates high variance in the support that sites give for the correct resolution of this branch).

#branch is short? probably most of the reason for the low gCF (short branches are hard to resolve particularly with short loci)
#the way boostrap is calculated, we indeed expect a high support if we have 131 vs 113 sites competing to resolve a single branch
#regardless, these sCF values tell us that despite the bootstrap support of 100%, the underlying data contain a lot of discordance (almost the maxmimim possible amound ) around hthis one branch





#------------------

#Using concordance factors to test the assumptions of an ILS model
#if the discrodance among gene tress and sites come from incomplete lineage sorting, we can maek a simple and testable prediction:
#that the number of gene trees or sites supporting the two discordant topologies should be roughly equal.
#Both of these ideas have been around for some time (for genes [link to Huson et al & Steel 2005 Recomb] and for sites [link to Greene et al 2010]),


# The basic idea is that we count up the genes or sites supporting the two discordant topologies, and use a chi-square test to see if they're significantly different.
# this requires a number of assumptions to hold, and if you're willing to make those assumptions, here is a simple way to calcularte the probability that the data can reject wqual frequencies for genes (gEFp, below) and sites (sEFp)

# test ILS assumptions

# first we use a slightly modified chisq function
# which behaves nicely when you feed it zeros
chisq = function(DF1, DF2, N){
    tryCatch({
        # converts percentages to counts, runs chisq, gets pvalue
        chisq.test(c(round(DF1*N)/100, round(DF2*N)/100))$p.value
    },
    error = function(err) {
        # errors come if you give chisq two zeros
        # but here we're sure that there's no difference
        return(1.0)
    })
}

e = d %>%
    group_by(ID) %>%
    mutate(gEF_p = chisq(gDF1, gDF2, gN)) %>%
    mutate(sEF_p = chisq(sDF1, sDF2, sN))


subset(data.frame(e), (gEF_p < 0.05 | sEF_p < 0.05))

#which shows that a lot of branches reject the ILS assumptions, assuming our chi-squared approach is accurate (which it won’t be among sites in a single gene because of linkage disequilibrium, so be careful with these p values):



# calculate internode certainty
IC = function(CF, DF1, DF2, N){
    # convert to counts
    X = CF * N / 100
    Y = max(DF1, DF2) * N / 100
    pX = X/(X+Y)
    pY = Y/(X+Y)
    IC = 1 + pX * log2(pX) +
             pY * log2(pY)
    return(IC)
}

e = e %>%
    group_by(ID) %>%
    mutate(gIC = IC(gCF, gDF1, gDF2, gN)) %>%
    mutate(sIC = IC(sCF, sDF1, sDF2, sN))

# entropy
ENT = function(CF, DF1, DF2, N){
    CF = CF * N / 100
    DF1 = DF1 * N / 100
    DF2 = DF2 * N / 100
    return(entropy::entropy(c(CF, DF1, DF2)))
}

ENTC = function(CF, DF1, DF2, N){
    maxent = 1.098612
    CF = CF * N / 100
    DF1 = DF1 * N / 100
    DF2 = DF2 * N / 100
    ent = entropy::entropy(c(CF, DF1, DF2))
    entc = 1 - (ent / maxent)
    return(entc)
}

e = e %>%
    group_by(ID) %>%
    mutate(sENT = ENT(sCF, sDF1, sDF2, sN)) %>%
    mutate(sENTC = ENTC(sCF, sDF1, sDF2, sN))


#=========================

#




# plot it
GGally::ggpairs(e, columns = c(2, 12, 18, 19,20,21,22,23,24))


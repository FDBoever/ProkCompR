#Installing and loading packages
#-------------------------------------------------------------#

#devtools
if(!require(devtools)) install.packages("devtools", repos = "http://cran.us.r-project.org",dependencies=TRUE)

#-------------------------------------------------------------#
#   Data manipulation
if(!require(plyr)) install.packages("plyr", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(tidyr)) install.packages("tidyr", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(tibble)) install.packages("tibble", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(gridExtra)) install.packages("gridExtra", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(reshape)) install.packages("reshape", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(reshape2)) install.packages("reshape2", repos = "http://cran.us.r-project.org",dependencies=TRUE)

#-------------------------------------------------------------#
#   Phylogenetics
#if(!require(future.apply)) devtools::install_github("HenrikBengtsson/future.apply")
if(!require(ape)) install.packages("ape", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(phyloseq)) install.packages("phyloseq", repos = "http://cran.us.r-project.org",dependencies=TRUE)

if(!require(phyloseq)) BiocManager::install("phyloseq")
if(!require(phytools)) install.packages("phytools", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(phangorn)) install.packages("phangorn", repos = "http://cran.us.r-project.org",dependencies=TRUE)

# -- since I am using R version 3.6.3 I had to install from source, and download an older version from CRAN
# first wget https://cran.r-project.org/src/contrib/Archive/phylolm/phylolm_2.5.tar.gz
# install.packages('phylolm_2.5.tar.gz', repos = NULL, type="source", dependencies=T)
if(!require(phylolm)) install.packages("phylolm", repos = "http://cran.us.r-project.org",dependencies=TRUE)

#-------------------------------------------------------------#
#   Analysis
if(!require(Rcpp)) devtools::install_github("RcppCore/Rcpp")
if(!require(Rtsne)) devtools::install_github("jkrijthe/Rtsne")

if(!require(Rtsne)) install.packages("Rtsne", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(micropan)) install.packages("micropan", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(randomForest)) install.packages("randomForest", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(tidymodels)) install.packages("tidymodels", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ranger)) install.packages("ranger", repos = "http://cran.us.r-project.org",dependencies=TRUE)

#-------------------------------------------------------------#
#For visualisation purposes
if(!require(gplots)) install.packages("gplots", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(RColorBrewer)) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ggrepel)) install.packages("ggrepel", repos = "http://cran.us.r-project.org",dependencies=TRUE)
#if(!require(ggtree)) install.packages("ggtree", repos = "http://cran.us.r-project.org",dependencies=TRUE)
#if(!require(ggtree)) devtools::install_github("YuLab-SMU/ggtree")
if(!require(ggtree)) BiocManager::install("ggtree")
if(!require(ggraph)) install.packages("ggraph", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(gggenes)) install.packages("gggenes", repos = "http://cran.us.r-project.org",dependencies=TRUE)
#if(!require(gggenes)) devtools::install_github("wilkox/gggenes")

#database related
suppressMessages(if(!require(RSQLite)) devtools::install_github("rstats-db/RSQLite"))
suppressMessages(if(!require(RMySQL)) install.packages("RMySQL", repos = "http://cran.us.r-project.org",dependencies=T))
suppressMessages(if(!require(pool)) install.packages("pool", repos = "http://cran.us.r-project.org",dependencies=T))

#-------------------------------------------------------------#
#other
#if(!require(ggmsa)) devtools::install_github('YuLab-SMU/ggmsa')
#if(!require(seqmagick)) devtools::install_github('YuLab-SMU/seqmagick')
#library("rfUtilities") # to test model significance
#library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function

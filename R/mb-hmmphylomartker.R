#hmm2fasta

#uses functions from mb-hmm.R

#load Hmm SCAN output
campbell.tbl = hmm2df(inPath="~/DATA/MarbGenomics/dbCAN_Campbell", method="tophit",evalue.cutoff=1.e-10, clean_genome_names=TRUE)
rinke.tbl = hmm2df(inPath="~/DATA/MarbGenomics/dbCAN_Rinke", method="tophit",evalue.cutoff=1.e-10, clean_genome_names=TRUE)


# convert to pan genome to inspect if these are SCO
campbell.pan <- tbl2pan(campbell.tbl)
rinke.pan <- tbl2pan(rinke.tbl)

#extract the SCO's only
campbell.sco <- pan2sco(campbell.pan)
rinke.sco <- pan2sco(rinke.pan)

#---------!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#--------!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Solve the issue of transposed pan tables in
dim(campbell.pan)
dim(campbell.sco)

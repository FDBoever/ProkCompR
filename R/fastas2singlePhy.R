make all.paml


#----

libary(chopper)

#copy the fasta folder to an phy outdir, so that when rerun, you can succesfully recompute

dir.create(out_dir)


fastaFiles <- dir('~/DATA/MarbGenomics/filtered_alignments/RibosomalProteins_mafft_enafpair_pal2nal_phy', pattern='*.fa', full.names = T)



for(f in fastaFiles){
	fas2phy(f, format = "sequential", overwrite = FALSE)
}


#in dir, cat *.phy > all.paml
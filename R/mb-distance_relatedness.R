# MANTEL TESTS
#---------------------------------------------------#


# Custom functions
#---------------------------------------------------#
#' mantle.on.list to calculate mantle test statistics from a list of distances
#' computes mantle test statistics for all combinations of distance matrices in list object
#' outputs a table
#' @param sim.list
#' @param method
#' @param perm
#'
#' @return
#' @export
#'
#' @examples
mantle.on.list <- function(sim.list, method='pearson',perm=999){
  out.df <- NULL
  for(i in names(sim.list)[1:length(names(sim.list))]){
    message('analysing ',i)
    for(j in names(sim.list)[1:length(names(sim.list))]){
      m.mantel = vegan::mantel(sim.list[[i]], sim.list[[j]],method='pearson',permutations=999)
      c.res  <-  c('i'=i,
                   'j'=j,
                   'method' = as.character(m.mantel[2]),
                   'statistic' = as.numeric(m.mantel[3]),
                   'signif' = as.numeric(m.mantel[4]))
      out.df <- rbind(out.df, c.res)
    }
  }

  #tidy up the output data frame
  out.df <- out.df %>%
    data.frame() %>%
    dplyr::mutate(statistic = as.numeric(as.character(statistic))) %>%
    dplyr::mutate(signif = as.numeric(as.character(signif))) %>%
    dplyr::mutate(i= factor(i,levels=names(sim.list))) %>%
    dplyr::mutate(j= factor(j,levels=rev(names(sim.list)))) %>%
    dplyr::mutate(sig = ifelse(signif<0.01,TRUE,FALSE)) %>%
    dplyr::mutate(comparison = paste0(i,'_',j))

  return(out.df)
}

#=================================================================#


#Set some variables
outdir = "~/DATA/MarbGenomics/Graphs/"

#create list object containing distance matrices
sim.list <- list(
  'COPH' = COPH[sel.genome,sel.genome],
  'ANIb' = 1-ANIb[sel.genome,sel.genome],
  'AAI' = 1-AAI[sel.genome,sel.genome],
  'OF' =  1-Ortho.Fraction[sel.genome,sel.genome],
  'GGR' = 1-pan.ggr[sel.genome,sel.genome],
  'OG'= OG.dist[sel.genome,sel.genome],
  'COG'= COG.dist[sel.genome,sel.genome],
  'KO'= kfm.dist[sel.genome,sel.genome],
  'CAZY'= cazy.dist[sel.genome,sel.genome],
  'TCDB'= tcdb.dist[sel.genome,sel.genome]
)


#calculate mantle test statistics from a list of distances
out.df <- mantle.on.list(sim.list,method='pearson',perm=999)


# Extract the lower triangle
out.matrix <- out.df %>%
  tidyr::pivot_wider(i, names_from = j, values_from = statistic)

lower.triangle <- out.matrix %>%
  dplyr::select(-i) %>%
  as.matrix()

rownames(lower.triangle) = colnames(lower.triangle)
lower.triangle[lower.tri(lower.triangle)] <- NA


# Generate table for visualisation and output
mantel.plot.df <- lower.triangle %>% data.frame() %>%
  tibble::rownames_to_column(var='o') %>%
  tidyr::pivot_longer(cols=COPH:TCDB,names_to='m') %>%
  dplyr::mutate(comparison = paste0(o,'_',m)) %>%
  dplyr::left_join(out.df,by='comparison')


# Save down supplementary data table (Mantel tests)
write.table(x=mantel.plot.df %>% data.frame(),file=paste0(outdir,'/SI_table_mantel','.txt'),quote=FALSE,row.names = FALSE, sep='\t')


#plot supplementary figure mantel tests
p.mantel_plot <- mantel.plot.df %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(i= factor(i,levels=names(sim.list))) %>%
  dplyr::mutate(j= factor(j,levels=rev(names(sim.list)))) %>%
  ggplot2::ggplot(aes(i,j))+
  ggplot2::geom_point(ggplot2::aes(fill=value),shape=21,size=6)+
  ggplot2::geom_text(aes(label=round(value,1)),color='white',size=2.8)+
  ggplot2::scale_fill_continuous(high = "navy", low = "skyblue")+
  fdb_style() +
  ggplot2::theme(axis.text.x = element_text(angle = 90,hjust=1))

ggplot2::ggsave(filename=paste0(outdir,'mantel_plot_all_v_all_numbered_navy','.pdf'),plot=p.mantel_plot, width = 4, height = 3,unit='in')

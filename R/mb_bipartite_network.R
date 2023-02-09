# bi partite network from pan genome tables

#inc <- matrix(sample(0:1, 15, repl=TRUE), 3, 5)
#colnames(inc) <- letters[1:5]
#rownames(inc) <- LETTERS[1:3]

ph2gen <- metadata$phylogroup
names(ph2gen) <- metadata$genome

q <- graph_from_incidence_matrix(OG.pan[sel.genome, allSig_hgm[1:100]],weighted = TRUE)
q <- graph_from_incidence_matrix(OG.pan[1:30, 1:100],weighted = TRUE)
q <- graph_from_incidence_matrix(OG.pan[sel.genome, allSig_hgm[1:500]],weighted = TRUE)


#q <- graph_from_incidence_matrix(OG.pan[sel.genome, sig_enriched_raw[1:500]],weighted = TRUE)
q <- graph_from_incidence_matrix(m.pan[],weighted = TRUE)

q <- igraph::simplify(q, remove.multiple=TRUE, remove.loops=TRUE)
q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < 2)])
q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight > 50)])

#q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < median(E(q)$weight))])

col <- c("steelblue", "orange", "green")
shape <- c(21,22)


qgn = ggnetwork::ggnetwork(intergraph::asNetwork(q),
                           layout="fruchtermanreingold",
                           weights="weight",
                           niter = 1000)
qgn <- qgn %>%
  mutate(v_type = ifelse(type==FALSE,'genome','protein')) %>%
  mutate(group = ph2gen[vertex.names])

ggplot2::ggplot(data = qgn,
                aes(x, y, xend = xend, yend = yend))+
  ggnetwork::geom_edges(color = "black",aes(alpha=weight)) + #,alpha=weight
  ggnetwork::geom_nodes(aes(shape=v_type,fill=group,size=v_type),color='black')+
  scale_fill_manual(values=phylogroup.colors)+
  scale_shape_manual(values=c(21,22))+
  scale_size_manual(values=c(2,0.3))


#q <- igraph::graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "upper", weighted = T, diag = F)
#q <- igraph::simplify(q, remove.multiple=TRUE, remove.loops=TRUE)
#q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < weight_threshold)])

#singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])==1])
#non_singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])>1])
#V(q)$comm = ifelse(cliques[names(V(q))] %in% singleton_cl, NA, cliques[names(V(q))] )

#n_vertices <- length(V(q))

#qgn = ggnetwork::ggnetwork(intergraph::asNetwork(q),
#                           layout="fruchtermanreingold",
#                           weights="weight",
#                           niter = 1000)

p.out <- ggplot2::ggplot(data = qgn,
                         aes(x, y, xend = xend, yend = yend)) +
  ggnetwork::geom_edges(color = "black") + #,alpha=weight
  ggnetwork::geom_nodes( size = 1.5, shape=21, alpha=0.7,aes(fill=as.factor(comm)))+
  ggtitle(paste0(weight_threshold,' (n=',n_vertices,')'))+
  theme_void()+theme(aspect.ratio=1)+
  scale_size(range=c(0.1,2))








plot.igraph <- function(similarity_matrix, weight_threshold, remove_singletons = TRUE,seed=04012021,colors=NULL,community.method='cluster_infomap'){
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  set.seed(seed)
  n_genomes <- nrow(similarity_matrix)
  q <- igraph::graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "upper", weighted = T, diag = F)
  q <- igraph::simplify(q, remove.multiple=TRUE, remove.loops=TRUE)
  q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < weight_threshold)])

  #Remove single nodes
  if(remove_singletons == TRUE){
    q = delete.vertices(q, which(igraph::degree(q)< 2))
  }

  if(!(community.method %in% c('cluster_walktrap','cluster_edge_betweenness','cluster_fast_greedy','cluster_infomap','cluster_label_prop','cluster_spinglass'))){
    stop('make sure you take the right method')
  }

  if(community.method == 'cluster_walktrap'){
    wc <- igraph::cluster_walktrap(q)
  }
  if(community.method == 'cluster_edge_betweenness'){
    wc <-igraph::cluster_edge_betweenness(q, weights = E(q)$weight, directed = FALSE)
  }
  if(community.method == 'cluster_fast_greedy'){
    wc <-igraph::cluster_fast_greedy(q, weights = E(q)$weight)
  }
  if(community.method == 'cluster_infomap'){
    wc <-igraph::cluster_infomap(q)
  }
  if(community.method == 'cluster_label_prop'){
    wc <-igraph::cluster_label_prop(q)
  }
  if(community.method == 'cluster_spinglass'){
    wc <-igraph::cluster_spinglass(q)
  }

  #wc <- cluster_infomap(q)
  cliques = membership(wc)

  V(q)$cliques = cliques[names(V(q))]

  #groups <- genome.tbl %>% filter(genome %in% names(V(q))) %>% select(SS) %>% pull() %>% as.character()
  #names(groups) <- names(V(q))
  #V(q)$groups <- groups

  #lat <- genome.tbl %>% filter(genome %in% names(V(q))) %>% select(lat) %>% pull() %>% abs() #%>% as.character()
  #names(lat) <- names(V(q))
  #V(q)$lat <- lat

  singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])==1])
  non_singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])>1])
  V(q)$comm = ifelse(cliques[names(V(q))] %in% singleton_cl, NA, cliques[names(V(q))] )

  n_vertices <- length(V(q))
  qgn = ggnetwork::ggnetwork(intergraph::asNetwork(q),
                             layout="fruchtermanreingold",
                             weights="weight",
                             niter = 1000)

  p.out <- ggplot2::ggplot(data = qgn,
                           aes(x, y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(color = "black") + #,alpha=weight
    ggnetwork::geom_nodes( size = 1.5, shape=21, alpha=0.7,aes(fill=as.factor(comm)))+
    ggtitle(paste0(weight_threshold,' (n=',n_vertices,')'))+
    theme_void()+theme(aspect.ratio=1)+
    scale_size(range=c(0.1,2))

  if(is.null(colors)){
    p.out <- p.out + scale_fill_manual(values=getPalette(length(non_singleton_cl)))
  }else{
    p.out <- p.out + scale_fill_manual(values=colors)
  }

  p.out
}



# SCREEN ANI VALUES AND BIODIVERSITY LOSS
#=================================================================================#
# CLIQUE detrction, ~ clique is maximal sibset of a network, where each node is a set of adjacent to all other nodes in the set
#cl <- cliques(q, min = 2, max = 106)
#largest_cliques(q)
#eb <- cluster_edge_betweenness(q)

#cliq <- largest_cliques(q)              # get the largest cliques of g
#plot(induced_subgraph(q, cliq[[1]]))

#dg <- decompose(q, mode="weak",
#                min.vertices = 5)
#plot.igraph(dg[[1]])   # Plot the largest one



#=================================================================================#
# FUNCTIONS
#=================================================================================#

#clique.igraph
#=================================================================================#
#' calculate igrapgh object and remove edges below threshold
#'
#' @param similarity_matrix
#' @param weight_threshold
#' @param remove_singletons
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
clique.igraph <- function(similarity_matrix, weight_threshold, remove_singletons = TRUE,seed=04012021){
  set.seed(seed)
  n_genomes <- nrow(similarity_matrix)
  q <- igraph::graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "upper", weighted = T, diag = F)
  q <- igraph::simplify(q, remove.multiple=TRUE, remove.loops=TRUE)
  q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < weight_threshold)])
  if(remove_singletons == TRUE){
    q = igraph::delete.vertices(q, which(igraph::degree(q)< 2))
  }
  wc <- igraph::cluster_walktrap(q)
  cliques = igraph::membership(wc)
  V(q)$cliques = cliques[names(V(q))]
  n_vertices <- length(V(q))
  return(q)
}

#range.igraph
#=================================================================================#
#' runs over clique.igraph at a range of thresholds (cutogffs).
#'
#' @param similarity_matrix
#' @param cutoff_range
#'
#' @return
#' @export
#'
#' @examples
range.igraph <- function(similarity_matrix, cutoff_range){
  out.df <- NULL
  for(ctf in cutoff_range){
    message(paste0('analysing cliques at >',ctf))
    q <- clique.igraph(similarity_matrix = similarity_matrix, weight_threshold = ctf, remove_singletons = TRUE)
    n_vertices <- length(V(q))
    message(n_vertices)
    out.df <- rbind(out.df,c('ctf'=ctf,'n'=n_vertices))
  }
  return(out.df)
}



#detect_cliques
#=================================================================================#
#' function that detects cliques in distance matrix at specified cutoff
#' define graph, delete edges < cutoff, asses cliques with igraph::cluster_walktrap and identifies genomes that belong to clusters
#' @param distance
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
#' cliques = detect_cliques(distance=ANIb, cutoff = 0.90)

detect_cliques <- function(distance, cutoff, community.method = 'cluster_walktrap'){
  if(!(community.method %in% c('cluster_walktrap','cluster_edge_betweenness','cluster_fast_greedy','cluster_infomap','cluster_label_prop','cluster_spinglass'))){
    stop('make sure you take the right method')
  }
  g <- igraph::graph_from_adjacency_matrix(as.matrix(distance), mode = "upper", weighted = T, diag = F)
  g2 <- igraph::delete.edges(g, which(E(g)$weight <cutoff))

  if(community.method == 'cluster_walktrap'){
    wc <- igraph::cluster_walktrap(g2)
  }
  if(community.method == 'cluster_edge_betweenness'){
    wc <-igraph::cluster_edge_betweenness(g2, weights = E(g2)$weight, directed = FALSE)
  }
  if(community.method == 'cluster_fast_greedy'){
    wc <-igraph::cluster_fast_greedy(g2, weights = E(g2)$weight)
  }
  if(community.method == 'cluster_infomap'){
    wc <-igraph::cluster_infomap(g2, e.weights = E(g2)$weight)
  }
  if(community.method == 'cluster_label_prop'){
    wc <-igraph::cluster_label_prop(g2)
  }
  if(community.method == 'cluster_spinglass'){
    wc <-igraph::cluster_spinglass(g2)
  }

  cliques = igraph::membership(wc)
  return(cliques)
}

#multiClique
#=================================================================================#
#' function that detects cliques in distance matrix at multiple specified cutoffs
#'
#' @param distance
#' @param cutoffs
#' @param no_singletons
#'
#' @return
#' @export
#'
#' @examples
#' multi_clique.ANIb <- multiClique(distance=ANIb, cutoffs = c(0.7,0.75,0.80,0.85,0.90,0.95,0.98),no_singletons = TRUE)

multiClique <- function(distance, cutoffs, no_singletons=TRUE, community.method = 'cluster_walktrap'){
  output <- NULL
  output <- data.frame(cbind(output,'genome'=rownames(distance)))
  for(ctf in cutoffs){
    cliques <- detect_cliques(distance=distance, cutoff = ctf, community.method = community.method)
    if(no_singletons==TRUE){
      cliques <- cliques[cliques %in% names(table(cliques)[table(cliques)>1])]
    }
    output <- cbind(output, unname(cliques[as.character(output$genome)]))
  }

  colnames(output) = c('genome',paste0('c',cutoffs))
  return(output)
}

#sortCliques
#=================================================================================#
#' Function that will use output of multiClique and classifies them hierarchically
#' it will change the numbers names of the clusters so that they can be viewed in a hierarchical manner
#'
#' @param multi.clique
#'
#' @return
#' @export
#'
#' @examples
#' multi.clique <- multi_clique.ANIb[,c('genome','c0.75','c0.8', 'c0.85', 'c0.9', 'c0.95', 'c0.98')]

sortCliques <- function(multi.clique, tree){
  message('sorting cliques hierarchically')
  if('genome' %in% colnames(multi.clique)){
    col.nms <- colnames(multi.clique)[!colnames(multi.clique)=='genome']
  }else{
    col.nms <- colnames(multi.clique)
  }

  rownames(multi.clique) <- as.character(multi.clique$genome)
  multi.clique <- multi.clique[tree$tip.label,]
  #---------- SORT ALONG rowname order (conveniently here, the order in tree)
  sorted <- NULL
  for(cln in col.nms){
    ccv <- as.factor(as.character(multi.clique[,cln]))
    detection.order <- unique(as.character(multi.clique[,cln]))
    detection.order <- detection.order[!is.na(detection.order)]
    convert.order = 1:length(detection.order)
    names(convert.order) = detection.order

    #levels(ccv)=c(detection.order)
    sorted <- cbind(sorted, as.numeric(unname(convert.order[as.character(ccv)])))
  }
  colnames(sorted)<-col.nms
  rownames(sorted)<-multi.clique$genome
  sorted<-sorted %>% data.frame() %>% tibble::rownames_to_column(var='genome')

  #find hierarchy
  reconstruct <-NULL
  for(i in 2:(length(col.nms))){
    deepgroup <- as.factor(as.character(sorted[,col.nms[i-1]]))
    names(deepgroup)<-sorted$genome
    dpgrp.uniq <-   unique(as.numeric(deepgroup))
    dpgrp.uniq <- dpgrp.uniq[!is.na(dpgrp.uniq)]

    sb <- NULL
    for(clique in dpgrp.uniq){
      sub.f <-  sorted[sorted[,col.nms[i-1]] == clique,]
      sub.f <- sub.f[!is.na(sub.f$genome),]
      detection.order <- unique(as.character(sub.f[,col.nms[i]]))
      detection.order <- detection.order[!is.na(detection.order)]
      if(length(detection.order)==0){
        ccv <- rep(NA,nrow(sub.f))
        names(ccv)<- sub.f$genome
      }
      if(length(detection.order)==1){
        ccv <- rep(1,nrow(sub.f))
        names(ccv)<- sub.f$genome
      }
      if(length(detection.order)>1){
        names(detection.order) <- 1:length(detection.order)
        detection.order = setNames(names( detection.order),  detection.order)
        ccv <- detection.order[as.character(sub.f[,col.nms[i]])]
        names(ccv)<- sub.f$genome
      }
      sb <- c(ccv, sb)
    }
    sb <- as.numeric(as.character(unname(sb[sorted$genome])))
    reconstruct <- cbind(reconstruct,unname(sb))
  }

  reconstruct <- cbind('c'=sorted[,col.nms[1]], reconstruct)
  colnames(reconstruct) = paste0('c',col.nms)
  reconstruct <- reconstruct %>% data.frame() %>%
    mutate_if(is.factor, ~as.numeric(as.character(.)))
  rownames(reconstruct) = sorted$genome

  #compile hrc (hierarchical clique)
  hrc <- NULL
  for(i in 1:nrow(reconstruct)){
    rw <- reconstruct[i,paste0('c',col.nms)]
    rw <- rw[!is.na(rw)]
    rw <- paste(rw,collapse='.')
    hrc <- c(hrc,rw)
  }

  reconstruct <- cbind(reconstruct, 'hrc'=hrc)
  reconstruct <- reconstruct %>% data.frame() %>% rownames_to_column(var='genome')
  colnames(sorted) <- c('genome',paste0('s',col.nms))
  combined.out <- multi.clique %>% left_join(sorted, by='genome') %>% left_join(reconstruct, by='genome')
  return(combined.out)
}



#plot.cliqueSort
#=================================================================================#
#' Title
#'
#' @param tree
#' @param sorted.cliques
#' @param circular
#'
#' @return
#' @export
#'
#' @examples
plot.cliqueSort <- function(tree, sorted.cliques,circular=TRUE,drop.outgroup=TRUE, anchor ='sc0.8',out.colors=NULL){

  if(drop.outgroup==TRUE){
    c.offset=0.07
    c.size=0.2 #if you want separating lines
    c.width=0.05
  }else{
    c.offset=0.01
    c.size=0.4
    c.width=0.2
  }

  #multi_clique_prep <- sorted.cl %>% mutate(ID=genome) %>% select(ID, everything()) %>% select(-genome)
  multi_clique_prep <- sorted.cliques %>%
    dplyr::mutate(ID=genome) %>%
    dplyr::select(ID, everything()) %>%
    dplyr::select(-genome)
  circ <- ggtree::ggtree(tr.tree, layout="fan", open.angle=10) +
    ggtree::geom_tiplab(align=TRUE,size=0,linetype='dashed',linesize=0.3)

  p1 <- circ %<+% multi_clique_prep
  p2 <- p1

  # NO COLORS PREDERINED
  if(is.null(out.colors)){
    sccols <- multi_clique_prep %>%
      dplyr::select(starts_with('sc')) %>%
      colnames()

    for(sccol in sccols){
      if(length(multi_clique_prep[,sccol]) != length(multi_clique_prep[,sccol][is.na(multi_clique_prep[,sccol])])){
        p2 <- p2 +
          ggnewscale::new_scale_fill() +
          ggtreeExtra::geom_fruit(geom=geom_tile,
                                  mapping=aes_string(fill=paste0("as.factor(",sccol,")"),alpha=paste0("-as.numeric(as.character(",sccol,"))")),
                                  color = "white", offset = c.offset,size = c.size,width=c.width)+
          ggplot2::scale_fill_discrete(na.value = 'white')+
          ggplot2::scale_alpha(na.value=0,range=c(0.5,1))
      }
    }

    #WITH PREDEFINED HIERARCHICAL COLORS from out.colors
  }else{
    sccols <- multi_clique_prep %>%
      dplyr::select(starts_with('sc')) %>%
      colnames()
    #sccols = sccols[10:20]
    for(sccol in sccols){

      #Check if all are NA
      if(length(multi_clique_prep[,sccol]) != length(multi_clique_prep[,sccol][is.na(multi_clique_prep[,sccol])])){
        #cls.scheme <- cliquelevel2colors(sorted.cliques = sorted.cl,
        #                   out.colors = out.colors,
        #                   sel.col=gsub('sc','cc',sccol))

        cls.df <- multi_clique_prep %>%
          dplyr::left_join(out.colors,by=c('ID'='genome')) %>%
          dplyr::select(sccol,paste0('col.',gsub('sc','cc',sccol))) %>%
          dplyr::filter(!is.na(!!as.name(sccol))) %>%
          unique()

        cls.scheme <- as.character(cls.df[,2])
        names(cls.scheme) = cls.df[,1]
        if(length(cls.scheme)==0){
          cls.scheme = 'white'
          message('attention, i manually assign')
        }else{
          cls.scheme <- c(cls.scheme,NA)
        }
        message(cls.scheme)
        p2 <- p2 + ggnewscale::new_scale_fill() +
          ggtreeExtra::geom_fruit(geom=geom_tile,
                                  mapping=aes_string(fill=paste0("as.factor(",sccol,")")),
                                  color = "white", offset = c.offset,size = c.size,width=c.width)+
          ggplot2::scale_fill_manual(values=unname(cls.scheme), na.value='white')
      }
    }
  }



  if(circular==FALSE){
    p2 <- p2 + layout_rectangular()

  }
  return(p2)
}



#blend_cluster_colors
#=================================================================================#
### IMPLEMENT EXTRACTION OF COLORS FROM THE BLOODY
#' Title
#'
#' @param sorted.cliques
#' @param columns
#' @param base.column
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
#' out.colors <- blend_cluster_colors(sorted.cliques, columns, base.column, colors )

blend_cluster_colors <- function(sorted.cliques, columns, base.column, colors ){
  base.colors <- NULL
  trans.colors <- NULL
  message("assigning colors, this shouldn't take too long ")
  for(row in 1:nrow(sorted.cliques)){
    #sorted.cliques[row,columns]
    base.color <- colors[sorted.cliques[row,columns[1]]]
    base.colors = c(base.colors, base.color)

    color.tans.df = NULL
    i<-0
    t.cols = NULL
    for(col in columns[1:length(columns)]){
      i <- i + 1
      if(is.na(sorted.cliques[row,col])){
        #if(is.na(sorted.cliques[row,columns2[i]])){
        t.col <- NA
        #message('assigned to NA as no group')
      }else{
        if(sorted.cliques[row,col] == 1){
          if(i==1){
            t.col <- base.color
            #message('assigned to base as i =1')
          }else{
            if(!is.null(t.cols)){
              #previsously assigned
              t.col <- t.cols[length(t.cols)]
              #message('assigned to previosu')

            }else{message('somthing went wrong')}
          }
        }else{
          if(!is.null(t.cols)){
            rgb.val <- col2rgb(t.cols[length(t.cols)])
            r <- rgb.val[1]-(sorted.cliques[row,col]*7)-(length(t.cols)*2)
            #r <- rgb.val[1]+(sorted.cliques[row,columns2[i]]*10)
            if(r>255){r=255}
            if(r<0){r=0}

            g <- rgb.val[2]-(sorted.cliques[row,col]*10)-(length(t.cols)*2)
            #g <- rgb.val[2]+(sorted.cliques[row,columns2[i]]*10)
            if(g>255){g=255}
            if(g<0){g=0}

            b <- rgb.val[3]+(sorted.cliques[row,col]*10)+(length(t.cols)*2)
            #b <- rgb.val[3]+(sorted.cliques[row,columns2[i]]*10)
            if(b>255){b=255}
            if(b<0){b=0}
          }else{
            rgb.val <- col2rgb(base.color)
            r <- rgb.val[1]-(sorted.cliques[row,col]*7)
            #r <- rgb.val[1]+(sorted.cliques[row,columns2[i]]*10)
            if(r>255){r=255}
            if(r<0){r=0}

            g <- rgb.val[2]-(sorted.cliques[row,col]*10)
            #g <- rgb.val[2]+(sorted.cliques[row,columns2[i]]*10)
            if(g>255){g=255}
            if(g<0){g=0}

            b <- rgb.val[3]+(sorted.cliques[row,col]*10)
            #b <- rgb.val[3]+(sorted.cliques[row,columns2[i]]*10)
            if(b>255){b=255}
            if(b<0){b=0}
          }



          t.col <- rgb(r, g, b, max = 255)
          #message('assigned to kinky color')
        }
      }
      t.cols = c(t.cols,t.col)
    }
    #message(t.cols)
    trans.colors <- rbind(trans.colors, t.cols)
  }
  out.colors <- cbind(trans.colors)
  colnames(out.colors) <- paste0('col.',columns)
  rownames(out.colors) <- sorted.cliques$genome

  out.colors <- out.colors %>%
    data.frame() %>%
    tibble::rownames_to_column(var='genome')

  message('done...')
  return(out.colors)
}


#cliquelevel2colors
#=================================================================================#
#' annotates colors respective to cutoff levels
#'
#' @param sorted.cliques
#' @param out.colors
#' @param sel.col
#'
#' @return
#' @export
#'
#' @examples
#' cliquelevel2colors(sorted.ANI.cliques,out.colors,'sc0.8')

cliquelevel2colors <- function(sorted.cliques,out.colors, sel.col){
  mtd <- genome.tbl %>%
    dplyr::left_join(sorted.cliques, by='genome') %>%
    dplyr::left_join(out.colors, by='genome') %>%
    data.frame()

  mtd$grp <- mtd[,gsub('cc','sc',sel.col)]

  mtd <- mtd %>%
    dplyr::mutate(group=paste0('cl',grp))

  colors <- mtd %>%
    dplyr::select(paste0('col.',sel.col), group) %>%
    unique()

  colors <-colors[order(colors$group),]
  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,2])
  return(grid.colors)
}


#' Title
#'
#' @param sorted.cliques
#' @param columns
#' @param base.column
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
#' out.colors <- blend_cluster_colors(sorted.cliques, columns=c('cc0.8','cc0.85','cc0.9','cc0.95','cc0.98'), base.column = 'cc0.8', colors = c("#39811D","#76AECF", "#A6D8D4", "#F2A968", "#F2EC70", "#E38FDD", "#898989", "#86DA81", "#83731B", "#B34D22") )
blend_cluster_colors <- function(sorted.cliques, columns, base.column, colors ){
  base.colors <- NULL
  trans.colors <- NULL
  for(row in 1:nrow(sorted.cliques)){
    #sorted.cliques[row,columns]
    base.color <- colors[sorted.cliques[row,columns[1]]]
    base.colors = c(base.colors, base.color)

    color.tans.df = NULL
    i<-0
    t.cols = NULL
    for(col in columns[1:length(columns)]){
      message(col)
      i <- i + 1
      if(is.na(sorted.cliques[row,col])){
        #if(is.na(sorted.cliques[row,columns2[i]])){
        t.col <- NA
        message('assigned to NA as no group')
      }else{
        if(sorted.cliques[row,col] == 1){
          if(i==1){
            t.col <- base.color
            message('assigned to base as i =1')
          }else{
            if(!is.null(t.cols)){
              #previsously assigned
              t.col <- t.cols[length(t.cols)]
              message('assigned to previosu')

            }else{message('somthing went wrong')}
          }
        }else{
          if(!is.null(t.cols)){
            rgb.val <- col2rgb(t.cols[length(t.cols)])
          }else{
            rgb.val <- col2rgb(base.color)
          }

          r <- rgb.val[1]-(sorted.cliques[row,col]*7)
          #r <- rgb.val[1]+(sorted.cliques[row,columns2[i]]*10)
          if(r>255){r=255}
          if(r<0){r=0}

          g <- rgb.val[2]-(sorted.cliques[row,col]*10)
          #g <- rgb.val[2]+(sorted.cliques[row,columns2[i]]*10)
          if(g>255){g=255}
          if(g<0){g=0}

          b <- rgb.val[3]+(sorted.cliques[row,col]*10)
          #b <- rgb.val[3]+(sorted.cliques[row,columns2[i]]*10)
          if(b>255){b=255}
          if(b<0){b=0}

          t.col <- rgb(r, g, b, max = 255)
          message('assigned to kinky color')
        }
      }
      t.cols = c(t.cols,t.col)
    }
    message(t.cols)
    trans.colors <- rbind(trans.colors, t.cols)
  }
  out.colors <- cbind(trans.colors)
  colnames(out.colors) <- paste0('col.',columns)
  rownames(out.colors) <- sorted.cliques$genome
  out.colors <- out.colors %>% data.frame %>% rownames_to_column(var='genome')


  return(out.colors)
}

#============================================================================================#
# Simulated data
# ---
# Chi Squared test linking habitat vs ani clique
#table() function that cross-classifies the number of rows that are in the categories specified by the two categorical variables.
#The null hypothesis with this test is that the two categories are independent. The alternative hypothesis is that there is some dependency between the two categories.
# low p-value from the Chi-Squared test. Thus indicates that there is a very low probability that the two categories are independent. Therefore we reject the null hypothesis.
# low p, fairly strong relationship between A and B
# how does the relationship between A and B changes with more stringent any clustering

#make this better by allowing selection of metadata, and metadata colums, anyway it works for the purpose now
#chisq_range
#============================================================================================#
#' Title
#'
#' @param sorted_cl
#' @param sel.col
#'
#' @return
#' @export
#'
#' @examples
#'
chisq_range <- function(sorted_cl, sel.col){
  output <- NULL
  for(selcol in sel.col){
    mtd <- genome.tbl %>%
      dplyr::left_join(sorted_cl, by='genome') %>%
      dplyr::filter(!is.na(!!as.name(selcol)))  %>%
      dplyr::filter(!is.na(SS)) %>%
      data.frame()

    mtd$grp <- mtd[,selcol]

    mtd <- mtd %>%
      dplyr::mutate(group=paste0('cl',grp))

    chording <- mtd %>%
      dplyr::filter(!is.na(!!as.name(selcol))) %>%
      dplyr::mutate() %>%
      dplyr::filter(!is.na(SS))

    if(length(unique(chording$group))>1){
      contingency.mat <- table(as.character(chording$group),as.character(chording$SS))
      chisq.out <- chisq.test(contingency.mat)

      names(chisq.out)
      chisq.out$p.value
      chisq.out$statistic

      c.out <- c( 'Xsquared'=as.numeric(unname(chisq.out$statistic)),
                  'df' = unname(chisq.out$parameter),
                  'p.value'=chisq.out$p.value,
                  'method'=chisq.out$method,
                  'n_A' = nrow(contingency.mat),
                  'n_B' = ncol(contingency.mat),
                  'grouping' = selcol
      )
      output <- rbind(output, c.out)
    }
  }
  output <- output %>% base::data.frame(stringsAsFactors = FALSE)  %>%
    dplyr::mutate_at(c('Xsquared','df','p.value','n_A','n_B'), as.numeric) %>%
    dplyr::mutate(depth=gsub('sc','',grouping)) %>%
    dplyr::mutate_at(c('depth'), as.numeric)
  return(output)
}


# hsv_colors
#=========================================================================
#' Title TREEMAP based HIERARCHICAL COLOR ASSIGNMENT
#'
#' @param sorted.cl
#' @param cutoff_range
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
#' out.colors<-hsv_colors(sorted.cl,cutoff_range,colors)

hsv_colors <- function(sorted.cl,cutoff_range,colors){
  column.names = paste0('cc',cutoff_range)
  s.dat = sorted.cl[column.names]

  # TWO METHODS EXIST, HSV or HCL
  hsv.colors <- treemap::treepalette(s.dat,method='HSV',palette= colors)
  #hsv.colors <- treemap::treepalette(s.dat,method='HCL', palette.HCL.options =list('luminance'=10))
  #cols(unique(as.character(test.colors$HCL.color)))

  #Assign hrc (hierarchical color range)
  df2hrc <- function(df){
    c.hrc <- NULL
    for(i in 1:nrow(df)){
      c.hrc = c(c.hrc, paste(df[i,][!is.na(df[i,])],collapse='.'))
    }
    return(c.hrc)
  }

  chrc =  df2hrc(hsv.colors[,column.names])
  hsv.colors <- hsv.colors %>% mutate(hrc=chrc)

  #assign named vector for color mapping
  color.mapping <- hsv.colors %>% select(hrc,HSV.color)
  color.map <- as.character(color.mapping$HSV.color)
  names(color.map) <- as.character(color.mapping$hrc)

  # generate dataframe with the respective hrc values
  df.hrc <- NULL
  for(j in 1:ncol(s.dat)){
    chrc =  df2hrc(data.frame(s.dat[,1:j]))
    df.hrc <- cbind(df.hrc,chrc )
  }
  colnames(df.hrc) <- column.names

  #out format
  out.colors <- df.hrc
  out.colors[] <- as.character(color.mapping$HSV.color[match(unlist(df.hrc), color.mapping$hrc)])
  rownames(out.colors) <- sorted.cl$genome
  colnames(out.colors) <- paste0('col.',column.names)
  out.colors <- out.colors %>% data.frame() %>% rownames_to_column(var='genome')
  return(out.colors)
}

#plot.igraph2
#=====================================================================================
#' Title
#'
#' @param similarity_matrix
#' @param weight_threshold
#' @param community
#' @param remove_singletons
#' @param seed
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
#'
plot.igraph2 <- function(similarity_matrix, weight_threshold,community=NULL, remove_singletons = TRUE,seed=04012021,colors=NULL){
  set.seed(seed)
  n_genomes <- nrow(similarity_matrix)
  q <- igraph::graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "upper", weighted = T, diag = F)
  q <- igraph::simplify(q, remove.multiple=TRUE, remove.loops=TRUE)
  q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < weight_threshold)])

  if(remove_singletons == TRUE){
    q = delete.vertices(q, which(igraph::degree(q)< 2))
  }


  if(is.null(community)){
    wc <- cluster_infomap(q, e.weights = E(q)$weight)
    cliques = membership(wc)
    V(q)$cliques = cliques[names(V(q))]
    singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])==1])
    non_singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])>1])
    V(q)$comm = ifelse(cliques[names(V(q))] %in% singleton_cl, NA, cliques[names(V(q))] )
  }else{
    V(q)$comm = community[names(V(q))]
    message('loading previous cluster identifcation')
  }

  n_vertices <- length(V(q))
  qgn = ggnetwork::ggnetwork(intergraph::asNetwork(q),
                             layout="fruchtermanreingold",
                             weights="weight",
                             niter = 1000)

  p.out <- ggplot2::ggplot(data = qgn,
                           aes(x, y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(color = "black") + #,alpha=weight
    ggnetwork::geom_nodes( size = 3, shape=21, aes(fill=as.factor(comm)))+
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


#cols
#=====================================================================================
#' Title
#'
#' @param a
#'
#' @return
#' @export
#'
#' @examples
#'cols(getPalette(15))

cols <- function(a){
  image(1:length(a), 1, as.matrix(1:length(a)), col=a, axes=FALSE , xlab="", ylab="")
}




#=================================================================================#
#=================================================================================#
#=================================================================================#
#=================================================================================#






#=================================================================================#
#!!!!!!!!! THIS COULD BE RUN BOTH ON FULL PHYLOGENY AS WELL AS OUTGROUP
#=================================================================================#
library(igraph)
#=================================================================================#weight_threshold = 0.825
similarity_matrix = ANIb[sel.genome,sel.genome]
set.seed(1234)
n_genomes <- nrow(similarity_matrix)
q <- igraph::graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "upper", weighted = T, diag = F)
q <- igraph::simplify(q, remove.multiple=TRUE, remove.loops=TRUE)
q <-  igraph::delete_edges(q, E(q)[which(E(q)$weight < weight_threshold)])

wc <- cluster_infomap(q,e.weights = E(q)$weight)
wc <- igraph::cluster_fast_greedy(q, weights = E(q)$weight)

cliques = membership(wc)


# make cor matrix, here i use a pan-genome table - needs, samples as rows, and taxa as columns (in my case is this genomes as rows and genes as columns)
cor.matrix <- as.matrix(cor(OG.pan,method = 'pearson'))
cor.cutoff <- 0.3
cor.adj <- ifelse(abs(cor.matrix) >= cor.cutoff, 1, 0) # Construct microbiome network from adjacency matrix cor.net <- graph.adjacency(cor.adj,
q <- igraph::graph_from_adjacency_matrix(cor.adj, mode='undirected',weighted = T, diag = F)

wc <- cluster_infomap(q,e.weights = E(q)$weight)
communities = membership(wc)

#here i use the detected communities to extract some intersting things, like singletons etc
singleton_cl <- names(table(communities[names(V(q))])[table(communities[names(V(q))])==1])
non_singleton_cl <- names(table(communities[names(V(q))])[table(communities[names(V(q))])>1])

#you can pass data to the vertices in the network using V(q) and the dolar sign
# so here I just add the community name to a vertex if it is not a singleton
V(q)$comm = ifelse(communities[names(V(q))] %in% singleton_cl, NA, communities[names(V(q))] )

#count nr of vertecis
n_vertices <- length(V(q))

#use ggnetwork and intergraph together, else it does not work for me
qgn = ggnetwork::ggnetwork(intergraph::asNetwork(q),
                           layout="fruchtermanreingold",
                           weights="weight",
                           niter = 10000)
head(qgn)

#this is now easily plotted in ggplot
# you could also start annotating this better, merging dataframes with qgn object
ggplot2::ggplot(data = qgn,
                aes(x, y, xend = xend, yend = yend)) +
  ggnetwork::geom_edges(color = "grey50") +
  ggnetwork::geom_nodes( size = 2, shape=21, aes(fill=as.factor(comm)))+
  theme_void()+theme(aspect.ratio=1)





#Remove single nodes
#if(remove_singletons == TRUE){
#  q = delete.vertices(q, which(igraph::degree(q)< 2))
#}

wc <- cluster_infomap(q)
cliques = membership(wc)

#V(q)$L2FC = data.frame(selectedRES)[names(V(q)),'log2FoldChange']
V(q)$cliques = cliques[names(V(q))]

#identify singletons
singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])==1])
non_singleton_cl <- names(table(cliques[names(V(q))])[table(cliques[names(V(q))])>1])
V(q)$comm = ifelse(cliques[names(V(q))] %in% singleton_cl, NA, cliques[names(V(q))] )


n_vertices <- length(V(q))
qgn = ggnetwork::ggnetwork(intergraph::asNetwork(q),
                           layout="fruchtermanreingold",
                           weights="weight",
                           niter = 10000)

p.out <- ggplot2::ggplot(data = qgn,
                         aes(x, y, xend = xend, yend = yend)) +
  ggnetwork::geom_edges(color = "grey50") +
  ggnetwork::geom_nodes( size = 2, shape=21, aes(fill=as.factor(comm)))+
  ggtitle(paste0(weight_threshold,' (n=',n_vertices,')'))+
  theme_void()+theme(aspect.ratio=1)+
  scale_fill_manual(values=getPalette(length(non_singleton_cl)))


p.out


habitat.colors = c(
  'Hydrothermal_Vent'="#E41A1C", 'Lake'="#F37912", 'Oil'="#FFD422", 'Photic'="#43997A", 'Phototroph'="#658E67", 'Polar'="#5D6795", 'Sediment'="#A35390", 'Terrestial'="#B85F49",'V1'='black','Deep'='grey','Saltern'='red')



#---------------------------
#ANIb bewteen genomes were calculated using pyANI, with default parameters.
#Cliques are a complete graph wghere each vertex represents a genome and every vcertex is linked with evry other vertex, links are formed between genomes at a specific ANI threshold such as >= 95% (See Parks, species cliques)

#Using only those that have phylogenetic data
#--- ANIb.f is the phylo filtered ANIb table, see earlier
#ANIb <- ANIb.f
#AAI <- data.frame(AAI.f/100)


outgroup

#--- drop the outgroup in a second iteration
drop.outgroup = TRUE
ANIb <- ANIb[!(rownames(ANIb)%in%outgroup),!(rownames(ANIb)%in%outgroup)]
AAI <- AAI[!(rownames(AAI)%in%outgroup),!(rownames(AAI)%in%outgroup)]



#full network isolates
#sim16s2 =  sim16s[!grepl('Marinobacter_sp_FDB33|Alteromonas_sp_FDB36',rownames(sim16s)),!grepl('Marinobacter_sp_FDB33|Alteromonas_sp_FDB36',colnames(sim16s))]
#those from bioassay
#sim16s2 =  sim16s[rownames(trait.binairy),rownames(trait.binairy)]




#====

dim(ANIb)

ANIb.sel <- ANIb[sel.genome,sel.genome]
AAI.sel <- AAI[sel.genome,sel.genome]



tip.order <- tip_order(lsTrees.94.rooted[[1]])

ANIb.sel <- ANIb.sel[rev(tip.order),rev(tip.order)]
AAI.sel <- AAI.sel[rev(tip.order),rev(tip.order)]

GRR.sel <- pan.ggr[rev(tip.order),rev(tip.order)]
#AAI.sel <- AAI.sel[rev(tip.order),rev(tip.order)]


#In case you want heatmap.2
# heatmap.2(as.matrix((AAI.sel)),
#           trace="none",
#           col = viridis::viridis(20),
#           srtCol=45,
#           cexRow=0.4,
#           cexCol=0.5,
#           main='',
#           #margins=c(10,20),
#           Rowv=FALSE,
#           Colv=FALSE,
#           dendrogram-'none')
#
# heatmap.2(as.matrix((ANIb.sel)),
#           trace="none",
#           #col = viridis::viridis(20),
#           col = colorRampPalette(brewer.pal(9, "Blues"))(100),
#           srtCol=45,
#           cexRow=0.4,
#           cexCol=0.5,
#           main='',
#           #margins=c(10,20),
#           Rowv=FALSE,
#           Colv=FALSE,
#           dendrogram-'none')




#=========
# using complexheatmap

col_fun = circlize::colorRamp2(c(0.6,0.8, 1), c("deepskyblue2", "black", "yellow"))


metadata <- genome.tbl %>% filter(genome %in% sco.94.concat.tree.rooted$tip.label)# %>% data.frame()

plotdat <-genome.tbl %>%
  left_join(sorted.cliques,by=c('genome'='genome')) %>%
  left_join(out.colors, by=c('genome'='genome')) %>%
  #mutate(phylogroup=phylogroup.x) %>%
  filter(genome %in% sco.94.concat.tree.rooted$tip.label) %>% data.frame

rownames(plotdat) = plotdat$genome

#metadata <- metadata %>% mutate(phylogroup=phylogroup.x)
ph2gen <- metadata$phylogroup
names(ph2gen) <- metadata$genome

type2gen <- metadata$TypeStrain
names(type2gen) <- metadata$genome



th2 <- ComplexHeatmap::HeatmapAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                        phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                        ANI0.8 = paste0('cl',plotdat[rev(tip.order),'sc0.8']),
                                        ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),
                                        ANI0.9 = paste0('cl',plotdat[rev(tip.order),'sc0.9']),
                                        ANI0.95 = paste0('cl',plotdat[rev(tip.order),'sc0.95']),
                                        ANI0.98 = paste0('cl',plotdat[rev(tip.order),'sc0.98']),

                                        col=list(
                                          phylogroup = phylogroup.colors,
                                          ANI0.8 = col.list[[1]],
                                          ANI0.85 = col.list[[2]],
                                          ANI0.9 = col.list[[3]],
                                          ANI0.95 = col.list[[4]],
                                          ANI0.98 = col.list[[5]],
                                          type=c('TRUE'='black','FALSE'='white')#,
                                )
)

th3 <- ComplexHeatmap::rowAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                     phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                     ANI0.8 = paste0('cl',plotdat[rev(tip.order),'sc0.8']),
                                     ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),
                                     ANI0.9 = paste0('cl',plotdat[rev(tip.order),'sc0.9']),
                                     ANI0.95 = paste0('cl',plotdat[rev(tip.order),'sc0.95']),
                                     ANI0.98 = paste0('cl',plotdat[rev(tip.order),'sc0.98']),

                                     col=list(
                                       phylogroup = phylogroup.colors,
                                       ANI0.8 = col.list[[1]],
                                       ANI0.85 = col.list[[2]],
                                       ANI0.9 = col.list[[3]],
                                       ANI0.95 = col.list[[4]],
                                       ANI0.98 = col.list[[5]],
                                       type=c('TRUE'='black','FALSE'='white')#,
                                     )
)

th <- ComplexHeatmap::HeatmapAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                         ANI0.98 = paste0('cl',plotdat[rev(tip.order),'sc0.98']),
                                         ANI0.95 = paste0('cl',plotdat[rev(tip.order),'sc0.95']),
                                         ANI0.9 = paste0('cl',plotdat[rev(tip.order),'sc0.9']),
                                         ANI0.85 = paste0('cl',plotdat[rev(tip.order),'sc0.85']),
                                         ANI0.8 = paste0('cl',plotdat[rev(tip.order),'sc0.8']),
                                         phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                        col=list(
                                          phylogroup = phylogroup.colors,
                                          ANI0.8 = col.list[[1]],
                                          ANI0.85 = col.list[[2]],
                                          ANI0.9 = col.list[[3]],
                                          ANI0.95 = col.list[[4]],
                                          ANI0.98 = col.list[[5]],
                                          type=c('TRUE'='black','FALSE'='white')#,
                                        )
)

rh <- ComplexHeatmap::rowAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                        phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                        col=list(phylogroup = phylogroup.colors,
                                                 type=c('TRUE'='black','FALSE'='white')#,
                                        )
)

# For AAI
col_fun2 = circlize::colorRamp2(c(0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1), viridis::viridis(8))
col_fun = circlize::colorRamp2(c(0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1), colorRampPalette(brewer.pal(9, "Blues"))(8))


pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_AAI_upsidedown.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(AAI.sel),
                        top_annotation = th2,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE
)
dev.off()

pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_AAI.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(AAI.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE
                        )
dev.off()

pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_AAI_viridis.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(AAI.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()


# For ANIb
col_fun2 = circlize::colorRamp2(c(0.7,0.75,0.8,0.85,0.9,0.95, 1), viridis::viridis(7))
col_fun = circlize::colorRamp2(c(0.7,0.75,0.8,0.85,0.9,0.95, 1), colorRampPalette(brewer.pal(9, "Blues"))(7))


pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_ANIb.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(ANIb.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()

pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_ANIb_viridis.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(ANIb.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()




pdf( paste("~/DATA/MarbGenomics/Graphs/","complex_heatmap_GRR_viridis.pdf",sep=''), width=6, height=5)
ComplexHeatmap::Heatmap(as.matrix(GRR.sel),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE)
dev.off()


#================================
pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_full.pdf",sep=''), width=14, height=4)

ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), ]),
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap.pdf",sep=''), width=14, height=4)

ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_hgm]),
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()


#Viridis
pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_v.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_hgm]),
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()



#Viridis
pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_cazy_heatmap_v.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pan_clean(cazy.pan[rev(tip.order), ])),
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_COG_heatmap_v.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pan_clean(COG.pan[rev(tip.order), ])),
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE
)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_kfm_heatmap_v.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pan_clean(kfm.pan[rev(tip.order), ])),
                        left_annotation= rh,
                        col = col_fun2,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()




pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam__heatmap_v.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pan_clean(pfam.pan[rev(tip.order), sigpfam ])),
                        left_annotation= rh,
                        #col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)
dev.off()

ComplexHeatmap::Heatmap(as.matrix(pan_clean(pfam.pan[rev(tip.order), sig.genes ])),
                        left_annotation= rh,
                        #col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE
)



pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('Cob|Dct|amylase|TPR|Ank|Transposase',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)

sig.genes[grepl('Transposase',sig.genes)]


pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_HTH_hmg_heatmap.pdf",sep=''), width=7, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('HTH',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()



pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_allegaartje_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('Nap|rhodopsin|lactamas|Cyt|toxin|MCP|resist|Phage|transposase|Tnp',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_DUF_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('DUF',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10), c(colorRampPalette(brewer.pal(9, "Blues"))(10))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


pdf( paste("~/DATA/MarbGenomics/Graphs/","pfam_NAD_ATP_hmg_heatmap.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('ATP|NAD',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)
dev.off()


ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(sig.genes[grepl('Gly|aldo|argin|alg|aldedh|pqq|malt|Mtase',sig.genes,ignore.case=T)]) ]),
                        left_annotation= th3,
                        col = circlize::colorRamp2(seq(0,20,by=1), c(colorRampPalette(brewer.pal(9, "Blues"))(21))),
                        show_row_names = FALSE,
                        #show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE)



groups <- ANI.hmg.grp.tbl %>%
  filter(level=='sc0.85') %>% select(group) %>% unique() %>% pull %>% as.character() %>% sort

for(sel.grp in groups){
  message(sel.grp)
  ssp <- ANI.hmg.grp.tbl %>%
    filter(level=='sc0.85') %>%
    filter(fdr.p.value_depletion_binary<0.05) %>%
    filter(group==sel.grp) %>% select(feature.id) %>%
    pull %>% as.character()


  pdf( paste0("~/DATA/MarbGenomics/Graphs/","pfam_hmg_heatmap_de",sel.grp,".pdf"), width=0.2*length(ssp), height=4)
  print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), unique(ssp)]),
                          left_annotation= th3,
                          col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                          show_row_names = FALSE,
                          #show_column_names = FALSE,
                          cluster_rows = FALSE,
                          cluster_columns = TRUE))
  dev.off()
}



sel.grp='cl1'

mda <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(MeanDecreaseAccuracy) %>% pull()
names(mda) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

penrbin <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_enriched_binary) %>% pull() %>%
names(penrbin) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

penrbin <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_enriched_binary) %>% pull()
names(penrbin) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

penrraw <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_enriched_rawcounts) %>% pull()
names(penrraw) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

pdepbin <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.p.value_depletion_binary) %>% pull()
names(pdepbin) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()

pdepraw <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(fdr.value_depletion_rawcounts) %>% pull()
names(pdepraw) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()


sel.features <- sig.genes

ch <- ComplexHeatmap::HeatmapAnnotation(MeanDecrAcc = mda[sel.features],
                                        penrbin=-log10(penrbin[sel.features]),
                                        penrraw=-log10(penrraw[sel.features]),
                                        pdepbin=-log10(pdepbin[sel.features]),
                                        pdepraw=-log10(pdepraw[sel.features]),
                                        col=list(
                                          mda = circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          penrraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          pdepbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                          pdepraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")))
                                                                                  )


print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order),sel.features ]),
                              left_annotation= th3,
                              top_annotation= ch ,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



mda1 <- HGM_RF_COMB %>% filter(group.x=='cl1') %>% select(MeanDecreaseAccuracy) %>% pull()
mda2 <- HGM_RF_COMB %>% filter(group.x=='cl2') %>% select(MeanDecreaseAccuracy) %>% pull()
mda3 <- HGM_RF_COMB %>% filter(group.x=='cl3') %>% select(MeanDecreaseAccuracy) %>% pull()
mda4 <- HGM_RF_COMB %>% filter(group.x=='cl4') %>% select(MeanDecreaseAccuracy) %>% pull()
mda5 <- HGM_RF_COMB %>% filter(group.x=='cl5') %>% select(MeanDecreaseAccuracy) %>% pull()
mda6 <- HGM_RF_COMB %>% filter(group.x=='cl6') %>% select(MeanDecreaseAccuracy) %>% pull()
mda7 <- HGM_RF_COMB %>% filter(group.x=='cl7') %>% select(MeanDecreaseAccuracy) %>% pull()
mda8 <- HGM_RF_COMB %>% filter(group.x=='cl8') %>% select(MeanDecreaseAccuracy) %>% pull()
mda9 <- HGM_RF_COMB %>% filter(group.x=='cl9') %>% select(MeanDecreaseAccuracy) %>% pull()
mda10 <- HGM_RF_COMB %>% filter(group.x=='cl10') %>% select(MeanDecreaseAccuracy) %>% pull()
mda11 <- HGM_RF_COMB %>% filter(group.x=='cl11') %>% select(MeanDecreaseAccuracy) %>% pull()

names(mda1) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda2) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda3) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda4) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda5) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda6) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda7) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda8) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda9) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda10) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()
names(mda11) <- HGM_RF_COMB %>% filter(group.x==sel.grp) %>% select(feature.id.x)%>% pull()



ch <- ComplexHeatmap::HeatmapAnnotation(cl1=mda1[sel.features],
                                        cl2=mda2[sel.features],
                                        cl3=mda3[sel.features],
                                        cl4=mda4[sel.features],
                                        cl5=mda5[sel.features],
                                        cl6=mda6[sel.features],
                                        cl7=mda7[sel.features],
                                        cl8=mda8[sel.features],
                                        cl9=mda9[sel.features],
                                        cl10=mda10[sel.features],
                                        cl11=mda11[sel.features]
                                        #col=list(
                                        #  mda = circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  penrbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  penrraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  pdepbin=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")),
                                        #  pdepraw=circlize::colorRamp2(c(0, 1.30103, 3), c("white", "black", "firebrick")))
)


print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order),sel.features ]),
                              left_annotation= th3,
                              top_annotation= ch ,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), grepl('RCC|NLR|LRR',colnames(pfam.pan))]),
                              left_annotation= rh,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))


print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), grepl('Arm-DNA-bind_1|recombinase|resolvase',colnames(pfam.pan),ignore.case=T)]),
                              left_annotation= rh,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



print(ComplexHeatmap::Heatmap(as.matrix(pfam.pan[rev(tip.order), grepl('Arm-DNA-bind_1|recombinase|resolvase',colnames(pfam.pan),ignore.case=T)]),
                              left_annotation= rh,
                              col = circlize::colorRamp2(c(0,1,2,3,4,5,6,8,9,10,20,30,40,50,60,70,80,90,100), c(colorRampPalette(brewer.pal(9, "Blues"))(10),colorRampPalette(brewer.pal(9, "Oranges"))(9))),
                              show_row_names = FALSE,
                              #show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE))



#===
n_k <- 15
th <- ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 2:4),
                                         labels = LETTERS[1:n_k],
                                         labels_gp = grid::gpar(col = "white", fontsize = 10)))

pdf( paste("~/DATA/MarbGenomics/Graphs/","phylogenomics_heatmap_vlocks.pdf",sep=''), width=14, height=4)
ComplexHeatmap::Heatmap(as.matrix(OG.pan[rev(tip.order), allSig_hgm]),
                        top_annotation = th,
                        left_annotation= rh,
                        col = col_fun,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE,
                        column_km = n_k
)
dev.off()





col_fun2 = circlize::colorRamp2(c(0.7,0.75,0.8,0.85,0.9,0.95, 1), viridis::viridis(7))
col_fun = circlize::colorRamp2(c(0.7,0.75,0.8,0.85,0.9,0.95, 1), colorRampPalette(brewer.pal(9, "Blues"))(7))

col_fun3 = circlize::colorRamp2(seq(50,100,5),
                                viridis::viridis(length(seq(50,100,5))))

#======== MODUEL
n.pan <- m.pan[tip.order,]
n.pan <- n.pan[,colSums(n.pan)!=0]

rowMax <- function (df) {
  apply(df, MARGIN=c(1), max)
}
rowMax(n.pan)


colMax <- function (df) {
  apply(df, MARGIN=c(2), max)
}

selected.features <- names(colMax(n.pan)[colMax(n.pan)>85])


th <- ComplexHeatmap::HeatmapAnnotation(type = unname(as.character(type2gen[rev(tip.order)])),
                                        phylogroup = unname(as.character(ph2gen[rev(tip.order)])),
                                        col=list(phylogroup = phylogroup.colors,
                                                 type=c('TRUE'='black','FALSE'='white')
                                        ))
heatmap.matrix <- t(n.pan[rev(tip.order),selected.features ])
htm <- ComplexHeatmap::Heatmap(as.matrix(heatmap.matrix),
                        top_annotation =  th,
                        col = col_fun3,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = TRUE,
                        cluster_columns = FALSE,
                        clustering_distance_rows = function(x) vegan::vegdist(heatmap.matrix,method = 'jaccard'),
                        clustering_method_rows = 'ward.D2'
                        )
htm


highlighted_features <- m2n2 %>% filter(grepl('transport',DESCRIPTION,ignore.case=TRUE)) %>% select(module, DESCRIPTION) %>% unique() %>% filter(module %in% rownames(heatmap.matrix)) %>% pull(module)

selDesc <- c('Glucose/mannose transport system',
'Trehalose/maltose transport system',
'Putative amino-acid transport system',
'Putative multiple sugar transport system',
#'Fructose transport system',
'Glycine betaine/proline transport system',
'Osmoprotectant transport system',
"Glycerol transport system",
"Maltose/maltodextrin transport system",
"Branched-chain amino acid transport system" ,
"Phosphonate transport system" ,
"Urea transport system",
"D-Methionine transport system",
"Spermidine/putrescine transport system",
"L-Arabinose/lactose transport system",
"Erythritol transport system",
"Oligopeptide transport system",
"Putative sorbitol/mannitol transport system",
"ABC-2 type transport system",
"Iron(III) transport system",
"Iron(III) transport system",
"Peptides/nickel transport system",
"General L-amino acid transport system",
#"Putative thiamine transport system",
"Cystine transport system",
#"Sulfonate transport system",
"Arginine transport system",
"Putative multiple sugar transport system")


highlighted_features <- m2n2 %>% filter(DESCRIPTION%in% selDesc ) %>% select(module, DESCRIPTION) %>% unique() %>% filter(module %in% rownames(heatmap.matrix)) %>% pull(module)

m2desc <- m2n2$DESCRIPTION
names(m2desc) <- m2n2$module

genelabels <- ComplexHeatmap::rowAnnotation(Genes= ComplexHeatmap::anno_mark(at=match( highlighted_features,rownames(heatmap.matrix)),
                                                                                labels= paste0(highlighted_features," - ",m2desc[highlighted_features]),
                                                                                labels_gp = grid::gpar(fontsize=10),
                                                                                padding=0.1
                                                                                ))
pdf( paste("~/DATA/MarbGenomics/Graphs/","modules_heatmap_v.pdf",sep=''), width=8, height=4)
ComplexHeatmap::draw(htm+genelabels)
dev.off()




#===
m2cat <- m2n2$Module.Category
names(m2cat) <-  m2n2$module

m2tpe <- m2n2$TYPE
names(m2tpe) <-  m2n2$module


rsplit <- as.character(m2tpe[rownames(heatmap.matrix)])
rsplit[is.na(rsplit)]='other'

heatmap.matrix <- t(n.pan[tip.order,selected.features ])


htm <- ComplexHeatmap::Heatmap(as.matrix(heatmap.matrix),
                               top_annotation =  th,
                               col = col_fun3,
                               show_row_names = FALSE,
                               show_column_names = FALSE,
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               #clustering_distance_rows = function(x) vegan::vegdist(heatmap.matrix,method = 'jaccard'),
                               clustering_method_rows = 'ward.D2',
                               row_split = rsplit
)
htm



#plot.igraph
#==========================================================================================
#' Title
#'
#' @param similarity_matrix
#' @param weight_threshold
#' @param remove_singletons
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#' plot.igraph(similarity_matrix = ANIb, weight_threshold = 0.9, remove_singletons = TRUE)

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
                            # layout="hall",
                              weights="weight",
                             niter = 1000)

  #p.out <- ggplot2::ggplot(#data = qgn,
  #                         data= qgn %>% left_join(genome.tbl, by=c('vertex.names'='genome')),
  #                         aes(x, y, xend = xend, yend = yend)) +
  #  ggnetwork::geom_edges(color = "black",aes(size=weight)) + #,alpha=weight
  #  ggnetwork::geom_nodes( size = 2, shape=21, alpha=0.7,aes(fill=as.factor(comm)))+
  #  ggtitle(paste0(weight_threshold,' (n=',n_vertices,')'))+
  #  theme_void()+theme(aspect.ratio=1)+
  #  scale_size(range=c(0.1,2))

  p.out <- ggplot2::ggplot(#data = qgn,
    data= qgn %>% left_join(genome.tbl, by=c('vertex.names'='genome')),
    aes(x, y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(color = "black",aes(alpha=weight,size=weight)) + #,alpha=weight
    ggnetwork::geom_nodes( size = 3, aes(shape=TypeStrain,fill=phylogroup))+
    ggtitle(paste0(weight_threshold,' (n=',n_vertices,')'))+
    theme_void()+theme(aspect.ratio=1)+
    scale_size(range=c(0.1,1)) + scale_shape_manual(values=c(21,23))#+scale_fill_manual(values=phylogroup.colors,na.value = 'white')

  if(is.null(colors)){
    p.out <- p.out + scale_fill_manual(values=getPalette(length(non_singleton_cl)),na.value='white')
  }else{
    p.out <- p.out + scale_fill_manual(values=colors,na.value='white')
  }

  p.out
}


p.igraph  <- plot.igraph(similarity_matrix = ANIb[sco.94.concat.tree$tip.label,sco.94.concat.tree$tip.label], weight_threshold = 0.9, remove_singletons = FALSE,colors=phylogroup.colors)


#run over several cutoffs and store plots in list
p.list <- list()
dist.name='OF'
for(ctf in c(0.7, 0.75,0.77,0.78,0.79, 0.8,0.85,0.9,0.95,0.98)){
  message(paste0('visualising cliques at >',ctf))
  p.list[[as.character(ctf)]] <- plot.igraph(similarity_matrix = pan.ggr, weight_threshold = ctf, remove_singletons = FALSE) + theme(legend.position = "none")
}
p.list <- list()
for(ctf in c(0.8, 0.81,0.82,0.83,0.84, 0.85,0.86,0.87,0.89,0.9)){
  message(paste0('visualising cliques at >',ctf))
  p.list[[as.character(ctf)]] <- plot.igraph(similarity_matrix = pan.ggr[sel.genome,sel.genome], weight_threshold = ctf, remove_singletons = FALSE) + theme(legend.position = "none")
}


p.list <- list()
dist.name='ANIb'
for(ctf in c(0.7, 0.75,0.77,0.78, 0.8,0.85,0.9,0.95,0.98)){
  message(paste0('visualising cliques at >',ctf))
  p.list[[as.character(ctf)]] <- plot.igraph(similarity_matrix = ANIb[sco.94.concat.tree$tip.label,sco.94.concat.tree$tip.label],
                                             weight_threshold = ctf,
                                             remove_singletons = FALSE,colors=phylogroup.colors) + fdb_style() +theme(legend.position = "none")
}

p.comb <- grid.arrange(grobs = p.list,nrow=3)
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph2','.pdf'),plot=p.comb, width = 9, height = 4,unit='in')
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph_test','.pdf'),plot=p.comb, width = 10, height = 10,unit='in')
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph_test2','.pdf'),plot=p.comb, width = 15, height = 15,unit='in')



p.list <- list()
dist.name='GRR'
for(ctf in c(0.7, 0.75,0.77,0.78,0.82, 0.8,0.85,0.9,0.95)){
  message(paste0('visualising cliques at >',ctf))
  p.list[[as.character(ctf)]] <- plot.igraph(similarity_matrix = GRR.sel[sco.94.concat.tree$tip.label,sco.94.concat.tree$tip.label],
                                             weight_threshold = ctf,
                                             remove_singletons = FALSE,colors=phylogroup.colors) + fdb_style() +theme(legend.position = "none")
}

p.comb <- grid.arrange(grobs = p.list,nrow=3)
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph2','.pdf'),plot=p.comb, width = 9, height = 4,unit='in')
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph_test','.pdf'),plot=p.comb, width = 10, height = 10,unit='in')
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph_test2','.pdf'),plot=p.comb, width = 15, height = 15,unit='in')



p.list <- list()
dist.name='AAI'
for(ctf in c(0.7, 0.75,0.77,0.78, 0.8,0.85,0.9,0.95,0.98)){
  message(paste0('visualising cliques at >',ctf))
  p.list[[as.character(ctf)]] <- plot.igraph(similarity_matrix = AAI[sco.94.concat.tree$tip.label,sco.94.concat.tree$tip.label],
                                             weight_threshold = ctf,
                                             remove_singletons = FALSE,colors=phylogroup.colors) + fdb_style() +theme(legend.position = "none")
}

p.comb <- grid.arrange(grobs = p.list,nrow=3)
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph2','.pdf'),plot=p.comb, width = 9, height = 4,unit='in')
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph_test','.pdf'),plot=p.comb, width = 10, height = 10,unit='in')
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph_test2','.pdf'),plot=p.comb, width = 15, height = 15,unit='in')





p.list <- list()
dist.name='COPH'
for(ctf in c(0.5, 0.55,0.6,0.65,0.75, 0.8,0.85,0.9,0.95,0.98)){
  message(paste0('visualising cliques at >',ctf))
  p.list[[as.character(ctf)]] <- plot.igraph(similarity_matrix = 1-COPH[sel.genome,sel.genome], weight_threshold = ctf, remove_singletons = FALSE) + theme(legend.position = "none")
}

p.comb <- grid.arrange(grobs = p.list,nrow=2)
ggplot2::ggsave(filename=paste0('~/DATA/MarbGenomics/Graphs/',dist.name,'_cliques_igraph2','.pdf'),plot=p.comb, width = 15, height = 7,unit='in')





#===================================================================
#===================================================================
#===================================================================
#!!!! REMOVE?  ?????
#!!!! REMOVE?  ?????
#!!!! REMOVE?  ?????
#!!!! REMOVE?  ?????

cutoff_range = seq(0.75, 1, by = 0.005)

multi.cl <- multiClique(distance=sim.list[['ANIb']],
                       cutoffs = cutoff_range,
                       no_singletons = TRUE,
                       community.method = 'cluster_infomap' )

sorted.cl <- sortCliques(multi.cl, tree=tr.tree)

out.colors <- blend_cluster_colors(sorted.cliques = sorted.cl ,
                                   columns = paste0('cc',cutoff_range),
                                   base.column = 'cc0.75',
                                   colors = c("#39811D","#76AECF", "#A6D8D4", "#F2A968", "#F2EC70", "#E38FDD", "#898989", "#86DA81", "#83731B", "#B34D22"))

p.anib <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.cl, circular=FALSE,out.colors = out.colors)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/testcolors_hierarchical3.pdf',plot=p.anib, width = 8, height = 8,unit='in')

cliquelevel2colors(sorted.cliques = sorted.cl,
                   out.colors = out.colors,
                   sel.col='cc0.85')



#=====================================================================================
#=====================================================================================
# USING HIERARCHICAL CLUSTERING TO OBTAIN COLORS, AND THEN BUILD TREE AND NETWORK PLOTS IN SAME SCHEME
# TREEMAP based HIERARCHICAL COLOR ASSIGNMENT

#cutoff_range = seq(0.785, 1, by = 0.005)
cutoff_range = seq(0.8, 1, by = 0.05)
cutoff_range = seq(0.8, 1, by = 0.05)

cuttoff_list <- list("COPH"=seq(0.75, 1, by = (1-0.75)/20),
                     "ANIb"=seq(0.78, 1, by = (1-0.78)/20),
                     "AAI"=seq(0.8, 1, by = (1-0.8)/20),
                     "OF"=seq(0.775, 1, by = (1-0.775)/20),
                     "GGR"=seq(0.835, 1, by = (1-0.835)/20))

for(sim.name in names(sim.list)){
  multi.cl <- multiClique(distance=sim.list[[sim.name]],
                          cutoffs = cuttoff_list[[sim.name]],
                          no_singletons = TRUE,
                          community.method = 'cluster_infomap' )

  sorted.cl <- sortCliques(multi.cl, tree=tr.tree)
  out.colors<-hsv_colors(sorted.cl,cuttoff_list[[sim.name]],colors)

  #plot trees
  p.cls <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.cl, circular=FALSE,out.colors = out.colors)
  ggplot2::ggsave(filename=paste0(outdir,sim.name,'_hierarchical_hsv_clustering_at_0.05inrements','.pdf'),plot=p.cls, width = 4, height = 8,unit='in')


  p.cls <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.cl[,c('genome',paste0('sc',cuttoff_list[[sim.name]][1]),paste0('sc',cuttoff_list[[sim.name]][1]))], circular=FALSE,out.colors = out.colors[,c('genome',paste0('col.cc',cuttoff_list[[sim.name]][1]))])
  ggplot2::ggsave(filename=paste0(outdir,sim.name,'_hierarchical_hsv_clustering_lowestlevel',cuttoff_list[[sim.name]][1],'.pdf'),plot=p.cls, width = 4, height = 8,unit='in')

  #plot the community clusters
  p.list <- list()
  for(ctf in cuttoff_list[[sim.name]]){
    if(ctf<1){
      tmp.col <- sorted.cl %>% left_join(out.colors,by='genome') %>% select(paste0('sc',ctf),paste0('col.cc',ctf))
      tmp.col <- tmp.col[!is.na(tmp.col[,1]),] %>% unique()
      tmp.col <- as.character(unname(tmp.col[,2]))
      names(tmp.col) <- 1:length(tmp.col)

      c.commnities <- sorted.cl %>% left_join(out.colors,by='genome') %>% select(paste0('sc',ctf)) %>% pull()
      names(c.commnities) = sorted.cl$genome
      p.list[[as.character(ctf)]] <- plot.igraph2(similarity_matrix = sim.list[[sim.name]], weight_threshold = ctf,community=c.commnities, remove_singletons = FALSE,colors=tmp.col)# + theme(legend.position = "none")
    }
  }
  p.comb <- grid.arrange(grobs = p.list,nrow=4)
  ggplot2::ggsave(filename=paste0(outdir,sim.name,'_hierarchical_hsv_clustering_at_0.05NET','.pdf'),plot=p.comb, width = 25, height = 25,unit='in')
}



#=============================================================================================================





#============================================================================================================================
# PARTIALLY DEVELOPED BUT OBSOLETE

#colors = getPalette(15)
#column.names <- paste0('cc',cutoff_range)
#sorted.cliques = sorted.cl
#colortracking = NULL
#i = 0
#firstsplit = FALSE
#for(colmn.name in column.names){
#  i <- i + 1
#  tmp <- sorted.cliques %>% select(genome,colmn.name)
#  t.cols=rep('lelijk',nrow(tmp))
#
#  if(i==1){
#    t.cols <- colors[tmp[,2]]
#    #tmp <- tmp %>% mutate(cols = colors[tmp[,2]])
#  }else{
#    if(firstsplit==FALSE){
#      firstsplit=TRUE
#      t.cols <- colors[tmp[,2]]
#    }else{
#      t.cols <- ifelse(tmp[,2] ==1, colortracking[,i-1],getPalette(30)[tmp[,2]])


#      }
    #tmp %>% mutate(cols = colors[tmp[,2]])
#  }
#  message(t.cols)
#  colortracking <- cbind(colortracking, t.cols)
#}

#colnames(colortracking) <- paste0('col.',column.names)
#rownames(colortracking) <- tmp$genome
#out.colors <- colortracking %>% data.frame() %>% rownames_to_column(var='genome')
#if first column, assign colors based on #original colors
#p.anib <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.cl, circular=FALSE,out.colors = out.colors)
#ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/testcolors_hierarchical3.pdf',plot=p.anib, width = 8, height = 8,unit='in')



#=====================================================================================
# BROADER SCALE SEARCH, 0.65 to 1 by 0.005
#
#===================================================================

cutoff_range <- seq(0.65, 1, by = 0.005)
similarity.list <- list('ANI'= ANIb,
                        'AAI'=AAI,
                        'COPH'=1-COPH[rownames(ANIb), rownames(ANIb)])

out.df<-NULL
for(i in 1:length(similarity.list)){
  int.df <- range.igraph(similarity_matrix = similarity.list[[i]], cutoff_range = cutoff_range )
  out.df = rbind(out.df, cbind(int.df,'matrix'=names(similarity.list)[i]))
  }

out.df <- out.df %>% data.frame()
out.df$ctf <- as.numeric(as.character(out.df$ctf ))
out.df$n <- as.numeric(as.character(out.df$n))
p.ANI_v_n <- out.df %>% ggplot(aes(ctf, n,color=matrix))+
  geom_hline(yintercept = 50, size = 0.25, colour = "#bdbdbd") +
  geom_vline(xintercept = 0.95, size = 0.25, colour = "#bdbdbd") +
  geom_vline(xintercept = 0.85, size = 0.25, colour = "#bdbdbd") +
  geom_line()+
  xlab('similarity cutoff')+
  ylab('n genomes after singleton removal')+
  fdb_style(aspect.ratio = 0.4)

out.df %>%
  filter(matrix=='ANI') %>%
  ggplot(aes(ctf, n))+
  geom_bar(stat='identity',aes(fill=matrix),color='white',size=.3)+
  xlab('similarity cutoff')+
  ylab('Genomes covered (%)')+
  fdb_style(aspect.ratio = 0.4)+facet_wrap(~matrix,ncol=1)

#Run this once from 0.7 to 1.
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANIb_cliques_v_n.pdf',plot=p.ANI_v_n, width = 4, height = 3,unit='in')

#and make new plot
p.ANI_v_n <- out.df %>%
  mutate(phylogenetic = ifelse(matrix=="COPH",'phylogenetic','phyloindependent')) %>%
  ggplot(aes(ctf, n,color=matrix))+
  geom_hline(yintercept = 50, size = 0.25, colour = "#bdbdbd") +
  geom_vline(xintercept = 0.95, size = 0.25, colour = "#bdbdbd") +
  geom_vline(xintercept = 0.85, size = 0.25, colour = "#bdbdbd") +
  geom_line()+
  xlab('similarity cutoff')+
  ylab('n genomes after singleton removal')+
  facet_wrap(~phylogenetic,scales = 'free')+
  fdb_style(aspect.ratio = 0.4)

ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/facet_wide_cliques_v_n.pdf',plot=p.ANI_v_n, width = 8, height = 3,unit='in')



#!!!!!===== COMMUNITY DETECTION THE INTENSIVE WAY ####
# In this section of code, we compile analysis for clique detection
#using different simmilarity measures, and different community detection methods

# RUNS OVER ALL SIMLIST SIMILARITY VALUES
#sim.list
dst.name = 'GGR'
cutoff_range <- seq(0.65, 1, by = 0.005)

tree=sco.94.concat.tree.rooted

sim.list <- list(
  'COPH' = 1-COPH[tree$tip.label,tree$tip.label],
  'ANIb' = ANIb[tree$tip.label,tree$tip.label],
  'AAI' = AAI[tree$tip.label,tree$tip.label],
  'OF' =  Ortho.Fraction[tree$tip.label,tree$tip.label],
  'GGR' = pan.ggr[tree$tip.label,tree$tip.label]
)


#Select genomes of interest from the genome.tbl
sel.genome <- genome.tbl %>%
  dplyr::filter(manualy_removed ==FALSE) %>%
  dplyr::filter(quality_class != 'low') %>%
  dplyr::filter(genus=='Marinobacter') %>%
  dplyr::select(genome) %>%
  dplyr::pull() %>%
  as.character()

#Reduce the original tree to only include the selected genomes
tr.tree <- ape::drop.tip(tree, tip = setdiff(tree$tip.label, sel.genome))


#Create a list of similarity measures, and filter to include only selected genomes
sim.list <- list(
  'COPH' = 1-COPH[sel.genome,sel.genome],
  'ANIb' = ANIb[sel.genome,sel.genome],
  'AAI' = AAI[sel.genome,sel.genome],
  'OF' =  Ortho.Fraction[sel.genome,sel.genome],
  'GGR' = pan.ggr[sel.genome,sel.genome]
)

#Set parameters of the analysis
#Names of the simmilarity values
set.seed(12345)
nms <- names(sim.list)
outdir = "~/DATA/MarbGenomics/Graphs/"
community_methods = c('cluster_edge_betweenness',
                      'cluster_fast_greedy',
                      'cluster_infomap',
                      'cluster_label_prop')
#'cluster_walktrap') #'cluster_spinglass'

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#getPalette = colorRampPalette(brewer.pal(9, "Set1")[1:length(brewer.pal(9, "Set1"))-3])

n_genomes <- length(rownames(sim.list[[1]]))


#initialising variables
out.out.tbl <- NULL
out.mono.checked.tbl <- NULL
chisq.out.tbl <- NULL
tr.tree <- tree

for(dst.name in nms){
  message('analysing ',dst.name)
  out.tbl <- NULL
  mono.checked.tbl <- NULL
  chisq.tbl <- NULL
    for(com.meth in community_methods){
    message('analysing ',com.meth)
    multi.cl <- multiClique(distance=sim.list[[dst.name]], cutoffs = cutoff_range,no_singletons = TRUE, community.method = com.meth )
    sorted.cl <- sortCliques(multi.cl, tree=tr.tree)
    sorted_m.cl <- sorted.cl %>% select(genome, sc0.65:sc1 )
    sorted_m.cl <- sorted_m.cl%>% melt(id="genome")
    sorted_m.cl$clique <- factor(sorted_m.cl$value, levels=c(NA,sort(unique(sorted_m.cl$value))), exclude = NULL)
    sorted_m.cl <- sorted_m.cl %>% mutate(method=com.meth)
    out.tbl <- rbind(out.tbl,sorted_m.cl)

    # STEP 2 Check monophyletic
    #convert to monophyly analysis input format
    monodat <- multi.cl %>% mutate(taxa = genome) %>%
      select(-genome) %>%
      select(taxa, c0.65:c1)
    mono.checked <- check.monoph(tr.tree ,monodat) #%>% mutate(tree=names(lsTrees.rooted)[1])
    mono.checked <- mono.checked %>% mutate(method=com.meth)
    mono.checked.tbl <- rbind(mono.checked.tbl,mono.checked)

    # STEP 3 analyse cluster overlaps with isolation sourcr
    # this step runs chisq_range function
    sorted_cl <- sorted.cl
    sel.col <- sorted_cl %>% select(sc0.65:sc0.995) %>% colnames()
    chisq.out <- chisq_range(sorted_cl, sel.col)

    #store output
    chisq.out <- chisq.out %>% mutate(com.method=com.meth)
    chisq.tbl <- rbind(chisq.tbl,chisq.out)
  }

  out.tbl$clique <- factor(out.tbl$value, levels=c(NA,sort(unique(out.tbl$value))), exclude = NULL)

  #Store dist name and make final output tibble
  out.tbl <- out.tbl %>% mutate(type=dst.name)
  out.out.tbl <- rbind(out.out.tbl, out.tbl)

  mono.checked.tbl <- mono.checked.tbl %>% mutate(type=dst.name)
  out.mono.checked.tbl <- rbind(out.mono.checked.tbl,mono.checked.tbl)

  chisq.tbl <- chisq.tbl %>% mutate(type=dst.name)
  chisq.out.tbl <- rbind(chisq.out.tbl,chisq.tbl)
}
out.out.tbl$clique <- factor(out.out.tbl$value, levels=c(NA,sort(unique(out.out.tbl$value))), exclude = NULL)
out.out.tbl$type <- factor(out.out.tbl$type, levels=c('ANIb','AAI','COPH','GGR'))

p.clique.dist <- out.out.tbl %>% mutate(simmilarity = as.numeric(gsub('sc','',variable))) %>%
  ggplot(aes(factor(simmilarity)))+
  geom_bar(aes(fill=clique), position='fill', width = 1)+
  scale_fill_manual(values=getPalette(length(sort(unique(out.out.tbl$value)))),na.value="white")+
  xlab('similarity cutoff')+
  ylab('Genomes covered (%)')+
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(method~type) +
  scale_x_discrete(breaks=c('0.7','0.75','0.8','0.85','0.9','0.95','c1'))+
  fdb_style(aspect.ratio = 0.4)
ggplot2::ggsave(filename=paste0(outdir,'_n=',n_genomes,'_range_communities_all_sorted','.pdf'),plot=p.clique.dist, width = 10, height = 8,unit='in')

p.clique.dist2 <- out.out.tbl %>% mutate(simmilarity = as.numeric(gsub('sc','',variable))) %>%
  filter(simmilarity %in% seq(0.75,0.95,by=0.05)) %>%
  ggplot(aes(factor(simmilarity)))+
  geom_bar(aes(fill=clique), position='fill', width = 1)+
  scale_fill_manual(values=getPalette(length(sort(unique(out.out.tbl$value)))),na.value="white")+
  xlab('similarity cutoff')+
  ylab('Genomes covered (%)')+
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(method~type) +
  fdb_style(aspect.ratio = 2)
ggplot2::ggsave(filename=paste0(outdir,'_n=',n_genomes,'_range_communities_selection_sorted','.pdf'),plot=p.clique.dist2, width = 8, height = 8,unit='in')

p.clique.dist3 <- out.out.tbl %>%
  group_by(variable,value,method,type) %>%
  dplyr::count() %>%
  ungroup() %>%
  group_by(variable,method,type) %>%
  dplyr::count() %>%
  ggplot(aes(variable,n))+
  geom_bar(stat='identity', fill='grey50',width = 1)+
  xlab('similarity cutoff')+
  ylab('# of cliques')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks=c('sc0.7','sc0.75','sc0.8','sc0.85','sc0.9','sc0.95','c1'))+
  facet_grid(method~type) +
  fdb_style(aspect.ratio = 0.4)
ggplot2::ggsave(filename=paste0(outdir,'_n=',n_genomes,'_range_N_communities','.pdf'),plot=p.clique.dist3, width = 10, height = 8,unit='in')

p.clique.dist4 <- out.out.tbl %>%
  group_by(variable,value,method,type) %>%
  dplyr::count() %>%
  ungroup() %>%
  group_by(variable,method,type) %>%
  dplyr::count() %>%
  mutate(simmilarity = as.numeric(gsub('sc','',variable))) %>%
  ggplot(aes(simmilarity,n,color=method))+
  geom_line()+
  xlab('similarity cutoff')+
  ylab('# of cliques')+
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~type,ncol=1) +
  scale_color_brewer(direction=1,palette="Dark2")+
  fdb_style(aspect.ratio = 0.4)
ggplot2::ggsave(filename=paste0(outdir,'_n=',n_genomes,'_range_overlay_N_communities','.pdf'),plot=p.clique.dist4, width = 6, height = 3,unit='in')

#Percentage of communities that is filles
p.monoph <- out.mono.checked.tbl %>% group_by(level,method,taxon,mono,type) %>% tally() %>%
  ggplot(aes(factor(level)))+
  geom_bar(aes(fill=mono), position='fill', width = 1)+
  xlab('similarity cutoff')+
  ylab('% of communities')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c('grey40','grey60'))+
  scale_x_discrete(breaks=c('c0.7','c0.75','c0.8','c0.85','c0.9','c0.95','c1'))+
  facet_grid(method~type) +
  fdb_style(aspect.ratio = 0.4)
ggplot2::ggsave(filename=paste0(outdir,'_n=',n_genomes,'_range_monophyletic_communities_all','.pdf'),plot=p.monoph, width = 10, height = 8,unit='in')

p.monoph <- out.mono.checked.tbl %>% group_by(level,method,type,taxon,mono) %>%
  tally()  %>%
  ggplot(aes(factor(level)))+
  geom_bar(aes(fill=mono), width = 1)+
  xlab('similarity cutoff')+
  ylab('# communities')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c('grey40','grey60'))+
  scale_x_discrete(breaks=c('c0.7','c0.75','c0.8','c0.85','c0.9','c0.95'))+
  facet_grid(method~type) +
  fdb_style(aspect.ratio = 0.4)
ggplot2::ggsave(filename=paste0(outdir,dst.name,'_n=',n_genomes,'_range_monophyletic_communities_count_all','.pdf'),plot=p.monoph, width = 10, height = 8,unit='in')

#Make sure we use this corretely, like faceting etc
p.chi <- chisq.out.tbl %>% ggplot(aes(x=as.factor(depth),y=-log10(p.value))) +
  geom_bar(stat = 'identity',aes(fill=df), width = 1)+
  ylab('-log10(p-value)')+
  xlab('ANI clique threshold')+
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05),color='grey60',size=0.5)+
  geom_hline(yintercept = -log10(0.01),color='grey60',size=0.5,linetype='dashed')+
  geom_hline(yintercept = -log10(0.001),color='grey60',size=0.5,linetype='dotted')+
  scale_x_discrete(breaks=c('sc0.7','0.75','0.8','0.85','0.9','0.95'))+
  facet_grid(com.method~type) +
  fdb_style(aspect.ratio = .5)
ggplot2::ggsave(filename=paste0(outdir,'chisq_NNIcliques_habitat_pvalue','.pdf'),plot=p.chi, width = 8, height = 3,unit='in')

p.chi <- chisq.out.tbl %>% ggplot(aes(x=as.factor(depth),y=Xsquared)) +
  geom_bar(stat = 'identity',aes(fill=df), width = 1)+
  ylab('Xsq')+
  xlab('ANI clique threshold')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks=c('sc0.7','0.75','0.8','0.85','0.9','0.95'))+
  facet_grid(com.method~type) +
  fdb_style(aspect.ratio = .5)
ggplot2::ggsave(filename=paste0(outdir,'chisq_NNIcliques_habitat_Chisq','.pdf'),plot=p.chi, width = 8, height = 3,unit='in')


#--------------------------------------------------------------------------------------------#
# Monophyly percentage plots

out.mono.checked.tbl %>%
  mutate(mono=as.logical(mono)) %>%
   group_by(level,method,type,taxon,mono) %>%
  tally()  %>%
    group_by(level,method,type) %>%
  filter(if(!all(mono)) TRUE else mono)

out.mono.checked.tbl %>%
  mutate(mono=as.logical(mono)) %>%
  group_by(level,method,type) %>%
  filter(!all(mono))

out.mono.checked.tbl %>% group_by(level,method,type,taxon,mono) %>%
  tally()  %>%
  ggplot(aes(factor(level)))+
  geom_bar(aes(fill=mono), width = 1)+
  xlab('similarity cutoff')+
  ylab('# communities')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c('grey40','grey60'))+
  scale_x_discrete(breaks=c('c0.7','c0.75','c0.8','c0.85','c0.9','c0.95'))+
  facet_grid(method~type) +
  fdb_style(aspect.ratio = 0.4)


#--------------------------------------------------------------------------------------------#
# condensed plot, total percentages of TRUE and FALSE monophely

p1 <- out.mono.checked.tbl %>%
  mutate(mono=as.logical(mono)) %>%
  mutate(method=gsub('cluster_','',method)) %>%
  ggplot(aes(method, fill=mono)) +
    geom_bar(position = 'fill',width = 0.9) +
    geom_text(data = . %>%
            group_by(method, type, mono) %>%
            group_by(method,type, mono) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::mutate(freq = n / sum(n)) %>% ungroup,
            aes(y = freq, label = scales::percent(freq)),
            position = position_fill(vjust = 0.5),
            show.legend = FALSE,size=1.5)+
  scale_fill_manual(values=c('grey40','grey60'))

#nrow=1
p1.1 <- p1 +
  facet_wrap(~type,nrow=1)+
  fdb_style(aspect.ratio = 1.5)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(outdir,'percentages_monophyletic_2','.pdf'),plot=p1.1, width = 8, height = 3,unit='in')

#ncol=1
p1.2 <- p1 + facet_wrap(~type,ncol=1)+
  fdb_style()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(outdir,'percentages_monophyletic2','.pdf'),plot=p1.2, width = 4, height = 8,unit='in')






#---------------------------------------------------
# building blocks for the composite image


getPalette_pl = colorRampPalette(brewer.pal(9, "Set1"))
getPalette_pl = colorRampPalette(brewer.pal(8, "Greys"))

getPalette_pl = colorRampPalette(brewer.pal(8, "BrBG"))
getPalette_pl2 = colorRampPalette(brewer.pal(8, "PRGn"))


#set aspect ratio the same
asp.rat <- 0.33

p1 <- intraPG %>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(value,fill=grouping)) +
  geom_density(alpha=0.5)+
  geom_hline(yintercept=0)+
  geom_point(aes(x=value,y=(as.numeric(as.factor(grouping))*5)-13,color=grouping),
             fill=NA,shape=3,alpha=0.5,
             position=position_jitter(w=0,h=1))+
  fdb_style(aspect.ratio = asp.rat) +
  scale_fill_manual(values=rev(getPalette_pl2(2))) +
  scale_color_manual(values=rev(getPalette_pl2(2)))+
  scale_x_continuous(expand=c(0,0),limits=c(0.7,1),breaks = c(.7,.75,0.8,0.85,0.9,0.95))+
  ylab('Density')

#plot with ridges per phylogroup
p2 <- intraGenus%>%
  dplyr::filter(itself==FALSE) %>%
  dplyr::filter(!is.na(grouping)) %>%
  ggplot(aes(x= value , y= group.x)) +
  ggridges::geom_density_ridges(mapping = aes(group=interaction(group.x, grouping), fill=grouping),alpha=0.5) +
  geom_point(aes(color=grouping),
             fill=NA,shape=3,alpha=0.5,
             position=position_jitter(w=0,h=0.1))+
  scale_color_manual(values=rev(getPalette_pl2(2)))+
  scale_fill_manual(values=rev(getPalette_pl2(2)))+
  fdb_style(aspect.ratio = asp.rat)+
  scale_y_discrete(limits=rev(levels(intraGenus$group.x)))+
  scale_x_continuous(expand=c(0,0),limits=c(0.7,1),breaks = c(.7,.75,0.8,0.85,0.9,0.95))+
  ylab('Phylogroup')


p3 <-  out.out.tbl %>% mutate(simmilarity = as.numeric(gsub('sc','',variable))) %>%
  filter(method=='cluster_infomap') %>%
  filter(simmilarity>=0.7) %>%
  ggplot(aes(factor(simmilarity)))+
  geom_bar(aes(fill=clique), position='fill', width = 1)+
  scale_fill_manual(values=rev(getPalette_pl(length(sort(unique(out.out.tbl$value))))),na.value="white")+
  xlab('similarity cutoff')+
  ylab('Genomes covered (%)')+
  scale_y_continuous(expand = c(0, 0)) +
  #facet_grid(method~type) +
  scale_x_discrete(breaks=c('0.7','0.75','0.8','0.85','0.9','0.95','c1'))+
  fdb_style(aspect.ratio = asp.rat)


p4 <- out.mono.checked.tbl %>%
  filter(method=='cluster_infomap') %>%
  mutate(simmilarity = as.numeric(gsub('c','',level))) %>%
  filter(simmilarity>=0.7) %>%
  group_by(simmilarity,method,type,taxon,mono) %>%
  tally()  %>%
  ggplot(aes(factor(simmilarity)))+
  geom_bar(aes(fill=mono), width = 1)+
  xlab('similarity cutoff')+
  ylab('# communities')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c('grey40','grey60'))+
  scale_x_discrete(breaks=c('0.7','0.75','0.8','0.85','0.9','0.95'))+
  #facet_grid(method~type) +
  fdb_style(aspect.ratio = asp.rat)



#cowplot::plot_grid(p1,p2,p3,p4,align='v',axis='rl',ncol=1,rel_heights = c(0.5,0.5,1,0.5,1))

pp.o <-ggpubr::ggarrange(p1,p2,p3,p4,labels=c('a','b','c','d'),ncol=1,align='v')#,heights = c(1,1.5))
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_breakup.pdf',plot=pp.o, width = 8, height = 8,unit='in')

#combine and save
#p.o <- ggpubr::ggarrange(p1,p2,labels=c('A','B'),ncol=1,align='v')#,heights = c(1,1.5))
#ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/inter_intra_non_homogeneous_1.pdf',plot=p.o, width = 8, height = 5,unit='in')




# to report percentage of monophyletic communities in infomap based analysis
out.mono.checked.tbl %>%
  filter(method=='cluster_infomap') %>%
  mutate(simmilarity = as.numeric(gsub('c','',level))) %>%
  filter(simmilarity>=0.7) %>% group_by(mono) %>% tally()


#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#
#   Older codes down here, could be interesting in many ways....
#
#--------------------------------------------------------------------------------------------#

#ANI CLIQUES
#============================================================================================#
#---- if you want to drop the outgroup from phylogeny here
tree <- sco.94.concat.tree.rooted
tr.tree <- tree
if(drop.outgroup==TRUE){
  tr.tree <- ape::drop.tip(phy=tr.tree,tip=outgroup)
}

#============================================================================================#

#AAI = AAI.f/100
multi_clique.ANIb <- multiClique(distance=ANIb[sel.genome,sel.genome], cutoffs = c(0.7,0.75,0.80,0.85,0.90,0.95,0.98),no_singletons = TRUE)
multi_clique.AAI <- multiClique(distance=AAI[sel.genome,sel.genome], cutoffs = c(0.7,0.75,0.80,0.85,0.90,0.95,0.98),no_singletons = TRUE)

#subset, ignore 0.7 when ingroup
multi.clique.ANI <- multi_clique.ANIb[,c('genome','c0.75','c0.8', 'c0.85', 'c0.9', 'c0.95', 'c0.98')]
sorted.ANI.cliques <- sortCliques(multi.clique.ANI, tree = tr.tree)

multi.clique.AAI <- multi_clique.AAI[,c('genome','c0.75','c0.8', 'c0.85', 'c0.9', 'c0.95', 'c0.98')]
sorted.AAI.cliques <- sortCliques(multi.clique.AAI, tree = tr.tree)



#============================================================================================#

p.anib <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.ANI.cliques, circular=TRUE)
p.aai <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.AAI.cliques, circular=TRUE)

p.combined <- gridExtra::arrangeGrob(p.anib, p.aai, nrow=1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/ANI_AAI_sortedCliques_phylogeny.pdf',plot=p.combined, width = 30, height = 15,unit='in')


#==============================================================================================
# One could ask, is this possible using different measures of distance?
#test on cophenetic distance, actually, this is real distance, 0 is similar, while in ANI and AAI this is similarity

COPH[multi_clique.AAI$genome, multi_clique.AAI$genome]
multi_clique.COPH <- multiClique(distance=1-COPH[multi_clique.AAI$genome, multi_clique.AAI$genome], cutoffs = c(0.7,0.75,0.80,0.85,0.90,0.95,0.98),no_singletons = TRUE)

#subset, ignore 0.7 when ingroup
multi.clique.COPH <- multi_clique.COPH[,c('genome','c0.75','c0.8', 'c0.85', 'c0.9', 'c0.95', 'c0.98')]
sorted.COPH.cliques <- sortCliques(multi.clique.COPH, tree=tr.tree)

p.coph <- plot.cliqueSort(tree=tr.tree, sorted.cliques = sorted.COPH.cliques, circular=TRUE)
p.combined <- gridExtra::arrangeGrob(p.coph, p.anib, p.aai, nrow=1)
ggplot2::ggsave(filename='~/DATA/MarbGenomics/Graphs/COPH_ANI_AAI_sortedCliques_phylogeny.pdf',plot=p.combined, width = 45, height = 15,unit='in')


#----
#============================================================================================#
# SET ANI CLIQUE COLORS
#============================================================================================#
outdir = '~/DATA/MarbGenomics/Graphs/'

sorted.ANI.cliques
sorted.cliques <- sorted.ANI.cliques
colnms <- colnames(sorted.cliques)[grepl('sc',colnames(sorted.cliques))]
base.column <- 'sc0.8'
#columns <- c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')

columns <- c('cc0.8','cc0.85','cc0.9','cc0.95','cc0.98')
base.column <- 'cc0.8'


sorted.cliques <- sorted.cliques %>%
  left_join(genome.tbl,by='genome') %>%
  select(genome, phylogroup,c0.75:cc0.98)
columns <- c('phylogroup','cc0.8','cc0.85','cc0.9','cc0.95','cc0.98')
base.column <- 'phylogroup'

sorted.cliques <- sorted.cliques %>%
  mutate(cphylo=as.numeric(phylogroup)) %>%
  mutate(scphylo=as.numeric(phylogroup))%>%
  mutate(ccphylo=as.numeric(phylogroup))

#colors = c("#39811D","#76AECF", "#A6D8D4", "#F2A968", "#F2EC70", "#E38FDD", "#898989", "#86DA81", "#83731B", "#B34D22")

### IMPLEMENT EXTRACTION OF COLORS FROM THE BLOODY
#out.colors <- blend_cluster_colors(sorted.cliques, columns, base.column, colors )

#out.colors<-hsv_colors(sorted.cliques,c(0.8,0.85,0.9,0.95,0.98),colors)
#out.colors<-hsv_colors(sorted.cliques,c(0.8,0.85,0.9,0.95,0.98),colors)

out.colors<-hsv_colors(sorted.cliques,c('phylo',0.8,0.85,0.9,0.95,0.98),phylogroup.colors)


#====





#single fixed one, but why bother, lets go multi
colcol <- 'scphylo'
columns = c('scphylo','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')
#Annotate trees
p.list <- list()
for(colcol in columns){
  p <- ggtree(tr.tree)
  p$data <- p$data %>% left_join(sorted.cliques,by=c('label'='genome')) %>%
    left_join(out.colors, by=c('label'='genome'))  %>% data.frame()

  #p$data$group <- p$data[,paste0("col.",colcol)]
  p$data$group <- paste0('cl',p$data[,colcol])
  p$data$group <- factor(p$data$group,levels = c(paste0('cl',1:(length(unique(p$data$group))-1)),'clNA'))
  p$data$group[p$data$group=='clNA']=NA

  colors <- p$data %>% filter(isTip==TRUE) %>% select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>% unique() %>% data.frame()
  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,'group'])
  grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]

  #grid.colors = factor(grid.colors,levels = c(paste0('cl',1:(length(unique(plotdat$group))-1)),'clNA'))
  p.tr <- p + geom_tippoint(data=p$data[!is.na(p$data$group),], aes(fill=group),shape=21,size=2) +scale_fill_manual(values=grid.colors)
  #ggplot2::ggsave(filename=paste0(outdir,'tree_ANIb_cliques',colcol,'.pdf'),plot=p.tr, width = 3, height = 6,unit='in')
  p.list[[colcol]]<-p.tr
}

p.comb.all <- gridExtra::grid.arrange(grobs = p.list, nrow=1)
p.comb.all
ggplot2::ggsave(filename=paste0(outdir,'tree_ANIb_cliques_all','.pdf'),plot=p.comb.all, width = 13, height = 6,unit='in')



#genome characteristics barchart
p.list <- list()
for(colcol in c('scphylo','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')){
  plotdat <-genome.tbl %>%
    left_join(sorted.cliques,by=c('genome'='genome')) %>%
    left_join(out.colors, by=c('genome'='genome')) %>%
    filter(genome %in% sel.genome) %>%
    data.frame()
  plotdat$group <- paste0('cl',plotdat[,colcol])
  plotdat$group <- factor(plotdat$group,levels = c(paste0('cl',1:(length(unique(plotdat$group))-1)),'clNA'))

  colors <- plotdat%>% select(paste0('col.',gsub('sc','cc',colcol)),colcol, group) %>%
    filter(group!='clNA') %>%
    unique() %>% data.frame()

  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,'group'])
  grid.colors <- grid.colors[names(grid.colors)[!is.na(names(grid.colors))]]
  grid.colors <- c(grid.colors,'clNA'='white')

  p.1 <- plotdat %>% ggplot(aes(x=group,y=Genome_size/1000000))+
    geom_hline(yintercept = mean(plotdat$Genome_size/1000000), size = 0.25, colour = '#bdbdbd') +
    geom_boxplot(aes(fill=group),alpha=0.3,outlier.size=-1)+
    geom_jitter(aes(fill=group),height=0,width=0.2,size=2,shape=21,alpha=0.7) +
    ylab('Genome size (Mbp)') +
    xlab('ANI cliques') +
    scale_fill_manual(values=grid.colors) + fdb_style(aspect.ratio=0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

  p.2 <- plotdat %>% ggplot(aes(x=group,y=GC))+
    geom_hline(yintercept = mean(plotdat$GC), size = 0.25, colour = '#bdbdbd') +
    geom_boxplot(aes(fill=group),alpha=0.3,outlier.size=-1)+
    geom_jitter(aes(fill=group),height=0,width=0.2,size=2,shape=21,alpha=0.7) +
    ylab('GC-content (%)') +
    xlab('ANI cliques') +
    scale_fill_manual(values=grid.colors) + fdb_style(aspect.ratio=0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))

  p.3 <- plotdat %>% ggplot(aes(x=group,y=Coding_density))+
    geom_hline(yintercept = mean(plotdat$Coding_density), size = 0.25, colour = '#bdbdbd') +
    geom_boxplot(aes(fill=group),alpha=0.3,outlier.size=-1)+
    geom_jitter(aes(fill=group),height=0,width=0.2,size=2,shape=21,alpha=0.7) +
    ylab('Coding density (%)') +
    xlab('ANI cliques') +
    scale_fill_manual(values=grid.colors) + fdb_style(aspect.ratio=0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))


  p.comb <- gridExtra::grid.arrange(p.1, p.2, p.3, ncol=1)
  ggplot2::ggsave(filename=paste0(outdir,'Genome_characteristics_ANIb_cliques',colcol,'.pdf'),plot=p.comb, width = 3, height = 6,unit='in')
  p.list[[colcol]]<-p.comb
}

p.comb.all <- gridExtra::grid.arrange(grobs = p.list, nrow=1)
ggplot2::ggsave(filename=paste0(outdir,'Genome_characteristics_ANIb_cliques_all','.pdf'),plot=p.comb.all, width = 13, height = 5,unit='in')



#------- BELOW IS IMPORTANT!!!!
# just some shuffeling into the genome table, and aims to start to think how to merge color scheme, and new grouping data


metadata <- genome.tbl %>% left_join(sorted.ANI.cliques, by='genome') %>%
  mutate(group=paste0('cl',as.character(c0.85)))

# CONSIDER SETTING IT AT DIFFERENT CUTOFFS

mtd <- genome.tbl %>% left_join(sorted.ANI.cliques, by='genome') %>%
  #filter(!is.na(c0.8)) %>%
  #filter(quality_class!='low') %>%
  mutate(group=paste0('cl',as.character(c0.85)))


#------------------------



#==
#sorted_cl <- sorted.ANI.cliques
#sel.col <- c('sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')

multi.cl <- multiClique(distance=sim.list[[dst.name]], cutoffs = cutoff_range,no_singletons = TRUE, community.method = com.meth )
sorted.cl <- sortCliques(multi.cl, tree=tr.tree)
sorted_cl <- sorted.cl

sel.col <- sorted_cl %>% select(sc0.65:sc0.995) %>% colnames()
output <- chisq_range(sorted_cl, sel.col)





p.chi <- output %>% ggplot(aes(x=as.factor(depth),y=-log10(p.value))) +
  geom_bar(stat = 'identity',aes(fill=df,alpha=method), width = 1)+
  ylab('-log10(p-value)')+
  xlab('ANI clique threshold')+
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05),color='grey60',size=0.5)+
  geom_hline(yintercept = -log10(0.01),color='grey60',size=0.5,linetype='dashed')+
  geom_hline(yintercept = -log10(0.001),color='grey60',size=0.5,linetype='dotted')+
  scale_x_discrete(breaks=c('sc0.7','0.75','0.8','0.85','0.9','0.95'))+
  fdb_style(aspect.ratio = .5)

ggplot2::ggsave(filename=paste0(outdir,'chisq_NNIcliques_habitat_pvalue','.pdf'),plot=p.chi, width = 8, height = 3,unit='in')

p.chi <- output %>% ggplot(aes(x=as.factor(depth),y=Xsquared)) +
  geom_bar(stat = 'identity',aes(fill=df,alpha=method), width = 1)+
  ylab('Xsq')+
  xlab('ANI clique threshold')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks=c('sc0.7','0.75','0.8','0.85','0.9','0.95'))+
  fdb_style(aspect.ratio = .5)

ggplot2::ggsave(filename=paste0(outdir,'chisq_NNIcliques_habitat_Chisq','.pdf'),plot=p.chi, width = 8, height = 3,unit='in')






#output %>% ggplot(aes(x=depth,y=-log10(p.value)))+geom_point()+geom_line()+fdb_style()
#output %>% ggplot(aes(x=depth,y=df))+geom_point()+geom_line()+fdb_style()
#output %>% ggplot(aes(x=depth,y=n_A))+geom_point()+geom_line()+fdb_style()
#output %>% ggplot(aes(x=depth,y=Xsquared))+geom_point()+geom_line()+fdb_style()
#output %>% ggplot(aes(x=df,y=n_A))+geom_point()+geom_line()+fdb_style()


#============================================================================================#
# Chord Diagrams
#============================================================================================#
library(circlize)

#range=c('sc0.75','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')
range=c('scphylo','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')
anglerot = 80
anglerot = 258

for(selcol in range){
  out.colors2 <-out.colors
  colnames(out.colors2) <- c('genome','col.scphylo','col.sc0.8','col.sc0.85','col.sc0.9','col.sc0.95','col.sc0.98')

  mtd <- genome.tbl %>% left_join(sorted.cliques, by='genome') %>%
    left_join(out.colors2, by='genome') %>%
    filter(!is.na(!!as.name(selcol)))  %>%
    filter(!is.na(SS)) %>%
    data.frame()

  mtd$grp <- mtd[,selcol]
  mtd <- mtd %>% mutate(group=paste0('cl',grp))
  chording <- mtd %>% filter(!is.na(!!as.name(selcol))) %>% mutate() %>% filter(!is.na(SS))

  n_genomes <- nrow(chording)

  chord.df = as.data.frame.matrix(table(as.character(chording$group),as.character(chording$SS)))

  minimal_chord.df = chord.df[,2:length(colnames(chord.df))]
  ordered_minimal_chord.df = minimal_chord.df

  mat = as.matrix(ordered_minimal_chord.df)


  df = data.frame(from =  rep(colnames(mat), each = nrow(mat)),
                  to = rep(rownames(mat), times = ncol(mat)),
                  value = as.vector(mat),
                  stringsAsFactors = FALSE)
  df <- df %>% filter(!is.na(to))

  colors <- mtd %>% filter(!is.na(!!as.name(selcol))) %>%  filter(!is.na(SS)) %>% select(paste0('col.',selcol), group) %>% unique()
  colors <-colors[order(colors$group),]
  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,2])

  ordered_minimal_chord.df <- ordered_minimal_chord.df[ paste0('cl',1:nrow(ordered_minimal_chord.df)),]

  circos.clear()
  circos.par(start.degree = anglerot)
  circos.par(gap.after = c(rep(5, length(unique(df[[2]]))-1), 20,
                           rep(5, length(unique(df[[1]]))-1), 20))
  habitat.colors = c(
    'Hydrothermal_Vent'="#E41A1C", 'Lake'="#F37912", 'Oil'="#FFD422", 'Photic'="#43997A", 'Phototroph'="#658E67", 'Polar'="#5D6795", 'Sediment'="#A35390", 'Terrestial'="#B85F49",'V1'='black'
  )

  grid.colors <- c(grid.colors, habitat.colors)
  # <- c(habitat.colors, grid.colors)

  pdf(paste0('~/DATA/MarbGenomics/Graphs/habitat_ANI_',selcol,'_',anglerot,'.pdf'), width=3, height=3)
  chordDiagram(as.matrix(ordered_minimal_chord.df), grid.col = grid.colors,
               annotationTrack = c("grid"), annotationTrackHeight = uh(4, "mm"),reduce=0)
  for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim),cex = 0.4, si, sector.index = si, track.index = 1,
                facing = "bending.inside", niceFacing = TRUE, col = "white")
  }
  title(main=paste0('ANI ',selcol,' (n=',n_genomes,')'))
  dev.off()
}

#highlight.sector(rownames(ordered_minimal_chord.df), track.index = 1, col = "grey",
#                 text = "ANI", cex = 0.8, text.col = "white", niceFacing = TRUE)
#highlight.sector(colnames(ordered_minimal_chord.df)[-2], track.index = 1, col = "darkgrey",
#                 text = "HABITAT", cex = 0.8, text.col = "white", niceFacing = TRUE)




#--------------------
# Alternative label positions

#range=c('sc0.75','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')
range=c('scphylo','sc0.8','sc0.85','sc0.9','sc0.95','sc0.98')
anglerot = 80
anglerot = 258
#anglerot = 270
for(selcol in range){
  out.colors2 <-out.colors
  colnames(out.colors2) <- c('genome','col.scphylo','col.sc0.8','col.sc0.85','col.sc0.9','col.sc0.95','col.sc0.98')

  mtd <- genome.tbl %>% left_join(sorted.cliques, by='genome') %>%
    left_join(out.colors2, by='genome') %>%
    filter(!is.na(!!as.name(selcol)))  %>%
    filter(!is.na(SS)) %>%
    data.frame()

  mtd$grp <- mtd[,selcol]
  mtd <- mtd %>% mutate(group=paste0('cl',grp))
  chording <- mtd %>% filter(!is.na(!!as.name(selcol))) %>% mutate() %>% filter(!is.na(SS))

  n_genomes <- nrow(chording)

  chord.df = as.data.frame.matrix(table(as.character(chording$group),as.character(chording$SS)))

  minimal_chord.df = chord.df[,2:length(colnames(chord.df))]
  ordered_minimal_chord.df = minimal_chord.df

  mat = as.matrix(ordered_minimal_chord.df)


  df = data.frame(from =  rep(colnames(mat), each = nrow(mat)),
                  to = rep(rownames(mat), times = ncol(mat)),
                  value = as.vector(mat),
                  stringsAsFactors = FALSE)
  df <- df %>% filter(!is.na(to))

  colors <- mtd %>% filter(!is.na(!!as.name(selcol))) %>%  filter(!is.na(SS)) %>% select(paste0('col.',selcol), group) %>% unique()
  colors <-colors[order(colors$group),]
  grid.colors <- as.character(colors[,1])
  names(grid.colors) <- as.character(colors[,2])

  ordered_minimal_chord.df <- ordered_minimal_chord.df[ paste0('cl',1:nrow(ordered_minimal_chord.df)),]

  circos.clear()
  circos.par(start.degree = anglerot)
  circos.par(circle.margin=0.25,gap.after = c(rep(5, length(unique(df[[2]]))-1), 20,
                           rep(5, length(unique(df[[1]]))-1), 20))
  #habitat.colors = c(
  #  'Hydrothermal_Vent'="#E41A1C", 'Lake'="#F37912", 'Oil'="#FFD422", 'Photic'="#43997A", 'Phototroph'="#658E67", 'Polar'="#5D6795", 'Sediment'="#A35390", 'Terrestial'="#B85F49",'V1'='black'
  #)

  grid.colors <- c(grid.colors, habitat.colors)
  # <- c(habitat.colors, grid.colors)

  pdf(paste0('~/DATA/MarbGenomics/Graphs/habitat_ANI_',selcol,'_',anglerot,'.pdf'), width=3, height=3)
  chordDiagram(as.matrix(ordered_minimal_chord.df), grid.col = grid.colors,
               order= c(rev(paste0('cl',1:length(rownames(mat)))),sort(colnames(mat))),
               annotationTrack = c("grid"),
               big.gap = 20,
               small.gap = 2,
               annotationTrackHeight = uh(2, "mm"),reduce=0)
  for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim),
                cex = 0.6,
                si,
                sector.index = si,
                track.index = 1,
                facing = "clockwise",
                niceFacing = TRUE,
                col = "black",
                adj=c(-0.7,0.5)#position labels outside
                )
  }
  title(main=paste0('ANI ',selcol,' (n=',n_genomes,')'))
  dev.off()
}

#highlight.sector(rownames(ordered_minimal_chord.df), track.index = 1, col = "grey",
#                 text = "ANI", cex = 0.8, text.col = "white", niceFacing = TRUE)
#highlight.sector(colnames(ordered_minimal_chord.df)[-2], track.index = 1, col = "darkgrey",
#                 text = "HABITAT", cex = 0.8, text.col = "white", niceFacing = TRUE)








g <- graph_from_adjacency_matrix(as.matrix(ANIb), mode = "upper", weighted = T, diag = F)
wc <- cluster_infomap(g, e.weights = E(g)$weight)
cliques(g)
modularity(wc)
membership(wc)


g2 <- delete.edges(g, which(E(g)$weight <0.95))
plot(g2)
wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)


cliques = membership(wc)



g <- graph_from_adjacency_matrix(as.matrix(AAI), mode = "upper", weighted = T, diag = F)
g2 <- delete.edges(g, which(E(g)$weight < 87))
plot(g2)
wc <- cluster_walktrap(g2)
modularity(wc)
membership(wc)
plot(wc, g2)
cliques = membership(wc)



p <- ggtree(tree)
large_cliques = cliques[cliques %in% names(table(cliques)[table(cliques)>1])]

p$data$group <- large_cliques[p$data$label]
p+geom_point(data=p$data[p$data$isTip ==TRUE,],aes(x,y,color=as.factor(group)))

#--------------------------------------------



cliques[cliques==12]

sortedCliqueSizes = sort(table(cliques),decreasing=TRUE)

ggplot(as.data.frame(sortedCliqueSizes), aes(x=cliques, y = Freq)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_bar(stat="identity",fill='#999999',color='white')+
  scale_color_manual(values='lightgrey')+
  ylab("Genomes")+
  xlab("ANI (95%) clique")+
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.position = "none")

ggplot(as.data.frame(sortedCliqueSizes[sortedCliqueSizes>1]), aes(x=cliques, y = Freq)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_bar(stat="identity",fill='#999999',color='white')+
  scale_color_manual(values='lightgrey')+
  ylab("Genomes")+
  xlab("ANI (95%) clique")+
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.position = "none")



#---------------------------------------

#none singleton cliques

non_signleton_cliques = sortedCliqueSizes[sortedCliqueSizes>1]
singleton_cliques = sortedCliqueSizes[sortedCliqueSizes==1]

length(sortedCliqueSizes)
length(non_signleton_cliques)
length(singleton_cliques)


#---------------------------------------

for(clq in names(non_signleton_cliques)){
  print('--------------------------')
  print(paste('---', clq))
  print(names(cliques[cliques == clq]))
  print('---')
  print(length(names(cliques[cliques == clq])))

  print(mergedEnv[names(cliques[cliques == clq]),c('TypeStrain', 'geo_loc_name','SS')])
  print('is monophyletic?')
  print(is.monophyletic(tree, tips= names(cliques[cliques == clq])))
}





#------------------------------------------------




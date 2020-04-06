library(qgraph)
library(igraph)
library(mclust)
library(clues)
library(shiny)
library(visNetwork)

#' Concordance Network
#'
#' This function receives a vector and returns the corresponding NxN concordance matrix.
#' @param vec vector.
#' @keywords vec to matrix
#' @return the corresponding NxN concordance matrix.
#' @examples
#'
#' vec <- sample(c(1,0), 4, replace=TRUE)
#' vec2concord(vec)

vec2concord <- function(vec){
  vec %*% t(vec)
}

#' Concordance Network
#'
#' This function receives an NXN matrix and returns the corresponding list of NxN concordance matrices.
#' @param mat NxN matrix.
#' @keywords vecs to matrices
#' @return the corresponding list of NxN concordance matrices.
#' @examples
#'
#' mat <- matrix(sample(c(1,0), 16, replace=TRUE), 3, 3)
#' vecs2concords(mat)

vecs2concords <- function(mat) {
  # aggregate cocncordance matrices in a list
  lapply(1:nrow(mat), function(i) vec2concord(mat[i,]))
}

#' Similarity Matrix
#'
#' This function receives a matrix where each row is a binary vector (e.g. symptoms' presence vector)
#' and returns the corresponding NxN similarity matrix (pair-wise similarities) using the Adjusted Hubert-Arabie RAND index (ARI).
#' @param mat matrix
#' @param concords list of concordance matrices (defaults to vecs2concords(mat)) to be compared.
#' @param simil_measure string value, similarity function used, can be either 'ARI' (default) or 'Euclidean'.
#' @param reactive_dom reactive domain in case used in a shiny app (will be used by incProgress()).
#' @param inProg_amount real number in [0,1], the amount of status bar to fill up (default is 1), in case used in a shiny app (will be used by incProgress()).
#' @keywords similarity matrix pair-wise comparisons
#' @return the corresponding NxN similarity matrix (pair-wise similarities) using the Adjusted Hubert-Arabie RAND index (ARI).
#' @examples 
#' mat <- matrix(sample(c(1,0), 16, replace=TRUE), 3, 3)
#' similarity_matrix(mat)

similarity_matrix <- function(mat, concords=vecs2concords(mat), simil_measure='ARI', reactive_dom=NULL, incProg_amount=1) {

  # validate incProg_amount
  if (incProg_amount <= 0 || incProg_amount > 1){ incProg_amount <- 1 }

  # Compute similarities
  sims <- NULL

  if (simil_measure == 'Euclidean') {
    a.conc <- t(apply(mat, 1, function(x) (x%*%t(x))[upper.tri(x%*%t(x), diag=FALSE)])) # taking just upper triangle w/o diagonal
    dists <- as.matrix(dist(a.conc))
    sims <- sqrt(1/(1+dists))
    if (!is.null(reactive_dom)){
      incProgress(0.5*incProg_amount, session=reactive_dom)
    }
  } else {
    sims <- matrix(0, nrow=length(concords), ncol=length(concords))

    # vectorize concords
    vConcords <- lapply(1:length(concords), function(i) as.vector(concords[[i]][upper.tri(concords[[i]], diag = F)]))

    # compute similarity matrix by calculting pairwise ARI
    # (upper triangular part - not including diagonal)
    for (i in 1:(nrow(sims)-1)) {
      for (j in (i+1):ncol(sims)) {
        sims[i,j] <- adjustedRand(vConcords[[i]], vConcords[[j]], randMethod='HA')
      }
      if (!is.null(reactive_dom)){
        incProgress((1/nrow(sims))*incProg_amount, session=reactive_dom)
      }
    }

    # diagonal (all 1s)
    for (i in 1:nrow(mat)) {
      sims[i,i] <- 1
    }

    # reflect upper triangle accross diagonal (cause symmetric)
    sims = sims + t(sims)
  }

  # set negative entries to 0
  sims[sims < 0] <- 0

  sims
}

#' Normalized Community Matrix
#'
#' This function receives a list of concordnce NxN matrices (a patient community)
#' and returns the corresponding normalized patient community concordance matrix.
#'
#' It does so by summing the concordance matrices of the given community,
#' and normalizing offdiagonal entries by dividing each entry.
#' by the sqrt of the product of the associated diagonal entries,
#' (giving it a value in [0,1]).
#'
#' @param mats list of NxN matrices.
#' @keywords list, matrices
#' @return the corresponding normalized patient community concordance matrix.
#' @examples
#'
#' normalized_community(mats)

normalize_community <- function(mats){

  # sum all the matrices in mats
  res <- Reduce(`+`, mats)

  # normalize entries
  for (i in 1:nrow(res)){
    for (j in 1:ncol(res)){
      # dividing off diagonal entries by sqrt of product of associated diagonal entries as the normalization mehtod (for now)
      # (consider manipulating normalizaiton method)
      # but if the entry is 0 than it's left as 0 as there's no point in dividing it by anything
      if (i != j) {
        res[i,j] <- ifelse(res[i,j] == 0, 0,  res[i,j]/sqrt(res[i,i]*res[j,j]))
      }
    }
  }

  res
}

#' Concordance Network Clustering
#'
#' This function receives an MxN matrix with named columns (will be used to name vertices),
#' and returns a list of igrpah graphs representing the resulting clusters (e.g. symptom clusters of "patient communities").
#'
#' It does so using the method of concordance networks clustering as described in Henry et al. 
#' "Concordance networks and application to clustering cancer symptomology - PLOS." 14 Mar. 2018. 
#'
#' The detection algorithm (from the igraph library) can be chosen out of the following list:
#'
#' 'WT' (default) for cluster_walktrap.	Community strucure via short random walks.
#' 'FG' for cluster_fast_greedy(),	Community structure via greedy optimization of modularity.
#' 'IM' for cluster_infomap(),	Infomap community finding.
#' 'LP' for cluster_label_prop(),	Finding communities based on propagating labels.
#' 'LE' for cluster_leading_eigen(),	Community structure detecting based on the leading eigenvector of the community matrix.
#' 'LV' for cluster_louvain(),	Finding community structure by multi-level optimization of modularity.
#'
#' @param mat MxN matrix
#' @param thresh discriminatory threshold to determine present vs. not present
#' @param nrows number of rows to analyze (defaults to all rows)
#' @param ncols number of columns to analyze (defaults to all columns)
#' @param thresh integer threshld below which entries (e.g. symptom scores) will be set to 0 (i.e. considered as non-present)
#' @param detectAlgo string value, the type of network community detection algorithm to use (defaults to 'WT')
#' @param simil_measure string value, similarity function used (string value), can be either 'ARI' (default) or 'Euclidean'.
#' @param centrality string value, centrality measure to use for indicating central nodes (via node size in the igraph networks)
#' @param cluster_colors string value, 
#' @param sparsify number in [0,100] percentage of the weakest network edges to remove from the centrality computation, to allow central nodes to be detected more easily.
#' @param simplify_graphs
#' @param reactive_dom reactive domain in case used in a shiny app (will be used by incProgress()).
#' @param inProg_amount real number in [0,1], the amount of status bar to fill up (default is 1), in case used in a shiny app (will be used by incProgress()).
#' @keywords matrix list graphs communities
#' @return list consisting of (1) a list of the normalized community networks, and (2) the associated communities (e.g. as generated by igraph::groups(cluster_walktrap(g)))
#' @export
#' @examples
#' data = rbinom(300*20, 1, .2)
#' data = matrix(data,300,20)
#' results = concord_cluster(data)

cluster_concordance <- function(mat, nrows=nrow(mat), ncols=ncol(mat), thresh=1, detectAlgo='WT', simil_measure='ARI',
                                centrality='Betweenness', cluster_colors=NULL, sparsify=0, simplify_graphs=TRUE,
                                reactive_dom=NULL, incProg_amount=1) {
    # validate incProg_amount
    if (incProg_amount <= 0 || incProg_amount > 1){ incProg_amount <- 1 }

    # if matrix is empty throw error
    if (nrow(mat) == 0){
      stop(safeError('Matrix was empty.'))
    }

    # ensure reproducibility
    set.seed(42)

    # fix row names
    rownames(mat) <- NULL

    # transform to present/not present (binary matrix)
    mat[mat < thresh] <- 0
    mat[mat >= thresh] <- 1

    # compute between patient similarity (generate similarity matrix)
    concords <- vecs2concords(mat)
    sims <- similarity_matrix(mat, concords, simil_measure=simil_measure, reactive_dom=reactive_dom, incProg_amount)

    # convert sims to a graph to prepare for random walk
    g <- graph_from_adjacency_matrix(sims, mode='undirected', weighted=TRUE, diag=FALSE)

    # remove multi-edges and loops
    if (simplify_graphs) {
      g <- simplify(g)
    }

    # increment progress
    if (!is.null(reactive_dom)){
      incProgress(0.2*incProg_amount, session=reactive_dom)
    }

    # run community detection algorithm to get communities (the default is Random Walk (walktrap))
    comms <- NULL

    if (detectAlgo == 'FG'){
      comms <- cluster_fast_greedy(g)
    } else if (detectAlgo == 'IM'){
      comms <- cluster_infomap(g)
    } else if (detectAlgo == 'LP'){
      comms <- cluster_label_prop(g)
    } else if (detectAlgo == 'LE'){
      comms <- cluster_leading_eigen(g)
    } else if (detectAlgo == 'LV'){
      comms <- cluster_louvain(g)
    } else {
      comms <- cluster_walktrap(g)
    }

    # convert comms to regular r list
    comms <- lapply(igraph::groups(comms), function(x) as.numeric(x))

    # create normalized community matrices
    # (based on now gained community information)
    norm_comms <- lapply(1:length(comms), FUN=function(i){
      # for each community
      # normalize and add community matrix to list of normalized communities
      normalize_community(concords[comms[[i]]])
    })

    # "inter-community" occurrence info for tooltip
    # (count number of occurrences of symptom greater than or equal to input$thresh in current category among all communities)
    inter_comm_occs <- lapply(1:ncol(norm_comms[[1]]), function(i) sum(vapply(norm_comms, FUN.VALUE=0, FUN=function(nc) diag(nc)[[i]])))

    # ADD EACH NORMALIZED COMMUNITY (network graph) TO GRAPHS LIST
    nets <- lapply(1:length(norm_comms), FUN=function(i) {

      # store matrix locally
      m <- norm_comms[[i]]
      diag(m) <- 0

      # Sparsify matrix to get more meaningful results
      sparsify_pctg <- sparsify/100
      m[m <= quantile(m, probs=sparsify_pctg)] <- 0

      if (simplify_graphs){
        # store grapah with original co-occurrences so can display that as edge thickness
        g_orig_co_occ <- simplify(graph_from_adjacency_matrix(m, mode='undirected', weighted=TRUE, diag=FALSE))
      }

      # Crucial data transformaiton to properly compute betweenness
      m[m != 0] <- 1/m[m != 0]
      if (max(m) != 0) {
        m <- m/max(m)
      }

      # skip network if no edges exist (return empty graph)
      if (all(m == 0)){
        comms[[i]] <<- NULL # get rid of community
        return(graph.empty())
      }

      # convert community to a graph
      g <- graph_from_adjacency_matrix(m, mode='undirected', weighted=TRUE, diag=FALSE)

      # simplify graph (remove multi-edges and loops)
      if(simplify_graphs){
        g <- simplify(g)
      }

      # vertices are implicitly labeled by symptom name (embedded in the matrix object)
      # (labels are stored in V(g)$label)

      # Compute centrality of nodes based on chosen centrality measure
      # and set size of each vertex by chosen centrality
      centralities <- NULL
      if (centrality == 'Strength'){
        centralities <- (strength(g_orig_co_occ, loops=FALSE)/max(strength(g_orig_co_occ, loops=FALSE)))
        V(g)$size <- 20 + 15*(centralities^2)
      } else if (centrality == 'Closeness'){
        centralities <- closeness(g, normalized = TRUE, mode='all')
        V(g)$size <- 20+50*(centralities^0.5)
      } else {
        centralities <- betweenness(g, normalized = TRUE, directed = FALSE)
        V(g)$size <- 20+40*(centralities^0.3)
      }

      # Add node centralities to nodes
      V(g)$centralities <- round(centralities, digits=3)

      # detect symptom clusters in graph and add to graph
      V(g)$cluster <- membership(cluster_walktrap(g))

      # color the nodes
      if (is.null(cluster_colors) || cluster_colors == 'DO NOT CLUSTER SYMPTOMS') {
        # based solely on beteenness
        V(g)$color <- rgb(.2,.2, (V(g)$size / max(V(g)$size))^1.5)
      } else {
        # differnt cluster colors
        color_pal <- hcl.colors(max(V(g)$cluster), palette = cluster_colors)
        V(g)$color <- lapply(V(g)$cluster, function(i) colorspace::darken(color_pal[[i]], amount = (1 - V(g)$size[[i]]/max(V(g)$size))))
      }

      # Add title  (neede for tooltip for vertices if a visNetwork is to be generated)
      V(g)$title <- paste('<b>', V(g)$name, '</b>',
                          '<ul>',
                          '<li>', 'Community:', i, '</li>',
                          '<li>', 'Centrality =', V(g)$centralities, '</li>',
                          '<li>', 'Community Occurrence =', diag(norm_comms[[i]]), '</li>',
                          '<li>', 'Inter-Comm. Occurrence =', inter_comm_occs, '</li>',
                          '</ul>')


      # set color & opacity of edges by their weight
      edge_color <- ifelse(is.null(cluster_colors) || cluster_colors == 'DO NOT CLUSTER SYMPTOMS', "dark green", "dark gray")
      E(g)$color <- apply(as.matrix((E(g_orig_co_occ)$weight/max(E(g_orig_co_occ)$weight))), c(1,2),
                          FUN=function(x) adjustcolor(edge_color, alpha.f=x^4 + 0.2))

      # Add title  (tooltip to edges)
      E(g)$title <- round(E(g_orig_co_occ)$weight, digits=2)

      # set width of edges by weight (to indicate the community's co-occurrence of pairs of symptomps)
      E(g)$width <- 2 + 12 * (E(g_orig_co_occ)$weight - min(E(g_orig_co_occ)$weight))/(max(E(g_orig_co_occ)$weight) - min(E(g_orig_co_occ)$weight))

      return(g)
    })

    # remove empty networks from 'nets'
    nets <- Filter(function(n) length(V(n))!=0, nets)

    # increment progress
    if (!is.null(reactive_dom)){
      incProgress(0.2*incProg_amount, session=reactive_dom)
    }

    return(list(nets, comms))
}

#' visNetwork Rendering
#'
#' This function receives a list of igrpah network objects, and optionally a title and a pre-title
#' and returns a visNetwork representing the given networks (using the visNetwork package)
#'
#' @param nets a list of igrpah network objects
#' @param title main title
#' @param pre_title title before the main title, as in "<pre-titile>: <title>"
#' @keywords list graphs networks visNetwork
#' @return a visNetwork representing the given networks
#' @examples
#' 
#' data = rbinom(300*20, 1, .2)
#' data = matrix(data,300,20)
#' res <- cluster_concordance(mat)
#' nets <- res[[1]]
#' 
#' getVisNet(nets, title, pre_title)

getVisNet <- function(nets, title=NULL, pre_title=NULL) {
  # initialize nodes df and edges df to be populated with all networks
  nodes <- data.frame()
  edges <- data.frame()

  # remove isolated vertices
  nets <- lapply(nets, FUN=function(net) delete.vertices(net, igraph::degree(net)==0))

  # remove empty networks from 'nets'
  nets <- Filter(function(n) length(V(n))!=0, nets)

  # concatenate the different netwroks/communities
  for (i in 1:length(nets)) {
    V(nets[[i]])$name <- paste0(V(nets[[i]])$name, i) # change names of nodes to be added
    node_df <- igraph::get.data.frame(nets[[i]], what='vertices') # get data frame of nodes from graph
    node_df['comm'] <- i # add node community number
    edge_df <- igraph::get.data.frame(nets[[i]], what='edges')    # get data frame of edges from graph
    edge_df['comm'] <- i # add edge community number
    nodes <- bind_rows(nodes, node_df)                # add nodes
    edges <- bind_rows(edges, edge_df)                # add nodes
  }

  # remove index from node labels
  nodes$label<- substr(nodes$name, 1, nchar(nodes$name)-1)

  # add necessary columns
  nodes$id <- nodes$name
  nodes$font.size <- nodes$size

  # title for visNetwork
  title <- ifelse(is.null(title), "", title)
  visNet_title <- ifelse(is.null(pre_title), title, paste0(pre_title, ': ', title))

  # return visNetwork output
  visNetwork(nodes, edges, main=visNet_title) %>%
    visEdges(smooth = FALSE) %>%
    visNodes(font=list(strokeWidth=4, strokeColor='white')) %>%
    visPhysics(timestep=0.2, barnesHut=list(avoidOverlap=1, gravitationalConstant=-10000, springLength=200, springConstant=0.01, centralGravity=0.3)) %>%
    visInteraction(hover=TRUE, navigationButtons=TRUE) %>%
    visLayout(randomSeed = 42) %>%
    visExport(name = ifelse(is.null(title), 'Network', title), label = 'Save as PNG', float='left')
}

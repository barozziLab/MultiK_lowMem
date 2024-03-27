# Leiden clustering functions
# This workflow was inspired by the Seurat package (see comments in the code below, describing the origin of the different parts more specifically)

##################################################################
##                        Find Neighbors                        ##
##################################################################
find_neighbors <- function (
    data,
    metric = "euclidean",
    n.trees = 50,
    k = 20,
    include.distance = T,
    search.k = -1,
    prune.SNN = 1/15,
    index = NULL
) {
  
  query <- data
  k.param <- k
  #~~~~~~~~~~~~~~~~~~
  # Annoy Build Index
  ###################
  
  # (inspired by the AnnoyBuildIndex function from the seurat package in https://github.com/satijalab/seurat/blob/HEAD/R/clustering.R)
  f <- ncol(data)
  
  index <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  
  for (ii in seq(nrow(x = data))) {
    index$addItem(ii - 1, data[ii, ])
  }
  index$build(n.trees)
  
  #~~~~~~~~~~~~
  # AnnoySearch
  #############
  
  # (inspired by the AnnoySearch function from the seurat package in https://github.com/satijalab/seurat/blob/HEAD/R/clustering.R)
  
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  
  res <- lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    #print(i)
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  nn <- list(nn.idx = idx, nn.dists = dist)
  
  #~~~~~~~~~~~~ (end of adapted AnnoySearch function part)
  
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  
  # create a new neighbor object (which is the output of the NNHelper function in the seurat clustering workflow)
  nn.ranked <- new(
    Class = 'Neighbor',
    nn.idx = nn$nn.idx,
    nn.dist = nn$nn.dists,
    alg.info = nn$alg.info %||% list(),
    cell.names = rownames(x = query)
  )
  
  nn.ranked <- Indices(object = nn.ranked)
  
  # Now the workflow continues according to the FindNeighbors function
  # convert nn.ranked into a Graph
  j <- as.numeric(x = Matrix::t(x = nn.ranked))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(object = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = data), nrow(x = data))), Class = "Graph")
  rownames(x = nn.matrix) <- rownames(x = data)
  colnames(x = nn.matrix) <- rownames(x = data)
  neighbor.graphs <- list(nn = nn.matrix)
  
  # compute SNN
  # Here, we need the Seurat implementation of ComputeSNN which is cpp code that can be found here: https://github.com/satijalab/seurat/blob/ca4c48b6300d7c10857e2d59a8ddee2858c7e2fc/src/RcppExports.cpp
  
  message("Computing SNN")
  
  snn.matrix <- Seurat:::ComputeSNN(
    nn_ranked = nn.ranked,
    prune = prune.SNN
  )
  
  rownames(x = snn.matrix) <- rownames(x = data)
  colnames(x = snn.matrix) <- rownames(x = data)
  snn.matrix <- as.Graph(x = snn.matrix)
  neighbor.graphs[["snn"]] <- snn.matrix
  
  return(neighbor.graphs)
  
}

#################################################################
##                        Find Clusters                        ##
#################################################################
leiden_clustering <- function (
    neighbor.graph,
    method = "igraph",
    partition.type = c(
      'RBConfigurationVertexPartition',
      'ModularityVertexPartition',
      'RBERVertexPartition',
      'CPMVertexPartition',
      'MutableVertexPartition',
      'SignificanceVertexPartition',
      'SurpriseVertexPartition'
    ),
    initial.membership = NULL,
    node.sizes = NULL,
    resolution.parameter = 1,
    n.iter = 10,
    random.seed = 0
) {
  
  object <- neighbor.graph
  clustering.results <- data.frame(row.names = colnames(x = neighbor.graph))
  
  #~~~~~~~~~~
  # RunLeiden
  ###########
  
  # inspired from here: https://github.com/satijalab/seurat/blob/HEAD/R/clustering.R
  
  # by default this is going to be the "matrix" switch, but maybe it pays off to use igraph for large datasets
  switch(
    EXPR = method,
    "matrix" = {
      input <- as(object = object, Class = "matrix")
    },
    "igraph" = {
      input <- if (inherits(x = object, what = 'list')) {
        graph_from_adj_list(adjlist = object)
      } else if (inherits(x = object, what = c('dgCMatrix', 'matrix', 'Matrix'))) {
        if (inherits(x = object, what = 'Graph')) {
          object <- as.sparse(x = object)
        }
        graph_from_adjacency_matrix(adjmatrix = object, weighted = TRUE)
      } else if (inherits(x = object, what = 'igraph')) {
        object
      } else {
        stop(
          "Method for Leiden not found for class", class(x = object),
          call. = FALSE
        )
      }
    },
    stop("Method for Leiden must be either 'matrix' or igraph'")
  )
  
  partition <- leiden(
    object = input,
    partition_type = partition.type,
    initial_membership = initial.membership,
    weights = NULL,
    node_sizes = node.sizes,
    resolution_parameter = resolution.parameter,
    seed = random.seed,
    n_iterations = n.iter
  )
  
  # There is a "group singletons" step. Include if needed
  
  # create a table/matrix of cell-ids and cluster-ids
  clustering.results[, paste("leiden", method, "res", resolution.parameter, sep = "_")] <- factor(x = partition)
  
  return(clustering.results)
  
}

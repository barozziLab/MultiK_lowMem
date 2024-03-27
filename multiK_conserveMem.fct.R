library(foreach)
library(HistogramTools)
library(progress)
library(doSNOW)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MultiK parallel conserve mem - main function
##############################

#' MultiKParallel main algorithm
#'
#' MultiKParallel main algorithm: takes a preprocessed gene expression matrix as input. Then subsamples 80\% of the cells and applies standard Seurat pipeline on the subsampled data matrix 100 times over multiple resolution parameters.
#' Adapted by Stephan Gruener (Nov 2023) to take a matrix of dimensionality reduced or non-negative matrix factorized counts (cells in rows, factors in columns) instead of a seurat object as input. Custom functions are used for computing the NN/SNN graphs and for leiden clustering. Cells are simply subsampled, but no additional normalization or dimensionality reduction is performed for subsets. This adapted workflow is intended to be used for clustering of cells based on NMF factors, but can in principle used on any dimensionality reduced matrix.
#' @param seu A Seurat object with normalized count
#' @param resolution A vector Seurat resolution parameters. Default is from 0.05 to 2 with step size of 0.05
#' @param nPC Number of principal components to use in clustering
#' @param reps Number of subsampling runs. Integer value. Default is 100
#' @param pSample Proportion of cells to sample. Numerical value. Default is 0.8
#' @param seed Optional numerical value. This sets a random seed for generating reproducible results
#' @return A list with components: k is a vector of number of runs for each K. clusters is a list containing the clustering labels for each subsampling run at each resolution parameter. consensus is a list containing a consensus matrix for each K.
#' @export
#' @importFrom foreach "%dopar%"
MultiK_par_conserveMem <- function(input.matr, 
                                   resolution = seq(0.05, 2, 0.05), 
                                   nFactors = ncol(input.matr), # by default, all factors are used
                                   reps = 100, 
                                   pSample = 0.8, 
                                   seed = NULL, 
                                   numCores_first = NULL,
                                   numCores_second = NULL,
                                   numChunks = NULL,
                                   leiden_method = "igraph",
                                   subrow = NULL
) {
  # setting seed for reproducibility
  if (is.null(seed) == TRUE) {
    seed <- timeSeed <- as.numeric(Sys.time())
  }
  set.seed(seed)
  
  # setting up doParallel
  if (is.null(numCores_first)) {
    numCores_first <- parallel::detectCores() / 2
  }
  cl <- parallel::makeCluster(numCores_first)
  doSNOW::registerDoSNOW(cl)
  #doParallel::registerDoParallel(numCores_first)
  print(paste0("doSNOW: numCores = ", numCores_first, "."))
  
  # step 1: subsampling
  if (is.null(subrow)) {
    subrow <- list()
    for (i in 1: reps) {
      subrow[[i]] <- sample(x=nrow(input.matr), size=round(nrow(input.matr) * pSample), replace=FALSE)
    }
  }
  
  # step 2: loop over subsampling runs, with each run subsampling 80% of cells, reselect genes for clustering
  print("Starting subsampling and clustering...")
  
  # setting up progress bar
  #progress <- function(n) cat(sprintf("-- rep %d is complete\n", n))
  #opts <- list(progress=progress)
  pb <- progress::progress_bar$new(format = "[:bar] :current/:total reps (:percent) |  :elapsed | eta: :eta", total = reps)
  
  progress <- function(n) {
    pb$tick()
  }
  opts <- list(progress = progress)
  
  results <- foreach::foreach (i = 1:reps, .options.snow = opts, .export = c("find_neighbors", "leiden_clustering"), .packages = c("Seurat", "Matrix", "leiden", "rlang", "igraph")) %dopar% {
    pb$tick()
    
    #print(paste("Rep: ", i, sep=""))
    
    # subsample the columns (the cells) from the full matrix
    subX <- input.matr[subrow[[i]], ]
    
    # # normalizing the data
    # #subX <- NormalizeData(object = subX, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)
    # 
    # # Find HVG genes ~ 2000 genes
    # subX <- Seurat::FindVariableFeatures(object = subX, selection.method = "vst", 
    #                                      nfeatures = 2000, loess.span = 0.3, 
    #                                      clip.max = "auto", num.bin = 20, 
    #                                      binning.method = "equal_width", verbose = F)
    # 
    # # Scaling unwanted variation
    # all.genes <- rownames(x = subX)
    # subX <- Seurat::ScaleData(object = subX, features = all.genes, verbose = F)
    # # Run PCA to reduce dimensions
    # subX <- Seurat::RunPCA(object = subX, features = Seurat::VariableFeatures(object = subX), npcs = nPC, verbose=F)
    
    # Run Clustering
    neighbors <- find_neighbors(data = subX[, 1:nFactors])
    
    clusters <- list()
    messages <- c()
    ks <- c()
    
    for (res in resolution) {
      print(paste("Rep", i, "Res", res, sep=" "))
      partition <- leiden_clustering(neighbor.graph = neighbors[["snn"]], 
                                     resolution.parameter = res, 
                                     method = leiden_method)
      subX.clusters <- partition[, 1]
      names(subX.clusters) <- rownames(partition)
      message <- paste("Rep_", i, "_res_", res, sep = "")
      messages <- c(messages, message)
      clusters[[message]] <- subX.clusters
      ks <- c(ks, length(unique(subX.clusters)))
    }
    
    # send results
    results_temp <- list("messages" = messages, "clusters" = clusters, "ks" = ks)
    results_temp
    
  }
  parallel::stopCluster(cl)
  
  results_combined <- rep_combine(results)
  rm(results)
  gc()
  
  messages <- results_combined$messages
  clusters <- results_combined$clusters
  ks <- results_combined$ks
  
  print("DONE!")
  print("Starting to count co-occurences across subsampling runs...")
  
  # create a list (across unique ks) of cell x run matrices, NAs where cells were not sampled
  res_by_k <- list()
  
  # TODO: parallelize this for loop
  for (k in sort(unique(results_combined$ks))) {
    
    k_names <- results_combined$messages[results_combined$ks == k]
    m <- matrix(NA, nrow = nrow(input.matr), ncol = sum(results_combined$ks == k))
    rownames(m) <- rownames(input.matr)
    colnames(m) <- k_names
    
    m <- sapply(k_names, function(k_name) {
      m[names(results_combined$clusters[[k_name]]), k_name] <- results_combined$clusters[[k_name]]
      m[,k_name]
    })
    
    res_by_k[[paste(k)]] <- m
    
  }
  
  # run scoring
  count_vecs <- list()
  
  for (k in names(res_by_k)) {
    # only consider ks that were found more than once (otherwise the scoring doesn't work)
    if (dim(res_by_k[[k]])[2] < 2) {
      warning(paste0("Skipping comparisons for k = ", k, " since this result was identified in less than 2 runs! It will therefore be removed and not show up in the results."))
      next
    }
    
    print(paste0("Running comparisons for k = ", k))
    count_vecs[[k]] <- multik_scoring_wrapper(k_res = res_by_k[[k]], num_cores = numCores_second, n_chunks = numChunks)
  }
  
  # dignostic plots
  plt <- DiagMultiKPlot_par(ks = ks, res_by_k = res_by_k, count_vecs = count_vecs)
  
  # results
  combined_ret_list <- list("messages" = messages, "clusters" = clusters, "ks" = ks, 
                            "res_by_k" = res_by_k, "count_vecs" = count_vecs, "plots" = plt)
  
  print("ALL DONE!")
  
  return(combined_ret_list)
  
}

#~~~~~~~~~~~~~~~~~~~
# combine replicates
####################

rep_combine <- function(ret_list) {
  
  messages <- c()
  clusters <- list()
  ks <- c()
  for (i in 1:length(ret_list)) {
    rep_results <- ret_list[[i]]
    messages <- c(messages, rep_results$messages)
    clusters <- append(clusters, rep_results$clusters)
    ks <- c(ks, rep_results$ks)
    
  }
  
  combined_ret_list = list("messages" = messages, "clusters" = clusters, "ks" = ks)
  
  return(combined_ret_list)
  
}

#~~~~~~~~~~~~~~~~~~
#multik scoring par - inner function
###################
multik_scoring_par <- function(cellA_index, k_res_chunks, mode = "self", ind_kresA = NULL, ind_kresB = NULL) {
  
  # Make sure, k_res is a list and ind_kresA and ind_kresB are set
  if(is.null(ind_kresA) | !(class(k_res_chunks)[1] == "list")) {
    stop("Mode 'versus' requires a list of chunks as k_res input! The indices for the chunks to be compared have to be passed to ind_kresA and ind_kresB!")
    }
  
  count_vector <- rep(0, 1001) # initiate count_vector
  
  if (mode == "self") { # comparison between combinations of cells from the same chunk
    if (!is.null(ind_kresB)) warning("ind_kresB will not be used, because mode is set to 'self'!")
    
    for (cellB_index in (cellA_index + 1):nrow(k_res_chunks[[ind_kresA]])) { # loop through all "other" cells, not yet compared to cell A
      # same cluster assignments / co-occurences in the same sample, rounded to 3 digits after the comma, * 1000 - 1 to assign it to the index in the count vector
      result <- round(sum(k_res_chunks[[ind_kresA]][cellA_index, ] == k_res_chunks[[ind_kresA]][cellB_index, ], na.rm = TRUE) / sum(!is.na(k_res_chunks[[ind_kresA]][cellA_index, ]) ** !is.na(k_res_chunks[[ind_kresA]][cellB_index, ])), 3) * 1000 + 1
      count_vector[result] <- count_vector[result] + 1 # update count vector
    }
    return(count_vector)
    
  } else if (mode == "versus") { # comparisons between all combinations of cells between chunks
    
    for (cellB_index in 1:nrow(k_res_chunks[[ind_kresB]])) { # loop through all "other" cells
      # same cluster assignments / co-occurences in the same sample, rounded to 3 digits after the comma, * 1000 - 1 to assign it to the index in the count vector
      result <- round(sum(k_res_chunks[[ind_kresA]][cellA_index, ] == k_res_chunks[[ind_kresB]][cellB_index, ], na.rm = TRUE) / sum(!is.na(k_res_chunks[[ind_kresA]][cellA_index, ]) ** !is.na(k_res_chunks[[ind_kresB]][cellB_index, ])), 3) * 1000 + 1
      count_vector[result] <- count_vector[result] + 1 # update count vector
    }
    return(count_vector)
    
  } else {
    stop ("Value for mode not supported. Please set mode to 'self' or 'versus'!")
  }
}
  
#~~~~~~~~~~~~~~~
#split in chunks
################
# Produce the chunks of k_res for parallelization
split_chunks <- function(k_res, n_chunks) {
  
  k_res_chunks <- list()
  chunk_length <- floor(nrow(k_res) / n_chunks)
  first_row <- 1
  last_row <- chunk_length

  
  for (i in 1 : (n_chunks - 1)) {
    k_res_chunks[[i]] <- k_res[first_row:last_row, ]
    first_row <- last_row + 1
    last_row <- last_row + chunk_length
  }

  k_res_chunks[[n_chunks]] <- k_res[first_row:nrow(k_res), ] # in case of a remainder in nrow(kres)/n_chunks
  
  message(paste0("Created ", n_chunks, " chunks of size ", chunk_length, "."))
  
  return(k_res_chunks)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~
#multik scoring par helper
##########################

multik_score_chunks <- function(k_res_chunks, num_cores) {
  
  # set up comparisons
  comparisons <- list()
  count <- 1
  for (cA in 1:length(k_res_chunks)) {
    for (cB in cA:length(k_res_chunks)) {
      comparisons[[count]] <- c("cA" = cA, "cB" = cB)
      count = count + 1
    }
  }
  
  # prepare cluster
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)
  print(paste0("doSNOW: numCores = ", num_cores, "."))
  
  # progress <- function(n) cat(sprintf("-- comparison %d out of %d is complete\n", n, length(comparisons)))
  # opts <- list(progress=progress)
  pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent) |  :elapsed | eta: :eta", total = length(comparisons))
  
  progress <- function(n) {
    pb$tick()
  }
  opts <- list(progress = progress)
  
  # run comparisons
  res_matr <- foreach::foreach(comp = 1:length(comparisons), .combine = "cbind", .options.snow=opts, .export = c("multik_scoring_par")) %dopar% {
    pb$tick()
    
    cA <- comparisons[[comp]]["cA"]
    cB <- comparisons[[comp]]["cB"]
    
    if (cA == cB) {
      # SELF comparison (all vs. all others)
      count_matr <- sapply(1:(nrow(k_res_chunks[[cA]]) - 1),
                           FUN = multik_scoring_par,
                           k_res_chunks = k_res_chunks,
                           mode = "self",
                           ind_kresA = cA,
                           ind_kresB = cA,
                           USE.NAMES = FALSE,
                           simplify = "matrix")
      count_vec <- apply(count_matr, 1, sum)
    } else {
      # VERSUS comparison (all vs. all)
      count_matr <- sapply(1:(nrow(k_res_chunks[[cA]])),
                           FUN = multik_scoring_par,
                           k_res_chunks = k_res_chunks,
                           mode = "versus",
                           ind_kresA = cA,
                           ind_kresB = cB,
                           USE.NAMES = FALSE,
                           simplify = "matrix")
      count_vec <- apply(count_matr, 1, sum)
    }
    
    count_vec
    
  }
  parallel::stopCluster(cl)
  
  # sum up res_matr into final count vector
  count_vector <- apply(res_matr, MARGIN = 1, sum)
  return(count_vector)
  
}

#~~~~~~~~~~~~~~~~~~~~~~
#multik scoring wrapper
#######################
multik_scoring_wrapper <- function(k_res, num_cores, n_chunks = 20) {
  if (num_cores == 1) {
    
    #run_sequential
    k_res_chunks <- list(k_res)
    count_matr <- sapply(1:(nrow(k_res_chunks[[cA]]) - 1),
                         FUN = multik_scoring_par,
                         k_res_chunks = k_res_chunks,
                         mode = "self",
                         ind_kresA = cA,
                         ind_kresB = cA,
                         USE.NAMES = FALSE,
                         simplify = "matrix")
    count_vec <- apply(count_matr, 1, sum)
    return(count_vec)
    
  } else if (num_cores > 1 & num_cores <= 125) {
    
    # produce chunks and run parallel
    k_res_chunks <- split_chunks(k_res = k_res, n_chunks = n_chunks)
    count_vec <- multik_score_chunks(k_res_chunks = k_res_chunks, num_cores = num_cores)
    return(count_vec)
    
  } else {
    stop("num_cores must be an integer between 1 and 125!")
  }
}

#~~~~~~~~~~~~~~~~~~~~~
# calculate rPAC score
######################
CalcPAC_par <- function(x1 = 0.1, x2 = 0.9, # threshold defining the intermediate sub-interval
                    k_vec, count_vecs) {
  
  PAC <- rep(NA, length(k_vec))
  names(PAC) <- k_vec #
  
  prop_zeroes <- rep(NA, length(k_vec))
  names(prop_zeroes) <- k_vec
  
  rPAC <- rep(NA, length(k_vec))
  names(rPAC) <- k_vec #
  
  for(k in k_vec) {
    count_vec <- count_vecs[[k]]
    hist <- list(breaks = seq(0, 1, 0.001), 
                 counts = count_vec[1:length(count_vec) - 1], 
                 density = count_vec[1:length(count_vec) - 1] / (1:1000),
                 xname="count_vec_hist")
    
    class(hist) <- "histogram"
    Fn <- HistToEcdf(hist)
    
    # calculate PAC
    PAC[k] <- Fn(x2) - Fn(x1)
    # calculate proportion of zeroes in count_vec
    prop_zeroes[k] <- count_vec[1] / sum(count_vec)
    # calculate relative PAC
    rPAC[k] <- PAC[k]/(1-prop_zeroes[k])
  }
  return(list("min PAC" = which.min(PAC),
              "min rPAC" = which.min(rPAC),
              "PAC" = PAC,
              "rPAC" = rPAC))
}

#~~~~~~~~~~~~~~~~~~~~~~~
# MultiK diagnostic plot
########################
DiagMultiKPlot_par = function(ks, res_by_k, count_vecs) {
  
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))
  suppressPackageStartupMessages(library(cowplot))
  
  # set a data frame for the data to plot
  tog <- as.data.frame(table(ks)[table(ks) > 1])
  
  # calculate rpac
  pacobj <- CalcPAC_par(x1 = 0.1, x2 = 0.9, k_vec = names(count_vecs), count_vecs = count_vecs)
  tog$rpac <- pacobj$rPAC
  tog$one_minus_rpac  <- 1-tog$rpac
  
  # Plot Freq of # of runs
  freqPlot <- ggplot(data=tog, aes(x=ks, y=Freq)) +
    geom_bar(stat="identity") +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    geom_text(aes(label=Freq),  vjust=-0.3, size=3.5) +
    scale_x_discrete("K") +
    scale_y_continuous("Number of clustering runs") +
    geom_hline(yintercept=100, linetype="dashed", color = "black")
  
  # Plot rPAC for each K
  rpacPlot <- ggplot(data=tog, aes(x=ks, y=rpac,group=1)) +
    geom_point(shape=21, color="black", fill="black", size=2) +
    geom_line() +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    scale_x_discrete("K") +
    scale_y_continuous("rPAC")
  
  # Plot (1-rPAC) Vs freq for each K
  # first find optimal K using the convex hull
  optK <- findOptK(tog)
  cat("Optimal K: ", optK)
  
  scatPlot <- ggplot(data=tog, aes(x=one_minus_rpac, y=Freq)) +
    geom_point(shape=21, color="black", fill="black", size=1.5) +
    geom_path(color="grey", alpha=0.75, linetype=2) +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    scale_x_continuous("1 - rPAC") +
    scale_y_continuous("Number of clustering runs") +
    geom_hline(yintercept=100, linetype="dashed", color = "black") +
    geom_label_repel(aes(label = ks), segment.color = 'grey50', size=3)  +
    geom_path(data=tog[match(findOptK(tog), tog$ks), ])
  
  plot_grid(freqPlot, rpacPlot, scatPlot, ncol=3)
  
}


#########################
# Find optimal K #
#########################
findOptK = function(tog) {
  
  # FAILED TEST FOR GETTING RID OF NON-SIGNIFICANT Ks
  #tog.f <- tog[tog$Freq > 100 | tog$Freq ==100, ]
  tog.f = tog
  hpts <- chull(tog.f[, c("one_minus_rpac", "Freq")]) # in clockwise order
  hpts <- c(hpts, hpts[1])
  ch.df <- tog.f[hpts, ]
  
  df <- ch.df[ , c("ks", "one_minus_rpac", "Freq")]
  colnames(df) <- c("k", "x", "y")
  b <- c()
  end_points <- c()
  for (i in 1: (nrow(df)-1)) {
    end_points[i] <- paste(as.character(df[i, ]$k), as.character(df[(i+1),]$k), sep="-")
    b[i] <- (df[(i+1), ]$y - df[i, ]$y)/(df[(i+1), ]$x - df[i, ]$x)
  }
  
  # put in data frame for each line segment
  lineseg.df <- data.frame("end_points"=end_points, "slope"=b)
  lineseg.df$p1 <- do.call("rbind", strsplit(lineseg.df$end_points, "-"))[, 1]
  lineseg.df$p2 <- do.call("rbind", strsplit(lineseg.df$end_points, "-"))[, 2]
  
  # step 1: find k with largest # of runs
  which.k <- as.character(ch.df[which.max(ch.df$Freq), ]$ks)
  
  # step 2: see if a line segment with negative slope coming out
  if ( all(lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ]$slope > 0) ) {
    optK <- which.k
  }
  else {
    
    # follow the line segment with the negative slope
    tmp <- which(lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ]$slope < 0)
    tmp <- lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ][tmp, ]
    which.k2 <- as.character(c(tmp$p1, tmp$p2)[which(c(tmp$p1, tmp$p2)!=which.k)])
    
    # check if the slope becomes more negative
    lineseg.df.sub <- lineseg.df[lineseg.df$p1!=which.k & lineseg.df$p2 !=which.k, ]
    
    if ( #any(lineseg.df[lineseg.df$p1==which.k2 | lineseg.df$p2==which.k2, ]$slope > 0)
      lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2 == which.k2, ]$slope > tmp$slope ) {
      optK <- c(which.k, which.k2)
    }
    
    else {
      tmp <- which(lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2==which.k2, ]$slope < 0)
      tmp <- lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2==which.k2, ][tmp, ]
      which.k3 <- as.character(c(tmp$p1, tmp$p2)[ which(c(tmp$p1, tmp$p2)!=which.k & c(tmp$p1, tmp$p2)!=which.k2)])
      
      lineseg.df.sub <- lineseg.df[lineseg.df$p1!=which.k & lineseg.df$p2 !=which.k
                                   & lineseg.df$p1!=which.k2 & lineseg.df$p2 !=which.k2, ]
      
      if ( lineseg.df.sub[lineseg.df.sub$p1==which.k3 | lineseg.df.sub$p2 == which.k3, ]$slope > tmp$slope ) {
        optK <- c(which.k, which.k2, which.k3)
      }
      else {
        optK <- c(which.k, which.k2, which.k3)
      }
    }
  }
  return(optK)
}



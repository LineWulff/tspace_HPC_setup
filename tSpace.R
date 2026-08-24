tSpace <- function (df, K = 20, L = NULL, D = "pearson_correlation", graph = 5, 
    trajectories = 200, wp = 20, ground_truth = F, weights = "exponential", 
    dr = "pca", seed = NULL, core_no = 1, env_path=.libPaths(), ...) 
{
    if (any(sapply(df, function(x) is.numeric(x)) == F)) {
        stop("Check that all values in your data frame are numeric")
    }
    if (!is.numeric(core_no)) {
        stop("Number of cores is not numeric")
    }
    if (!is.numeric(K) | !is.numeric(graph) | !is.numeric(wp) | 
        !is.numeric(trajectories)) {
        stop("K, graph, waypoints or trajectories variables are not numbers")
    }
    if (!(D %in% c("euclidean", "manhattan", "chebyshev", "canberra", 
        "braycurtis", "pearson_correlation", "simple_matching_coefficient", 
        "minkowski", "hamming", "mahalanobis", "jaccard_coefficient", 
        "Rao_coefficient"))) {
        stop("distance can be any of 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient', 'minkowski',  'hamming', 'mahalanobis', 'jaccard_coefficient', 'Rao_coefficient'")
    }
    if (!(weights %in% c("uniform", "linear", "quadratic", "exponential"))) {
        stop("weights can be any of 'uniform', 'linear', 'quadratic', 'exponential'")
    }
    if (!(dr %in% c("pca", "umap", "both"))) {
        stop("dimensionality reduction can be any of 'pca', 'umap', 'both'")
    }
    if (is.null(seed)) {
        seed <- 1111
    }
    if (is.null(L)) {
        L = as.numeric(round(0.75 * K, digits = 0))
    }
    if (ground_truth == T) {
        numPop <- 1:nrow(df)
        tspacem <- matrix(data = NA, nrow = dim(df)[1], ncol = length(unique(numPop)))
        graph_panel <- list()
        s <- numPop
        trajectories <- length(numPop)
        Index <- numPop
    }
    else {
        numPop <- kmeans(df, centers = trajectories, iter.max = 10000)
        Index <- seq(1, nrow(df), by = 1)
        tspacem <- matrix(data = NA, nrow = dim(df)[1], ncol = length(unique(numPop$cluster)))
        graph_panel <- list()
        s <- unlist(lapply(split(Index, as.factor(numPop$cluster)), 
            function(x) {
                as.numeric(x[1])
            }))
    }
    cat(paste0("Step 0: debug - ",env_path),"\n")
    cat(paste0("Step 1:Finding graph"))
    knn <- graphfinder(x = df, k = K, distance = D, core_n = core_no, env_path=env_path)
    if (min(knn) < 0) {
        cat(paste0("\nSome cell-cell pairs have negative distances. \nIn order for this analysis to proceed, these distances will be aproximated to zero. \nIf substantial number of cells exhibit negative distances\nplease check the original data and examine which cells are causing the issue. \nSome of these cells may be just noise and should be removed"))
        negative.distance <- knn[which(knn[, 3] < 0), ]
        knn[which(knn[, 3] < 0), 3] <- 0
    }
    knn <- igraph::get.adjacency(igraph::graph.adjacency(Matrix::sparseMatrix(i = knn[, 
        "I"], j = knn[, "J"], x = knn[, "D"]), mode = "max", 
        weighted = TRUE), attr = "weight")
    cat(paste0("\nStep 2: Finding trajectories in sub-graphs \nCalculation may take time, don't close R\n"))
    graph_panel <- list()
    percentage <- seq(100/graph, 100, by = 100/graph)
    tictoc::tic("graphs_loop\n")
    for (graph_iter in 1:graph) {
        svMisc::progress(percentage[graph_iter], progress.bar = T)
        if (K != L) {
            l.knn = find_lknn(knn, l = L, core_n = core_no, env_path=env_path)
        }
        else {
            l.knn = knn
        }
        cl <- parallel::makeCluster(core_no)
        #clusterEvalQ(cl, .libPaths(lib_path))
        e <- new.env()
        e$libs <- c(env_path, .libPaths())
        source("/work/mccoy_lab/0000_generalscripts/tSp_helperfunctions.R")
        clusterExport(cl, "libs", envir=e)
        clusterEvalQ(cl, .libPaths(libs))
        clusterExport(cl, "pathfinder")
        doParallel::registerDoParallel(cl)

        tspacem <- foreach::foreach(i = 1:trajectories, .combine = cbind, 
            .packages = c("igraph", "Matrix", "KernelKnn", "pracma")) %dopar% 
            {
                s_c = as.numeric(s[i])
                tspacem <- pathfinder(data = df, lknn = l.knn, 
                  s = s_c, waypoints = wp, voting_scheme = weights, 
                  distance = D)$final_trajectory
            }
        parallel::stopCluster(cl)
        graph_panel[[graph_iter]] <- tspacem
        if (graph_iter == graph) {
            cat("\nLast sub-graph is done!")
        }
    }
    time <- tictoc::toc()
    rm(tspacem)
    tspace_mat <- array(unlist(graph_panel), c(nrow(graph_panel[[1]]), 
        ncol(graph_panel[[1]]), graph))
    rm(graph_panel)
    tspace_mat <- rowMeans(tspace_mat, dims = 2)
    colnames(tspace_mat) <- paste0("T_", s)
    cat(paste0("\nStep 3: Low dimensionality embbeding for visualization step\n"))
    if (dr == "pca") {
        if (ncol(tspace_mat) > 20) {
            set.seed(seed)
            pca_tspace <- prcomp(t(tspace_mat), rank. = 20)
        }
        else {
            set.seed(seed)
            pca_tspace <- prcomp(t(tspace_mat), center = T, scale. = T)
        }
        pca_out <- pca_tspace$rotation
        colnames(pca_out) <- paste0("tPC", seq(1, ncol(pca_out), 
            1))
        data.out <- as.data.frame(cbind(Index = Index, pca_out, 
            df))
        if (exists("negative.distance") == T) {
            tspace_obj <- list(ts_file = data.out, pca_embbeding = pca_tspace, 
                tspace_matrix = tspace_mat, negative_distances = negative.distance)
        }
        else {
            tspace_obj <- list(ts_file = data.out, pca_embbeding = pca_tspace, 
                tspace_matrix = tspace_mat)
        }
    }
    if (dr == "umap") {
        config_tspace <- umap::umap.defaults
        config_tspace$n_neighbors <- 7
        config_tspace$min_dist <- 0.2
        config_tspace$metric <- "manhattan"
        set.seed(seed)
        umap_tspace <- umap::umap(tspace_mat, config = config_tspace)
        umap_out <- as.data.frame(umap_tspace$layout)
        colnames(umap_out) <- paste0("umap", seq(1, ncol(umap_tspace$layout), 
            1))
        data.out <- as.data.frame(cbind(Index = Index, umap_out, 
            df))
        if (exists("negative.distance") == T) {
            tspace_obj <- list(ts_file = data.out, umap_embbeding = umap_tspace, 
                tspace_matrix = tspace_mat, negative_distances = negative.distance)
        }
        else {
            tspace_obj <- list(ts_file = data.out, umap_embbeding = umap_tspace, 
                tspace_matrix = tspace_mat)
        }
    }
    if (dr == "both") {
        if (ncol(tspace_mat) > 20) {
            set.seed(seed)
            pca_tspace <- prcomp(t(tspace_mat), rank. = 20)
        }
        else {
            set.seed(seed)
            pca_tspace <- prcomp(t(tspace_mat), center = T, scale. = T)
        }
        pca_out <- as.data.frame(pca_tspace$rotation)
        colnames(pca_out) <- paste0("tPC", seq(1, ncol(pca_out), 
            1))
        config_tspace <- umap::umap.defaults
        config_tspace$n_neighbors <- 7
        config_tspace$min_dist <- 0.2
        config_tspace$metric <- "manhattan"
        # 2D vs 3D umap outputs
        config_tspace$n_components <- 3
        set.seed(seed)
        umap_tspace <- umap::umap(tspace_mat, config = config_tspace)
        umap_out <- as.data.frame(umap_tspace$layout)
        colnames(umap_out) <- paste0("umap", seq(1, ncol(umap_tspace$layout), 
            1))
        data.out <- as.data.frame(cbind(Index = Index, pca_out, 
            umap_out, df))
        if (exists("negative.distance") == T) {
            tspace_obj <- list(ts_file = data.out, pca = pca_tspace, 
                umap_embbeding = umap_tspace, tspace_matrix = tspace_mat, 
                negative_distances = negative.distance)
        }
        else {
            tspace_obj <- list(ts_file = data.out, pca = pca_tspace, 
                umap_embbeding = umap_tspace, tspace_matrix = tspace_mat)
        }
    }
    return(tspace_obj)
}

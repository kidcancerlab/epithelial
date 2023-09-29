recurluster <- function(sobj,
                        level = 1,
                        max_level = 3,
                        min_cells = 0,
                        data_type = "RNA",
                        do_plots = TRUE,
                        parallel = FALSE) {
    message("Clustering level ", level)

    # can't do plots while parallelizing
    if (parallel && do_plots) {
        do_plots <- FALSE
        warning("Cannot plot while operating in parallel")
    }

    # number of PCs needs to be less than number of cells
    pc_num <- min(50, length(Seurat::Cells(sobj)) - 1)

    # Recluster at course resolution
    if (data_type == "RNA") {
        temp_sobj <-
            sobj %>%
            SCTransform(vars.to.regress = "percent.mt",
                        verbose = FALSE,
                        return.only.var.genes = FALSE) %>%
            RunPCA(npcs = pc_num,
                   approx = FALSE,
                   verbose = FALSE,
                   assay = "SCT") %>%
            FindNeighbors(verbose = FALSE) %>%
            RunUMAP(dims = 1:min(10, pc_num),
                    n.neighbors = min(30, pc_num / 2),
                    verbose = FALSE)

        opt_res <-
            optimize_silhouette(temp_sobj,
                                test_res = seq(0.05, 0.3, by = 0.05),
                                summary_plot = FALSE) %>%
            arrange(sil_vals * -1) %>%
            pull(res_vals) %>%
            head(n = 1)

        temp_sobj <- temp_sobj %>%
            Seurat::FindClusters(resolution = opt_res, verbose = FALSE)
            message("Optimal resolution: ", opt_res)
    } else if (data_type == "ATAC") {
        temp_sobj <- sobj %>%
            RunTFIDF(verbose = FALSE) %>%
            FindTopFeatures(min.cutoff = "q25", verbose = FALSE) %>%
            RunSVD(verbose = FALSE) %>%
            Seurat::RunUMAP(reduction = "lsi",
                            dims = 2:min(30, pc_num),
                            verbose = FALSE) %>%
            FindNeighbors(reduction = "lsi", verbose = FALSE) %>%
            FindClusters(algorithm = 3, resolution = res, verbose = FALSE)
    }

    # Add new clusters to the Seurat object, at approprate level
    if (level == 1) {
        temp_sobj[[paste0("clust_", level)]] <-
            paste0("clust_", temp_sobj$seurat_clusters)
    } else {
        temp_sobj[[paste0("clust_", level)]] <-
            paste0(temp_sobj[[paste0("clust_", level - 1)]][, 1],
                   temp_sobj$seurat_clusters)
    }

    # Plot the clusters at this level
    if (do_plots) {
        Seurat::Idents(temp_sobj) <- temp_sobj[[paste0("clust_", level)]]
        print(Seurat::DimPlot(temp_sobj, pt.size = 1) +
            ggtitle(paste0("Level ", level))) +
            theme(plot.title = element_text(hjust = 0.5))
    }

    # If more than one cluster isolate and recurluster each cluster
    # Also stop if we've reached the max level
    if (temp_sobj$seurat_clusters %>%
            as.character() %>%
            as.numeric() %>%
            max() > 0 &&
        level < max_level) {

        if (!parallel) {
            sobj_list <-
                lapply(unique(temp_sobj$seurat_clusters),
                       function(x) re_recurluster(sobject = temp_sobj,
                                                  cluster_num = x,
                                                  level = level,
                                                  min_cells = min_cells,
                                                  data_type = data_type,
                                                  do_plots = do_plots,
                                                  parallel = parallel))
        } else {
            num_cores <- parallel::detectCores()
            sobj_list <-
                parallel::mclapply(unique(temp_sobj$seurat_clusters),
                                   function(x) re_recurluster(sobject = temp_sobj,
                                                              cluster_num = x,
                                                              level = level,
                                                              min_cells = min_cells,
                                                              data_type = data_type,
                                                              do_plots = do_plots,
                                                              parallel = parallel),
                                   mc.cores = num_cores)
        }
        # Merge the subclusters into a single Seurat object
        temp_sobj <- merge(sobj_list[[1]],
                        sobj_list[2:length(sobj_list)])
    }

    # Fill in NAs in the clustser columns with parent cluster plus "_0"
    for (clust_col in grep("clust_",
                           names(temp_sobj@meta.data),
                           value = TRUE) %>%
            grep("clust_1", ., value = TRUE, invert = TRUE)) {
        clust_num <- str_remove(clust_col, "clust_") %>%
            as.numeric()

        temp_sobj[[clust_col]] <-
            temp_sobj@meta.data %>%
            as.data.frame() %>%
            mutate({{ clust_col }} := if_else(is.na(get(clust_col)),
                                              paste0(get(paste0("clust_",
                                                                clust_num - 1)),
                                                     "0"),
                                              get(clust_col))) %>%
            pull(get(clust_col))
    }
    return(temp_sobj)
}

re_recurluster <- function(sobject,
                           cluster_num,
                           level,
                           min_cells,
                           data_type,
                           do_plots,
                           parallel) {
    if (length(which(sobject$seurat_clusters == cluster_num)) > min_cells) {
                message("Processing cluster ", cluster_num, " at level ", level)
        sobj_out <-
            recurluster(subset(sobject,
                                subset = seurat_clusters == cluster_num),
                        level = level + 1,
                        min_cells = min_cells,
                        data_type = data_type,
                        do_plots = do_plots,
                        parallel = parallel)
    } else {
        sobj_out <-
            subset(sobject,
                   subset = seurat_clusters == cluster_num)
    }
}

recurluster_fill_nas <- function(sobject) {
    for (clust_col in grep("clust_",
                        names(sobject@meta.data),
                        value = TRUE) %>%
        grep("clust_1", ., value = TRUE, invert = TRUE)) {
        clust_num <- str_remove(clust_col, "clust_") %>%
            as.numeric()

        sobject[[clust_col]] <-
            sobject@meta.data %>%
            as.data.frame() %>%
            mutate({{ clust_col }} := if_else(is.na(get(clust_col)),
                                                paste0(get(paste0("clust_",
                                                                clust_num - 1)),
                                                        "0"),
                                                get(clust_col))) %>%
            pull(get(clust_col))
    }
}

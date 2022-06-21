UMAP_SCT <- function(seurat_data, do.SCT = TRUE, assay = "RNA", umap.n.neighbors = 25,
    umap.min.dist = 0.4, leiden.res = 0.5, algorithm = 1, umap.components = 2, umap.metric = "correlation") {
    if (do.SCT) {
        seurat_data <- SCTransform(seurat_data, assay = assay, verbose = FALSE, method = "glmGamPoi")  # vars.to.regress = c('Rep'))
    }

    # These are now standard steps in the Seurat workflow for visualization and
    # clustering
    seurat_data <- RunPCA(seurat_data, verbose = FALSE)
    seurat_data <- RunUMAP(seurat_data, dims = 1:30, verbose = FALSE, n.neighbors = umap.n.neighbors,
        min.dist = umap.min.dist, n.components = umap.components, metric = umap.metric)

    seurat_data <- FindNeighbors(seurat_data, dims = 1:30, , k.param = 20, verbose = FALSE)
    seurat_data <- FindClusters(seurat_data, resolution = leiden.res, algorithm = algorithm, 
        verbose = FALSE)

    seurat_data
}

UMAP_LogScale <- function(seurat_data, do.NORM = TRUE, assay = "RNA", n.high.var.features = 2000,
    umap.n.neighbors = 25, umap.min.dist = 0.4, leiden.res = 0.5) {
    if (do.NORM) {
        seurat_data <- NormalizeData(seurat_data, assay = assay, verbose = FALSE)  # vars.to.regress = c('Rep'))
        seurat_data <- FindVariableFeatures(seurat_data, selection.method = "vst",
            nfeatures = n.high.var.features)
        seurat_data <- ScaleData(seurat_data, verbose = FALSE)
    }

    # These are now standard steps in the Seurat workflow for visualization and
    # clustering
    seurat_data <- RunPCA(seurat_data, verbose = FALSE)
    seurat_data <- RunUMAP(seurat_data, dims = 1:30, verbose = FALSE, n.neighbors = umap.n.neighbors,
        min.dist = umap.min.dist, metric = "correlation")

    seurat_data <- FindNeighbors(seurat_data, dims = 1:30, verbose = FALSE)
    seurat_data <- FindClusters(seurat_data, resolution = leiden.res, verbose = FALSE)

    seurat_data
}

IntegrationbyGroup <- function(to_combine, split.by = "group", do.SCT = FALSE, n.features.sct = 5000,
    n.features.integration = 5000, k.weight = 200, reference = NULL, ...) {
    list_to_combine <- SplitObject(to_combine, split.by = split.by)
    if (do.SCT) {
        cat("Performing SCT transform\n")
        list_to_combine <- lapply(X = list_to_combine, FUN = function(obj) SCTransform(obj,
            verbose = FALSE, method = "glmGamPoi", variable.features.n = n.features.sct))
    }

    cat("Finding anchors for integration\n")
    cat(names(list_to_combine))
    cat(subset(names(list_to_combine), !grepl("E4.5", names(list_to_combine))))
    features <- SelectIntegrationFeatures(object.list = list_to_combine, nfeatures = n.features.integration)
    list_to_combine <- PrepSCTIntegration(object.list = list_to_combine, anchor.features = features)
    integrate.anchors <- FindIntegrationAnchors(object.list = list_to_combine, normalization.method = "SCT",
        anchor.features = features, scale = FALSE, reference = reference, verbose = FALSE
    )

    # determine k weight
    k.weight.sample <- min(unlist(lapply(list_to_combine, ncol)))
    if (k.weight.sample < k.weight) {
        k.weight = k.weight.sample
    }
    cat("Integrating\n")
    combined.sct <- IntegrateData(anchorset = integrate.anchors, normalization.method = "SCT",
        k.weight = k.weight, verbose = FALSE, ...)

    combined.sct
}


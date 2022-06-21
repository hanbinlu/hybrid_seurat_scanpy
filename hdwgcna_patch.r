# Modify code from hdWGCNA to adjust to sct integrated workflow See comments at
# places of modifications for more details


#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_net$colors
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are 'linear', 'poisson', or 'negbinom'
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' ComputeModuleEigengene(pbmc)
ComputeModuleEigengene.1 <- function(seurat_obj, cur_mod, modules, group.by.vars = NULL,
                                     verbose = TRUE, vars.to.regress = NULL, scale.model.use = "linear", wgcna_name = NULL,
                                     assay = "RNA", slot = "data", do.scale.data = FALSE, ...) {
    ##### Add control for assay to use for compute ME, you can opt out
    ##### ScaleData as well #####

    # set as active assay if wgcna_name is not given
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }

    # get genes in this module:
    cur_genes <- modules %>%
        subset(module == cur_mod) %>%
        .$gene_name

    # subset seurat object by these genes only:
    cur_seurat <- seurat_obj[cur_genes, ]

    # MODIFIED: opt out ScaleData
    if (do.scale.data) {
        # scale the subsetted expression dataset:
        if (is.null(vars.to.regress)) {
            cur_seurat <- ScaleData(cur_seurat,
                features = rownames(cur_seurat),
                model.use = scale.model.use
            )
        } else if (all(vars.to.regress %in% colnames(seurat_obj@meta.data))) {
            cur_seurat <- ScaleData(cur_seurat,
                features = rownames(cur_seurat),
                model.use = scale.model.use, vars.to.regress = vars.to.regress
            )
        } else {
            stop(paste0("Some variables specified in vars.to.regress are not found in seurat_obj@meta.data"))
        }
    }

    # MODIFIED: Assay to use. We want to use integrated data to compute ME
    # compute average expression of each gene
    DefaultAssay(cur_seurat) <- assay
    cur_expr <- GetAssayData(cur_seurat, slot = slot)
    expr <- t(as.matrix(cur_expr))
    averExpr <- rowSums(expr) / ncol(expr)

    # run PCA with Seurat function This should be run on selected assay, e.g.
    # 'integrated'
    cur_pca <- Seurat::RunPCA(cur_seurat, features = cur_genes, reduction.key = paste0(
        "pca",
        cur_mod
    ), verbose = verbose, ...)@reductions$pca
    pc <- cur_pca@cell.embeddings[, 1]
    pc_loadings <- cur_pca@feature.loadings[, 1]

    # correlate average expression with eigengene
    pca_cor <- cor(averExpr, pc)

    # run harmony ATTENTION: whether the chose assay justified to perform
    # harmony, if not, keep `group.by.vars` NULL
    if (!is.null(group.by.vars)) {

        # add this PCA as its own reduction in the seurat object
        seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
            embeddings = cur_pca@cell.embeddings,
            assay = Seurat::DefaultAssay(cur_seurat)
        )

        cur_harmony <- harmony::RunHarmony(seurat_obj,
            group.by.vars = group.by.vars,
            reduction = "ME", verbose = verbose, ...
        )@reductions$harmony
        ha <- cur_harmony@cell.embeddings[, 1]
        ha_loadings <- cur_pca@feature.loadings[, 1]

        if (pca_cor < 0) {
            cur_harmony@cell.embeddings[, 1] <- -ha
            ha_loadings <- -ha_loadings
        }

        # add harmonized PCA as its own reduction in the seurat object
        seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(
            embeddings = cur_harmony@cell.embeddings,
            assay = Seurat::DefaultAssay(seurat_obj)
        )

        seurat_obj <- SetMELoadings(seurat_obj,
            loadings = ha_loadings, harmonized = TRUE,
            wgcna_name = wgcna_name
        )
    }

    if (pca_cor < 0) {
        cur_pca@cell.embeddings[, 1] <- -pc
        pc_loadings <- -pc_loadings
    }

    # add this PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
        embeddings = cur_pca@cell.embeddings,
        assay = Seurat::DefaultAssay(cur_seurat)
    )

    seurat_obj <- SetMELoadings(seurat_obj,
        loadings = pc_loadings, harmonized = FALSE,
        wgcna_name = wgcna_name
    )

    # return seurat object
    seurat_obj
}

#' ModuleEigengenes
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are 'linear', 'poisson', or 'negbinom'
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleEigengenes(pbmc)
ModuleEigengenes.1 <- function(seurat_obj, group.by.vars = NULL, modules = NULL,
                               vars.to.regress = NULL, scale.model.use = "linear", verbose = TRUE, wgcna_name = NULL,
                               assay = "RNA", slot = "data", do.scale.data = FALSE, ...) {
    ##### Add control for assay to use for compute ME, you can opt out
    ##### ScaleData as well #####

    # set as active assay if wgcna_name is not given
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }

    # are we going to run Harmony?
    harmonized <- !is.null(group.by.vars)

    me_list <- list()
    harmonized_me_list <- list()

    # re-set feature loadings:
    seurat_obj <- SetMELoadings(seurat_obj,
        loadings = c(""), harmonized = FALSE,
        wgcna_name = wgcna_name
    )
    if (harmonized) {
        seurat_obj <- SetMELoadings(seurat_obj,
            loadings = c(""), harmonized = TRUE,
            wgcna_name = wgcna_name
        )
    }

    # get modules from Seurat object, else use provided modules
    if (is.null(modules)) {
        modules <- GetModules(seurat_obj, wgcna_name)
        projected <- FALSE
    } else {
        projected <- TRUE
    }

    # get list of modules:
    mods <- levels(modules$module)

    # loop over modules:
    for (cur_mod in mods) {
        print(cur_mod)

        # compute module eigengenes for this module MODIFIED: add assay options
        seurat_obj <- ComputeModuleEigengene.1(
            seurat_obj = seurat_obj, cur_mod = cur_mod,
            modules = modules, group.by.vars = group.by.vars, vars.to.regress = vars.to.regress,
            scale.model.use = scale.model.use, verbose = verbose, wgcna_name = wgcna_name,
            assay = assay, slot = slot, do.scale.data = do.scale.data, ...
        )

        # add module eigengene to ongoing list
        cur_me <- seurat_obj@reductions$ME@cell.embeddings[, 1]
        me_list[[cur_mod]] <- cur_me

        # run harmony
        if (harmonized) {
            # add module eigengene to ongoing list
            cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[
                ,
                1
            ]
            harmonized_me_list[[cur_mod]] <- cur_harmonized_me
        }
    }

    # merge module eigengene lists into a dataframe, order modules, add to
    # Seurat obj
    me_df <- do.call(cbind, me_list)
    if (!projected) {
        me_df <- WGCNA::orderMEs(me_df)
    }
    seurat_obj <- SetMEs(seurat_obj, me_df, harmonized = FALSE, wgcna_name)

    # merge harmonized module eigengene lists into a dataframe, add to Seurat
    # obj
    if (!is.null(group.by.vars)) {
        hme_df <- do.call(cbind, harmonized_me_list)
        if (!projected) {
            hme_df <- WGCNA::orderMEs(hme_df)
        }
        seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized = TRUE, wgcna_name)
    }

    # set module factor levels based on order
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
    modules$module <- factor(as.character(modules$module), levels = colnames(MEs))
    seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

    # remove temp dim reductions by setting to NULL
    seurat_obj@reductions$ME <- NULL
    seurat_obj@reductions$ME_harmony <- NULL

    # return seurat object
    seurat_obj
}

###### Based on Meta Cell ######

# Aggregate expression for meta cell
MetaCellDatExpr <- function(seurat_obj, wgcna_name, group.by = NULL, group.names = NULL) {
    s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
    if (!(group.by %in% colnames(s_obj@meta.data))) {
        m_cell_message <- "metacell"
        stop(paste0(
            group.by, " not found in the meta data of the ", m_cell_message,
            " Seurat object"
        ))
    }

    # get the metadata from the seurat object:
    seurat_meta <- s_obj@meta.data
    # columns to group by for cluster/celltype
    if (!is.null(group.by)) {
        seurat_meta <- seurat_meta %>%
            subset(get(group.by) %in% group.names)
    }

    # get list of cells to use
    cells <- rownames(seurat_meta)

    # get expression data from seurat obj Notice that this is a meta cell obj
    datExpr <- as.data.frame(Seurat::GetAssayData(s_obj, assay = "RNA", slot = "data")[
        gene_list,
        cells
    ])
    # transpose data
    datExpr <- as.data.frame(t(datExpr))

    gene_list <- gene_list[WGCNA::goodGenes(datExpr)]
    datExpr <- datExpr[, gene_list]

    # update the WGCNA gene list:
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)
    seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
    seurat_obj
}

# Modify code from hdWGCNA to adjust to sct integrated workflow See comments at
# places of modifications for more details


#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_net$colors
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are 'linear', 'poisson', or 'negbinom'
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' ComputeModuleEigengene(pbmc)
ComputeModuleEigengene.Meta <- function(seurat_obj, cur_mod, modules, group.by.vars = NULL,
                                        verbose = TRUE, vars.to.regress = NULL, scale.model.use = "linear", wgcna_name = NULL,
                                        assay = "RNA", slot = "counts", do.scale.data = FALSE, ...) {
    ##### Add control for assay to use for compute ME, you can opt out
    ##### ScaleData as well #####

    # set as active assay if wgcna_name is not given
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }

    # get genes in this module:
    cur_genes <- modules %>%
        subset(module == cur_mod) %>%
        .$gene_name

    # subset seurat object by these genes only:
    s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
    cur_seurat <- s_obj[cur_genes, ]

    # MODIFIED: opt out ScaleData
    if (do.scale.data) {
        # scale the subsetted expression dataset:
        if (is.null(vars.to.regress)) {
            cur_seurat <- ScaleData(cur_seurat,
                features = rownames(cur_seurat),
                model.use = scale.model.use
            )
        } else if (all(vars.to.regress %in% colnames(seurat_obj@meta.data))) {
            cur_seurat <- ScaleData(cur_seurat,
                features = rownames(cur_seurat),
                model.use = scale.model.use, vars.to.regress = vars.to.regress
            )
        } else {
            stop(paste0("Some variables specified in vars.to.regress are not found in seurat_obj@meta.data"))
        }
    }

    # MODIFIED: Assay to use. We want to use integrated data to compute ME
    # compute average expression of each gene
    DefaultAssay(cur_seurat) <- assay
    cur_expr <- GetAssayData(cur_seurat, slot = slot)
    expr <- t(as.matrix(cur_expr))
    averExpr <- rowSums(expr) / ncol(expr)

    # run PCA with Seurat function This should be run on selected assay, e.g.
    # 'integrated'
    cur_pca <- Seurat::RunPCA(cur_seurat, features = cur_genes, reduction.key = paste0(
        "pca",
        cur_mod
    ), verbose = verbose, ...)@reductions$pca
    pc <- cur_pca@cell.embeddings[, 1]
    pc_loadings <- cur_pca@feature.loadings[, 1]

    # correlate average expression with eigengene
    pca_cor <- cor(averExpr, pc)

    # run harmony ATTENTION: whether the chose assay justified to perform
    # harmony, if not, keep `group.by.vars` NULL
    # if (!is.null(group.by.vars)) {

    #    # add this PCA as its own reduction in the seurat object
    #    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
    #        embeddings = cur_pca@cell.embeddings,
    #        assay = Seurat::DefaultAssay(cur_seurat)
    #    )

    #    cur_harmony <- harmony::RunHarmony(seurat_obj,
    #        group.by.vars = group.by.vars,
    #        reduction = "ME", verbose = verbose, ...
    #    )@reductions$harmony
    #    ha <- cur_harmony@cell.embeddings[, 1]
    #    ha_loadings <- cur_pca@feature.loadings[, 1]

    #    if (pca_cor < 0) {
    #        cur_harmony@cell.embeddings[, 1] <- -ha
    #        ha_loadings <- -ha_loadings
    #    }

    #    # add harmonized PCA as its own reduction in the seurat object
    #    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(
    #        embeddings = cur_harmony@cell.embeddings,
    #        assay = Seurat::DefaultAssay(seurat_obj)
    #    )

    #    seurat_obj <- SetMELoadings(seurat_obj,
    #        loadings = ha_loadings, harmonized = TRUE,
    #        wgcna_name = wgcna_name
    #    )
    # }

    if (pca_cor < 0) {
        cur_pca@cell.embeddings[, 1] <- -pc
        pc_loadings <- -pc_loadings
    }

    # add this PCA as its own reduction in the seurat object
    cur_pca@cell.embeddings[, 1]
}

#' ModuleEigengenes
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are 'linear', 'poisson', or 'negbinom'
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleEigengenes(pbmc)
ModuleEigengenes.Meta <- function(seurat_obj, group.by.vars = NULL, modules = NULL,
                                  vars.to.regress = NULL, scale.model.use = "linear", verbose = TRUE, wgcna_name = NULL,
                                  assay = "RNA", slot = "counts", do.scale.data = FALSE, ...) {
    ##### Add control for assay to use for compute ME, you can opt out
    ##### ScaleData as well #####

    # set as active assay if wgcna_name is not given
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }

    # are we going to run Harmony?
    harmonized <- !is.null(group.by.vars)

    me_list <- list()
    harmonized_me_list <- list()

    # re-set feature loadings:
    seurat_obj <- SetMELoadings(seurat_obj,
        loadings = c(""), harmonized = FALSE,
        wgcna_name = wgcna_name
    )
    if (harmonized) {
        seurat_obj <- SetMELoadings(seurat_obj,
            loadings = c(""), harmonized = TRUE,
            wgcna_name = wgcna_name
        )
    }

    # get modules from Seurat object, else use provided modules
    if (is.null(modules)) {
        modules <- GetModules(seurat_obj, wgcna_name)
        projected <- FALSE
    } else {
        projected <- TRUE
    }

    # get list of modules:
    mods <- levels(modules$module)

    # loop over modules:
    for (cur_mod in mods) {
        print(cur_mod)

        # compute module eigengenes for this module MODIFIED: add assay options
        cur_me <- ComputeModuleEigengene.Meta(
            seurat_obj = seurat_obj, cur_mod = cur_mod,
            modules = modules, group.by.vars = group.by.vars, vars.to.regress = vars.to.regress,
            scale.model.use = scale.model.use, verbose = verbose, wgcna_name = wgcna_name,
            assay = assay, slot = slot, do.scale.data = do.scale.data, ...
        )

        me_list[[cur_mod]] <- cur_me

        # run harmony
        # if (harmonized) {
        #     # add module eigengene to ongoing list
        #     cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[
        #         ,
        #         1
        #     ]
        #     harmonized_me_list[[cur_mod]] <- cur_harmonized_me
        # }
    }

    # merge module eigengene lists into a dataframe, order modules, add to
    # Seurat obj
    me_df <- do.call(cbind, me_list)
    # if (!projected) {
    #     me_df <- WGCNA::orderMEs(me_df)
    # }
    seurat_obj <- SetMEs(seurat_obj, me_df, harmonized = FALSE, wgcna_name)

    # set module factor levels based on order
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
    modules$module <- factor(as.character(modules$module), levels = colnames(MEs))
    seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

    # return seurat object
    seurat_obj
}

ModuleConnectivity.Meta <- function(seurat_obj,
                                    harmonized = FALSE,
                                    group.by = NULL,
                                    group_name = NULL,
                                    wgcna_name = NULL,
                                    ...) {

    # set as active assay if wgcna_name is not given
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }

    # get module df, wgcna genes, and wgcna params:
    modules <- GetModules(seurat_obj, wgcna_name)
    # MODIFIED: use the scoring scheme as module signature
    # MEs <- seurat_obj@misc[[wgcna_name]]$module_scores
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
    genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
    params <- GetWGCNAParams(seurat_obj, wgcna_name)

    # exclude the grey module:
    # modules <- subset(modules, module != "grey")
    # MEs <- MEs[,colnames(MEs) != 'grey']

    if (!is.null(group.by)) {
        cells.use <- GetMetacellObject(seurat_obj, wgcna_name)@meta.data %>%
            subset(get(group.by) %in% group_name) %>%
            rownames()
        MEs <- MEs[cells.use, ]
    } else {
        cells.use <- colnames(seurat_obj)
    }

    # datExpr <- t(as.matrix(exp_mat))
    datExpr <- seurat_obj@misc[[wgcna_name]]$datExpr[cells.use, genes_use]
    print(dim(datExpr))
    print("running signedKME...")

    kMEs <- WGCNA::signedKME(
        datExpr,
        MEs,
        outputColumnName = "kME",
        corFnc = "bicor",
        ...
    )

    # add module color to the kMEs table
    modules <- modules[, 1:3]
    kMEs <- cbind(modules, kMEs)
    colnames(kMEs) <- c(colnames(modules), paste0("kME_", colnames(MEs)))

    # update the modules table in the Seurat object:
    seurat_obj <- SetModules(seurat_obj, kMEs, wgcna_name)

    seurat_obj
}
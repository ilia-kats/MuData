setGeneric("WriteH5MU", function(object, file, overwrite = TRUE) standardGeneric("WriteH5MU"))
# setGeneric("WriteH5AD", function(object, assay, file, overwrite = TRUE) standardGeneric("WriteH5AD"))

#' Save a MultiAssayExperimet to an .h5mu file
#'
#' Note than the primary key is used as obs_names
#' so the behaviour of WriteH5MU when there are multiple samples
#' for one primary key is not guaranteed.
#'
#' @importFrom rhdf5 H5Gcreate H5Gclose
#' @importMethodsFrom MultiAssayExperiment colData experiments sampleMap
#' @importMethodsFrom SummarizedExperiment colData assay
#' @importMethodsFrom SingleCellExperiment reducedDims
#'
#' @export
setMethod("WriteH5MU", "MultiAssayExperiment", function(object, file, overwrite) {
    h5 <- open_h5(file)

    obs <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    write_data_frame(h5, "obs", obs)

    modalities <- names(experiments(object))

    mods <- H5Gcreate(h5, "mod")
    vars <- lapply(modalities, function(mod) {
        mod_group <- H5Gcreate(mods, mod)

        # .obs
        meta <- sampleMap(object)[sampleMap(object)$assay == mod,]
        obs_global <- as.data.frame(colData(object)[meta$primary,], stringsAsFactors = FALSE)
        obs <- data.frame(row.names=rownames(obs_global))
        if (is(object[[mod]], "SummarizedExperiment")) {
            obs_local <- colData(object[[mod]])
            obs_local <- obs_local[meta$colname,,drop=FALSE]
        if (ncol(obs_local) > 0)
            obs <- cbind(obs, obs_local)
        }
        write_data_frame(mod_group, "obs", obs)

        # .obsm
        if (is(object[[mod]], "SingleCellExperiment")) {
            obsm <- reducedDims(object[[mod]])
            if (length(obsm) > 0) {
                obsmgrp <- H5Gcreate(mod_group, "obsm")
                mapply(function(name, data) {
                    if (is.data.frame(data)) {
                        rownames(data) <- rownames(obs_global)
                        write_data_frame(obsmgrp, name, data)
                    } else {
                        if (length(dim(data)) == 1)
                            data <- as.vector(data)
                        write_matrix(obsmgrp, name, data)
                    }
                }, names(obsm), obsm)
                H5Gclose(obsmgrp)
            }
        }

        # X
        x <- object[[mod]]
        x <- x[,meta$colname]
        write_matrix(mod_group, "X", t(assay(x)))

        # .var
        var <- data.frame("mod" = rep(mod, nrow(x)), row.names = rownames(x), stringsAsFactors = FALSE)
        write_data_frame(mod_group, "var", var)

        finalize_anndata_internal(mod_group)
        H5Gclose(mod_group)
        var
    })
    H5Gclose(mods)

    var <- do.call(rbind, vars)
    write_data_frame(h5, "var", var)

    finalize_mudata(h5)
    invisible(NULL)
})

#' @importFrom rhdf5 h5writeDataset h5writeAttribute H5Gcreate H5Gclose
write_matrix <- function(parent, key, mat) {
    if (is.matrix(mat) || is.vector(mat)) {
        h5writeDataset(mat, parent, key)
    } else if (is(mat, "dgCMatrix") || is(mat, "dgRMatrix")) {
        grp <- H5Gcreate(parent, key)
        h5writeDataset(mat@p, grp, "indptr")
        h5writeDataset(mat@x, grp, "data")
        h5writeAttribute(dim(mat), grp, "shape")
        h5writeAttribute("0.1.0", grp, "encoding-version", variableLengthString=TRUE, asScalar=TRUE)
        if (is(mat, "dgCMatrix")) {
            h5writeDataset(mat@i, grp, "indices")
            h5writeAttribute("csc_matrix", grp, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
        } else {
            h5writeDataset(mat@j, grp, "indices")
            h5writeAttribute("csr_matrix", grp, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
        }
    } else {
        stop(paste0("Writing matrices of type ", class(mat), " is not implemented."))
    }
}

#' @importFrom rhdf5 H5Gcreate H5Gclose h5writeDataset h5writeAttribute h5createAttribute H5Dclose
write_data_frame <- function(parent, key, df) {
    group <- H5Gcreate(parent, key)

    columns <- colnames(df)
    df["_index"] <- rownames(df)
    categories <- list()
    for (col in colnames(df)) {
        v <- df[[col]]
        if (is.factor(v)) {
            categories[[col]] <- levels(v)
            h5writeDataset(as.integer(v) - 1L, group, col)
        } else {
            h5writeDataset(v, group, col)
        }
    }
    if (length(categories) > 0) {
        catgrp <- H5Gcreate(group, "__categories")

        mapply(function(colname, categories) {
            h5writeDataset(categories, catgrp, colname, variableLengthString=TRUE)
            cat_dset <- catgrp & colname
            h5writeAttribute(as.integer(is.ordered(df[[colname]])), cat_dset, "ordered", asScalar=TRUE)
            h5writeAttribute(cat_dset, group & colname, "categories", asScalar=TRUE)
            H5Dclose(cat_dset)
        }, names(categories), categories)
        H5Gclose(catgrp)
    }

    # Write attributes
    h5writeAttribute("_index", group, "_index", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute("dataframe", group, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute("0.1.0", group, "encoding-version", variableLengthString=TRUE, asScalar=TRUE)
    if (length(columns) > 0) {
        h5writeAttribute(columns, group, "column-order", variableLengthString=TRUE)
    } else {
        # When there are no columns, null buffer can't be written to a file.
        h5createAttribute(group, "column-order", dims=0)
    }
    H5Gclose(group)
    invisible(NULL)
}

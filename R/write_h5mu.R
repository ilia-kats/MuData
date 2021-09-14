setGeneric("WriteH5AD", function(object, file, overwrite=TRUE) standardGeneric("WriteH5AD"), signature=c("object", "file"))
setGeneric("WriteH5MU", function(object, file, overwrite=TRUE) standardGeneric("WriteH5MU"), signature=c("object", "file"))

#' @importFrom rhdf5 H5Iget_type
setMethod("WriteH5AD", c(object="mMatrix", file="H5IdComponent"), function(object, file, overwrite) {
    if (!(H5Iget_type(file) %in% c("H5I_FILE", "H5I_GROUP")))
        stop("object must be a file or group")
    write_matrix(file, "X", object)
    finalize_anndata_internal(file)
})

#' @importFrom rhdf5 H5Iget_type
#' @importMethodsFrom SummarizedExperiment colData assay
setMethod("WriteH5AD", c(object="SummarizedExperiment", file="H5IdComponent"), function(object, file, overwrite) {
    if (!(H5Iget_type(file) %in% c("H5I_FILE", "H5I_GROUP")))
        stop("object must be a file or group")
    write_data_frame(file, "obs", colData(object))
    rdata <- rowData(object)
    if (ncol(rdata) > 0 || !is.null(rownames(rdata)))
        write_data_frame(file, "var", rdata)
    WriteH5AD(assay(object), file, overwrite)
})

#' @importFrom rhdf5 H5Iget_type H5Gcreate H5Gclose
#' @importFrom SingleCellExperiment reducedDims
setMethod("WriteH5AD", c(object="SingleCellExperiment", file="H5IdComponent"), function(object, file, overwrite) {
    if (!(H5Iget_type(file) %in% c("H5I_FILE", "H5I_GROUP")))
        stop("object must be a file or group")

    write_data_frame(file, "var", rowData(object))
    obsm <- reducedDims(object)
    if (length(obsm) > 0) {
        obsmgrp <- H5Gcreate(file, "obsm")
        mapply(function(name, data) {
            if (is.data.frame(data)) {
                rownames(data) <- rownames(colData(object))
                write_data_frame(obsmgrp, name, data)
            } else {
                if (length(dim(data)) == 1)
                    data <- as.vector(data)
                else
                    data <- t(data)
                write_matrix(obsmgrp, name, data)
            }
        }, names(obsm), obsm)
        H5Gclose(obsmgrp)
    }
    WriteH5AD(as(object, "SummarizedExperiment"), file, overwrite)
})

#' @export
setMethod("WriteH5AD", c(object="ANY", file="character"), function(object, file, overwrite) {
    h5 <- open_h5(file)
    WriteH5AD(object, file, overwrite)
    finalize_anndata(h5)
    invisible(NULL)
})

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
setMethod("WriteH5MU", c(object="MultiAssayExperiment", file="character"), function(object, file, overwrite) {
    h5 <- open_h5(file)

    obs <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    write_data_frame(h5, "obs", obs)

    modalities <- names(experiments(object))

    mods <- H5Gcreate(h5, "mod")
    vars <- lapply(modalities, function(mod) {
        mod_group <- H5Gcreate(mods, mod)
        WriteH5AD(object[[mod]], mod_group)
        H5Gclose(mod_group)

        data.frame(row.names = rownames(object[[mod]]))
    })
    H5Gclose(mods)

    var <- do.call(rbind, vars)
    write_data_frame(h5, "var", var)

    finalize_mudata(h5)
    invisible(NULL)
})

#' @importFrom rhdf5 h5writeDataset h5writeAttribute H5Gcreate H5Gclose H5Fget_name H5Iget_name
write_matrix <- function(parent, key, mat) {
    if (is.matrix(mat) || is.vector(mat)) {
        h5writeDataset(mat, parent, key)
    } else if (is(mat, "dgCMatrix") || is(mat, "dgRMatrix") || is(mat, "DelayedArray") && DelayedArray::is_sparse(mat)) {
        if (is(mat, "DelayedArray"))
            mat <- as(mat, "dgRMatrix")

        grp <- H5Gcreate(parent, key)
        h5writeDataset(mat@p, grp, "indptr")
        h5writeDataset(mat@x, grp, "data")
        h5writeAttribute(rev(dim(mat)), grp, "shape")
        h5writeAttribute("0.1.0", grp, "encoding-version", variableLengthString=TRUE, asScalar=TRUE)
        if (is(mat, "dgCMatrix")) {
            h5writeDataset(mat@i, grp, "indices")
            h5writeAttribute("csr_matrix", grp, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
        } else {
            h5writeDataset(mat@j, grp, "indices")
            h5writeAttribute("csc_matrix", grp, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
        }
        H5Gclose(grp)
    } else if (is(mat, "DelayedArray") && requireNamespace("HDF5Array", quietly=TRUE)) {
        writeArrayToMuData(mat, parent, key)
    } else {
        stop(paste0("Writing matrices of type ", class(mat), " is not implemented."))
    }
}

#' @importFrom rhdf5 H5Gcreate H5Gclose h5writeDataset h5writeAttribute h5createAttribute H5Dclose
write_data_frame <- function(parent, key, df) {
    group <- H5Gcreate(parent, key)

    columns <- colnames(df)
    df[["_index"]] <- rownames(df)
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
            dset <- group & colname
            h5writeAttribute(as.integer(is.ordered(df[[colname]])), cat_dset, "ordered", asScalar=TRUE)
            h5writeAttribute(cat_dset, dset, "categories", asScalar=TRUE)
            H5Dclose(cat_dset)
            H5Dclose(dset)
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

.registeredDelayedArrayMethods <- FALSE
#' @importFrom rhdf5 h5write H5Dclose
registerDelayedArrayMethods <- function() {
    if (!.registeredDelayedArrayMethods) {
        haveDelayedArray <- requireNamespace("DelayedArray", quietly=TRUE)
        haveHDF5Array <- requireNamespace("HDF5Array", quietly=TRUE)
        if (!haveDelayedArray || !haveHDF5Array)
            return(FALSE)
        setClass("MuDataFileRealizationSink",
                 contains="HDF5RealizationSink",
                 slots=c(parent="H5IdComponent",
                         datasetname="character"))


        setMethod(DelayedArray::write_block, "MuDataFileRealizationSink", function(sink, viewport, block) {
            if (!is.array(block))
                block <- as.array(block)
            h5write(block, sink@parent, sink@datasetname, start=DelayedArray::start(viewport), count=DelayedArray::width(viewport))
            sink
        })

        .registeredDelayedArrayMethods <<- TRUE
    }
    .registeredDelayedArrayMethods
}

#' @importFrom rhdf5 h5createDataset H5Fget_name H5Iget_name
MuDataFileRealizationSink <- function(dim, type, parent, key, dimnames=NULL, as.sparse=FALSE) {
    chunkdim <- HDF5Array::getHDF5DumpChunkDim(dim)
    h5createDataset(parent, key, dim, storage.mode=type, chunk=chunkdim, level=9, shuffle=FALSE)
    file <- H5Fget_name(parent)
    path <- paste(H5Iget_name(parent), key, sep="/")
    new("MuDataFileRealizationSink", dim=dim, dimnames=dimnames, type=type, as_sparse=as.sparse,
                                     filepath=file, name=path, chunkdim=chunkdim, parent=parent, datasetname=key)
}

writeArrayToMuData <- function(x, parent, key, verbose=NA) {
    if (!registerDelayedArrayMethods())
        stop("The HDF5Array and DelayedArray packages must be installed to save DelayedArrays.")
    as.sparse <- DelayedArray::is_sparse(x)
    sink_dimnames <- dimnames(x)
    sink <- MuDataFileRealizationSink(dim(x), DelayedArray::type(x), parent, key, sink_dimnames, as.sparse)

    verbose <- DelayedArray:::normarg_verbose(verbose)
    sink <- DelayedArray::BLOCK_write_to_sink(sink, x, verbose=verbose)
    invisible(NULL)
}

setGeneric("WriteH5AD", function(object, file, overwrite=TRUE, ...) standardGeneric("WriteH5AD"), signature=c("object", "file"))
setGeneric("WriteH5MU", function(object, file, overwrite=TRUE) standardGeneric("WriteH5MU"), signature=c("object", "file"))

#' @importClassesFrom Matrix Matrix
#' @importClassesFrom DelayedArray DelayedMatrix
setClassUnion("Matrix_OR_DelayedMatrix", c("matrix", "Matrix", "DelayedMatrix"))

#' @importFrom rhdf5 H5Iget_type
setMethod("WriteH5AD", c(object="Matrix_OR_DelayedMatrix", file="H5IdComponent"), function(object, file, overwrite, write_dimnames=TRUE) {
    if (!(H5Iget_type(file) %in% c("H5I_FILE", "H5I_GROUP")))
        stop("object must be a file or group")
    write_matrix(file, "X", object)
    if (write_dimnames) {
        rownames <- rownames(object)
        colnames <- colnames(object)
        if (is.null(rownames))
            rownames <- as.character(1:nrow(object))
        if (is.null(colnames))
            colnames <- as.character(1:ncol(object))
        var <- data.frame(row.names=rownames)
        obs <- data.frame(row.names=colnames)
        write_data_frame(file, "obs", obs)
        write_data_frame(file, "var", var)
    }
    finalize_anndata_internal(file)
})

#' @importFrom rhdf5 H5Iget_type H5Gcreate H5Gclose
#' @importFrom SummarizedExperiment colData assays
setMethod("WriteH5AD", c(object="SummarizedExperiment", file="H5IdComponent"), function(object, file, overwrite) {
    if (!(H5Iget_type(file) %in% c("H5I_FILE", "H5I_GROUP")))
        stop("object must be a file or group")
    write_data_frame(file, "obs", colData(object))
    rdata <- rowData(object)
    if (ncol(rdata) > 0 || !is.null(rownames(rdata)))
        write_data_frame(file, "var", rdata)

    assays <- assays(object)
    nassays <- length(assays)
    if (nassays > 1) {
        layersgrp <- H5Gcreate(file, "layers")
        mapply(function(name, mat) {
            write_matrix(layersgrp, name, mat)
        }, names(assays[2:nassays]), assays[2:nassays])
        H5Gclose(layersgrp)
    }
    WriteH5AD(assays[[1]], file, overwrite, write_dimnames=FALSE)
})

#' @importFrom rhdf5 H5Iget_type H5Gcreate H5Gclose
#' @importFrom SingleCellExperiment rowData colData colPairNames colPair rowPairNames rowPair reducedDims
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

    lapply(list(list(name="obsp", names=colPairNames, getter=colPair), list(name="varp", names=rowPairNames, getter=rowPair)), function(cp) {
        names <- cp$names(object)
        if (length(names) > 0) {
            pairgrp <- H5Gcreate(file, cp$name)
            lapply(names, function(cname) {
                write_matrix(pairgrp, cname, cp$getter(object, cname, asSparse=TRUE))
            })
            H5Gclose(pairgrp)
        }
    })
    WriteH5AD(as(object, "SummarizedExperiment"), file, overwrite)
})

setMethod("WriteH5AD", c(object="ANY", file="H5IdComponent"), function(object, file, overwrite) {
    warning(paste("Objects of class", class(object), "are currently unsupported, skipping..."))
})

#' Save an experiment to an .h5ad file.
#'
#' @param object The object to save.
#' @param file Name of the file to save to.
#' @param overwrite Currently unused.
#'
#' @rdname WriteH5AD
#' @export
setMethod("WriteH5AD", c(object="ANY", file="character"), function(object, file, overwrite) {
    h5 <- open_h5(file)
    WriteH5AD(object, h5, overwrite)
    finalize_anndata(h5)
    invisible(NULL)
})

#' Save a \code{\linkS4class{MultiAssayExperiment}} to an .h5mu file.
#'
#' @param object A \code{\linkS4class{MultiAssayExperiment}}.
#' @param file Name of the file to save to.
#' @param overwrite Currently unused.
#'
#' @rdname WriteH5MU
#'
#' @importFrom rhdf5 H5Gcreate H5Gclose h5writeDataset
#' @importFrom MultiAssayExperiment colData experiments sampleMap
#'
#' @export
setMethod("WriteH5MU", c(object="MultiAssayExperiment", file="character"), function(object, file, overwrite) {
    h5 <- open_h5(file)

    obs <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    write_data_frame(h5, "obs", obs)

    mods <- H5Gcreate(h5, "mod")
    obsmgrp <- H5Gcreate(h5, "obsm")
    obsmapgrp <- H5Gcreate(h5, "obsmap")
    samplemap <- sampleMap(object)
    globalrownames <- rownames(obs)
    vars <- mapply(function(mname, mod) {
        mod_group <- H5Gcreate(mods, mname)
        WriteH5AD(mod, mod_group)
        H5Gclose(mod_group)

        cmap <- samplemap[samplemap$assay == mname,]
        cmaporder <- match(globalrownames, cmap$primary)
        localorder <- match(cmap$colname, colnames(mod))
        obsmap <- sapply(cmaporder, function(o)ifelse(is.na(o), 0L, localorder[o]))
        h5writeDataset(obsmap, obsmapgrp, mname)
        h5writeDataset(as.integer(!is.na(cmaporder)), obsmgrp, mname)

        data.frame(row.names = rownames(mod))
    }, names(object), object)
    H5Gclose(mods)
    H5Gclose(obsmgrp)
    H5Gclose(obsmapgrp)

    var <- do.call(rbind, vars)
    write_data_frame(h5, "var", var)

    finalize_mudata(h5)
    invisible(NULL)
})

#' @importFrom rhdf5 h5writeDataset h5writeAttribute H5Gcreate H5Gclose H5Fget_name H5Iget_name
write_matrix <- function(parent, key, mat) {
    if (is.matrix(mat) || is.vector(mat)) {
        writeDataset(parent, key, mat)
    } else if (is(mat, "dgCMatrix") || is(mat, "dgRMatrix") || is(mat, "DelayedArray") && DelayedArray::is_sparse(mat)) {
        if (is(mat, "DelayedArray"))
            mat <- as(mat, "dgRMatrix")

        grp <- H5Gcreate(parent, key)
        writeDataset(grp, "indptr", mat@p)
        writeDataset(grp, "data", mat@x)
        h5writeAttribute(rev(dim(mat)), grp, "shape")
        h5writeAttribute("0.1.0", grp, "encoding-version", variableLengthString=TRUE, asScalar=TRUE)
        if (is(mat, "dgCMatrix")) {
            writeDataset(grp, "indices", mat@i)
            h5writeAttribute("csr_matrix", grp, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
        } else {
            writeDataset(grp, "indices", mat@j)
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
            writeDataset(group, col, as.integer(v) - 1L)
        } else {
            writeDataset(group, col, v)
        }
    }
    if (length(categories) > 0) {
        catgrp <- H5Gcreate(group, "__categories")

        mapply(function(colname, categories) {
            writeDataset(categories, catgrp, colname, variableLengthString=TRUE)
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

.registeredHDF5ArrayMethods <- FALSE
#' @importFrom rhdf5 h5write H5Dclose
#' @importFrom DelayedArray write_block start width
registerHDF5ArrayMethods <- function() {
    if (!.registeredHDF5ArrayMethods) {
        haveHDF5Array <- requireNamespace("HDF5Array", quietly=TRUE)
        if (!haveHDF5Array)
            return(FALSE)
        setClass("MuDataFileRealizationSink",
                 contains="HDF5RealizationSink",
                 slots=c(parent="H5IdComponent",
                         datasetname="character"))


        setMethod(write_block, "MuDataFileRealizationSink", function(sink, viewport, block) {
            if (!is.array(block))
                block <- as.array(block)
            h5write(block, sink@parent, sink@datasetname, start=start(viewport), count=width(viewport))
            sink
        })

        .registeredHDF5ArrayMethods <<- TRUE
    }
    .registeredHDF5ArrayMethods
}

#' @importFrom rhdf5 h5createDataset H5Fget_name H5Iget_name
MuDataFileRealizationSink <- function(dim, type, parent, key, dimnames=NULL, as.sparse=FALSE) {
    chunkdim <- HDF5Array::getHDF5DumpChunkDim(dim)
    h5createDataset(parent, key, dim, storage.mode=type, chunk=chunkdim, level=9, shuffle=TRUE)
    file <- H5Fget_name(parent)
    path <- paste(H5Iget_name(parent), key, sep="/")
    new("MuDataFileRealizationSink", dim=dim, dimnames=dimnames, type=type, as_sparse=as.sparse,
                                     filepath=file, name=path, chunkdim=chunkdim, parent=parent, datasetname=key)
}

#' @importFrom DelayedArray is_sparse type BLOCK_write_to_sink
writeArrayToMuData <- function(x, parent, key, verbose=NA) {
    if (!registerHDF5ArrayMethods())
        stop("The HDF5Array packages must be installed to save DelayedArrays.")
    as.sparse <- is_sparse(x)
    sink_dimnames <- dimnames(x)
    sink <- MuDataFileRealizationSink(dim(x), type(x), parent, key, sink_dimnames, as.sparse)

    verbose <- DelayedArray:::normarg_verbose(verbose)
    sink <- BLOCK_write_to_sink(sink, x, verbose=verbose)
    invisible(NULL)
}

# this is a straight port of the h5py guess_chunk function
chunk_base <- 16 * 1024
chunk_min <- 8 * 1024
chunk_max <- 1024 * 1024
guess_chunk <- function(shape, storage.mode) {
    nbytes <- switch(storage.mode, double=8, integer=4, logical=1, character=8)
    ndims <- length(shape)

    dset_size <- prod(shape) * nbytes
    target_size <- chunk_base * (2^log10(dset_size / 1024^2))
    if (target_size > chunk_max)
        target_size <- chunk_max
    else if (target_size < chunk_min)
        target_size <- chunk_min

   idx <- 0
   while(TRUE) {
       chunk_bytes <- prod(shape) * nbytes
       if ((chunk_bytes < target_size || abs(chunk_bytes - target_size) / target_size < 0.5) && chunk_bytes < chunk_max)
           break
       if (prod(shape) == 1)
           break

       shape[idx %% ndims + 1] <- ceiling(shape[idx %% ndims + 1] / 2)
       idx <- idx + 1
   }
   as.integer(shape)
}

writeDataset <- function(parent, key, data) {
    shape <- dim(data)
    if (is.null(shape))
        shape <- length(data)
    chunksize <- guess_chunk(shape, storage.mode(data))
    h5createDataset(parent, key, shape, storage.mode=storage.mode(data), chunk=chunksize, level=9, shuffle=TRUE)
    h5writeDataset(data, parent, key)
}

#' @importClassesFrom Matrix Matrix
#' @importClassesFrom DelayedArray DelayedMatrix
setClassUnion("Matrix_OR_DelayedMatrix", c("matrix", "Matrix", "DelayedMatrix"))

#' @importClassesFrom S4Vectors DataFrame
setClassUnion("data.frame_OR_DataFrame", c("data.frame", "DataFrame"))

#' Save an experiment to an .h5ad file.
#'
#' Note that NA values are not supported by HDF5, and therefore by h5ad. The behavior of this
#' function if NAs are present is undefined.
#'
#' @param object The object to save.
#' @param file Name of the file to save to.
#' @param overwrite Currently unused.
#'
#' @returns NULL, invisibly
#'
#' @examples
#' data(miniACC, package="MultiAssayExperiment")
#' writeH5AD(miniACC[[1]], "miniacc.h5ad")
#'
#' @importFrom rhdf5 H5Iget_type H5Gcreate H5Gclose
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors metadata
#' @importFrom methods hasMethod as
#' @importFrom SingleCellExperiment altExps rowData colData colPairNames colPair rowPairNames rowPair reducedDims
#'
#' @export
writeH5AD <- function(object, file, overwrite) {
    need_finalize <- FALSE
    written_object <- FALSE
    if (is.character(file)) {
        file <- open_h5(file)
        need_finalize <- TRUE
    } else if (!(H5Iget_type(file) %in% c("H5I_FILE", "H5I_GROUP")))
        stop("file must be a character, file or group")

    cls <- class(object)
    if (is(object, "RangedSummarizedExperiment") && !is(object, "SingleCellExperiment")) {
        warning("Ranged data is currently unsupported. Coercing to SummarizedExperiment...")
        object <- as(object, "SummarizedExperiment")
    }
    if (is(object, "SingleCellExperiment")) {
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
        object <- as(object, "SummarizedExperiment")
        written_object <- TRUE
    }
    if (is(object, "SummarizedExperiment")) {
        write_data_frame(file, "obs", colData(object))
        rdata <- rowData(object)
        if (ncol(rdata) > 0 || !is.null(rownames(rdata)))
            write_data_frame(file, "var", rdata)

        assays <- assays(object)
        nassays <- length(assays)
        write_matrix(file, "X", assays[[1]])
        if (nassays > 1) {
            layersgrp <- H5Gcreate(file, "layers")
            mapply(function(name, mat) {
                write_matrix(layersgrp, name, mat)
            }, names(assays[2:nassays]), assays[2:nassays])
            H5Gclose(layersgrp)
        }

        writeList(file, "uns", metadata(object))

        if (hasMethod("altExps", class(object))) {
            naltexps <- length(altExps(object))
            if (naltexps > 1) {
                warning("Alternative experiments are currently unsupported. Construct a MultiAssayExperiment object or write them as individual H5AD files.")
            }
        }
        written_object <- TRUE
    }
    if (is(object, "Matrix_OR_DelayedMatrix")) {
        write_matrix(file, "X", object)
        rownames <- rownames(object)
        colnames <- colnames(object)
        if (is.null(rownames))
            rownames <- as.character(seq_len(nrow(object)))
        if (is.null(colnames))
            colnames <- as.character(seq_len(ncol(object)))
        var <- data.frame(row.names=rownames)
        obs <- data.frame(row.names=colnames)
        write_data_frame(file, "obs", obs)
        write_data_frame(file, "var", var)
        written_object <- TRUE
    }

    if (!written_object) {
        warning("Objects of class ", class(object), " are currently unsupported, skipping...")
    } else {
        write_object_class(file, cls[1])
        finalize_anndata_internal(file)
        if (need_finalize)
            finalize_anndata(file)
    }
    invisible(NULL)
}

#' Save a \code{\linkS4class{MultiAssayExperiment}} to an .h5mu file.
#'
#' Note that NA values are not supported by HDF5, and therefore by h5mu. The behavior of this
#' function if NAs are present is undefined.
#'
#' @param object A \code{\linkS4class{MultiAssayExperiment}}.
#' @param file Name of the file to save to.
#' @param overwrite Currently unused.
#'
#' @returns NULL, invisibly
#'
#' @examples
#' data(miniACC, package="MultiAssayExperiment")
#' writeH5MU(miniACC, "miniacc.h5mu")
#'
#' @importFrom rhdf5 H5Gcreate H5Gclose
#' @importFrom MultiAssayExperiment colData experiments sampleMap metadata
#'
#' @export
writeH5MU <- function(object, file, overwrite) {
    if (!is(object, "MultiAssayExperiment"))
        stop("Only MultiAssayExperiment objects are currently supported.")
    if (is.character(file)) {
        h5 <- open_h5(file)
        need_finalize <- TRUE
    } else
        stop("file must be a character")

    obs <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    write_data_frame(h5, "obs", obs)

    mods <- H5Gcreate(h5, "mod")
    obsmgrp <- H5Gcreate(h5, "obsm")
    obsmapgrp <- H5Gcreate(h5, "obsmap")
    samplemap <- sampleMap(object)
    globalrownames <- rownames(obs)
    vars <- mapply(function(mname, mod) {
        mod_group <- H5Gcreate(mods, mname)
        writeH5AD(mod, mod_group)
        H5Gclose(mod_group)

        cmap <- samplemap[samplemap$assay == mname,]
        cmaporder <- match(globalrownames, cmap$primary)
        localorder <- match(cmap$colname, colnames(mod))
        obsmap <- vapply(cmaporder, function(o)ifelse(is.na(o), 0L, localorder[o]), 0L)
        writeDataset(obsmapgrp, mname, obsmap)
        writeDataset(obsmgrp, mname, as.integer(!is.na(cmaporder)))

        data.frame(row.names = rownames(mod))
    }, names(object), object)
    writeAttribute(mods, "order", names(object))
    H5Gclose(mods)
    H5Gclose(obsmgrp)
    H5Gclose(obsmapgrp)

    var <- do.call(rbind, vars)
    write_data_frame(h5, "var", var)

    writeList(h5, "uns", metadata(object))

    write_object_class(h5, class(object)[1])
    finalize_mudata(h5)
    invisible(NULL)
}

#' @importFrom rhdf5 H5Gcreate H5Gclose
#' @importFrom methods is as
write_matrix <- function(parent, key, mat) {
    if (is(mat, "dgeMatrix"))
        mat <- as.matrix(mat)
    if (is.matrix(mat) || is.vector(mat) || is.array(mat) || is.numeric(mat) || is.integer(mat) || is.logical(mat) || is.character(mat)) { # is.vector returns false for vectors with attributes
        isscalar <- length(mat) == 1 & !is.null(attr(mat, "encoding-scalar"))
        hasna <- anyNA(mat)
        if (hasna && is.double(mat)) {
            # FIXME: extend anndata spec to handle double NAs?
            mat[is.na(mat)] <- NaN
            hasna <- FALSE
        }

        if (isscalar || !hasna) {
            writeDataset(parent, key, mat, scalar=isscalar)
            dset <- h5autoclose(parent & key)

            if (!isscalar)
                writeAttribute(dset, "encoding-type", ifelse(is.character(mat), "string-array", "array"))
            else
                writeAttribute(dset, "encoding-type", ifelse(is.character(mat), "string", "numeric-scalar"))
            writeAttribute(dset, "encoding-version", "0.2.0")
        } else {
            grp <- H5Gcreate(parent, key)
            write_matrix(grp, "values", mat)
            write_matrix(grp, "mask", is.na(mat))
            writeAttribute(grp, "encoding-type", ifelse(is.logical(mat), "nullable-boolean", "nullable-integer"))
            writeAttribute(grp, "encoding-version", "0.1.0")
            H5Gclose(grp)
        }
    } else if (is.factor(mat)) {
        grp <- H5Gcreate(parent, key)
        codes <- as.integer(mat)
        codes[is.na(mat)] <- 0L
        write_matrix(grp, "codes", codes - 1L)
        write_matrix(grp, "categories", levels(mat))
        writeAttribute(grp, "ordered", is.ordered(mat))
        writeAttribute(grp, "encoding-type", "categorical")
        writeAttribute(grp, "encoding-version", "0.2.0")
        H5Gclose(grp)
    } else if (is(mat, "dgCMatrix") || is(mat, "dgRMatrix") || is(mat, "DelayedArray") && DelayedArray::is_sparse(mat)) {
        if (is(mat, "DelayedArray"))
            mat <- as(mat, "RsparseMatrix")

        grp <- H5Gcreate(parent, key)
        writeDataset(grp, "indptr", mat@p)
        writeDataset(grp, "data", mat@x)
        writeAttribute(grp, "shape", rev(dim(mat)))
        writeAttribute(grp, "encoding-version", "0.1.0")
        if (is(mat, "dgCMatrix")) {
            writeDataset(grp, "indices", mat@i)
            writeAttribute(grp, "encoding-type", "csr_matrix")
        } else {
            writeDataset(grp, "indices", mat@j)
            writeAttribute(grp, "encoding-type", "csc_matrix")
        }
        H5Gclose(grp)
    } else if (is(mat, "DelayedArray") && requireNamespace("HDF5Array", quietly=TRUE)) {
        writeArrayToMuData(mat, parent, key)
    } else {
        stop("Writing matrices of type ", class(mat), " is not implemented.")
    }
}

#' @importFrom rhdf5 H5Gcreate H5Gclose h5createAttribute
write_data_frame <- function(parent, key, df) {
    group <- H5Gcreate(parent, key)

    columns <- colnames(df)
    df[["_index"]] <- rownames(df)
    for (col in colnames(df)) {
        write_matrix(group, col, df[[col]])
    }

    # Write attributes
    writeAttribute(group, "_index", "_index")
    writeAttribute(group, "encoding-type", "dataframe")
    writeAttribute(group, "encoding-version", "0.2.0")
    if (length(columns) > 0) {
        writeAttribute(group, "column-order", columns)
    } else {
        # When there are no columns, null buffer can't be written to a file.
        h5createAttribute(group, "column-order", dims=0)
    }
    H5Gclose(group)
    invisible(NULL)
}

.registeredHDF5ArrayMethods <- new.env()
.registeredHDF5ArrayMethods$registered <- FALSE
#' @importFrom rhdf5 h5write H5Dclose
#' @importFrom DelayedArray write_block start width
#' @importFrom methods setClass
registerHDF5ArrayMethods <- function() {
    if (!.registeredHDF5ArrayMethods$registered) {
        haveHDF5Array <- requireNamespace("HDF5Array", quietly=TRUE)
        if (!haveHDF5Array)
            return(FALSE)

        setClass("MuDataFileRealizationSink",
                 contains="HDF5RealizationSink",
                 slots=c(parent="H5IdComponent",
                         datasetname="character"),
                 where=.registeredHDF5ArrayMethods)


        setMethod(write_block, "MuDataFileRealizationSink", function(sink, viewport, block) {
            if (!is.array(block))
                block <- as.array(block)
            h5write(block, sink@parent, sink@datasetname, start=start(viewport), count=width(viewport))
            sink
        }, where=.registeredHDF5ArrayMethods)

        .registeredHDF5ArrayMethods$registered <- TRUE
    }
    .registeredHDF5ArrayMethods$registered
}

#' @importFrom rhdf5 h5createDataset H5Fget_name H5Iget_name
#' @importFrom methods new
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

#' @importFrom rhdf5 h5createDataset h5writeDataset
writeDataset <- function(parent, key, data, scalar=FALSE) {
    shape <- dim(data)
    if (is.null(shape))
        shape <- length(data)
    if (length(shape) == 1 && shape == 1 && scalar) {
        shape <- chunksize <- NULL
        level <- 0
    } else {
        chunksize <- guess_chunk(shape, storage.mode(data))
        level <- 9
    }
    h5createDataset(parent, key, shape, storage.mode=storage.mode(data), chunk=chunksize, level=level, shuffle=TRUE, encoding="UTF-8")
    h5writeDataset(data, parent, key, variableLengthString=TRUE, encoding="UTF-8")
}

#' @importFrom rhdf5 h5writeAttribute
writeAttribute <- function(obj, name, value) {
    if (is.logical(value))
        value <- as.integer(value) # rhdf5 hasn't implemented logical attributes yet
    args <- list(attr=value, h5obj=obj, name=name)
    if (is.character(value))
        args$variableLengthString <- TRUE
    if (length(value) == 1)
        args$asScalar <- TRUE
    do.call(h5writeAttribute, args)
}

#' @importFrom rhdf5 H5Gcreate H5Gclose
#' @importFrom methods slotNames slot
write_elem <- function(parent, key, data) {
    if (is(data, "list_OR_List"))
        writeList(parent, key, data)
    else if (is(data, "Matrix_OR_DelayedMatrix") || is(data, "vector_OR_Vector"))
        write_matrix(parent, key, data)
    else if (is(data, "data.frame_OR_DataFrame"))
        write_data_frame(parent, key, data)
    else if (isS4(data)) {
        grp <- H5Gcreate(parent, key)
        for (slotnm in slotNames(data)) {
            write_elem(grp, slotnm, slot(data, slotnm))
        }
        writeattribute(grp, "encoding-type", "dict")
        writeAttribute(grp, "encoding-version", "0.1.0")
        H5Gclose(grp)
    } else
        warning(paste("Cannot write object of class", class(data), "skipping..."))
}

#' @importFrom rhdf5 H5Gcreate H5Gclose
writeList <- function(parent, key, data) {
    if (length(data) > 0) {
        nms <- names(data)
        if (is.null(nms))
            nms <- as.character(seq_along(data))
        grp <- H5Gcreate(parent, key)
        mapply(function(name, data) {
            write_elem(grp, name, data)
        }, nms, data)
        writeAttribute(grp, "encoding-type", "dict")
        writeAttribute(grp, "encoding-version", "0.1.0")
        H5Gclose(grp)
    }
}

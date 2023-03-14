#' @importFrom rhdf5 H5Aexists H5Aopen H5Aread H5Aclose H5Dread H5Dclose H5Rdereference
#' @importMethodsFrom rhdf5 &
#' @importFrom methods is
read_dataframe <- function(group) {
    indexcol <- "_index"
    if (H5Aexists(group, "_index")) {
        indexattr <- H5Aopen(group, "_index")
        indexcol <- H5Aread(indexattr)
        H5Aclose(indexattr)
    }

    orderedattr <- H5Aopen(group, "column-order")
    columnorder <- H5Aread(orderedattr)
    H5Aclose(orderedattr)

    col_list <- lapply(columnorder, function(name) {
        col <- group & name
        values <- read_attribute(col)
        if (H5Aexists(col, "categories")) {
            attr <- H5Aopen(col, "categories")
            labels <- H5Aread(attr)
            if (!is(labels, "H5Ref")) {
                warning("found categories attribute for column ",
                        name, ", but it is not a reference")
            } else {
                labels <- H5Rdereference(labels, h5loc=col)
                labels_items <- H5Dread(labels)
                n_labels <- length(unique(values))
                if (length(labels_items) > n_labels) {
                    labels_items <- labels_items[seq_len(n_labels)]
                }
                values <- factor(as.integer(values), labels=labels_items)
                H5Dclose(labels)
            }
            H5Aclose(attr)
        }
        H5Dclose(col)
        values
    })
    names(col_list) <- columnorder
    index <- group & indexcol
    col_list[["row.names"]] <- H5Dread(index)
    H5Dclose(index)
    do.call(data.frame, args=col_list)
}

#' @importFrom rhdf5 H5Dread H5Aexists H5Aopen H5Aread H5Aclose
read_dataframe_legacy <- function(dataset) {
    table <- H5Dread(dataset)

    indexcol <- "_index"
    if (H5Aexists(dataset, "_index")) {
        index <- H5Aopen(dataset, "_index")
        indexcol <- H5Aread(index)
        H5Aclose(index)
    }

    if (indexcol %in% colnames(table)) {
        rownames(table) <- table[,indexcol,drop=TRUE]
        table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }
    table
}

#' @importFrom rhdf5 H5Iget_type H5Aexists H5Aopen H5Aread H5Aclose
read_with_index <- function(dataset) {
    cls <- H5Iget_type(dataset)
    if (cls == "H5I_GROUP" && H5Aexists(dataset, "encoding-type")) {
        encattr <- H5Aopen(dataset, "encoding-type")
        encoding <- H5Aread(encattr, "encoding-type")
        H5Aclose(encattr)
        if (encoding != "dataframe") {
            warning("Unknown encoding ", encoding,
                    " when attempting to read data frame")
            return(data.frame())
        }
        read_dataframe(dataset)
    } else {
        read_dataframe_legacy(dataset)
    }
}

#' @importFrom rhdf5 H5Dread H5Dclose H5Aopen H5Aread H5Aclose
#' @importMethodsFrom rhdf5 &
read_sparse_matrix <- function(group, encoding, backed=FALSE) {
    indices <- group & "indices"
    indptr <- group & "indptr"
    data <- group & "data"
    i <- as.vector(H5Dread(indices))
    p <- as.vector(H5Dread(indptr))
    x <- as.vector(H5Dread(data))
    H5Dclose(indices)
    H5Dclose(indptr)
    H5Dclose(data)
    shapeattr <- H5Aopen(group, "shape")
    shape <- H5Aread(shapeattr)
    H5Aclose(shapeattr)
    if (encoding == "csr_matrix") {
        sparseMatrix(i=i, p=p, x=x, dims=rev(shape), repr="C", index1=FALSE)
    } else {
        t(sparseMatrix(i=i, p=p, x=x, dims=shape, repr="C", index1=FALSE))
    }
}

#' @importFrom rhdf5 H5Iget_type H5Iget_name H5Aexists H5Aopen H5Aread H5Aclose H5Dread H5Dclose H5Gclose H5Fget_name
#' @importFrom methods new
read_matrix <- function(dataset, backed=FALSE) {
    if (backed) {
        have_delayedarray <- requireNamespace("HDF5Array", quietly=TRUE)
        if (!have_delayedarray) {
            stop("Could not load the HDF5Array package. HDF5Array is required for backed matrices.")
        }
    }

    if (H5Iget_type(dataset) == "H5I_GROUP" && H5Aexists(dataset, "encoding-type")) {
        encattr <- H5Aopen(dataset, "encoding-type")
        encoding <- H5Aread(encattr)
        H5Aclose(encattr)
        if (encoding %in% c("csr_matrix", "csc_matrix")) {
            if (backed) {
                cls <- ifelse(encoding == "csr_matrix", "CSC_H5ADMatrixSeed", "CSR_H5ADMatrixSeed")
                seed <- HDF5Array::H5SparseMatrixSeed
            } else {
                return(read_sparse_matrix(dataset, encoding))
            }
        } else {
            warning("Unknown encoding ", encoding, "when attempting to read matrix")
            return(matrix())
        }
    } else {
        if (backed) {
            cls <- "Dense_H5ADMatrixSeed"
            seed <- HDF5Array::HDF5ArraySeed
        } else {
            return(H5Dread(dataset))
        }
    }
    if (backed) {
        file <- H5Fget_name(dataset)
        name <- H5Iget_name(dataset)
        suppressWarnings({
            seed <- seed(file, name)
            seed <- new(cls, seed)
            HDF5Array::H5ADMatrix(seed)
        })
    }
}

#' @importFrom rhdf5 H5Aopen H5Aread H5Aclose H5Aexists h5ls
#' @importFrom stats setNames
read_group <- function(group, read_uns=FALSE) {
    if (H5Aexists(group, "encoding-type")) {
        encattr <- H5Aopen(group, "encoding-type")
        encoding <- H5Aread(encattr)
        H5Aclose(encattr)

        if (encoding == "dataframe") {
            return(read_dataframe(group))
        } else if (endsWith(encoding, "matrix")) {
            return(read_sparse_matrix(group))
        } else if (!encoding %in% c("dict")) {
            # "dict" encoding is just a regular group of attributes and should 
            # be read as though there is no defined encoding
            warning("Unknown encoding ", encoding)
            if (!read_uns)
                return(invisible(NULL))
        }
    }

    objects <- h5ls(group, recursive=FALSE, datasetinfo=FALSE)$name
    lapply(setNames(nm=objects), function(x)read_attribute(h5autoclose(group & x)))
}

#' @importFrom rhdf5 H5Iget_type H5Dread
read_attribute <- function(attr) {
    if (H5Iget_type(attr) == "H5I_GROUP")
        read_group(attr)
    else {
        values <- H5Dread(attr)
        # h5py saves boolean arrays as HDF5 enums
        if (is.factor(values) && all(levels(values) == c("FALSE", "TRUE"))) {
            values <- as.logical(values)
        }
        values
    }
}

#' @importFrom stats setNames
#' @importFrom rhdf5 h5ls
#' @importMethodsFrom rhdf5 &
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment colPair<- rowPair<-
#' @importFrom SummarizedExperiment SummarizedExperiment metadata<-
#' @importFrom methods is
read_modality <- function(view, backed=FALSE) {
    X <- read_matrix(h5autoclose(view & "X"), backed=backed)
    var <- read_with_index(h5autoclose(view & "var"))
    obs <- read_with_index(h5autoclose(view  & "obs"))
    rownames(X) <- rownames(var)
    colnames(X) <- rownames(obs)

    viewnames <- h5ls(view, recursive=FALSE)$name
    layers <- list()
    if ("layers" %in% viewnames) {
        layers <- lapply(setNames(nm=h5ls(h5autoclose(view & "layers"), recursive=FALSE)$name), function(layer) {
            read_matrix(h5autoclose(view & paste("layers", layer, sep="/")), backed=backed)
        })
    }

    args <- list(assays=c(list(X), layers), rowData=var, colData=obs)

    if ("obsm" %in% viewnames) {
        obsmnames <- h5ls(h5autoclose(view & "obsm"), recursive=FALSE)$name
        obsm <- lapply(obsmnames, function(space) {
            elem <- read_attribute(h5autoclose(view & paste("obsm", space, sep="/")))
            if (!is.data.frame(elem) && length(dim(elem)) > 1)
                elem <- t(elem)
            rownames(elem) <- rownames(obs)
            elem
        })
        names(obsm) <- obsmnames
        args$reducedDims <- obsm
    }

    if (H5Aexists(view, "origin-class")) {
        originattr <- H5Aopen(view, "origin-class")
        objectclass <- H5Aread(originattr)
        H5Aclose(originattr)

	if (objectclass == "SummarizedExperiment") {
		se <- do.call(SummarizedExperiment, args)
	} else {
		if (objectclass != "SingleCellExperiment") {
            		message("Reading as SingleCellExperiment where the original object class is ", objectclass)
		}
		se <- do.call(SingleCellExperiment, args)
	}
    } else {
	se <- do.call(SingleCellExperiment, args)
    }

    for (cp in list(list(name="obsp", setter=`colPair<-`), list(name="varp", setter=`rowPair<-`))) {
        if (cp$name %in% viewnames) {
            names <- h5ls(h5autoclose(view & cp$name), recursive=FALSE)$name
            for (name in names) {
                cpair <- read_matrix(h5autoclose(view & paste(cp$name, name, sep="/")))
                if (!is(cpair, "dsparseMatrix")) {
                    warning("Pairwise ", cp$name, " matrix ", name, " is not a sparse matrix. Only sparse matrices are currently supported, skipping...")
                } else {
                    se <- cp$setter(se, name, value=cpair)
                }
            }
        }
    }

    if ("uns" %in% viewnames)
        metadata(se) <- read_group(h5autoclose(view & "uns"))
    se
}


#' Read an .h5ad file and create a \code{\linkS4class{SingleCellExperiment}}.
#'
#' In file-backed mode, the main \code{X} matrix is not read into memory,
#' but references the HDF5 file and its required parts are read on demand.
#' This requires the HDF5Array package to be installed.
#'
#' @param file Path to the .h5ad file.
#' @param backed Whether to use file-backed mode.
#'
#' @return A \code{\linkS4class{SingleCellExperiment}}.
#'
#' @examples
#' data(miniACC, package="MultiAssayExperiment")
#' writeH5AD(miniACC[[1]], "miniacc.h5ad")
#' sce <- readH5AD("miniacc.h5ad")
#'
#' @importFrom rhdf5 H5Fclose
#' @export
readH5AD <- function(file, backed=FALSE) {
    h5 <- H5Fopen(file, flags="H5F_ACC_RDONLY", native=FALSE)
    res <- read_modality(h5, backed)
    H5Fclose(h5)
    res
}

#' Read an .h5mu file and create a \code{\link{MultiAssayExperiment}}.
#'
#' In file-backed mode, the main \code{X} matrices are not read into memory,
#' but reference the HDF5 file and their required parts are read on demand.
#' This requires the HDF5Array package to be installed.
#'
#' @param file Path to the .h5mu file.
#' @param backed Whether to use file-backed mode.
#'
#' @return A \code{\linkS4class{MultiAssayExperiment}}
#'
#' @examples
#' data(miniACC, package="MultiAssayExperiment")
#' writeH5MU(miniACC, "miniacc.h5mu")
#' mae <- readH5MU("miniacc.h5mu")
#'
#' @importFrom stats setNames
#' @importFrom rhdf5 h5ls H5Fclose H5Lexists H5Dread H5Dclose
#' @importMethodsFrom rhdf5 &
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
readH5MU <- function(file, backed=FALSE) {
    # Connect to the the file
    h5 <- open_and_check_mudata(file)

    # Check all the assays are written
    assays <- setNames(nm=h5ls(h5autoclose(h5 & "mod"), recursive=FALSE)$name)

    # Read modality order if available and matches available assays
    mod_order <- check_mod_order(h5)

    # Create global colData
    metadata <- read_with_index(h5autoclose(h5 & "obs"))

    # Create an experiments list
    modalities <- lapply(assays, function(mod) {
        read_modality(h5autoclose(h5 & paste("mod", mod, sep="/")), backed)
    })
    modalities <- modalities[mod_order]

    args <- list(experiments=modalities, colData=metadata)

    # create a sampleMap. This is needed for a round-trip MAE -> .h5mu -> MAE
    # if colData(MAE) has different row names than the experiments
    if (H5Lexists(h5, "obsmap")) {
        samplemaps <- lapply(assays, function(mod) {
            cmapdset <- h5 & paste("obsmap", mod, sep="/")
            cmap <- H5Dread(cmapdset)
            H5Dclose(cmapdset)

            idx <- which(cmap > 0)
            data.frame(assay=mod, primary=rownames(metadata)[idx], colname=colnames(modalities[[mod]])[cmap[idx]])
        })
        sampleMap <- do.call(rbind, samplemaps)
        rownames(sampleMap) <- NULL
        args$sampleMap <- sampleMap
    }

    if (H5Lexists(h5, "uns"))
        args$metadata <- read_group(h5autoclose(h5 & "uns"))

    # Close the connection
    H5Fclose(h5)

    # Create a MAE object
    do.call(MultiAssayExperiment, args)
}


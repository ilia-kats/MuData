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
        values <- H5Dread(col)
        # h5py saves boolean arrays as HDF5 enums
        if (is.factor(values) && levels(values) == c("FALSE", "TRUE")) {
            values <- as.logical(values)
        }
        if (H5Aexists(col, "categories")) {
            attr <- H5Aopen(col, "categories")
            labels <- H5Aread(attr)
            if (!is(labels, "H5Ref")) {
                warning("found categories attribute for column ", name, ", but it is not a reference")
            } else {
                labels <- H5Rdereference(labels, h5loc=col)
                values <- factor(as.integer(values), labels=H5Dread(labels))
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
            warning("Unknown encoding ", encoding, " when attempting to read data frame")
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
            warning("Could not load the HDF5Array package. HDF5Array is required for backed matrices. Loading matrix into memory...")
            backed <- FALSE
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

#' @importFrom rhdf5 H5Aopen H5Aread H5Aclose
read_group <- function(group) {
    encattr <- H5Aopen(group, "encoding-type")
    encoding <- H5Aread(encattr)
    H5Aclose(encattr)
    if (encoding == "dataframe") {
        read_dataframe(group)
    } else if (endsWith(encoding, "matrix")) {
        read_sparse_matrix(group)
    } else {
        warning("Unknown encoding ", encoding)
        invisible(NULL)
    }
}

#' @importFrom rhdf5 H5Iget_type H5Dread
read_attribute <- function(attr) {
    if (H5Iget_type(attr) == "H5I_GROUP")
        read_group(attr)
    else
        H5Dread(attr)
}

#' @importFrom stats setNames
#' @importFrom rhdf5 h5ls
#' @importMethodsFrom rhdf5 &
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment colPair<- rowPair<-
#' @importFrom SummarizedExperiment SummarizedExperiment
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

    se <- do.call(SingleCellExperiment, args)

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
#' WriteH5AD(miniACC[[1]], "miniacc.h5ad")
#' sce <- ReadH5AD("miniacc.h5ad")
#'
#' @importFrom rhdf5 H5Fclose
#' @export
ReadH5AD <- function(file, backed=FALSE) {
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
#' WriteH5MU(miniACC, "miniacc.h5mu")
#' mae <- ReadH5MU("miniacc.h5mu")
#'
#' @importFrom stats setNames
#' @importFrom rhdf5 h5ls H5Fclose H5Lexists H5Dread H5Dclose
#' @importMethodsFrom rhdf5 &
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
ReadH5MU <- function(file, backed=FALSE) {
    # Connect to the the file
    h5 <- open_and_check_mudata(file)

    # Check all the assays are written
    assays <- setNames(nm=h5ls(h5autoclose(h5 & "mod"), recursive=FALSE)$name)

    # Create global colData
    metadata <- read_with_index(h5autoclose(h5 & "obs"))

    # Create an experiments list
    modalities <- lapply(assays, function(mod) {
        read_modality(h5autoclose(h5 & paste("mod", mod, sep="/")), backed)
    })

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

    # Close the connection
    H5Fclose(h5)

    # Create a MAE object
    do.call(MultiAssayExperiment, args)
}


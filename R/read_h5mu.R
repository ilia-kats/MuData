#' @importFrom rhdf5 H5Aexists H5Aopen H5Aread H5Aclose H5Dread H5Dclose H5Rdereference
#' @importFrom S4Vectors DataFrame
#' @importMethodsFrom rhdf5 &
#' @importFrom methods is
read_dataframe_020 <- function(group, encoding, strict=TRUE) {
    if (strict)
        check_encodingversion(group, encoding, "0.2.0")
    indexattr <- H5Aopen(group, "_index")
    indexcol <- H5Aread(indexattr)
    H5Aclose(indexattr)

    orderedattr <- H5Aopen(group, "column-order")
    columnorder <- H5Aread(orderedattr)
    H5Aclose(orderedattr)

    col_list <- lapply(columnorder, function(name) {
        read_attribute(h5autoclose(group & name))
    })
    names(col_list) <- columnorder
    index <- h5autoclose(group & indexcol)
    col_list[["row.names"]] <- H5Dread(index)
    do.call(DataFrame, args=col_list)
}

#' @importFrom rhdf5 H5Aexists H5Aopen H5Aread H5Aclose H5Dread H5Dclose H5Rdereference
#' @importFrom S4Vectors DataFrame
#' @importMethodsFrom rhdf5 &
#' @importFrom methods is
read_dataframe_010 <- function(group, encoding, strict=TRUE) {
    if (strict)
        check_encodingversion(group, encoding, "0.1.0")
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
                values <- convert_categoricals(values, labels_items)
                H5Dclose(labels)
            }
            H5Aclose(attr)
        }
        H5Dclose(col)
        values
    })
    names(col_list) <- columnorder
    col_list[["row.names"]] <- H5Dread(h5autoclose(group & indexcol))
    do.call(DataFrame, args=col_list)
}

read_dataframe <- function(group, encoding) {
    version <- check_encodingversion(group, encoding, c("0.1.0", "0.2.0"))
    if (version == "0.1.0")
        read_dataframe_010(group, encoding, strict=FALSE)
    else
        read_dataframe_020(group, encoding, strict=FALSE)
}

#' Helper function to convert values + labels into factors
#'
#' @description  A helper function to convert categories into factors.
#'  Assumptions:
#'      - values correspond to the zero indexed categories
#'          (i.e. value 0 is the first category)
#'      - NA are encoded with a value -1
#'  Categories not uses will be dropped.
#'
#' @param values Vector of integer level numbers (zero indexed). -1 indicate NA
#' @param categories Labels for level numbers (zero indexed).
#'
#' @returns factor with categorical values
#'
#' @keywords internal
#' @noRd
convert_categoricals <- function(values, categories) {
    # The levels are 0 indexed integers
    levels <- seq_len(length(categories))-1
    value_factor <- factor(as.integer(values), levels, labels=categories)
    # Drop unused levels
    droplevels(value_factor)
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

#' @importFrom rhdf5 H5Dread H5Aopen H5Aread H5Aclose
#' @importMethodsFrom rhdf5 &
read_sparse_matrix <- function(group, encoding, backed=FALSE) {
    check_encodingversion(group, encoding, "0.1.0")
    indices <- h5autoclose(group & "indices")
    indptr <- h5autoclose(group & "indptr")
    data <- h5autoclose(group & "data")
    i <- as.vector(H5Dread(indices))
    p <- as.vector(H5Dread(indptr))
    x <- as.vector(H5Dread(data))
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
            dset <- H5Dread(dataset)
            if (length(dim(dset)) == 1)
                dset <- as.vector(dset)
            else if (length(dim(dset)) == 2)
                dset <- as.matrix(dset)
            return(dset)
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

#' @importFrom rhdf5 h5ls
#' @importMethodsFrom rhdf5 &
#' @importFrom stats setNames
read_dict <- function(group, encoding, strict=TRUE) {
    if (strict)
        check_encodingversion(group, encoding, "0.1.0")

    objects <- h5ls(group, recursive=FALSE, datasetinfo=FALSE)$name
    lapply(setNames(nm=objects), function(x)read_attribute(h5autoclose(group & x)))
}

#' @importFrom rhdf5 H5Dread
read_array <- function(attr, encoding, strict=TRUE) {
    if (strict)
        check_encodingversion(attr, encoding, "0.2.0")

    ret <- H5Dread(attr)
    if (length(dim(ret)) == 1)
        ret <- as.vector(ret)
    else if (length(dim(ret)) == 2)
        ret <- as.matrix(ret)
    # h5py saves boolean arrays as HDF5 enums
    if (is.factor(ret) && all(levels(ret) == c("FALSE", "TRUE"))) {
        ret <- as.logical(ret)
    }

    if (!is.null(encoding) && (endsWith(encoding, "-scalar") || encoding == "string")) {
        attr(ret, "encoding-scalar") <- TRUE
    }
    if (length(dim(ret)) > 1)
        ret <- t(ret)
    ret
}

#' @importFrom rhdf5 H5Aopen H5Aread H5Aclose H5Dread
#' @importMethodsFrom rhdf5 &
read_categorical <- function(group, encoding) {
    check_encodingversion(group, encoding, "0.2.0")

    orderedattr <- H5Aopen(group, "ordered")
    ordered <- as.logical(H5Aread(orderedattr))
    H5Aclose(orderedattr)

    values <- as.integer(H5Dread(h5autoclose(group & "codes")))
    labels <- H5Dread(h5autoclose(group & "categories"))
    convert_categoricals(values, labels)
}

#' @importFrom rhdf5 H5Dread
#' @importMethodsFrom rhdf5 &
read_nullable <- function(group, encoding) {
    check_encodingversion(group, encoding, "0.1.0")
    values <- H5Dread(h5autoclose(group & "values"))
    mask <- as.logical(H5Dread(h5autoclose(group & "mask")))
    values[mask] <- NA
    if (endsWith(encoding, "-boolean")) {
        values <- as.logical(values)
    }
    values
}

.read_funcs <- list("array"=read_array,
                    "csr_matrix"=read_sparse_matrix,
                    "csc_matrix"=read_sparse_matrix,
                    "dataframe"=read_dataframe,
                    "dict"=read_dict,
                    "numeric-scalar"=read_array,
                    "string"=read_array,
                    "categorical"=read_categorical,
                    "string-array"=read_array,
                    "nullable-integer"=read_nullable,
                    "nullable-boolean"=read_nullable)

#' @importFrom rhdf5 H5Aexists H5Aopen H5Aread H5Aclose H5Iget_type H5Iget_name
read_attribute <- function(attr) {
    ret <- NULL
    if (H5Aexists(attr, "encoding-type")) {
        encattr <- H5Aopen(attr, "encoding-type")
        encoding <- H5Aread(encattr)
        H5Aclose(encattr)

        func <- switch(encoding,
                      "array"=read_array,
                      "csr_matrix"=read_sparse_matrix,
                      "csc_matrix"=read_sparse_matrix,
                      "dataframe"=read_dataframe,
                      "dict"=read_dict,
                      "numeric-scalar"=read_array,
                      "string"=read_array,
                      "categorical"=read_categorical,
                      "string-array"=read_array,
                      "nullable-integer"=read_nullable,
                      "nullable-boolean"=read_nullable)
        if (is.null(func)) {
            warning("Unknown encoding ", encoding, " for element ", H5Iget_name(attr))
        } else {
            ret <- func(attr, encoding)
        }
    }

    if (is.null(ret)) {
        if (H5Iget_type(attr) == "H5I_GROUP")
            ret <- read_dict(attr, NULL, strict=FALSE)
        else
            ret <- read_array(attr, NULL, strict=FALSE)
    }
    ret
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
        metadata(se) <- read_dict(h5autoclose(view & "uns"), NULL, strict=FALSE)
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
            cmapdset <- h5autoclose(h5 & paste("obsmap", mod, sep="/"))
            cmap <- H5Dread(cmapdset)

            idx <- which(cmap > 0)
            data.frame(assay=mod, primary=rownames(metadata)[idx], colname=colnames(modalities[[mod]])[cmap[idx]])
        })
        sampleMap <- do.call(rbind, samplemaps)
        rownames(sampleMap) <- NULL
        args$sampleMap <- sampleMap
    }

    if (H5Lexists(h5, "uns"))
        args$metadata <- read_dict(h5autoclose(h5 & "uns"), NULL, strict=FALSE)

    # Close the connection
    H5Fclose(h5)

    # Create a MAE object
    do.call(MultiAssayExperiment, args)
}


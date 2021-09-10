#' @importFrom rhdf5 H5Iget_type h5readAttributes h5ls H5Aexists H5Aopen H5Aread H5Aclose
read_with_index <- function(dataset) {
    cls <- H5Iget_type(dataset)
    if (cls == "H5I_GROUP" && H5Aexists(dataset, "encoding-type")) {
        encattr <- H5Aopen(dataset, "encoding-type")
        encoding <- H5Aread(encattr, "encoding-type")
        H5Aclose(encattr)
        if (encoding != "dataframe") {
            warning(paste0("Unknown encoding ", encoding, " when attempting to read data frame"))
            return(data.frame())
        }
        # Table is saved as a group rather than a dataset
        indexcol <- "_index"
        if (H5Aexists(dataset, "_index")) {
            indexattr <- H5Aopen(dataset, "_index")
            indexcol <- H5Aread(indexattr)
            H5Aclose(indexattr)
        }

        columns <- h5ls(dataset, recursive=FALSE, datasetinfo=FALSE)$name
        columns <- columns[columns != "__categories"]

        col_list <- lapply(columns, function(name) {
            values <- H5Dread(dataset & name)
            col <- dataset & name
            if (H5Aexists(col, "categories")) {
                attr <- H5Aopen(col, "categories")
                if (H5is_attr_reference(attr)) {
                    labels <- H5deref_attr_reference(attr)
                    values <- factor(as.integer(values), labels=H5Dread(labels))
                } else {
                    warning(paste0("found categories attribute for column ", name, ", but it is not a reference"))
                }
                H5Aclose(attr)
            }
            values
        })
        table <- data.frame(Reduce(cbind, col_list))
        colnames(table) <- columns

        if (indexcol %in% colnames(table)) {
            rownames(table) <- table[,indexcol,drop=TRUE]
            table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
        }

        # Fix column order
        if (H5Aexists(dataset, "column-order")) {
            orderedattr <- H5Aopen(dataset, "column-order")
            ordered_columns <- H5Aread(orderedattr)
            H5Aclose(orderedattr)

            ordered_columns <- ordered_columns[ordered_columns != indexcol]
            table <- table[,ordered_columns[ordered_columns %in% columns],drop=FALSE]
        }
    } else {
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
    }
        table
}

#' @importFrom rhdf5 H5Iget_type H5Aexists H5Aopen H5Aread H5Aclose H5Dread
read_matrix <- function(dataset) {
    if (H5Iget_type(dataset) == "H5I_GROUP" && H5Aexists(dataset, "encoding-type")) {
        encattr <- H5Aopen(dataset, "encoding-type")
        encoding <- H5Aread(dataset, "encoding")
        H5Aclose(encattr)
        if (encoding %in% c("csr_matrix", "csc_matrix")) {
            i <- H5Dread(dataset & "indices")
            p <- H5Dread(dataset & "indptr")
            x <- H5Dread(dataset & "data")
            shapeattr <- H5Aopen(dataset, "shape")
            shape <- H5Aread(shapeattr)
            H5Aclose(shapeattr)
            if (encoding == "csr_matrix") {
                sparseMatrix(j=i, p=p, x=x, dims=shape, repr="R", index1=FALSE)
            } else {
                sparseMatrix(i=i, p=p, x=x, dims=shape, repr="C", index1=FALSE)
            }
        } else {
            warning(paste0("Unknown encoding ", encoding, "when attempting to read matrix"))
            matrix()
        }
    } else {
        H5Dread(dataset)
    }
}

#' @details MultiAssayExperiment-helpers
#'
#' @description Create a MultiAssayExperiment or a Seurat object from the .h5mu file
#'
#' @importFrom rhdf5 h5ls H5Dread
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
ReadH5MU <- function(file) {
    # Connect to the the file
    h5 <- open_and_check_mudata(file)

    # Check all the assays are written
    assays <- h5ls(h5 & "mod", recursive=FALSE)$name

    # Create global colData
    metadata <- read_with_index(h5 & "obs")

    # Create an experiments list
    modalities <- lapply(assays, function(mod) {
        view <- h5 & paste("mod", mod, sep="/")
        X <- read_matrix(view & "X")

        var <- read_with_index(view & "var")

        obs <- read_with_index(view  & "obs")
        if (is("obs", "data.frame"))
            rownames(obs) <- paste(mod, rownames(obs), sep="-")

        viewnames <- h5ls(view, recursive=FALSE)$name
        if ("obsm" %in% viewnames) {
            obsmnames <- h5ls(view & "obsm", recursive=FALSE)$name
            obsm <- lapply(obsnames, function(space) {
                H5Dread(view & paste("obsm", space, sep="/"))
            })
            names(obsm) <- obsmnames
            se <- SingleCellExperiment(assays=SimpleList(counts=X), rowData=var, colData=obs, reducedDims=obsm)
        } else {
            se <- SummarizedExperiment(assays=SimpleList(counts=X), rowData=var, colData=obs)
        }

        se
    })
    names(modalities) <- assays

    # Create sampleMap
    mapping <- lapply(assays, function(mod) {
        primary <- colnames(modalities[[mod]])

        view <- h5 & paste("mod", mod, sep="/")
        view_attr <- h5attributes(view[["obs"]])
        indexcol <- "_index"
        if ("_index" %in% names(view_attr)) {
          indexcol <- view_attr$`_index`
        }
        obs_names <- view[['obs']]$read()[,indexcol,drop=TRUE]
        sm <- data.frame(primary = obs_names,
                         colname = rownames(colData(modalities[[mod]])),
                         stringsAsFactors = FALSE)
        sm
    })

    obsmap <- do.call(rbind, mapping)
    obsmap["assay"] <- rep(assays, times=vapply(mapping, nrow, 1))
    obsmap <- obsmap[,c("assay", "primary", "colname")]

    # Close the connection
    h5$close_all()

    # Create a MAE object
    MultiAssayExperiment(modalities, metadata, obsmap)
}


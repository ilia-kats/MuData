setGeneric("WriteH5MU", function(object, file, overwrite = TRUE) standardGeneric("WriteH5MU"))
setGeneric("WriteH5AD", function(object, assay, file, overwrite = TRUE) standardGeneric("WriteH5AD"))

#' @details MultiAssayExperiment-helpers
#'
#' @description Save MultiAssayExperiment object to .h5mu file
#' Please note than the primary key is used as obs_names 
#' so the behaviour of WriteH5MU when there are multiple samples
#' for one primary key is not guaranteed.
#'
#' @import hdf5r
#'
#' @exportMethod WriteH5MU
setMethod("WriteH5MU", "MultiAssayExperiment", function(object, file, overwrite) {
  h5 <- open_h5(file)

  obs <- as.data.frame(colData(object), stringsAsFactors = FALSE)
  obs[["_index"]] <- rownames(obs)
  obs_columns <- colnames(obs)
  obs <- obs[,c("_index", obs_columns)]

  h5[["obs"]] <- obs
  h5attr(h5[["obs"]], "_index") <- "_index"
  h5attr(h5[["obs"]], "column-order") <- colnames(obs)

  modalities <- names(experiments(object))

  h5$create_group("mod")
  vars <- lapply(modalities, function(mod) {
    mod_group <- h5$create_group(paste0("mod/", mod))

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
    obs <- data.frame(obs, stringsAsFactors = FALSE)
    obs_columns <- colnames(obs)
    obs[["_index"]] <- rownames(obs)
    obs <- obs[,c("_index", obs_columns),drop=FALSE]

    obs_dataset <- mod_group$create_dataset("obs", obs)
    h5attr(obs_dataset, "_index") <- "_index"
    if (length(obs_columns) > 0)
      h5attr(obs_dataset, "column-order") <- obs_columns

    # .obsm
    if (is(object[[mod]], "SingleCellExperiment")) {
      obsm <- reducedDims(object[[mod]])
      if (length(obsm) > 0) {
        mod_obsm <- mod_group$create_group("obsm")
        lapply(names(obsm), function(space) {
          mod_obsm$create_dataset(space, t(obsm[[space]]))
        })
      }
    }

    # X
    x <- object[[mod]]
    x <- x[,meta$colname]
    h5[[paste0("mod/", mod, "/X")]] <- assay(x)

    # .var
    var <- data.frame("mod" = rep(mod, nrow(x)), row.names = rownames(x), stringsAsFactors = FALSE)
    var_columns <- colnames(var)
    var[["_index"]] <- rownames(var)
    var <- var[,c("_index", var_columns)]
    h5[[paste0("mod/", mod, "/var")]] <- var
    h5attr(h5[[paste0("mod/", mod, "/var")]], "_index") <- "_index"

    var
  })

  var <- do.call(rbind, vars)
  h5[["var"]] <- var
  h5attr(h5[["var"]], "_index") <- "_index"
  h5attr(h5[["var"]], "column-order") <- colnames(var)


  finalize_mudata(h5)

  TRUE
})

# setMethod("WriteH5AD", "MultiAssayExperiment", function(object, assay, file, overwrite) {
#   h5 <- open_h5(file)



#' @description Save an assay to .h5ad / AnnData object
#'
#' @import hdf5r
#' @importFrom Matrix t
#'
WriteH5ADHelper <- function(object, assay, root) {

  mod_object <- Seurat::GetAssay(object, assay)

  # .obs
  obs_group <- root$create_group("obs")
  # There is no local metadata in Seurat objects
  obs_names <- colnames(object)
  obs <- data.frame(row.names = obs_names)
  write_data_frame(obs_group, obs)

  # .var
  var <- mod_object@meta.features

  # Define highly variable features, if any
  if ('var.features' %in% slotNames(mod_object)) {
    message("Defining highly variable features...")
    var$highly_variable <- rownames(var) %in% mod_object@var.features
  }

  var_group <- root$create_group("var")
  write_data_frame(var_group, var)

  # .X, .layers['counts']. .raw.X
  if ('counts' %in% slotNames(mod_object)) {
    x_counts <- Seurat::GetAssayData(mod_object, 'counts')
    sparse_type <- ifelse(class(x_counts) == "dgCMatrix", "csc_matrix", "csr_matrix")
    # case 1: only counts available
    if (!(('data' %in% slotNames(mod_object)) || ('scale.data' %in% slotNames(mod_object)))) {
      if ("i" %in% slotNames(x_counts)) {
        # sparse matrix
        x_counts <- Matrix::t(x_counts)
        counts_group <- root$create_group("X")
        counts_group$create_dataset("indices", x_counts@i)
        counts_group$create_dataset("indptr", x_counts@p)
        counts_group$create_dataset("data", x_counts@x)
        h5attr(counts_group, "shape") <- dim(x_counts)
        counts_group$create_attr("encoding-type", sparse_type, space=H5S$new("scalar"))
        counts_group$create_attr("encoding-version", "0.1.0", space=H5S$new("scalar"))
      } else {
        # dense matrix
        root$create_dataset("X", x_counts)
      }
    } else {
      layers_group <- root$create_group("layers")
      if ("i" %in% slotNames(x_counts)) {
        # sparse matrix
        print("Writing sparse counts...")
        x_counts <- Matrix::t(x_counts)
        counts_group <- layers_group$create_group("counts")
        # counts_group$create_dataset("indices", x_counts@i)
        # counts_group$create_dataset("indptr", x_counts@p)
        # counts_group$create_dataset("data", x_counts@x)
        counts_group[["indices"]] <- x_counts@i
        counts_group[["indptr"]] <- x_counts@p
        counts_group[["data"]] <- x_counts@x
        h5attr(counts_group, "shape") <- dim(x_counts)
        counts_group$create_attr("encoding-type", sparse_type, space=H5S$new("scalar"))
        counts_group$create_attr("encoding-version", "0.1.0", space=H5S$new("scalar"))
      } else {
        # dense matrix
        layers_group$create_dataset("counts", x_counts)
      }
      if ('data' %in% slotNames(mod_object)) {
        x_data <- Seurat::GetAssayData(mod_object, 'data')
        sparse_type <- ifelse(class(x_data) == "dgCMatrix", "csc_matrix", "csr_matrix")
        if ('scale.data' %in% slotNames(mod_object) && length(mod_object@scale.data) > 0) {
          # case 2: counts, data, and scale.data are available
          # .X
          x_scaled <- t(Seurat::GetAssayData(mod_object, 'scale.data'))
          root$create_dataset("X", x_scaled)
          # .raw
          raw_group <- root$create_group("raw")
          if ("i" %in% slotNames(x_data)) {
            # sparse matrix
            x_data <- Matrix::t(x_data)
            data_group <- raw_group$create_group("X")
            data_group$create_dataset("indices", x_data@i)
            data_group$create_dataset("indptr", x_data@p)
            data_group$create_dataset("data", x_data@x)
            h5attr(data_group, "shape") <- dim(x_data)
            data_group$create_attr("encoding-type", sparse_type, space=H5S$new("scalar"))
            data_group$create_attr("encoding-version", "0.1.0", space=H5S$new("scalar"))
          } else {
            # dense matrix
            raw_group$create_dataset("X", t(x_data))
          }
        } else {
          # case 3: counts and data are available but not scale.data
          if ("i" %in% slotNames(x_data)) {
            # sparse matrix
            x_data <- Matrix::t(x_data)
            print("Writing sparse X...")
            data_group <- root$create_group("X")
            data_group$create_dataset("indices", x_data@i)
            data_group$create_dataset("indptr", x_data@p)
            data_group$create_dataset("data", x_data@x)
            h5attr(data_group, "shape") <- dim(x_data)
            data_group$create_attr("encoding-type", sparse_type, space=H5S$new("scalar"))
            data_group$create_attr("encoding-version", "0.1.0", space=H5S$new("scalar"))
          } else {
            # dense matrix
            root$create_dataset("X", x_data)
          }
        }
      }
      # 'data' should to be available when 'scale.data' is available
    }
  }

  finalize_anndata_internal(root)


  TRUE
}


#' @details Seurat-helpers
#'
#' @description Save Seurat object to .h5mu file
#'
#' @import hdf5r
#'
#' @exportMethod WriteH5MU
setMethod("WriteH5MU", "Seurat", function(object, file, overwrite) {
  h5 <- open_h5(file)

  # .obs
  obs_group <- h5$create_group("obs")
  obs <- object@meta.data

  write_data_frame(obs_group, obs)

  modalities <- Seurat::Assays(object)

  h5$create_group("mod")
  var_names <- lapply(modalities, function(mod) {
    mod_group <- h5$create_group(paste0("mod/", mod))

    WriteH5ADHelper(object, mod, mod_group)

    mod_object <- object[[mod]]
    rownames(mod_object)
  })

  # global .var will only contain rownames
  var <- data.frame(row.names = do.call(c, var_names))
  var_group <- h5$create_group("var")
  write_data_frame(var_group, var)

  # .obsm
  if ('reductions' %in% slotNames(object)) {
    for (red_name in names(object@reductions)) {
      red <- object@reductions[[red_name]]
      emb <- t(red@cell.embeddings)
      assay <- red@assay.used
      if (!is.null(assay) && assay != "") {
        obsm <- h5$create_group(paste0("mod/", assay, "/obsm"))
      } else {
        obsm <- h5$create_group("obsm")
      }
      obsm$create_dataset(paste0("X_", red_name), emb)
    }
  }

  # TODO: .varm
  # object@reductions$...@feature.loadings
  # also: # object@reductions$...@stdev

  # TODO: .obsp
  # object@graphs (sparse)
  # object@neighbors (k nearest neighbours)

  finalize_mudata(h5)

  TRUE
})


write_data_frame <- function(attr_group, attr_df) {
  attr_columns <- colnames(attr_df)

  attr_df["_index"] <- rownames(attr_df)

  attr_df <- attr_df[,c("_index", attr_columns),drop=FALSE]

  categories <- list()
  for (col in colnames(attr_df)) {
    v <- attr_df[[col]]
    if ("factor" %in% class(v)) {
      # Write a factor
      categories[[col]] <- levels(v)
      attr_group$create_dataset(col, as.integer(v) - 1, dtype = h5types$H5T_NATIVE_INT)
    } else {
      attr_group$create_dataset(col, v)
    }
  }
  if (length(categories) > 0) {
    cats <- attr_group$create_group("__categories")
    for (cat in names(categories)) {
      cat_dataset <- cats$create_dataset(cat, categories[[cat]])
      cat_dataset$create_attr("ordered", FALSE, space = H5S$new("scalar"))
      attr_group[[cat]]$create_attr("categories", 
                                    cats$create_reference(cat), 
                                    space = H5S$new("scalar"))
    }
  }

  # Write attributes
  attr_group$create_attr("_index", "_index", space = H5S$new("scalar"))
  attr_group$create_attr("encoding-type", "dataframe", space = H5S$new("scalar"))
  attr_group$create_attr("encoding-version", "0.1.0", space = H5S$new("scalar"))
  if (length(attr_columns) > 0) {
    attr_group$create_attr("column-order", attr_columns)
  } else {
    # When there are no columns, null buffer can't be written to a file.
    attr_group$create_attr("column-order", dtype=h5types$H5T_NATIVE_DOUBLE, space=H5S$new("simple", 0, 0))
  }
  
}
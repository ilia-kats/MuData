.mudataversion <- "0.1.0"
.anndataversion <- "0.1.0"
.name <- paste0(getPackageName(), ".r")
.version <- as.character(packageVersion(getPackageName()))

#' @importFrom methods setGeneric setMethod
#' @import Matrix
NULL

#' @importFrom rhdf5 H5Pcreate H5Pset_userblock H5Fcreate
open_h5 <- function(filename) {
    h5p_create <- H5Pcreate("H5P_FILE_CREATE")
    res <- H5Pset_userblock(h5p_create, 512)
    if (res < 0) {
        stop("could not set HDF5 user block")
    }
    file <- H5Fcreate(filename, fcpl=h5p_create, native=FALSE)
    H5Pclose(h5p_create)
    file
}

#' @importFrom rhdf5 h5writeAttribute H5Fget_name H5Fclose
finalize_mudata <- function(h5) {
    h5writeAttribute("MuData", h5, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute(.mudataversion, h5, "encoding-version", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute(.name, h5, "encoder", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute(.version, h5, "encoder-version", variableLengthString=TRUE, asScalar=TRUE)

    filename <- H5Fget_name(h5)
    H5Fclose(h5)
    h5 <- file(filename, "r+b")
    writeChar(paste0("MuData (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
    close(h5)
}

#' @importFrom rhdf5 h5writeAttribute
finalize_anndata_internal <- function(h5) {
    h5writeAttribute("AnnData", h5, "encoding-type", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute(.anndataversion, h5, "encoding-version", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute(.name, h5, "encoder", variableLengthString=TRUE, asScalar=TRUE)
    h5writeAttribute(.version, h5, "encoder-version", variableLengthString=TRUE, asScalar=TRUE)
}

#' @importFrom rhdf5 H5Fis_hdf5 H5Fopen
open_and_check_mudata <- function(filename) {
    if (readChar(filename, 6) != "MuData") {
        if (H5Fis_hdf5(filename)) {
            warning("The HDF5 file was not created by muon, we can't guarantee that everything will work correctly", call.=FALSE)
        } else (
            stop("The file is not an HDF5 file", call.=FALSE)
        )
    }

    H5Fopen(filename, flags="H5F_ACC_RDONLY", native=FALSE)
}

#' @importFrom rhdf5 H5iget_type H5Iis_valid H5Dclose H5Gclose H5Aclose H5Fclose
h5autoclose <- function(obj) {
    obj <- force(obj)
    do.call(on.exit, list(substitute({
        if (H5Iis_valid(obj)) {
            switch(H5Iget_type(obj),
                H5I_FILE=H5Fclose(obj),
                H5I_GROUP=H5Gclose(obj),
                H5I_DATASET=H5Dclose(obj),
                H5I_ATTR=H5Aclose(obj)
            )
        }
    }), add=TRUE),
    envir=parent.frame())
    obj
}

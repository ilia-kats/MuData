# Read an .h5ad file and create a `SingleCellExperiment`.

In file-backed mode, the main `X` matrix is not read into memory, but
references the HDF5 file and its required parts are read on demand. This
requires the HDF5Array package to be installed.

## Usage

``` r
readH5AD(file, backed = FALSE)
```

## Arguments

- file:

  Path to the .h5ad file.

- backed:

  Whether to use file-backed mode.

## Value

A `SingleCellExperiment`.

## Examples

``` r
data(miniACC, package="MultiAssayExperiment")
writeH5AD(miniACC[[1]], "miniacc.h5ad")
sce <- readH5AD("miniacc.h5ad")
```

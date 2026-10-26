# Read an .h5mu file and create a [`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).

In file-backed mode, the main `X` matrices are not read into memory, but
reference the HDF5 file and their required parts are read on demand.
This requires the HDF5Array package to be installed.

## Usage

``` r
readH5MU(file, backed = FALSE)
```

## Arguments

- file:

  Path to the .h5mu file.

- backed:

  Whether to use file-backed mode.

## Value

A
[`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-class.html)

## Examples

``` r
data(miniACC, package="MultiAssayExperiment")
writeH5MU(miniACC, "miniacc.h5mu")
mae <- readH5MU("miniacc.h5mu")
#> Reading as SingleCellExperiment where the original object class is matrix
#> Warning: sampleMap[['assay']] coerced with as.factor()
```

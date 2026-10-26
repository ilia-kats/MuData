# Save an experiment to an .h5ad file.

Note that NA values are not supported by HDF5, and therefore by h5ad.
The behavior of this function if NAs are present is undefined.

## Usage

``` r
writeH5AD(object, file, overwrite)
```

## Arguments

- object:

  The object to save.

- file:

  Name of the file to save to.

- overwrite:

  Currently unused.

## Value

NULL, invisibly

## Examples

``` r
data(miniACC, package="MultiAssayExperiment")
writeH5AD(miniACC[[1]], "miniacc.h5ad")
```

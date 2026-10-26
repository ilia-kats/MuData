# Save a [`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-class.html) to an .h5mu file.

Note that NA values are not supported by HDF5, and therefore by h5mu.
The behavior of this function if NAs are present is undefined.

## Usage

``` r
writeH5MU(object, file, overwrite)
```

## Arguments

- object:

  A
  [`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-class.html).

- file:

  Name of the file to save to.

- overwrite:

  Currently unused.

## Value

NULL, invisibly

## Examples

``` r
data(miniACC, package="MultiAssayExperiment")
writeH5MU(miniACC, "miniacc.h5mu")
```

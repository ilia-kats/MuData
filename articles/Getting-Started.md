# Getting started with MuData for MultiAssayExperiment

## Introduction

Multimodal data format — [MuData](https://mudata.readthedocs.io/) — [has
been
introduced](https://www.biorxiv.org/content/10.1101/2021.06.01.445670v1)
to address the need for cross-platform standard for sharing large-scale
multimodal omics data. Importantly, it develops ideas of and is
compatible with [AnnData](https://anndata.readthedocs.io/) standard for
storing raw and derived data for unimodal datasets.

In Bioconductor, multimodal datasets have been stored in
[MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment)
(MAE) objects. This `MuData` package provides functionality to read data
from MuData files into MAE objects as well as to save MAE objects into
H5MU files.

## Installation

The most recent dev build can be installed from GitHub:

``` r

library(remotes)
remotes::install_github("ilia-kats/MuData")
```

Stable version of `MuData` will be available in future bioconductor
versions.

## Loading libraries

``` r

library(MuData)
library(MultiAssayExperiment)

library(rhdf5)
```

## Writing H5MU files

We’ll use a simple MAE object from the `MultiAssayExperiment` package
that we’ll then save in a H5MU file.

``` r

data(miniACC)
miniACC
#> A MultiAssayExperiment object of 5 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 5:
#>  [1] RNASeq2GeneNorm: SummarizedExperiment with 198 rows and 79 columns
#>  [2] gistict: SummarizedExperiment with 198 rows and 90 columns
#>  [3] RPPAArray: SummarizedExperiment with 33 rows and 46 columns
#>  [4] Mutations: matrix with 97 rows and 90 columns
#>  [5] miRNASeqGene: SummarizedExperiment with 471 rows and 80 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

We will now write its contents into an H5MU file with `WriteH5MU`:

``` r

writeH5MU(miniACC, "miniacc.h5mu")
```

## Reading H5MU files

We can manually check the top level of the structure of the file:

``` r

rhdf5::h5ls("miniacc.h5mu", recursive = FALSE)
#>   group   name     otype dclass dim
#> 0     /    mod H5I_GROUP           
#> 1     /    obs H5I_GROUP           
#> 2     /   obsm H5I_GROUP           
#> 3     / obsmap H5I_GROUP           
#> 4     /    uns H5I_GROUP           
#> 5     /    var H5I_GROUP
```

Or dig deeper into the file:

``` r

h5 <- rhdf5::H5Fopen("miniacc.h5mu")
h5&'mod'
#> HDF5 GROUP 
#>         name /mod
#>     filename 
#> 
#>              name     otype dclass dim
#> 0 Mutations       H5I_GROUP           
#> 1 RNASeq2GeneNorm H5I_GROUP           
#> 2 RPPAArray       H5I_GROUP           
#> 3 gistict         H5I_GROUP           
#> 4 miRNASeqGene    H5I_GROUP
rhdf5::H5close()
```

### Creating MAE objects from H5MU files

This package provides `ReadH5MU` to create an object with data from an
H5MU file. Since H5MU structure has been designed to accommodate more
structured information than MAE, only some data will be read. For
instance, MAE has no support for loading multimodal embeddings or
pairwise graphs.

``` r

acc <- readH5MU("miniacc.h5mu")
#> Reading as SingleCellExperiment where the original object class is matrix
#> Warning: sampleMap[['assay']] coerced with as.factor()
acc
#> A MultiAssayExperiment object of 5 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 5:
#>  [1] RNASeq2GeneNorm: SummarizedExperiment with 198 rows and 79 columns
#>  [2] gistict: SummarizedExperiment with 198 rows and 90 columns
#>  [3] RPPAArray: SummarizedExperiment with 33 rows and 46 columns
#>  [4] Mutations: SingleCellExperiment with 97 rows and 90 columns
#>  [5] miRNASeqGene: SummarizedExperiment with 471 rows and 80 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

Importantly, we recover the information from the original MAE object:

``` r

head(colData(miniACC)[,1:4])
#> DataFrame with 6 rows and 4 columns
#>                 patientID years_to_birth vital_status days_to_death
#>               <character>      <integer>    <integer>     <integer>
#> TCGA-OR-A5J1 TCGA-OR-A5J1             58            1          1355
#> TCGA-OR-A5J2 TCGA-OR-A5J2             44            1          1677
#> TCGA-OR-A5J3 TCGA-OR-A5J3             23            0            NA
#> TCGA-OR-A5J4 TCGA-OR-A5J4             23            1           423
#> TCGA-OR-A5J5 TCGA-OR-A5J5             30            1           365
#> TCGA-OR-A5J6 TCGA-OR-A5J6             29            0            NA
head(colData(acc)[,1:4])
#> DataFrame with 6 rows and 4 columns
#>                 patientID years_to_birth vital_status days_to_death
#>               <character>      <integer>    <integer>     <integer>
#> TCGA-OR-A5J1 TCGA-OR-A5J1             58            1          1355
#> TCGA-OR-A5J2 TCGA-OR-A5J2             44            1          1677
#> TCGA-OR-A5J3 TCGA-OR-A5J3             23            0            NA
#> TCGA-OR-A5J4 TCGA-OR-A5J4             23            1           423
#> TCGA-OR-A5J5 TCGA-OR-A5J5             30            1           365
#> TCGA-OR-A5J6 TCGA-OR-A5J6             29            0            NA
```

Features metadata is also recovered:

``` r

head(rowData(miniACC[["gistict"]]))
#> DataFrame with 6 rows and 3 columns
#>        Gene.Symbol    Locus.ID    Cytoband
#>        <character> <character> <character>
#> DIRAS3      DIRAS3        9077      1p31.3
#> MAPK14      MAPK14        1432     6p21.31
#> YAP1          YAP1       10413     11q22.1
#> CDKN1B      CDKN1B        1027     12p13.1
#> ERBB2        ERBB2        2064       17q12
#> G6PD          G6PD        2539        Xq28
head(rowData(acc[["gistict"]]))
#> DataFrame with 6 rows and 3 columns
#>        Gene.Symbol    Locus.ID    Cytoband
#>        <character> <character> <character>
#> DIRAS3      DIRAS3        9077      1p31.3
#> MAPK14      MAPK14        1432     6p21.31
#> YAP1          YAP1       10413     11q22.1
#> CDKN1B      CDKN1B        1027     12p13.1
#> ERBB2        ERBB2        2064       17q12
#> G6PD          G6PD        2539        Xq28
```

#### Backed objects

It is possible to read H5MU files while keeping matrices (both `.X` and
`.layers`) on disk.

``` r

acc_b <- readH5MU("miniacc.h5mu", backed = TRUE)
#> Reading as SingleCellExperiment where the original object class is matrix
#> Warning: sampleMap[['assay']] coerced with as.factor()
assay(acc_b, "RNASeq2GeneNorm")[1:5,1:3]
#> <5 x 3> DelayedMatrix object of type "double":
#>        TCGA-OR-A5J1-01A-11R-A29S-07 TCGA-OR-A5J2-01A-11R-A29S-07
#> DIRAS3                    1487.0317                       9.6631
#> MAPK14                     778.5783                    2823.6469
#> YAP1                      1009.6061                    2305.0590
#> CDKN1B                    2101.3449                    4248.9584
#> ERBB2                      651.2968                     246.4098
#>        TCGA-OR-A5J3-01A-11R-A29S-07
#> DIRAS3                      18.9602
#> MAPK14                    1061.7686
#> YAP1                      1561.2502
#> CDKN1B                    1348.5410
#> ERBB2                       90.0607
```

The data in the assay is a `DelayedMatrix` object:

``` r

class(assay(acc_b, "RNASeq2GeneNorm"))
#> [1] "DelayedMatrix"
#> attr(,"package")
#> [1] "DelayedArray"
```

This is in contrast to the `acc` object that has matrices in memory:

``` r

assay(acc, "RNASeq2GeneNorm")[1:5,1:3]
#>        TCGA-OR-A5J1-01A-11R-A29S-07 TCGA-OR-A5J2-01A-11R-A29S-07
#> DIRAS3                    1487.0317                       9.6631
#> MAPK14                     778.5783                    2823.6469
#> YAP1                      1009.6061                    2305.0590
#> CDKN1B                    2101.3449                    4248.9584
#> ERBB2                      651.2968                     246.4098
#>        TCGA-OR-A5J3-01A-11R-A29S-07
#> DIRAS3                      18.9602
#> MAPK14                    1061.7686
#> YAP1                      1561.2502
#> CDKN1B                    1348.5410
#> ERBB2                       90.0607
class(assay(acc, "RNASeq2GeneNorm"))
#> [1] "matrix" "array"
```

## References

- [Muon: multimodal omics analysis
  framework](https://www.biorxiv.org/content/10.1101/2021.06.01.445670)
  preprint

- [mudata](https://mudata.readthedocs.io/) (Python) documentation

- muon [documentation](https://muon.readthedocs.io/) and [web
  page](https://gtca.github.io/muon/)

## Session Info

``` r

sessionInfo()
#> R Under development (unstable) (2026-10-26 r90766)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] MultiAssayExperiment_1.38.0 SummarizedExperiment_1.42.0
#>  [3] Biobase_2.72.0              GenomicRanges_1.64.0       
#>  [5] Seqinfo_1.2.0               IRanges_2.46.0             
#>  [7] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#>  [9] MuData_1.17.1               rhdf5_2.56.1               
#> [11] S4Vectors_0.50.3            BiocGenerics_0.58.1        
#> [13] generics_0.1.4              Matrix_1.7-6               
#> [15] BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10                 SparseArray_1.12.3         
#>  [3] lattice_0.23-1              h5mread_1.4.1              
#>  [5] digest_0.6.39               evaluate_1.0.5             
#>  [7] grid_4.7.0                  bookdown_0.48              
#>  [9] fastmap_1.2.0               jsonlite_2.0.0             
#> [11] BiocManager_1.30.27         SingleCellExperiment_1.34.0
#> [13] HDF5Array_1.40.0            textshaping_1.0.5          
#> [15] jquerylib_0.1.4             abind_1.4-8                
#> [17] cli_3.6.6                   rlang_1.3.0                
#> [19] XVector_0.52.0              cachem_1.1.0               
#> [21] DelayedArray_0.38.2         yaml_2.3.12                
#> [23] otel_0.2.0                  S4Arrays_1.12.1            
#> [25] tools_4.7.0                 Rhdf5lib_2.0.0             
#> [27] R6_2.6.1                    lifecycle_1.0.5            
#> [29] fs_2.1.0                    htmlwidgets_1.6.4          
#> [31] ragg_1.5.2                  desc_1.4.3                 
#> [33] pkgdown_2.2.1               bslib_0.12.0               
#> [35] systemfonts_1.3.2           xfun_0.62                  
#> [37] knitr_1.52                  rhdf5filters_1.24.1        
#> [39] htmltools_0.5.9             rmarkdown_2.32             
#> [41] compiler_4.7.0
```

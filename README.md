# MuDataSeurat

`MuDataSeurat` is a package that provides I/O funcitonality for `.h5mu` files and [MultiAssayExperiment](http://waldronlab.io/MultiAssayExperiment/) objects.

You can learn more about multimodal data containers in the [`mudata` documentation](https://mudata.readthedocs.io/en/latest/io/mudata.html).

## Installation

```R
remotes::install_github("pmbio/MuDataMAE")
```

## Quick start

`MuDataMAE` provides a set of I/O operations for multimodal data.

`MuDataMAE` implements `WriteH5MU()` that saves MultiAssayExperiment objects to `.h5mu` files that can be further integrated into workflows in multiple programming languages, including the [`muon` Python library](https://github.com/pmbio/muon) and the [`Muon.jl` Julia library](https://github.com/pmbio/Muon.jl). `ReadH5MU()` reads `.h5mu` files into MultiAssayExperiment objects.


### Writing files

Start with an existing dataset, e.g. a [MultiAssayExperiment](http://waldronlab.io/MultiAssayExperiment/) object with five distinct modalities:

```R
library(MultiAssayExperiment)
data(miniACC)
```

`WriteH5MU()` allows to save the object into a `.h5mu` file:

```R
library(MuDataMAE)
WriteH5MU(miniACC, "miniACC.h5mu")
```

### Reading files

```R
miniACC <- ReadH5MU("miniACC.h5mu")
```

## Relevant projects

Other R packages for multimodal I/O include:

- [MuData](https://github.com/PMBio/MuDataSeurat) for [Seurat](https://github.com/satijalab/seurat/) objects
- [SeuratDisk](https://github.com/mojaveazure/seurat-disk)


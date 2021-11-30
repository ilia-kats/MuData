# MuData

[Documentation](https://pmbio.github.io/MuDataMAE/) | [Preprint](https://www.biorxiv.org/content/10.1101/2021.06.01.445670v1) | [Discord](https://discord.com/invite/MMsgDhnSwQ)

[![R-CMD-check](https://github.com/PMBio/MuDataMAE/workflows/R-CMD-check/badge.svg)](https://github.com/PMBio/MuDataMAE/actions)

`MuData` is a package that provides I/O funcitonality for `.h5mu` files and [MultiAssayExperiment](http://waldronlab.io/MultiAssayExperiment/) objects.

You can learn more about multimodal data containers in the reference [`mudata` documentation](https://mudata.readthedocs.io/en/latest/io/mudata.html).

## Installation
`MuData` uses [`rhdf5`](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) to access `.h5mu` and `.h5ad` files.
We use `rhdf5` over `hdf5r` to stay compatible with the rest of the Bioconductor ecosystem.
In particular, using `hdf5r` would make integrating with other packages building on `rhdf5`, such as `HDF5Array`, much more difficult, if not impossible. We have implemented necessary HDF5 features that the `.h5ad` and consequently `.h5mu` formats make use of upstream, including [file creation properties](https://github.com/grimbough/rhdf5/pull/95) and [object references](https://github.com/grimbough/rhdf5/pull/96).

In the meantime, the most recent dev `rhdf5` version from GitHub must be used. `rhdf5` and `MuData` can be installed by running

```R
remotes::install_github("grimbough/rhdf5")
remotes::install_github("pmbio/MuDataMAE")
```

## Quick start

`MuData` provides a set of I/O operations for multimodal data.

`MuData` implements `WriteH5MU()` that saves MultiAssayExperiment objects to `.h5mu` files that can be further integrated into workflows in multiple programming languages, including the [`muon` Python library](https://github.com/pmbio/muon) and the [`Muon.jl` Julia library](https://github.com/pmbio/Muon.jl). `ReadH5MU()` reads `.h5mu` files into MultiAssayExperiment objects.


### Writing files

Start with an existing dataset, e.g. a [MultiAssayExperiment](http://waldronlab.io/MultiAssayExperiment/) object with five distinct modalities:

```R
library(MultiAssayExperiment)
data(miniACC)
```

`WriteH5MU()` allows to save the object into a `.h5mu` file:

```R
library(MuData)
WriteH5MU(miniACC, "miniACC.h5mu")
```

### Reading files

```R
miniACC <- ReadH5MU("miniACC.h5mu")
```

## Relevant projects

Other R packages for multimodal I/O include:

- [MuDataSeurat](https://github.com/PMBio/MuDataSeurat) for [Seurat](https://github.com/satijalab/seurat/) objects
- [SeuratDisk](https://github.com/mojaveazure/seurat-disk)


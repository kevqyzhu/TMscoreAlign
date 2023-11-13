
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TMscoreAlign

<!-- badges: start -->
<!-- badges: end -->

`TMscoreAlign` is an R package that implements the TM-align program in R
([Zhang & Skolnick, 2005](https://doi.org/10.1002/prot.20264)). TM-align
is a method for aligning proteins based on TM-score (Template Modelling
score), which calculates topological similarity between two protein
structures ([Zhang & Skolnick,
2004](https://doi.org/10.1093/nar/gki524)). TM-score improves upon
existing means of structural comparison such as RMSD (root-mean-square
deviation) as it is length-independent and more sensitive to global
similarities between structures rather than local deviations. This
package is targeted for structural biologists who may use this package
to investigate different conformations of the same protein, informing a
structural basis of protein functions. The tool includes functions that
performs structure alignment, calculates TM-scores and RMSD, and
visualizes aligned structures. `TMscoreAlign` was developed using
`R version 4.3.1 (2023-06-16)`,
`Platform: aarch64-apple-darwin20 (64-bit)` and Running under:
`macOS Ventura 13.2.1`.

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("kevqyzhu/TMscoreAlign", build_vignettes = TRUE)
library("TMscoreAlign")
```

## Overview

``` r
ls("package:TMscoreAlign")
data(package = "TMscoreAlign") 
browseVignettes("TMscoreAlign")
```

There are 5 functions in `TMscoreAlign`

`get_alignment`: This function performs structure alignment between two
protein structures specified by PDB files. First, it selects common
residues using a user-specified method (either based on ‘index’ or
‘alignment’). Then, the alignment parameters, including translation and
rotation values, are optimized for TM-score (Template Modeling Score)
and returned to the user as a list. Other values such as TM-score and
RMSD are also returned in the list.

`optimize`: This function performs optimization to improve the alignment
of two protein structures based on their atomic coordinates. The
optimization aims to find the best parameters that minimize the
objective function (TM-score). The optimization can be restarted with
default values or continue from a given set of values.

`get_tmscore`: This function calculates the TM-Score (Template Modeling
Score) between two protein structures based on their alignment. The
alignment is performed using the specified method (‘alignment’ or
‘index’). This function is a wrapper for `get_alignment` so that the
user can directly calculate the TM-score of two protein structures given
their PDB files.

`get_rmsd`: This function calculates the RMSD (root mean square
deviation) between two protein structures based on their alignment. The
alignment is performed using the specified method (‘alignment’ or
‘index’). This function is a wrapper for `get_alignment` so that the
user can directly calculate the RMSD of two protein structures given
their PDB files.

`visualize_alignment_pdb`: This function visualizes the structural
alignment of two protein structures by displaying the aligned
coordinates in a 3D viewer. The viewer is set up using `3Dmol`, and the
alignment is represented by color-coded cartoon-style structures for
each protein chain. The user can specify the alignment parameters and
chain coloration.

Below is a flowchart demonstart the work flow of this package.

<img src="./inst/extdata/flowchart.jpeg" style="width:75.0%" /> See
`help(package = "DeregGenes")` for further details and references
provided by `citation("DeregGenes")`.

## Contributions

The package is created by Kevin Zhu.

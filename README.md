
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

There are 6 functions in `TMscoreAlign`

`get_alignment`: This function performs structure alignment between two
protein structures specified by PDB files. First, it selects common
residues using a user-specified method (either based on ‘index’ or
‘alignment’). Then, the alignment parameters, including translation and
rotation values, are initialized. These values can be further optimized
for TM-score (Template Modeling score) before being returned to the
user.

`optimize_alignment`: This function performs optimization to improve the
alignment of two protein structures based on their atomic coordinates.
It takes as input the translation and rotation values and aims to find
the best parameters that minimize the objective function (TM-score).
Note that this process can be initialized with default values or
continue from a given set of values. Optimization is performed via the
limited-memory Broyden–Fletcher–Goldfarb–Shanno algorithm.

`get_tmscore`: This function calculates the TM-score from the alignment
parameters and coordinates obtained in a structural alignment between
two protein structures.

`get_tm_samples`:This function calculates the TM local scores from the
alignment parameters and coordinates obtained in a structural alignment
between two protein structures.

`get_rmsd`: This function calculates the RMSD (root mean square
deviation) between two protein structures based on their alignment
parameters.

`visualize_alignment_pdb`: This function visualizes the structural
alignment of two protein structures by displaying the aligned
coordinates in a 3D viewer. The viewer is set up using `3Dmol`, and the
alignment is represented by color-coded cartoon-style structures for
each protein chain. The user can specify the alignment parameters and
chain coloration.

## Contributions

The package is created by Kevin Zhu. `get_alignment` calls a helper
function that uses `read.pdb`, `clean.pdb`, `pdbseq`, and `atom.select`
from the `bio3d` package to read the input PDB files. The
`pairwiseAlignment` and `subject` functions from `Biostrings` were used
to generate the sequence alignment. The `optimize_alignment` function
uses the `optim` function from the `stats` package in order to optimize
the alignment parameters for maximum TM-score. Furthermore, when
generating the transformation matrix to represent the alignment
parameters, the `cross` function from the `pracma` package was used to
generate the vector cross product. The visualization from this package
is accomplished using `r3dmol`. All functions from this package were
written by the author. No generative AI tools were used in the
development of this package.

## References

- Berman, H., Henrick, K. and Nakamura, H. Announcing the worldwide
  Protein Data Bank. *Nature Structural & Molecular Biology*, 10, 980
  (2003). <https://doi.org/10.1038/nsb1203-980>
- Borchers H (2022). *pracma: Practical Numerical Math Functions*. R
  package version 2.4.2 <https://CRAN.R-project.org/package=pracma>
- Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory
  algorithm for bound constrained optimization. *SIAM Journal on
  Scientific Computing*, 16(5), 1190–1208.
  <https://epubs.siam.org/doi/10.1137/0916069>
- Grant, B.J. et al. (2006). *Bio3D: An R package for the comparative
  analysis of protein structures*. *Bioinformatics*, 22, 2695–2696.
  <https://doi.org/10.1093/bioinformatics/btl461> R package version
  2.4.4, <https://cran.r-project.org/web/packages/bio3d/index.html>
- Jordi Cruzado
  (<https://math.stackexchange.com/users/96872/jordi-cruzado>), Explain
  3d transformation matrix…, URL (version: 2023-05-30):
  <https://math.stackexchange.com/q/532974>
- Jur van den Berg
  (<https://math.stackexchange.com/users/91768/jur-van-den-berg>),
  Calculate Rotation Matrix to align Vector $A$ to Vector $B$ in $3D$?,
  URL (version: 2016-09-01): <https://math.stackexchange.com/q/476311>
- Pagès H, Aboyoun P, Gentleman R, DebRoy S (2023). *Biostrings:
  Efficient manipulation of biological strings*. R package version
  2.68.1, <https://doi.org/10.18129/B9.bioc.Biostrings>
- R Core Team (2023). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>
- Su W, Johnston B (2021). *r3dmol: Create Interactive 3D Visualizations
  of Molecular Data*. R package version 0.1.2,
  <https://CRAN.R-project.org/package=r3dmol>
- Zhang, Y., and Skolnick, J. (2004). Scoring function for automated
  assessment of protein structure template quality. *Proteins,
  Structure, Function, and Bioinformatics*, 57(4), 702–710.
  <https://doi.org/10.1002/prot.20264>

## Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `TMscoreAlign` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/kevqyzhu/TMscoreAlign/issues). Many thanks to
those who provided feedback to improve this package.

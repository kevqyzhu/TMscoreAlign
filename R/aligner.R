# Purpose: Instantiate and optimize alignment
# Author: Kevin Zhu
# Date: 11.14.2023
# Version: 1.0.0
# Bugs and Issues: N/A

#' Load Data and Perform Sequence Alignment or Index-based Analysis
#'
#' This function reads two Protein Data Bank (PDB) files, extracts the CA atom
#' coordinates for common residues between the specified chains, and returns the
#' coordinates along with the number of common residues.
#'
#' @param pdb_file1 Character string. Path to the first PDB file.
#' @param pdb_file2 Character string. Path to the second PDB file.
#' @param chain1 Character string. Chain identifier for the first PDB file
#'  (default is 'A').
#' @param chain2 Character string. Chain identifier for the second PDB file
#'  (default is 'A').
#' @param method Character string. Method for residue selection.
#'  Options: "alignment" or "index"\cr
#'   - "alignment": Perform sequence alignment using pairwiseAlignment from
#'      Bioconductor.\cr
#'   - "index": Use residue indices based on the specified chains to identify
#'      common residues.
#' @return A list containing the following components:\cr
#'   - coord1: 3D coordinates (matrix) of CA atoms for common residues in the
#'      first structure.\cr
#'   - coord2: 3D coordinates (matrix) of CA atoms for common residues in the
#'      second structure.\cr
#'   - N: Number of common residues.
#'
#' @examples
#' \dontrun{
#' # Example 1: Perform sequence alignment
#' result_alignment <- load_data_alignment("path/to/file1.pdb",
#'                                         "path/to/file2.pdb",
#'                                         method = "alignment")
#'
#' # Example 2: Use residue indices for common residue selection
#' result_index <- load_data_alignment("path/to/file1.pdb",
#'                                     "path/to/file2.pdb",
#'                                     method = "index")
#' }
#'
#' @references
#' Grant, B.J. et al. (2006). Bio3D: An R package for the comparative analysis
#' of protein structures. \emph{Bioinformatics}, 22, 2695--2696.
#' \href{https://doi.org/10.1093/bioinformatics/btl461}{Link}, R package version
#' 2.4.4, \href{https://cran.r-project.org/web/packages/bio3d/index.html}{Link}
#'
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2023). Biostrings: Efficient
#' manipulation of biological strings.
#' \href{https://doi.org/10.18129/B9.bioc.Biostrings}{Link}, R package version
#' 2.68.1, \href{https://bioconductor.org/packages/Biostrings}{Link}.
#'
#' @seealso
#' \code{\link{pairwiseAlignment}} for sequence alignment.
#' \code{\link{atom.select}} for atom selection in the \code{\link{bio3d}}
#'   package.
#'
#' @export
#' @importFrom bio3d read.pdb clean.pdb pdbseq atom.select
#' @importFrom Biostrings pairwiseAlignment subject
load_data_alignment <- function(pdb_file1, pdb_file2, chain1 = 'A',
                                chain2 = 'A', method = "alignment") {
  # Read PDB structures
  pdb_data1 <- clean.pdb(read.pdb(pdb_file1), fix.chain = TRUE)
  pdb_data2 <- clean.pdb(read.pdb(pdb_file2), fix.chain = TRUE)

  seq1 <- paste(pdbseq(pdb_data1, inds = NULL, aa1 = TRUE), collapse = " ")
  seq2 <-paste(pdbseq(pdb_data2, inds = NULL, aa1 = TRUE), collapse = " ")

  if (method == "alignment") {
    # Perform sequence alignment
    alignment <- pairwiseAlignment(seq1, seq2)
    aligned_seq <- subject(alignment)

    # Convert the aligned sequence to a character string
    aligned_str <- as.character(aligned_seq)
    common_residues <- which(strsplit(aligned_str, " ")[[1]]!="-")
  }  else if (method == "index") {
    common_residues <- intersect(unique(pdb_data2$atom[pdb_data2$atom$chain ==
                                                         chain1,]$resno),
                                 unique(pdb_data2$atom[pdb_data2$atom$chain ==
                                                         chain2,]$resno))
  }

  # Extract coordinates of CA atoms for common residues
  sele_1 <- atom.select(pdb_data1, 'calpha', resno=common_residues,
                        chain=chain1)
  coord1 <- matrix(pdb_data1$xyz[sele_1$xyz], nrow=3, byrow=FALSE)
  coord1 <- rbind(coord1, rep(1, ncol(coord1)))


  sele_2 <- atom.select(pdb_data2, 'calpha', resno=common_residues,
                        chain=chain2)
  coord2 <- matrix(pdb_data2$xyz[sele_2$xyz], nrow=3, byrow=FALSE)
  coord2 <- rbind(coord2, rep(1, ncol(coord2)))

  return(list(coord1 = coord1, coord2 = coord2, N = length(common_residues)))
}

#' Optimize Protein Structure Alignment Parameters
#'
#' This function optimizes the alignment parameters for better accuracy of the
#' structural alignment between two protein structures.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements:\cr
#'   - N: Numeric. Number of common residues in the alignment.\cr
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.\cr
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.\cr
#'   - values: Numeric vector. Initial alignment parameters.
#' @param restart Logical. If TRUE, the optimization starts from the default
#'   values. If FALSE, it continues from the current alignment parameters.
#'
#' @return The input alignment list updated with optimized parameters
#'
#' @examples
#' \dontrun{
#' # Example: Optimize alignment parameters
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment", optimize = FALSE)
#' optimized_results <- optimize(alignment_results, restart = TRUE)
#' print(optimized_results)
#' }
#'
#' @references
#' Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995).
#' A limited memory algorithm for bound constrained optimization.
#' \emph{SIAM Journal on Scientific Computing}, 16(5), 1190–1208.
#' \href{https://epubs.siam.org/doi/10.1137/0916069}{Link}
#'
#' @seealso
#' \code{\link{get_default_values}} for parameter initialization.
#' \code{\link{tm}} and \code{\link{rmsd}} for TM-score and RMSD calculation.
#' \code{\link{optim}} for optimization functions.
#'
#' @export
#' @importFrom stats optim
optimize <- function(alignment, restart = TRUE) {
  coord1 <- alignment$coord1
  coord2 <- alignment$coord2
  d0_values <- estimate_d0(alignment$N)
  d02 <- d0_values$d02

  if (restart) {
    default_values <- get_default_values(coord1, coord2)
  } else {
    default_values <- alignment$values
  }

  method <- "L-BFGS-B"
  result <- stats::optim(par = default_values, fn = tm, coord1 = coord1,
                         coord2 = coord2, d02 = d02, method = method,
                         control = list(fnscale = -1))
  alignment$values <- result$par
  return(alignment)
}

#' Get Protein Structure Alignment
#'
#' This function performs structural alignment between two protein structures
#' specified by PDB files. It can use either a specified alignment method or
#' index-based alignment. The alignment parameters, including translation and
#' rotation values, are initialized and can be further optimized.
#'
#' @param pdb1 Character. Path to the PDB file of the first protein structure.
#' @param pdb2 Character. Path to the PDB file of the second protein structure.
#' @param chain1 Character. Chain identifier for the first protein structure.
#'   Defaults to 'A'.
#' @param chain2 Character. Chain identifier for the second protein structure.
#'   Defaults to 'A'.
#' @param method Character. Alignment method to be used. Options are "alignment"
#'   for sequence-based alignment or "index" for index-based alignment.
#' @param optimize Logical. If TRUE, the alignment parameters are optimized for
#'   better accuracy. Defaults to TRUE.
#'
#' @return A list containing the following elements: \cr
#'   - N: Numeric. Number of common residues in the alignment. \cr
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure. \cr
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure. \cr
#'   - values: Numeric vector. Alignment parameters obtained from structure
#'     alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Get structural alignment
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment", optimize = TRUE)
#' print(alignment_results)
#' }
#'
#' @seealso
#' \code{\link{load_data_alignment}} for loading data for structure alignment.
#' \code{\link{optimize}} for optimizing alignment parameters.
#'
#' @export
get_alignment <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', method,
                          optimize = TRUE) {
  # Load data alignment
  data <- load_data_alignment(pdb1, pdb2, chain1, chain2, method)

  default_values <- get_default_values(data$coord1, data$coord2)

  alignment <- list(N = data$N,
                    coord1 = data$coord1,
                    coord2 = data$coord2,
                    values = default_values)
  if (optimize) {
    # Optimize the alignment
    alignment <- optimize(alignment)
  }

  return(alignment)
}


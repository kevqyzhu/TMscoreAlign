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
#'  Options: "alignment" or "index"
#'   - "alignment": Perform sequence alignment using pairwiseAlignment from
#'      Bioconductor.
#'   - "index": Use residue indices based on the specified chains to identify
#'      common residues.
#' @return A list containing the following components:
#'   - coord1: 3D coordinates (matrix) of CA atoms for common residues in the
#'      first structure.
#'   - coord2: 3D coordinates (matrix) of CA atoms for common residues in the
#'      second structure.
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
#'   The list should contain the following elements:
#'   - N: Numeric. Number of common residues in the alignment.
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.
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

#' Get Default Values for Structure Alignment Parameters
#'
#' This function calculates default values for the parameters used in structure
#' alignment. The parameters include translation along x, y, z axes
#' (dx, dy, dz), and rotation angles (theta, phi, psi) based on the given 3D
#' coordinates of two structures.
#'
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @return A numeric vector containing default values for the structure
#'   alignment parameters:
#'   - dx: Translation along the x-axis.
#'   - dy: Translation along the y-axis.
#'   - dz: Translation along the z-axis.
#'   - theta: Rotation angle around the x-axis.
#'   - phi: Rotation angle around the y-axis.
#'   - psi: Rotation angle around the z-axis.
#'
#' @examples
#' \dontrun{
#' # Example: Get default values for structure alignment parameters
#' default_values <- get_default_values(coord1, coord2)
#' }
#'
#' @references
#' Borchers H (2022). pracma: Practical Numerical Math Functions. R package
#' version 2.4.2, \href{https://CRAN.R-project.org/package=pracma}{Link}.
#' Mathematics Stack Exchange:
#' \href{https://math.stackexchange.com/questions/180418}{Link}
#'
#' @export
#' @importFrom pracma cross
get_default_values <- function(coord1, coord2) {
  # Initialize a list to store the default values
  values <- list()

  # Set initial translation values to zero
  values$dx <- 0
  values$dy <- 0
  values$dz <- 0

  # Calculate the mean displacement along each axis
  dist <- rowMeans(coord1 - coord2)
  values$dx <- dist[1]
  values$dy <- dist[2]
  values$dz <- dist[3]

  # Calculate normalized vectors and cross product
  vec1 <- coord1[-nrow(coord1), 2] - coord1[-nrow(coord1), ncol(coord1)]
  vec2 <- coord2[-nrow(coord2), 2] - coord2[-nrow(coord2), ncol(coord1)]
  vec1 <- vec1 / sqrt(sum(vec1^2))
  vec2 <- vec2 / sqrt(sum(vec2^2))
  v <- cross(vec1, vec2)

  # Calculate rotation parameters
  s <- sqrt(sum(v^2)) + .Machine$double.eps
  c <- sum(vec1 * vec2)
  vx <- matrix(c(0, -v[3], v[2], v[3], 0, -v[1], -v[2], v[1], 0), nrow = 3,
               byrow = TRUE)
  rotation_matrix <- diag(3) + vx + vx%*%vx * (1 - c) / (s * s)
  values$theta <- atan2(rotation_matrix[3, 2], rotation_matrix[3, 3])
  values$phi <- atan2(-rotation_matrix[3, 1],
                      sqrt(rotation_matrix[3, 2]^2 + rotation_matrix[3, 3]^2))
  values$psi <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1])

  # Convert the list to a numeric vector
  return(unlist(values))
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
#' @return A list containing the following elements:
#'   - N: Numeric. Number of common residues in the alignment.
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.
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

#' Get TM Score from Protein Structure Alignment
#'
#' This function calculates the TM score from the alignment parameters and
#' coordinates obtained in a structural alignment between two protein
#' structures.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements:
#'   - N: Numeric. Number of common residues in the alignment.
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.
#'   - values: Numeric vector. Alignment parameters.
#'
#' @return Numeric. TM score of the protein structure alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate TM score from alignment
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment")
#' optimized_alignment <- optimize(alignment_results)
#' tmscore <- get_tmscore(optimized_alignment)
#' print(tmscore)
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining alignment details.
#' \code{\link{optimize}} for optimizing alignment parameters.
#' \code{\link{estimate_d0}} for estimating initial distance parameters.
#' \code{\link{tm}} for calculating TM-score.
#'
#' @export
get_tmscore <- function(alignment) {
  # Return the TM-Score from the alignment
  tmscore <- tm(alignment$values,
                alignment$coord1,
                alignment$coord2,
                estimate_d0(alignment$N)$d02)
  return(tmscore)
}

#' Get TM local scores from Protein Structure Alignment
#'
#' This function calculates the TM local scores from the alignment parameters
#' and coordinates obtained in a structural alignment between two protein
#' structures.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements:
#'   - N: Numeric. Number of common residues in the alignment.
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.
#'   - values: Numeric vector. Alignment parameters.
#'
#' @return Numeric vector. TM samples of the protein structure alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate TM samples from alignment
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment")
#' tm_samples <- get_tm_samples(alignment_results)
#' print(tm_samples)
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining alignment details.
#' \code{\link{optimize}} for optimizing alignment parameters.
#' \code{\link{estimate_d0}} for estimating initial distance parameters.
#' \code{\link{tm_samples}} for calculating TM local scores
#'
#' @export
get_tm_samples <- function(alignment) {
  # Return the TM-samoles from the alignment
  tm_samples <- tm_samples(alignment$values,
                alignment$coord1,
                alignment$coord2,
                estimate_d0(alignment$N)$d02)
  return(tm_samples)
}

#' Get Root Mean Square Deviation (RMSD) from Protein Structure Alignment
#'
#' This function calculates the Root Mean Square Deviation (RMSD) from the
#' alignment parameters and coordinates obtained in a structural alignment
#' between two protein structures.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements:
#'   - N: Numeric. Number of common residues in the alignment.
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.
#'   - values: Numeric vector. Alignment parameters.
#'
#' @return Numeric. Root Mean Square Deviation (RMSD) of the protein structure
#'   alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate RMSD from alignment
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment")
#' rmsd_value <- get_rmsd(alignment_results)
#' print(rmsd_value)
#' }
#'
#' @references
#' RMSD calculation:
#' \href{https://zhanglab.ccmb.med.umich.edu/}{Zhang Lab}.
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining structural alignment between two
#'   protein structures.
#' \code{\link{optimize}} for optimizing alignment parameters.
#' \code{\link{estimate_d0}} for estimating initial distance parameters.
#' \code{\link{rmsd}} for calculating RMSD.
#'
#' @export
get_rmsd <- function(alignment) {
  rmsd <- rmsd(alignment$values,
                alignment$coord1,
                alignment$coord2)

  # Return the Root Mean Square Deviation (RMSD)
  return(rmsd)
}

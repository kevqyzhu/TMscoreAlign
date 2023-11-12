library(pracma)

#' Load Data and Perform Sequence Alignment or Index-based Analysis
#'
#' This function reads two Protein Data Bank (PDB) files, extracts the CA atom coordinates
#' for common residues between the specified chains, and returns the coordinates along with
#' the number of common residues.
#'
#' @param pdb_file1 Character string. Path to the first PDB file.
#' @param pdb_file2 Character string. Path to the second PDB file.
#' @param chain1 Character string. Chain identifier for the first PDB file (default is 'A').
#' @param chain2 Character string. Chain identifier for the second PDB file (default is 'A').
#' @param method Character string. Method for residue selection. Options: "alignment" or "index"
#'   - "alignment": Perform sequence alignment using pairwiseAlignment from Bioconductor.
#'   - "index": Use residue indices based on the specified chains to identify common residues.
#' @return A list containing the following components:
#'   - coord1: 3D coordinates (matrix) of CA atoms for common residues in the first structure.
#'   - coord2: 3D coordinates (matrix) of CA atoms for common residues in the second structure.
#'   - N: Number of common residues.
#'
#' @examples
#' \dontrun{
#' # Example 1: Perform sequence alignment
#' result_alignment <- load_data_alignment("path/to/file1.pdb", "path/to/file2.pdb", method = "alignment")
#'
#' # Example 2: Use residue indices for common residue selection
#' result_index <- load_data_alignment("path/to/file1.pdb", "path/to/file2.pdb", method = "index")
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
#' \code{\link{atom.select}} for atom selection in the \code{\link{bio3d}} package.
#'
#' @export
#' @importFrom bio3d read.pdb clean.pdb pdbseq atom.select
#' @importFrom Biostrings pairwiseAlignment
load_data_alignment <- function(pdb_file1, pdb_file2, chain1 = 'A', chain2 = 'A', method = "alignment") {
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
    common_residues <- intersect(unique(pdb_data2$atom[pdb_data2$atom$chain == chain1,]$resno),
                                 unique(pdb_data2$atom[pdb_data2$atom$chain == chain2,]$resno))
  }

  # Extract coordinates of CA atoms for common residues
  sele_1 <- atom.select(pdb_data1, 'calpha', resno=common_residues, chain=chain1)
  coord1 <- matrix(pdb_data1$xyz[sele_1$xyz], nrow=3, byrow=FALSE)
  coord1 <- rbind(coord1, rep(1, ncol(coord1)))


  sele_2 <- atom.select(pdb_data2, 'calpha', resno=common_residues, chain=chain2)
  coord2 <- matrix(pdb_data2$xyz[sele_2$xyz], nrow=3, byrow=FALSE)
  coord2 <- rbind(coord2, rep(1, ncol(coord2)))

  return(list(coord1 = coord1, coord2 = coord2, N = length(common_residues)))
}

#' Optimize Protein Structure Alignment
#'
#' This function performs optimization to improve the alignment of two protein structures based
#' on their atomic coordinates. The optimization aims to find the best parameters that minimize
#' the objective function, which includes the TM-score and RMSD calculations. The optimization
#' can be restarted with default values or continue from a given set of values.
#'
#' @param coord1 Numeric matrix. Coordinates of the reference protein structure.
#' @param coord2 Numeric matrix. Coordinates of the target protein structure.
#' @param d02 Numeric. Reference distance for the TM-score calculation.
#' @param values Numeric vector. Initial values for the optimization parameters.
#'   If provided, the optimization continues from these values; otherwise, default values are used.
#' @param restart Logical. If TRUE, the optimization is restarted with default values; otherwise, it continues.
#' @return A list containing the optimized values, TM-score, and RMSD.
#'   - values: Numeric vector. Optimized parameters for the alignment.
#'   - tmscore: Numeric. TM-score of the aligned structures.
#'   - rmsd: Numeric. Root Mean Square Deviation (RMSD) between the aligned structures.
#'
#' @examples
#' \dontrun{
#' # Example: Optimize protein structure alignment
#' reference_structure <- matrix(rnorm(300), ncol = 3)
#' target_structure <- matrix(rnorm(300), ncol = 3)
#' reference_d0 <- 10
#' optimization_results <- optimize(coord1 = reference_structure, coord2 = target_structure, d02 = reference_d0)
#' }
#'
#' @references
#' Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995).
#' A limited memory algorithm for bound constrained optimization.
#' \emph{SIAM Journal on Scientific Computing}, 16(5), 1190–1208.
#' \href{https://doi.org/10.1137/091606}{Link}
#'
#' @seealso
#' \code{\link{get_default_values}} and \code{\link{get_current_values}} for parameter initialization.
#' \code{\link{tm}} and \code{\link{rmsd}} for TM-score and RMSD calculation.
#' \code{\link{stats::optim}} for optimization functions.
#'
#' @export
#' @importFrom stats optim
optimize <- function(coord1, coord2, d02, values = NULL, restart = TRUE) {
  if (restart) {
    default_values <- get_default_values(coord1, coord2)
  } else {
    default_values <- values
  }
  method <- "L-BFGS-B"
  result <- stats::optim(par = default_values, fn = tm, coord1 = coord1, coord2 = coord2, d02 = d02, method = method, control = list(fnscale = -1))
  values <- result$par
  return(list(values = values, tmscore = tm(values, coord1, coord2, d02), rmsd = rmsd(values, coord1, coord2)))
}

#' Get Default Values for Structure Alignment Parameters
#'
#' This function calculates default values for the parameters used in structure alignment.
#' The parameters include translation along x, y, z axes (dx, dy, dz), and rotation angles
#' (theta, phi, psi) based on the given 3D coordinates of two structures.
#'
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms (rows: dimensions, columns: atoms).
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms (rows: dimensions, columns: atoms).
#' @return A numeric vector containing default values for the structure alignment parameters:
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
#' Mathematics Stack Exchange: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/897677#897677
#'
#' @export
#' @importfrom pracma cross
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
  vx <- matrix(c(0, -v[3], v[2], v[3], 0, -v[1], -v[2], v[1], 0), nrow = 3, byrow = TRUE)
  rotation_matrix <- diag(3) + vx + vx%*%vx * (1 - c) / (s * s)
  values$theta <- atan2(rotation_matrix[3, 2], rotation_matrix[3, 3])
  values$phi <- atan2(-rotation_matrix[3, 1], sqrt(rotation_matrix[3, 2]^2 + rotation_matrix[3, 3]^2))
  values$psi <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1])

  # Convert the list to a numeric vector
  return(unlist(values))
}

#' Get Structure Alignment Parameters
#'
#' This function performs structure alignment between two protein structures specified by PDB files.
#' It uses the given method to select common residues and optimizes the alignment parameters based
#' on 3D coordinates. The resulting alignment parameters, including translation and rotation values,
#' are returned.
#'
#' @param pdb1 Character string. Path to the first PDB file.
#' @param pdb2 Character string. Path to the second PDB file.
#' @param chain1 Character string. Chain identifier for the first PDB file (default is 'A').
#' @param chain2 Character string. Chain identifier for the second PDB file (default is 'A').
#' @param method Character string. Method for residue selection. Options: "alignment" or "index"
#'   - "alignment": Perform sequence alignment using pairwiseAlignment from Bioconductor.
#'   - "index": Use residue indices based on the specified chains to identify common residues.
#' @return A list containing the following components:
#'   - values: Optimized parameters for structure alignment.
#'   - tmscore: TM-score value calculated with the optimized parameters.
#'   - rmsd: RMSD value calculated with the optimized parameters.
#'
#' @examples
#' \dontrun{
#' # Example: Get structure alignment parameters using sequence alignment method
#' alignment_seq <- get_alignment("path/to/file1.pdb", "path/to/file2.pdb", chain1 = 'A', chain2 = 'B', method = "alignment")
#'
#' # Example: Get structure alignment parameters using index-based method
#' alignment_index <- get_alignment("path/to/file1.pdb", "path/to/file2.pdb", chain1 = 'A', chain2 = 'B', method = "index")
#' }
#'
#' @references
#' - For residue selection: Please refer to the documentation of the load_data_alignment function.
#' - For parameter optimization: Please refer to the documentation of the optimize function.
#'
#' @seealso
#' \code{\link{load_data_alignment}} for loading and aligning PDB data.
#' \code{\link{optimize}} for optimizing structure alignment parameters.
#'
#' @export
get_alignment <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', method) {
  # Load data alignment
  data <- load_data_alignment(pdb1, pdb2, chain1, chain2, method)

  # Estimate d0
  d0_values <- estimate_d0(data$N)
  d02 <- d0_values$d02

  # Optimize the alignment
  alignment <- optimize(data$coord1, data$coord2, d02)
  return(alignment)
}

#' Get TM-Score Between Two Protein Structures
#'
#' This function calculates the TM-Score (Template Modeling Score) between two protein structures
#' based on their alignment. The alignment is performed using the specified method ('alignment' or 'index').
#' The TM-Score is a measure of structural similarity and is commonly used in structural bioinformatics.
#'
#' @param pdb1 Character. Path to the first PDB file containing the coordinates of the first protein structure.
#' @param pdb2 Character. Path to the second PDB file containing the coordinates of the second protein structure.
#' @param chain1 Character. Chain identifier for the first protein structure. Default is 'A'.
#' @param chain2 Character. Chain identifier for the second protein structure. Default is 'A'.
#' @param method Character. Alignment method to be used. Options are 'alignment' for sequence-based alignment
#'   or 'index' for index-based alignment.
#' @return The TM-Score between the two protein structures based on the specified alignment method.
#'
#' @examples
#' \dontrun{
#' # Example 1: Calculate TM-Score using sequence-based alignment
#' tm_score_seq <- get_tmscore("structure1.pdb", "structure2.pdb", chain1 = 'A', chain2 = 'A', method = 'alignment')
#'
#' # Example 2: Calculate TM-Score using index-based alignment
#' tm_score_index <- get_tmscore("structure1.pdb", "structure2.pdb", chain1 = 'A', chain2 = 'A', method = 'index')
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function, and
#' Bioinformatics}, 57(4), 702–710. \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining alignment details.
#'
#' @export
get_tmscore <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', method) {
  # Get alignment details
  alignment <- get_alignment(pdb1, pdb2, chain1, chain2, method)

  # Return the TM-Score from the alignment
  return(alignment$tmscore)
}

#' Calculate Root Mean Square Deviation (RMSD) Between Two Protein Structures
#'
#' This function calculates the Root Mean Square Deviation (RMSD) between two protein structures
#' based on their atomic coordinates. The RMSD is a measure of the average distance between
#' corresponding atoms in the two structures, providing an indication of structural similarity.
#'
#' @param pdb1 Character. Path to the PDB file of the first protein structure.
#' @param pdb2 Character. Path to the PDB file of the second protein structure.
#' @param chain1 Character. Chain identifier for the first protein structure (default: 'A').
#' @param chain2 Character. Chain identifier for the second protein structure (default: 'A').
#' @param method Character. The method used for structure alignment, either "alignment" or "index".
#' @return The Root Mean Square Deviation (RMSD) between the two protein structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate RMSD between two protein structures
#' pdb_file1 <- "path/to/structure1.pdb"
#' pdb_file2 <- "path/to/structure2.pdb"
#' chain_id1 <- 'A'
#' chain_id2 <- 'B'
#' alignment_method <- "alignment"
#' rmsd_value <- get_rmsd(pdb_file1, pdb_file2, chain_id1, chain_id2, alignment_method)
#' }
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining structural alignment between two protein structures.
#'
#' @export
get_rmsd <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', method) {
  # Obtain structural alignment between two protein structures
  alignment <- get_alignment(pdb1, pdb2, chain1, chain2, method)

  # Return the Root Mean Square Deviation (RMSD)
  return(alignment$rmsd)
}

#' Write Transformed Coordinates to PDB File
#'
#' This function writes the transformed coordinates obtained from a structure alignment to a new
#' PDB file. The transformed coordinates are generated by applying the alignment parameters to
#' the original coordinates of the second protein structure. The resulting PDB file includes the
#' transformed coordinates alongside the original coordinates of the first structure.
#'
#' @param alignment List. Structure alignment results, including alignment parameters and transformed coordinates.
#'   The list should contain the following elements:
#'   - values: Numeric vector. Alignment parameters obtained from structure alignment.
#'   - rmsd: Numeric. Root Mean Square Deviation (RMSD) between the two protein structures.
#' @param outputfile Character. The name of the output PDB file to be created (default: "out.pdb").
#' @param appended Logical. If TRUE, the transformed coordinates are appended to the original coordinates
#'   in the output PDB file. If FALSE, the output PDB file contains only the transformed coordinates.
#' @param pdb1 Character. Path to the PDB file of the first protein structure.
#' @param pdb2 Character. Path to the PDB file of the second protein structure.
#' @param chain_1 Character. Chain identifier for the first protein structure in output file.
#' @param chain_2 Character. Chain identifier for the second protein structure in output file.
#'
#' @examples
#' \dontrun{
#' # Example: Write transformed coordinates to a new PDB file
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb", chain1 = 'A', chain2 = 'A', method = "alignment")
#' write_pdb(alignment_results, outputfile = "aligned_structure.pdb", appended = TRUE, pdb1 = "structure1.pdb", pdb2 = "structure2.pdb", chain_1 = 'A', chain_2 = 'A')
#' }
#'
#' @references
#' PDB formatting: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining structural alignment between two protein structures.
#' \code{\link{get_matrix}} for obtaining the transformation matrix from alignment parameters.
#'
#' @export
write_pdb <- function(alignment, outputfile = "out.pdb", appended = TRUE, pdb1, pdb2, chain_1, chain_2) {

  values <- alignment$values
  matrix <- get_matrix(values)
  # browser()

  out <- file(outputfile, "w")
  atomid <- 1

  if (appended) {
    # Process the first PDB file
    lines <- readLines(pdb1)
    for (line in lines) {
      if (!grepl("^ATOM", line) || (substring(line, 22, 22) != " " && substring(line, 22, 22) != chain_1)) {
        next
      }
      cat(substr(line, 1, 6), sprintf("%4d", atomid), substr(line, 13, 20), "A", substr(line, 24, nchar(line)), "\n", file = out)
      atomid <- atomid + 1
    }
  }

  # Process the second PDB file
  lines <- readLines(pdb2)
  for (line in lines) {
    if (!grepl("^ATOM", line) || (substring(line, 22, 22) != " " && substring(line, 22, 22) != chain_2)) {
      next
    }

    x <- as.numeric(substr(line, 32, 38))
    y <- as.numeric(substr(line, 39, 46))
    z <- as.numeric(substr(line, 48, 54))

    vec <- c(x, y, z, 1)
    transformed_vec <- matrix %*% vec

    cat(substr(line, 1, 6), sprintf("%4d", atomid), substr(line, 13, 20), "B", substr(line, 24, 29),
        sprintf("%8.3f%8.3f%8.3f", transformed_vec[1], transformed_vec[2], transformed_vec[3]), substr(line, 56, nchar(line)), "\n", file = out)
    atomid <- atomid + 1
  }

  close(out)
}

#' Visualize Protein Structure Alignment Using 3Dmol
#'
#' This function visualizes the structural alignment of two protein structures by displaying
#' the aligned coordinates in a 3D viewer. The viewer is set up using 3Dmol, and the alignment
#' is represented by color-coded cartoon-style structures for each protein chain.
#'
#' @param alignment_pdb Character. Path to the PDB file containing the aligned coordinates.
#'   The PDB file should include both the original and transformed coordinates of the protein structures.
#' @param chain1 Character. Color code (hexadecimal) for the first protein structure (default: "#636efa").
#' @param chain2 Character. Color code (hexadecimal) for the second protein structure (default: "#ff7f0e").
#' @return A 3Dmol viewer displaying the aligned protein structures.
#'
#' @examples
#' \dontrun{
#' # Example: Visualize protein structure alignment
#' alignment_pdb_file <- "aligned_structure.pdb"
#' chain1_color <- "#636efa"  # Blue
#' chain2_color <- "#ff7f0e"  # Orange
#' visualize_alignment_pdb(alignment_pdb_file, chain1 = chain1_color, chain2 = chain2_color)
#' }
#'
#' @references
#' Su W, Johnston B (2021). r3dmol: Create Interactive 3D Visualizations of
#' Molecular Data. R package version 0.1.2.
#' \href{https://CRAN.R-project.org/package=r3dmol}{Link}.
#'
#' @seealso
#' \code{\link{write_pdb}} for generating a PDB file with transformed coordinates.
#' \code{\link{get_alignment}} for obtaining structural alignment between two protein structures.
#'
#' @export
#' @import r3dmol
visualize_alignment_pdb <- function(alignment_pdb = "out.pdb", chain1 = "#636efa", chain2 = "#ff7f0e") {
  r3dmol(                         # Set up the initial viewer
    viewer_spec = m_viewer_spec(
      cartoonQuality = 40,
      lowerZoomLimit = 50,
      upperZoomLimit = 350
    )
  ) %>%
    m_add_model(                  # Add model to scene
      data = alignment_pdb,
      format = "pdb"
    ) %>%
    m_zoom_to() %>%               # Zoom to encompass the whole scene
    m_set_style(                  # Set style of specific selection
      sel = m_sel(chain = "A"),      # (selecting by secondary)
      style = m_style_cartoon(
        color = chain1,
        arrows = TRUE
      )
    ) %>%
    m_set_style(                  # Style the alpha helix
      sel = m_sel(chain = "B"),      # (selecting by alpha helix)
      style = m_style_cartoon(
        color = chain2,
        arrows = TRUE
      )
    ) %>%
    m_rotate(                     # Rotate the scene by given angle on given axis
      angle = 90,
      axis = "y"
    ) %>%
    m_spin()                      # Animate the scene by spinning it
}


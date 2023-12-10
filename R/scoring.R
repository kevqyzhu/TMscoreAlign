# Purpose: Score alignment
# Author: Kevin Zhu
# Date: 12.10.2023
# Version: 1.0.0
# Bugs and Issues: N/A

#' Get TM Score from Protein Structure Alignment
#'
#' This function calculates the TM score from the alignment parameters and
#' coordinates obtained in a structural alignment between two protein
#' structures. Note that unlike \code{\link{calculate_tmscore}}, this function
#' is a wrapper that is directly accessible to the user. This allows the
#' TM-Score to be easily obtained from the alignment parameters.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements: \cr
#'   - N: Numeric. Number of common residues in the alignment. \cr
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.\cr
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.\cr
#'   - values: Numeric vector. Alignment parameters.
#'
#' @return Numeric. TM score of the protein structure alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate TM score from alignment
#' pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' alignment_results <- get_alignment(pdb_file1, pdb_file2,
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment"
#'                                   )
#' optimized_alignment <- optimize_alignment(alignment_results)
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
#' \code{\link{optimize_alignment}} for optimizing alignment parameters.
#' \code{\link{estimate_d0}} for estimating initial distance parameters.
#' \code{\link{calculate_tmscore}} for calculating TM-score.
#'
#' @export
get_tmscore <- function(alignment) {
  if (typeof(alignment) != "list") {
    stop("Alignment type must be List.")
    }

  if (!setequal(names(alignment), c("N", "coord1", "coord2", "values"))) {
    stop("Alignment does not have the correct elements.")
    }

  if (typeof(alignment$N) != "integer") {
    stop("The N in alignment must be an integer.")
    }

  if (length(dim(alignment$coord1)) != 2 |
      length(dim(alignment$coord2)) != 2) {
    stop("The coord1 and coord2 matrices in alignment must be 2D matrices.")
    }

  if (dim(alignment$coord1)[1] != 4 |
      dim(alignment$coord2)[1] != 4) {
    stop("The first dimension of coord1 and coord2 matrices must be 4.")
    }

  if (dim(alignment$coord1)[2] != alignment$N |
      dim(alignment$coord2)[2] != alignment$N) {
    stop("The second dimension of coord1 and coord2 matrices must be equal to
         N.")
    }

  if (!is.vector(alignment$values)) {
    stop("The values in alignment must be a vector.")
    }

  if (!setequal(names(alignment$values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("The values in alignment does not have the correct elements.")
    }

  # Return the TM-Score from the alignment
  tmscore <- calculate_tmscore(alignment$values,
                               alignment$coord1,
                               alignment$coord2,
                               estimate_d0(alignment$N)$d02
                               )
  return(tmscore)
}

#' Get TM local scores from Protein Structure Alignment
#'
#' This function calculates the TM local scores from the alignment parameters
#' and coordinates obtained in a structural alignment between two protein
#' structures. Note that unlike \code{\link{calculate_tm_samples}}, this
#' function is a wrapper that is directly accessible to the user. This allows
#' the TM local scores to be easily obtained from the alignment parameters.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements:\cr
#'   - N: Numeric. Number of common residues in the alignment. \cr
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.\cr
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.\cr
#'   - values: Numeric vector. Alignment parameters.
#'
#' @return Numeric vector. TM samples of the protein structure alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate TM samples from alignment
#' pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' alignment_results <- get_alignment(pdb_file1, pdb_file2,
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment"
#'                                   )
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
#' \code{\link{optimize_alignment}} for optimizing alignment parameters.
#' \code{\link{estimate_d0}} for estimating initial distance parameters.
#' \code{\link{calculate_tm_samples}} for calculating TM local scores
#'
#' @export
get_tm_samples <- function(alignment) {
  if (typeof(alignment) != "list") {
    stop("Alignment type must be List.")
    }

  if (!setequal(names(alignment), c("N", "coord1", "coord2", "values"))) {
    stop("Alignment does not have the correct elements.")
    }

  if (typeof(alignment$N) != "integer") {
    stop("The N in alignment must be an integer.")
    }

  if (length(dim(alignment$coord1)) != 2 |
      length(dim(alignment$coord2)) != 2) {
    stop("The coord1 and coord2 matrices in alignment must be 2D matrices.")
    }

  if (dim(alignment$coord1)[1] != 4 |
      dim(alignment$coord2)[1] != 4) {
    stop("The first dimension of coord1 and coord2 matrices must be 4.")
    }

  if (dim(alignment$coord1)[2] != alignment$N |
      dim(alignment$coord2)[2] != alignment$N) {
    stop("The second dimension of coord1 and coord2 matrices must be equal to
         N.")
    }

  if (!is.vector(alignment$values)) {
    stop("The values in alignment must be a vector.")
    }

  if (!setequal(names(alignment$values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("The values in alignment does not have the correct elements.")
    }

  # Return the TM-samoles from the alignment
  tm_samples <- calculate_tm_samples(alignment$values,
                                     alignment$coord1,
                                     alignment$coord2,
                                     estimate_d0(alignment$N)$d02
                                     )
  return(tm_samples)
}

#' Get Root Mean Square Deviation (RMSD) from Protein Structure Alignment
#'
#' This function calculates the Root Mean Square Deviation (RMSD) from the
#' alignment parameters and coordinates obtained in a structural alignment
#' between two protein structures. Note that unlike
#' \code{\link{calculate_rmsd}}, this function is a wrapper that is directly
#' accessible to the user. This allows the RMSD to be easily obtained from the
#' alignment parameters.
#'
#' @param alignment List. Structure alignment results, including alignment
#'   parameters, coordinates, and other information.
#'   The list should contain the following elements:\cr
#'   - N: Numeric. Number of common residues in the alignment.\cr
#'   - coord1: Matrix. 3D coordinates of CA atoms for common residues in the
#'     first structure.\cr
#'   - coord2: Matrix. 3D coordinates of CA atoms for common residues in the
#'     second structure.\cr
#'   - values: Numeric vector. Alignment parameters.
#'
#' @return Numeric. Root Mean Square Deviation (RMSD) of the protein structure
#'   alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate RMSD from alignment
#' pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' alignment_results <- get_alignment(pdb_file1, pdb_file2,
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment"
#'                                   )
#' rmsd_value <- get_rmsd(alignment_results)
#' print(rmsd_value)
#' }
#'
#' @seealso
#' \code{\link{get_alignment}} for obtaining structural alignment between two
#'   protein structures.
#' \code{\link{optimize_alignment}} for optimizing alignment parameters.
#' \code{\link{estimate_d0}} for estimating initial distance parameters.
#' \code{\link{calculate_rmsd}} for calculating RMSD.
#'
#' @export
get_rmsd <- function(alignment) {
  if (typeof(alignment) != "list") {
    stop("Alignment type must be List.")
    }

  if (!setequal(names(alignment), c("N", "coord1", "coord2", "values"))) {
    stop("Alignment does not have the correct elements.")
    }

  if (typeof(alignment$N) != "integer") {
    stop("The N in alignment must be an integer.")
    }

  if (length(dim(alignment$coord1)) != 2 |
      length(dim(alignment$coord2)) != 2) {
    stop("The coord1 and coord2 matrices in alignment must be 2D matrices.")
    }

  if (dim(alignment$coord1)[1] != 4 |
      dim(alignment$coord2)[1] != 4) {
    stop("The first dimension of coord1 and coord2 matrices must be 4.")
    }

  if (dim(alignment$coord1)[2] != alignment$N |
      dim(alignment$coord2)[2] != alignment$N) {
    stop("The second dimension of coord1 and coord2 matrices must be equal to
         N.")
    }

  if (!is.vector(alignment$values)) {
    stop("The values in alignment must be a vector.")
    }

  if (!setequal(names(alignment$values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("The values in alignment does not have the correct elements.")
    }

  rmsd <- calculate_rmsd(alignment$values,
                         alignment$coord1,
                         alignment$coord2
                         )

  # Return the Root Mean Square Deviation (RMSD)
  return(rmsd)
}

#' Estimate d0 for Structure Alignment
#'
#' This function estimates the d0 parameter for structure alignment based on the
#' number of common residues. The formula used for estimation is derived from
#' empirical observations and is commonly used in structure alignment methods.
#' The square of the estimated d0 value (d02) is returned.
#'
#' @param N Numeric. Number of common residues between two structures.
#'
#' @return A list containing the following component:\cr
#'   - d02: The square of the estimated d0 value for structure alignment.
#'
#' @examples
#' \dontrun{
#' # Example: Estimate d0 for structure alignment with 50 common residues
#' d0_values <- estimate_d0(N = 50)
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @export
estimate_d0 <- function(N) {
  if (typeof(N) != "integer") {
    stop("N must be an integer.")
    }

  d0 <- 1.24 * abs(N - 15)^(1/3) - 1.8
  d02 <- d0^2

  return(list(d0 = d0, d02 = d02))
}

#' Calculate Distances Between Transformed Coordinates
#'
#' This function calculates the Euclidean distances between transformed
#' coordinates based on a 4x4 transformation matrix obtained from the alignment
#' parameters. The alignment parameters include translation along x, y, z axes
#' (dx, dy, dz), and rotation angles (theta, phi, psi). The transformed
#' coordinates are calculated by applying the transformation matrix to the
#' original coordinates of the second structure (coord2).
#'
#' @param values Numeric vector. Alignment parameters obtained from structure
#'  alignment.
#'   The vector includes the following elements:\cr
#'   - dx: Translation along the x-axis.\cr
#'   - dy: Translation along the y-axis.\cr
#'   - dz: Translation along the z-axis.\cr
#'   - theta: Rotation angle around the x-axis.\cr
#'   - phi: Rotation angle around the y-axis.\cr
#'   - psi: Rotation angle around the z-axis.
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms
#'  (rows: dimensions, columns: atoms).
#'
#' @return A matrix of Euclidean distances between the transformed coordinates
#'  (coord) and the original coordinates of the first structure (coord1).
#'
#' @examples
#' \dontrun{
#' # Example: Calculate distances between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4,
#'                                byrow = TRUE
#'                                )
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7, 1, 8), nrow = 4,
#'                                   byrow = TRUE
#'                                   )
#' distances <- calculate_dist_samples(alignment_params, original_coordinates,
#'                                     transformed_coordinates
#'                                    )
#' }
#'
#' @references
#' Jordi Cruzado (https://math.stackexchange.com/users/96872/jordi-cruzado),
#' Explain 3d transformation matrix..., URL (version: 2023-05-30):
#' \href{https://math.stackexchange.com/q/532974}{Link}
#'
#' @seealso
#' \code{\link{get_matrix}} for obtaining the transformation matrix from
#'   alignment parameters.
#'
#' @export
calculate_dist_samples <- function(values, coord1, coord2) {
  if (!is.vector(values)) {
    stop("Values must be a vector.")
    }

  if (!setequal(names(values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("Values in alignment does not have the correct elements.")
    }

  if (length(dim(coord1)) != 2) {
    stop("coord1 must be a 2D matrix.")
    }

  if (dim(coord1)[1] != 4) {
    stop("The first dimension of coord1 must be 4.")
    }

  if (length(dim(coord2)) != 2) {
    stop("coord2 must be a 2D matrix.")
    }

  if (dim(coord2)[1] != 4) {
    stop("The first dimension of coord2 must be 4.")
    }

  # Get the 4x4 transformation matrix from alignment parameters
  matrix <- get_matrix(values)

  # Apply the transformation matrix to the coordinates of the second structure
  coord <- matrix %*% coord2

  # Calculate the distances between the transformed coordinates and coord1
  dist <- coord - coord1
  return(dist)
}

#' Calculate TM-Score Between Transformed Coordinates
#'
#' This function calculates the TM-Score (Template Modeling Score) between the
#' transformed coordinates and the original coordinates of two structures. The
#' TM-Score is a measure of structural similarity and is commonly used in
#' structural bioinformatics. The alignment parameters, original coordinates of
#' the first structure (coord1), original coordinates of the second structure
#' (coord2), and the square of the d0 parameter (d02) are required. Note that
#' unlike \code{\link{get_tm_samples}}, this function is a helper function that
#' performs the actual calculations to derive the TM local scores.
#'
#' @param values Numeric vector. Alignment parameters obtained from structure
#'  alignment.
#'   The vector includes the following elements:\cr
#'   - dx: Translation along the x-axis.\cr
#'   - dy: Translation along the y-axis.\cr
#'   - dz: Translation along the z-axis.\cr
#'   - theta: Rotation angle around the x-axis.\cr
#'   - phi: Rotation angle around the y-axis.\cr
#'   - psi: Rotation angle around the z-axis.
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param d02 Numeric. The square of the d0 parameter, a critical parameter in
#'  the TM-Score calculation. It represents an effective distance cutoff for
#'  defining structural similarity.
#'
#' @return The TM-Score between the transformed coordinates and the original
#'  coordinates of the structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate TM-Score between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4,
#'                                byrow = TRUE
#'                                )
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7, 1, 8), nrow = 4,
#'                                   byrow = TRUE
#'                                   )
#' d0_squared <- 5.0
#' tm_score <- calculate_tm_samples(alignment_params, original_coordinates,
#'                                  transformed_coordinates, d0_squared
#'                                  )
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{calculate_dist_samples}} for calculating distances between
#' transformed coordinates.
#'
#' @export
calculate_tm_samples <- function(values, coord1, coord2, d02) {
  if (!is.vector(values)) {
    stop("Values must be a vector.")
    }

  if (!setequal(names(values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("Values in alignment does not have the correct elements.")
    }

  if (length(dim(coord1)) != 2) {
    stop("coord1 must be a 2D matrix.")
    }

  if (dim(coord1)[1] != 4) {
    stop("The first dimension of coord1 must be 4.")
    }

  if (length(dim(coord2)) != 2) {
    stop("coord2 must be a 2D matrix.")
    }

  if (dim(coord2)[1] != 4) {
    stop("The first dimension of coord2 must be 4.")
    }

  if (!is.numeric(d02)) {
    stop("d02 must be a number.")
    }

  # Calculate distances between transformed coordinates
  dist <- calculate_dist_samples(values, coord1, coord2)

  # Calculate the sum of squared distances
  d_i2 <- colSums(dist^2)

  # Calculate the TM-Score
  tm <- 1 / (1 + d_i2 / d02)
  return(tm)
}

#' Calculate Average TM-Score Between Transformed Coordinates
#'
#' This function calculates the average TM-Score (Template Modeling Score)
#' between the transformed coordinates and the original coordinates of two
#' structures. The TM-Score is a measure of structural similarity and is
#' commonly used in structural bioinformatics. The alignment parameters,
#' original coordinates of the first structure (coord1), original coordinates of
#' the second structure (coord2), and the square of the d0 parameter (d02) are
#' required. Note that unlike \code{\link{get_tmscore}}, this function is a
#' helper function that performs the actual TM-Score calculation.
#'
#' @param values Numeric vector. Alignment parameters obtained from structure
#'  alignment.
#'   The vector includes the following elements:\cr
#'   - dx: Translation along the x-axis.\cr
#'   - dy: Translation along the y-axis.\cr
#'   - dz: Translation along the z-axis.\cr
#'   - theta: Rotation angle around the x-axis.\cr
#'   - phi: Rotation angle around the y-axis.\cr
#'   - psi: Rotation angle around the z-axis.
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param d02 Numeric. The square of the d0 parameter, a critical parameter in
#'  the TM-Score calculation. It represents an effective distance cutoff for
#'  defining structural similarity.
#'
#' @return The average TM-Score between the transformed coordinates and the
#'  original coordinates of the structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate average TM-Score between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4,
#'                                byrow = TRUE
#'                                )
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7, 1, 8), nrow = 4,
#'                                   byrow = TRUE
#'                                   )
#' d0_squared <- 5.0
#' avg_tm_score <- calculate_tmscore(alignment_params, original_coordinates,
#'                                   transformed_coordinates, d0_squared
#'                                   )
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{calculate_tm_samples}} for calculating individual TM-Scores.
#'
#' @export
calculate_tmscore <- function(values, coord1, coord2, d02) {
  if (!is.vector(values)) {
    stop("Values must be a vector.")
    }

  if (!setequal(names(values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("Values in alignment does not have the correct elements.")
    }

  if (length(dim(coord1)) != 2) {
    stop("coord1 must be a 2D matrix.")
    }

  if (dim(coord1)[1] != 4) {
    stop("The first dimension of coord1 must be 4.")
    }

  if (length(dim(coord2)) != 2) {
    stop("coord2 must be a 2D matrix.")
    }

  if (dim(coord2)[1] != 4) {
    stop("The first dimension of coord2 must be 4.")
    }

  if (!is.numeric(d02)) {
    stop("d02 must be a number.")
    }

  # Calculate individual TM-Scores
  tm_scores <- calculate_tm_samples(values, coord1, coord2, d02)

  # Calculate the average TM-Score
  avg_tm_score <- mean(tm_scores)

  return(avg_tm_score)
}

#' Calculate Root Mean Square Deviation (RMSD) Between Transformed Coordinates
#'
#' This function calculates the Root Mean Square Deviation (RMSD) between the
#' transformed coordinates and the original coordinates of two structures. The
#' RMSD is a measure of the average distance between corresponding atoms in the
#' two structures after optimal superposition. The alignment parameters,
#' original coordinates of the first structure (coord1), and original
#' coordinates of the second structure (coord2) are required. Note that unlike
#' \code{\link{get_rmsd}}, this function is a helper function that performs
#' the actual RMSD calculation.
#'
#' @param values Numeric vector. Alignment parameters obtained from structure
#'  alignment.
#'   The vector includes the following elements:\cr
#'   - dx: Translation along the x-axis.\cr
#'   - dy: Translation along the y-axis.\cr
#'   - dz: Translation along the z-axis.\cr
#'   - theta: Rotation angle around the x-axis.\cr
#'   - phi: Rotation angle around the y-axis.\cr
#'   - psi: Rotation angle around the z-axis.
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms
#'  (rows: dimensions, columns: atoms).
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms
#'  (rows: dimensions, columns: atoms).
#'
#' @return The Root Mean Square Deviation (RMSD) between the transformed
#' coordinates and the original coordinates of the structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate RMSD between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4,
#'                                byrow = TRUE
#'                                )
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7, 1, 8), nrow = 4,
#'                                   byrow = TRUE
#'                                   )
#' rmsd_value <- calculate_rmsd(alignment_params, original_coordinates,
#'                              transformed_coordinates
#'                              )
#' }
#'
#' @seealso
#' \code{\link{calculate_dist_samples}} for calculating distances between
#' transformed coordinates.
#'
#' @export
calculate_rmsd <- function(values, coord1, coord2) {
  if (!is.vector(values)) {
    stop("Values must be a vector.")
    }

  if (!setequal(names(values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("Values in alignment does not have the correct elements.")
    }

  if (length(dim(coord1)) != 2) {
    stop("coord1 must be a 2D matrix.")
    }

  if (dim(coord1)[1] != 4) {
    stop("The first dimension of coord1 must be 4.")
    }

  if (length(dim(coord2)) != 2) {
    stop("coord2 must be a 2D matrix.")
    }

  if (dim(coord2)[1] != 4) {
    stop("The first dimension of coord2 must be 4.")
    }

  # Calculate distances between transformed coordinates
  dist <- calculate_dist_samples(values, coord1, coord2)

  # Calculate the Root Mean Square Deviation (RMSD)
  rmsd_value <- sqrt(mean(dist^2))

  return(rmsd_value)
}

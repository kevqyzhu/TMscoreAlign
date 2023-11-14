# Purpose: Score alignment
# Author: Kevin Zhu
# Date: 11.14.2023
# Version: 1.0.0
# Bugs and Issues: N/A

# Define constants
DTYPE <- 'numeric'

#' Get TM Score from Protein Structure Alignment
#'
#' This function calculates the TM score from the alignment parameters and
#' coordinates obtained in a structural alignment between two protein
#' structures.
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
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
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
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
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
#' alignment_results <- get_alignment("structure1.pdb", "structure2.pdb",
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

#' Estimate d0 for Structure Alignment
#'
#' This function estimates the d0 parameter for structure alignment based on the
#' number of common residues.The formula used for estimation is derived from
#' empirical observations and is commonly used in structure alignment methods.
#' The square of the estimated d0 value (d02) is returned.
#'
#' @param N Numeric. Number of common residues between two structures.
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
  d0 <- 1.24 * (N - 15)^(1/3) - 1.8
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
#' @return A matrix of Euclidean distances between the transformed coordinates
#'  (coord) and the original
#'   coordinates of the first structure (coord1).
#'
#' @examples
#' \dontrun{
#' # Example: Calculate distances between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE
#'                                   )
#' distances <- dist_samples(alignment_params, original_coordinates,
#'                           transformed_coordinates
#'                           )
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
dist_samples <- function(values, coord1, coord2) {
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
#' (coord2), and the square of the d0 parameter (d02) are required.
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
#'  the TM-Score calculation.
#'   It represents an effective distance cutoff for defining structural
#'   similarity.
#' @return The TM-Score between the transformed coordinates and the original
#'  coordinates of the structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate TM-Score between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE
#'                                   )
#' d0_squared <- 5.0
#' tm_score <- tm_samples(alignment_params, original_coordinates,
#'                        transformed_coordinates, d0_squared)
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{dist_samples}} for calculating distances between transformed
#'   coordinates.
#'
#' @export
tm_samples <- function(values, coord1, coord2, d02) {
  # Calculate distances between transformed coordinates
  dist <- dist_samples(values, coord1, coord2)

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
#' required.
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
#'  the TM-Score calculation.
#'   It represents a distance cutoff for defining structural similarity.
#' @return The average TM-Score between the transformed coordinates and the
#'  original coordinates of the structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate average TM-Score between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE
#'                                   )
#' d0_squared <- 5.0
#' avg_tm_score <- tm(alignment_params, original_coordinates,
#'                    transformed_coordinates, d0_squared
#'                    )
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' @seealso
#' \code{\link{tm_samples}} for calculating individual TM-Scores.
#'
#' @export
tm <- function(values, coord1, coord2, d02) {
  # Calculate individual TM-Scores
  tm_scores <- tm_samples(values, coord1, coord2, d02)

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
#' coordinates of the second structure (coord2) are required.
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
#' @return The Root Mean Square Deviation (RMSD) between the transformed
#'  coordinates and the original
#'   coordinates of the structures.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate RMSD between transformed coordinates
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3
#'                       )
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE
#'                                   )
#' rmsd_value <- rmsd(alignment_params, original_coordinates,
#'                    transformed_coordinates
#'                    )
#' }
#'
#' @seealso
#' \code{\link{dist_samples}} for calculating distances between transformed
#'   coordinates.
#'
#' @export
rmsd <- function(values, coord1, coord2) {
  # Calculate distances between transformed coordinates
  dist <- dist_samples(values, coord1, coord2)

  # Calculate the Root Mean Square Deviation (RMSD)
  rmsd_value <- sqrt(mean(dist^2))

  return(rmsd_value)
}

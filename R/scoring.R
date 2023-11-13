# Define constants
DTYPE <- 'numeric'

#' Estimate d0 for Structure Alignment
#'
#' This function estimates the d0 parameter for structure alignment based on the
#' number of common residues.The formula used for estimation is derived from
#' empirical observations and is commonly used in structure alignment methods.
#' The square of the estimated d0 value (d02) is returned.
#'
#' @param N Numeric. Number of common residues between two structures.
#' @return A list containing the following component:
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


#' Get Transformation Matrix from Alignment Parameters
#'
#' This function constructs a 4x4 transformation matrix based on the alignment
#' parameters obtained from the structure alignment. The alignment parameters
#' include translation along x, y, z axes (dx, dy, dz), and rotation angles
#' (theta, phi, psi). The resulting transformation matrix can be used to
#' transform the coordinates of one structure to align with the other.
#'
#' @param values Numeric vector. Alignment parameters obtained from structure
#'  alignment.
#'   The vector includes the following elements:
#'   - dx: Translation along the x-axis.
#'   - dy: Translation along the y-axis.
#'   - dz: Translation along the z-axis.
#'   - theta: Rotation angle around the x-axis.
#'   - phi: Rotation angle around the y-axis.
#'   - psi: Rotation angle around the z-axis.
#' @return A 4x4 transformation matrix for aligning structures based on the
#'  input alignment parameters.
#'
#' @examples
#' \dontrun{
#' # Example: Get transformation matrix for structure alignment
#' alignment_params <- c(dx = 1.0, dy = 2.0, dz = 0.5, theta = 0.1, phi = 0.2,
#'                       psi = 0.3)
#' transformation_matrix <- get_matrix(alignment_params)
#' }
#'
#' @references
#' Jordi Cruzado (https://math.stackexchange.com/users/96872/jordi-cruzado),
#' Explain 3d transformation matrix..., URL (version: 2023-05-30):
#' \href{https://math.stackexchange.com/q/532974}{Link}
#'
#' @export
get_matrix <- function(values) {
  ctheta <- cos(values["theta"])
  stheta <- sin(values["theta"])
  cphi <- cos(values["phi"])
  sphi <- sin(values["phi"])
  cpsi <- cos(values["psi"])
  spsi <- sin(values["psi"])

  # Create the rotation matrix
  rotation <- matrix(0, nrow = 3, ncol = 3)
  rotation[1, 1] <- ctheta * cpsi - stheta * cphi * spsi
  rotation[1, 2] <- ctheta * spsi + stheta * cphi * cpsi
  rotation[1, 3] <- stheta * sphi
  rotation[2, 1] <- -stheta * cpsi - ctheta * cphi * spsi
  rotation[2, 2] <- -stheta * spsi + ctheta * cphi * cpsi
  rotation[2, 3] <- ctheta * sphi
  rotation[3, 1] <- sphi * spsi
  rotation[3, 2] <- -sphi * cpsi
  rotation[3, 3] <- cphi

  # Create the translation vector
  translation <- c(values["dx"], values["dy"], values["dz"])

  # Combine rotation and translation into the transformation matrix
  matrix <- cbind(rotation, translation)
  matrix <- rbind(matrix, c(0, 0, 0, 1))
  # browse()
  return(matrix)
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
#'   The vector includes the following elements:
#'   - dx: Translation along the x-axis.
#'   - dy: Translation along the y-axis.
#'   - dz: Translation along the z-axis.
#'   - theta: Rotation angle around the x-axis.
#'   - phi: Rotation angle around the y-axis.
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
#'                       psi = 0.3)
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE)
#' distances <- dist_samples(alignment_params, original_coordinates,
#'                           transformed_coordinates)
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
#'   The vector includes the following elements:
#'   - dx: Translation along the x-axis.
#'   - dy: Translation along the y-axis.
#'   - dz: Translation along the z-axis.
#'   - theta: Rotation angle around the x-axis.
#'   - phi: Rotation angle around the y-axis.
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
#'                       psi = 0.3)
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE)
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
#'   The vector includes the following elements:
#'   - dx: Translation along the x-axis.
#'   - dy: Translation along the y-axis.
#'   - dz: Translation along the z-axis.
#'   - theta: Rotation angle around the x-axis.
#'   - phi: Rotation angle around the y-axis.
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
#'                       psi = 0.3)
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE)
#' d0_squared <- 5.0
#' avg_tm_score <- tm(alignment_params, original_coordinates,
#'                    transformed_coordinates, d0_squared)
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
#'   The vector includes the following elements:
#'   - dx: Translation along the x-axis.
#'   - dy: Translation along the y-axis.
#'   - dz: Translation along the z-axis.
#'   - theta: Rotation angle around the x-axis.
#'   - phi: Rotation angle around the y-axis.
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
#'                       psi = 0.3)
#' original_coordinates <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' transformed_coordinates <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 3,
#'                                   byrow = TRUE)
#' rmsd_value <- rmsd(alignment_params, original_coordinates,
#'                    transformed_coordinates)
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

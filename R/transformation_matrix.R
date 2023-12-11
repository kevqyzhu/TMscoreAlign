# Purpose: Instantiate alignment parameters and transformation matrix
# Author: Kevin Zhu
# Date: 12.10.2023
# Version: 1.0.0
# Bugs and Issues: N/A

#' Get Default Values for Structure Alignment Parameters
#'
#' This function calculates default values for the parameters used in structure
#' alignment. The parameters include translation along x, y, z axes
#' (dx, dy, dz), and rotation angles (theta, phi, psi) based on the given 3D
#' coordinates of two structures.
#'
#' @param coord1 Matrix. 3D coordinates of the first structure's atoms
#'   (rows: dimensions, columns: atoms). This parameter should be a 4xN matrix,
#'   where N is the number of atoms. The matrix is structured such that each
#'   column corresponds to the 3D coordinates of an alpha-carbon atom in the
#'   first structure, and each row represents a dimension (x, y, z, and a
#'   placeholder coordinate). The placeholder coordinate is set to 1 for each
#'   atom, facilitating matrix transformations.
#' @param coord2 Matrix. 3D coordinates of the second structure's atoms
#'   (rows: dimensions, columns: atoms). This parameter should be a 4xN matrix,
#'   where N is the number of atoms. The matrix is structured such that each
#'   column corresponds to the 3D coordinates of an alpha-carbon atom in the
#'   second structure, and each row represents a dimension (x, y, z, and a
#'   placeholder coordinate). The placeholder coordinate is set to 1 for each
#'   atom, facilitating matrix transformations.
#'
#' @return A numeric vector containing default values for the structure
#'   alignment parameters:\cr
#'   - dx: Translation along the x-axis.\cr
#'   - dy: Translation along the y-axis.\cr
#'   - dz: Translation along the z-axis.\cr
#'   - theta: Rotation angle around the x-axis.\cr
#'   - phi: Rotation angle around the y-axis.\cr
#'   - psi: Rotation angle around the z-axis.
#'
#' @examples
#' \dontrun{
#' # Example: Get default values for structure alignment parameters
#' pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
#'                           package="TMscoreAlign"
#'                           )
#' alignment_results <- get_alignment(pdb_file1, pdb_file2,
#'                                   chain1 = 'A', chain2 = 'A',
#'                                   method = "alignment", optimize = FALSE
#'                                   )
#' coord1 <- alignment_results$coord1
#' coord2 <- alignment_results$coord2
#' default_values <- get_default_values(coord1, coord2)
#' }
#'
#' @references
#' Borchers H (2022). pracma: Practical Numerical Math Functions. R package
#' version 2.4.2, \href{https://CRAN.R-project.org/package=pracma}{Link}.
#'
#' Jur van den Berg
#' (https://math.stackexchange.com/users/91768/jur-van-den-berg),
#' Calculate Rotation Matrix to align Vector $A$ to Vector $B$ in $3D$?,
#' URL (version: 2016-09-01):
#' \href{https://math.stackexchange.com/q/476311}{Link}
#'
#' @export
#' @importFrom pracma cross
get_default_values <- function(coord1, coord2) {
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
  v <- pracma::cross(vec1, vec2)

  # Calculate rotation parameters
  s <- sqrt(sum(v^2)) + .Machine$double.eps
  c <- sum(vec1 * vec2)
  vx <- matrix(c(0, -v[3], v[2], v[3], 0, -v[1], -v[2], v[1], 0),
               nrow = 3, byrow = TRUE
               )
  rotation_matrix <- diag(3) + vx + vx%*%vx * (1 - c) / (s * s)
  values$theta <- atan2(rotation_matrix[3, 2], rotation_matrix[3, 3])
  values$phi <- atan2(-rotation_matrix[3, 1],
                      sqrt(rotation_matrix[3, 2]^2 + rotation_matrix[3, 3]^2)
                      )
  values$psi <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1])

  # Convert the list to a numeric vector
  return(unlist(values))
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
#'   The vector includes the following elements:\cr
#'   - dx: Translation along the x-axis.\cr
#'   - dy: Translation along the y-axis.\cr
#'   - dz: Translation along the z-axis.\cr
#'   - theta: Rotation angle around the x-axis.\cr
#'   - phi: Rotation angle around the y-axis.\cr
#'   - psi: Rotation angle around the z-axis.
#'
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
  if (!is.vector(values)) {
    stop("Values must be a vector.")
  }

  if (!setequal(names(values),
                c("dx", "dy", "dz", "theta", "phi", "psi"))) {
    stop("Values in alignment does not have the correct elements.")
  }

  ctheta <- cos(values["theta"])
  stheta <- sin(values["theta"])
  cphi <- cos(values["phi"])
  sphi <- sin(values["phi"])
  cpsi <- cos(values["psi"])
  spsi <- sin(values["psi"])

  # Create the rotation matrix
  rotation <- matrix(0, nrow = 3, ncol = 3)
  rotation[1, 1] <- (ctheta * cpsi) - (stheta * cphi * spsi)
  rotation[1, 2] <- (ctheta * spsi) + (stheta * cphi * cpsi)
  rotation[1, 3] <- stheta * sphi
  rotation[2, 1] <- (-stheta * cpsi) - (ctheta * cphi * spsi)
  rotation[2, 2] <- (-stheta * spsi) + (ctheta * cphi * cpsi)
  rotation[2, 3] <- ctheta * sphi
  rotation[3, 1] <- sphi * spsi
  rotation[3, 2] <- -sphi * cpsi
  rotation[3, 3] <- cphi

  # Create the translation vector
  translation <- c(values["dx"], values["dy"], values["dz"])

  # Combine rotation and translation into the transformation matrix
  matrix <- cbind(rotation, translation)
  matrix <- rbind(matrix, c(0, 0, 0, 1))

  return(matrix)
}

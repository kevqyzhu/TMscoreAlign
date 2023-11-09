library(bio3d)
library(geigen)
library(Biostrings)

# Define constants
DTYPE <- 'numeric'

# Function to estimate d0
estimate_d0 <- function(N, d0s = 5) {
  d0 <- 1.24 * (N - 15)^(1/3) - 1.8
  d02 <- d0^2
  d0s2 <- d0s^2
  return(list(d02 = d02, d0s2 = d0s2))
}

# Function to calculate the rotation-translation matrix
get_matrix <- function(theta, phi, psi, dx, dy, dz) {
  # Convert angles from degrees to radians
  theta <- theta * pi / 180
  phi <- phi * pi / 180
  psi <- psi * pi / 180

  # Calculate trigonometric values
  ctheta <- cos(theta)
  stheta <- sin(theta)
  cphi <- cos(phi)
  sphi <- sin(phi)
  cpsi <- cos(psi)
  spsi <- sin(psi)

  # Create the rotation matrix
  rotation <- matrix(0, nrow = 3, ncol = 3)
  rotation[1, 1] <- cphi * cpsi
  rotation[1, 2] <- -ctheta * spsi + stheta * sphi * cpsi
  rotation[1, 3] <- stheta * spsi + ctheta * sphi * cpsi
  rotation[2, 1] <- cphi * spsi
  rotation[2, 2] <- ctheta * cpsi + stheta * sphi * spsi
  rotation[2, 3] <- -stheta * cpsi + ctheta * sphi * spsi
  rotation[3, 1] <- -sphi
  rotation[3, 2] <- stheta * cphi
  rotation[3, 3] <- ctheta * cphi

  # Create the translation vector
  translation <- c(dx, dy, dz)

  # Combine rotation and translation into the transformation matrix
  matrix <- cbind(rotation, translation)
  matrix <- rbind(matrix, c(0, 0, 0, 1))

  return(matrix)
}

# Function to calculate the TM score
tm <- function(theta, phi, psi, dx, dy, dz, coord1, coord2, d02) {
  matrix <- get_matrix(theta, phi, psi, dx, dy, dz)
  coord <- matrix %*% coord2
  dist <- coord - coord1
  d_i2 <- rowSums(dist^2)
  tm <- -1 / (1 + d_i2 / d02)
  return(sum(tm))
}

# Function to calculate the S score
s <- function(theta, phi, psi, dx, dy, dz, coord1, coord2, d0s2) {
  matrix <- get_matrix(theta, phi, psi, dx, dy, dz)
  coord <- matrix %*% coord2
  dist <- coord - coord1
  d_i2 <- rowSums(dist^2)
  tm <- -1 / (1 + d_i2 / d0s2)
  return(sum(tm))
}

# Function to calculate the RMSD
rmsd <- function(theta, phi, psi, dx, dy, dz, coord1, coord2) {
  matrix <- get_matrix(theta, phi, psi, dx, dy, dz)
  coord <- matrix %*% coord2
  dist <- coord - coord1
  return(sum(dist^2))
}

# Function to optimize the alignment
optimise <- function(coord1, coord2, d02, d0s2, restart = TRUE) {
  if (restart) {
    default_values <- get_default_values(coord1, coord2)
  } else {
    default_values <- get_current_values()
  }

  result <- stats::optim(par = default_values, fn = tm, coord1 = coord1, coord2 = coord2, d02 = d02, method = "Nelder-Mead")
  values <- result$par
  return(list(values = values, tmscore = tm(values, coord1, coord2, d02), rmsd = sqrt(rmsd(values, coord1, coord2))))
}

# Function to get the default alignment values
get_default_values <- function(coord1, coord2) {
  values <- list()
  values$dx <- 0
  values$dy <- 0
  values$dz <- 0
  dist <- rowMeans(coord1 - coord2)
  values$dx <- dist[1]
  values$dy <- dist[2]
  values$dz <- dist[3]
  vec1 <- coord1[, 2] - coord1[, 4]
  vec2 <- coord2[, 2] - coord2[, 4]
  vec1 <- vec1 / sqrt(sum(vec1^2))
  vec2 <- vec2 / sqrt(sum(vec2^2))
  v <- crossprod(vec1, vec2)
  s <- sqrt(sum(v^2)) + .Machine$double.eps
  c <- sum(vec1 * vec2)
  rotation_matrix <- diag(3) + tcrossprod(v) + tcrossprod(tcrossprod(v), tcrossprod(v)) * (1 - c) / (s * s)
  values$theta <- atan2(rotation_matrix[3, 2], rotation_matrix[3, 3])
  values$phi <- atan2(-rotation_matrix[3, 1], sqrt(rotation_matrix[3, 2]^2 + rotation_matrix[3, 3]^2))
  values$psi <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1])
  return(values)
}

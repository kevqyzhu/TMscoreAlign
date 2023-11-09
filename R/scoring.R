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





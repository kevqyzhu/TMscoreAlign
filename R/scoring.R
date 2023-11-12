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

dist_samples <- function(values, coord1, coord2) {
  matrix <- get_matrix(values)
  coord <- matrix %*% coord2
  dist <- coord - coord1
  return(dist)
}

tm_samples <- function(values, coord1, coord2, d02) {
  dist <- dist_samples(values, coord1, coord2)
  d_i2 <- colSums(dist^2)
  tm <- 1 / (1 + d_i2 / d02)
  return(tm)
}

# Function to calculate the TM score
tm <- function(values, coord1, coord2, d02) {
  tm <- tm_samples(values, coord1, coord2, d02)
  return(mean(tm))
}

# Function to calculate the RMSD
rmsd <- function(values, coord1, coord2) {
  dist <- dist_samples(values, coord1, coord2)
  return(sqrt((mean(dist^2))))
}





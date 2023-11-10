source("R/scoring.R")
library(pracma)

# Function to load data from PDB files
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

# Function to optimize the alignment
optimise <- function(coord1, coord2, d02, d0s2, restart = TRUE) {
  if (restart) {
    default_values <- get_default_values(coord1, coord2)
  } else {
    default_values <- get_current_values()
  }
  result <- stats::optim(par = default_values, fn = tm, coord1 = coord1, coord2 = coord2, d02 = d02, method = "Nelder-Mead")
  values <- result$par
  return(list(values = values, tmscore = -tm(values, coord1, coord2, d02), rmsd = sqrt(rmsd(values, coord1, coord2))))
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
  vec1 <- coord1[-nrow(coord1), 2] - coord1[-nrow(coord1), ncol(coord1)]
  vec2 <- coord2[-nrow(coord2), 2] - coord2[-nrow(coord2), ncol(coord1)]
  vec1 <- vec1 / sqrt(sum(vec1^2))
  vec2 <- vec2 / sqrt(sum(vec2^2))
  v <- cross(vec1, vec2)
  s <- sqrt(sum(v^2)) + .Machine$double.eps
  c <- sum(vec1 * vec2)
  vx <- matrix(c(0, -v[3], v[2], v[3], 0, -v[1], -v[2], v[1], 0), nrow = 3, byrow = TRUE)
  rotation_matrix <- diag(3) + vx + vx%*%vx * (1 - c) / (s * s)
  values$theta <- atan2(rotation_matrix[3, 2], rotation_matrix[3, 3])
  values$phi <- atan2(-rotation_matrix[3, 1], sqrt(rotation_matrix[3, 2]^2 + rotation_matrix[3, 3]^2))
  values$psi <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1])
  return(unlist(values))
}

# Function to calculate the TM score between two PDB files
get_alignment <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', d0s = 5, method) {
  # Load data alignment
  data <- load_data_alignment(pdb1, pdb2, chain1, chain2, method)

  # Estimate d0
  d0_values <- estimate_d0(data$N, d0s)
  d02 <- d0_values$d02
  d0s2 <- d0_values$d0s2

  # Optimize the alignment
  alignment <- optimise(data$coord1, data$coord2, d02, d0s2)
  return(alignment)
}

get_tmscore <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', d0s = 5, method) {
  # Load data alignment
  alignment <- get_alignment(pdb1, pdb2, chain1, chain2, d0s, method)
  return(alignment$tmscore)
}

get_rmsd <- function(pdb1, pdb2, chain1 = 'A', chain2 = 'A', d0s = 5, method) {
  # Load data alignment
  alignment <- get_alignment(pdb1, pdb2, chain1, chain2, d0s, method)
  return(alignment$rmsd)
}

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

      cat(substr(line, 1, 7), sprintf("%4d", atomid), substr(line, 12, 21), "A", substr(line, 23), "\n", file = out)
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

    cat(substr(line, 1, 7), sprintf("%4d", atomid), substr(line, 12, 21), "B", substr(line, 23, 30),
        sprintf("%8.3f%8.3f%8.3f", transformed_vec[1], transformed_vec[2], transformed_vec[3]), substr(line, 55), "\n", file = out)
    atomid <- atomid + 1
  }

  close(out)
}


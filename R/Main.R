# File: Main.R

source("Aligning.R")
source("Utilities.R")

# Function to calculate TM score
# Args:
#   pdb_1: The path to the first PDB file.
#   pdb_2: The path to the second PDB file.
# Returns:
#   The TM score between the two input PDB files.
get_tm <- function(pdb_1, pdb_2) {
  align_obj <- Aligning(pdb_1, pdb_2)
  result <- align_obj$tm(0, 0, 0, 0, 0, 0)
  return(result)
}

# Function to calculate RMSD
# Args:
#   pdb_1: The path to the first PDB file.
#   pdb_2: The path to the second PDB file.
# Returns:
#   The RMSD between the two input PDB files.
get_rmsd <- function(pdb_1, pdb_2) {
  align_obj <- Aligning(pdb_1, pdb_2)
  result <- align_obj$rmsd(0, 0, 0, 0, 0, 0)
  return(result)
}

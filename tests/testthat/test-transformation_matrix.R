set.seed(1) # for exact reproducibility

# Define test cases for get_default_values function
test_that("get_default_values calculates default alignment parameters", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE)

  coord1 <- alignment$coord1
  coord2 <- alignment$coord2
  default_values <- get_default_values(coord1, coord2)

  expect_equal(as.double(default_values["dx"]), 66.96311, tolerance = 1e-6)
  expect_equal(as.double(default_values["dy"]), 14.65792, tolerance = 1e-6)
  expect_equal(as.double(default_values["dz"]), -5.640854, tolerance = 1e-6)
  expect_equal(as.double(default_values["theta"]), 1.588594, tolerance = 1e-6)
  expect_equal(as.double(default_values["phi"]), 0.7806635, tolerance = 1e-6)
  expect_equal(as.double(default_values["psi"]), 0.8305835, tolerance = 1e-6)
})

# Define test cases for get_matrix function
test_that("get_matrix constructs correct transformation matrix", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE)

  coord1 <- alignment$coord1
  coord2 <- alignment$coord2
  default_values <- get_default_values(coord1, coord2)
  matrix <- get_matrix(default_values)

  #check matrix values
  expect_equal(matrix[1, 1], -0.5364607, tolerance = 1e-6)
  expect_equal(matrix[1, 2], 0.4659413, tolerance = 1e-6)
  expect_equal(matrix[1, 3], 0.7036395, tolerance = 1e-6)
  expect_equal(matrix[1, 4], 66.96311, tolerance = 1e-6)

  expect_equal(matrix[2, 1], -0.66500293, tolerance = 1e-6)
  expect_equal(matrix[2, 2], -0.74673572, tolerance = 1e-6)
  expect_equal(matrix[2, 3], -0.01252476, tolerance = 1e-6)
  expect_equal(matrix[2, 4], 14.65791667, tolerance = 1e-6)

  expect_equal(matrix[3, 1], 0.5195970, tolerance = 1e-6)
  expect_equal(matrix[3, 2], -0.4746414, tolerance = 1e-6)
  expect_equal(matrix[3, 3], 0.7104467, tolerance = 1e-6)
  expect_equal(matrix[3, 4], -5.6408542, tolerance = 1e-6)

  expect_equal(matrix[4, 1], 0)
  expect_equal(matrix[4, 2], 0)
  expect_equal(matrix[4, 3], 0)
  expect_equal(matrix[4, 4], 1)
})

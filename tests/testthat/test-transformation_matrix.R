set.seed(1) # for exact reproducibility

# Define test cases for get_default_values function
test_that("get_default_values calculates default alignment parameters", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign"
  )
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign"
  )

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE
  )

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

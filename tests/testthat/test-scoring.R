set.seed(1) # for exact reproducibility

# Test alignment scoring functions
test_that("scoring functions return the correct values", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE)

  expect_equal(get_tmscore(alignment), 0.0133905, tolerance = 1e-6)
  expect_equal(get_rmsd(alignment), 16.67055, tolerance = 1e-6)

  optimized_alignment <- optimize_alignment(alignment, restart = TRUE)

  expect_equal(get_tmscore(optimized_alignment), 0.6054935, tolerance = 1e-6)
  expect_equal(get_rmsd(optimized_alignment), 2.639272, tolerance = 1e-6)
})

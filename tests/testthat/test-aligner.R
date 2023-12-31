set.seed(1) # for exact reproducibility

# Test get_alignment function
test_that("get_alignment returns a list", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE)

  expect_type(alignment, "list")
  expect_named(alignment, c("N", "coord1", "coord2", "values"))
  expect_type(alignment$N, "integer")
  expect_equal(length(dim(alignment$coord1)), 2)
  expect_equal(length(dim(alignment$coord2)), 2)
  expect_vector(alignment$values)
})

# Test get_alignment function shows error message
test_that("get_alignment returns error message with incorrect file input", {
  pdb_file1 <- ""
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  expect_error(alignment <- get_alignment(pdb_file1, pdb_file2,
                                          chain1 = 'A', chain2 = 'A',
                                          method = "alignment",
                                          optimize = FALSE))

  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- ""
  expect_error(alignment <- get_alignment(pdb_file1, pdb_file2,
                                          chain1 = 'A', chain2 = 'A',
                                          method = "alignment",
                                          optimize = FALSE))
})

# Test load_data_alignment function
test_that("load_data_alignment returns a list", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  # use "alignment" method for sequence alignment
  alignment_data <- load_data_alignment(pdb_file1, pdb_file2,
                                        chain1 = "A", chain2 = "A",
                                        method = "alignment")
  expect_type(alignment_data, "list")
  expect_named(alignment_data, c("coord1", "coord2", "N"))
  expect_equal(length(dim(alignment_data$coord1)), 2)
  expect_equal(length(dim(alignment_data$coord2)), 2)
  expect_type(alignment_data$N, "integer")
  expect_equal(alignment_data$N, 96)

  # use "index" method for sequence alignment
  alignment_data <- load_data_alignment(pdb_file1, pdb_file2,
                                        chain1 = "A", chain2 = "A",
                                        method = "alignment")
  expect_type(alignment_data, "list")
  expect_named(alignment_data, c("coord1", "coord2", "N"))
  expect_equal(length(dim(alignment_data$coord1)), 2)
  expect_equal(length(dim(alignment_data$coord2)), 2)
  expect_type(alignment_data$N, "integer")
  expect_equal(alignment_data$N, 96)
})

# Test optimize_alignment function
test_that("optimize_alignment returns a list", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE)
  optimized_alignment <- optimize_alignment(alignment, restart = TRUE)

  expect_type(optimized_alignment, "list")
  expect_named(optimized_alignment, c("N", "coord1", "coord2", "values"))
  expect_type(optimized_alignment$N, "integer")
  expect_equal(length(dim(optimized_alignment$coord1)), 2)
  expect_equal(length(dim(optimized_alignment$coord2)), 2)
  expect_vector(alignment$values)
})

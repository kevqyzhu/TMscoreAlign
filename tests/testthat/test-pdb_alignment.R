set.seed(1) # for exact reproducibility

# Test write_pdb function
test_that("write_pdb returns error message with incorrect file input", {
  pdb_file1 <- system.file("extdata", "1LNIA_decoy1_4.pdb",
                           package = "TMscoreAlign")
  pdb_file2 <- system.file("extdata", "1LNIA_decoy2_180.pdb",
                           package = "TMscoreAlign")

  alignment <- get_alignment(pdb_file1, pdb_file2,
                             chain1 = 'A', chain2 = 'A', method = "alignment",
                             optimize = FALSE)

  # should not be able to write to output PDB file if the input is incorrect
  expect_error(write_pdb(alignment, outputfile = "out.pdb", appended = TRUE,
                         pdb_file1, "", chain1, chain2))
  # should not be able to write to output PDB file if the input is incorrect
  expect_error(write_pdb("", outputfile = "out.pdb", appended = TRUE,
                         pdb_file1, pdb_file2, chain1, chain2))
})

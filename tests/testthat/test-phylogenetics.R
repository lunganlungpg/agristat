test_that("agristat package loads for phylogenetics tests", {
  expect_true(requireNamespace("agristat", quietly = TRUE))
})

# Placeholder tests for Module 5 — Phylogenetics & Evolutionary Mapping.
# Full tests will be added in Phase 6 when functions are implemented.

# test_that("build_phylo_tree() returns an ape phylo object", {
#   seqs <- c(
#     taxon_A = "ACGTACGT",
#     taxon_B = "ACGTATGT",
#     taxon_C = "ACGCACGT"
#   )
#   tree <- build_phylo_tree(seqs, method = "NJ")
#   expect_s3_class(tree, "phylo")
# })

# test_that("compute_distances() returns a dist matrix", {
#   seqs <- c(taxon_A = "ACGTACGT", taxon_B = "ACGTATGT")
#   d <- compute_distances(seqs)
#   expect_true(inherits(d, "dist"))
# })

# tests/testthat/test-phylogenetics.R
# Comprehensive tests for Module 5: Phylogenetics & Evolutionary Mapping
# (Steps 42-49)

# ---------------------------------------------------------------------------
# Shared test data helpers
# ---------------------------------------------------------------------------

make_dna_fasta <- function(n = 4, len = 40, file = tempfile(fileext = ".fasta")) {
  set.seed(123)
  bases   <- c("A", "T", "G", "C")
  headers <- paste0(">Seq", seq_len(n))
  seqs    <- vapply(seq_len(n), function(i) {
    paste(sample(bases, len, replace = TRUE), collapse = "")
  }, character(1L))
  lines <- unlist(mapply(c, headers, seqs, SIMPLIFY = FALSE))
  writeLines(lines, file)
  file
}

make_protein_fasta <- function(n = 3, len = 30,
                               file = tempfile(fileext = ".fasta")) {
  aa <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y")
  set.seed(42)
  headers <- paste0(">Prot", seq_len(n))
  seqs    <- vapply(seq_len(n), function(i) {
    paste(sample(aa, len, replace = TRUE), collapse = "")
  }, character(1L))
  lines <- unlist(mapply(c, headers, seqs, SIMPLIFY = FALSE))
  writeLines(lines, file)
  file
}

# ---------------------------------------------------------------------------
# Step 42: read_fasta()
# ---------------------------------------------------------------------------

test_that("read_fasta() returns correct class for DNA", {
  fasta <- make_dna_fasta(4, 40)
  seqs  <- read_fasta(fasta)
  expect_s3_class(seqs, "dna_sequences")
  expect_s3_class(seqs, "fasta_sequences")
  expect_equal(seqs$type, "nucleotide")
})

test_that("read_fasta() returns correct class for protein", {
  fasta <- make_protein_fasta(3, 30)
  seqs  <- read_fasta(fasta, type = "protein")
  expect_s3_class(seqs, "protein_sequences")
  expect_equal(seqs$type, "protein")
})

test_that("read_fasta() auto-detects protein sequences", {
  fasta <- make_protein_fasta(3, 30)
  seqs  <- read_fasta(fasta, type = "auto")
  # Protein sequences contain D, E, F, I, L, P, Q which are not in DNA set
  expect_equal(seqs$type, "protein")
})

test_that("read_fasta() parses correct number of sequences", {
  fasta <- make_dna_fasta(6, 20)
  seqs  <- read_fasta(fasta)
  expect_equal(length(seqs$sequences), 6L)
  expect_equal(nrow(seqs$stats), 6L)
})

test_that("read_fasta() has correct sequence names", {
  fasta <- make_dna_fasta(3, 20)
  seqs  <- read_fasta(fasta)
  expect_equal(names(seqs$sequences), paste0("Seq", 1:3))
})

test_that("read_fasta() reports correct lengths in stats", {
  len   <- 30L
  fasta <- make_dna_fasta(3, len)
  seqs  <- read_fasta(fasta, min_length = 5L)
  expect_true(all(seqs$stats$length == len))
})

test_that("read_fasta() computes GC content for DNA", {
  fasta_lines <- c(">SeqA", "GGGGCCCC",  # 100% GC
                   ">SeqB", "AAAATTTT")  # 0% GC
  tmp <- tempfile(fileext = ".fasta")
  writeLines(fasta_lines, tmp)
  seqs <- read_fasta(tmp, min_length = 5L)
  expect_equal(seqs$stats$gc_content[seqs$stats$name == "SeqA"], 1.0)
  expect_equal(seqs$stats$gc_content[seqs$stats$name == "SeqB"], 0.0)
})

test_that("read_fasta() GC content is NA for protein", {
  fasta <- make_protein_fasta(2, 20)
  seqs  <- read_fasta(fasta, type = "protein")
  expect_true(all(is.na(seqs$stats$gc_content)))
})

test_that("read_fasta() warns on duplicate names", {
  lines <- c(">Seq1", "ATGCATGCATGCATGCATGC",
             ">Seq1", "ATGCTTGCATGCATGCTTGC")
  tmp <- tempfile(fileext = ".fasta")
  writeLines(lines, tmp)
  expect_warning(seqs <- read_fasta(tmp), "Duplicate")
  expect_equal(length(seqs$sequences), 2L)
  expect_false(anyDuplicated(names(seqs$sequences)) > 0L)
})

test_that("read_fasta() warns on short sequences", {
  lines <- c(">Short", "ATGC", ">Long", "ATGCATGCATGCATGCATGC")
  tmp <- tempfile(fileext = ".fasta")
  writeLines(lines, tmp)
  expect_warning(seqs <- read_fasta(tmp, min_length = 10L), "min_length")
})

test_that("read_fasta() errors on missing file", {
  expect_error(read_fasta("/tmp/nonexistent_file_xyz.fasta"), "not found")
})

test_that("read_fasta() errors on empty file", {
  tmp <- tempfile()
  writeLines("", tmp)
  expect_error(read_fasta(tmp), "empty|header")
})

test_that("read_fasta() handles IUPAC ambiguity codes", {
  lines <- c(">Amb", "ATGCNRYWSMKBDHV")
  tmp <- tempfile(fileext = ".fasta")
  writeLines(lines, tmp)
  expect_no_warning(seqs <- read_fasta(tmp, min_length = 5L))
  expect_equal(seqs$type, "nucleotide")
})

test_that("read_fasta() handles multi-line sequences", {
  lines <- c(">MultiLine",
             "ATGCATGC",
             "ATGCATGC",
             "ATGCATGC")
  tmp <- tempfile(fileext = ".fasta")
  writeLines(lines, tmp)
  seqs <- read_fasta(tmp, min_length = 5L)
  expect_equal(nchar(seqs$sequences[["MultiLine"]]), 24L)
})

test_that("read_fasta() metadata is populated", {
  fasta <- make_dna_fasta(2, 20)
  seqs  <- read_fasta(fasta)
  expect_equal(seqs$metadata$n_seqs, 2L)
  expect_true(file.exists(seqs$metadata$file))
  expect_s3_class(seqs$metadata$parse_date, "POSIXct")
})

test_that("print.fasta_sequences works", {
  seqs <- read_fasta(make_dna_fasta(2, 20))
  expect_output(print(seqs), "FASTA")
})

test_that("summary.fasta_sequences works", {
  seqs <- read_fasta(make_dna_fasta(3, 30))
  expect_output(summary(seqs), "Summary")
})

# ---------------------------------------------------------------------------
# Step 43: align_sequences()
# ---------------------------------------------------------------------------

test_that("align_sequences() returns aligned_sequences class", {
  seqs <- read_fasta(make_dna_fasta(4, 20))
  aln  <- align_sequences(seqs)
  expect_s3_class(aln, "aligned_sequences")
})

test_that("align_sequences() produces correct matrix dimensions", {
  n <- 5L
  seqs <- read_fasta(make_dna_fasta(n, 20))
  aln  <- align_sequences(seqs)
  expect_equal(nrow(aln$alignment_matrix), n)
  expect_true(ncol(aln$alignment_matrix) >= 20L)
})

test_that("align_sequences() row names match sequence names", {
  seqs <- read_fasta(make_dna_fasta(4, 20))
  aln  <- align_sequences(seqs)
  expect_equal(rownames(aln$alignment_matrix), names(seqs$sequences))
})

test_that("align_sequences() stats data frame is populated", {
  seqs <- read_fasta(make_dna_fasta(4, 20))
  aln  <- align_sequences(seqs)
  expect_true(all(c("n_seqs","n_pos","gap_pct","coverage") %in%
                    names(aln$stats)))
  expect_equal(aln$stats$n_seqs, 4L)
  expect_true(aln$stats$gap_pct >= 0 & aln$stats$gap_pct <= 1)
})

test_that("align_sequences() errors on single sequence", {
  lines <- c(">Only", "ATGCATGCATGCATGCATGC")
  tmp <- tempfile(fileext = ".fasta")
  writeLines(lines, tmp)
  seqs <- read_fasta(tmp, min_length = 5L)
  expect_error(align_sequences(seqs), "At least 2")
})

test_that("align_sequences() errors on non-fasta input", {
  expect_error(align_sequences(list(a = 1)), "fasta_sequences")
})

test_that("print.aligned_sequences works", {
  seqs <- read_fasta(make_dna_fasta(3, 20))
  aln  <- align_sequences(seqs)
  expect_output(print(aln), "Aligned")
})

# ---------------------------------------------------------------------------
# Step 44: build_tree_ml()
# ---------------------------------------------------------------------------

test_that("build_tree_ml() returns phylo_tree_ml class", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_ml(aln)
  expect_s3_class(tree, "phylo_tree_ml")
  expect_s3_class(tree, "phylo_tree")
})

test_that("build_tree_ml() tree has correct number of tips", {
  skip_if_not_installed("ape")
  n    <- 5L
  seqs <- read_fasta(make_dna_fasta(n, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_ml(aln)
  expect_equal(length(tree$tree$tip.label), n)
})

test_that("build_tree_ml() tip labels match sequence names", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_ml(aln)
  expect_setequal(tree$tree$tip.label, names(seqs$sequences))
})

test_that("build_tree_ml() metadata is populated", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_ml(aln)
  expect_equal(tree$metadata$n_seqs, 4L)
  expect_equal(tree$metadata$n_positions, ncol(aln$alignment_matrix))
})

test_that("build_tree_ml() default model is GTR for DNA", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_ml(aln)
  expect_equal(tree$model, "GTR")
})

test_that("build_tree_ml() errors on non-aligned input", {
  expect_error(build_tree_ml(list(x = 1)), "aligned_sequences")
})

test_that("print.phylo_tree_ml works", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_ml(aln)
  expect_output(print(tree), "ML")
})

# ---------------------------------------------------------------------------
# Step 45: build_tree_nj()
# ---------------------------------------------------------------------------

test_that("build_tree_nj() returns phylo_tree_nj class", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  expect_s3_class(tree, "phylo_tree_nj")
  expect_s3_class(tree, "phylo_tree")
})

test_that("build_tree_nj() tree has correct number of tips", {
  skip_if_not_installed("ape")
  n    <- 6L
  seqs <- read_fasta(make_dna_fasta(n, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  expect_equal(length(tree$tree$tip.label), n)
})

test_that("build_tree_nj() distance matrix is symmetric", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  dm   <- tree$distance_matrix
  expect_true(isSymmetric(dm))
})

test_that("build_tree_nj() distance matrix diagonal is zero", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  expect_true(all(diag(tree$distance_matrix) == 0))
})

test_that("build_tree_nj() distance matrix is non-negative", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  dm   <- tree$distance_matrix
  dm[is.na(dm)] <- 0
  expect_true(all(dm >= 0))
})

test_that("build_tree_nj() accepts a distance matrix directly", {
  skip_if_not_installed("ape")
  n  <- 4L
  nm <- paste0("S", seq_len(n))
  dm <- matrix(c(0, 0.1, 0.2, 0.3,
                 0.1, 0,  0.15, 0.25,
                 0.2, 0.15, 0,  0.2,
                 0.3, 0.25, 0.2, 0),
               n, n, dimnames = list(nm, nm))
  tree <- build_tree_nj(dm)
  expect_s3_class(tree, "phylo_tree_nj")
  expect_equal(length(tree$tree$tip.label), n)
})

test_that("build_tree_nj() method is neighbor-joining", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  expect_equal(tree$method, "neighbor-joining")
})

test_that("print.phylo_tree_nj works", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  expect_output(print(tree), "NJ")
})

# ---------------------------------------------------------------------------
# Step 46: bootstrap_tree()
# ---------------------------------------------------------------------------

test_that("bootstrap_tree() returns bootstrap_result class", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  expect_s3_class(bs, "bootstrap_result")
})

test_that("bootstrap_tree() support values in 0-100 range", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 10L, seed = 42L)
  sv   <- bs$support_values$support_pct
  sv   <- sv[!is.na(sv)]
  expect_true(all(sv >= 0 & sv <= 100))
})

test_that("bootstrap_tree() consensus tree has correct tips", {
  skip_if_not_installed("ape")
  n    <- 5L
  seqs <- read_fasta(make_dna_fasta(n, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  expect_equal(length(bs$consensus_tree$tip.label), n)
})

test_that("bootstrap_tree() num_replicates is stored correctly", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 15L, seed = 1L)
  expect_equal(bs$num_replicates, 15L)
})

test_that("bootstrap_tree() with keep_replicates stores trees", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L,
                         keep_replicates = TRUE)
  expect_true(!is.null(bs$replicate_trees))
  expect_true(length(bs$replicate_trees) > 0L)
})

test_that("bootstrap_tree() replicate_trees is NULL by default", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  expect_null(bs$replicate_trees)
})

test_that("bootstrap_tree() is reproducible with same seed", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs1  <- bootstrap_tree(aln, num_replicates = 10L, seed = 99L)
  bs2  <- bootstrap_tree(aln, num_replicates = 10L, seed = 99L)
  expect_equal(bs1$support_values$support_pct,
               bs2$support_values$support_pct)
})

test_that("bootstrap_tree() errors on invalid num_replicates", {
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  expect_error(bootstrap_tree(aln, num_replicates = 0L), "num_replicates")
})

test_that("print.bootstrap_result works", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  expect_output(print(bs), "Bootstrap")
})

test_that("summary.bootstrap_result works", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  expect_output(summary(bs), "Bootstrap")
})

# ---------------------------------------------------------------------------
# Step 47: plot_phylo() and plot_tree_unrooted()
# ---------------------------------------------------------------------------

test_that("plot_phylo() runs without error (base graphics)", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  # Suppress plot output in test
  pdf(tempfile())
  expect_no_error(plot_phylo(tree))
  dev.off()
})

test_that("plot_phylo() works with bootstrap_result", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  pdf(tempfile())
  expect_no_error(plot_phylo(bs))
  dev.off()
})

test_that("plot_tree_unrooted() runs without error", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  pdf(tempfile())
  expect_no_error(plot_tree_unrooted(tree))
  dev.off()
})

# ---------------------------------------------------------------------------
# Step 48: write_tree() and reroot_tree()
# ---------------------------------------------------------------------------

test_that("write_tree() writes Newick file", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  tmp  <- tempfile(fileext = ".nwk")
  write_tree(tree, file = tmp, format = "newick")
  expect_true(file.exists(tmp))
  content <- readLines(tmp, warn = FALSE)
  expect_true(any(grepl("Seq", content)))
  # Valid Newick ends with semicolon
  combined <- paste(content, collapse = "")
  expect_true(grepl(";", combined))
})

test_that("write_tree() writes Nexus file", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  tmp  <- tempfile(fileext = ".nex")
  write_tree(tree, file = tmp, format = "nexus")
  expect_true(file.exists(tmp))
  content <- readLines(tmp, warn = FALSE)
  expect_true(any(grepl("#NEXUS", content, ignore.case = TRUE)))
})

test_that("write_tree() returns tree string invisibly", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  tmp  <- tempfile(fileext = ".nwk")
  str  <- write_tree(tree, file = tmp, format = "newick")
  expect_type(str, "character")
  expect_true(nchar(str) > 0L)
})

test_that("reroot_tree() returns an ape::phylo object", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  re   <- reroot_tree(tree, outgroup = "Seq1")
  expect_s3_class(re, "phylo")
})

test_that("reroot_tree() outgroup becomes a tip", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  re   <- reroot_tree(tree, outgroup = "Seq2")
  expect_true("Seq2" %in% re$tip.label)
})

test_that("reroot_tree() errors on unknown outgroup", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  expect_error(reroot_tree(tree, outgroup = "UnknownSeq"), "not found")
})

test_that("reroot_tree() accepts multiple outgroup tips (MRCA)", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 30))
  aln  <- align_sequences(seqs)
  tree <- build_tree_nj(aln)
  re   <- reroot_tree(tree, outgroup = c("Seq1", "Seq2"))
  expect_s3_class(re, "phylo")
})

# ---------------------------------------------------------------------------
# Step 49: validate_sequence_format(), compute_distance_matrix(),
#          extract_bootstrap_values()
# ---------------------------------------------------------------------------

test_that("validate_sequence_format() passes for valid DNA", {
  seqs <- c(S1 = "ATGCATGCATGCNRYW", S2 = "GCATGCATGCATGCNN")
  res  <- validate_sequence_format(seqs, type = "nucleotide")
  expect_true(all(res$valid))
})

test_that("validate_sequence_format() fails for illegal DNA chars", {
  seqs <- c(OK = "ATGCATGC", Bad = "ATGZATGC")
  res  <- validate_sequence_format(seqs, type = "nucleotide")
  expect_true(res$valid[res$name == "OK"])
  expect_false(res$valid[res$name == "Bad"])
  expect_true(grepl("Z", res$illegal_chars[res$name == "Bad"]))
})

test_that("validate_sequence_format() passes for valid protein", {
  seqs <- c(P1 = "ACDEFGHIKLMNPQRSTVWY", P2 = "ACDEFGHIKLM")
  res  <- validate_sequence_format(seqs, type = "protein")
  expect_true(all(res$valid))
})

test_that("validate_sequence_format() returns data frame with correct cols", {
  seqs <- c(A = "ATGC")
  res  <- validate_sequence_format(seqs)
  expect_true(all(c("name", "valid", "illegal_chars") %in% names(res)))
})

test_that("compute_distance_matrix() returns symmetric matrix", {
  fasta <- make_dna_fasta(4, 30)
  seqs  <- read_fasta(fasta)
  aln   <- align_sequences(seqs)
  dm    <- compute_distance_matrix(aln$alignment_matrix, "nucleotide", "JC69")
  expect_true(isSymmetric(dm))
})

test_that("compute_distance_matrix() diagonal is zero", {
  fasta <- make_dna_fasta(4, 30)
  seqs  <- read_fasta(fasta)
  aln   <- align_sequences(seqs)
  dm    <- compute_distance_matrix(aln$alignment_matrix, "nucleotide", "JC69")
  expect_true(all(diag(dm) == 0))
})

test_that("compute_distance_matrix() non-negative off-diagonal", {
  fasta <- make_dna_fasta(4, 30)
  seqs  <- read_fasta(fasta)
  aln   <- align_sequences(seqs)
  dm    <- compute_distance_matrix(aln$alignment_matrix, "nucleotide", "JC69")
  off_diag <- dm[lower.tri(dm)]
  off_diag[is.na(off_diag)] <- 0
  expect_true(all(off_diag >= 0))
})

test_that("compute_distance_matrix() identical seqs give distance 0", {
  mat <- matrix(rep(c("A","T","G","C"), 2), nrow = 2, byrow = TRUE,
                dimnames = list(c("s1","s2"), NULL))
  dm  <- compute_distance_matrix(mat, "nucleotide", "JC69")
  expect_equal(dm["s1","s2"], 0)
})

test_that("extract_bootstrap_values() returns named numeric vector", {
  skip_if_not_installed("ape")
  seqs <- read_fasta(make_dna_fasta(4, 40))
  aln  <- align_sequences(seqs)
  bs   <- bootstrap_tree(aln, num_replicates = 5L, seed = 1L)
  sv   <- extract_bootstrap_values(bs)
  expect_type(sv, "double")
  expect_true(length(sv) > 0L)
  expect_false(is.null(names(sv)))
})

test_that("extract_bootstrap_values() errors on non-bootstrap input", {
  expect_error(extract_bootstrap_values(list(a = 1)), "bootstrap_result")
})

# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

test_that("read_fasta() handles gaps in sequences", {
  lines <- c(">GappedSeq", "ATGC-ATGC-ATGC-ATGC")
  tmp <- tempfile(fileext = ".fasta")
  writeLines(lines, tmp)
  seqs <- read_fasta(tmp, min_length = 5L)
  expect_equal(seqs$type, "nucleotide")
  expect_true(nchar(seqs$sequences[["GappedSeq"]]) > 0L)
})

test_that("Full pipeline: FASTA -> align -> NJ -> bootstrap -> export", {
  skip_if_not_installed("ape")
  fasta <- make_dna_fasta(4, 50)
  seqs  <- read_fasta(fasta)
  aln   <- align_sequences(seqs)
  tree  <- build_tree_nj(aln)
  bs    <- bootstrap_tree(aln, num_replicates = 5L, seed = 7L)
  tmp   <- tempfile(fileext = ".nwk")
  str   <- write_tree(bs, file = tmp, format = "newick")

  expect_s3_class(seqs, "fasta_sequences")
  expect_s3_class(aln,  "aligned_sequences")
  expect_s3_class(tree, "phylo_tree_nj")
  expect_s3_class(bs,   "bootstrap_result")
  expect_true(file.exists(tmp))
  expect_type(str, "character")
})

test_that("Full pipeline: FASTA -> align -> ML -> reroot", {
  skip_if_not_installed("ape")
  fasta <- make_dna_fasta(4, 40)
  seqs  <- read_fasta(fasta)
  aln   <- align_sequences(seqs)
  tree  <- build_tree_ml(aln)
  re    <- reroot_tree(tree, outgroup = "Seq1")

  expect_s3_class(tree, "phylo_tree_ml")
  expect_s3_class(re,   "phylo")
  expect_equal(length(re$tip.label), 4L)
})

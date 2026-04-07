#' @title Module 5: Phylogenetics & Evolutionary Mapping
#'
#' @description
#' Functions for phylogenetic analysis of DNA and protein sequence data,
#' including FASTA parsing, multiple sequence alignment, tree construction
#' (Maximum Likelihood and Neighbor-Joining), bootstrap resampling,
#' tree visualization, and export to standard formats (Newick/Nexus).
#'
#' @name phylogenetics
#' @keywords internal
"_PACKAGE"

# ---------------------------------------------------------------------------
# Step 42: read_fasta()
# ---------------------------------------------------------------------------

#' Read and Validate FASTA Sequences
#'
#' Parses a FASTA-format file containing nucleotide or protein sequences.
#' Automatically detects sequence type, validates characters, and returns
#' an S3 object with per-sequence statistics.
#'
#' @param file Character. Path to a FASTA file.
#' @param type Character. One of `"auto"` (default), `"nucleotide"`, or
#'   `"protein"`.
#' @param min_length Integer. Minimum acceptable sequence length (default 10).
#'
#' @return S3 object of class `dna_sequences` or `protein_sequences`:
#'   \describe{
#'     \item{sequences}{Named character vector.}
#'     \item{type}{`"nucleotide"` or `"protein"`.}
#'     \item{stats}{Data frame: `name`, `length`, `gc_content`.}
#'     \item{metadata}{List: `file`, `n_seqs`, `parse_date`.}
#'   }
#'
#' @details
#' **DNA alphabet**: A C G T plus IUPAC codes N R Y W S K M B D H V `-`.
#' **Protein alphabet**: standard 20 amino acids plus X Z B J U O `-`.
#'
#' @references
#' IUPAC codes: \url{https://www.bioinformatics.org/sms/iupac.html}
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1 Arabidopsis thaliana",
#'   "ATGCATGCATGCATGCATGC",
#'   ">Seq2 Oryza sativa",
#'   "ATGCTTGCATGCATGCTTGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' print(seqs)
#'
#' @export
read_fasta <- function(file, type = "auto", min_length = 10L) {
  if (!file.exists(file)) rlang::abort(paste0("File not found: ", file))
  type <- match.arg(type, c("auto", "nucleotide", "protein"))

  lines <- readLines(file, warn = FALSE)
  lines <- trimws(gsub("\r", "", lines))
  lines <- lines[nchar(lines) > 0L]

  if (length(lines) == 0L) rlang::abort("FASTA file is empty.")

  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0L) rlang::abort("No FASTA headers ('>') found.")

  n_seqs    <- length(header_idx)
  seq_names <- sub("^>\\s*", "", lines[header_idx])

  if (anyDuplicated(seq_names)) {
    dups <- seq_names[duplicated(seq_names)]
    rlang::warn(paste0(
      "Duplicate sequence names: ",
      paste(unique(dups), collapse = ", "),
      ". Names made unique."
    ))
    seq_names <- make.unique(seq_names, sep = "_")
  }

  end_idx <- c(header_idx[-1L] - 1L, length(lines))
  sequences <- vapply(seq_len(n_seqs), function(i) {
    body <- lines[(header_idx[i] + 1L):end_idx[i]]
    body <- body[!grepl("^>", body)]
    toupper(paste(body, collapse = ""))
  }, character(1L))
  names(sequences) <- seq_names

  seq_lengths <- nchar(sequences)
  if (any(seq_lengths < min_length)) {
    rlang::warn(paste0(
      "Sequences shorter than min_length (", min_length, "): ",
      paste(seq_names[seq_lengths < min_length], collapse = ", ")
    ))
  }

  detected_type <- .detect_sequence_type(sequences)
  if (type == "auto") {
    type <- detected_type
  } else if (type != detected_type) {
    rlang::warn(paste0(
      "Specified type '", type, "' may not match detected type '",
      detected_type, "'. Proceeding with specified type."
    ))
  }

  .validate_sequence_chars(sequences, type)

  gc_content <- if (type == "nucleotide") {
    vapply(sequences, function(s) {
      gc <- nchar(gsub("[^GCgc]", "", s))
      if (nchar(s) == 0L) NA_real_ else round(gc / nchar(s), 4L)
    }, numeric(1L))
  } else {
    rep(NA_real_, n_seqs)
  }

  stats_df <- data.frame(
    name       = seq_names,
    length     = seq_lengths,
    gc_content = gc_content,
    stringsAsFactors = FALSE
  )

  obj <- list(
    sequences = sequences,
    type      = type,
    stats     = stats_df,
    metadata  = list(
      file       = normalizePath(file),
      n_seqs     = n_seqs,
      parse_date = Sys.time()
    )
  )
  class(obj) <- c(
    if (type == "nucleotide") "dna_sequences" else "protein_sequences",
    "fasta_sequences"
  )
  obj
}

#' @export
print.fasta_sequences <- function(x, ...) {
  cat(sprintf(
    "FASTA Sequences (%s)\n  %d sequence(s)  |  type: %s\n",
    basename(x$metadata$file), length(x$sequences), x$type
  ))
  cat("  Lengths:", paste(x$stats$length, collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.fasta_sequences <- function(object, ...) {
  cat("FASTA Sequences Summary\n")
  cat(sprintf("  File   : %s\n", object$metadata$file))
  cat(sprintf("  Type   : %s\n", object$type))
  cat(sprintf("  Seqs   : %d\n", object$metadata$n_seqs))
  cat(sprintf("  Length : %d - %d (mean %.1f)\n",
              min(object$stats$length), max(object$stats$length),
              mean(object$stats$length)))
  if (object$type == "nucleotide") {
    cat(sprintf("  GC%%    : %.1f%% - %.1f%%\n",
                min(object$stats$gc_content, na.rm = TRUE) * 100,
                max(object$stats$gc_content, na.rm = TRUE) * 100))
  }
  invisible(object)
}

# ---------------------------------------------------------------------------
# Step 43: align_sequences()
# ---------------------------------------------------------------------------

#' Align Multiple Sequences
#'
#' Performs multiple sequence alignment using ClustalOmega or MUSCLE via the
#' **msa** package when available, otherwise pads to maximum length.
#'
#' @param sequences A `dna_sequences` or `protein_sequences` object.
#' @param method Character. `"clustal"` (default) or `"muscle"`.
#' @param ... Extra arguments forwarded to `msa::msa()`.
#'
#' @return S3 object of class `aligned_sequences`:
#'   \describe{
#'     \item{alignment}{Raw alignment from `msa::msa()` or character matrix.}
#'     \item{alignment_matrix}{Character matrix (sequences x positions).}
#'     \item{method}{Alignment method.}
#'     \item{stats}{Data frame: `n_seqs`, `n_pos`, `gap_pct`, `coverage`.}
#'     \item{type}{Sequence type.}
#'   }
#'
#' @references
#' Sievers & Higgins (2018). Clustal Omega for making accurate alignments of
#' many protein sequences. *Protein Science*, 27(1), 135-145.
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' dim(aln$alignment_matrix)
#'
#' @export
align_sequences <- function(sequences, method = "clustal", ...) {
  if (!inherits(sequences, "fasta_sequences")) {
    rlang::abort("'sequences' must be a fasta_sequences object from read_fasta().")
  }
  method <- match.arg(method, c("clustal", "muscle"))
  seqs   <- sequences$sequences
  n      <- length(seqs)
  if (n < 2L) rlang::abort("At least 2 sequences required for alignment.")

  aln_mat <- NULL
  aln_obj <- NULL

  if (requireNamespace("msa", quietly = TRUE)) {
    seq_type   <- if (sequences$type == "nucleotide") "dna" else "protein"
    msa_method <- if (method == "clustal") "ClustalOmega" else "Muscle"
    ss <- msa::msa(inputSeqs = seqs, method = msa_method,
                   type = seq_type, ...)
    aln_obj <- ss
    aln_strs <- as.character(ss)
    aln_mat  <- do.call(rbind, strsplit(aln_strs, ""))
    rownames(aln_mat) <- names(seqs)
  } else {
    max_len <- max(nchar(seqs))
    padded  <- vapply(seqs, function(s) {
      g <- max_len - nchar(s)
      if (g > 0L) paste0(s, strrep("-", g)) else s
    }, character(1L))
    aln_mat <- do.call(rbind, strsplit(padded, ""))
    rownames(aln_mat) <- names(seqs)
    aln_obj <- aln_mat
  }

  n_pos   <- ncol(aln_mat)
  gap_pct <- round(sum(aln_mat == "-") / (n * n_pos), 4L)
  coverage <- round(mean(colMeans(aln_mat == "-") < 0.5), 4L)

  obj <- list(
    alignment        = aln_obj,
    alignment_matrix = aln_mat,
    method           = method,
    stats            = data.frame(n_seqs = n, n_pos = n_pos,
                                  gap_pct = gap_pct, coverage = coverage),
    type             = sequences$type
  )
  class(obj) <- "aligned_sequences"
  obj
}

#' @export
print.aligned_sequences <- function(x, ...) {
  cat(sprintf(
    "Aligned Sequences\n  Method : %s\n  Seqs   : %d  |  Positions: %d\n  Gaps   : %.1f%%\n",
    x$method, x$stats$n_seqs, x$stats$n_pos, x$stats$gap_pct * 100
  ))
  invisible(x)
}

#' @export
summary.aligned_sequences <- function(object, ...) {
  cat("Multiple Sequence Alignment\n")
  cat(sprintf("  Method    : %s\n", object$method))
  cat(sprintf("  Seqs      : %d\n", object$stats$n_seqs))
  cat(sprintf("  Positions : %d\n", object$stats$n_pos))
  cat(sprintf("  Gap %%     : %.1f%%\n", object$stats$gap_pct * 100))
  cat(sprintf("  Coverage  : %.1f%%\n", object$stats$coverage * 100))
  invisible(object)
}

# ---------------------------------------------------------------------------
# Step 44: build_tree_ml()
# ---------------------------------------------------------------------------

#' Construct a Maximum Likelihood Phylogenetic Tree
#'
#' Builds an ML tree from aligned sequences using **phangorn**, starting from
#' a Neighbor-Joining tree and optimising likelihood.
#'
#' @param alignment An `aligned_sequences` object.
#' @param model Character. Substitution model (default `"GTR"` for DNA,
#'   `"LG"` for protein).
#' @param optimize Logical. Run `phangorn::optim.pml()`? Default `TRUE`.
#' @param ... Additional arguments to `phangorn::pml()`.
#'
#' @return S3 object of class `phylo_tree_ml`:
#'   \describe{
#'     \item{tree}{`ape::phylo` object.}
#'     \item{likelihood}{Log-likelihood value.}
#'     \item{model}{Substitution model.}
#'     \item{edge_lengths}{Named numeric branch lengths.}
#'     \item{pml_fit}{Full `pml` object from phangorn.}
#'     \item{metadata}{List: `n_seqs`, `n_positions`, `alignment_method`.}
#'   }
#'
#' @references
#' Schliep (2011). phangorn: phylogenetic analysis in R.
#' *Bioinformatics*, 27(4), 592-593.
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGCTTGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGCATGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' tree <- build_tree_ml(aln)
#' print(tree)
#'
#' @export
build_tree_ml <- function(alignment, model = NULL, optimize = TRUE, ...) {
  if (!inherits(alignment, "aligned_sequences")) {
    rlang::abort("'alignment' must be an aligned_sequences object.")
  }
  if (is.null(model)) {
    model <- if (alignment$type == "nucleotide") "GTR" else "LG"
  }
  aln_mat  <- alignment$alignment_matrix
  n_seqs   <- nrow(aln_mat)
  n_pos    <- ncol(aln_mat)
  seq_type <- alignment$type

  pml_fit <- NULL
  tree    <- NULL

  if (requireNamespace("phangorn", quietly = TRUE) &&
      requireNamespace("ape", quietly = TRUE)) {
    phydat   <- .aln_mat_to_phydat(aln_mat, seq_type)
    dist_mat <- .compute_distance_matrix(aln_mat, seq_type, "JC69")
    nj_tree  <- ape::nj(dist_mat)
    pml_fit  <- phangorn::pml(nj_tree, data = phydat, model = model, ...)
    if (optimize) {
      pml_fit <- phangorn::optim.pml(
        pml_fit,
        optNni        = TRUE,
        optBf         = TRUE,
        optQ          = TRUE,
        optGamma      = TRUE,
        rearrangement = "NNI",
        control       = phangorn::pml.control(trace = 0)
      )
    }
    tree <- pml_fit$tree
  } else {
    rlang::warn(paste0(
      "Package 'phangorn' not available; using Neighbor-Joining as fallback. ",
      "Install phangorn for ML analysis."
    ))
    if (!requireNamespace("ape", quietly = TRUE)) {
      rlang::abort("Package 'ape' is required.")
    }
    dist_mat <- .compute_distance_matrix(aln_mat, seq_type, "JC69")
    tree     <- ape::nj(dist_mat)
  }

  edge_lengths <- stats::setNames(
    tree$edge.length, paste0("e", seq_along(tree$edge.length))
  )

  obj <- list(
    tree         = tree,
    likelihood   = if (!is.null(pml_fit)) pml_fit$logLik else NA_real_,
    model        = model,
    edge_lengths = edge_lengths,
    pml_fit      = pml_fit,
    metadata     = list(
      n_seqs           = n_seqs,
      n_positions      = n_pos,
      alignment_method = alignment$method
    )
  )
  class(obj) <- c("phylo_tree_ml", "phylo_tree")
  obj
}

#' @export
print.phylo_tree_ml <- function(x, ...) {
  cat(sprintf(
    "ML Phylogenetic Tree\n  Model  : %s\n  Tips   : %d\n  logLik : %.4f\n",
    x$model,
    length(x$tree$tip.label),
    if (!is.na(x$likelihood)) x$likelihood else 0
  ))
  invisible(x)
}

#' @export
summary.phylo_tree_ml <- function(object, ...) {
  cat("Maximum Likelihood Phylogenetic Tree\n")
  cat(sprintf("  Model       : %s\n", object$model))
  cat(sprintf("  Tips        : %d\n", length(object$tree$tip.label)))
  cat(sprintf("  log-Lik     : %.4f\n",
              if (!is.na(object$likelihood)) object$likelihood else 0))
  cat(sprintf("  Tip labels  : %s\n",
              paste(object$tree$tip.label, collapse = ", ")))
  invisible(object)
}

# ---------------------------------------------------------------------------
# Step 45: build_tree_nj()
# ---------------------------------------------------------------------------

#' Construct a Neighbor-Joining Phylogenetic Tree
#'
#' Builds an NJ tree from an aligned sequence object or distance matrix using
#' `ape::nj()`.
#'
#' @param alignment An `aligned_sequences` object **or** a symmetric numeric
#'   distance matrix with row/column names.
#' @param dist_method Character. `"JC69"` (default) or `"JTT"`.
#'
#' @return S3 object of class `phylo_tree_nj`:
#'   \describe{
#'     \item{tree}{`ape::phylo` object.}
#'     \item{distance_matrix}{Pairwise distance matrix.}
#'     \item{method}{`"neighbor-joining"`.}
#'     \item{metadata}{List: `dist_model`, `n_seqs`.}
#'   }
#'
#' @references
#' Saitou & Nei (1987). The neighbor-joining method.
#' *Mol. Biol. Evol.*, 4(4), 406-425.
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' tree <- build_tree_nj(aln)
#' print(tree)
#'
#' @export
build_tree_nj <- function(alignment, dist_method = "JC69") {
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for build_tree_nj().")
  }

  if (inherits(alignment, "aligned_sequences")) {
    aln_mat  <- alignment$alignment_matrix
    seq_type <- alignment$type
    dm       <- .compute_distance_matrix(aln_mat, seq_type, dist_method)
  } else if (is.matrix(alignment) && is.numeric(alignment)) {
    dm          <- alignment
    dist_method <- "user-supplied"
  } else {
    rlang::abort(
      "'alignment' must be an aligned_sequences object or a numeric matrix."
    )
  }

  tree <- ape::nj(dm)

  obj <- list(
    tree            = tree,
    distance_matrix = dm,
    method          = "neighbor-joining",
    metadata        = list(dist_model = dist_method, n_seqs = nrow(dm))
  )
  class(obj) <- c("phylo_tree_nj", "phylo_tree")
  obj
}

#' @export
print.phylo_tree_nj <- function(x, ...) {
  cat(sprintf(
    "NJ Phylogenetic Tree\n  Method : %s\n  Tips   : %d\n  Dist   : %s\n",
    x$method, length(x$tree$tip.label), x$metadata$dist_model
  ))
  invisible(x)
}

#' @export
summary.phylo_tree_nj <- function(object, ...) {
  cat("Neighbor-Joining Phylogenetic Tree\n")
  cat(sprintf("  Method      : %s\n", object$method))
  cat(sprintf("  Dist model  : %s\n", object$metadata$dist_model))
  cat(sprintf("  Tips        : %d\n", length(object$tree$tip.label)))
  cat(sprintf("  Tip labels  : %s\n",
              paste(object$tree$tip.label, collapse = ", ")))
  diag_vals <- diag(object$distance_matrix)
  cat(sprintf("  Diag (all 0): %s\n", all(diag_vals == 0)))
  invisible(object)
}

# ---------------------------------------------------------------------------
# Step 46: bootstrap_tree()
# ---------------------------------------------------------------------------

#' Bootstrap Support Values for a Phylogenetic Tree
#'
#' Resamples alignment columns with replacement, reconstructs a tree for each
#' bootstrap replicate, and computes clade support values (0-100 %).
#'
#' @param alignment An `aligned_sequences` object.
#' @param tree_method Character. `"nj"` (default) or `"ml"`.
#' @param num_replicates Integer. Bootstrap replicates (default 100).
#' @param seed Integer. Random seed for reproducibility (default 42).
#' @param keep_replicates Logical. Store replicate trees? Default `FALSE`.
#' @param tree_args Named list of extra arguments to the tree builder.
#'
#' @return S3 object of class `bootstrap_result`:
#'   \describe{
#'     \item{consensus_tree}{`ape::phylo` with bootstrap support in
#'       `node.label`.}
#'     \item{support_values}{Data frame: `node`, `support_pct`.}
#'     \item{num_replicates}{Integer.}
#'     \item{tree_method}{Character.}
#'     \item{replicate_trees}{List or `NULL`.}
#'   }
#'
#' @details
#' Values >= 70 % indicate moderate support; >= 95 % strong support.
#' Uses `ape::prop.clades()` to map support counts onto the original topology.
#'
#' @references
#' Felsenstein (1985). Confidence limits on phylogenies: an approach using the
#' bootstrap. *Evolution*, 39(4), 783-791.
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGCATGCTTGCATGCTTGCATGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGCTTGCATGCATGCTTGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGCATGCATGCATGCTTGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' bs   <- bootstrap_tree(aln, num_replicates = 10, seed = 1)
#' print(bs)
#'
#' @export
bootstrap_tree <- function(alignment,
                           tree_method     = "nj",
                           num_replicates  = 100L,
                           seed            = 42L,
                           keep_replicates = FALSE,
                           tree_args       = list()) {
  if (!inherits(alignment, "aligned_sequences")) {
    rlang::abort("'alignment' must be an aligned_sequences object.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for bootstrap_tree().")
  }
  tree_method    <- match.arg(tree_method, c("nj", "ml"))
  num_replicates <- as.integer(num_replicates)
  if (num_replicates < 1L) rlang::abort("num_replicates must be >= 1.")

  set.seed(seed)
  aln_mat <- alignment$alignment_matrix
  n_pos   <- ncol(aln_mat)

  orig_tree <- .build_one_tree(alignment, tree_method, tree_args)

  rep_trees <- vector("list", num_replicates)
  for (i in seq_len(num_replicates)) {
    col_idx      <- sample(seq_len(n_pos), n_pos, replace = TRUE)
    boot_aln     <- alignment
    boot_aln$alignment_matrix <- aln_mat[, col_idx, drop = FALSE]
    rep_trees[[i]] <- tryCatch(
      .build_one_tree(boot_aln, tree_method, tree_args),
      error = function(e) NULL
    )
  }
  rep_trees <- Filter(Negate(is.null), rep_trees)
  n_ok      <- length(rep_trees)

  support_counts <- ape::prop.clades(orig_tree, rep_trees, rooted = FALSE)
  support_pct    <- round(support_counts / n_ok * 100, 1)

  consensus_tree            <- orig_tree
  consensus_tree$node.label <- as.character(support_pct)

  support_df <- data.frame(
    node        = seq_along(support_pct) + length(orig_tree$tip.label),
    support_pct = support_pct,
    stringsAsFactors = FALSE
  )

  obj <- list(
    consensus_tree  = consensus_tree,
    support_values  = support_df,
    num_replicates  = num_replicates,
    tree_method     = tree_method,
    replicate_trees = if (keep_replicates) rep_trees else NULL
  )
  class(obj) <- "bootstrap_result"
  obj
}

#' @export
print.bootstrap_result <- function(x, ...) {
  sv <- x$support_values$support_pct
  sv <- sv[!is.na(sv)]
  cat(sprintf(
    "Bootstrap Result\n  Method : %s\n  Reps   : %d\n  Tips   : %d\n",
    x$tree_method, x$num_replicates,
    length(x$consensus_tree$tip.label)
  ))
  if (length(sv) > 0) {
    cat(sprintf("  Support: %.1f%% - %.1f%% (mean %.1f%%)\n",
                min(sv), max(sv), mean(sv)))
  }
  invisible(x)
}

#' @export
summary.bootstrap_result <- function(object, ...) {
  cat("Bootstrap Resampling Result\n")
  cat(sprintf("  Tree method   : %s\n", object$tree_method))
  cat(sprintf("  Replicates    : %d\n", object$num_replicates))
  cat(sprintf("  Tips          : %d\n",
              length(object$consensus_tree$tip.label)))
  sv <- object$support_values$support_pct
  sv <- sv[!is.na(sv)]
  if (length(sv) > 0) {
    cat(sprintf("  Support range : %.1f%% - %.1f%%\n", min(sv), max(sv)))
    cat(sprintf("  Mean support  : %.1f%%\n", mean(sv)))
    cat(sprintf("  Nodes >= 70%%  : %d / %d\n",
                sum(sv >= 70, na.rm = TRUE), length(sv)))
  }
  invisible(object)
}

# ---------------------------------------------------------------------------
# Step 47: plot_phylo() and plot_tree_unrooted()
# ---------------------------------------------------------------------------

#' Plot a Phylogenetic Tree
#'
#' Produces a phylogenetic tree plot. Uses **ggtree** when available;
#' otherwise uses `ape::plot.phylo()` base graphics.
#'
#' @param tree A `phylo_tree_ml`, `phylo_tree_nj`, `bootstrap_result`, or
#'   `ape::phylo` object.
#' @param main Character. Plot title. Default `"Phylogenetic Tree"`.
#' @param support_threshold Numeric. Bootstrap support threshold for
#'   highlighting clades (default 70).
#' @param tip_labels Logical. Show tip labels? Default `TRUE`.
#' @param ... Extra arguments forwarded to `ape::plot.phylo()`.
#'
#' @return ggplot2 object (ggtree) or `NULL` invisibly (base graphics).
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' tree <- build_tree_nj(aln)
#' plot_phylo(tree)
#'
#' @export
plot_phylo <- function(tree,
                       main              = "Phylogenetic Tree",
                       support_threshold = 70,
                       tip_labels        = TRUE,
                       ...) {
  phylo_obj      <- .extract_phylo(tree)
  support_labels <- NULL
  if (inherits(tree, "bootstrap_result")) {
    support_labels <- tree$consensus_tree$node.label
  } else if (!is.null(phylo_obj$node.label)) {
    support_labels <- phylo_obj$node.label
  }

  if (requireNamespace("ggtree", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggtree::ggtree(phylo_obj) +
      ggplot2::labs(title = main) +
      ggplot2::theme_bw()
    if (tip_labels) p <- p + ggtree::geom_tiplab(size = 3)
    return(p)
  }

  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for plot_phylo().")
  }

  ape::plot.phylo(phylo_obj, main = main,
                  show.tip.label = tip_labels, ...)

  if (!is.null(support_labels)) {
    support_num <- suppressWarnings(as.numeric(support_labels))
    n_tips   <- length(phylo_obj$tip.label)
    node_ids <- seq_along(support_labels) + n_tips
    high     <- !is.na(support_num) & support_num >= support_threshold
    if (any(high)) {
      ape::nodelabels(
        text  = support_labels[high], node = node_ids[high],
        frame = "circle", cex = 0.6, bg = "lightblue"
      )
    }
    low <- !is.na(support_num) & support_num < support_threshold
    if (any(low)) {
      ape::nodelabels(
        text  = support_labels[low], node = node_ids[low],
        frame = "none", cex = 0.5, col = "grey50"
      )
    }
  }
  ape::add.scale.bar()
  invisible(NULL)
}

#' Plot an Unrooted (Radial) Phylogenetic Tree
#'
#' Displays a tree in unrooted/radial layout using `ape` or ggtree.
#'
#' @param tree A `phylo_tree_ml`, `phylo_tree_nj`, `bootstrap_result`, or
#'   `ape::phylo` object.
#' @param main Character. Plot title. Default `"Unrooted Tree"`.
#' @param tip_labels Logical. Show tip labels? Default `TRUE`.
#' @param ... Extra arguments forwarded to `ape::plot.phylo()`.
#'
#' @return ggplot2 object (ggtree) or `NULL` invisibly (base graphics).
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' tree <- build_tree_nj(aln)
#' plot_tree_unrooted(tree)
#'
#' @export
plot_tree_unrooted <- function(tree,
                               main       = "Unrooted Tree",
                               tip_labels = TRUE,
                               ...) {
  phylo_obj <- .extract_phylo(tree)
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for plot_tree_unrooted().")
  }

  if (requireNamespace("ggtree", quietly = TRUE)) {
    p <- ggtree::ggtree(phylo_obj, layout = "unrooted") +
      ggplot2::labs(title = main) +
      ggplot2::theme_bw()
    if (tip_labels) p <- p + ggtree::geom_tiplab(size = 3)
    return(p)
  }

  ape::plot.phylo(phylo_obj, type = "unrooted", main = main,
                  show.tip.label = tip_labels, ...)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Step 48: write_tree() and reroot_tree()
# ---------------------------------------------------------------------------

#' Export a Phylogenetic Tree to File
#'
#' Writes a tree to Newick or Nexus format (compatible with FigTree, iTOL,
#' Mesquite).
#'
#' @param tree A `phylo_tree_ml`, `phylo_tree_nj`, `bootstrap_result`, or
#'   `ape::phylo` object.
#' @param file Character. Output file path.
#' @param format Character. `"newick"` (default) or `"nexus"`.
#' @param append Logical. Append to existing file? Default `FALSE`.
#'
#' @return The tree string (invisibly).
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGC"
#' )
#' tmp_fasta <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp_fasta)
#' seqs <- read_fasta(tmp_fasta)
#' aln  <- align_sequences(seqs)
#' tree <- build_tree_nj(aln)
#' tmp_nwk <- tempfile(fileext = ".nwk")
#' write_tree(tree, file = tmp_nwk, format = "newick")
#' cat(readLines(tmp_nwk), sep = "\n")
#'
#' @export
write_tree <- function(tree, file, format = "newick", append = FALSE) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for write_tree().")
  }
  format    <- match.arg(format, c("newick", "nexus"))
  phylo_obj <- .extract_phylo(tree)

  if (format == "newick") {
    ape::write.tree(phylo_obj, file = file, append = append)
    tree_str <- ape::write.tree(phylo_obj, file = "")
  } else {
    ape::write.nexus(phylo_obj, file = file, translate = TRUE)
    tree_str <- paste(readLines(file, warn = FALSE), collapse = "\n")
  }

  invisible(tree_str)
}

#' Reroot a Phylogenetic Tree
#'
#' Reroots a tree on a specified outgroup taxon or internal node.
#'
#' @param tree A `phylo_tree_ml`, `phylo_tree_nj`, `bootstrap_result`, or
#'   `ape::phylo` object.
#' @param outgroup Character vector of one or more tip names, **or** an
#'   integer node number.
#'
#' @return A rerooted `ape::phylo` object.
#'
#' @examples
#' fasta_lines <- c(
#'   ">Seq1", "ATGCATGCATGCATGCATGC",
#'   ">Seq2", "ATGCTTGCATGCATGCTTGC",
#'   ">Seq3", "ATGCATGCTTGCATGCATGC",
#'   ">Seq4", "ATGCATGCATGCTTGCATGC"
#' )
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(fasta_lines, tmp)
#' seqs <- read_fasta(tmp)
#' aln  <- align_sequences(seqs)
#' tree <- build_tree_nj(aln)
#' rerooted <- reroot_tree(tree, outgroup = "Seq4")
#' rerooted$tip.label
#'
#' @export
reroot_tree <- function(tree, outgroup) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for reroot_tree().")
  }
  phylo_obj  <- .extract_phylo(tree)
  tip_labels <- phylo_obj$tip.label

  if (is.numeric(outgroup)) {
    node <- as.integer(outgroup)
  } else {
    missing_tips <- setdiff(outgroup, tip_labels)
    if (length(missing_tips) > 0L) {
      rlang::abort(paste0(
        "Outgroup tip(s) not found: ",
        paste(missing_tips, collapse = ", ")
      ))
    }
    node <- if (length(outgroup) == 1L) {
      which(tip_labels == outgroup)
    } else {
      ape::getMRCA(phylo_obj, outgroup)
    }
  }

  ape::root(phylo_obj, node = node, resolve.root = TRUE)
}

# ---------------------------------------------------------------------------
# Exported utility functions (Step 49 additions)
# ---------------------------------------------------------------------------

#' Validate Sequence Format
#'
#' Checks sequences for illegal IUPAC characters and returns a summary.
#'
#' @param sequences Named character vector of sequences (upper-cased
#'   internally).
#' @param type Character. `"nucleotide"` (default) or `"protein"`.
#'
#' @return Data frame: `name`, `valid` (logical), `illegal_chars`.
#'
#' @examples
#' validate_sequence_format(
#'   c(Seq1 = "ATGCATGC", Seq2 = "ATGZATGC"),
#'   type = "nucleotide"
#' )
#'
#' @export
validate_sequence_format <- function(sequences, type = "nucleotide") {
  type <- match.arg(type, c("nucleotide", "protein"))
  valid_chars <- if (type == "nucleotide") {
    c("A","C","G","T","N","R","Y","W","S","K","M","B","D","H","V","-")
  } else {
    c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R",
      "S","T","V","W","Y","X","Z","B","J","U","O","-")
  }
  results <- lapply(names(sequences), function(nm) {
    chars   <- unique(strsplit(toupper(sequences[[nm]]), "")[[1L]])
    illegal <- setdiff(chars, valid_chars)
    data.frame(
      name          = nm,
      valid         = length(illegal) == 0L,
      illegal_chars = paste(illegal, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#' Compute a Pairwise Distance Matrix
#'
#' Wrapper around the internal distance computation, exposed for direct use.
#'
#' @param aln_mat Character matrix (sequences x positions) as returned in
#'   `aligned_sequences$alignment_matrix`.
#' @param seq_type Character. `"nucleotide"` or `"protein"`.
#' @param method Character. Distance model: `"JC69"` (default), `"raw"`, or
#'   any model supported by `phangorn::dist.ml()`.
#'
#' @return A symmetric numeric distance matrix with row and column names.
#'
#' @examples
#' mat <- matrix(c("A","T","G","C","A","T","G","C",
#'                  "A","C","G","T","A","C","G","T"),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("s1","s2"), NULL))
#' compute_distance_matrix(mat, "nucleotide", "JC69")
#'
#' @export
compute_distance_matrix <- function(aln_mat, seq_type, method = "JC69") {
  if (!is.matrix(aln_mat) || !is.character(aln_mat)) {
    rlang::abort("'aln_mat' must be a character matrix.")
  }
  seq_type <- match.arg(seq_type, c("nucleotide", "protein"))
  .compute_distance_matrix(aln_mat, seq_type, method)
}

#' Extract Bootstrap Support Values
#'
#' Returns support percentages for all internal nodes as a named numeric
#' vector.
#'
#' @param bootstrap A `bootstrap_result` object from [bootstrap_tree()].
#'
#' @return Named numeric vector (names are node numbers as characters).
#'
#' @examples
#' \dontrun{
#'   bs <- bootstrap_tree(aln, num_replicates = 10, seed = 1)
#'   extract_bootstrap_values(bs)
#' }
#'
#' @export
extract_bootstrap_values <- function(bootstrap) {
  if (!inherits(bootstrap, "bootstrap_result")) {
    rlang::abort("'bootstrap' must be a bootstrap_result object.")
  }
  sv <- bootstrap$support_values
  stats::setNames(sv$support_pct, as.character(sv$node))
}

# ---------------------------------------------------------------------------
# Internal helpers (not exported)
# ---------------------------------------------------------------------------

#' @noRd
.detect_sequence_type <- function(sequences) {
  all_chars <- toupper(paste(sequences, collapse = ""))
  all_chars <- gsub("[-[:digit:][:space:]]", "", all_chars)
  chars     <- unique(strsplit(all_chars, "")[[1L]])
  chars     <- chars[nchar(chars) > 0L]

  dna_only  <- c("A","C","G","T","N","R","Y","W","S","K","M","B","D","H","V")
  prot_only <- c("D","E","F","I","L","P","Q")

  if (all(chars %in% dna_only) && !any(chars %in% prot_only)) {
    "nucleotide"
  } else {
    "protein"
  }
}

#' @noRd
.validate_sequence_chars <- function(sequences, type) {
  valid <- if (type == "nucleotide") {
    c("A","C","G","T","N","R","Y","W","S","K","M","B","D","H","V","-")
  } else {
    c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R",
      "S","T","V","W","Y","X","Z","B","J","U","O","-")
  }
  for (nm in names(sequences)) {
    chars   <- unique(strsplit(sequences[[nm]], "")[[1L]])
    illegal <- setdiff(chars, valid)
    if (length(illegal) > 0L) {
      rlang::warn(paste0(
        "Sequence '", nm, "' has potentially invalid characters: ",
        paste(illegal, collapse = ", ")
      ))
    }
  }
}

#' @noRd
.compute_distance_matrix <- function(aln_mat, seq_type, method = "JC69") {
  n  <- nrow(aln_mat)
  nm <- rownames(aln_mat)
  dm <- matrix(0, n, n, dimnames = list(nm, nm))

  if (seq_type == "nucleotide" && method %in% c("JC69", "raw")) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        si   <- aln_mat[i, ]
        sj   <- aln_mat[j, ]
        keep <- si != "-" & sj != "-"
        if (!any(keep)) { dm[i, j] <- dm[j, i] <- NA_real_; next }
        p <- mean(si[keep] != sj[keep])
        d <- if (method == "JC69" && p < 0.75) {
          -3/4 * log(1 - 4/3 * p)
        } else if (method == "JC69") {
          Inf
        } else {
          p
        }
        dm[i, j] <- dm[j, i] <- max(d, 0, na.rm = TRUE)
      }
    }
  } else if (requireNamespace("phangorn", quietly = TRUE)) {
    phydat <- .aln_mat_to_phydat(aln_mat, seq_type)
    dm_ph  <- phangorn::dist.ml(phydat, model = method)
    dm     <- as.matrix(dm_ph)
    rownames(dm) <- colnames(dm) <- nm
  } else {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        si   <- aln_mat[i, ]
        sj   <- aln_mat[j, ]
        keep <- si != "-" & sj != "-"
        p    <- if (!any(keep)) NA_real_ else mean(si[keep] != sj[keep])
        dm[i, j] <- dm[j, i] <- p
      }
    }
  }
  diag(dm) <- 0
  dm
}

#' @noRd
.aln_mat_to_phydat <- function(aln_mat, seq_type) {
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    rlang::abort("Package 'phangorn' is required.")
  }
  type_str  <- if (seq_type == "nucleotide") "DNA" else "AA"
  seqs_list <- lapply(seq_len(nrow(aln_mat)), function(i) aln_mat[i, ])
  names(seqs_list) <- rownames(aln_mat)
  phangorn::phyDat(seqs_list, type = type_str)
}

#' @noRd
.build_one_tree <- function(alignment, method, tree_args) {
  if (method == "nj") {
    obj <- do.call(build_tree_nj, c(list(alignment = alignment), tree_args))
  } else {
    obj <- do.call(build_tree_ml, c(list(alignment = alignment), tree_args))
  }
  obj$tree
}

#' @noRd
.extract_phylo <- function(tree) {
  if (inherits(tree, "bootstrap_result")) return(tree$consensus_tree)
  if (inherits(tree, "phylo_tree"))      return(tree$tree)
  if (inherits(tree, "phylo"))           return(tree)
  rlang::abort(paste0(
    "Unsupported tree type: ", paste(class(tree), collapse = ", ")
  ))
}

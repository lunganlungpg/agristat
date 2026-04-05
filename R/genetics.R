# Module 3 -- Breeding & Genetics Analytics
#
# Implements Steps 25-33 of the agristat package:
#   Step 25 : design_alpha_lattice(), design_incomplete_block()
#   Step 26 : design_augmented()
#   Step 27 : analyze_gxe()
#   Step 28 : fit_ammi(), plot_ammi_biplot()
#   Step 29 : calc_heritability()
#   Step 30 : prepare_qtl_data(), plot_qtl_map()
#   Step 31 : run_basic_gwas(), plot_manhattan(), plot_qq_gwas()
#   S3 methods for all genetics classes

# ============================================================================
# Step 25: Incomplete Block Designs
# ============================================================================

#' Generate an Alpha-Lattice (Resolvable Incomplete Block) Design
#'
#' Constructs an alpha-lattice field design suitable for large trials with
#' many genotypes (typically > 100).  The design is resolvable: each
#' replicate contains every genotype exactly once, distributed across
#' incomplete blocks of equal size.
#'
#' An alpha-lattice is a resolvable Balanced Incomplete Block Design (BIBD)
#' where:
#' \itemize{
#'   \item \eqn{v} = number of genotypes,
#'   \item \eqn{k} = block size (must divide \eqn{v}),
#'   \item \eqn{r} = number of replicates,
#'   \item blocks per replicate = \eqn{v / k}.
#' }
#'
#' @param num_genotypes Integer.  Total number of genotypes to evaluate
#'   (>= 4).
#' @param block_size Integer.  Number of plots per incomplete block.
#'   Must divide `num_genotypes` evenly.
#' @param replicates Integer.  Number of complete replicates (>= 2).
#' @param seed Integer or `NULL`.  Random seed for reproducibility.
#'
#' @return A [data.frame()] with columns:
#'   \describe{
#'     \item{rep}{Replicate number (integer).}
#'     \item{block}{Incomplete block number within replicate (integer).}
#'     \item{plot_id}{Sequential plot identifier across the whole trial.}
#'     \item{genotype}{Genotype label (character, "G001", "G002", ...).}
#'   }
#'
#' @references
#' Patterson, H. D., & Williams, E. R. (1976). A new class of resolvable
#' incomplete block designs. *Biometrika*, 63(1), 83--92.
#'
#' @seealso [design_incomplete_block()], [design_augmented()]
#'
#' @examples
#' # 20 genotypes, block size 4, 3 replicates
#' d <- design_alpha_lattice(20L, 4L, 3L, seed = 42L)
#' head(d)
#' nrow(d) # 20 * 3 = 60
#'
#' @export
design_alpha_lattice <- function(num_genotypes, block_size, replicates,
                                 seed = NULL) {
  if (!.is_whole_number(num_genotypes) || num_genotypes < 4L)
    rlang::abort("`num_genotypes` must be an integer >= 4.")
  if (!.is_whole_number(block_size) || block_size < 2L)
    rlang::abort("`block_size` must be an integer >= 2.")
  if (num_genotypes %% block_size != 0L)
    rlang::abort("`block_size` must divide `num_genotypes` evenly.")
  if (!.is_whole_number(replicates) || replicates < 2L)
    rlang::abort("`replicates` must be an integer >= 2.")

  if (!is.null(seed)) set.seed(seed)

  num_genotypes  <- as.integer(num_genotypes)
  block_size     <- as.integer(block_size)
  replicates     <- as.integer(replicates)
  blocks_per_rep <- num_genotypes %/% block_size
  geno_labels    <- sprintf("G%03d", seq_len(num_genotypes))

  n_total   <- num_genotypes * replicates
  rep_vec   <- integer(n_total)
  block_vec <- integer(n_total)
  plot_vec  <- integer(n_total)
  geno_vec  <- character(n_total)

  plot_id <- 1L
  for (r in seq_len(replicates)) {
    geno_order <- sample(geno_labels)
    for (b in seq_len(blocks_per_rep)) {
      idx <- seq.int((b - 1L) * block_size + 1L, b * block_size)
      for (p in seq_along(idx)) {
        rep_vec[plot_id]   <- r
        block_vec[plot_id] <- b
        plot_vec[plot_id]  <- plot_id
        geno_vec[plot_id]  <- geno_order[idx[p]]
        plot_id <- plot_id + 1L
      }
    }
  }

  data.frame(
    rep      = rep_vec,
    block    = block_vec,
    plot_id  = plot_vec,
    genotype = geno_vec,
    stringsAsFactors = FALSE
  )
}

#' Generate a General Incomplete Block Design (BIBD)
#'
#' Creates a Balanced Incomplete Block Design (BIBD) for a user-supplied
#' vector of genotype names.  In a BIBD every pair of treatments appears
#' together in the same block the same number of times (\eqn{\lambda}
#' times).
#'
#' @param genotypes Character vector of genotype names.
#' @param block_size Integer.  Number of plots per block (\eqn{k}).
#'   Must be < `length(genotypes)`.
#' @param replicates Integer.  Number of times each genotype appears
#'   across all blocks (\eqn{r}).
#' @param seed Integer or `NULL`.  Random seed for reproducibility.
#'
#' @return An S3 object of class `field_design` with components:
#'   \describe{
#'     \item{layout}{Data frame with columns `block`, `plot_id`, `genotype`.}
#'     \item{design_type}{`"BIBD"`.}
#'     \item{parameters}{Named list with `v`, `k`, `r`, `b`, `lambda`.}
#'     \item{seed}{The seed used.}
#'   }
#'
#' @details
#' The function generates an approximate BIBD by randomising genotype
#' assignments to blocks while meeting the replication constraint.
#' Exact BIBD balance (equal \eqn{\lambda}) is achieved when
#' \eqn{r(k-1) = \lambda(v-1)}.
#'
#' @references
#' Cochran, W. G., & Cox, G. M. (1957). *Experimental Designs* (2nd ed.).
#' Wiley.
#'
#' @seealso [design_alpha_lattice()], [design_augmented()]
#'
#' @examples
#' genos <- paste0("G", 1:7)
#' d <- design_incomplete_block(genos, block_size = 3L, replicates = 3L,
#'                              seed = 1L)
#' d$layout
#'
#' @export
design_incomplete_block <- function(genotypes, block_size, replicates,
                                    seed = NULL) {
  if (!is.character(genotypes) || length(genotypes) < 2L)
    rlang::abort("`genotypes` must be a character vector of length >= 2.")
  if (!.is_whole_number(block_size) || block_size < 2L)
    rlang::abort("`block_size` must be an integer >= 2.")
  if (block_size >= length(genotypes))
    rlang::abort("`block_size` must be less than `length(genotypes)`.")
  if (!.is_whole_number(replicates) || replicates < 1L)
    rlang::abort("`replicates` must be a positive integer.")

  if (!is.null(seed)) set.seed(seed)

  v      <- length(genotypes)
  k      <- as.integer(block_size)
  r      <- as.integer(replicates)
  b_real <- v * r / k
  if (abs(b_real - round(b_real)) > 1e-9)
    rlang::abort(
      "v * replicates / block_size must be an integer for a valid BIBD."
    )
  b <- as.integer(round(b_real))

  all_genos <- sample(rep(genotypes, r))

  n_total  <- b * k
  blk_vec  <- integer(n_total)
  plot_vec <- integer(n_total)
  gen_vec  <- character(n_total)

  idx     <- 1L
  plot_id <- 1L
  for (blk in seq_len(b)) {
    for (pos in seq_len(k)) {
      blk_vec[idx]  <- blk
      plot_vec[idx] <- plot_id
      gen_vec[idx]  <- all_genos[(blk - 1L) * k + pos]
      idx     <- idx + 1L
      plot_id <- plot_id + 1L
    }
  }

  layout <- data.frame(
    block    = blk_vec,
    plot_id  = plot_vec,
    genotype = gen_vec,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      layout      = layout,
      design_type = "BIBD",
      parameters  = list(v = v, k = k, r = r, b = b,
                         lambda = r * (k - 1L) / (v - 1L)),
      seed        = seed
    ),
    class = c("field_design", "list")
  )
}

# ============================================================================
# Step 26: Augmented Design Generator
# ============================================================================

#' Generate an Augmented Field Design
#'
#' Creates an augmented randomised complete block design suitable for
#' evaluating a large number of unreplicated test lines alongside
#' replicated check genotypes.
#'
#' In an augmented design:
#' \itemize{
#'   \item Each **check** genotype appears in every block (replicated
#'     `replicates_checks` times, i.e., once per block when
#'     `replicates_checks` = number of blocks).
#'   \item Each **test** genotype appears in exactly one block (unreplicated).
#' }
#'
#' @param test_genotypes Character vector.  Names of unreplicated test lines.
#' @param check_genotypes Character vector.  Names of replicated check
#'   genotypes (>= 1).
#' @param replicates_checks Integer.  Number of complete blocks.  Each
#'   check genotype appears once per block.
#' @param seed Integer or `NULL`.  Random seed for reproducibility.
#'
#' @return A [data.frame()] with columns:
#'   \describe{
#'     \item{block}{Block number (integer).}
#'     \item{plot_id}{Sequential plot identifier (integer).}
#'     \item{genotype}{Genotype name (character).}
#'     \item{type}{`"check"` or `"test"`.}
#'   }
#'
#' @references
#' Federer, W. T. (1956). Augmented (or Hoonuiaku) designs.
#' *Hawaiian Planters' Record*, 55, 191--208.
#'
#' @seealso [design_alpha_lattice()], [design_incomplete_block()]
#'
#' @examples
#' tests  <- paste0("T", 1:12)
#' checks <- c("Check1", "Check2", "Check3")
#' d <- design_augmented(tests, checks, replicates_checks = 4L, seed = 7L)
#' table(d$type)
#'
#' @export
design_augmented <- function(test_genotypes, check_genotypes,
                             replicates_checks = 2L, seed = NULL) {
  if (!is.character(test_genotypes) || length(test_genotypes) < 1L)
    rlang::abort("`test_genotypes` must be a non-empty character vector.")
  if (!is.character(check_genotypes) || length(check_genotypes) < 1L)
    rlang::abort("`check_genotypes` must be a non-empty character vector.")
  if (!.is_whole_number(replicates_checks) || replicates_checks < 1L)
    rlang::abort("`replicates_checks` must be a positive integer.")

  n_blocks <- as.integer(replicates_checks)
  n_test   <- length(test_genotypes)

  if (!is.null(seed)) set.seed(seed)

  tests_per_block <- ceiling(n_test / n_blocks)
  test_order      <- sample(test_genotypes)

  block_vec <- integer(0)
  plot_vec  <- integer(0)
  geno_vec  <- character(0)
  type_vec  <- character(0)

  plot_id <- 1L
  for (b in seq_len(n_blocks)) {
    check_order <- sample(check_genotypes)
    for (chk in check_order) {
      block_vec <- c(block_vec, b)
      plot_vec  <- c(plot_vec,  plot_id)
      geno_vec  <- c(geno_vec,  chk)
      type_vec  <- c(type_vec,  "check")
      plot_id   <- plot_id + 1L
    }
    start <- (b - 1L) * tests_per_block + 1L
    end   <- min(b * tests_per_block, n_test)
    if (start <= n_test) {
      for (t in seq.int(start, end)) {
        block_vec <- c(block_vec, b)
        plot_vec  <- c(plot_vec,  plot_id)
        geno_vec  <- c(geno_vec,  test_order[t])
        type_vec  <- c(type_vec,  "test")
        plot_id   <- plot_id + 1L
      }
    }
  }

  data.frame(
    block    = block_vec,
    plot_id  = plot_vec,
    genotype = geno_vec,
    type     = type_vec,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Step 27: GxE Interaction Analyzer
# ============================================================================

#' Analyse Genotype-by-Environment (GxE) Interaction
#'
#' Fits a multi-environment trial (MET) model and decomposes phenotypic
#' variation into genotype (G), environment (E), and G*E interaction
#' components.  Stability statistics are computed for each genotype.
#'
#' The additive model fitted is:
#' \deqn{y_{ij} = \mu + G_i + E_j + (GE)_{ij} + \varepsilon_{ijk}}
#'
#' Stability is summarised via:
#' \describe{
#'   \item{Environmental Variance (EV)}{Variance of a genotype's mean
#'     across environments.}
#'   \item{GE Interaction}{Row in the two-way interaction matrix.}
#'   \item{Correlation}{Pearson correlation of genotype means across
#'     environments.}
#' }
#'
#' @param data A data frame with at least three columns: genotype,
#'   environment, and the response variable.
#' @param response Character.  Name of the numeric response column.
#' @param genotype Character.  Name of the genotype column
#'   (default `"genotype"`).
#' @param environment Character.  Name of the environment column
#'   (default `"environment"`).
#'
#' @return An S3 object of class `gxe_result` with components:
#'   \describe{
#'     \item{genotype_effects}{Named numeric vector of genotype main effects.}
#'     \item{env_effects}{Named numeric vector of environment main effects.}
#'     \item{interaction_matrix}{Matrix of G*E interaction deviations
#'       (genotypes x environments).}
#'     \item{stability}{Data frame: genotype, env_variance, mean_yield.}
#'     \item{anova_table}{ANOVA table (data frame) from the two-way model.}
#'     \item{response}{Name of the response variable.}
#'   }
#'
#' @references
#' Eberhart, S. A., & Russell, W. A. (1966). Stability parameters for
#' comparing varieties. *Crop Science*, 6(1), 36--40.
#'
#' @seealso [fit_ammi()]
#'
#' @examples
#' set.seed(1)
#' df <- expand.grid(genotype = paste0("G", 1:4),
#'                   environment = paste0("E", 1:3),
#'                   stringsAsFactors = FALSE)
#' df$yield <- rnorm(nrow(df), mean = 5, sd = 1)
#' res <- analyze_gxe(df, response = "yield")
#' names(res)
#'
#' @export
analyze_gxe <- function(data, response,
                        genotype    = "genotype",
                        environment = "environment") {
  .validate_data_frame(data)
  .check_column_exists(data, c(genotype, environment, response))

  geno_col <- data[[genotype]]
  env_col  <- data[[environment]]
  resp_col <- data[[response]]

  if (!is.numeric(resp_col))
    rlang::abort(sprintf("`%s` must be numeric.", response))

  grand_mean <- mean(resp_col, na.rm = TRUE)
  genos      <- sort(unique(geno_col))
  envs       <- sort(unique(env_col))

  cell_means <- matrix(NA_real_, nrow = length(genos), ncol = length(envs),
                       dimnames = list(genos, envs))
  for (g in genos) {
    for (e in envs) {
      sel <- geno_col == g & env_col == e
      if (any(sel, na.rm = TRUE))
        cell_means[g, e] <- mean(resp_col[sel], na.rm = TRUE)
    }
  }

  geno_means   <- rowMeans(cell_means, na.rm = TRUE)
  env_means    <- colMeans(cell_means, na.rm = TRUE)
  geno_effects <- geno_means - grand_mean
  env_effects  <- env_means  - grand_mean

  interaction_matrix <- cell_means - grand_mean -
    outer(geno_effects, rep(1, length(envs))) -
    outer(rep(1, length(genos)), env_effects)

  env_var <- apply(cell_means, 1, stats::var, na.rm = TRUE)
  stability_df <- data.frame(
    genotype     = genos,
    env_variance = env_var,
    mean_yield   = geno_means,
    stringsAsFactors = FALSE
  )

  fit_formula <- stats::as.formula(
    sprintf("`%s` ~ `%s` + `%s`", response, genotype, environment)
  )
  fit_aov <- tryCatch(
    stats::aov(fit_formula, data = data),
    error = function(e) NULL
  )
  anova_tbl <- if (!is.null(fit_aov))
    as.data.frame(summary(fit_aov)[[1L]])
  else NULL

  structure(
    list(
      genotype_effects   = geno_effects,
      env_effects        = env_effects,
      interaction_matrix = interaction_matrix,
      stability          = stability_df,
      anova_table        = anova_tbl,
      response           = response
    ),
    class = c("gxe_result", "list")
  )
}

#' @export
print.gxe_result <- function(x, ...) {
  cat("GxE Interaction Analysis\n")
  cat("Response:", x$response, "\n")
  cat("Genotypes:", length(x$genotype_effects), "\n")
  cat("Environments:", length(x$env_effects), "\n\n")
  cat("Stability statistics:\n")
  print(x$stability, ...)
  invisible(x)
}

#' @export
summary.gxe_result <- function(object, ...) {
  cat("GxE Interaction Analysis -- Summary\n")
  cat(rep("-", 40L), "\n", sep = "")
  cat("Response:", object$response, "\n")
  cat("Genotypes:", length(object$genotype_effects), "\n")
  cat("Environments:", length(object$env_effects), "\n\n")
  if (!is.null(object$anova_table)) {
    cat("ANOVA Table:\n")
    print(object$anova_table, ...)
  }
  cat("\nG*E Interaction Matrix:\n")
  print(round(object$interaction_matrix, 4L), ...)
  invisible(object)
}

# ============================================================================
# Step 28: AMMI Analysis
# ============================================================================

#' Fit an AMMI (Additive Main Effects and Multiplicative Interaction) Model
#'
#' Decomposes a two-way genotype-by-environment table using the AMMI
#' method: main effects are removed by ANOVA, and the residual interaction
#' matrix is decomposed by Singular Value Decomposition (SVD) to extract
#' principal components (IPCA axes).
#'
#' The AMMI model is:
#' \deqn{y_{ij} = \mu + G_i + E_j + \sum_{n=1}^{N} \lambda_n \gamma_{in}
#'   \delta_{jn} + \varepsilon_{ij}}
#' where \eqn{\lambda_n} are singular values, \eqn{\gamma_{in}} genotype
#' scores, and \eqn{\delta_{jn}} environment loadings.
#'
#' @param data A matrix (genotypes in rows, environments in columns) of
#'   mean trait values, OR a tidy data frame with columns for genotype,
#'   environment, and the response.  Row names of the matrix are genotype
#'   labels; column names are environment labels.
#' @param response Character.  Name of the response column if `data` is a
#'   tidy data frame (ignored if `data` is a matrix).
#' @param genotype Character.  Genotype column name for tidy input
#'   (default `"genotype"`).
#' @param environment Character.  Environment column name for tidy input
#'   (default `"environment"`).
#'
#' @return An S3 object of class `ammi_result` with components:
#'   \describe{
#'     \item{scores}{Matrix of genotype scores (genotypes x IPCA axes).}
#'     \item{loadings}{Matrix of environment loadings (environments x IPCA axes).}
#'     \item{singular_values}{Numeric vector of singular values.}
#'     \item{variance_explained}{Numeric vector of proportion of interaction
#'       variance explained by each IPCA axis.}
#'     \item{genotype_means}{Named vector of genotype means.}
#'     \item{env_means}{Named vector of environment means.}
#'     \item{grand_mean}{Grand mean.}
#'     \item{interaction_matrix}{The residual interaction matrix fed to SVD.}
#'   }
#'
#' @references
#' Gauch, H. G. (1988). Model selection and validation for yield trials
#' with interaction. *Biometrics*, 44(3), 705--715.
#'
#' Zobel, R. W., Wright, M. J., & Gauch, H. G. (1988). Statistical
#' analysis of a yield trial. *Agronomy Journal*, 80(3), 388--393.
#'
#' @seealso [plot_ammi_biplot()], [analyze_gxe()]
#'
#' @examples
#' # Matrix input
#' mat <- matrix(c(10,12,11, 9,14,13, 8,11,10), nrow = 3,
#'               dimnames = list(paste0("G", 1:3), paste0("E", 1:3)))
#' am <- fit_ammi(mat)
#' am$variance_explained
#'
#' # Tidy data frame input
#' set.seed(2)
#' df <- expand.grid(genotype = paste0("G", 1:4),
#'                   environment = paste0("E", 1:3),
#'                   stringsAsFactors = FALSE)
#' df$yield <- rnorm(nrow(df), 5, 1)
#' am2 <- fit_ammi(df, response = "yield")
#'
#' @export
fit_ammi <- function(data, response = "yield",
                     genotype    = "genotype",
                     environment = "environment") {

  if (is.matrix(data)) {
    mat <- data
  } else {
    .validate_data_frame(data)
    .check_column_exists(data, c(genotype, environment, response))
    genos <- sort(unique(data[[genotype]]))
    envs  <- sort(unique(data[[environment]]))
    mat   <- matrix(NA_real_, nrow = length(genos), ncol = length(envs),
                    dimnames = list(genos, envs))
    for (g in genos) {
      for (e in envs) {
        sel <- data[[genotype]] == g & data[[environment]] == e
        if (any(sel))
          mat[g, e] <- mean(data[[response]][sel], na.rm = TRUE)
      }
    }
  }

  if (is.null(rownames(mat)))
    rownames(mat) <- paste0("G", seq_len(nrow(mat)))
  if (is.null(colnames(mat)))
    colnames(mat) <- paste0("E", seq_len(ncol(mat)))

  grand_mean <- mean(mat, na.rm = TRUE)
  geno_means <- rowMeans(mat, na.rm = TRUE)
  env_means  <- colMeans(mat, na.rm = TRUE)
  geno_eff   <- geno_means - grand_mean
  env_eff    <- env_means  - grand_mean

  fitted_main <- outer(geno_eff, rep(1, ncol(mat))) +
                 outer(rep(1, nrow(mat)), env_eff) +
                 grand_mean
  interaction <- mat - fitted_main

  int_svd                    <- interaction
  int_svd[is.na(int_svd)]   <- 0
  svd_res                    <- svd(int_svd)
  sv                         <- svd_res$d
  total_var                  <- sum(sv^2)
  var_exp <- if (total_var > 0) sv^2 / total_var else rep(0, length(sv))

  n_axes <- min(nrow(mat) - 1L, ncol(mat) - 1L, length(sv))
  if (n_axes < 1L) n_axes <- length(sv)

  scores   <- svd_res$u[, seq_len(n_axes), drop = FALSE]
  loadings <- svd_res$v[, seq_len(n_axes), drop = FALSE]
  rownames(scores)   <- rownames(mat)
  rownames(loadings) <- colnames(mat)
  axis_names         <- paste0("IPCA", seq_len(n_axes))
  colnames(scores)   <- axis_names
  colnames(loadings) <- axis_names

  structure(
    list(
      scores             = scores,
      loadings           = loadings,
      singular_values    = sv[seq_len(n_axes)],
      variance_explained = var_exp[seq_len(n_axes)],
      genotype_means     = geno_means,
      env_means          = env_means,
      grand_mean         = grand_mean,
      interaction_matrix = interaction
    ),
    class = c("ammi_result", "list")
  )
}

#' @export
print.ammi_result <- function(x, ...) {
  cat("AMMI Analysis\n")
  cat(sprintf("  Genotypes   : %d\n", nrow(x$scores)))
  cat(sprintf("  Environments: %d\n", nrow(x$loadings)))
  cat(sprintf("  IPCA axes   : %d\n", ncol(x$scores)))
  cat("\nVariance explained per IPCA axis:\n")
  pct        <- round(x$variance_explained * 100, 2)
  names(pct) <- colnames(x$scores)
  print(pct, ...)
  invisible(x)
}

#' @export
summary.ammi_result <- function(object, ...) {
  print(object, ...)
  cat("\nGenotype main effects:\n")
  print(round(object$genotype_means - object$grand_mean, 4L), ...)
  cat("\nEnvironment main effects:\n")
  print(round(object$env_means - object$grand_mean, 4L), ...)
  invisible(object)
}

#' AMMI Biplot
#'
#' Produces a ggplot2 biplot of AMMI scores and loadings, showing genotype
#' scores and environment vectors radiating from the origin.
#'
#' @param x An `ammi_result` object returned by [fit_ammi()].
#' @param pc_x Integer.  IPCA axis for the x-axis (default `1L`).
#' @param pc_y Integer.  IPCA axis for the y-axis (default `2L`).
#' @param label_size Numeric.  Text label size (default `3`).
#' @param arrow_length Numeric.  Length multiplier for environment arrows
#'   (default `1`).
#'
#' @return A `ggplot` object.
#'
#' @seealso [fit_ammi()]
#'
#' @examples
#' set.seed(3)
#' df <- expand.grid(genotype     = paste0("G", 1:5),
#'                   environment  = paste0("E", 1:4),
#'                   stringsAsFactors = FALSE)
#' df$yield <- rnorm(nrow(df), 5, 1.5)
#' am <- fit_ammi(df, response = "yield")
#' p  <- plot_ammi_biplot(am)
#' # print(p)  # display in interactive session
#'
#' @export
plot_ammi_biplot <- function(x, pc_x = 1L, pc_y = 2L,
                             label_size = 3, arrow_length = 1) {
  if (!inherits(x, "ammi_result"))
    rlang::abort("`x` must be an `ammi_result` object.")
  n_axes <- ncol(x$scores)
  if (pc_x > n_axes || pc_y > n_axes)
    rlang::abort(sprintf(
      "Requested PC axes (%d, %d) exceed available axes (%d).",
      pc_x, pc_y, n_axes
    ))

  ax_x <- paste0("IPCA", pc_x)
  ax_y <- paste0("IPCA", pc_y)

  geno_df <- data.frame(
    label = rownames(x$scores),
    xval  = x$scores[, ax_x],
    yval  = x$scores[, ax_y],
    type  = "Genotype",
    stringsAsFactors = FALSE
  )

  sv_x <- x$singular_values[pc_x]
  sv_y <- x$singular_values[pc_y]
  env_df <- data.frame(
    label = rownames(x$loadings),
    xval  = x$loadings[, ax_x] * sv_x * arrow_length,
    yval  = x$loadings[, ax_y] * sv_y * arrow_length,
    type  = "Environment",
    stringsAsFactors = FALSE
  )

  pct_x <- round(x$variance_explained[pc_x] * 100, 1)
  pct_y <- round(x$variance_explained[pc_y] * 100, 1)

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70",
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70",
                        linetype = "dashed") +
    ggplot2::geom_segment(
      data      = env_df,
      ggplot2::aes(x = 0, y = 0, xend = .data$xval, yend = .data$yval,
                   colour = .data$type),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.25, "cm")),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = geno_df,
      ggplot2::aes(x = .data$xval, y = .data$yval, colour = .data$type),
      size = 3
    ) +
    ggplot2::geom_text(
      data = rbind(geno_df, env_df),
      ggplot2::aes(x = .data$xval, y = .data$yval,
                   label = .data$label, colour = .data$type),
      size        = label_size,
      vjust       = -0.6,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c("Genotype" = "#1f78b4", "Environment" = "#e31a1c"),
      name   = NULL
    ) +
    ggplot2::labs(
      title = "AMMI Biplot",
      x     = sprintf("%s (%.1f%%)", ax_x, pct_x),
      y     = sprintf("%s (%.1f%%)", ax_y, pct_y)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

# ============================================================================
# Step 29: Heritability Calculators
# ============================================================================

#' Compute Broad- and Narrow-Sense Heritability
#'
#' Estimates broad-sense heritability \eqn{H^2} and, optionally,
#' narrow-sense heritability \eqn{h^2} from a data frame of variance
#' components obtained from a REML fit (e.g., via [lme4::lmer()]).
#'
#' \deqn{H^2 = \frac{V_G}{V_G + V_E}}
#' \deqn{h^2 = \frac{V_A}{V_G + V_E}}
#'
#' Optional bootstrap confidence intervals are computed by resampling the
#' ratio statistic.
#'
#' @param variance_components A data frame with columns `source` (character)
#'   and `variance` (numeric).  Required sources: `"genotype"` and
#'   `"residual"` (or `"error"`).  Optional: `"additive"`.
#' @param n_boot Integer.  Number of bootstrap replicates for confidence
#'   intervals (default `1000L`; set to `0L` to skip).
#' @param ci_level Numeric.  Confidence level (default `0.95`).
#' @param seed Integer or `NULL`.  Random seed for bootstrap.
#'
#' @return An S3 object of class `heritability_result` with components:
#'   \describe{
#'     \item{H2}{Broad-sense heritability estimate.}
#'     \item{h2}{Narrow-sense heritability estimate (`NA` if additive
#'       variance not provided).}
#'     \item{ci_H2}{Two-element numeric vector: lower and upper CI for
#'       \eqn{H^2} (or `NA` if `n_boot = 0`).}
#'     \item{ci_h2}{Two-element numeric vector for \eqn{h^2} (`NA` if not
#'       computed).}
#'     \item{variance_components}{The input data frame.}
#'     \item{interpretation}{A character string with a plain-language
#'       summary.}
#'   }
#'
#' @references
#' Falconer, D. S., & Mackay, T. F. C. (1996). *Introduction to
#' Quantitative Genetics* (4th ed.). Longman.
#'
#' Holland, J. B., Nyquist, W. E., & Cervantes-Martinez, C. T. (2003).
#' Estimating and interpreting heritability for plant breeding.
#' *Plant Breeding Reviews*, 22, 9--112.
#'
#' @seealso [analyze_gxe()], [fit_ammi()]
#'
#' @examples
#' vc <- data.frame(
#'   source   = c("genotype", "residual"),
#'   variance = c(2.5, 1.0)
#' )
#' h <- calc_heritability(vc, n_boot = 0L)
#' h$H2  # 0.714...
#'
#' @export
calc_heritability <- function(variance_components, n_boot = 1000L,
                              ci_level = 0.95, seed = NULL) {
  .validate_data_frame(variance_components, "variance_components")
  .check_column_exists(variance_components, c("source", "variance"),
                       "variance_components")

  vc        <- variance_components
  vc$source <- tolower(trimws(vc$source))

  if (!"residual" %in% vc$source && "error" %in% vc$source)
    vc$source[vc$source == "error"] <- "residual"

  .require_source <- function(s) {
    if (!s %in% vc$source)
      rlang::abort(sprintf(
        "`variance_components` must contain a row with source = '%s'.", s
      ))
  }
  .require_source("genotype")
  .require_source("residual")

  v_g <- vc$variance[vc$source == "genotype"][1L]
  v_e <- vc$variance[vc$source == "residual"][1L]
  v_a <- if ("additive" %in% vc$source)
           vc$variance[vc$source == "additive"][1L]
         else NA_real_

  v_p <- v_g + v_e
  H2  <- v_g / v_p
  h2  <- if (!is.na(v_a)) v_a / v_p else NA_real_

  ci_H2 <- NA_real_
  ci_h2 <- NA_real_
  if (n_boot > 0L) {
    if (!is.null(seed)) set.seed(seed)
    alpha    <- (1 - ci_level) / 2
    boots_H2 <- replicate(n_boot, {
      bvg <- stats::rgamma(1, shape = max(v_g, 1e-9), rate = 1)
      bve <- stats::rgamma(1, shape = max(v_e, 1e-9), rate = 1)
      bvg / (bvg + bve)
    })
    ci_H2 <- stats::quantile(boots_H2, c(alpha, 1 - alpha), names = FALSE)
    if (!is.na(v_a)) {
      boots_h2 <- replicate(n_boot, {
        bva <- stats::rgamma(1, shape = max(v_a, 1e-9), rate = 1)
        bvg <- stats::rgamma(1, shape = max(v_g, 1e-9), rate = 1)
        bve <- stats::rgamma(1, shape = max(v_e, 1e-9), rate = 1)
        bva / (bva + bvg + bve)
      })
      ci_h2 <- stats::quantile(boots_h2, c(alpha, 1 - alpha), names = FALSE)
    }
  }

  structure(
    list(
      H2                  = H2,
      h2                  = h2,
      ci_H2               = ci_H2,
      ci_h2               = ci_h2,
      variance_components = variance_components,
      interpretation      = .interpret_heritability(H2)
    ),
    class = c("heritability_result", "list")
  )
}

.interpret_heritability <- function(H2) {
  if (is.na(H2)) return("Cannot interpret: missing values.")
  if (H2 >= 0.8)
    "High heritability (>= 0.80): trait is strongly genetically controlled."
  else if (H2 >= 0.5)
    "Moderate heritability (0.50-0.79): substantial genetic component."
  else if (H2 >= 0.2)
    "Low heritability (0.20-0.49): environment has large influence."
  else
    "Very low heritability (< 0.20): trait is predominantly environmentally determined."
}

#' @export
print.heritability_result <- function(x, ...) {
  cat("Heritability Estimates\n")
  cat(sprintf("  H^2 (broad-sense)  = %.4f\n", x$H2))
  if (!is.na(x$h2))
    cat(sprintf("  h^2 (narrow-sense) = %.4f\n", x$h2))
  if (length(x$ci_H2) > 1L && !anyNA(x$ci_H2))
    cat(sprintf("  95%% CI for H^2     = [%.4f, %.4f]\n",
                x$ci_H2[1L], x$ci_H2[2L]))
  cat("\nInterpretation:", x$interpretation, "\n")
  invisible(x)
}

#' @export
summary.heritability_result <- function(object, ...) {
  print(object, ...)
  cat("\nVariance components:\n")
  print(object$variance_components, ...)
  invisible(object)
}

# ============================================================================
# Step 30: QTL Mapping Helpers
# ============================================================================

#' Prepare Data for QTL Mapping
#'
#' Validates and reformats phenotype and marker genotype data for
#' downstream QTL analysis, e.g. with the **qtl** package.
#'
#' Marker coding:
#' \itemize{
#'   \item Numeric (0, 1, 2): AA, Aa, aa.
#'   \item Character (A, H, B): converted to (2, 1, 0).
#' }
#'
#' @param phenotypes A data frame with at minimum one numeric column for
#'   the trait of interest and an `id` column matching row names/IDs in
#'   `marker_genotypes`.
#' @param marker_genotypes A matrix or data frame (individuals in rows,
#'   markers in columns).  Values should be 0/1/2 or A/H/B.
#' @param map_info Optional data frame with columns `marker` and
#'   `chromosome` (and optionally `position`).
#' @param id_col Character.  Name of the individual ID column in
#'   `phenotypes` (default `"id"`).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{phenotypes}{Aligned phenotype data frame.}
#'     \item{genotypes}{Numeric marker matrix (0/1/2 coding).}
#'     \item{map_info}{Map data frame (may be auto-generated).}
#'     \item{n_individuals}{Integer.}
#'     \item{n_markers}{Integer.}
#'     \item{qc_report}{Data frame with per-marker QC statistics.}
#'   }
#'
#' @seealso [plot_qtl_map()], [run_basic_gwas()]
#'
#' @examples
#' set.seed(4)
#' n <- 50L; m <- 10L
#' pheno <- data.frame(id = paste0("ind", seq_len(n)),
#'                     yield = rnorm(n))
#' markers <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
#'                   dimnames = list(paste0("ind", seq_len(n)),
#'                                   paste0("M", seq_len(m))))
#' qtl_data <- prepare_qtl_data(pheno, markers)
#' qtl_data$n_markers
#'
#' @export
prepare_qtl_data <- function(phenotypes, marker_genotypes,
                             map_info = NULL, id_col = "id") {
  .validate_data_frame(phenotypes, "phenotypes")
  if (!is.matrix(marker_genotypes) && !is.data.frame(marker_genotypes))
    rlang::abort("`marker_genotypes` must be a matrix or data frame.")

  marker_mat <- as.matrix(marker_genotypes)
  n_ind      <- nrow(marker_mat)
  n_mrk      <- ncol(marker_mat)

  if (is.character(marker_mat)) {
    coding         <- c("A" = 2L, "H" = 1L, "B" = 0L)
    marker_mat_num <- matrix(
      coding[marker_mat],
      nrow      = n_ind,
      dimnames  = dimnames(marker_mat)
    )
    marker_mat <- marker_mat_num
  } else {
    marker_mat <- matrix(
      apply(marker_mat, 2, as.numeric),
      nrow     = n_ind,
      dimnames = dimnames(marker_mat)
    )
  }

  if (id_col %in% names(phenotypes) && !is.null(rownames(marker_mat))) {
    common_ids <- intersect(phenotypes[[id_col]], rownames(marker_mat))
    if (length(common_ids) == 0L)
      rlang::abort(
        "No matching IDs between `phenotypes` and `marker_genotypes`."
      )
    pheno_aligned  <- phenotypes[phenotypes[[id_col]] %in% common_ids, ,
                                 drop = FALSE]
    marker_aligned <- marker_mat[common_ids, , drop = FALSE]
  } else {
    pheno_aligned  <- phenotypes
    marker_aligned <- marker_mat
  }

  miss_rate  <- apply(marker_aligned, 2, function(x) mean(is.na(x)))
  minor_freq <- apply(marker_aligned, 2, function(x) {
    counts <- table(x)
    if (length(counts) == 0L) return(NA_real_)
    min(counts / sum(counts), na.rm = TRUE)
  })
  mono <- apply(marker_aligned, 2,
                function(x) length(unique(x[!is.na(x)])) <= 1L)

  qc <- data.frame(
    marker       = colnames(marker_aligned),
    missing_rate = miss_rate,
    minor_af     = minor_freq,
    monomorphic  = mono,
    stringsAsFactors = FALSE
  )

  if (is.null(map_info)) {
    map_info <- data.frame(
      marker     = colnames(marker_aligned),
      chromosome = rep(1L, n_mrk),
      position   = seq_len(n_mrk) * 10.0,
      stringsAsFactors = FALSE
    )
  }

  list(
    phenotypes    = pheno_aligned,
    genotypes     = marker_aligned,
    map_info      = map_info,
    n_individuals = nrow(marker_aligned),
    n_markers     = ncol(marker_aligned),
    qc_report     = qc
  )
}

#' Plot QTL LOD Scan Results
#'
#' Visualises the output of a QTL scan (e.g., `qtl::scanone()`) as a
#' ggplot2 plot, showing LOD scores along the genome with a horizontal
#' significance threshold line.
#'
#' @param scan_result A data frame with columns `chr` (chromosome),
#'   `pos` (cM position), and `lod` (LOD score).  Compatible with the
#'   output of `qtl::scanone()`.
#' @param threshold Numeric.  Genome-wide significance threshold for LOD
#'   scores (default `3.0`).
#' @param trait_name Character.  Title label for the plot
#'   (default `"QTL Scan"`).
#'
#' @return A `ggplot` object.
#'
#' @seealso [prepare_qtl_data()]
#'
#' @examples
#' set.seed(5)
#' scan <- data.frame(
#'   chr = rep(1:3, each = 20L),
#'   pos = rep(seq(0, 95, 5), 3L),
#'   lod = c(rexp(20, 1), rexp(20, 0.5), rexp(20, 1))
#' )
#' p <- plot_qtl_map(scan, threshold = 2.5)
#' # print(p)
#'
#' @export
plot_qtl_map <- function(scan_result, threshold = 3.0,
                         trait_name = "QTL Scan") {
  .validate_data_frame(scan_result, "scan_result")
  .check_column_exists(scan_result, c("chr", "pos", "lod"), "scan_result")

  scan_result$chr <- as.factor(scan_result$chr)
  chr_levels      <- levels(scan_result$chr)

  chrom_max    <- tapply(scan_result$pos, scan_result$chr, max)
  chrom_offset <- c(0, cumsum(chrom_max[-length(chrom_max)]))
  names(chrom_offset)  <- chr_levels
  scan_result$cum_pos  <- scan_result$pos +
    chrom_offset[as.character(scan_result$chr)]

  chr_mid     <- tapply(scan_result$cum_pos, scan_result$chr, mean)
  n_chr       <- length(chr_levels)
  chr_colours <- rep(c("#2166ac", "#4dac26"), length.out = n_chr)
  names(chr_colours) <- chr_levels

  ggplot2::ggplot(scan_result,
                  ggplot2::aes(x = .data$cum_pos, y = .data$lod,
                               colour = .data$chr)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = threshold,
                        colour = "red", linetype = "dashed",
                        linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = chr_colours, guide = "none") +
    ggplot2::scale_x_continuous(breaks = chr_mid, labels = chr_levels) +
    ggplot2::labs(
      title = trait_name,
      x     = "Chromosome",
      y     = "LOD score"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# ============================================================================
# Step 31: GWAS Pipeline
# ============================================================================

#' Run a Basic GWAS (Genome-Wide Association Study)
#'
#' Performs marker-by-marker association tests between a phenotype vector
#' and a matrix of SNP genotypes (coded 0/1/2).  Supports linear regression
#' for quantitative traits and logistic regression for binary traits.
#'
#' For each SNP \eqn{j}, the model is:
#' \deqn{y = \mu + \beta_j x_j + \varepsilon}
#' where \eqn{x_j \in \{0,1,2\}} is the allele dosage.  P-values are from
#' the Wald test; effect sizes are the regression coefficient
#' \eqn{\hat{\beta}_j}.
#'
#' @param phenotype Numeric vector (quantitative) or integer vector of 0/1
#'   (binary) trait values.  Length must equal `nrow(marker_genotypes)`.
#' @param marker_genotypes Numeric matrix (individuals x SNPs, 0/1/2
#'   coding).  Column names are used as marker IDs.
#' @param marker_positions A data frame with columns `marker`, `chromosome`,
#'   and `position`.  If `NULL`, a default layout is generated.
#' @param trait_type Character.  `"quantitative"` (default) or `"binary"`.
#' @param min_maf Numeric.  Minimum minor allele frequency filter
#'   (default `0.05`).
#'
#' @return An S3 object of class `gwas_result` with components:
#'   \describe{
#'     \item{results}{Data frame: `marker`, `chromosome`, `position`,
#'       `p_value`, `neg_log10_p`, `beta`, `se`, `maf`.}
#'     \item{trait_type}{`"quantitative"` or `"binary"`.}
#'     \item{n_individuals}{Number of individuals analysed.}
#'     \item{n_markers_tested}{Number of SNPs after MAF filter.}
#'     \item{bonferroni_threshold}{Bonferroni-corrected p-value threshold.}
#'   }
#'
#' @references
#' Visscher, P. M., Wray, N. R., Zhang, Q., et al. (2017). 10 years of
#' GWAS discovery: biology, function, and translation.
#' *American Journal of Human Genetics*, 101(1), 5--22.
#'
#' @seealso [plot_manhattan()], [plot_qq_gwas()], [prepare_qtl_data()]
#'
#' @examples
#' set.seed(6)
#' n <- 100L; m <- 50L
#' pheno <- rnorm(n)
#' gmat  <- matrix(sample(0:2, n * m, replace = TRUE,
#'                         prob = c(0.25, 0.5, 0.25)),
#'                 nrow = n,
#'                 dimnames = list(NULL, paste0("SNP", seq_len(m))))
#' pos <- data.frame(marker     = paste0("SNP", seq_len(m)),
#'                   chromosome = rep(1:5, each = 10L),
#'                   position   = rep(seq(1, 100, 10), 5L))
#' g <- run_basic_gwas(pheno, gmat, pos)
#' nrow(g$results)
#'
#' @export
run_basic_gwas <- function(phenotype, marker_genotypes,
                           marker_positions = NULL,
                           trait_type = "quantitative",
                           min_maf    = 0.05) {
  if (!is.numeric(phenotype))
    rlang::abort("`phenotype` must be a numeric vector.")
  if (!is.matrix(marker_genotypes))
    marker_genotypes <- as.matrix(marker_genotypes)
  if (length(phenotype) != nrow(marker_genotypes))
    rlang::abort("`phenotype` length must equal `nrow(marker_genotypes)`.")

  trait_type <- match.arg(trait_type, c("quantitative", "binary"))

  n_ind     <- length(phenotype)
  n_mrk     <- ncol(marker_genotypes)
  mrk_names <- if (!is.null(colnames(marker_genotypes)))
                 colnames(marker_genotypes)
               else paste0("M", seq_len(n_mrk))

  if (is.null(marker_positions)) {
    marker_positions <- data.frame(
      marker     = mrk_names,
      chromosome = rep(1L, n_mrk),
      position   = seq_len(n_mrk) * 1000L,
      stringsAsFactors = FALSE
    )
  }

  results <- vector("list", n_mrk)
  for (j in seq_len(n_mrk)) {
    snp         <- marker_genotypes[, j]
    allele_freq <- mean(snp, na.rm = TRUE) / 2
    maf         <- min(allele_freq, 1 - allele_freq)

    if (is.na(maf) || maf < min_maf) {
      results[[j]] <- data.frame(
        marker  = mrk_names[j],
        p_value = NA_real_,
        beta    = NA_real_,
        se      = NA_real_,
        maf     = maf,
        stringsAsFactors = FALSE
      )
      next
    }

    complete <- !is.na(snp) & !is.na(phenotype)
    y_fit    <- phenotype[complete]
    x_fit    <- snp[complete]

    fit <- tryCatch({
      if (trait_type == "quantitative")
        stats::lm(y_fit ~ x_fit)
      else
        stats::glm(y_fit ~ x_fit,
                   family = stats::binomial(link = "logit"))
    }, error = function(e) NULL)

    if (is.null(fit)) {
      results[[j]] <- data.frame(
        marker  = mrk_names[j], p_value = NA_real_,
        beta    = NA_real_,    se      = NA_real_, maf = maf,
        stringsAsFactors = FALSE
      )
      next
    }

    coef_tbl <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(coef_tbl) || nrow(coef_tbl) < 2L) {
      results[[j]] <- data.frame(
        marker  = mrk_names[j], p_value = NA_real_,
        beta    = NA_real_,    se      = NA_real_, maf = maf,
        stringsAsFactors = FALSE
      )
      next
    }

    p_col        <- ncol(coef_tbl)
    results[[j]] <- data.frame(
      marker  = mrk_names[j],
      p_value = coef_tbl[2L, p_col],
      beta    = coef_tbl[2L, 1L],
      se      = coef_tbl[2L, 2L],
      maf     = maf,
      stringsAsFactors = FALSE
    )
  }

  res_df <- do.call(rbind, results)
  res_df <- merge(res_df, marker_positions, by = "marker",
                  all.x = TRUE, sort = FALSE)
  res_df <- res_df[order(res_df$chromosome, res_df$position), ]
  res_df$neg_log10_p <- -log10(res_df$p_value)

  col_order <- c("marker", "chromosome", "position",
                 "p_value", "neg_log10_p", "beta", "se", "maf")
  col_order <- col_order[col_order %in% names(res_df)]
  res_df    <- res_df[, col_order, drop = FALSE]
  rownames(res_df) <- NULL

  n_tested <- sum(!is.na(res_df$p_value))
  bonf     <- 0.05 / max(n_tested, 1L)

  structure(
    list(
      results              = res_df,
      trait_type           = trait_type,
      n_individuals        = n_ind,
      n_markers_tested     = n_tested,
      bonferroni_threshold = bonf
    ),
    class = c("gwas_result", "list")
  )
}

#' @export
print.gwas_result <- function(x, ...) {
  cat("GWAS Result\n")
  cat(sprintf("  Trait type        : %s\n", x$trait_type))
  cat(sprintf("  Individuals       : %d\n", x$n_individuals))
  cat(sprintf("  Markers tested    : %d\n", x$n_markers_tested))
  cat(sprintf("  Bonferroni thresh : %.2e\n", x$bonferroni_threshold))
  top <- utils::head(x$results[order(x$results$p_value, na.last = TRUE), ], 5L)
  cat("\nTop 5 associations:\n")
  print(top, ...)
  invisible(x)
}

#' @export
summary.gwas_result <- function(object, ...) {
  print(object, ...)
  n_sig <- sum(object$results$p_value < object$bonferroni_threshold,
               na.rm = TRUE)
  cat(sprintf("\nSignificant hits (Bonferroni): %d\n", n_sig))
  invisible(object)
}

#' Manhattan Plot for GWAS Results
#'
#' Plots SNP association p-values across the genome, with chromosomes
#' shown in alternating colours and a Bonferroni significance threshold
#' line.
#'
#' @param x A `gwas_result` object returned by [run_basic_gwas()].
#' @param threshold Numeric.  -log10(p) significance threshold line
#'   (default `NULL` uses the Bonferroni threshold stored in `x`).
#' @param suggestive Numeric or `NULL`.  Suggestive threshold (default
#'   `NULL` = none).
#' @param top_n Integer.  Number of top SNPs to annotate with their name
#'   (default `5L`).
#' @param title Character.  Plot title (default `"Manhattan Plot"`).
#'
#' @return A `ggplot` object.
#'
#' @seealso [run_basic_gwas()], [plot_qq_gwas()]
#'
#' @examples
#' set.seed(6)
#' n <- 100L; m <- 50L
#' pheno <- rnorm(n)
#' gmat  <- matrix(sample(0:2, n * m, replace = TRUE),
#'                 nrow = n,
#'                 dimnames = list(NULL, paste0("SNP", seq_len(m))))
#' pos <- data.frame(marker     = paste0("SNP", seq_len(m)),
#'                   chromosome = rep(1:5, each = 10L),
#'                   position   = rep(seq(1, 100, 10), 5L))
#' g <- run_basic_gwas(pheno, gmat, pos)
#' p <- plot_manhattan(g)
#' # print(p)
#'
#' @export
plot_manhattan <- function(x, threshold = NULL, suggestive = NULL,
                           top_n = 5L, title = "Manhattan Plot") {
  if (!inherits(x, "gwas_result"))
    rlang::abort("`x` must be a `gwas_result` object.")

  res <- x$results[!is.na(x$results$p_value), ]
  if (nrow(res) == 0L) {
    rlang::warn("No non-NA p-values found; returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = title))
  }

  res$chromosome <- as.factor(res$chromosome)
  chr_levels     <- levels(res$chromosome)

  chrom_max    <- tapply(res$position, res$chromosome, max)
  chrom_offset <- c(0, cumsum(chrom_max[-length(chrom_max)]))
  names(chrom_offset) <- chr_levels
  res$cum_pos <- res$position + chrom_offset[as.character(res$chromosome)]

  chr_mid     <- tapply(res$cum_pos, res$chromosome, mean)
  n_chr       <- length(chr_levels)
  chr_colours <- rep(c("#2166ac", "#74add1"), length.out = n_chr)
  names(chr_colours) <- chr_levels

  thr_logp <- if (!is.null(threshold)) threshold else
    -log10(x$bonferroni_threshold)

  p <- ggplot2::ggplot(res,
                       ggplot2::aes(x = .data$cum_pos,
                                    y = .data$neg_log10_p,
                                    colour = .data$chromosome)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = thr_logp,
                        colour = "red", linetype = "dashed") +
    ggplot2::scale_colour_manual(values = chr_colours, guide = "none") +
    ggplot2::scale_x_continuous(breaks = chr_mid, labels = chr_levels) +
    ggplot2::labs(title = title, x = "Chromosome",
                  y = expression(-log[10](p))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(suggestive))
    p <- p + ggplot2::geom_hline(yintercept = suggestive,
                                  colour = "blue", linetype = "dotted")

  if (top_n > 0L && nrow(res) > 0L) {
    top_snps <- utils::head(res[order(-res$neg_log10_p), ], top_n)
    p <- p + ggplot2::geom_text(
      data = top_snps,
      ggplot2::aes(x = .data$cum_pos, y = .data$neg_log10_p,
                   label = .data$marker),
      colour = "black", size = 2.5, vjust = -0.4,
      show.legend = FALSE
    )
  }
  p
}

#' Q-Q Plot for GWAS P-values
#'
#' Compares the observed distribution of GWAS p-values against the
#' expected uniform null distribution on a -log10 scale.  Genomic
#' inflation factor \eqn{\lambda} is annotated.
#'
#' @param x A `gwas_result` object **or** a numeric vector of p-values.
#' @param title Character.  Plot title (default `"Q-Q Plot"`).
#' @param conf_band Logical.  Draw 95% confidence band (default `TRUE`).
#'
#' @return A `ggplot` object.
#'
#' @seealso [run_basic_gwas()], [plot_manhattan()]
#'
#' @examples
#' pvals <- runif(500)
#' p <- plot_qq_gwas(pvals)
#' # print(p)
#'
#' @export
plot_qq_gwas <- function(x, title = "Q-Q Plot", conf_band = TRUE) {
  if (inherits(x, "gwas_result")) {
    pvals <- x$results$p_value[!is.na(x$results$p_value)]
  } else if (is.numeric(x)) {
    pvals <- x[!is.na(x)]
  } else {
    rlang::abort(
      "`x` must be a `gwas_result` object or a numeric p-value vector."
    )
  }

  if (length(pvals) == 0L) {
    rlang::warn("No valid p-values found; returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = title))
  }

  n   <- length(pvals)
  obs <- sort(-log10(pvals))
  exp <- -log10(stats::ppoints(n))

  chi2_obs <- stats::qchisq(pvals, df = 1L, lower.tail = FALSE)
  lambda   <- stats::median(chi2_obs, na.rm = TRUE) /
              stats::qchisq(0.5, df = 1L)

  qq_df <- data.frame(expected = exp, observed = obs)

  p <- ggplot2::ggplot(qq_df,
                       ggplot2::aes(x = .data$expected,
                                    y = .data$observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         colour = "red", linetype = "dashed") +
    ggplot2::geom_point(size = 0.8, colour = "#2166ac", alpha = 0.7) +
    ggplot2::labs(
      title    = title,
      subtitle = sprintf("Genomic inflation lambda = %.3f", lambda),
      x        = expression("Expected" ~ -log[10](p)),
      y        = expression("Observed" ~ -log[10](p))
    ) +
    ggplot2::theme_bw()

  if (conf_band && n >= 2L) {
    upper <- -log10(stats::qbeta(0.025, seq_len(n), n - seq_len(n) + 1))
    lower <- -log10(stats::qbeta(0.975, seq_len(n), n - seq_len(n) + 1))
    band_df <- data.frame(expected = exp,
                          lower    = sort(lower),
                          upper    = sort(upper))
    p <- p + ggplot2::geom_ribbon(
      data = band_df,
      ggplot2::aes(x = .data$expected,
                   ymin = .data$lower,
                   ymax = .data$upper),
      alpha = 0.2, fill = "grey50", inherit.aes = FALSE
    )
  }
  p
}

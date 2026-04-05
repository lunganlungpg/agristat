# Module 2 — Field Layout Engine
#
# Implements experimental design generation, ANOVA analysis dispatch,
# error-term computation, and field layout visualisation.
#
# Steps 15–24 of the agristat Phase 3 implementation.

# ============================================================================
# Step 15: Experimental Design Framework — internal constructor
# ============================================================================

#' Internal Constructor for field_design Objects
#'
#' Creates the canonical S3 `"field_design"` list structure shared by all
#' design generator functions.  End users should call the public functions
#' such as [design_crd()], [design_rcbd()], etc.
#'
#' @param design_type Character.  One of `"CRD"`, `"RCBD"`, `"factorial"`,
#'   `"split_plot"`, `"split_split_plot"`, `"strip_plot"`,
#'   `"latin_square"`.
#' @param treatments Character or numeric vector of treatment labels.
#' @param blocks Integer.  Number of blocks (NA for CRD).
#' @param replications Integer.  Number of replications.
#' @param randomization_seed Integer or NULL.  Seed used for randomisation.
#' @param layout A data frame containing the plot-level layout.
#' @param metadata Named list of additional design metadata.
#' @param error_structure Named list describing the error strata.
#'
#' @return A list of class `c("<design_type>_design", "field_design")`.
#'
#' @keywords internal
new_field_design <- function(design_type,
                             treatments,
                             blocks       = NA_integer_,
                             replications = NA_integer_,
                             randomization_seed = NULL,
                             layout       = NULL,
                             metadata     = list(),
                             error_structure = list()) {
  structure(
    list(
      design_type        = design_type,
      treatments         = treatments,
      blocks             = blocks,
      replications       = replications,
      randomization_seed = randomization_seed,
      layout             = layout,
      metadata           = metadata,
      error_structure    = error_structure
    ),
    class = c(paste0(tolower(design_type), "_design"), "field_design")
  )
}

# ============================================================================
# Step 16: CRD and RCBD Generators
# ============================================================================

#' Generate a Completely Randomised Design (CRD)
#'
#' Creates a randomised plot layout for a **Completely Randomised Design**
#' (CRD) in which treatments are assigned to experimental units with no
#' blocking.  This is the simplest experimental design and assumes homogeneous
#' experimental material.
#'
#' The statistical model for CRD is:
#' \deqn{y_{ij} = \mu + \tau_i + \varepsilon_{ij}}
#' where \eqn{\mu} is the overall mean, \eqn{\tau_i} is the \eqn{i}-th
#' treatment effect, and \eqn{\varepsilon_{ij} \sim N(0, \sigma^2)}.
#'
#' @param treatments Character or numeric vector of treatment labels.  Must
#'   have at least 2 unique values.
#' @param replicates Positive integer.  Number of replicates (experimental
#'   units) per treatment.
#' @param seed Integer or `NULL`.  Random seed for reproducibility.  Default
#'   `NULL` uses the current RNG state.
#'
#' @return An S3 object of class `c("crd_design", "field_design")` with a
#'   `layout` data frame containing columns:
#'   \describe{
#'     \item{`plot_id`}{Integer plot number (sequential).}
#'     \item{`treatment`}{Treatment label.}
#'     \item{`block`}{Always `NA` for CRD.}
#'     \item{`row`}{Integer row position in the field.}
#'     \item{`column`}{Integer column position in the field.}
#'   }
#'
#' @examples
#' # Three treatments, four replicates each
#' d <- design_crd(
#'   treatments = c("Control", "Low-N", "High-N"),
#'   replicates = 4L,
#'   seed       = 42L
#' )
#' print(d)
#' head(d$layout)
#'
#' @seealso [design_rcbd()], [design_factorial()], [analyze_design()],
#'   [plot_field_layout()]
#' @export
design_crd <- function(treatments, replicates, seed = NULL) {
  # --- input validation ---
  if (missing(treatments) || length(treatments) < 2L) {
    rlang::abort("`treatments` must have at least 2 elements.")
  }
  treatments <- as.character(treatments)
  unique_trts <- unique(treatments)
  if (length(unique_trts) < 2L) {
    rlang::abort("`treatments` must contain at least 2 distinct values.")
  }
  if (!.is_whole_number(replicates) || length(replicates) != 1L ||
      replicates < 1L) {
    rlang::abort("`replicates` must be a positive integer.")
  }
  replicates <- as.integer(replicates)

  # --- randomisation ---
  if (!is.null(seed)) set.seed(seed)
  n_plots <- length(unique_trts) * replicates
  trt_vec  <- rep(unique_trts, times = replicates)
  trt_vec  <- sample(trt_vec)

  # --- grid positions ---
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)
  row_vec <- ceiling(seq_len(n_plots) / n_cols)
  col_vec <- ((seq_len(n_plots) - 1L) %% n_cols) + 1L

  layout <- data.frame(
    plot_id   = seq_len(n_plots),
    treatment = trt_vec,
    block     = NA_character_,
    row       = row_vec,
    column    = col_vec,
    stringsAsFactors = FALSE
  )

  new_field_design(
    design_type        = "CRD",
    treatments         = unique_trts,
    blocks             = NA_integer_,
    replications       = replicates,
    randomization_seed = seed,
    layout             = layout,
    metadata           = list(
      n_treatments = length(unique_trts),
      n_plots      = n_plots,
      n_rows       = n_rows,
      n_cols       = n_cols
    ),
    error_structure = list(
      residual = list(df = n_plots - length(unique_trts))
    )
  )
}

#' Generate a Randomised Complete Block Design (RCBD)
#'
#' Creates a randomised plot layout for a **Randomised Complete Block Design**
#' (RCBD) in which each block contains exactly one complete set of treatments,
#' randomised independently within each block.  Blocks account for a known
#' source of field variation (e.g., soil gradient, topography).
#'
#' The statistical model for RCBD is:
#' \deqn{y_{ij} = \mu + \tau_i + \beta_j + \varepsilon_{ij}}
#' where \eqn{\beta_j} is the \eqn{j}-th block effect.
#'
#' @param treatments Character or numeric vector of treatment labels.
#' @param blocks Positive integer.  Number of complete blocks.
#' @param seed Integer or `NULL`.  Random seed.
#'
#' @return An S3 object of class `c("rcbd_design", "field_design")` with a
#'   `layout` data frame with columns: `plot_id`, `treatment`, `block`,
#'   `row`, `column`.
#'
#' @examples
#' d <- design_rcbd(
#'   treatments = c("V1", "V2", "V3", "V4"),
#'   blocks     = 3L,
#'   seed       = 7L
#' )
#' print(d)
#' head(d$layout, 8)
#'
#' @seealso [design_crd()], [analyze_design()], [plot_field_layout()]
#' @export
design_rcbd <- function(treatments, blocks, seed = NULL) {
  if (missing(treatments) || length(treatments) < 2L) {
    rlang::abort("`treatments` must have at least 2 elements.")
  }
  treatments <- as.character(treatments)
  unique_trts <- unique(treatments)
  if (!.is_whole_number(blocks) || length(blocks) != 1L || blocks < 1L) {
    rlang::abort("`blocks` must be a positive integer.")
  }
  blocks <- as.integer(blocks)

  n_trts  <- length(unique_trts)
  n_plots <- n_trts * blocks

  if (!is.null(seed)) set.seed(seed)

  # Randomise treatments independently within each block
  layout_list <- vector("list", blocks)
  for (b in seq_len(blocks)) {
    trt_rand <- sample(unique_trts)
    layout_list[[b]] <- data.frame(
      treatment = trt_rand,
      block     = as.character(b),
      row       = b,
      column    = seq_len(n_trts),
      stringsAsFactors = FALSE
    )
  }
  layout_df <- do.call(rbind, layout_list)
  layout_df$plot_id <- seq_len(n_plots)
  layout_df <- layout_df[, c("plot_id", "treatment", "block", "row", "column")]

  new_field_design(
    design_type        = "RCBD",
    treatments         = unique_trts,
    blocks             = blocks,
    replications       = blocks,
    randomization_seed = seed,
    layout             = layout_df,
    metadata           = list(
      n_treatments = n_trts,
      n_blocks     = blocks,
      n_plots      = n_plots
    ),
    error_structure = list(
      block    = list(df = blocks - 1L),
      residual = list(df = (blocks - 1L) * (n_trts - 1L))
    )
  )
}

# ============================================================================
# Step 17: Factorial Design Handler
# ============================================================================

#' Generate a Factorial Experimental Design
#'
#' Creates a complete factorial design crossing all combinations of two or
#' more factors.  Supports 2-way, 3-way, and higher-order factorials.  Each
#' unique treatment combination (cell) is replicated `replicates` times.
#'
#' For a two-factor factorial with factors A (a levels) and B (b levels):
#' \deqn{y_{ijk} = \mu + \alpha_i + \beta_j + (\alpha\beta)_{ij} +
#'   \varepsilon_{ijk}}
#'
#' @param factor_names Character vector of factor names (length \eqn{\geq 2}).
#' @param factor_levels Named list where each element is a character (or
#'   numeric) vector of levels for the corresponding factor.  The list names
#'   must match `factor_names`.
#' @param replicates Positive integer.  Number of replicates per cell.
#' @param seed Integer or `NULL`.  Random seed.
#'
#' @return An S3 object of class `c("factorial_design", "field_design")` whose
#'   `layout` data frame contains columns for each factor plus `plot_id`,
#'   `treatment` (concatenated factor-level string), `block` (`NA`), `row`,
#'   and `column`.  The `metadata` element includes `contrast_coding` with
#'   all main-effect and interaction term names.
#'
#' @examples
#' d <- design_factorial(
#'   factor_names  = c("Nitrogen", "Variety"),
#'   factor_levels = list(
#'     Nitrogen = c("Low", "High"),
#'     Variety  = c("V1", "V2", "V3")
#'   ),
#'   replicates = 3L,
#'   seed       = 1L
#' )
#' print(d)
#' nrow(d$layout)   # 2 * 3 * 3 = 18 plots
#'
#' @seealso [design_crd()], [design_rcbd()], [analyze_design()]
#' @export
design_factorial <- function(factor_names, factor_levels, replicates,
                             seed = NULL) {
  if (missing(factor_names) || length(factor_names) < 2L) {
    rlang::abort("`factor_names` must have at least 2 elements.")
  }
  if (!is.list(factor_levels) ||
      !all(factor_names %in% names(factor_levels))) {
    rlang::abort(
      "`factor_levels` must be a named list with names matching `factor_names`."
    )
  }
  if (!.is_whole_number(replicates) || length(replicates) != 1L ||
      replicates < 1L) {
    rlang::abort("`replicates` must be a positive integer.")
  }
  replicates <- as.integer(replicates)

  # Build the complete factorial grid
  fl_ordered <- factor_levels[factor_names]
  grid <- do.call(expand.grid,
                  c(lapply(fl_ordered, as.character),
                    list(stringsAsFactors = FALSE)))
  # expand.grid puts first factor as fastest-cycling; rename to match factor_names
  names(grid) <- factor_names

  n_cells  <- nrow(grid)
  n_plots  <- n_cells * replicates

  # Replicate and randomise
  full_grid <- grid[rep(seq_len(n_cells), times = replicates), ,
                    drop = FALSE]
  if (!is.null(seed)) set.seed(seed)
  full_grid <- full_grid[sample(nrow(full_grid)), , drop = FALSE]
  rownames(full_grid) <- NULL

  # Build treatment label as concatenation of all factor levels
  full_grid$treatment <- apply(full_grid, 1L, paste, collapse = ":")
  full_grid$block     <- NA_character_
  full_grid$plot_id   <- seq_len(n_plots)

  # Grid positions
  n_cols <- ceiling(sqrt(n_plots))
  full_grid$row    <- ceiling(seq_len(n_plots) / n_cols)
  full_grid$column <- ((seq_len(n_plots) - 1L) %% n_cols) + 1L

  # Re-order columns: plot_id, factor cols, treatment, block, row, column
  fixed_cols <- c("plot_id", factor_names, "treatment", "block", "row",
                  "column")
  full_grid  <- full_grid[, fixed_cols, drop = FALSE]

  # Generate contrast coding terms (all main effects + interactions)
  contrast_terms <- unlist(lapply(seq_along(factor_names), function(k) {
    combs <- utils::combn(factor_names, k, simplify = FALSE)
    vapply(combs, paste, FUN.VALUE = character(1L), collapse = " * ")
  }))

  new_field_design(
    design_type        = "factorial",
    treatments         = unique(full_grid$treatment),
    blocks             = NA_integer_,
    replications       = replicates,
    randomization_seed = seed,
    layout             = full_grid,
    metadata           = list(
      factor_names   = factor_names,
      factor_levels  = fl_ordered,
      n_cells        = n_cells,
      n_plots        = n_plots,
      contrast_coding = contrast_terms
    ),
    error_structure = list(
      residual = list(df = n_plots - n_cells)
    )
  )
}

# ============================================================================
# Step 18: Split-Plot Designs
# ============================================================================

#' Generate a Split-Plot Design
#'
#' Creates a split-plot experimental design where one factor (the *whole-plot*
#' factor) is applied to large experimental units, and a second factor (the
#' *sub-plot* factor) is applied within each whole-plot.
#'
#' The split-plot model contains two error strata:
#' \deqn{y_{ijk} = \mu + \rho_k + \alpha_i + \delta_{ik} + \beta_j +
#'   (\alpha\beta)_{ij} + \varepsilon_{ijk}}
#' where \eqn{\delta_{ik}} is the whole-plot error and
#' \eqn{\varepsilon_{ijk}} is the sub-plot error.
#'
#' @param whole_plot_factor Character vector of whole-plot factor levels.
#' @param sub_plot_factor Character vector of sub-plot factor levels.
#' @param blocks Positive integer.  Number of blocks (replications of the
#'   whole-plot factor).
#' @param seed Integer or `NULL`.  Random seed.
#'
#' @return An S3 object of class `c("split_plot_design", "field_design")` with
#'   a `layout` data frame and an `error_structure` list containing
#'   `whole_plot_error` and `sub_plot_error` degrees of freedom.
#'
#' @examples
#' d <- design_split_plot(
#'   whole_plot_factor = c("Irrigation-On", "Irrigation-Off"),
#'   sub_plot_factor   = c("V1", "V2", "V3"),
#'   blocks            = 4L,
#'   seed              = 99L
#' )
#' print(d)
#'
#' @seealso [design_split_split_plot()], [design_strip_plot()],
#'   [analyze_design()]
#' @export
design_split_plot <- function(whole_plot_factor, sub_plot_factor,
                              blocks, seed = NULL) {
  if (missing(whole_plot_factor) || length(whole_plot_factor) < 2L) {
    rlang::abort("`whole_plot_factor` must have at least 2 levels.")
  }
  if (missing(sub_plot_factor) || length(sub_plot_factor) < 2L) {
    rlang::abort("`sub_plot_factor` must have at least 2 levels.")
  }
  if (!.is_whole_number(blocks) || length(blocks) != 1L || blocks < 2L) {
    rlang::abort("`blocks` must be an integer >= 2.")
  }
  whole_plot_factor <- as.character(whole_plot_factor)
  sub_plot_factor   <- as.character(sub_plot_factor)
  blocks            <- as.integer(blocks)

  if (!is.null(seed)) set.seed(seed)

  a <- length(whole_plot_factor)   # whole-plot levels
  b <- length(sub_plot_factor)     # sub-plot levels
  n_whole_plots <- a * blocks
  n_plots       <- n_whole_plots * b

  all_rows <- list()
  wp_id    <- 0L
  for (blk in seq_len(blocks)) {
    wp_order <- sample(whole_plot_factor)
    for (wp_level in wp_order) {
      wp_id    <- wp_id + 1L
      sp_order <- sample(sub_plot_factor)
      all_rows[[wp_id]] <- data.frame(
        whole_plot_id    = wp_id,
        whole_plot_level = wp_level,
        sub_plot_level   = sp_order,
        block            = as.character(blk),
        stringsAsFactors = FALSE
      )
    }
  }

  layout_df <- do.call(rbind, all_rows)
  layout_df$plot_id   <- seq_len(nrow(layout_df))
  layout_df$treatment <- paste(layout_df$whole_plot_level,
                               layout_df$sub_plot_level, sep = ":")
  layout_df$row       <- ceiling(layout_df$plot_id / b)
  layout_df$column    <- ((layout_df$plot_id - 1L) %% b) + 1L
  layout_df <- layout_df[, c("plot_id", "whole_plot_id", "whole_plot_level",
                              "sub_plot_level", "treatment", "block",
                              "row", "column")]

  new_field_design(
    design_type        = "split_plot",
    treatments         = unique(layout_df$treatment),
    blocks             = blocks,
    replications       = blocks,
    randomization_seed = seed,
    layout             = layout_df,
    metadata           = list(
      whole_plot_factor = whole_plot_factor,
      sub_plot_factor   = sub_plot_factor,
      n_whole_plots     = n_whole_plots,
      n_plots           = n_plots
    ),
    error_structure = list(
      whole_plot_error = list(
        df = blocks * (a - 1L),
        description = "Error(block:whole_plot_factor)"
      ),
      sub_plot_error = list(
        df = n_whole_plots * (b - 1L),
        description = "Error(whole_plot_id:sub_plot_factor)"
      )
    )
  )
}

#' Generate a Split-Split-Plot Design
#'
#' Extends the split-plot design to three factor levels: whole-plot, sub-plot,
#' and sub-sub-plot.  This design has three error strata.
#'
#' @param whole_plot_factor Character vector of whole-plot factor levels.
#' @param sub_plot_factor Character vector of sub-plot factor levels.
#' @param sub_sub_plot_factor Character vector of sub-sub-plot factor levels.
#' @param blocks Positive integer >= 2.  Number of blocks.
#' @param seed Integer or `NULL`.  Random seed.
#'
#' @return An S3 object of class `c("split_split_plot_design", "field_design")`
#'   with three-stratum `error_structure`.
#'
#' @examples
#' d <- design_split_split_plot(
#'   whole_plot_factor     = c("Tillage-Deep", "Tillage-Shallow"),
#'   sub_plot_factor       = c("Irrigation-On", "Irrigation-Off"),
#'   sub_sub_plot_factor   = c("V1", "V2"),
#'   blocks                = 3L,
#'   seed                  = 5L
#' )
#' print(d)
#'
#' @seealso [design_split_plot()], [analyze_design()]
#' @export
design_split_split_plot <- function(whole_plot_factor, sub_plot_factor,
                                    sub_sub_plot_factor, blocks,
                                    seed = NULL) {
  if (length(whole_plot_factor) < 2L)
    rlang::abort("`whole_plot_factor` must have at least 2 levels.")
  if (length(sub_plot_factor) < 2L)
    rlang::abort("`sub_plot_factor` must have at least 2 levels.")
  if (length(sub_sub_plot_factor) < 2L)
    rlang::abort("`sub_sub_plot_factor` must have at least 2 levels.")
  if (!.is_whole_number(blocks) || blocks < 2L)
    rlang::abort("`blocks` must be an integer >= 2.")

  whole_plot_factor   <- as.character(whole_plot_factor)
  sub_plot_factor     <- as.character(sub_plot_factor)
  sub_sub_plot_factor <- as.character(sub_sub_plot_factor)
  blocks              <- as.integer(blocks)

  if (!is.null(seed)) set.seed(seed)

  a <- length(whole_plot_factor)
  b <- length(sub_plot_factor)
  cc <- length(sub_sub_plot_factor)

  rows <- list()
  wp_id  <- 0L
  ssp_id <- 0L

  for (blk in seq_len(blocks)) {
    wp_order <- sample(whole_plot_factor)
    for (wp_level in wp_order) {
      wp_id <- wp_id + 1L
      sp_order <- sample(sub_plot_factor)
      for (sp_level in sp_order) {
        ssp_order <- sample(sub_sub_plot_factor)
        for (ssp_level in ssp_order) {
          ssp_id <- ssp_id + 1L
          rows[[ssp_id]] <- data.frame(
            whole_plot_id        = wp_id,
            whole_plot_level     = wp_level,
            sub_plot_level       = sp_level,
            sub_sub_plot_level   = ssp_level,
            block                = as.character(blk),
            stringsAsFactors     = FALSE
          )
        }
      }
    }
  }

  layout_df <- do.call(rbind, rows)
  layout_df$plot_id   <- seq_len(nrow(layout_df))
  layout_df$treatment <- paste(layout_df$whole_plot_level,
                               layout_df$sub_plot_level,
                               layout_df$sub_sub_plot_level, sep = ":")
  layout_df$row    <- layout_df$whole_plot_id
  layout_df$column <- seq_len(nrow(layout_df))
  layout_df <- layout_df[, c("plot_id", "whole_plot_id", "whole_plot_level",
                              "sub_plot_level", "sub_sub_plot_level",
                              "treatment", "block", "row", "column")]

  n_plots <- nrow(layout_df)

  new_field_design(
    design_type        = "split_split_plot",
    treatments         = unique(layout_df$treatment),
    blocks             = blocks,
    replications       = blocks,
    randomization_seed = seed,
    layout             = layout_df,
    metadata           = list(
      whole_plot_factor   = whole_plot_factor,
      sub_plot_factor     = sub_plot_factor,
      sub_sub_plot_factor = sub_sub_plot_factor,
      n_plots             = n_plots
    ),
    error_structure = list(
      whole_plot_error = list(
        df = blocks * (a - 1L),
        description = "Error A (whole-plot)"
      ),
      sub_plot_error = list(
        df = blocks * a * (b - 1L),
        description = "Error B (sub-plot)"
      ),
      sub_sub_plot_error = list(
        df = blocks * a * b * (cc - 1L),
        description = "Error C (sub-sub-plot)"
      )
    )
  )
}

#' Generate a Strip-Plot Design
#'
#' Creates a strip-plot (or split-block) design in which two treatment factors
#' are applied in orthogonal strips across a field, forming a grid of
#' intersecting plots.  Both main effects and their interaction can be
#' estimated.
#'
#' @param row_factor Character vector of factor levels applied to horizontal
#'   strips (rows).
#' @param col_factor Character vector of factor levels applied to vertical
#'   strips (columns).
#' @param blocks Positive integer.  Number of complete blocks.
#' @param seed Integer or `NULL`.  Random seed.
#'
#' @return An S3 object of class `c("strip_plot_design", "field_design")`.
#'
#' @examples
#' d <- design_strip_plot(
#'   row_factor = c("Herbicide-A", "Herbicide-B", "Herbicide-C"),
#'   col_factor = c("Low-N", "Mid-N", "High-N"),
#'   blocks     = 3L,
#'   seed       = 22L
#' )
#' print(d)
#'
#' @seealso [design_split_plot()], [analyze_design()]
#' @export
design_strip_plot <- function(row_factor, col_factor, blocks, seed = NULL) {
  if (length(row_factor) < 2L)
    rlang::abort("`row_factor` must have at least 2 levels.")
  if (length(col_factor) < 2L)
    rlang::abort("`col_factor` must have at least 2 levels.")
  if (!.is_whole_number(blocks) || blocks < 2L)
    rlang::abort("`blocks` must be an integer >= 2.")

  row_factor <- as.character(row_factor)
  col_factor <- as.character(col_factor)
  blocks     <- as.integer(blocks)

  if (!is.null(seed)) set.seed(seed)

  a <- length(row_factor)
  b <- length(col_factor)

  rows <- list()
  plot_id <- 0L
  for (blk in seq_len(blocks)) {
    rf_order <- sample(row_factor)
    cf_order <- sample(col_factor)
    for (ri in seq_along(rf_order)) {
      for (ci in seq_along(cf_order)) {
        plot_id <- plot_id + 1L
        rows[[plot_id]] <- data.frame(
          row_level  = rf_order[ri],
          col_level  = cf_order[ci],
          block      = as.character(blk),
          row        = (blk - 1L) * a + ri,
          column     = ci,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  layout_df           <- do.call(rbind, rows)
  layout_df$plot_id   <- seq_len(nrow(layout_df))
  layout_df$treatment <- paste(layout_df$row_level, layout_df$col_level,
                               sep = ":")
  layout_df <- layout_df[, c("plot_id", "row_level", "col_level",
                              "treatment", "block", "row", "column")]

  n_plots <- nrow(layout_df)

  new_field_design(
    design_type        = "strip_plot",
    treatments         = unique(layout_df$treatment),
    blocks             = blocks,
    replications       = blocks,
    randomization_seed = seed,
    layout             = layout_df,
    metadata           = list(
      row_factor = row_factor,
      col_factor = col_factor,
      n_plots    = n_plots
    ),
    error_structure = list(
      row_error = list(
        df = blocks * (a - 1L),
        description = "Error for row factor"
      ),
      col_error = list(
        df = blocks * (b - 1L),
        description = "Error for column factor"
      ),
      interaction_error = list(
        df = (a - 1L) * (b - 1L) * blocks,
        description = "Error for row x col interaction"
      )
    )
  )
}

# ============================================================================
# Step 19: Latin Square Support
# ============================================================================

#' Generate a Latin Square Design
#'
#' Creates an \eqn{n \times n} Latin Square in which each treatment appears
#' exactly once in each row and exactly once in each column, providing
#' two-way blocking against known sources of variation.
#'
#' The Latin Square model is:
#' \deqn{y_{ijk} = \mu + \rho_i + \kappa_j + \tau_k + \varepsilon_{ijk}}
#' where \eqn{\rho_i}, \eqn{\kappa_j} are the row and column blocking effects
#' respectively, and each \eqn{(i, j, k)} combination appears once.
#'
#' @param factor_levels Character or numeric vector of treatment labels.
#'   The number of treatments equals the number of rows and columns
#'   (\eqn{n = }\code{length(factor_levels)}).
#' @param seed Integer or `NULL`.  Random seed for the random Latin Square
#'   generation.
#'
#' @return An S3 object of class `c("latin_square_design", "field_design")`
#'   with a `layout` data frame containing columns: `plot_id`, `treatment`,
#'   `block` (`NA`), `row`, `column`.
#'
#' @examples
#' d <- design_latin_square(
#'   factor_levels = c("A", "B", "C", "D"),
#'   seed          = 2024L
#' )
#' print(d)
#' d$layout
#'
#' @seealso [design_crd()], [design_rcbd()], [analyze_design()]
#' @export
design_latin_square <- function(factor_levels, seed = NULL) {
  if (missing(factor_levels) || length(factor_levels) < 2L) {
    rlang::abort("`factor_levels` must have at least 2 elements.")
  }
  factor_levels <- as.character(factor_levels)
  n <- length(factor_levels)

  if (!is.null(seed)) set.seed(seed)

  # Build a random Latin square using row-by-row permutation
  # Start with a standard cyclic Latin square, then randomly permute rows
  # and columns
  base_square <- matrix(NA_character_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    base_square[i, ] <- factor_levels[((seq_len(n) + i - 2L) %% n) + 1L]
  }

  # Random permutation of rows and columns
  row_perm <- sample(n)
  col_perm <- sample(n)
  ls       <- base_square[row_perm, col_perm]

  # Build data frame
  layout_df <- data.frame(
    plot_id   = seq_len(n * n),
    treatment = as.vector(t(ls)),
    block     = NA_character_,
    row       = rep(seq_len(n), each = n),
    column    = rep(seq_len(n), times = n),
    stringsAsFactors = FALSE
  )

  new_field_design(
    design_type        = "latin_square",
    treatments         = factor_levels,
    blocks             = NA_integer_,
    replications       = 1L,
    randomization_seed = seed,
    layout             = layout_df,
    metadata           = list(
      n        = n,
      n_plots  = n * n,
      square   = ls
    ),
    error_structure = list(
      residual = list(df = (n - 1L) * (n - 2L))
    )
  )
}

# ============================================================================
# Step 20: ANOVA Engine
# ============================================================================

#' Analyse an Experimental Design
#'
#' Dispatches to the appropriate statistical analysis based on the design type
#' stored in a `field_design` object.  Supports CRD, RCBD, factorial, Latin
#' Square, split-plot, split-split-plot, and strip-plot analyses.
#'
#' | Design | Model | Function |
#' |---|---|---|
#' | CRD | `response ~ treatment` | `stats::aov()` |
#' | RCBD | `response ~ treatment + block` | `stats::aov()` |
#' | Factorial | `response ~ A * B * ...` | `stats::aov()` |
#' | Split-plot | mixed model | `lme4::lmer()` |
#' | Latin Square | `response ~ row + column + treatment` | `stats::aov()` |
#'
#' @param design A `field_design` object (created by one of the
#'   `design_*()` functions).
#' @param response Character string.  Name of the response variable column
#'   in `data`.
#' @param data A data frame containing the response variable and all design
#'   columns (must include `plot_id` to join with `design$layout`).
#'
#' @return An S3 object of class `c("agristat_aov", "aov")` or
#'   `c("agristat_aov", "lmerMod")` (for split-plot designs), depending on
#'   the design type.  The object carries a `design` attribute with the
#'   original `field_design`.
#'
#' @examples
#' set.seed(1L)
#' d <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
#' df <- d$layout
#' df$yield <- c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 11))
#' result <- analyze_design(d, response = "yield", data = df)
#' summary(result)
#'
#' @seealso [design_crd()], [design_rcbd()], [design_factorial()],
#'   [plot_anova_diagnostics()]
#' @export
analyze_design <- function(design, response, data) {
  if (!inherits(design, "field_design")) {
    rlang::abort("`design` must be a `field_design` object.")
  }
  if (!is.character(response) || length(response) != 1L ||
      nchar(response) == 0L) {
    rlang::abort("`response` must be a non-empty character string.")
  }
  .validate_data_frame(data, "data")
  if (!response %in% names(data)) {
    rlang::abort(sprintf("Column '%s' not found in `data`.", response))
  }

  dtype <- design$design_type

  result <- switch(dtype,
    CRD = {
      if (!"treatment" %in% names(data)) {
        data <- merge(data, design$layout[, c("plot_id", "treatment")],
                      by = "plot_id", all.x = TRUE)
      }
      data$treatment <- as.factor(data$treatment)
      fm <- stats::as.formula(paste(response, "~ treatment"))
      stats::aov(fm, data = data)
    },
    RCBD = {
      needed <- c("treatment", "block")
      miss   <- needed[!needed %in% names(data)]
      if (length(miss) > 0L) {
        data <- merge(data,
                      design$layout[, c("plot_id", miss)],
                      by = "plot_id", all.x = TRUE)
      }
      data$treatment <- as.factor(data$treatment)
      data$block     <- as.factor(data$block)
      fm <- stats::as.formula(paste(response, "~ treatment + block"))
      stats::aov(fm, data = data)
    },
    factorial = {
      fn <- design$metadata$factor_names
      miss <- fn[!fn %in% names(data)]
      if (length(miss) > 0L) {
        data <- merge(data, design$layout[, c("plot_id", fn)],
                      by = "plot_id", all.x = TRUE)
      }
      for (f in fn) data[[f]] <- as.factor(data[[f]])
      fm_rhs <- paste(fn, collapse = " * ")
      fm     <- stats::as.formula(paste(response, "~", fm_rhs))
      stats::aov(fm, data = data)
    },
    split_plot = {
      needed <- c("whole_plot_level", "sub_plot_level", "block",
                  "whole_plot_id")
      miss <- needed[!needed %in% names(data)]
      if (length(miss) > 0L) {
        data <- merge(data, design$layout[, c("plot_id", miss)],
                      by = "plot_id", all.x = TRUE)
      }
      data$whole_plot_level <- as.factor(data$whole_plot_level)
      data$sub_plot_level   <- as.factor(data$sub_plot_level)
      data$block            <- as.factor(data$block)
      data$whole_plot_id    <- as.factor(data$whole_plot_id)
      fm <- stats::as.formula(
        paste(response,
              "~ whole_plot_level * sub_plot_level + (1 | whole_plot_id)")
      )
      if (!requireNamespace("lme4", quietly = TRUE)) {
        rlang::abort("Package 'lme4' is required for split-plot analysis.")
      }
      lme4::lmer(fm, data = data, REML = TRUE)
    },
    split_split_plot = {
      needed <- c("whole_plot_level", "sub_plot_level",
                  "sub_sub_plot_level", "block", "whole_plot_id")
      miss <- needed[!needed %in% names(data)]
      if (length(miss) > 0L) {
        data <- merge(data, design$layout[, c("plot_id", miss)],
                      by = "plot_id", all.x = TRUE)
      }
      data$whole_plot_level   <- as.factor(data$whole_plot_level)
      data$sub_plot_level     <- as.factor(data$sub_plot_level)
      data$sub_sub_plot_level <- as.factor(data$sub_sub_plot_level)
      data$whole_plot_id      <- as.factor(data$whole_plot_id)
      if (!requireNamespace("lme4", quietly = TRUE)) {
        rlang::abort("Package 'lme4' is required for split-split-plot analysis.")
      }
      fm <- stats::as.formula(paste(
        response,
        "~ whole_plot_level * sub_plot_level * sub_sub_plot_level",
        "+ (1 | whole_plot_id)"
      ))
      lme4::lmer(fm, data = data, REML = TRUE)
    },
    latin_square = {
      needed <- c("treatment", "row", "column")
      miss <- needed[!needed %in% names(data)]
      if (length(miss) > 0L) {
        data <- merge(data, design$layout[, c("plot_id", miss)],
                      by = "plot_id", all.x = TRUE)
      }
      data$treatment <- as.factor(data$treatment)
      data$row       <- as.factor(data$row)
      data$column    <- as.factor(data$column)
      fm <- stats::as.formula(paste(response, "~ row + column + treatment"))
      stats::aov(fm, data = data)
    },
    strip_plot = {
      needed <- c("row_level", "col_level", "block")
      miss <- needed[!needed %in% names(data)]
      if (length(miss) > 0L) {
        data <- merge(data, design$layout[, c("plot_id", miss)],
                      by = "plot_id", all.x = TRUE)
      }
      data$row_level <- as.factor(data$row_level)
      data$col_level <- as.factor(data$col_level)
      data$block     <- as.factor(data$block)
      fm <- stats::as.formula(
        paste(response, "~ block + row_level + col_level + row_level:col_level")
      )
      stats::aov(fm, data = data)
    },
    rlang::abort(sprintf("Unsupported design type: '%s'.", dtype))
  )

  # Wrap the model in an agristat_aov list to unify S3 and S4 results
  structure(
    list(
      model    = result,
      design   = design,
      response = response
    ),
    class = "agristat_aov"
  )
}
# ============================================================================

#' Compute Error Term Degrees of Freedom for a Design
#'
#' A convenience wrapper around [get_error_structure()] that returns a named
#' numeric vector of degrees of freedom for each error stratum in the
#' `field_design` object.
#'
#' @param design A `field_design` object.
#'
#' @return Named numeric vector of degrees of freedom.
#'
#' @examples
#' d <- design_split_plot(
#'   whole_plot_factor = c("I+", "I-"),
#'   sub_plot_factor   = c("V1", "V2", "V3"),
#'   blocks = 4L, seed = 1L
#' )
#' compute_error_terms(d)
#'
#' @seealso [get_error_structure()], [extract_variance_components()]
#' @export
compute_error_terms <- function(design) {
  if (!inherits(design, "field_design")) {
    rlang::abort("`design` must be a `field_design` object.")
  }
  es <- design$error_structure
  if (length(es) == 0L) return(numeric(0L))
  vapply(es, function(x) as.numeric(x$df), FUN.VALUE = numeric(1L))
}

# ============================================================================
# Step 22: Visualization Functions
# ============================================================================

#' Plot a Field Layout
#'
#' Produces a **ggplot2** tile plot showing the spatial arrangement of
#' treatments across rows and columns of the field.  Plots are coloured by
#' treatment and labelled with their `plot_id`.
#'
#' @param design A `field_design` object.
#' @param title Character.  Optional plot title.  Default uses the design type.
#' @param label_col Character.  Column from `design$layout` to use as tile
#'   labels (default `"plot_id"`).
#' @param fill_col Character.  Column to use for fill colour
#'   (default `"treatment"`).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- design_rcbd(c("A", "B", "C", "D"), blocks = 3L, seed = 1L)
#' p <- plot_field_layout(d)
#' print(p)
#'
#' @seealso [design_crd()], [design_rcbd()], [plot_anova_diagnostics()]
#' @export
plot_field_layout <- function(design,
                              title     = NULL,
                              label_col = "plot_id",
                              fill_col  = "treatment") {
  if (!inherits(design, "field_design")) {
    rlang::abort("`design` must be a `field_design` object.")
  }
  layout <- design$layout
  if (!all(c("row", "column") %in% names(layout))) {
    rlang::abort(
      "`design$layout` must contain 'row' and 'column' columns."
    )
  }
  if (!fill_col %in% names(layout)) {
    rlang::abort(sprintf(
      "Column '%s' not found in `design$layout`.", fill_col
    ))
  }
  if (!label_col %in% names(layout)) {
    rlang::abort(sprintf(
      "Column '%s' not found in `design$layout`.", label_col
    ))
  }

  if (is.null(title)) {
    title <- sprintf("Field Layout: %s", design$design_type)
  }

  layout[[fill_col]]  <- as.factor(layout[[fill_col]])
  layout[[label_col]] <- as.character(layout[[label_col]])

  # ggplot2 aesthetics use string-based column references via aes()
  ggplot2::ggplot(
    layout,
    ggplot2::aes(
      x    = .data[["column"]],
      y    = .data[["row"]],
      fill = .data[[fill_col]]
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[[label_col]]),
      size = 3
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_fill_brewer(palette = "Set3", name = fill_col) +
    ggplot2::labs(
      title = title,
      x     = "Column",
      y     = "Row"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(size = 9)
    )
}

#' Plot ANOVA Diagnostic Plots
#'
#' Produces a 4-panel diagnostic plot for a fitted model object.  The panels
#' show:
#' 1. Residuals vs Fitted
#' 2. Normal Q-Q plot of residuals
#' 3. Scale-Location (sqrt(|residuals|) vs Fitted)
#' 4. Residuals vs Leverage (Cook's distance contours)
#'
#' @param model An `aov`, `lm`, or `lmerMod` object, or an `agristat_aov`
#'   object.
#' @param which Integer vector (1–4) selecting which panels to draw.  Default
#'   `1:4` (all four panels).
#'
#' @return Called for its side effect (base graphics plot).  Returns the model
#'   invisibly.
#'
#' @examples
#' set.seed(1L)
#' d   <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
#' df  <- d$layout
#' df$yield <- c(rnorm(5, 10), rnorm(5, 12), rnorm(5, 11))
#' res <- analyze_design(d, response = "yield", data = df)
#' plot_anova_diagnostics(res)
#'
#' @seealso [analyze_design()], [plot_field_layout()]
#' @export
plot_anova_diagnostics <- function(model, which = 1:4) {
  # Unwrap agristat_aov list wrapper
  inner <- if (inherits(model, "agristat_aov")) model$model else model

  if (inherits(inner, "lmerMod")) {
    # For mixed models use base plots of residuals/fitted
    resid_vals  <- stats::residuals(inner)
    fitted_vals <- stats::fitted(inner)
    old_par <- graphics::par(mfrow = c(2L, 2L))
    on.exit(graphics::par(old_par), add = TRUE)

    if (1L %in% which) {
      graphics::plot(fitted_vals, resid_vals,
                     xlab = "Fitted values", ylab = "Residuals",
                     main = "Residuals vs Fitted")
      graphics::abline(h = 0, lty = 2)
    }
    if (2L %in% which) {
      stats::qqnorm(resid_vals, main = "Normal Q-Q")
      stats::qqline(resid_vals, col = "steelblue")
    }
    if (3L %in% which) {
      graphics::plot(fitted_vals, sqrt(abs(resid_vals)),
                     xlab = "Fitted values",
                     ylab = expression(sqrt(abs(Residuals))),
                     main = "Scale-Location")
    }
    if (4L %in% which) {
      graphics::plot(stats::residuals(inner, type = "pearson"),
                     xlab = "Observation index",
                     ylab = "Standardised residuals",
                     main = "Standardised Residuals")
      graphics::abline(h = c(-2, 0, 2), lty = c(2, 1, 2))
    }
    return(invisible(model))
  }
  # For aov/lm objects use the built-in plot method
  graphics::plot(inner, which = which)
  invisible(model)
}

# ============================================================================
# S3 Methods: field_design
# ============================================================================

#' Print a field_design Object
#'
#' @param x A `field_design` object.
#' @param ... Ignored.
#'
#' @return `x` invisibly.
#' @export
print.field_design <- function(x, ...) {
  cat(sprintf(
    "Field Design: %s\n  Treatments  : %d (%s)\n  Blocks      : %s\n  Replications: %s\n  Plots       : %d\n  Seed        : %s\n",
    x$design_type,
    length(x$treatments),
    paste(x$treatments, collapse = ", "),
    if (is.na(x$blocks)) "NA" else as.character(x$blocks),
    if (is.na(x$replications)) "NA" else as.character(x$replications),
    if (!is.null(x$layout)) nrow(x$layout) else 0L,
    if (is.null(x$randomization_seed)) "NULL" else
      as.character(x$randomization_seed)
  ))
  invisible(x)
}

#' Summarise a field_design Object
#'
#' @param object A `field_design` object.
#' @param ... Ignored.
#'
#' @return `object` invisibly.
#' @export
summary.field_design <- function(object, ...) {
  print(object)
  cat("\nError structure:\n")
  es <- object$error_structure
  if (length(es) == 0L) {
    cat("  (none)\n")
  } else {
    for (nm in names(es)) {
      cat(sprintf("  %s: df = %s\n", nm,
                  if (is.null(es[[nm]]$df)) "?" else
                    as.character(es[[nm]]$df)))
    }
  }
  if (!is.null(object$layout)) {
    cat(sprintf("\nLayout (first 6 rows):\n"))
    print(utils::head(object$layout, 6L))
  }
  invisible(object)
}

# ============================================================================
# S3 Methods: agristat_aov
# ============================================================================

#' Print an agristat_aov Object
#'
#' @param x An `agristat_aov` object.
#' @param ... Passed to the underlying model's print method.
#'
#' @return `x` invisibly.
#' @export
print.agristat_aov <- function(x, ...) {
  cat(sprintf("agristat ANOVA -- Design: %s\n", x$design$design_type))
  print(x$model, ...)
  invisible(x)
}

#' Summarise an agristat_aov Object
#'
#' @param object An `agristat_aov` object.
#' @param ... Passed to the underlying model's summary method.
#'
#' @return The summary object returned by the underlying model class.
#' @export
summary.agristat_aov <- function(object, ...) {
  cat(sprintf("agristat ANOVA Summary -- Design: %s\n\n",
              object$design$design_type))
  summary(object$model, ...)
}

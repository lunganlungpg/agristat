# Module 6 — Thesis Assembly & Export Suite
#
# Implements Steps 43-50 of the agristat Phase 7 implementation:
#   Step 43 : generate_methods_text()  — auto-generate methods paragraph
#   Step 44 : export_results_table()   — publication-ready tables
#   Step 45 : export_figure()          — save plots in publication formats
#   Step 46 : generate_appendix()      — appendix with assumption tests
#   Step 47 : compile_report()         — assemble Rmd/md report
#   Step 48 : thesis_template_methods()   — methods section template
#   Step 49 : thesis_template_results()   — results section template
#   Step 50 : thesis_template_discussion() / fill_template()

# ============================================================================
# Internal helpers
# ============================================================================

.p_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return(".")
  return("ns")
}

.detect_design_info <- function(design_obj) {
  cls <- class(design_obj)

  if (inherits(design_obj, "field_design") ||
      inherits(design_obj, "crd_design") ||
      inherits(design_obj, "rcbd_design") ||
      inherits(design_obj, "factorial_design") ||
      inherits(design_obj, "split_plot_design") ||
      inherits(design_obj, "latin_square_design")) {
    dtype  <- if (!is.null(design_obj$design_type)) design_obj$design_type else cls[1L]
    n_trt  <- length(unique(design_obj$treatments))
    n_rep  <- design_obj$replications
    n_blk  <- if (is.na(design_obj$blocks)) NULL else design_obj$blocks
    return(list(type = dtype, n_treatments = n_trt,
                n_replicates = n_rep, n_blocks = n_blk))
  }

  if (inherits(design_obj, "agristat_aov")) {
    d <- design_obj$design
    if (!is.null(d)) return(.detect_design_info(d))
    return(list(type = "ANOVA", n_treatments = NA, n_replicates = NA, n_blocks = NULL))
  }

  if (inherits(design_obj, "gxe_result") || inherits(design_obj, "ammi_result")) {
    return(list(type = "GxE / AMMI", n_treatments = NA, n_replicates = NA, n_blocks = NULL))
  }

  if (inherits(design_obj, "gwas_result")) {
    return(list(type = "GWAS", n_treatments = NA, n_replicates = NA, n_blocks = NULL))
  }

  # Bare named list with $design slot
  if (is.list(design_obj) && !is.null(design_obj$design)) {
    return(.detect_design_info(design_obj$design))
  }

  list(type = paste(cls, collapse = "/"), n_treatments = NA,
       n_replicates = NA, n_blocks = NULL)
}

.design_model_formula <- function(design_type) {
  switch(
    toupper(as.character(design_type)),
    "CRD"          = "y_ij = mu + tau_i + epsilon_ij",
    "RCBD"         = "y_ij = mu + tau_i + beta_j + epsilon_ij",
    "FACTORIAL"    = "y_ijk = mu + alpha_i + beta_j + (alpha*beta)_ij + epsilon_ijk",
    "SPLIT_PLOT"   = "y_ijk = mu + alpha_i + delta_ij + beta_k + (alpha*beta)_ik + epsilon_ijk",
    "LATIN_SQUARE" = "y_ijk = mu + tau_i + rho_j + gamma_k + epsilon_ijk",
    "GXE / AMMI"   = "Y_ger = mu + g_g + e_e + (ge)_ge + epsilon_ger",
    "GWAS"         = "y = Xb + Zu + epsilon",
    "y = mu + treatment + error"
  )
}

.design_human_name <- function(design_type) {
  switch(
    toupper(as.character(design_type)),
    "CRD"                  = "Completely Randomised Design (CRD)",
    "RCBD"                 = "Randomised Complete Block Design (RCBD)",
    "FACTORIAL"            = "Factorial Design",
    "SPLIT_PLOT"           = "Split-Plot Design",
    "SPLIT_SPLIT_PLOT"     = "Split-Split-Plot Design",
    "STRIP_PLOT"           = "Strip-Plot Design",
    "LATIN_SQUARE"         = "Latin Square Design",
    "ALPHA_LATTICE"        = "Alpha-Lattice Design",
    "INCOMPLETE_BLOCK"     = "Incomplete Block Design",
    "GXE / AMMI"           = "Genotype-by-Environment Interaction (GxE/AMMI) Analysis",
    "GWAS"                 = "Genome-Wide Association Study (GWAS)",
    as.character(design_type)
  )
}

.r_version_string <- function() {
  rv <- R.Version()
  paste0("R version ", rv$major, ".", rv$minor,
         " (", rv$nickname, ")")
}

.key_packages_string <- function() {
  pkgs <- c("agristat", "ggplot2", "stats")
  versions <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) "?")
  }, character(1))
  paste(sprintf("%s (v%s)", pkgs, versions), collapse = ", ")
}

# ============================================================================
# Step 43: generate_methods_text()
# ============================================================================

#' Generate a Methods-Section Paragraph
#'
#' Automatically constructs a publication-quality methods paragraph from a
#' design object and optional analysis/post-hoc objects.  The text is suitable
#' for direct inclusion in a thesis or journal article methods section.
#'
#' @param design_obj An experimental design or analysis object.  Accepted
#'   classes include `field_design` (and its subclasses such as `crd_design`,
#'   `rcbd_design`, etc.), `agristat_aov`, `gxe_result`, `ammi_result`, and
#'   `gwas_result`.  A bare named list containing a `$design` slot is also
#'   accepted.
#' @param analysis_obj Optional. An `agristat_aov` object produced by
#'   [analyze_design()].  When provided, the analysis model is included in
#'   the text.
#' @param posthoc_obj Optional. A post-hoc result object (`dmrt_result`,
#'   `tukey_hsd_result`, `dunnett_result`, or `mean_separation_result`).
#'   When provided, the post-hoc procedure is described.
#' @param software Logical.  If `TRUE` (default) a software note citing the R
#'   version and key packages is appended.
#'
#' @return An S3 object of class `"methods_text"` with components:
#'   \describe{
#'     \item{`text`}{Character string. Full methods paragraph ready for pasting
#'       into a document.}
#'     \item{`sections`}{Named list with elements `design_description`,
#'       `analysis_description`, and `software_note`.}
#'     \item{`references`}{Character vector of suggested citations.}
#'   }
#'
#' @examples
#' d <- design_rcbd(
#'   treatments = c("Control", "T1", "T2", "T3"),
#'   blocks     = 4L,
#'   seed       = 1L
#' )
#' mt <- generate_methods_text(d)
#' cat(mt$text)
#'
#' @seealso [thesis_template_methods()], [fill_template()]
#' @export
generate_methods_text <- function(design_obj,
                                  analysis_obj = NULL,
                                  posthoc_obj  = NULL,
                                  software     = TRUE) {
  if (is.null(design_obj)) {
    rlang::abort("`design_obj` must not be NULL.")
  }

  info <- .detect_design_info(design_obj)
  design_name <- .design_human_name(info$type)
  formula_str <- .design_model_formula(info$type)

  # ---- design description --------------------------------------------------
  n_trt_str <- if (!is.na(info$n_treatments)) {
    paste0(info$n_treatments, " treatment", if (info$n_treatments != 1) "s" else "")
  } else {
    "multiple treatments"
  }
  n_rep_str <- if (!is.null(info$n_replicates) && !is.na(info$n_replicates)) {
    paste0(info$n_replicates, " replication", if (info$n_replicates != 1) "s" else "")
  } else {
    "multiple replications"
  }
  blk_str <- if (!is.null(info$n_blocks) && !is.na(info$n_blocks)) {
    paste0(" arranged in ", info$n_blocks, " blocks")
  } else {
    ""
  }

  design_desc <- paste0(
    "The experiment was laid out as a ", design_name,
    " with ", n_trt_str, " and ", n_rep_str, blk_str, ". ",
    "The statistical model was: ", formula_str, "."
  )

  # ---- analysis description ------------------------------------------------
  analysis_desc <- if (!is.null(analysis_obj)) {
    if (!inherits(analysis_obj, "agristat_aov")) {
      rlang::warn("`analysis_obj` is not an `agristat_aov` object; skipping analysis description.")
      ""
    } else {
      resp <- if (!is.null(analysis_obj$response)) analysis_obj$response else "the response variable"
      paste0(
        "Analysis of variance (ANOVA) was performed on ", resp,
        " using the `analyze_design()` function in agristat. ",
        "Treatment differences were considered significant at P <= 0.05."
      )
    }
  } else {
    "Data were analysed by analysis of variance (ANOVA) at P <= 0.05."
  }

  # ---- post-hoc description ------------------------------------------------
  posthoc_desc <- if (!is.null(posthoc_obj)) {
    if (inherits(posthoc_obj, "dmrt_result")) {
      alpha_val <- if (!is.null(posthoc_obj$alpha)) posthoc_obj$alpha else 0.05
      paste0(
        "Significant treatment means were separated using Duncan's Multiple Range Test (DMRT) ",
        "at the ", alpha_val * 100, "% significance level (Duncan, 1955)."
      )
    } else if (inherits(posthoc_obj, "tukey_result") ||
               inherits(posthoc_obj, "tukey_hsd_result")) {
      alpha_val <- if (!is.null(posthoc_obj$alpha)) posthoc_obj$alpha else 0.05
      paste0(
        "Significant treatment means were separated using Tukey's Honestly Significant Difference (HSD) ",
        "test at the ", alpha_val * 100, "% significance level (Tukey, 1949)."
      )
    } else if (inherits(posthoc_obj, "dunnett_result")) {
      paste0(
        "Treatments were compared to the control using Dunnett's test ",
        "(Dunnett, 1955) at P <= 0.05."
      )
    } else if (inherits(posthoc_obj, "mean_separation_result")) {
      method <- if (!is.null(posthoc_obj$method)) posthoc_obj$method else "mean separation"
      paste0(
        "Treatment means were separated using ", method, " at P <= 0.05."
      )
    } else {
      "Post-hoc mean separation was performed at P <= 0.05."
    }
  } else {
    ""
  }

  # ---- software note -------------------------------------------------------
  software_note <- if (isTRUE(software)) {
    paste0(
      "All statistical analyses were carried out in ", .r_version_string(),
      " using the packages: ", .key_packages_string(), ". ",
      "Figures were produced with ggplot2 (Wickham, 2016)."
    )
  } else {
    ""
  }

  # ---- assemble full text --------------------------------------------------
  parts <- c(design_desc, analysis_desc,
             if (nchar(posthoc_desc) > 0) posthoc_desc,
             if (nchar(software_note) > 0) software_note)
  full_text <- paste(parts, collapse = " ")

  refs <- c(
    "R Core Team (2024). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.",
    "Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York."
  )
  if (!is.null(posthoc_obj)) {
    if (inherits(posthoc_obj, "dmrt_result"))
      refs <- c(refs, "Duncan, D.B. (1955). Multiple range and multiple F tests. Biometrics, 11, 1-42.")
    if (inherits(posthoc_obj, "tukey_result") || inherits(posthoc_obj, "tukey_hsd_result"))
      refs <- c(refs, "Tukey, J.W. (1949). Comparing individual means in the analysis of variance. Biometrics, 5(2), 99-114.")
    if (inherits(posthoc_obj, "dunnett_result"))
      refs <- c(refs, "Dunnett, C.W. (1955). A multiple comparison procedure for comparing several treatments with a control. JASA, 50(272), 1096-1121.")
  }

  structure(
    list(
      text     = full_text,
      sections = list(
        design_description   = design_desc,
        analysis_description = paste(c(analysis_desc, posthoc_desc), collapse = " "),
        software_note        = software_note
      ),
      references = refs
    ),
    class = "methods_text"
  )
}

#' @export
print.methods_text <- function(x, ...) {
  cat("=== Methods Text ===\n\n")
  cat(x$text, "\n\n")
  if (length(x$references) > 0L) {
    cat("--- References ---\n")
    cat(paste0("[", seq_along(x$references), "] ", x$references), sep = "\n")
  }
  invisible(x)
}

# ============================================================================
# Step 44: export_results_table()
# ============================================================================

#' Export a Results Table in Various Formats
#'
#' Converts an agristat result object into a formatted table suitable for
#' inclusion in a thesis, report, or publication.  Supported output formats
#' are Markdown, HTML, LaTeX, and CSV.
#'
#' @param results_obj A result object.  Accepted classes:
#'   \describe{
#'     \item{`agristat_aov`}{ANOVA table (Source, Df, SS, MS, F, p).}
#'     \item{`ammi_result`}{AMMI variance partitioning table.}
#'     \item{`gwas_result`}{Top-SNP summary table.}
#'     \item{`dmrt_result`, `tukey_result`, `tukey_hsd_result`,
#'       `dunnett_result`, `mean_separation_result`}{Means and letter
#'       groupings.}
#'   }
#' @param format Character.  One of `"markdown"` (default), `"html"`,
#'   `"latex"`, or `"csv"`.
#' @param file Character or `NULL` (default).  If a file path is supplied the
#'   table is written to that file; otherwise the formatted string is returned.
#' @param digits Integer.  Number of significant digits for numeric columns
#'   (default `4`).
#' @param stars Logical.  If `TRUE` (default) significance stars are appended
#'   to p-value columns where applicable.
#'
#' @return A character string containing the formatted table, or `NULL`
#'   (invisibly) if `file` is specified.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   yield = c(rnorm(5, 4), rnorm(5, 5), rnorm(5, 3)),
#'   trt   = rep(c("A", "B", "C"), each = 5)
#' )
#' d   <- design_crd(c("A", "B", "C"), replicates = 5L, seed = 1L)
#' res <- analyze_design(d, response = "yield", data = d$layout)
#' cat(export_results_table(res))
#'
#' @seealso [generate_methods_text()], [compile_report()]
#' @export
export_results_table <- function(results_obj,
                                 format = "markdown",
                                 file   = NULL,
                                 digits = 4L,
                                 stars  = TRUE) {
  format <- match.arg(format, c("markdown", "html", "latex", "csv"))
  if (!is.null(file) && (!is.character(file) || length(file) != 1L)) {
    rlang::abort("`file` must be a single character string or NULL.")
  }

  df <- .results_to_df(results_obj, digits = digits, stars = stars)
  out <- .format_table(df, format = format)

  if (!is.null(file)) {
    writeLines(out, con = file)
    return(invisible(NULL))
  }
  out
}

.results_to_df <- function(obj, digits, stars) {
  fmt_num <- function(x) {
    if (is.numeric(x)) signif(x, digits) else x
  }

  if (inherits(obj, "agristat_aov")) {
    sm <- tryCatch(summary(obj$model)[[1L]], error = function(e) NULL)
    if (is.null(sm)) {
      rlang::warn("Could not extract ANOVA table from `agristat_aov` object.")
      return(data.frame(Note = "ANOVA table unavailable"))
    }
    p_vals <- sm[["Pr(>F)"]]
    f_vals <- sm[["F value"]]
    df_out <- data.frame(
      Source  = rownames(sm),
      Df      = sm[["Df"]],
      SS      = fmt_num(sm[["Sum Sq"]]),
      MS      = fmt_num(sm[["Mean Sq"]]),
      `F value` = fmt_num(f_vals),
      `p value` = fmt_num(p_vals),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (isTRUE(stars)) {
      df_out[["Sig."]] <- vapply(p_vals, .p_stars, character(1))
    }
    return(df_out)
  }

  if (inherits(obj, "ammi_result")) {
    if (!is.null(obj$anova_table)) {
      return(obj$anova_table)
    }
    if (!is.null(obj$variance_components)) {
      return(obj$variance_components)
    }
    return(data.frame(Note = "AMMI table not extractable from this object"))
  }

  if (inherits(obj, "gwas_result")) {
    gwas_df <- if (!is.null(obj$results)) obj$results else NULL
    if (is.null(gwas_df)) return(data.frame(Note = "GWAS results not available"))
    top_n <- min(20L, nrow(gwas_df))
    keep_cols <- intersect(c("snp", "chr", "pos", "p_value", "effect"),
                           names(gwas_df))
    gwas_sub <- gwas_df[order(gwas_df[["p_value"]])[seq_len(top_n)], keep_cols,
                        drop = FALSE]
    gwas_sub[["p_value"]] <- fmt_num(gwas_sub[["p_value"]])
    if (isTRUE(stars) && "p_value" %in% names(gwas_sub)) {
      raw_p <- obj$results[order(obj$results[["p_value"]])[seq_len(top_n)], "p_value"]
      gwas_sub[["Sig."]] <- vapply(raw_p, .p_stars, character(1))
    }
    return(gwas_sub)
  }

  # Post-hoc / means tables
  mt <- NULL
  if (inherits(obj, c("dmrt_result", "tukey_result", "tukey_hsd_result",
                       "dunnett_result", "mean_separation_result"))) {
    mt <- obj$means_table
  }
  if (!is.null(mt) && is.data.frame(mt)) {
    mt_out <- mt
    num_cols <- sapply(mt_out, is.numeric)
    mt_out[num_cols] <- lapply(mt_out[num_cols], fmt_num)
    return(mt_out)
  }

  rlang::warn("Unrecognised result object type; returning coerced data frame.")
  as.data.frame(unclass(obj))
}

.format_table <- function(df, format) {
  nc   <- ncol(df)
  nr   <- nrow(df)
  cnames <- names(df)

  col_w <- vapply(seq_len(nc), function(j) {
    max(nchar(cnames[j]), max(nchar(as.character(df[[j]])), na.rm = TRUE))
  }, integer(1))

  pad <- function(x, w) formatC(as.character(x), width = w, flag = "-")

  if (format == "csv") {
    lines <- c(paste(cnames, collapse = ","),
               apply(df, 1L, function(r) paste(r, collapse = ",")))
    return(paste(lines, collapse = "\n"))
  }

  if (format == "markdown") {
    header <- paste0("| ", paste(mapply(pad, cnames, col_w), collapse = " | "), " |")
    sep    <- paste0("| ", paste(vapply(col_w, function(w) paste(rep("-", w), collapse = ""), character(1)), collapse = " | "), " |")
    rows   <- apply(df, 1L, function(r) {
      paste0("| ", paste(mapply(pad, r, col_w), collapse = " | "), " |")
    })
    return(paste(c(header, sep, rows), collapse = "\n"))
  }

  if (format == "html") {
    th <- paste0("<th>", cnames, "</th>", collapse = "")
    trs <- apply(df, 1L, function(r) {
      paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")
    })
    return(paste0(
      "<table>\n<thead><tr>", th, "</tr></thead>\n<tbody>\n",
      paste(trs, collapse = "\n"),
      "\n</tbody>\n</table>"
    ))
  }

  if (format == "latex") {
    col_spec <- paste(rep("l", nc), collapse = " ")
    header   <- paste(cnames, collapse = " & ")
    rows_l   <- apply(df, 1L, function(r) paste(r, collapse = " & "))
    return(paste(
      paste0("\\begin{tabular}{", col_spec, "}"),
      "\\hline",
      paste0(header, " \\\\"),
      "\\hline",
      paste(paste0(rows_l, " \\\\"), collapse = "\n"),
      "\\hline",
      "\\end{tabular}",
      sep = "\n"
    ))
  }
}

# ============================================================================
# Step 45: export_figure()
# ============================================================================

#' Export a Plot to a File
#'
#' Saves a ggplot2 object (or the current base-graphics device) to a file in
#' one of several publication-quality formats.
#'
#' @param plot_obj A `ggplot` object, or `NULL` to capture the current
#'   base-graphics device.
#' @param file Character.  Output file path (including extension or just stem
#'   — the correct extension is appended automatically if missing).
#' @param format Character.  One of `"png"` (default), `"pdf"`, `"svg"`,
#'   or `"eps"`.
#' @param width Numeric.  Figure width in inches (default `7`).
#' @param height Numeric.  Figure height in inches (default `5`).
#' @param dpi Integer.  Dots per inch for raster formats (default `300`).
#' @param title Character or `NULL`.  Optional title prepended via
#'   `ggplot2::ggtitle()` before saving (ggplot2 objects only).
#' @param experiment_id Character or `NULL`.  Optional identifier appended to
#'   the figure metadata comment.
#'
#' @return Invisibly, a named list with elements:
#'   \describe{
#'     \item{`file_path`}{Absolute path to the saved file.}
#'     \item{`size_kb`}{File size in kilobytes.}
#'     \item{`resolution`}{DPI value used.}
#'     \item{`format`}{Format string used.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' info <- export_figure(p, file = "scatter.png", format = "png", dpi = 150)
#' info$size_kb
#' }
#'
#' @seealso [compile_report()], [generate_appendix()]
#' @export
export_figure <- function(plot_obj      = NULL,
                          file,
                          format        = "png",
                          width         = 7,
                          height        = 5,
                          dpi           = 300,
                          title         = NULL,
                          experiment_id = NULL) {
  if (missing(file) || !is.character(file) || length(file) != 1L) {
    rlang::abort("`file` must be a single character string.")
  }
  format <- match.arg(format, c("png", "pdf", "svg", "eps"))
  if (!is.numeric(width) || width <= 0)  rlang::abort("`width` must be a positive number.")
  if (!is.numeric(height) || height <= 0) rlang::abort("`height` must be a positive number.")
  if (!is.numeric(dpi) || dpi <= 0)      rlang::abort("`dpi` must be a positive number.")

  # Ensure correct extension
  ext_map <- c(png = "png", pdf = "pdf", svg = "svg", eps = "eps")
  ext     <- ext_map[[format]]
  if (!grepl(paste0("\\.", ext, "$"), file, ignore.case = TRUE)) {
    file <- paste0(file, ".", ext)
  }

  if (inherits(plot_obj, "ggplot")) {
    if (!is.null(title)) {
      plot_obj <- plot_obj + ggplot2::ggtitle(title)
    }
    ggplot2::ggsave(
      filename = file,
      plot     = plot_obj,
      device   = format,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  } else {
    # Base-graphics capture
    dev_fn <- switch(format,
      png  = function(f) grDevices::png(f,  width = width * dpi, height = height * dpi, res = dpi),
      pdf  = function(f) grDevices::pdf(f,  width = width, height = height),
      svg  = function(f) grDevices::svg(f,  width = width, height = height),
      eps  = function(f) grDevices::postscript(f, width = width, height = height,
                                               horizontal = FALSE, paper = "special")
    )
    dev_fn(file)
    if (!is.null(plot_obj)) {
      tryCatch(print(plot_obj), error = function(e) NULL)
    }
    grDevices::dev.off()
  }

  file_path <- normalizePath(file, mustWork = FALSE)
  size_kb   <- if (file.exists(file_path)) {
    round(file.info(file_path)$size / 1024, 2)
  } else {
    NA_real_
  }

  invisible(list(
    file_path  = file_path,
    size_kb    = size_kb,
    resolution = dpi,
    format     = format
  ))
}

# ============================================================================
# Step 46: generate_appendix()
# ============================================================================

#' Generate an Appendix Section
#'
#' Produces a structured appendix comprising a data summary, optional
#' assumption tests (Shapiro-Wilk normality test on residuals), and a raw-data
#' table summary.  All text is formatted in Markdown.
#'
#' @param data A `data.frame` of raw field/experimental data.
#' @param design_obj Optional. A field design object used to contextualise the
#'   appendix header.
#' @param analysis_obj Optional. An `agristat_aov` object.  When provided,
#'   residuals are extracted and tested for normality.
#'
#' @return An S3 object of class `"appendix"` with components:
#'   \describe{
#'     \item{`sections`}{Named list of Markdown-formatted character strings:
#'       `data_summary`, `assumption_tests`, and `raw_data_table`.}
#'     \item{`tables`}{Named list of data frames: `descriptive_stats` and
#'       (if `analysis_obj` provided) `normality_test`.}
#'     \item{`n_obs`}{Integer. Number of observations in `data`.}
#'     \item{`n_missing`}{Integer. Total count of `NA` values in `data`.}
#'   }
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   yield = rnorm(20, 5, 1),
#'   trt   = rep(c("A", "B", "C", "D"), each = 5)
#' )
#' app <- generate_appendix(dat)
#' cat(app$sections$data_summary)
#'
#' @seealso [compile_report()], [generate_methods_text()]
#' @export
generate_appendix <- function(data,
                               design_obj   = NULL,
                               analysis_obj = NULL) {
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }

  n_obs     <- nrow(data)
  n_missing <- sum(is.na(data))

  # ---- data summary --------------------------------------------------------
  num_cols <- names(data)[sapply(data, is.numeric)]
  desc_rows <- lapply(num_cols, function(cn) {
    x <- data[[cn]]
    x_clean <- x[!is.na(x)]
    data.frame(
      Variable = cn,
      N        = length(x_clean),
      Missing  = sum(is.na(x)),
      Mean     = round(mean(x_clean), 4),
      SD       = round(stats::sd(x_clean), 4),
      Min      = round(min(x_clean), 4),
      Max      = round(max(x_clean), 4),
      stringsAsFactors = FALSE
    )
  })
  desc_df <- if (length(desc_rows) > 0L) {
    do.call(rbind, desc_rows)
  } else {
    data.frame(Note = "No numeric columns found")
  }

  data_summary_text <- paste0(
    "## Appendix A: Data Summary\n\n",
    "**Total observations:** ", n_obs, "  \n",
    "**Total missing values:** ", n_missing, "  \n\n",
    "### Descriptive Statistics\n\n",
    export_results_table(
      structure(list(means_table = desc_df), class = "mean_separation_result"),
      format = "markdown"
    ),
    "\n"
  )

  # ---- assumption tests ----------------------------------------------------
  norm_df   <- NULL
  assump_md <- "## Appendix B: Assumption Tests\n\n"

  if (!is.null(analysis_obj)) {
    if (!inherits(analysis_obj, "agristat_aov")) {
      rlang::warn("`analysis_obj` is not an `agristat_aov`; skipping assumption tests.")
    } else {
      resids <- tryCatch(stats::residuals(analysis_obj$model), error = function(e) NULL)
      if (!is.null(resids) && length(resids) >= 3L) {
        sw_test <- stats::shapiro.test(resids)
        norm_df <- data.frame(
          Test      = "Shapiro-Wilk",
          Statistic = round(sw_test$statistic, 4),
          `p value` = round(sw_test$p.value, 4),
          Result    = if (sw_test$p.value >= 0.05) "Normal (fail to reject H0)"
                      else "Non-normal (reject H0 at 5%)",
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        assump_md <- paste0(
          assump_md,
          "### Normality of Residuals (Shapiro-Wilk)\n\n",
          export_results_table(norm_df, format = "markdown"),
          "\n\n",
          if (sw_test$p.value >= 0.05)
            "_Residuals are consistent with normality (p >= 0.05)._\n"
          else
            "_Residuals may depart from normality (p < 0.05). Consider transformations._\n"
        )
      } else {
        assump_md <- paste0(assump_md,
          "_Insufficient residuals for Shapiro-Wilk test (n < 3)._\n")
      }
    }
  } else {
    assump_md <- paste0(assump_md,
      "_No `analysis_obj` provided; assumption tests were not performed._\n")
  }

  # ---- raw data table ------------------------------------------------------
  n_show <- min(20L, n_obs)
  raw_md <- paste0(
    "## Appendix C: Raw Data (first ", n_show, " rows)\n\n",
    export_results_table(
      structure(list(means_table = utils::head(data, n_show)),
                class = "mean_separation_result"),
      format = "markdown"
    ),
    "\n"
  )

  tables <- list(descriptive_stats = desc_df)
  if (!is.null(norm_df)) tables$normality_test <- norm_df

  structure(
    list(
      sections  = list(
        data_summary      = data_summary_text,
        assumption_tests  = assump_md,
        raw_data_table    = raw_md
      ),
      tables    = tables,
      n_obs     = n_obs,
      n_missing = n_missing
    ),
    class = "appendix"
  )
}

#' @export
print.appendix <- function(x, ...) {
  cat("=== Appendix ===\n")
  cat("Observations:", x$n_obs, "| Missing:", x$n_missing, "\n")
  cat("Sections:", paste(names(x$sections), collapse = ", "), "\n")
  invisible(x)
}

# ============================================================================
# Step 47: compile_report()
# ============================================================================

#' Compile a Complete Statistical Analysis Report
#'
#' Assembles a collection of report components (methods text, results tables,
#' figures, appendix) into a single document.  If the **rmarkdown** package is
#' available the document is rendered to the requested output format; otherwise
#' a plain Markdown (`.md`) file is written.
#'
#' @param methods_text_obj Optional. A `methods_text` object from
#'   [generate_methods_text()].
#' @param results_tables Optional. A named list of character strings (formatted
#'   tables) or result objects.  Each element is included as a section.
#' @param figures Optional. A character vector of file paths to figures that
#'   have been saved with [export_figure()].  Each path is embedded as an
#'   image link.
#' @param appendix_obj Optional. An `appendix` object from
#'   [generate_appendix()].
#' @param output_file Character.  Path for the output file (without extension;
#'   the appropriate extension is appended).
#' @param template Character.  One of `"html"` (default), `"pdf"`, or
#'   `"word"`.
#' @param title Character.  Report title (default `"Statistical Analysis
#'   Report"`).
#' @param author Character or `NULL`.  Author name for the YAML front-matter.
#' @param date Character or `NULL`.  Date string; defaults to today's date.
#'
#' @return Invisibly, a named list with elements:
#'   \describe{
#'     \item{`output_file`}{Path to the rendered/written file.}
#'     \item{`format`}{Format string used (`"html"`, `"pdf"`, or `"word"`).}
#'     \item{`sections_included`}{Character vector of sections written.}
#'   }
#'
#' @examples
#' \dontrun{
#' mt  <- generate_methods_text(design_rcbd(c("A","B","C"), blocks = 3L))
#' compile_report(methods_text_obj = mt, output_file = "my_report",
#'                template = "html", title = "Rice Yield Trial")
#' }
#'
#' @seealso [generate_methods_text()], [export_results_table()],
#'   [generate_appendix()]
#' @export
compile_report <- function(methods_text_obj = NULL,
                           results_tables   = NULL,
                           figures          = NULL,
                           appendix_obj     = NULL,
                           output_file,
                           template         = "html",
                           title            = "Statistical Analysis Report",
                           author           = NULL,
                           date             = NULL) {
  if (missing(output_file) || !is.character(output_file) || length(output_file) != 1L) {
    rlang::abort("`output_file` must be a single character string.")
  }
  template <- match.arg(template, c("html", "pdf", "word"))

  date_str   <- if (is.null(date)) format(Sys.Date(), "%d %B %Y") else date
  author_str <- if (is.null(author)) "" else author

  # ---- build YAML front-matter ---------------------------------------------
  yaml_output <- switch(template,
    html = "rmarkdown::html_document",
    pdf  = "rmarkdown::pdf_document",
    word = "rmarkdown::word_document"
  )
  yaml_block <- paste0(
    "---\n",
    "title: \"", title, "\"\n",
    if (nchar(author_str) > 0) paste0("author: \"", author_str, "\"\n") else "",
    "date: \"", date_str, "\"\n",
    "output: ", yaml_output, "\n",
    "---\n\n"
  )

  # ---- sections ------------------------------------------------------------
  sections_included <- character(0)
  body_parts        <- character(0)

  if (!is.null(methods_text_obj)) {
    if (inherits(methods_text_obj, "methods_text")) {
      body_parts <- c(body_parts,
        "## Methods\n\n",
        methods_text_obj$text, "\n\n")
      if (length(methods_text_obj$references) > 0L) {
        body_parts <- c(body_parts,
          "### References\n\n",
          paste0(seq_along(methods_text_obj$references), ". ",
                 methods_text_obj$references, collapse = "\n"), "\n\n")
      }
      sections_included <- c(sections_included, "methods")
    }
  }

  if (!is.null(results_tables)) {
    body_parts <- c(body_parts, "## Results\n\n")
    tbl_list <- if (is.list(results_tables)) results_tables else list(results_tables)
    for (nm in names(tbl_list)) {
      if (!is.null(nm) && nchar(nm) > 0)
        body_parts <- c(body_parts, paste0("### ", nm, "\n\n"))
      tbl_item <- tbl_list[[nm]]
      tbl_text <- if (is.character(tbl_item)) {
        tbl_item
      } else {
        tryCatch(export_results_table(tbl_item, format = "markdown"),
                 error = function(e) paste("*Table could not be formatted:*", e$message))
      }
      body_parts <- c(body_parts, tbl_text, "\n\n")
    }
    sections_included <- c(sections_included, "results")
  }

  if (!is.null(figures) && length(figures) > 0L) {
    body_parts <- c(body_parts, "## Figures\n\n")
    for (fig_path in figures) {
      body_parts <- c(body_parts,
        paste0("![Figure](", fig_path, ")\n\n"))
    }
    sections_included <- c(sections_included, "figures")
  }

  if (!is.null(appendix_obj) && inherits(appendix_obj, "appendix")) {
    body_parts <- c(body_parts,
      "## Appendix\n\n",
      appendix_obj$sections$data_summary, "\n\n",
      appendix_obj$sections$assumption_tests, "\n\n",
      appendix_obj$sections$raw_data_table, "\n\n")
    sections_included <- c(sections_included, "appendix")
  }

  rmd_content <- paste0(yaml_block, paste(body_parts, collapse = ""))

  # ---- write and optionally render -----------------------------------------
  rmd_file <- paste0(output_file, ".Rmd")
  writeLines(rmd_content, con = rmd_file)

  final_file <- rmd_file
  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    tryCatch({
      final_file <- rmarkdown::render(
        input       = rmd_file,
        output_file = basename(output_file),
        output_dir  = dirname(output_file),
        quiet       = TRUE
      )
    }, error = function(e) {
      rlang::warn(paste0("rmarkdown::render() failed: ", e$message,
                         ". Keeping .Rmd file only."))
      final_file <<- rmd_file
    })
  } else {
    md_file <- paste0(output_file, ".md")
    md_content <- sub("^---.*?---\n\n", "", rmd_content, perl = TRUE)
    writeLines(md_content, con = md_file)
    final_file <- md_file
    rlang::inform("rmarkdown not available; written plain Markdown file.")
  }

  invisible(list(
    output_file       = final_file,
    format            = template,
    sections_included = sections_included
  ))
}

# ============================================================================
# Step 48-50: Thesis Templates & fill_template()
# ============================================================================

#' Methods-Section Thesis Template
#'
#' Returns a professional academic **methods section** template for a given
#' experimental design and analysis type.  The template contains placeholders
#' of the form `{placeholder_name}` that can be filled using [fill_template()].
#'
#' @param design_type Character.  One of `"crd"`, `"rcbd"`, `"factorial"`,
#'   `"split_plot"`, or `"latin_square"`.  Default `"rcbd"`.
#' @param analysis_type Character.  One of `"anova"`, `"mixed"`, or `"gwas"`.
#'   Default `"anova"`.
#' @param experiment_type Character.  Context label, e.g. `"plant"` (default)
#'   or `"animal"`.  Affects generic language in the template.
#'
#' @return An S3 object of class `"thesis_template"` with components:
#'   \describe{
#'     \item{`text`}{Character string. Markdown template with placeholders.}
#'     \item{`placeholders`}{Character vector of placeholder names (without
#'       braces).}
#'     \item{`section`}{`"methods"`.}
#'     \item{`design_type`}{The `design_type` argument.}
#'   }
#'
#' @examples
#' tmpl <- thesis_template_methods(design_type = "rcbd")
#' cat(tmpl$text)
#'
#' # Fill placeholders
#' filled <- fill_template(tmpl, list(
#'   n_treatments      = 5,
#'   n_replicates      = 4,
#'   design_description = "The experiment was arranged in an RCBD.",
#'   model_formula      = "y_ij = mu + tau_i + beta_j + epsilon_ij",
#'   software_versions  = "R 4.3.0",
#'   citations          = "R Core Team (2024)."
#' ))
#' cat(filled)
#'
#' @seealso [fill_template()], [generate_methods_text()]
#' @export
thesis_template_methods <- function(design_type    = "rcbd",
                                    analysis_type  = "anova",
                                    experiment_type = "plant") {
  design_type   <- match.arg(tolower(design_type),
                             c("crd", "rcbd", "factorial", "split_plot", "latin_square"))
  analysis_type <- match.arg(tolower(analysis_type),
                             c("anova", "mixed", "gwas"))

  design_human <- switch(design_type,
    crd          = "Completely Randomised Design (CRD)",
    rcbd         = "Randomised Complete Block Design (RCBD)",
    factorial    = "Factorial Design",
    split_plot   = "Split-Plot Design",
    latin_square = "Latin Square Design"
  )

  trt_lang <- if (experiment_type == "animal") "animals" else "plots"

  text <- paste0(
    "## 3. Materials and Methods\n\n",
    "### 3.1 Experimental Design\n\n",
    "The experiment was conducted as a **", design_human, "** with ",
    "{n_treatments} treatments and {n_replicates} replications, ",
    "giving a total of {n_total} experimental ", trt_lang, ". ",
    "{design_description}\n\n",
    "### 3.2 Statistical Analysis\n\n",
    "Data were subjected to analysis of variance (ANOVA) using the ",
    "statistical model:\n\n",
    "$$\n{model_formula}\n$$\n\n",
    "where $\\mu$ is the overall mean, $\\tau_i$ is the treatment effect, ",
    "and $\\varepsilon_{ij}$ is the random error term assumed to be ",
    "independently and normally distributed with mean zero and constant ",
    "variance.\n\n",
    if (analysis_type == "mixed")
      "Mixed-effects models were fitted using the `lme4` package.\n\n"
    else "",
    "Significance was declared at $P \\leq 0.05$. ",
    "Mean separation was performed using an appropriate post-hoc test.\n\n",
    "### 3.3 Software\n\n",
    "All analyses were performed using {software_versions}. {citations}\n"
  )

  placeholders <- c("n_treatments", "n_replicates", "n_total",
                    "design_description", "model_formula",
                    "software_versions", "citations")

  structure(
    list(
      text         = text,
      placeholders = placeholders,
      section      = "methods",
      design_type  = design_type
    ),
    class = "thesis_template"
  )
}

#' Results-Section Thesis Template
#'
#' Returns a professional academic **results section** template for a given
#' analysis type.  Use [fill_template()] to substitute the placeholders.
#'
#' @param analysis_type Character.  One of `"anova"` (default), `"mixed"`,
#'   or `"gwas"`.
#' @param include_posthoc Logical.  If `TRUE` (default) the template includes a
#'   paragraph for mean-separation results.
#'
#' @return An S3 object of class `"thesis_template"` with components:
#'   `text`, `placeholders`, `section` (`"results"`).
#'
#' @examples
#' tmpl <- thesis_template_results()
#' tmpl$placeholders
#'
#' @seealso [thesis_template_methods()], [fill_template()]
#' @export
thesis_template_results <- function(analysis_type  = "anova",
                                    include_posthoc = TRUE) {
  analysis_type <- match.arg(tolower(analysis_type), c("anova", "mixed", "gwas"))

  anova_section <- paste0(
    "### 4.1 Analysis of Variance\n\n",
    "The ANOVA revealed {anova_significance} differences among treatments ",
    "for {response_variable} (Table {table_number}).\n\n",
    "{anova_table}\n\n"
  )

  posthoc_section <- if (include_posthoc) {
    paste0(
      "### 4.2 Mean Separation\n\n",
      "Treatment means ranged from {min_mean} to {max_mean}. ",
      "{means_letters}\n\n",
      "_Figure {figure_number}: {figure_caption}_\n\n"
    )
  } else ""

  interp_section <- paste0(
    "### 4.3 Summary\n\n",
    "{interpretation}\n"
  )

  text <- paste0("## 4. Results\n\n", anova_section, posthoc_section, interp_section)

  placeholders <- c("anova_significance", "response_variable", "table_number",
                    "anova_table", "min_mean", "max_mean", "means_letters",
                    "figure_number", "figure_caption", "interpretation")
  if (!include_posthoc) {
    placeholders <- setdiff(placeholders,
      c("min_mean", "max_mean", "means_letters", "figure_number", "figure_caption"))
  }

  structure(
    list(
      text         = text,
      placeholders = placeholders,
      section      = "results"
    ),
    class = "thesis_template"
  )
}

#' Discussion-Section Thesis Template
#'
#' Returns a professional academic **discussion section** template.
#' Use [fill_template()] to substitute placeholders.
#'
#' @param topic Character.  Topic label used in contextual language.  Default
#'   `"crop_yield"`.
#' @param include_limitations Logical.  If `TRUE` (default) a limitations
#'   subsection is included.
#'
#' @return An S3 object of class `"thesis_template"` with components:
#'   `text`, `placeholders`, `section` (`"discussion"`).
#'
#' @examples
#' tmpl <- thesis_template_discussion(topic = "crop_yield")
#' tmpl$placeholders
#'
#' @seealso [thesis_template_results()], [fill_template()]
#' @export
thesis_template_discussion <- function(topic               = "crop_yield",
                                       include_limitations = TRUE) {
  if (!is.character(topic) || length(topic) != 1L) {
    rlang::abort("`topic` must be a single character string.")
  }

  topic_phrase <- gsub("_", " ", topic)

  text <- paste0(
    "## 5. Discussion\n\n",
    "### 5.1 Interpretation of Key Findings\n\n",
    "The present study investigated ", topic_phrase,
    " under the experimental conditions described. ",
    "{key_findings}\n\n",
    "### 5.2 Comparison with the Literature\n\n",
    "{comparison_literature}\n\n",
    "### 5.3 Practical Implications\n\n",
    "{implications}\n\n",
    if (include_limitations) {
      paste0(
        "### 5.4 Limitations\n\n",
        "{limitations}\n\n"
      )
    } else "",
    "### 5.5 Conclusion\n\n",
    "{conclusion}\n"
  )

  placeholders <- c("key_findings", "comparison_literature",
                    "implications", "conclusion")
  if (include_limitations) placeholders <- c(placeholders, "limitations")

  structure(
    list(
      text         = text,
      placeholders = placeholders,
      section      = "discussion"
    ),
    class = "thesis_template"
  )
}

#' @export
print.thesis_template <- function(x, ...) {
  cat("=== Thesis Template: ", x$section, " ===\n\n", sep = "")
  cat("Placeholders:", paste0("{", x$placeholders, "}", collapse = ", "), "\n\n")
  cat(x$text, "\n")
  invisible(x)
}

#' Fill Placeholders in a Thesis Template
#'
#' Replaces `{placeholder}` tokens in a `thesis_template` object with the
#' supplied values, returning a single character string ready for inclusion in
#' a report or document.
#'
#' @param template_obj A `thesis_template` object produced by
#'   [thesis_template_methods()], [thesis_template_results()], or
#'   [thesis_template_discussion()].
#' @param values A named list mapping placeholder names (without braces) to
#'   replacement character strings (or values coercible to character).
#'
#' @return A character string with all matching placeholders replaced.
#'   Unmatched placeholders are left as-is.
#'
#' @examples
#' tmpl <- thesis_template_methods(design_type = "crd")
#' filled <- fill_template(tmpl, list(
#'   n_treatments       = 4,
#'   n_replicates       = 5,
#'   n_total            = 20,
#'   design_description = "Treatments were randomly assigned to pots.",
#'   model_formula      = "y_{ij} = mu + tau_i + epsilon_{ij}",
#'   software_versions  = "R 4.3.0 (agristat v1.0.0)",
#'   citations          = "R Core Team (2024)."
#' ))
#' cat(filled)
#'
#' @seealso [thesis_template_methods()], [thesis_template_results()]
#' @export
fill_template <- function(template_obj, values) {
  if (!inherits(template_obj, "thesis_template")) {
    rlang::abort("`template_obj` must be a `thesis_template` object.")
  }
  if (!is.list(values)) {
    rlang::abort("`values` must be a named list.")
  }
  if (is.null(names(values)) && length(values) > 0L) {
    rlang::abort("`values` must be a named list.")
  }

  text <- template_obj$text
  for (nm in names(values)) {
    placeholder <- paste0("{", nm, "}")
    text <- gsub(placeholder, as.character(values[[nm]]),
                 text, fixed = TRUE)
  }
  text
}

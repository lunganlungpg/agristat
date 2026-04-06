# Tests for Module 6 — Thesis Assembly & Export Suite

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_crd <- function() {
  design_crd(
    treatments = c("Control", "T1", "T2", "T3"),
    replicates = 4L,
    seed       = 42L
  )
}

make_rcbd <- function() {
  design_rcbd(
    treatments = c("A", "B", "C", "D"),
    blocks     = 4L,
    seed       = 1L
  )
}

make_aov <- function() {
  set.seed(99)
  dat <- data.frame(
    yield = c(rnorm(5, 4.2, 0.5), rnorm(5, 5.1, 0.5), rnorm(5, 3.8, 0.5)),
    trt   = rep(c("A", "B", "C"), each = 5)
  )
  analyze_design(dat, response = "yield", treatment = "trt")
}

make_dmrt <- function() {
  means <- c(A = 14.8, B = 10.2, C = 12.5, D = 11.0)
  dmrt(means, mse = 2.0, df = 16, alpha = 0.05, n = 5)
}

make_tukey <- function() {
  means <- c(A = 14.8, B = 10.2, C = 12.5)
  tukey_hsd(means, mse = 2.0, df = 12, alpha = 0.05, n = 5)
}

make_df <- function() {
  set.seed(7)
  data.frame(
    yield = rnorm(20, 5, 1),
    trt   = rep(c("A", "B", "C", "D"), each = 5),
    block = rep(1:4, times = 5)
  )
}

# ---------------------------------------------------------------------------
# generate_methods_text()
# ---------------------------------------------------------------------------

test_that("generate_methods_text: returns methods_text class from crd_design", {
  d  <- make_crd()
  mt <- generate_methods_text(d)
  expect_s3_class(mt, "methods_text")
})

test_that("generate_methods_text: $text is non-empty character", {
  d  <- make_rcbd()
  mt <- generate_methods_text(d)
  expect_true(is.character(mt$text))
  expect_true(nchar(mt$text) > 50)
})

test_that("generate_methods_text: $sections contains required names", {
  d  <- make_rcbd()
  mt <- generate_methods_text(d)
  expect_named(mt$sections,
               c("design_description", "analysis_description", "software_note"),
               ignore.order = TRUE)
})

test_that("generate_methods_text: $references is a character vector", {
  d  <- make_crd()
  mt <- generate_methods_text(d, software = FALSE)
  expect_true(is.character(mt$references))
  expect_true(length(mt$references) >= 1L)
})

test_that("generate_methods_text: software = FALSE omits software note", {
  d  <- make_crd()
  mt <- generate_methods_text(d, software = FALSE)
  expect_equal(nchar(mt$sections$software_note), 0L)
})

test_that("generate_methods_text: includes posthoc description when provided", {
  d  <- make_rcbd()
  ph <- make_dmrt()
  mt <- generate_methods_text(d, posthoc_obj = ph)
  expect_true(grepl("Duncan", mt$text))
})

test_that("generate_methods_text: includes Tukey description when tukey_hsd provided", {
  d  <- make_rcbd()
  ph <- make_tukey()
  mt <- generate_methods_text(d, posthoc_obj = ph)
  expect_true(grepl("Tukey", mt$text, ignore.case = TRUE))
})

test_that("generate_methods_text: NULL design_obj errors", {
  expect_error(generate_methods_text(NULL))
})

test_that("generate_methods_text: print method works", {
  d  <- make_crd()
  mt <- generate_methods_text(d)
  expect_output(print(mt), "Methods")
})

# ---------------------------------------------------------------------------
# export_results_table()
# ---------------------------------------------------------------------------

test_that("export_results_table: markdown from agristat_aov", {
  res <- make_aov()
  out <- export_results_table(res, format = "markdown")
  expect_true(is.character(out))
  expect_true(grepl("\\|", out))
})

test_that("export_results_table: html format produces <table>", {
  res <- make_aov()
  out <- export_results_table(res, format = "html")
  expect_true(grepl("<table>", out))
  expect_true(grepl("</table>", out))
})

test_that("export_results_table: latex format produces \\begin{tabular}", {
  res <- make_aov()
  out <- export_results_table(res, format = "latex")
  expect_true(grepl("\\\\begin\\{tabular\\}", out))
  expect_true(grepl("\\\\end\\{tabular\\}", out))
})

test_that("export_results_table: csv format does not contain pipes", {
  res <- make_aov()
  out <- export_results_table(res, format = "csv")
  expect_false(grepl("\\|", out))
  expect_true(grepl(",", out))
})

test_that("export_results_table: dmrt_result returns means_table", {
  ph  <- make_dmrt()
  out <- export_results_table(ph, format = "markdown")
  expect_true(is.character(out))
  expect_true(grepl("\\|", out))
})

test_that("export_results_table: tukey_hsd_result returns table", {
  ph  <- make_tukey()
  out <- export_results_table(ph, format = "markdown")
  expect_true(is.character(out))
})

test_that("export_results_table: file argument writes file", {
  res      <- make_aov()
  out_file <- file.path(getwd(), "test_tbl_output.md")
  on.exit(unlink(out_file), add = TRUE)
  result <- export_results_table(res, format = "markdown", file = out_file)
  expect_null(result)
  expect_true(file.exists(out_file))
})

test_that("export_results_table: invalid format errors", {
  res <- make_aov()
  expect_error(export_results_table(res, format = "docx"))
})

# ---------------------------------------------------------------------------
# export_figure()
# ---------------------------------------------------------------------------

test_that("export_figure: saves ggplot2 png and returns list", {
  skip_if_not_installed("ggplot2")
  p        <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  out_file <- file.path(getwd(), "test_fig_export")
  on.exit(unlink(paste0(out_file, ".png")), add = TRUE)
  info <- export_figure(p, file = out_file, format = "png", dpi = 72)
  expect_true(file.exists(paste0(out_file, ".png")))
  expect_named(info, c("file_path", "size_kb", "resolution", "format"),
               ignore.order = TRUE)
  expect_equal(info$format, "png")
  expect_true(info$size_kb > 0)
})

test_that("export_figure: adds extension if missing", {
  skip_if_not_installed("ggplot2")
  p        <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  out_file <- file.path(getwd(), "test_fig_noext")
  on.exit(unlink(paste0(out_file, ".png")), add = TRUE)
  info <- export_figure(p, file = out_file, format = "png", dpi = 72)
  expect_true(endsWith(info$file_path, ".png"))
})

test_that("export_figure: title argument prepended", {
  skip_if_not_installed("ggplot2")
  p        <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  out_file <- file.path(getwd(), "test_fig_title")
  on.exit(unlink(paste0(out_file, ".png")), add = TRUE)
  expect_no_error(
    export_figure(p, file = out_file, format = "png", dpi = 72,
                  title = "My Figure")
  )
})

test_that("export_figure: missing file errors", {
  skip_if_not_installed("ggplot2")
  expect_error(export_figure(NULL))
})

test_that("export_figure: bad format errors", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  expect_error(export_figure(p, file = "out", format = "bmp"))
})

# ---------------------------------------------------------------------------
# generate_appendix()
# ---------------------------------------------------------------------------

test_that("generate_appendix: returns appendix class", {
  app <- generate_appendix(make_df())
  expect_s3_class(app, "appendix")
})

test_that("generate_appendix: $n_obs is correct", {
  dat <- make_df()
  app <- generate_appendix(dat)
  expect_equal(app$n_obs, nrow(dat))
})

test_that("generate_appendix: $n_missing counts NAs", {
  dat        <- make_df()
  dat$yield[3] <- NA
  app <- generate_appendix(dat)
  expect_equal(app$n_missing, 1L)
})

test_that("generate_appendix: $sections has required names", {
  app <- generate_appendix(make_df())
  expect_named(app$sections,
               c("data_summary", "assumption_tests", "raw_data_table"),
               ignore.order = TRUE)
})

test_that("generate_appendix: $tables contains descriptive_stats", {
  app <- generate_appendix(make_df())
  expect_true("descriptive_stats" %in% names(app$tables))
  expect_true(is.data.frame(app$tables$descriptive_stats))
})

test_that("generate_appendix: with analysis_obj adds normality_test", {
  dat <- make_df()
  aov_res <- analyze_design(dat, response = "yield", treatment = "trt")
  app <- generate_appendix(dat, analysis_obj = aov_res)
  expect_true("normality_test" %in% names(app$tables))
})

test_that("generate_appendix: sections are markdown strings", {
  app <- generate_appendix(make_df())
  expect_true(is.character(app$sections$data_summary))
  expect_true(grepl("#", app$sections$data_summary))
})

test_that("generate_appendix: non-data.frame errors", {
  expect_error(generate_appendix(list(a = 1:5)))
})

test_that("generate_appendix: print method works", {
  app <- generate_appendix(make_df())
  expect_output(print(app), "Appendix")
})

# ---------------------------------------------------------------------------
# compile_report()
# ---------------------------------------------------------------------------

test_that("compile_report: creates .Rmd file", {
  d   <- make_rcbd()
  mt  <- generate_methods_text(d)
  out <- file.path(getwd(), "test_report_out")
  on.exit(unlink(c(paste0(out, ".Rmd"), paste0(out, ".md"),
                   paste0(out, ".html")), force = TRUE),
          add = TRUE)
  res <- compile_report(methods_text_obj = mt, output_file = out,
                        template = "html", title = "Test Report")
  expect_true(is.list(res))
  expect_true("output_file" %in% names(res))
  expect_true(file.exists(paste0(out, ".Rmd")) ||
              file.exists(paste0(out, ".html")) ||
              file.exists(paste0(out, ".md")))
})

test_that("compile_report: sections_included records methods", {
  d   <- make_rcbd()
  mt  <- generate_methods_text(d)
  out <- file.path(getwd(), "test_report_sec")
  on.exit(unlink(c(paste0(out, ".Rmd"), paste0(out, ".md"),
                   paste0(out, ".html")), force = TRUE),
          add = TRUE)
  res <- compile_report(methods_text_obj = mt, output_file = out)
  expect_true("methods" %in% res$sections_included)
})

test_that("compile_report: results_tables section recorded", {
  aov_res <- make_aov()
  tbl     <- export_results_table(aov_res, format = "markdown")
  out     <- file.path(getwd(), "test_report_tbl")
  on.exit(unlink(c(paste0(out, ".Rmd"), paste0(out, ".md"),
                   paste0(out, ".html")), force = TRUE),
          add = TRUE)
  res <- compile_report(results_tables = list(ANOVA = tbl),
                        output_file = out)
  expect_true("results" %in% res$sections_included)
})

test_that("compile_report: missing output_file errors", {
  expect_error(compile_report())
})

test_that("compile_report: invalid template errors", {
  expect_error(compile_report(output_file = "x", template = "docx"))
})

# ---------------------------------------------------------------------------
# thesis_template_methods()
# ---------------------------------------------------------------------------

test_that("thesis_template_methods: returns thesis_template class", {
  tmpl <- thesis_template_methods()
  expect_s3_class(tmpl, "thesis_template")
})

test_that("thesis_template_methods: section is 'methods'", {
  tmpl <- thesis_template_methods()
  expect_equal(tmpl$section, "methods")
})

test_that("thesis_template_methods: $text is character", {
  tmpl <- thesis_template_methods(design_type = "crd")
  expect_true(is.character(tmpl$text))
  expect_true(nchar(tmpl$text) > 20)
})

test_that("thesis_template_methods: $placeholders is character vector", {
  tmpl <- thesis_template_methods()
  expect_true(is.character(tmpl$placeholders))
  expect_true(length(tmpl$placeholders) >= 3L)
})

test_that("thesis_template_methods: text contains placeholder tokens", {
  tmpl <- thesis_template_methods()
  for (ph in tmpl$placeholders) {
    expect_true(grepl(paste0("{", ph, "}"), tmpl$text, fixed = TRUE),
                label = paste("placeholder", ph, "in template text"))
  }
})

test_that("thesis_template_methods: all design types accepted", {
  for (dt in c("crd", "rcbd", "factorial", "split_plot", "latin_square")) {
    expect_no_error(thesis_template_methods(design_type = dt))
  }
})

test_that("thesis_template_methods: invalid design_type errors", {
  expect_error(thesis_template_methods(design_type = "invalid_design"))
})

test_that("thesis_template_methods: print method works", {
  tmpl <- thesis_template_methods()
  expect_output(print(tmpl), "methods")
})

# ---------------------------------------------------------------------------
# thesis_template_results()
# ---------------------------------------------------------------------------

test_that("thesis_template_results: returns thesis_template with section 'results'", {
  tmpl <- thesis_template_results()
  expect_s3_class(tmpl, "thesis_template")
  expect_equal(tmpl$section, "results")
})

test_that("thesis_template_results: include_posthoc = FALSE reduces placeholders", {
  t1 <- thesis_template_results(include_posthoc = TRUE)
  t2 <- thesis_template_results(include_posthoc = FALSE)
  expect_true(length(t1$placeholders) > length(t2$placeholders))
})

test_that("thesis_template_results: text contains anova_table placeholder", {
  tmpl <- thesis_template_results()
  expect_true(grepl("{anova_table}", tmpl$text, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# thesis_template_discussion()
# ---------------------------------------------------------------------------

test_that("thesis_template_discussion: returns thesis_template with section 'discussion'", {
  tmpl <- thesis_template_discussion()
  expect_s3_class(tmpl, "thesis_template")
  expect_equal(tmpl$section, "discussion")
})

test_that("thesis_template_discussion: include_limitations = FALSE omits placeholder", {
  t1 <- thesis_template_discussion(include_limitations = TRUE)
  t2 <- thesis_template_discussion(include_limitations = FALSE)
  expect_true("limitations" %in% t1$placeholders)
  expect_false("limitations" %in% t2$placeholders)
})

test_that("thesis_template_discussion: key_findings placeholder present", {
  tmpl <- thesis_template_discussion()
  expect_true("key_findings" %in% tmpl$placeholders)
  expect_true(grepl("{key_findings}", tmpl$text, fixed = TRUE))
})

test_that("thesis_template_discussion: non-character topic errors", {
  expect_error(thesis_template_discussion(topic = 123))
})

# ---------------------------------------------------------------------------
# fill_template()
# ---------------------------------------------------------------------------

test_that("fill_template: replaces a single placeholder", {
  tmpl   <- thesis_template_methods(design_type = "crd")
  filled <- fill_template(tmpl, list(n_treatments = 4))
  expect_false(grepl("{n_treatments}", filled, fixed = TRUE))
  expect_true(grepl("4", filled))
})

test_that("fill_template: replaces all placeholders", {
  tmpl <- thesis_template_methods(design_type = "rcbd")
  vals <- stats::setNames(
    as.list(rep("REPLACED", length(tmpl$placeholders))),
    tmpl$placeholders
  )
  filled <- fill_template(tmpl, vals)
  for (ph in tmpl$placeholders) {
    expect_false(grepl(paste0("{", ph, "}"), filled, fixed = TRUE),
                 label = paste("placeholder", ph, "removed"))
  }
})

test_that("fill_template: unmatched placeholders left as-is", {
  tmpl   <- thesis_template_methods()
  filled <- fill_template(tmpl, list(n_treatments = 3))
  expect_true(grepl("{n_replicates}", filled, fixed = TRUE))
})

test_that("fill_template: non-thesis_template errors", {
  expect_error(fill_template(list(text = "hello {x}"), list(x = "world")))
})

test_that("fill_template: unnamed values list errors", {
  tmpl <- thesis_template_methods()
  expect_error(fill_template(tmpl, list("a", "b")))
})

test_that("fill_template: numeric values coerced to character", {
  tmpl   <- thesis_template_methods()
  filled <- fill_template(tmpl, list(n_treatments = 5L))
  expect_true(grepl("5", filled))
})

test_that("fill_template: returns character string", {
  tmpl   <- thesis_template_results()
  filled <- fill_template(tmpl, list(anova_table = "| A | B |"))
  expect_true(is.character(filled))
  expect_length(filled, 1L)
})

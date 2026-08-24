
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)


latex_result_table <- function(
  data,
  settings_cols,
  metric_groups,
  column_names = names(data),
  caption = NULL,
  label = NULL,
  digits = 3,
  align = NULL,
  scale_down = FALSE,
  font_size = NULL,
  latex_hline = "major",
  custom_latex_hline = NULL,
  collapse_cols = seq_along(settings_cols),
  metric_subgroups = NULL,
  nested_multirow_vmove = NULL
) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop(
      "Package 'kableExtra' is required to render LaTeX result tables.",
      call. = FALSE
    )
  }
  if (!nrow(data)) return(NULL)
  stopifnot(length(settings_cols) >= 1L)
  stopifnot(length(column_names) == ncol(data))
  stopifnot(sum(unlist(metric_groups)) + length(settings_cols) == ncol(data))
  if (!is.null(metric_subgroups)) {
    stopifnot(
      sum(unlist(metric_subgroups)) + length(settings_cols) == ncol(data)
    )
  }
  stopifnot(
    is.numeric(collapse_cols),
    all(collapse_cols >= 1L & collapse_cols <= ncol(data))
  )
  if (!is.null(nested_multirow_vmove)) {
    stopifnot(
      is.numeric(nested_multirow_vmove),
      length(nested_multirow_vmove) == 1L,
      is.finite(nested_multirow_vmove),
      nested_multirow_vmove >= 0
    )
  }

  if (is.null(align)) align <- rep("c", ncol(data))
  latex_options <- "hold_position"
  if (scale_down) latex_options <- c(latex_options, "scale_down")
  if (!is.null(label)) label <- sub("^tab:", "", label)

  output <- kableExtra::kbl(
    data,
    format = "latex",
    booktabs = TRUE,
    digits = digits,
    align = align,
    col.names = column_names,
    caption = caption,
    label = label,
    escape = FALSE,
    linesep = ""
  )

  # Header rows are prepended, so add the lower grouping level first.
  if (!is.null(metric_subgroups)) {
    output <- kableExtra::add_header_above(
      output,
      c(" " = length(settings_cols), metric_subgroups),
      escape = FALSE
    )
  }
  output <- kableExtra::add_header_above(
    output,
    c(" " = length(settings_cols), metric_groups),
    escape = FALSE
  ) |>
    kableExtra::collapse_rows(
      columns = collapse_cols,
      valign = "middle",
      latex_hline = latex_hline,
      custom_latex_hline = custom_latex_hline
    )

  # Remove kableExtra's automatic offsets first; they are too aggressive for
  # nested booktabs groups.  A smaller higher-level adjustment is added below
  # when requested by the caller.
  multirow_vmove <- paste0(
    "\\[[+-]?[0-9]+(?:\\.[0-9]+)?",
    "\\\\dimexpr\\\\aboverulesep\\+\\\\belowrulesep\\+\\\\cmidrulewidth\\]"
  )
  output[1] <- gsub(multirow_vmove, "", output[1], perl = TRUE)

  # When multiple collapsed levels end on the same row, the earlier
  # \multirow commands are the higher-level groups.  Give only those groups
  # an optional, modest upward adjustment; the lowest-level group remains at
  # LaTeX's natural vertical center.
  if (!is.null(nested_multirow_vmove) && nested_multirow_vmove > 0) {
    modifier <- paste0(
      "[", format(nested_multirow_vmove, trim = TRUE, scientific = FALSE),
      "\\dimexpr\\aboverulesep+\\belowrulesep+\\cmidrulewidth]"
    )
    output_lines <- strsplit(output[1], "\n", fixed = TRUE)[[1]]
    for (line_idx in seq_along(output_lines)) {
      line <- output_lines[line_idx]
      multirow_starts <- gregexpr("\\multirow{", line, fixed = TRUE)[[1]]
      multirow_starts <- multirow_starts[multirow_starts > 0]
      if (length(multirow_starts) < 2L) next

      for (start in rev(head(multirow_starts, -1L))) {
        command_tail <- substring(line, start)
        width_end <- regexpr("}{*}", command_tail, fixed = TRUE)[1]
        if (width_end < 1L) next
        insert_after <- start + width_end + nchar("}{*}") - 2L
        line <- paste0(
          substr(line, 1L, insert_after),
          modifier,
          substring(line, insert_after + 1L)
        )
      }
      output_lines[line_idx] <- line
    }
    output[1] <- paste(output_lines, collapse = "\n")
  }

  output <- kableExtra::kable_styling(
    output,
    latex_options = latex_options,
    position = "center",
    font_size = font_size
  )

  # kableExtra's hold_position option emits [!h], which raises a LaTeX
  # warning because `h` is being tightened without offering an alternative.
  # Keep the table near its source while allowing placement at the page top.
  output <- sub(
    "\\begin{table}[!h]",
    "\\begin{table}[!ht]",
    as.character(output),
    fixed = TRUE
  )
  structure(output, format = "latex", class = "knitr_kable")
}

print_latex_result_table <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  if (is.null(table)) {
    cat("No matching results.\n")
  } else {
    cat(as.character(table), "\n", sep = "")
  }
  invisible(table)
}

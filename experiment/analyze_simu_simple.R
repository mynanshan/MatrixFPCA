library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

fix_lines = function(filename, start_key = "MatFPCA,") {
  # 1) Read all lines from your original file:
  orig_lines <- readLines(filename, warn = FALSE)
  
  # 2) Prepare a place to hold the “fixed” lines:
  fixed_lines <- character(0)
  
  # 3) Process each line in turn:
  for (lin in orig_lines) {
    # Trim whitespace at the ends:
    lin_trim <- trimws(lin)
    
    # If the line is completely blank, skip it:
    if (lin_trim == "") next
    
    # If it contains ">>simple,", we may need to split it:
    if (grepl(start_key, lin_trim, fixed = TRUE)) {
      # Insert a newline before any ">>simple," that does NOT already 
      # appear at the very start of the string.
      # The (?<!^) look‐behind ensures we don’t insert a newline before 
      # the first ">>simple," if it’s already at position 1.
      split_pattern <- paste0("(?<!^)(",start_key,")")
      with_newlines <- gsub(split_pattern, "\n\\1", lin_trim, perl = TRUE)
      
      # Now split on those newlines.  This yields one element per record.
      parts <- strsplit(with_newlines, "\n", fixed = TRUE)[[1]]
      
      # Drop any elements that are empty or only whitespace
      parts <- parts[nzchar(trimws(parts))]
      
      # Append them to our output vector
      fixed_lines <- c(fixed_lines, parts)
      
    } else {
      # If there is no ">>simple," but it’s nonblank, just keep it as is.
      fixed_lines <- c(fixed_lines, lin_trim)
    }
  }
  
  # 4) At this point, `fixed_lines` is a character vector where every 
  #    element begins with ">>simple," (and no line contains more than one).
  #    You can inspect the first few to double-check:
  # print(head(fixed_lines, 10))
  
  # 5) Finally, write the cleaned‐up lines back to disk:
  writeLines(fixed_lines, filename)
}

n = 30
pq = c(10, 20, 50, 100, seq(200, 2000, 200))
rt = 3
cm = 0

dir_path <- file.path("experiment", "results")

res <- do.call(
  rbind, lapply(
    pq, \(p) {
      print(p)
      file_id <- paste0(
        "_n",n,"_p",p,"_q",p,
        "_rt",rt,
        "_cm",as.integer(cm)
      )
      metafile_name <- paste0("meta-simu-simple", file_id, ".txt")
      file_name <- paste0("res-simu-simple", file_id, ".txt")
      if (!file.exists(file.path(dir_path, file_name))) return(NULL)
      header <- read.table(file.path(dir_path, metafile_name))[['V1']]
      fix_lines(file.path(dir_path, file_name))
      read.csv(file.path(dir_path, file_name), col.names = header, header=FALSE)
    }
  )
)

dim_counts <- res |> group_by(p) |> count(name = "num") |> ungroup()
# dims_rm <- dim_counts |> filter(num < median(num)) |> pull(p)
# res <- res |> filter(!(p %in% dims_rm))


metric_name = "mse"

# ggplot(res) +
#   geom_line(aes(x = p, y = normV^2))
# ggplot(res) +
#   geom_point(aes(x = p, y = sumU))

res_long <- res |>
  pivot_longer(
    cols = paste0(metric_name, c("U", "V", "K")),
    names_to = "metric",
    values_to = "value") |> 
  mutate(
    metric = factor(
      metric,
      levels=paste0(metric_name, c("U", "V", "K"))
    )
  )

ggplot(res_long |>
    filter(metric %in% paste0(metric_name, c("U", "V")))) +
  geom_point(aes(x = p, y = value * p)) +
  facet_wrap(~metric, scales = "free") + 
  scale_y_log10() +
  # scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()

ggplot(res_long |>
    filter(metric %in% paste0(metric_name, c("K")))) +
  geom_point(aes(x = p, y = value * log(p))) +
  facet_wrap(~metric, scales = "free") + 
  scale_y_log10() + 
  scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()

ggplot(res_long |>
    filter(metric %in% paste0(metric_name, c("K"))) |> 
  filter(p<=100)) +
  geom_point(aes(x = p, y = value * p^2)) +
  facet_wrap(~metric, scales = "free") + 
  scale_y_log10() + 
  scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()

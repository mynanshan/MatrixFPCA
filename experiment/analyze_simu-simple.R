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

n = 20
pq = c(10, 20, 50, 100, seq(200, 2000, 200))
s1vec = 1:4
rt = 3

dir_path <- file.path("experiment", "results")

res <- list()
for (p in pq) {
  for (s1 in s1vec) {
    message("Dimension:", p, " Sampling level:", s1)
    file_id <- paste0(
      "_n",n,"_p",p,"_q",p,
      "_s",s1,"s",s1+1,"_rt",rt
    )
    metafile_name <- paste0("meta-simu-simple.txt")
    file_name <- paste0("res-simu-simple", file_id, ".txt")
    if (!file.exists(file.path(dir_path, file_name))) return(NULL)
    header <- read.table(file.path(dir_path, metafile_name))[['V1']]
    fix_lines(file.path(dir_path, file_name))
    res[[length(res)+1]] <- read.csv(
      file.path(dir_path, file_name),
      col.names = header, header=FALSE
    )
  }
}
res <- do.call(rbind, res)

res |> group_by(p, s1) |> count(name = "num") |> ungroup()


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
  geom_point(aes(x = p, y = value * p), shape = 1) +
  facet_grid(s1 ~ metric, scales = "free") + 
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()

ggplot(res_long |>
    filter(metric %in% paste0(metric_name, c("K")))) +
  geom_point(aes(x = p, y = value * log(p)), shape = 1) +
  facet_wrap( ~ s1, scales = "free") + 
  scale_y_log10() + 
  # scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()

ggplot(
  res_long |>
    filter(metric %in% paste0(metric_name, c("K")))
      # |> filter(p<=400)
  ) +
  geom_point(aes(x = p, y = value * p^2), shape = 1) +
  facet_wrap(~ s1, scales = "free") + 
  scale_y_log10() + 
  scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()

ggplot(
  res_long |>
    filter(metric %in% paste0(metric_name, c("K")))
      |> filter(p>500)
  ) +
  geom_point(aes(x = p, y = value), shape = 1) +
  facet_wrap(~ s1, scales = "free") + 
  scale_y_log10() + 
  # scale_x_log10() +
  labs(x = "Image Dimension", y = "Estimation Error") +
  theme_bw()



# New version starts here ------------------------------------------------

# analyze_simu-simple.R
# ------------------------------------------------------------
# Interactive script to summarize simulation results from
# experiment/results/simu_r{npc_row}-{npc_col}-{npc_time}_set{eigset}_sd{seed}.txt
#
# User should set:
#   npc_row, npc_col, npc_time, eigset, max_seed
#
# Output:
#   - combined raw table
#   - summary table with one row per (n, p, sig)
#   - columns grouped by weighting scheme and metric
# ------------------------------------------------------------

library(data.table)

# =========================
# User-specified parameters
# =========================
npc_row  <- 30
npc_col  <- 30
npc_time <- 5
eigset   <- 2
max_seed <- 100

results_dir <- file.path("experiment", "results")

# =========================
# Helper
# =========================
safe_read_one <- function(seed, npc_row, npc_col, npc_time, eigset,
                          results_dir = file.path("experiment", "results")) {
  file_id <- paste0(
    "_r", npc_row, "-", npc_col, "-", npc_time,
    "_set", eigset, "_sd", seed
  )
  file_name <- paste0("simu", file_id, ".txt")
  file_path <- file.path(results_dir, file_name)

  if (!file.exists(file_path)) {
    message("Missing file: ", file_path)
    return(NULL)
  }

  dt <- tryCatch(
    fread(file_path),
    error = function(e) {
      message("Failed to read: ", file_path)
      message("  ", e$message)
      return(NULL)
    }
  )

  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  dt[, file_seed := seed]
  dt[, file_name := file_name]
  dt
}

# =========================
# Read and combine
# =========================
dt_list <- lapply(
  X = 1:max_seed,
  FUN = safe_read_one,
  npc_row = npc_row,
  npc_col = npc_col,
  npc_time = npc_time,
  eigset = eigset,
  results_dir = results_dir
)

dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) {
  stop("No result files were successfully read.")
}

res <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# =========================
# Basic checks
# =========================
required_cols <- c(
  "seed", "weight", "n", "p", "q", "sig", "s1", "s2", "nt",
  "relsseU", "relsseV", "reliseK", "relsiseC"
)

missing_cols <- setdiff(required_cols, names(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Coerce numeric columns just in case fread inferred oddly
num_cols <- c(
  "seed", "n", "p", "q", "sig", "s1", "s2", "nt",
  "relsseU", "relsseV", "reliseK", "relsiseC"
)
for (cc in intersect(num_cols, names(res))) {
  res[, (cc) := as.numeric(get(cc))]
}

res[, weight := as.character(weight)]

# Check q = p assumption
bad_q <- res[!is.na(q) & q != p]
if (nrow(bad_q) > 0) {
  warning("Some rows have q != p. The summary still groups only by p.")
}

# Check s1, s2, nt fixed assumption
check_fixed <- unique(res[, .(s1, s2, nt)])
if (nrow(check_fixed) > 1) {
  warning("s1, s2, nt are not fixed across files. ",
          "The summary below still groups only by (n, p, sig, weight).")
}

# =========================
# Convert logged metrics back
# =========================
metric_log_cols <- c("relsseU", "relsseV", "reliseK", "relsiseC")
metric_raw_cols <- paste0(metric_log_cols, "_raw")

for (j in seq_along(metric_log_cols)) {
  log_col <- metric_log_cols[j]
  raw_col <- metric_raw_cols[j]
  res[, (raw_col) := exp(get(log_col))]
}

# =========================
# Summarize mean and sd
# =========================
summary_long <- res[
  ,
  .(
    n_rep = .N,
    relsseU_mean = mean(relsseU_raw, na.rm = TRUE),
    relsseU_sd   = sd(relsseU_raw, na.rm = TRUE),
    relsseV_mean = mean(relsseV_raw, na.rm = TRUE),
    relsseV_sd   = sd(relsseV_raw, na.rm = TRUE),
    reliseK_mean = mean(reliseK_raw, na.rm = TRUE),
    reliseK_sd   = sd(reliseK_raw, na.rm = TRUE),
    relsiseC_mean = mean(relsiseC_raw, na.rm = TRUE),
    relsiseC_sd   = sd(relsiseC_raw, na.rm = TRUE)
  ),
  by = .(n, p, sig, weight)
]

setorder(summary_long, n, p, sig, weight)

# =========================
# Make wide table:
# rows = (n, p, sig)
# cols = weight scheme × metric
# =========================
summary_wide <- dcast(
  summary_long,
  n + p + sig ~ weight,
  value.var = c(
    "relsseU_mean", "relsseU_sd",
    "relsseV_mean", "relsseV_sd",
    "reliseK_mean", "reliseK_sd",
    "relsiseC_mean", "relsiseC_sd",
    "n_rep"
  )
)

setorder(summary_wide, n, p, sig)

# Optional: reorder columns to group by weight first, then metric
id_cols <- c("n", "p", "sig")
other_cols <- setdiff(names(summary_wide), id_cols)

split_name <- function(x) {
  pos <- regexpr("_", x)[1]
  if (pos < 0) return(c(metric = x, weight = ""))
  c(metric = substr(x, 1, pos - 1), weight = substr(x, pos + 1, nchar(x)))
}

meta <- data.table(col = other_cols)
meta[, c("metric", "weight") := transpose(lapply(col, split_name))]
metric_order <- c(
  "n_rep",
  "relsseU_mean", "relsseU_sd",
  "relsseV_mean", "relsseV_sd",
  "reliseK_mean", "reliseK_sd",
  "relsiseC_mean", "relsiseC_sd"
)
weight_order <- c("identity", "invcov", "optimal")

meta[, metric_rank := match(metric, metric_order)]
meta[, weight_rank := match(weight, weight_order)]
setorder(meta, weight_rank, metric_rank, col)

summary_wide <- summary_wide[, c(id_cols, meta$col), with = FALSE]

# =========================
# Pretty rounded version
# =========================
summary_wide_round <- copy(summary_wide)

num_cols_wide <- setdiff(names(summary_wide_round), id_cols)
summary_wide_round[, (num_cols_wide) := lapply(.SD, round, 6), .SDcols = num_cols_wide]

# =========================
# Objects available to user
# =========================
# res                : combined raw table
# summary_long       : grouped long summary
# summary_wide       : grouped wide summary
# summary_wide_round : rounded display table

cat("Loaded ", nrow(res), " rows from ", length(dt_list), " files.\n", sep = "")
cat("Objects created:\n")
cat("  res\n")
cat("  summary_long\n")
cat("  summary_wide\n")
cat("  summary_wide_round\n")

print(summary_wide_round)



# analyze_simu-simple.R
# ------------------------------------------------------------
# Interactive script to summarize simulation results from
# experiment/results/simu_r{npc_row}-{npc_col}-{npc_time}_set{eigset}_sd{seed}.txt
#
# User should set:
#   npc_row, npc_col, npc_time, eigset, max_seed
#
# Output:
#   - combined raw table
#   - summary table with one row per (n, p, sig)
#   - LaTeX table string: latex_table
# ------------------------------------------------------------

library(data.table)
library(xtable)

# =========================
# User-specified parameters
# =========================
npc_row  <- 30
npc_col  <- 30
npc_time <- 5
eigset   <- 2
max_seed <- 100

results_dir <- file.path("experiment", "results")

# =========================
# Helper
# =========================
safe_read_one <- function(seed, npc_row, npc_col, npc_time, eigset,
                          results_dir = file.path("experiment", "results")) {
  file_id <- paste0(
    "_r", npc_row, "-", npc_col, "-", npc_time,
    "_set", eigset, "_sd", seed
  )
  file_name <- paste0("simu", file_id, ".txt")
  file_path <- file.path(results_dir, file_name)

  if (!file.exists(file_path)) {
    message("Missing file: ", file_path)
    return(NULL)
  }

  dt <- tryCatch(
    fread(file_path),
    error = function(e) {
      message("Failed to read: ", file_path)
      message("  ", e$message)
      return(NULL)
    }
  )

  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  dt[, file_seed := seed]
  dt[, file_name := file_name]
  dt
}

# =========================
# Read and combine
# =========================
dt_list <- lapply(
  X = 1:max_seed,
  FUN = safe_read_one,
  npc_row = npc_row,
  npc_col = npc_col,
  npc_time = npc_time,
  eigset = eigset,
  results_dir = results_dir
)

dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) {
  stop("No result files were successfully read.")
}

res <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# =========================
# Basic checks
# =========================
required_cols <- c(
  "seed", "weight", "n", "p", "q", "sig", "s1", "s2", "nt",
  "relsseU", "relsseV", "reliseK", "relsiseC"
)

missing_cols <- setdiff(required_cols, names(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

num_cols <- c(
  "seed", "n", "p", "q", "sig", "s1", "s2", "nt",
  "relsseU", "relsseV", "reliseK", "relsiseC"
)
for (cc in intersect(num_cols, names(res))) {
  res[, (cc) := as.numeric(get(cc))]
}

res[, weight := as.character(weight)]

bad_q <- res[!is.na(q) & q != p]
if (nrow(bad_q) > 0) {
  warning("Some rows have q != p. The summary still groups only by p.")
}

check_fixed <- unique(res[, .(s1, s2, nt)])
if (nrow(check_fixed) > 1) {
  warning("s1, s2, nt are not fixed across files.")
}

# =========================
# Convert logged metrics back
# =========================
res[, relsseU_raw := exp(relsseU)]
res[, relsseV_raw := exp(relsseV)]
res[, reliseK_raw := exp(reliseK)]
res[, relsiseC_raw := exp(relsiseC)]

# =========================
# Summarize means only
# =========================
summary_long <- res[
  ,
  .(
    U = mean(relsseU_raw, na.rm = TRUE),
    V = mean(relsseV_raw, na.rm = TRUE),
    K = mean(reliseK_raw, na.rm = TRUE),
    C = mean(relsiseC_raw, na.rm = TRUE),
    n_rep = .N
  ),
  by = .(n, p, sig, weight)
]

setorder(summary_long, n, p, sig, weight)

# =========================
# Wide table with desired order
# =========================
summary_wide <- dcast(
  summary_long,
  n + p + sig ~ weight,
  value.var = c("U", "V", "K", "C")
)

# Rename to desired style
old_names <- names(summary_wide)
new_names <- old_names

rename_map <- c(
  "U_identity" = "idt_U",
  "V_identity" = "idt_V",
  "K_identity" = "idt_K",
  "C_identity" = "idt_C",
  "U_invcov"   = "inv_U",
  "V_invcov"   = "inv_V",
  "K_invcov"   = "inv_K",
  "C_invcov"   = "inv_C",
  "U_optimal"  = "opt_U",
  "V_optimal"  = "opt_V",
  "K_optimal"  = "opt_K",
  "C_optimal"  = "opt_C"
)

for (nm in names(rename_map)) {
  new_names[old_names == nm] <- rename_map[[nm]]
}
setnames(summary_wide, old = old_names, new = new_names)

desired_cols <- c(
  "n", "p", "sig",
  "idt_U", "idt_V", "idt_K", "idt_C",
  "inv_U", "inv_V", "inv_K", "inv_C",
  "opt_U", "opt_V", "opt_K", "opt_C"
)

missing_desired <- setdiff(desired_cols, names(summary_wide))
if (length(missing_desired) > 0) {
  stop("Missing expected columns after reshaping: ",
       paste(missing_desired, collapse = ", "))
}

summary_wide <- summary_wide[, ..desired_cols]
setorder(summary_wide, n, p, sig)

# Rounded display version
summary_wide_round <- copy(summary_wide)
num_cols_wide <- setdiff(names(summary_wide_round), c("n", "p", "sig"))
summary_wide_round[, (num_cols_wide) := lapply(.SD, round, 6), .SDcols = num_cols_wide]

# =========================
# Build LaTeX table
# =========================
latex_df <- as.data.frame(summary_wide)

colnames(latex_df) <- c(
  "n", "p", "$\\sigma$",
  "U", "V", "K", "C",
  "U", "V", "K", "C",
  "U", "V", "K", "C"
)

align_vec <- c(
  "r",  # rownames
  "r", "r", "r",
  rep("r", 12)
)

# digits format for xtable:
# first entry is for row names (ignored here),
# then one entry per displayed column
digits_vec <- c(
  0,  # row names
  0,  # n as integer
  0,  # p as integer
  2,  # sigma
  4, 4, 4, 4,  # identity: U V K C
  4, 4, 4, 4,  # invcov:   U V K C
  4, 4, 4, 4   # optimal:  U V K C
)

add_to_row <- list()
add_to_row$pos <- list(-1, 0)
add_to_row$command <- c(
  paste0(
    "\\hline\n",
    " &  &  & \\multicolumn{4}{c}{identity}",
    " & \\multicolumn{4}{c}{invcov}",
    " & \\multicolumn{4}{c}{optimal} \\\\\n",
    "\\cline{4-7} \\cline{8-11} \\cline{12-15}\n"
  ),
  "\\hline\n"
)

latex_table <- print(
  xtable(
    latex_df,
    align = align_vec,
    digits = digits_vec,
    caption = paste0(
      "Average relative errors for ",
      "$(r_U,r_V,r_K)=(",
      npc_row, ",", npc_col, ",", npc_time, ")$",
      ", eigset=", eigset, "."
    ),
    label = paste0("tab:simu_r", npc_row, "_", npc_col, "_", npc_time, "_set", eigset)
  ),
  include.rownames = FALSE,
  sanitize.colnames.function = identity,
  add.to.row = add_to_row,
  hline.after = NULL,
  print.results = FALSE
)

cat(latex_table)


library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

n = 50
pq = c(10, 20, 50, 100, 200, 400, 600, 800, 1000)
rt = 3
cm = 0

dir_path <- file.path("experiment", "results")

res <- do.call(
  rbind, lapply(
    pq, \(p) {
      file_id <- paste0(
        "_n",n,"_p",p,"_q",p,
        "_rt",rt,
        "_cm",as.integer(cm)
      )
      metafile_name <- paste0("meta-simu-simple", file_id, ".txt")
      file_name <- paste0("res-simu-simple", file_id, ".txt")
      if (!file.exists(file.path(dir_path, file_name))) return(NULL)
      header <- read.table(file.path(dir_path, metafile_name))[['V1']]
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
  scale_x_log10() +
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

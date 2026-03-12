library(tidyverse)


res <- do.call(
  bind_rows,
  lapply(
    list.files("test/results", full.names = TRUE),
    \(f) read_csv(f, col_names = c("seed", "n", "m", "nb", "s", "mse"))
  )
)

res <- res |> 
  filter(s <= 1) |> 
  mutate(
    seed = as.integer(seed),
    n = as.integer(n),
    m = as.integer(m),
    nb = as.integer(nb),
    p = 1 / s
  ) |> 
  arrange(seed, s, n, m, nb)

mList <- unique(res$m)
nbList <- unique(res$nb)
sList <- unique(res$s)

ggplot(res) +
  geom_point(aes(x = p, y = mse), shape = 1) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(m~nb, labeller = label_both) + 
  theme_bw()

source(file.path("R", "matfpca.R"))

Nis_large <- rep(8L, 1000)
nb_small_space <- choose.nbasis(
  nbasis = "theory",
  ndeg = 3,
  ndiag = sum(Nis_large),
  noffdiag = sum(Nis_large * (Nis_large - 1L)),
  nloc = 51,
  Nis = Nis_large,
  p = 10,
  q = 10
)
nb_large_space <- choose.nbasis(
  nbasis = "theory",
  ndeg = 3,
  ndiag = sum(Nis_large),
  noffdiag = sum(Nis_large * (Nis_large - 1L)),
  nloc = 51,
  Nis = Nis_large,
  p = 100,
  q = 100
)

stopifnot(
  nb_large_space[["K"]] > nb_small_space[["K"]],
  nb_large_space[["G"]] > nb_small_space[["G"]]
)

Nis_run_like <- rep(8L, 50)
nb_run_p30 <- choose.nbasis(
  nbasis = "theory",
  ndeg = 3,
  ndiag = sum(Nis_run_like),
  noffdiag = sum(Nis_run_like * (Nis_run_like - 1L)),
  nloc = 51,
  Nis = Nis_run_like,
  p = 30,
  q = 30
)
nb_run_p200 <- choose.nbasis(
  nbasis = "theory",
  ndeg = 3,
  ndiag = sum(Nis_run_like),
  noffdiag = sum(Nis_run_like * (Nis_run_like - 1L)),
  nloc = 51,
  Nis = Nis_run_like,
  p = 200,
  q = 200
)

stopifnot(
  nb_run_p30[["K"]] >= 15,
  nb_run_p200[["K"]] >= nb_run_p30[["K"]],
  nb_run_p200[["K"]] <= floor(sqrt(sum(Nis_run_like * (Nis_run_like - 1L)) / 8))
)

U_low_rank <- matrix(0, 100, 100)
V_low_rank <- matrix(0, 100, 100)
U_low_rank[1, 1] <- 100
V_low_rank[1, 1] <- 100
nb_low_rank_pilot <- choose.nbasis(
  nbasis = "theory",
  ndeg = 3,
  ndiag = sum(Nis_large),
  noffdiag = sum(Nis_large * (Nis_large - 1L)),
  nloc = 51,
  Nis = Nis_large,
  p = 100,
  q = 100,
  U_tilde = U_low_rank,
  V_tilde = V_low_rank
)

stopifnot(
  nb_low_rank_pilot[["K"]] < nb_large_space[["K"]],
  nb_low_rank_pilot[["G"]] < nb_large_space[["G"]]
)

cat("Temporal basis choice checks passed.\n")

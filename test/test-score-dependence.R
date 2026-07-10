source(file.path("experiment", "gendata-simple.R"))

score_matrix <- function(dependence, dependence_param, seed = 1L) {
  set.seed(seed)
  generated <- gen_data_simple(
    n = 5000,
    p = 2,
    q = 2,
    npc_row = 2,
    npc_col = 2,
    npc_time = 3,
    nt = 9,
    sampling_ratio = NULL,
    sig = 0,
    eigvecSet = 1,
    scores_distr = "normal",
    score_dependence = dependence,
    dependence_param = dependence_param
  )
  matrix(generated$scores@data, nrow = 5000)
}

common_amplitude <- score_matrix("common_amplitude", 1)
common_signed <- score_matrix("common_signed", 0)
heavy_amplitude <- score_matrix("heavy_amplitude", 1.5)

set.seed(17)
ar1_generated <- gen_data_simple(
  n = 3,
  p = 6,
  q = 5,
  npc_row = 6,
  npc_col = 5,
  npc_time = 3,
  nt = 9,
  sampling_ratio = NULL,
  sig = 0,
  eigvalSet = 1,
  spatial_model = "ar1"
)
ar1_U <- 0.7^abs(outer(seq_len(6), seq_len(6), "-"))
ar1_V <- 0.7^abs(outer(seq_len(5), seq_len(5), "-"))

offdiag_max <- function(C) max(abs(C[row(C) != col(C)]))
offdiag_mean <- function(C) mean(C[row(C) != col(C)])

stopifnot(
  # Pairwise covariance remains close to zero in the finite-moment examples.
  offdiag_max(cor(common_amplitude)) < 0.10,
  offdiag_max(cor(common_signed)) < 0.10,

  # Squared scores are dependent, maximally so in the common-signed model.
  offdiag_mean(cor(common_amplitude^2)) > 0.10,
  max(abs(cor(common_signed^2) - 1)) < 1e-10,

  # The AR(1) option is trace-normalized and reproduces its target covariance.
  isTRUE(all.equal(ar1_generated$cov$U, ar1_U, tolerance = 1e-10)),
  isTRUE(all.equal(ar1_generated$cov$V, ar1_V, tolerance = 1e-10)),
  identical(ar1_generated$spatial_model, "ar1"),

  # The infinite-fourth-moment stress test generates visibly more extreme scores.
  max(abs(heavy_amplitude)) > stats::quantile(abs(common_amplitude), 0.999)
)

cat("Dependent-score construction checks passed.\n")

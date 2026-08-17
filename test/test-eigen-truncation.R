source(file.path("R", "matfdata.R"))

make_test_eigensystem <- function(
  lam_phi = c(100, 0.1),
  lam_psi = c(100, 0.1),
  lam_theta = c(6, 3, 1)
) {
  create_matrix_eigobj(
    Phi = diag(length(lam_phi)),
    Psi = diag(length(lam_psi)),
    Theta = diag(length(lam_theta)),
    lam_phi = lam_phi,
    lam_psi = lam_psi,
    lam_theta = lam_theta,
    t_grid = seq(0, 1, length.out = length(lam_theta))
  )
}

retained_ranks <- function(eig) {
  c(
    Row = ncol(eig$funcs$Phi),
    Col = ncol(eig$funcs$Psi),
    Time = ncol(eig$funcs$Theta)
  )
}

eig <- make_test_eigensystem()

# The manuscript rule is the default: keep the full spatial eigensystems and
# apply the PVE threshold only to the temporal eigenvalues.
default_trunc <- trunc.matfeigen(eig, pve = 0.95, verbose = FALSE)
stopifnot(
  identical(retained_ranks(default_trunc), c(Row = 2L, Col = 2L, Time = 3L))
)

# The former greedy product-PVE rule remains available explicitly.
joint_trunc <- trunc.matfeigen(
  eig, pve = 0.95, verbose = FALSE, pve_rule = "joint"
)
stopifnot(
  identical(retained_ranks(joint_trunc), c(Row = 1L, Col = 1L, Time = 3L))
)

# Explicit ranks continue to bypass PVE-based selection.
fixed_trunc <- trunc.matfeigen(
  eig, I = 1, J = 2, K = 1, verbose = FALSE
)
stopifnot(
  identical(retained_ranks(fixed_trunc), c(Row = 1L, Col = 2L, Time = 1L))
)

# Temporal PVE uses the represented nonnegative spectrum, as in the manuscript.
indefinite_eig <- make_test_eigensystem(lam_theta = c(3, 1, -0.5))
nonnegative_trunc <- trunc.matfeigen(
  indefinite_eig, pve = 0.8, verbose = FALSE
)
stopifnot(
  identical(retained_ranks(nonnegative_trunc), c(Row = 2L, Col = 2L, Time = 2L))
)

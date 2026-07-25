# Matrix-valued Functional Data -------------------------------------------

create_matrix_funcdata_grid = function(X, t_grid) {
  # validate t
  stopifnot(is.vector(t_grid) && !is.unsorted(t_grid))
  
  # validate X
  if (!is(X, "Tensor"))
    stop("X has to be an {rTensor} tensor")
  stopifnot(length(dim(X)) == 4L)
  stopifnot(dim(X)[4] == length(t_grid))
  
  # formulate object
  obj = list(X = X, t_grid = t_grid, vdim = dim(X)[2:3])
  class(obj) <- c("matfdata", "matfdata.grid")
  obj
}

create_matrix_funcdata_irreg = function(
    X_list, t_list=NULL, t_grid=NULL, t_idx_list=NULL) {
  # validate time points
  t_idx_obj = make_t_idx_list(t_list, t_grid, t_idx_list)
  t_grid = t_idx_obj$t_grid
  t_idx_list = t_idx_obj$t_idx_list
  nt = length(t_grid)
  ms = sapply(t_idx_list, length)
  
  # validate X
  # X can be
  # (1) a list of {rTensor} tensors
  # (2) a vector of paths to RDS files, each stores an rTensor for one subject
  in_mem = TRUE # flag for whether X are read into the list
  if (is.list(X_list)) {
    if (!all(sapply(X_list, function(x) is(x, "Tensor"))))
      stop("Each X_list[[i]] has to be an {rTensor} tensor")
    if (!all(sapply(X_list, function(x) length(dim(x))) == 3L))
      stop("Each X_list[[i]] has to be a 3d")
    p = dim(X_list[[1]])[1]
    q = dim(X_list[[1]])[2]
    if (!all(ms == sapply(X_list, function(x) dim(x)[3])))
      stop("Not all t_idx_list[[i]] and X_list[[i]] have matched length.")
  } else if (is.vector(X_list) && is.character(X_list)) {
    X1 = readRDS(X_list[1])
    if (!is(X1, "Tensor"))
      stop("Each X_list[[i]] has to be an {rTensor} tensor")
    if (!length(dim(X1)) == 3L)
      stop("Each X_list[[i]] has to be a 3d")
    p = dim(X1[[1]])[1]
    q = dim(X1[[1]])[2]
    in_mem = FALSE
  } else stop("X_list has to be a list")
  
  # formulate object
  obj = list(
    X_list = X_list,
    t_grid = t_grid,
    t_idx_list = t_idx_list,
    vdim = c(p, q),
    in_mem = in_mem
  )
  class(obj) <- c("matfdata", "matfdata.irreg")
  obj
}

# TODO: complete this
make_t_idx_list = function(
    t_list=NULL, t_grid=NULL, t_idx_list=NULL, max_nt=101) {
  if (is.null(t_list)) {
    if (is.null(t_grid) || is.null(t_idx_list)) {
      stop("Both t_grid and t_idx_list cannot be NULL if t_list is NULL.")
    }
    if (!is.vector(t_grid) || is.unsorted(t_grid)) {
      stop("t_grid is not a vector or is not sorted.")
    }
    nt = length(t_grid)
    if (!is.list(t_idx_list)) stop("t_idx_list is not a list.")
    if (!all(is.integer(unlist(t_idx_list))))
      stop("All elements in t_idx_list have to be integers.")
    if (any(unlist(t_idx_list) > nt | unlist(t_idx_list) < 1))
      stop("All elements in t_idx_list have to be in [1,length(t_grid)].")
  } else {
    stop("NOT IMPLEMENTED")
  }
  list(t_grid = t_grid, t_idx_list = t_idx_list)
}

get_t_list <- function(data) {
  stopifnot(is(data, "matfdata.irreg"))
  lapply(data$t_idx_list, \(idx) data$t_grid[idx])
}


# Count, Retrieve, Subset Subjects ------------------------------------------------

count_subj = function(data) UseMethod("count_subj")

count_subj.matfdata.grid = function(data) dim(data$X)[1]
count_subj.matfdata.irreg = function(data) length(data$X_list)

subset_subj = function(data, i) UseMethod("subset_subj")

subset_subj.matfdata.grid = function(data, i) {
  i = as.integer(as.vector(i))
  create_matrix_funcdata_grid(
    X = rTensor::as.tensor(data$X@data[i, , , , drop = FALSE]),
    t_grid = data$t_grid)
}

subset_subj.matfdata.irreg = function(
  data, i, force_in_mem = FALSE) {
  i = as.integer(as.vector(i))
  if ((!data$in_mem) && force_in_mem) {
    X_list = lapply(data$X_list[i], get_subj_tnsr.matfdata.irreg)
  } else {
    X_list = data$X_list[i]
  }
  create_matrix_funcdata_irreg(
    X_list = X_list,
    t_grid = data$t_grid,
    t_idx_list = data$t_idx_list[i])
}

get_subj_tnsr = function(data, i) UseMethod("get_subj_tnsr")

get_subj_tnsr.matfdata.grid = function(data, i) {
  i = as.integer(i)
  stopifnot(length(i) == 1)
  shape = dim(data$X)[2:length(dim(data$X))]
  X_i <- data$X@data[i, , , , drop = FALSE]
  dim(X_i) <- shape
  rTensor::as.tensor(X_i)
}

get_subj_tnsr.matfdata.irreg = function(data, i, shape_check=TRUE) {
  i = as.integer(i)
  stopifnot(length(i) == 1)
  if (data$in_mem) {
    return(data$X_list[[i]])
  } else {
    X_i = readRDS(data$X_list[[i]])
    if (shape_check) {
      # shape check when reading the data
      stopifnot(length(dim(X_i)==3))
      stopifnot(dim(X_i)[1] == data$vdim[1])
      stopifnot(dim(X_i)[2] == data$vdim[2])
    }
    return(X_i)
  }
}

get_subj_times = function(data, i) UseMethod("get_subj_times")

get_subj_times.matfdata.grid = function(data, i) {
  data$t_grid
}

get_subj_times.matfdata.irreg = function(data, i) {
  i = as.integer(i)
  stopifnot(length(i) == 1)
  data$t_grid[data$t_idx_list[[i]]]
}

# Cov. and Eigen. of Matrix-valued FuncData -----------------------------------

create_matrix_covfunc = function(U, V, K, t_grid, sig2 = 0) {
  # validate t
  stopifnot(is.vector(t_grid) && !is.unsorted(t_grid))
  
  # validate U, V, K
  # stopifnot(isSymmetric(U))
  # stopifnot(isSymmetric(V))
  # stopifnot(isSymmetric(K))
  stopifnot(nrow(K) == length(t_grid))
  
  # formulate object
  obj = list(
    U = U, V = V, K = K, t_grid = t_grid,
    vdim = c(nrow(U), nrow(V)),
    sig2 = sig2
  )
  class(obj) <- c("matfcov")
  obj
}

create_matrix_eigfunc = function(Phi, Psi, Theta, t_grid) {
  # validate t
  stopifnot(is.vector(t_grid) && !is.unsorted(t_grid))

  # validate Phi, Psi, Theta
  stopifnot(nrow(Theta) == length(t_grid))

  # formulate object
  obj = list(Phi = Phi, Psi = Psi, Theta = Theta, t_grid = t_grid,
    vdim = c(nrow(Phi), nrow(Psi)))
  class(obj) <- c("matfeigf")
  obj
}

create_matrix_eigval = function(lam_phi, lam_psi, lam_theta) {

  # validate lam_phi, lam_psi, lam_theta
  lam_phi = as.vector(lam_phi)
  lam_psi = as.vector(lam_psi)
  lam_theta = as.vector(lam_theta)

  # formulate object
  obj = list(lam_phi = lam_phi, lam_psi = lam_psi, lam_theta = lam_theta)
  class(obj) <- c("matfeigval")
  obj
}

create_matrix_eigobj = function(
    Phi, Psi, Theta, lam_phi, lam_psi, lam_theta, t_grid) {
  # validate t
  stopifnot(is.vector(t_grid) && !is.unsorted(t_grid))

  # validate Phi, Psi, Theta
  stopifnot(nrow(Theta) == length(t_grid))
  
  # validate Lam
  stopifnot(ncol(Phi) == length(lam_phi))
  stopifnot(ncol(Psi) == length(lam_psi))
  stopifnot(ncol(Theta) == length(lam_theta))

  # formulate object
  obj = list(
    funcs = create_matrix_eigfunc(Phi, Psi, Theta, t_grid),
    values = create_matrix_eigval(lam_phi, lam_psi, lam_theta) 
  )
  class(obj) <- c("matfeigen")
  obj
}


# Evaluate Functions --------------------------------------------------------

evalf = function(f, t, ...) UseMethod("evalf")

evalf.matfcov = function(f, t=NULL, t_idx=NULL, tol=1e-15) {
  stopifnot(is(f, "matfcov"))
  if (is.null(t) && is.null(t_idx)) {
    K = f$K
  } else {
    if (!is.null(t)) {
      t = sort(as.vector(t))
      diffs = abs(outer(f$t_grid, t, "-"))
      t_idx = apply(diffs < tol, 2, function(x) which(x)[1])
      if (any(is.na(t_idx))) {
        stop("Not all elements in t are in the t_grid.")
      }
    } else if (!is.null(t_idx)) {
      stopifnot(is.integer(t_idx))
      t_idx = sort(t_idx)
    }
    K = f$K[t_idx, t_idx, drop=FALSE]
  }
  K %x% f$V %x% f$U
}

evalf.matfeigf = function(f, t=NULL, t_idx=NULL, tol=1e-15) {
  stopifnot(is(f, "matfeigf"))
  if (is.null(t) && is.null(t_idx)) {
    Theta = f$Theta
  } else {
    if (!is.null(t)) {
      t = sort(as.vector(t))
      diffs = abs(outer(f$t_grid, t, "-"))
      t_idx = apply(diffs < tol, 2, function(x) which(x)[1])
      if (any(is.na(t_idx))) {
        stop("Not all elements in t are in the t_grid.")
      }
    } else if (!is.null(t_idx)) {
      stopifnot(is.integer(t_idx))
      t_idx = sort(t_idx)
    }
    Theta = f$Theta[t_idx, , drop=FALSE]
  }
  Theta %x% f$Psi %x% f$Phi
}


# Wrangling Eigensystem --------------------------------------------------

get_eigval_tnsr.matfeigval = function(eigval) {
  stopifnot(is(eigval, "matfeigval"))
  eigval$lam_phi %o% eigval$lam_psi %o% eigval$lam_theta
}

eigen.matfcov = function(C) {
  stopifnot(is(C, "matfcov"))
  eig_U = eigen(C$U)
  eig_V = eigen(C$V)
  eig_K = eigen(C$K)
  nt = length(C$t_grid)
  t_span = diff(range(C$t_grid))
  int_factor = t_span / nt
  obj = create_matrix_eigobj(
    eig_U$vectors, eig_V$vectors, eig_K$vectors / sqrt(int_factor),
    eig_U$values, eig_V$values, eig_K$values * int_factor, C$t_grid
  )
  obj
}

match_signs.matfeigf = function(eigf, eigf0) {
  stopifnot(is(eigf, "matfeigf"))
  stopifnot(is(eigf0, "matfeigf"))
  for (d in c("Phi", "Psi", "Theta")) {
    if (ncol(eigf[[d]]) > ncol(eigf0[[d]])) {
      message(
        "More eigenfunctions in the target than in the reference.\n",
        "Only matching the first ", ncol(eigf0[[d]]), " eigenfunctions.")
    }
    r = min(ncol(eigf[[d]]), ncol(eigf0[[d]]))
    idx = colMeans(eigf[[d]][,1:r,drop=F] * eigf0[[d]][,1:r,drop=F]) < 0
    if (sum(idx) > 0) {
      eigf[[d]][,idx] = -eigf[[d]][,idx]
    }
  }
  eigf
}

match_signs.matfeigen = function(eig, eig0) {
  stopifnot(is(eig, "matfeigen"))
  stopifnot(is(eig0, "matfeigen"))
  eig$funcs = match_signs.matfeigf(eig$funcs, eig0$funcs)
  eig
}

.trunc.matfeigen = function(eig, I, J, K) {
  obj = create_matrix_eigobj(
    eig$funcs$Phi[, 1:I, drop=FALSE],
    eig$funcs$Psi[, 1:J, drop=FALSE],
    eig$funcs$Theta[, 1:K, drop=FALSE],
    eig$values$lam_phi[1:I],
    eig$values$lam_psi[1:J],
    eig$values$lam_theta[1:K],
    eig$funcs$t_grid
  )
  obj
}

.select_temporal_rank = function(lam_theta, pve) {
  # Estimated covariance operators can have small negative eigenvalues.
  # PVE is defined using the represented nonnegative temporal spectrum.
  lam_theta = pmax(as.vector(lam_theta), 0)
  total_var = sum(lam_theta)
  if (!is.finite(total_var) || total_var <= 0) {
    stop("The temporal eigensystem has no positive variance.")
  }
  which(cumsum(lam_theta) / total_var >= pve)[1]
}

.select_joint_ranks = function(eigval, rmax, pve) {
  # Legacy rule: greedily grow a row-column-time rank box according to the
  # largest marginal increase in explained product variance.
  r = rep(1L, 3)
  eigval_tnsr = get_eigval_tnsr.matfeigval(eigval)
  total_var = sum(eigval_tnsr)
  if (!is.finite(total_var) || total_var <= 0) {
    stop("The product eigensystem has no positive total variance.")
  }
  curr_pve = eigval_tnsr[1, 1, 1] / total_var

  while (curr_pve < pve) {
    if (all(r >= rmax)) {
      stop("Maximum rank reached before attaining the requested PVE.")
    }
    var_incrs = c(
      ifelse(
        r[1] >= rmax[1], 0,
        sum(eigval_tnsr[r[1] + 1, 1:r[2], 1:r[3]])
      ),
      ifelse(
        r[2] >= rmax[2], 0,
        sum(eigval_tnsr[1:r[1], r[2] + 1, 1:r[3]])
      ),
      ifelse(
        r[3] >= rmax[3], 0,
        sum(eigval_tnsr[1:r[1], 1:r[2], r[3] + 1])
      )
    )
    idx = which.max(var_incrs)
    r[idx] = r[idx] + 1L
    curr_pve = sum(
      eigval_tnsr[1:r[1], 1:r[2], 1:r[3], drop = FALSE]
    ) / total_var
  }
  r
}

trunc.matfeigen = function(
    eig, pve=0.95, I=NULL, J=NULL, K=NULL, verbose = TRUE,
    pve_rule = c("temporal", "joint")
  ) {
  stopifnot(is(eig, "matfeigen"))
  if (!is.numeric(pve) || length(pve) != 1L || !is.finite(pve) ||
      pve <= 0 || pve > 1) {
    stop("pve must be one finite number in (0, 1].")
  }
  pve_rule = match.arg(pve_rule)
  rmax_phi = ncol(eig$funcs$Phi)
  rmax_psi = ncol(eig$funcs$Psi)
  rmax_theta = ncol(eig$funcs$Theta)
  rmax = c(rmax_phi, rmax_psi, rmax_theta)

  ranks = list(I = I, J = J, K = K)
  valid_rank = function(rank, upper) {
    is.numeric(rank) && length(rank) == 1L && is.finite(rank) &&
      rank == as.integer(rank) && rank >= 1L && rank <= upper
  }
  for (d in seq_along(ranks)) {
    if (!is.null(ranks[[d]]) && !valid_rank(ranks[[d]], rmax[d])) {
      stop(
        names(ranks)[d], " must be an integer between 1 and ", rmax[d], "."
      )
    }
  }

  # Explicit ranks take precedence over PVE selection.
  supplied = !vapply(ranks, is.null, logical(1))
  if (all(supplied)) {
    return(.trunc.matfeigen(eig, as.integer(I), as.integer(J), as.integer(K)))
  }
  if (any(supplied)) {
    stop("I, J, and K must either all be supplied or all be NULL.")
  }

  if (pve_rule == "temporal") {
    r = c(
      rmax_phi,
      rmax_psi,
      .select_temporal_rank(eig$values$lam_theta, pve)
    )
  } else {
    r = .select_joint_ranks(eig$values, rmax, pve)
  }

  if (verbose) {
    cat(
      paste0("Selected rank (", pve_rule, " PVE rule):"),
      "Row =", r[1], "Col =", r[2], "Time =", r[3], "\n"
    )
  }

  .trunc.matfeigen(eig, r[1], r[2], r[3])
}

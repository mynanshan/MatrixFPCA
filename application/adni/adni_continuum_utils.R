project_to_polyline <- function(points, curve_points) {
  points <- as.matrix(points)
  curve_points <- as.matrix(curve_points)
  segment_start <- curve_points[-nrow(curve_points), , drop = FALSE]
  segment_vector <- curve_points[-1L, , drop = FALSE] - segment_start
  segment_length_squared <- rowSums(segment_vector^2)
  valid <- is.finite(segment_length_squared) & segment_length_squared > 1e-14
  segment_start <- segment_start[valid, , drop = FALSE]
  segment_vector <- segment_vector[valid, , drop = FALSE]
  segment_length_squared <- segment_length_squared[valid]
  if (!nrow(segment_start)) stop("The fitted principal curve is degenerate.")

  segment_lengths <- sqrt(segment_length_squared)
  cumulative_length <- c(0, cumsum(segment_lengths))
  projected <- matrix(NA_real_, nrow(points), ncol(points))
  arc_length <- numeric(nrow(points))
  squared_distance <- rep(Inf, nrow(points))
  segment_index <- integer(nrow(points))
  segment_fraction <- numeric(nrow(points))

  for (segment in seq_len(nrow(segment_start))) {
    residual <- sweep(points, 2L, segment_start[segment, ], "-")
    fraction <- as.numeric(
      residual %*% segment_vector[segment, ] / segment_length_squared[[segment]]
    )
    fraction <- pmin(1, pmax(0, fraction))
    candidate <- sweep(
      fraction %o% segment_vector[segment, ],
      2L, segment_start[segment, ], "+"
    )
    candidate_distance <- rowSums((points - candidate)^2)
    improved <- candidate_distance < squared_distance
    if (any(improved)) {
      squared_distance[improved] <- candidate_distance[improved]
      projected[improved, ] <- candidate[improved, , drop = FALSE]
      segment_index[improved] <- segment
      segment_fraction[improved] <- fraction[improved]
      arc_length[improved] <- cumulative_length[[segment]] +
        fraction[improved] * segment_lengths[[segment]]
    }
  }
  list(
    projected = projected,
    arc_length = arc_length,
    squared_distance = squared_distance,
    segment_index = segment_index,
    segment_fraction = segment_fraction,
    curve_arc_length = cumulative_length
  )
}

fit_score_principal_curve <- function(
  score_input, spline_df = 5, max_iterations = 20L,
  tolerance = 1e-4, grid_points = 800L
) {
  score_input <- as.matrix(score_input)
  if (nrow(score_input) < 8L || ncol(score_input) < 2L) {
    stop("Principal-curve fitting requires at least eight rows and two columns.")
  }
  if (any(!is.finite(score_input))) stop("Principal-curve scores must be finite.")
  spline_df <- min(as.numeric(spline_df), nrow(score_input) - 1L)
  if (spline_df < 2) stop("spline_df must be at least two.")
  grid_points <- max(as.integer(grid_points), 4L * nrow(score_input), 200L)

  centered <- scale(score_input, center = TRUE, scale = FALSE)
  linear_fit <- stats::prcomp(centered, center = FALSE, scale. = FALSE)
  coordinate <- as.numeric(linear_fit$x[, 1L])
  coordinate <- (coordinate - min(coordinate)) / diff(range(coordinate))
  total_variation <- sum(centered^2)
  linear_residual <- centered - tcrossprod(
    linear_fit$x[, 1L], linear_fit$rotation[, 1L]
  )
  linear_variance_explained <- 1 - sum(linear_residual^2) / total_variation

  previous_loss <- Inf
  loss_history <- numeric()
  final_projection <- NULL
  curve_scores <- NULL
  curve_coordinate <- NULL
  for (iteration in seq_len(as.integer(max_iterations))) {
    unique_coordinates <- length(unique(signif(coordinate, 12L)))
    effective_df <- min(spline_df, unique_coordinates - 1L)
    if (effective_df < 2) stop("Principal-curve coordinate collapsed.")
    coordinate_grid <- seq(0, 1, length.out = grid_points)
    curve_scores <- vapply(seq_len(ncol(score_input)), function(dimension) {
      smoother <- stats::smooth.spline(
        x = coordinate,
        y = score_input[, dimension],
        df = effective_df
      )
      stats::predict(smoother, coordinate_grid)$y
    }, numeric(grid_points))
    projection <- project_to_polyline(score_input, curve_scores)
    loss <- sum(projection$squared_distance)
    loss_history <- c(loss_history, loss)
    coordinate <- projection$arc_length / max(projection$curve_arc_length)
    final_projection <- projection
    if (is.finite(previous_loss) &&
        abs(previous_loss - loss) / max(previous_loss, 1e-12) < tolerance) {
      break
    }
    previous_loss <- loss
  }

  curve_segment_lengths <- sqrt(rowSums(diff(curve_scores)^2))
  curve_coordinate <- c(0, cumsum(curve_segment_lengths))
  curve_coordinate <- curve_coordinate / max(curve_coordinate)
  final_projection <- project_to_polyline(score_input, curve_scores)
  coordinate <- final_projection$arc_length /
    max(final_projection$curve_arc_length)
  curve_variance_explained <- 1 -
    sum(final_projection$squared_distance) / total_variation

  list(
    coordinate = coordinate,
    projected_scores = final_projection$projected,
    residual_distance = sqrt(final_projection$squared_distance),
    curve_coordinate = curve_coordinate,
    curve_scores = curve_scores,
    variance_explained = curve_variance_explained,
    linear_variance_explained = linear_variance_explained,
    linear_coordinate = as.numeric(linear_fit$x[, 1L]),
    iterations = length(loss_history),
    loss_history = loss_history,
    spline_df = spline_df
  )
}

orient_principal_curve <- function(curve_fit, reference) {
  reference <- as.numeric(reference)
  association <- suppressWarnings(stats::cor(
    curve_fit$coordinate, reference, method = "spearman", use = "complete.obs"
  ))
  if (is.finite(association) && association < 0) {
    curve_fit$coordinate <- 1 - curve_fit$coordinate
    curve_fit$curve_coordinate <- 1 - rev(curve_fit$curve_coordinate)
    curve_fit$curve_scores <- curve_fit$curve_scores[
      nrow(curve_fit$curve_scores):1L, , drop = FALSE
    ]
    curve_fit$projected_scores <- curve_fit$projected_scores
    curve_fit$linear_coordinate <- -curve_fit$linear_coordinate
    association <- -association
  }
  curve_fit$orientation_reference_spearman <- association
  curve_fit
}

interpolate_curve_scores <- function(curve_fit, positions) {
  positions <- pmin(1, pmax(0, as.numeric(positions)))
  vapply(seq_len(ncol(curve_fit$curve_scores)), function(dimension) {
    stats::approx(
      x = curve_fit$curve_coordinate,
      y = curve_fit$curve_scores[, dimension],
      xout = positions,
      ties = "ordered",
      rule = 2
    )$y
  }, numeric(length(positions)))
}

embedding_trustworthiness <- function(high_dimensional, embedding, neighbors) {
  high_dimensional <- as.matrix(high_dimensional)
  embedding <- as.matrix(embedding)
  n_subjects <- nrow(high_dimensional)
  neighbors <- min(as.integer(neighbors), n_subjects - 1L)
  high_distance <- as.matrix(stats::dist(high_dimensional))
  low_distance <- as.matrix(stats::dist(embedding))
  diag(high_distance) <- Inf
  diag(low_distance) <- Inf
  high_ranks <- t(apply(high_distance, 1L, rank, ties.method = "average"))
  low_neighbors <- t(apply(low_distance, 1L, order))[
    , seq_len(neighbors), drop = FALSE
  ]
  penalty <- sum(vapply(seq_len(n_subjects), function(index) {
    ranks <- high_ranks[index, low_neighbors[index, ]]
    sum(pmax(ranks - neighbors, 0))
  }, numeric(1L)))
  1 - 2 * penalty /
    (n_subjects * neighbors * (2 * n_subjects - 3 * neighbors - 1))
}

embedding_neighbor_overlap <- function(high_dimensional, embedding, neighbors) {
  high_distance <- as.matrix(stats::dist(high_dimensional))
  low_distance <- as.matrix(stats::dist(embedding))
  diag(high_distance) <- Inf
  diag(low_distance) <- Inf
  neighbors <- min(as.integer(neighbors), nrow(high_dimensional) - 1L)
  high_neighbors <- t(apply(high_distance, 1L, order))[
    , seq_len(neighbors), drop = FALSE
  ]
  low_neighbors <- t(apply(low_distance, 1L, order))[
    , seq_len(neighbors), drop = FALSE
  ]
  mean(vapply(seq_len(nrow(high_dimensional)), function(index) {
    length(intersect(high_neighbors[index, ], low_neighbors[index, ])) / neighbors
  }, numeric(1L)))
}

compare_score_embeddings <- function(
  score_input, neighbors = 30L, min_dist = 0.20,
  tsne_perplexity = 30, seed = 1L
) {
  score_input <- as.matrix(score_input)
  n_subjects <- nrow(score_input)
  neighbors <- min(as.integer(neighbors), n_subjects - 1L)
  evaluation_neighbors <- min(
    15L, neighbors, max(1L, floor((2L * n_subjects - 2L) / 3L))
  )

  pca_embedding <- stats::prcomp(
    score_input, center = TRUE, scale. = FALSE, rank. = 2L
  )$x[, 1:2, drop = FALSE]
  set.seed(seed)
  umap_embedding <- uwot::umap(
    score_input,
    n_neighbors = neighbors,
    n_components = 2L,
    metric = "euclidean",
    min_dist = min_dist,
    init = "spectral",
    n_threads = 1L,
    n_sgd_threads = 1L,
    fast_sgd = FALSE,
    ret_model = FALSE,
    verbose = FALSE,
    seed = seed
  )
  perplexity <- min(
    as.numeric(tsne_perplexity),
    max(1, floor((n_subjects - 2L) / 3L))
  )
  set.seed(seed + 1L)
  tsne_embedding <- Rtsne::Rtsne(
    score_input,
    dims = 2L,
    pca = FALSE,
    check_duplicates = FALSE,
    perplexity = perplexity,
    theta = 0.5,
    max_iter = 1000L,
    verbose = FALSE
  )$Y
  embeddings <- list(
    PCA = pca_embedding, UMAP = umap_embedding, `t-SNE` = tsne_embedding
  )
  embeddings <- lapply(embeddings, function(value) {
    value <- as.matrix(value)
    scale(value, center = TRUE, scale = FALSE)
  })
  high_distance <- stats::dist(score_input)
  metrics <- data.table::rbindlist(lapply(names(embeddings), function(method) {
    embedding <- embeddings[[method]]
    data.table::data.table(
      method = method,
      neighborhood_size = evaluation_neighbors,
      trustworthiness = embedding_trustworthiness(
        score_input, embedding, evaluation_neighbors
      ),
      neighbor_overlap = embedding_neighbor_overlap(
        score_input, embedding, evaluation_neighbors
      ),
      global_distance_spearman = suppressWarnings(stats::cor(
        high_distance, stats::dist(embedding), method = "spearman"
      ))
    )
  }))
  metrics[, recommended_use := data.table::fcase(
    method == "PCA", "global-geometry check",
    method == "UMAP", "primary local-structure display",
    default = "local-neighborhood sensitivity check"
  )]
  list(
    embeddings = embeddings,
    metrics = metrics,
    umap_neighbors = neighbors,
    tsne_perplexity = perplexity
  )
}

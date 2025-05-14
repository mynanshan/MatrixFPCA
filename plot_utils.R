
tile_plot.matrix = function(matrix_data, axes=FALSE, flip_y_axis=TRUE, ...) {
  matrix_data = t(matrix_data)
  n = ncol(matrix_data)
  if (flip_y_axis) matrix_data = matrix_data[,n:1]
  image(matrix_data, axes=axes, ...)
}

tile_plot.tensor = function(
    tensor_data, row_idx, col_idx, axes=FALSE, flip_y_axis=TRUE, ...) {
  tile_plot.matrix(
    rTensor::unfold(tensor_data, row_idx, col_idx)@data,
    axes = axes, flip_y_axis = flip_y_axis, ...
  )
}

tile_plot.matf = function(X, Xlabels = NULL, label_cex = 1.5, ...) {
  p = dim(X)[1]
  q = dim(X)[2]
  nt = dim(X)[3]
  
  X_unfold = rTensor::unfold(X, 1, 2:3)
  tile_plot.matrix(X_unfold@data, axes=FALSE, ...)
  for (i in seq_len(nt)) {
    x_start <- ((i - 1) * q - 0.5) / (q * nt - 1)
    x_end <- (i * q - 0.5) / (q * nt - 1)
    rect(x_start, -0.5 / (p-1), x_end, 1 + 0.5 / (p-1), border = "black", lwd = 2)
  }
  
  if (is.null(Xlabels)) Xlabels = as.character(1:nt)
  for (i in seq_len(nt)) {
    x_mid <- ((i - 1) * q - 0.5 + 0.1 * q) / (q * nt - 1)
    y_mid <- 0.9 
    text(x_mid, y_mid, labels = Xlabels[i], col = "black", cex = label_cex)
  }
}

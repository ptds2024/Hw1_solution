hankel_matrix_fast <- function(dim_mat, h_seq) {
  # input validation and informative messages
  if (!is.numeric(dim_mat) || dim_mat <= 0 || dim_mat %% 1 != 0) stop("`dim_mat` must be a positive integer.")
  if (!is.numeric(h_seq)) stop("`h_seq` must be a numeric vector.")
  if (length(h_seq) < 2*dim_mat - 1) stop("Length of `h_seq` must be at least 2 * `dim_mat` - 1.")
  message("Computing a Hankel matrix of size ", dim_mat, "x", dim_mat, " without for loops")
  
  x <- rep(seq_len(dim_mat), times=dim_mat) + rep(seq_len(dim_mat)-1L, each=dim_mat)
  hankel <- matrix(h_seq[x], nrow = dim_mat, ncol = dim_mat)
  
  return(hankel)
}

# Construct the Denoised Time Series
denoised_time_series <- function(hankel, rank) {
  # Tests for Hankel matrix and rank
  if (!is.matrix(hankel)) stop("Input must be a matrix.")
  if (rank < 1 || rank > min(nrow(hankel), ncol(hankel))) stop("Rank must be between 1 and the minimum of the number of rows and columns of the Hankel matrix.")
  
  # Compute SVD
  svd_result <- svd(hankel)
  u <- svd_result$u[, 1:rank]
  d <- diag(svd_result$d[1:rank])
  v <- svd_result$v[, 1:rank]
  
  # construct lower ranked Hankel matrix
  lower_rank_hankel <- u %*% d %*% t(v)
  
  # Convert the lower ranked Hankel matrix back to time series by averaging the anti-diagonals
  A <- lower_rank_hankel[nrow(lower_rank_hankel):1, ] # so we can take the diagonal elements of A
  n <- 2 * nrow(A) - 1
  denoised_time_series <- rep(NA_real_, n)
  
  # create an indicator for all diagonals in the matrix
  d <- row(A) - col(A)
  
  # use split to group on these values
  denoised_time_series <- rev(sapply(split(A, d), FUN = mean))
  
  return(denoised_time_series)
}

# Task 5: Evaluate the Denoising
evaluate_denoising <- function(original, denoised) {
  rmse <- sqrt(mean((original - denoised)^2))
  return(rmse)
}

# Example Usage:

# Generate a noisy time series
dim_mat <- 50
n <- 2 * dim_mat - 1 # sample size
set.seed(123)
signal <- sin(seq(0, 2*pi, length.out = n))
ts0 <- signal + rnorm(n, sd = 0.25)

# Construct Hankel matrix from the noisy time series
hankel <- hankel_matrix_fast(dim_mat, ts0)

# Reconstruct the time series using the first few singular values (e.g., rank = 2)
denoised_ts0 <- denoised_time_series(hankel, rank = 2)

# Evaluate the denoising
mse <- evaluate_denoising(signal, denoised_ts0)
print(paste("Root Mean Squared Error (RMSE):", mse))

# Plot the original, noisy, and denoised time series
plot(original, type = "l", col = "blue", ylim = range(c(original, noisy, denoised)), 
     main = "Denoising Time Series", xlab = "Time", ylab = "Value")
lines(noisy, col = "red")
lines(denoised, col = "green")
legend("topright", legend = c("Signal", "Time series", "Denoised time series"), col = c("blue", "red", "green"), lty = 1)


# Application to BJsales
data <- BJsales
data <- data[1:149]
dim_mat <- (length(data) - 1) / 2

# Construct Hankel matrix from the noisy time series
hankel <- hankel_matrix_fast(dim_mat, data)

# Reconstruct the time series using the first few singular values (e.g., rank = 2)
denoised_geyser <- denoised_time_series(hankel, rank = 50)
plot(data, type = "l", col = "blue",  
     main = "Sales Data with Leading Indicator", xlab = "Time", ylab = "Value")
lines(denoised_geyser, col = "green")
legend("topleft", legend = c("Time series", "Denoised time series"), col = c("blue", "green"), lty = 1)


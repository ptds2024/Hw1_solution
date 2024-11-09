# Solution
# 1.
hankel_matrix <- function(dim_mat, h_seq) {
  hankel <- matrix(0, nrow = dim_mat, ncol = dim_mat)
  for (i in 1:dim_mat) {
    for (j in 1:dim_mat) {
      hankel[i, j] = h_seq[i + j - 1]
    }
  }
  
  return(hankel)
}

# Example usage
h_seq <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
hankel_matrix(5, h_seq)


# 2.
hankel_matrix <- function(dim_mat, h_seq) {
  # input validation and informative messages
  if (!is.numeric(dim_mat) || dim_mat <= 0 || dim_mat %% 1 != 0) stop("`dim_mat` must be a positive integer.")
  if (!is.numeric(h_seq)) stop("`h_seq` must be a numeric vector.")
  if (length(h_seq) < 2*dim_mat - 1) stop("Length of `h_seq` must be at least 2 * `dim_mat` - 1.")
  message("Computing a Hankel matrix of size ", dim_mat, "x", dim_mat)
  
  # compute Hankel matrix
  hankel <- matrix(0, nrow = dim_mat, ncol = dim_mat)
  for (i in 1:dim_mat) {
    for (j in 1:dim_mat) {
      hankel[i, j] = h_seq[i + j - 1]
    }
  }
  
  return(hankel)
}

# Example usage
hankel_matrix(-1, 1:9)
hankel_matrix(1.1, 1:9)
hankel_matrix(5, 1:9)
hankel_matrix(6, 1:9)
hankel_matrix(5, letters[1:9])

# 3.
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

# Example usage
hankel_matrix_fast(5, h_seq)
all.equal(hankel_matrix(5, h_seq), hankel_matrix_fast(5, h_seq))

# 4.
library(microbenchmark)

microbenchmark(
  loop_version = hankel_matrix(100, h_seq),
  vectorized_version = hankel_matrix_fast(100, h_seq),
  times = 1000
)

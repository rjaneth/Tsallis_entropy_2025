al_omari_1_estimator <- function(data) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # m-spacing
  data_sorted <- sort(data)
  sum_term <- 0
  
  for (i in 1:n) {
    if (i <= m) {
      omega_i <- 3/2
      diff_term <- data_sorted[i + m] - data_sorted[1]
    } else if (i >= n - m + 1) {
      omega_i <- 3/2
      diff_term <- data_sorted[n] - data_sorted[i - m]
    } else {
      omega_i <- 2
      diff_term <- data_sorted[i + m] - data_sorted[i - m]
    }
    
    sum_term <- sum_term + log((n / (omega_i * m)) * diff_term)
  }
  
  return(sum_term / n)
}

# set.seed(123)
# test_data <- rnorm(10)
# result <- al_omari_1_estimator(test_data)
# print(result)

# Example usage:
# data <- c(1.2, 1.5,3.6, 2.0, 2.1, 2.5, 2.6, 3.0, 3.1, 3.5)
# result <- al_omari_1_estimator(data)
# cat("Al-Omari-1 Estimator:", result, "\n")

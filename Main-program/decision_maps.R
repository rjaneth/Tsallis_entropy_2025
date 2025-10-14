
# decision map

load("./Data_paper/envi_munich1024adaptive_tsallis_FIXED_ETA_intensity_L12_eta_3.00_5x11_B100_lambda_085.Rdata")
calculate_p_values_matrix <- function(data_matrix, mu, sigma) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  p_values_matrix <- matrix(NA, nrow = rows, ncol = cols)
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      test_difference <- data_matrix[i, j]
      
      epsilon <-(test_difference - 0) / (sigma)
      p_value <-  2*pnorm(-abs(epsilon))# 
      
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(p_values_matrix )
}




mean_difference_values <- mean(difference_values, na.rm = TRUE)
sd_difference_values <- sd(difference_values, na.rm = TRUE)


p_values_tsallis_munich_eta3 <- calculate_p_values_matrix(difference_values, mean_difference_values, sd_difference_values)
save(p_values_tsallis_munich_eta3, file = "./Data/p_values_tsallis_munich_eta3.Rdata")


alpha <- 0.05
decisions_tsallis_munich_eta3 <- ifelse(p_values_tsallis_munich_eta3 < alpha, 1L, 0L)



save(decisions_tsallis_munich_eta3,
     file = "./Data/decisions_tsallis_munich_eta3_alpha005.Rdata")
source("imagematrix_visualizer_v2.R")
plot(imagematrix(decisions_tsallis_munich_eta3))

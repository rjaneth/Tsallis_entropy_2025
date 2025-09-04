# =============================================================================
# This script loads the .Rdata file with precomputed statistics,
# calculates p-values from the difference_values matrix,
# and generates visualizations of the p-value maps using the
# imagematrix_visualizer_v2.
# =============================================================================


load("./Data/dublin1100_munich_eta5_B100_eta_5_5x11_085.Rdata")
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


p_values <- calculate_p_values_matrix(difference_values, mean_difference_values, sd_difference_values)
#save(p_values_renyi, file = "./Data/p_values_.Rdata")


source("imagematrix_visualizer_v2.R")

#plot(imagematrix(p_values))
#plot(imagematrix(equalize(difference_values)))

plot(imagematrix(equalize(x))) # visualize the SAR image

plot(imagematrix(p_values< 0.05)) # binary map


# save images in grayscale

# imagematrixPNG(imagematrix(equalize(x)), name = "image1.png")
# imagematrixPNG(imagematrix(equalize(difference_values)), name = "test.png")
#imagematrixPNG(imagematrix(p_values), name="H_pvalue.png")
imagematrixPNG(imagematrix(p_values<0.05), name="H_005_.png")


# save images with color viridis-H gradient

imagematrix_colorPNG(
  imagematrix_color(p_values),
  name = "p-values_image.png",
  palette_option = "viridis-H",
  legend_width_px = 540,
  scale_factor = 2,
  direction = -1
)



     
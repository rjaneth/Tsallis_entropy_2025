# p-values and plot


#load("./Data/L16_envi_mexico_512_renyi_lambda09_B200.Rdata")
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

source("../imagematrix_visualizer.R")

#plot(imagematrix(p_values< 0.05))

#plot(imagematrix(equalize(x)))
# imagematrixPNG(imagematrix(equalize(x)), name = "image.png")

#imagematrixPNG(imagematrix(p_values), name="H_pvalue_image.png")
imagematrixPNG(imagematrix(p_values<0.05), name="H_005_tsallis_Foggia_1900_09_B100_eta_4_5x11.png")

imagematrix_colorPNG(
  imagematrix_color(p_values),
  name = "p-values_tsallis_Foggia_1900_09_B100_eta_4_5x11_n.png",
  palette_option = "viridis-H",
  legend_width_px = 540,
  scale_factor = 2
)


# imagematrix_colorPNG(
#   imagematrix_color(p_values),
#   name = "p-values_shannon_munich01-H.png",
#   palette_option = "viridis-H",
#   legend_width_px = 220,
#   scale_factor = 1
# )

#dublin
# imagematrix_colorPNG(
#   imagematrix_color(p_values),
#   name = "p-values_nyi-dublin_B.png",
#   palette_option = "viridis-B",
#   legend_width_px = 540,
#   scale_factor = 2
# )

#Munich final p-values_renyi-munich-G, p-values_shannon-munich-H
imagematrix_colorPNG(
  imagematrix_color(p_values),
  name = "p-values_tsallis_Foggia_1900_09_B100_eta_4_5x11_n.png",
  palette_option = "viridis-H",
  legend_width_px = 540,
  scale_factor = 2
)


#dublin p-values_renyi-dublin-H1 p-values_tsallis_1100-dublin-H1

# imagematrix_colorPNG(
#   imagematrix_color(p_values),
#   name = "p-values_tsallis_Foggia_1900_085_B100_eta_4_5x11_bar.png",
#   palette_option = "viridis-H",
#   legend_width_px = 280,
#   scale_factor = 2
# )
# 
# # #london
# imagematrix_colorPNG(
#   imagematrix_color(p_values),
#   name = "p-values_tsallis_Foggia_1900_085_B100_eta_4_5x11_H.png",
#   palette_option = "viridis-H",
#   legend_width_px = 1000,
#   scale_factor = 2
# )


#new_orlr

# imagematrix_colorPNG(
#   imagematrix_color(p_values),
#   name = "p-values_renyi-dublin_b10_L16-Hss11s1.png",
#   palette_option = "viridis-H",
#   legend_width_px = 540,
#   scale_factor = 2
# )


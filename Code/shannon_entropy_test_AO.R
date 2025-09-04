# Start the timer
start_time <- Sys.time()

#Estimators
source("./Code/al_omari_1_estimator.R")
source("./Code/bootstrap_al_omari_1_estimator.R")


source("./Code/renyi_entropy_estimator_v1.R")
source("./Code/bootstrap_renyi_entropy_estimator_v1.R")



source("./Code/read_ENVI_images.R")


#Load and read data from SAR images

 x <- myread.ENVI(file='./Code/SAR/envi_dublin_1100_HH/Intensity_HH.img', 
                  headerfile='./Code/SAR/envi_dublin_1100_HH/Intensity_HH.hdr')#dublin


rows <- nrow(x)
cols <- ncol(x)



L <- 16 # Number of looks
B <- 100 # Replications bootstrap
window_size <- 7 # sliding window size


difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
   
    difference_values[i, j] <- bootstrap_al_omari_1_estimator(window_data,B) - (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
                                                             
  }
}



save(difference_values, x, file = "./Data/results_dublin_AO_w7_L16_b100.Rdata")
# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time


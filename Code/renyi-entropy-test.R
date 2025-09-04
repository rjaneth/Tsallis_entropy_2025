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
lambda<- 0.9
window_size <- 7 # sliding window size


difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
    
    difference_values[i, j] <-bootstrap_renyi_entropy_estimator_v1(window_data, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                                                (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
    
    
    #difference_values[i, j]  <-renyi_entropy_estimator_v1(window_data, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    #                                                                   (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
  }
}


# Save the results

save(difference_values, x, file = "./Data/results_dublin_renyi_w7_09_L16_b100.Rdata")#
# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time


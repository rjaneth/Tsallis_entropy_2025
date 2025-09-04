# Function for Generating Samples from the GI0 Distribution

# Input:
#   - L: Looks.
#   - alpha: Roughness parameter.
#   - mu: The mean parameter 
#   - n: Sample size
#
# Output:
#   - Random samples following the GI0 distribution.

#library(invgamma)

gi0_sample <- function(mu, alpha, L,  n) {
  
  X_samples <- rinvgamma(n, shape = -alpha, rate =mu*(-alpha-1))
  
  Y_samples <- rgamma(n, shape = L, rate = L)
  Z_samples <- X_samples * Y_samples
  return(Z_samples)
}



# gi0_sample2 <- function(mu, alpha, L, n) {
#   X_samples <- rgamma(n, shape = -alpha, rate = mu * (-alpha - 1))
#   Y_samples <- rgamma(n, shape = L, rate = L)
#   Z_samples <-  Y_samples/X_samples
#   return(Z_samples)
# }
# 
# 
# set.seed(1234567890, kind = "Mersenne-Twister")
# 
# #c<-gi0_sample(1, -8, 5, 10)
# m<-gi0_sample2(1, -8, 5, 10)

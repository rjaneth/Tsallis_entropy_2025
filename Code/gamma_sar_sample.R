# Function for Gamma SAR Sample Generation
# Description: This function generates random samples from the Gamma SAR distribution
# based on the provided parameters L, mu, and sample size (n).
#
# Input:
#   - L: Shape parameter of the distribution.
#   - mu: Mean parameter of the distribution.
#   - n: Sample size.
#
# Output:
#   - Random samples following the Gamma SAR distribution.

gamma_sar_sample <- function(L, mu, n) {
  samples <- rgamma(n, shape = L, rate = L / mu)
  return(samples)
}


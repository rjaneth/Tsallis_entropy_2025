entropy_renyi_gamma_sar <- function(L, mu, alpha) {
  entropy <- (alpha * lgamma(L) - lgamma(alpha * (L - 1) + 1) +
                (alpha * (L - 1) + 1) * log(alpha)) / (alpha - 1) + log(mu / L)
  return(entropy)
}

entropy_renyi_gamma_sar <- function(L, mu, lambda) {
  entropy <- (lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mu / L)
  return(entropy)
}

# alpha <- 2
# mu <- 1
# L <- 3
# 
# # Cálculo de la entropía de Rényi
# H <- entropy_renyi_gamma_sar( L, mu, alpha)
# print(H)

# entropy_renyi_gamma_sar_m <- function(L, mu, alpha) {
#   entropy <- (alpha * lgamma(L) - lgamma(alpha * (L - 1) + 1) +
#                 (alpha * (L - 1) + 1) * log(alpha)) / (alpha - 1) + log(mu / L)
#   return(entropy)
# }
# 
# entropy_renyi_gamma_sar_m(5, 10, 2000)

# entropy_renyi_gamma_sar <- function(L, mu, lambda) {
#   
#   # Verificar las condiciones: lambda > 0 y lambda*L - lambda + 1 > 0
#   if (lambda <= 0) {
#     stop("The lambda parameter must be greater than 0")
#   }
#   if (lambda * L - lambda + 1 <= 0) {
#     stop("The expression lambda * L - lambda + 1 must be greater than 0")
#   }
#   
#   term1 <- -log(L)
#   term2 <- (1 / (1 - lambda)) * (
#     lgamma(lambda * L - lambda + 1) - 
#       ((lambda * L - lambda + 1) * log(lambda)) -  # Corrección aquí
#       lambda * lgamma(L)
#   )
#   term3 <- log(mu)
#   
#   # Entropía de Rényi
#   return(term1 + term2 + term3)
# }

 # (-log(L)+ (1 / (1 - lambda)) * (lgamma(lambda * L - lambda + 1) -  ((lambda * L - lambda + 1) * log(lambda)) -   lambda * lgamma(L))+log(mu))
#entropy_renyi_gamma_sar(36, 1, 0.999)
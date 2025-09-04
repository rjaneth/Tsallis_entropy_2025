# entropia exacta formula correcta 
rm(list = ls())
entropy_GI0_renyi <- function(alpha, mu, L, lambda) {
  if (lambda <= 0 || lambda == 1) {
    stop("Lambda debe ser mayor que 0 y diferente de 1.")
  }
  
  # Calcular gamma
  gamma <- -mu * (alpha + 1)
  if (gamma <= 0) {
    stop("Gamma debe ser positiva. Verifica los valores de mu y alpha.")
  }
  
  # Calcular a, b y suma
  a <- lambda * (L - 1) + 1
  b <- lambda * (-alpha + 1) - 1
  ab_sum <- lambda * (L - alpha)
  
  # Verificar que los argumentos de las funciones gamma sean positivos
  if (a <= 0 || b <= 0 || ab_sum <= 0) {
    stop("Los argumentos de las funciones gamma deben ser positivos. Verifica los valores de lambda, L y alpha.")
  }
  
  # Calcular los tC)rminos
  term1 <- log(gamma / L)
  
  term2 <- lambda * (lgamma(L - alpha) - lgamma(-alpha) - lgamma(L))
  
  term3 <- lgamma(a)
  term4 <- lgamma(b)
  term5 <- lgamma(ab_sum)
  
  numerator <- term2 + term3 + term4 - term5
  
  entropy <- term1 + numerator / (1 - lambda)
  
  return(entropy)
}

alpha <- -3
mu <- 1.0
L <- 5
lambda <- 0.8

entropy_value <- entropy_GI0_renyi(alpha, mu, L, lambda)
print(entropy_value)

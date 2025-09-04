#optimizado
renyi_entropy_estimator_v1 <- function(data, a) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # Espaciado m
  data_sorted <- sort(data)   # Ordenar los datos
  
  # Calcular c_i para cada valor de i usando vectores
  c_i <- numeric(n)
  c_i[1:m] <- 1 + (0:(m-1)) / m
  c_i[(m+1):(n-m)] <- 2
  c_i[(n-m+1):n] <- 1 + (n - (n-m+1):n) / m
  
  # Manejar los límites y calcular diff_term para cada i
  diff_term <- numeric(n)
  for (i in 1:n) {
    if (i <= m) {
      diff_term[i] <- data_sorted[i + m] - data_sorted[1]
    } else if (i >= n - m + 1) {
      diff_term[i] <- data_sorted[n] - data_sorted[i - m]
    } else {
      diff_term[i] <- data_sorted[i + m] - data_sorted[i - m]
    }
  }
  
  # Calcular la suma para el estimador
  sum_term <- sum(((n / (c_i * m)) * diff_term)^(1 - a))
  
  # Calcular la entropía de Rényi usando la fórmula general
  renyi_entropy <- (1 / (1 - a)) * log(sum_term / n)
  
  return(renyi_entropy)
}
# set.seed(1234567890, kind = "Mersenne-Twister")
# #set.seed(123)
# data <- rnorm(1000)  # Generar una muestra aleatoria de 100 datos
# alpha <- 0.999       # Parámetro de orden de la entropía de Rényi
# renyi_entropy <- renyi_entropy_estimator_v1(data, alpha)
# print(renyi_entropy)
# # set.seed(1234567890, kind = "Mersenne-Twister")
# # # Ejemplo de uso
# # data <- rnorm(100)  # Datos de ejemplo
# # a <- 0.5            # Parámetro de entropía de Rényi
# # entropy <- renyi_entropy_estimator(data, a)
# # print(entropy)

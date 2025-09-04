bootstrap_renyi_entropy_estimator_v1 <- function(x, B, a) {
  
  n <- length(x)
  v.Bootstrap <- numeric(B)
  
  for (b in 1:B) {
    # Muestra bootstrap con reemplazo
    bx <- sample(x, size = n, replace = TRUE)
    
    # Verificar si todos los valores no son iguales para evitar un resultado no informativo
    if (length(unique(bx)) > 1) {
      entropy_result <- renyi_entropy_estimator_v1(bx, a)
      if (is.finite(entropy_result)) {
        v.Bootstrap[b] <- entropy_result
      } else {
        v.Bootstrap[b] <- NA
      }
    } else {
      v.Bootstrap[b] <- NA  # Asignar NA si todos los valores son iguales
    }
  }
  
  # Calcular la entropía de Rényi para la muestra original
  t <- renyi_entropy_estimator_v1(x, a)
  
  # Calcular el promedio de los valores bootstrap válidos
  estimated_mean <- mean(v.Bootstrap, na.rm = TRUE)
  
  # Estimador mejorado
  return(2 * t - estimated_mean)
}


# # Ejemplo de uso
# data <- rnorm(1000)  # Datos de ejemplo
# B <- 200             # Número de iteraciones de bootstrap
# a <- 0.5             # Parámetro de entropía de Rényi
# entropy <- bootstrap_renyi_entropy_estimator(data, B, a)
# print(entropy)

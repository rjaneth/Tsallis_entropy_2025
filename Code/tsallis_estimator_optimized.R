#------------------------------------------------------------
# 1. Estimador optimizado de Tsallis (spacing-based, orden = λ)
#------------------------------------------------------------
tsallis_estimator_optimized <- function(x, lambda, m = NULL) {
  if (lambda <= 0) stop("lambda debe ser > 0 (orden Tsallis)")
  n <- length(x)
  if (n < 2) return(NA_real_)
  
  if (is.null(m))
    m <- if (n < 10) max(1, floor(sqrt(n))) else round(sqrt(n) + 0.5)
  m <- max(1, min(m, floor((n - 1)/2)))
  
  xs <- sort(x, method = "quick")
  i  <- seq_len(n)
  li <- pmax(i - m, 1L)
  ri <- pmin(i + m, n)
  diff <- xs[ri] - xs[li]
  
  ci <- ifelse(i <= m, (m + i - 1)/m,
               ifelse(i >= n - m + 1, (n + m - i)/m, 2))
  
  r <- n * diff / (ci * m)
  valid <- diff > .Machine$double.eps
  if (!any(valid)) return(NA_real_)
  
  s <- sum(r[valid]^(1 - lambda))
  (1 / (lambda - 1)) * (1 - s / n)
}

# #------------------------------------------------------------
# # 1. Estimador optimizado de Tsallis (spacing-based)
# #------------------------------------------------------------
# tsallis_estimator_optimized <- function(x, alpha, m = NULL) {
#   # Validación de parámetros
#   if (alpha <= 0) stop("alpha debe ser > 0 (Tsallis)")
#   n <- length(x)
#   if (n < 2) {
#     warning("Muestra demasiado pequeña (n<2). Devolviendo NA.")
#     return(NA)
#   }
#   
#   # Ajuste automático de m con límites seguros m = round(sqrt(length(x)) + 0.5
#   if (is.null(m)) {
#     m <- if (n < 10) max(1, floor(sqrt(n))) else round(sqrt(n) + 0.5)
#   }
#   m <- max(1, min(m, floor((n - 1)/2)))  # 1 ≤ m ≤ (n-1)/2
#   
#   xs <- sort(x, method = "quick")
#   i <- seq_len(n)
#   
#   # Manejo seguro de índices
#   li <- pmax(i - m, 1L)
#   ri <- pmin(i + m, n)
#   diff <- xs[ri] - xs[li]
#   
#   # Coeficientes con corrección de bordes
#   ci <- ifelse(
#     i <= m, (m + i - 1)/m,
#     ifelse(i >= n - m + 1, (n + m - i)/m, 2)
#   )
#   
#   # Término central con protección numérica
#   r <- n * diff / (ci * m)
#   valid <- diff > .Machine$double.eps
#   if (sum(valid) == 0) return(NA)  # Todos los espaciamientos son cero
#   
#   # Cálculo final con manejo de alpha = 1
#   s <- sum(r[valid]^(1 - alpha), na.rm = TRUE)
#   (1 / (alpha - 1)) * (1 - s / n)
# }
# 

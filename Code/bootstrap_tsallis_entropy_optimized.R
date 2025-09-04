
#------------------------------------------------------------
# 2. Bootstrap bias-corrected para Tsallis (λ)
#------------------------------------------------------------
bootstrap_tsallis_entropy_optimized <- function(x, B = 200L, lambda,
                                                m = round(sqrt(length(x)) + 0.5),
                                                parallel = FALSE)
{
  n  <- length(x)
  t0 <- tsallis_estimator_optimized(x, lambda, m)
  if (is.na(t0)) return(NA_real_)
  
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  f   <- function(col) tsallis_estimator_optimized(x[col], lambda, m)
  
  boot_vals <- if (parallel &&
                   requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_apply(idx, 2, f)
  } else {
    apply(idx, 2, f)
  }
  
  valid <- boot_vals[!is.na(boot_vals)]
  if (length(valid) < B/2) return(NA_real_)
  
  2 * t0 - mean(valid)
}

# #------------------------------------------------------------
# # 2. Bootstrap optimizado para Tsallis (sin cambios)
# #------------------------------------------------------------
# bootstrap_tsallis_entropy_optimized <- function(x, B = 200L, alpha,
#                                       m = round(sqrt(length(x)) + 0.5),
#                                       parallel = FALSE) 
# {
#   n <- length(x)
#   
#   # Estimación original con nuevo tsallis_estimator
#   t0 <- tsallis_estimator_optimized (x, alpha, m)
#   if (is.na(t0)) {
#     warning("Estimación original falló. Abortando bootstrap.")
#     return(NA)
#   }
#   
#   # Matriz de índices bootstrap
#   idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
#   
#   # Función auxiliar
#   f_boot <- function(col) tsallis_estimator_optimized (x[col], alpha, m)
#   
#   # Paralelización segura
#   if (parallel) {
#     if (!requireNamespace("future.apply", quietly = TRUE)) {
#       message("Usando modo secuencial (instala future.apply para paralelizar)")
#       boot_vals <- apply(idx, 2, f_boot)
#     } else {
#       boot_vals <- future.apply::future_apply(idx, 2, f_boot)
#     }
#   } else {
#     boot_vals <- apply(idx, 2, f_boot)
#   }
#   
#   # Manejo de réplicas fallidas
#   if (any(is.na(boot_vals))) {
#     valid_vals <- boot_vals[!is.na(boot_vals)]
#     if (length(valid_vals) < B/2) {
#       warning(">50% de réplicas bootstrap fallaron. Resultado no confiable")
#       return(NA)
#     }
#     boot_vals <- valid_vals
#   }
#   
#   # Corrección de sesgo
#   2 * t0 - mean(boot_vals, na.rm = TRUE)
# }
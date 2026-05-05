# =============================================================================
# run_shannon_null_diagnostic_sim.R
#
# Diagnóstico nulo para Shannon con imagen homogénea simulada bajo H0.
# Ejecuta:
#   1) ventanas fijas W = 5, 7, 9, 11
#   2) ventana adaptativa (Park-type)
# y produce:
#   - resultados .RData
#   - tablas de calibración nula
#   - frecuencias de W_map (global e interior)
#   - figuras comparativas de distribuciones
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Librerías
# -----------------------------------------------------------------------------
req_pkgs <- c("future.apply", "future", "ggplot2")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(future.apply)
  library(future)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Plan paralelo
# -----------------------------------------------------------------------------
opt_workers <- getOption("future.workers", default = NA_real_)
env_workers <- suppressWarnings(as.numeric(Sys.getenv("FUTURE_WORKERS", "")))
auto_workers <- tryCatch({
  ac <- future::availableCores()
  if (is.numeric(ac) && length(ac) == 1L) max(1, ac - 1L) else 1L
}, error = function(e) 1L)

workers <- if (!is.na(opt_workers) && opt_workers >= 1) {
  as.integer(opt_workers)
} else if (!is.na(env_workers) && env_workers >= 1) {
  as.integer(env_workers)
} else {
  as.integer(auto_workers)
}

parallel <- TRUE
try({
  if (workers > 1L) {
    plan(multisession, workers = workers)
  } else {
    plan(sequential)
    parallel <- FALSE
  }
}, silent = TRUE)

set.seed(1234567890, kind = "Mersenne-Twister")
message(sprintf("Plan: %s con %d worker(s).",
                if (parallel) "multisession" else "sequential",
                future::nbrOfWorkers()))

# -----------------------------------------------------------------------------
# 2. Parámetros
# -----------------------------------------------------------------------------
sim_rdata  <- "L16_simulated_image_1000_gamma_H0.Rdata"
sim_var    <- "Z"

data_scale    <- "intensity"
L             <- 16
use_bootstrap <- TRUE
B             <- 100

fixed_windows <- c(5, 7, 9, 11)

W_min <- 5
W_max <- 11
eta   <- 2

outdir <- "./Data/null_diagnostic_shannon"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 3. Funciones externas
# -----------------------------------------------------------------------------
source("./Code/gamma_sar_sample.R")
source("./Code/al_omari_1_estimator.R")

# -----------------------------------------------------------------------------
# 4. Helper: cargar imagen simulada
# -----------------------------------------------------------------------------
load_simulated_matrix <- function(rdata_path, var_name = NULL) {
  stopifnot(file.exists(rdata_path))
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  
  pick_object <- function(obj) {
    if (is.data.frame(obj)) obj <- as.matrix(obj)
    if (is.array(obj) && length(dim(obj)) > 2) obj <- obj[, , 1, drop = FALSE]
    if (!is.matrix(obj)) obj <- tryCatch(as.matrix(obj), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    storage.mode(obj) <- "double"
    obj
  }
  
  if (!is.null(var_name)) {
    if (!exists(var_name, envir = env, inherits = FALSE)) {
      stop(sprintf("Objeto '%s' no encontrado en '%s'.", var_name, rdata_path))
    }
    obj <- get(var_name, envir = env, inherits = FALSE)
    obj <- pick_object(obj)
    if (is.null(obj)) stop(sprintf("Objeto '%s' no es una matriz 2D numérica.", var_name))
    return(obj)
  }
  
  nms <- ls(envir = env, all.names = TRUE)
  for (nm in nms) {
    cand <- get(nm, envir = env, inherits = FALSE)
    obj  <- pick_object(cand)
    if (!is.null(obj) && is.matrix(obj) && all(dim(obj) >= 3)) {
      return(obj)
    }
  }
  
  stop(sprintf("No se encontró una matriz 2D numérica en '%s'.", rdata_path))
}

# -----------------------------------------------------------------------------
# 5. Cargar imagen
# -----------------------------------------------------------------------------
x <- load_simulated_matrix(sim_rdata, sim_var)

if (!is.matrix(x) || any(dim(x) < 3)) {
  stop("La imagen 'x' debe ser una matriz 2D numérica con al menos 3x3 píxeles.")
}

rows <- nrow(x)
cols <- ncol(x)

cat(sprintf("Imagen simulada homogénea cargada: %d x %d\n", rows, cols))

sigma_n <- if (identical(tolower(data_scale), "amplitude")) {
  0.523 / sqrt(L)
} else if (identical(tolower(data_scale), "intensity")) {
  1 / sqrt(L)
} else {
  stop("data_scale debe ser 'intensity' o 'amplitude'.")
}

cat(sprintf("Escala: %s | L = %d | bootstrap = %s | B = %d\n",
            data_scale, L, ifelse(use_bootstrap, "TRUE", "FALSE"), B))

# -----------------------------------------------------------------------------
# 6. Utilidades
# -----------------------------------------------------------------------------
clamp <- function(v, n) pmin(pmax(v, 1), n)

get_window <- function(mat, i, j, N) {
  r0 <- (i - N):(i + N)
  c0 <- (j - N):(j + N)
  r  <- clamp(r0, nrow(mat))
  c  <- clamp(c0, ncol(mat))
  mat[r, c, drop = FALSE]
}

border_stats <- function(win) {
  Lw <- nrow(win)
  is_edge <- (row(win) == 1) | (row(win) == Lw) | (col(win) == 1) | (col(win) == Lw)
  vals <- win[is_edge]
  mu_b <- mean(vals)
  var_b_raw <- mean(vals^2) - mu_b^2
  var_b <- if (var_b_raw < 0) 0 else var_b_raw
  list(mu = mu_b, var = var_b)
}

extract_interior <- function(mat, margin) {
  if (margin <= 0) return(as.vector(mat))
  nr <- nrow(mat)
  nc <- ncol(mat)
  if ((2 * margin + 1) >= nr || (2 * margin + 1) >= nc) {
    stop("El margen es demasiado grande para la matriz.")
  }
  interior <- mat[(margin + 1):(nr - margin), (margin + 1):(nc - margin), drop = FALSE]
  as.vector(interior)
}

extract_interior_matrix <- function(mat, margin) {
  if (margin <= 0) return(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  if ((2 * margin + 1) >= nr || (2 * margin + 1) >= nc) {
    stop("El margen es demasiado grande para la matriz.")
  }
  mat[(margin + 1):(nr - margin), (margin + 1):(nc - margin), drop = FALSE]
}

summary_stats <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  c(
    n       = length(x),
    mean    = mean(x),
    sd      = sd(x),
    median  = median(x),
    q01     = as.numeric(quantile(x, 0.01)),
    q05     = as.numeric(quantile(x, 0.05)),
    q95     = as.numeric(quantile(x, 0.95)),
    q99     = as.numeric(quantile(x, 0.99))
  )
}

# -----------------------------------------------------------------------------
# 7. Shannon helpers
# -----------------------------------------------------------------------------
shannon_theoretical <- function(L, mu) {
  log(mu) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L))
}

bootstrap_al_omari_1_estimator_local <- function(x, B) {
  x <- x[is.finite(x) & x > 0]
  
  n <- length(x)
  if (n < 2L) return(NA_real_)
  
  v_bootstrap <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    same_values <- TRUE
    tries <- 0L
    
    while (same_values) {
      tries <- tries + 1L
      if (tries > 100L) break
      
      bx <- sample(x, size = n, replace = TRUE)
      
      if (!all(bx == bx[1])) {
        same_values <- FALSE
        
        entropy_result <- tryCatch(
          al_omari_1_estimator(bx),
          error = function(e) NA_real_
        )
        
        if (is.finite(entropy_result)) {
          v_bootstrap[b] <- entropy_result
        } else {
          v_bootstrap[b] <- NA_real_
        }
      }
    }
  }
  
  t0 <- tryCatch(al_omari_1_estimator(x), error = function(e) NA_real_)
  if (!is.finite(t0)) return(NA_real_)
  
  valid_vals <- v_bootstrap[is.finite(v_bootstrap)]
  if (length(valid_vals) < B / 2) return(NA_real_)
  
  estimated_mean <- mean(valid_vals)
  2 * t0 - estimated_mean
}

stat_shannon_window <- function(z_win, L, use_bootstrap, B) {
  z <- as.vector(z_win)
  z <- z[is.finite(z) & z > 0]
  
  if (length(z) < 2L) return(NA_real_)
  
  mu_hat <- mean(z)
  
  est <- if (use_bootstrap) {
    bootstrap_al_omari_1_estimator_local(z, B = B)
  } else {
    tryCatch(al_omari_1_estimator(z), error = function(e) NA_real_)
  }
  
  if (!is.finite(est)) return(NA_real_)
  
  theo <- shannon_theoretical(L = L, mu = mu_hat)
  
  est - theo
}

# -----------------------------------------------------------------------------
# 8. Runner FIXED
# -----------------------------------------------------------------------------
run_fixed_shannon <- function(mat, window_size, L, use_bootstrap, B, parallel = TRUE) {
  if (window_size %% 2 == 0 || window_size < 3) {
    stop("window_size debe ser impar y >= 3.")
  }
  
  rows <- nrow(mat)
  cols <- ncol(mat)
  out_rows <- rows - window_size + 1
  out_cols <- cols - window_size + 1
  
  idx <- expand.grid(i = seq_len(out_rows), j = seq_len(out_cols))
  
  compute_one <- function(k) {
    i <- idx$i[k]
    j <- idx$j[k]
    win <- mat[i:(i + window_size - 1), j:(j + window_size - 1), drop = FALSE]
    stat_shannon_window(win, L, use_bootstrap, B)
  }
  
  start_time <- Sys.time()
  
  res <- if (parallel) {
    future_sapply(seq_len(nrow(idx)), compute_one, future.seed = TRUE)
  } else {
    sapply(seq_len(nrow(idx)), compute_one)
  }
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  diff_map <- matrix(res, nrow = out_rows, byrow = FALSE)
  W_map <- matrix(window_size, nrow = out_rows, ncol = out_cols)
  
  list(
    diff_map      = diff_map,
    W_map         = W_map,
    window_mode   = "fixed",
    window_size   = window_size,
    crop_margin   = (window_size - 1) / 2,
    L             = L,
    bootstrap     = use_bootstrap,
    B             = if (use_bootstrap) B else 0L,
    exec_seconds  = elapsed
  )
}

# -----------------------------------------------------------------------------
# 9. Runner ADAPTIVE
# -----------------------------------------------------------------------------
run_adaptive_shannon <- function(mat, W_min, W_max, eta, sigma_n, L,
                                 use_bootstrap, B, parallel = TRUE) {
  if (W_min %% 2 == 0 || W_max %% 2 == 0) stop("W_min y W_max deben ser impares.")
  if (W_min < 3) stop("W_min debe ser >= 3.")
  if (W_min > W_max) stop("W_min no puede ser > W_max.")
  
  rows <- nrow(mat)
  cols <- ncol(mat)
  
  adaptive_row <- function(i) {
    N_prev   <- (W_min - 1) / 2
    row_stat <- numeric(cols)
    row_W    <- integer(cols)
    
    eps_mu <- .Machine$double.eps
    
    for (j in seq_len(cols)) {
      win_big <- get_window(mat, i, j, N_prev)
      bs      <- border_stats(win_big)
      mu_b    <- bs$mu
      sd_b    <- sqrt(bs$var)
      Cij     <- if (mu_b > eps_mu) sd_b / mu_b else Inf
      
      W_prev <- 2 * N_prev + 1
      Tij <- eta * (1 + sqrt((1 + 2 * sigma_n^2) / (8 * (W_prev - 1)))) * sigma_n
      
      if (Cij <= Tij) {
        N_curr <- min(N_prev + 1, (W_max - 1) / 2)
      } else {
        N_curr <- max(N_prev - 1, (W_min - 1) / 2)
      }
      
      W_curr <- 2 * N_curr + 1
      row_W[j] <- W_curr
      
      win_full <- get_window(mat, i, j, N_curr)
      row_stat[j] <- stat_shannon_window(win_full, L, use_bootstrap, B)
      
      N_prev <- N_curr
    }
    
    list(stat = row_stat, Ws = row_W)
  }
  
  start_time <- Sys.time()
  
  res_list <- if (parallel) {
    future_lapply(seq_len(rows), adaptive_row, future.seed = TRUE)
  } else {
    lapply(seq_len(rows), adaptive_row)
  }
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  diff_map <- t(sapply(res_list, `[[`, "stat"))
  W_map    <- t(sapply(res_list, `[[`, "Ws"))
  
  list(
    diff_map      = diff_map,
    W_map         = W_map,
    window_mode   = "adaptive",
    W_min         = W_min,
    W_max         = W_max,
    crop_margin   = (W_max - 1) / 2,
    eta           = eta,
    sigma_n       = sigma_n,
    L             = L,
    bootstrap     = use_bootstrap,
    B             = if (use_bootstrap) B else 0L,
    exec_seconds  = elapsed
  )
}

# -----------------------------------------------------------------------------
# 10. Ejecutar FIXED
# -----------------------------------------------------------------------------
fixed_results <- list()
fixed_summary_rows <- list()

for (w in fixed_windows) {
  cat(sprintf("\n================ FIXED W = %d ================\n", w))
  res_fix <- run_fixed_shannon(
    mat = x,
    window_size = w,
    L = L,
    use_bootstrap = use_bootstrap,
    B = B,
    parallel = parallel
  )
  
  fixed_results[[paste0("W", w)]] <- res_fix
  
  vals <- as.vector(res_fix$diff_map)
  stats <- summary_stats(vals)
  
  fixed_summary_rows[[paste0("W", w)]] <- data.frame(
    method       = "fixed",
    W            = w,
    n            = stats["n"],
    mu0_hat      = stats["mean"],
    sigma0_hat   = stats["sd"],
    median       = stats["median"],
    q01          = stats["q01"],
    q05          = stats["q05"],
    q95          = stats["q95"],
    q99          = stats["q99"],
    exec_seconds = res_fix$exec_seconds,
    stringsAsFactors = FALSE
  )
  
  save(
    res_fix,
    file = file.path(outdir, sprintf("shannon_null_fixed_w%d_L%d%s.RData",
                                     w, L,
                                     if (use_bootstrap) sprintf("_B%d", B) else ""))
  )
}

fixed_summary <- do.call(rbind, fixed_summary_rows)

# -----------------------------------------------------------------------------
# 11. Ejecutar ADAPTIVE
# -----------------------------------------------------------------------------
cat("\n================ ADAPTIVE ================\n")

adaptive_result <- run_adaptive_shannon(
  mat = x,
  W_min = W_min,
  W_max = W_max,
  eta = eta,
  sigma_n = sigma_n,
  L = L,
  use_bootstrap = use_bootstrap,
  B = B,
  parallel = parallel
)

save(
  adaptive_result,
  file = file.path(outdir, sprintf("shannon_null_adaptive_w%d-%d_eta_%s_L%d%s.RData",
                                   W_min, W_max, gsub("\\.", "", format(eta)),
                                   L,
                                   if (use_bootstrap) sprintf("_B%d", B) else ""))
)

# -----------------------------------------------------------------------------
# 12. Resumen ADAPTIVE
# -----------------------------------------------------------------------------
adaptive_stats <- summary_stats(as.vector(adaptive_result$diff_map))

adaptive_summary <- data.frame(
  method       = "adaptive",
  W            = NA_integer_,
  n            = adaptive_stats["n"],
  mu0_hat      = adaptive_stats["mean"],
  sigma0_hat   = adaptive_stats["sd"],
  median       = adaptive_stats["median"],
  q01          = adaptive_stats["q01"],
  q05          = adaptive_stats["q05"],
  q95          = adaptive_stats["q95"],
  q99          = adaptive_stats["q99"],
  exec_seconds = adaptive_result$exec_seconds,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 13. Frecuencias de W_map
# -----------------------------------------------------------------------------
margin_ad <- adaptive_result$crop_margin

W_map_all <- adaptive_result$W_map
W_map_interior <- extract_interior_matrix(W_map_all, margin_ad)

freq_all <- as.data.frame(table(W_map_all))
names(freq_all) <- c("W", "count")
freq_all$W <- as.integer(as.character(freq_all$W))
freq_all$prop <- freq_all$count / sum(freq_all$count)
freq_all$region <- "all"

freq_interior <- as.data.frame(table(W_map_interior))
names(freq_interior) <- c("W", "count")
freq_interior$W <- as.integer(as.character(freq_interior$W))
freq_interior$prop <- freq_interior$count / sum(freq_interior$count)
freq_interior$region <- "interior"

W_freq <- rbind(freq_all, freq_interior)

# -----------------------------------------------------------------------------
# 14. Resumen adaptive condicionado por W
# -----------------------------------------------------------------------------
adaptive_diff_interior <- extract_interior_matrix(adaptive_result$diff_map, margin_ad)
adaptive_W_interior    <- extract_interior_matrix(adaptive_result$W_map, margin_ad)

cond_rows <- list()
for (w in sort(unique(as.vector(adaptive_W_interior)))) {
  vals_w <- adaptive_diff_interior[adaptive_W_interior == w]
  st_w <- summary_stats(vals_w)
  
  cond_rows[[paste0("W", w)]] <- data.frame(
    method       = "adaptive_conditional",
    W            = as.integer(w),
    n            = st_w["n"],
    mu0_hat      = st_w["mean"],
    sigma0_hat   = st_w["sd"],
    median       = st_w["median"],
    q01          = st_w["q01"],
    q05          = st_w["q05"],
    q95          = st_w["q95"],
    q99          = st_w["q99"],
    exec_seconds = NA_real_,
    stringsAsFactors = FALSE
  )
}
adaptive_conditional_summary <- do.call(rbind, cond_rows)

# -----------------------------------------------------------------------------
# 15. Tablas
# -----------------------------------------------------------------------------
null_calibration_table <- rbind(
  fixed_summary[, c("method", "W", "n", "mu0_hat", "sigma0_hat", "median", "q01", "q05", "q95", "q99")],
  adaptive_summary[, c("method", "W", "n", "mu0_hat", "sigma0_hat", "median", "q01", "q05", "q95", "q99")],
  adaptive_conditional_summary[, c("method", "W", "n", "mu0_hat", "sigma0_hat", "median", "q01", "q05", "q95", "q99")]
)

write.csv(fixed_summary,
          file = file.path(outdir, "shannon_null_fixed_summary.csv"),
          row.names = FALSE)

write.csv(adaptive_summary,
          file = file.path(outdir, "shannon_null_adaptive_summary.csv"),
          row.names = FALSE)

write.csv(adaptive_conditional_summary,
          file = file.path(outdir, "shannon_null_adaptive_conditional_summary.csv"),
          row.names = FALSE)

write.csv(W_freq,
          file = file.path(outdir, "shannon_null_adaptive_Wmap_frequencies.csv"),
          row.names = FALSE)

write.csv(null_calibration_table,
          file = file.path(outdir, "shannon_null_calibration_master_table.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 16. Datos para gráficas
# -----------------------------------------------------------------------------
density_df <- data.frame(
  value = c(
    extract_interior(fixed_results$W5$diff_map, 0),
    extract_interior(fixed_results$W7$diff_map, 0),
    extract_interior(fixed_results$W9$diff_map, 0),
    extract_interior(fixed_results$W11$diff_map, 0),
    extract_interior(adaptive_result$diff_map, margin_ad)
  ),
  method = c(
    rep("Fixed W=5",  length(extract_interior(fixed_results$W5$diff_map, 0))),
    rep("Fixed W=7",  length(extract_interior(fixed_results$W7$diff_map, 0))),
    rep("Fixed W=9",  length(extract_interior(fixed_results$W9$diff_map, 0))),
    rep("Fixed W=11", length(extract_interior(fixed_results$W11$diff_map, 0))),
    rep("Adaptive (interior)", length(extract_interior(adaptive_result$diff_map, margin_ad)))
  ),
  stringsAsFactors = FALSE
)

density_cond_df <- data.frame(
  value = as.vector(adaptive_diff_interior),
  W     = factor(as.vector(adaptive_W_interior),
                 levels = sort(unique(as.vector(adaptive_W_interior))))
)

# -----------------------------------------------------------------------------
# 17. Gráficas
# -----------------------------------------------------------------------------
p1 <- ggplot(density_df, aes(x = value, color = method)) +
  geom_density(linewidth = 1) +
  labs(
    title = "Distribuciones del estadístico bajo H0",
    subtitle = "Comparación entre ventanas fijas y adaptive (interior)",
    x = "diff_map",
    y = "Densidad",
    color = "Método"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "density_fixed_vs_adaptive.png"),
  plot = p1, width = 10, height = 6, dpi = 300
)

p2 <- ggplot(W_freq, aes(x = factor(W), y = prop, fill = region)) +
  geom_col(position = "dodge") +
  labs(
    title = "Frecuencia de tamaños de ventana seleccionados por adaptive",
    subtitle = "Comparación entre toda la imagen e interior",
    x = "Tamaño de ventana W",
    y = "Proporción",
    fill = "Zona"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "adaptive_Wmap_frequencies.png"),
  plot = p2, width = 9, height = 5.5, dpi = 300
)

p3 <- ggplot(density_cond_df, aes(x = value, color = W)) +
  geom_density(linewidth = 1) +
  labs(
    title = "Distribución adaptive condicionada por W seleccionado",
    subtitle = "Solo región interior",
    x = "diff_map",
    y = "Densidad",
    color = "W"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "adaptive_conditional_densities.png"),
  plot = p3, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 18. Guardar objeto maestro
# -----------------------------------------------------------------------------
master_result <- list(
  image = x,
  L = L,
  use_bootstrap = use_bootstrap,
  B = B,
  data_scale = data_scale,
  sigma_n = sigma_n,
  fixed_windows = fixed_windows,
  W_min = W_min,
  W_max = W_max,
  eta = eta,
  fixed_results = fixed_results,
  adaptive_result = adaptive_result,
  fixed_summary = fixed_summary,
  adaptive_summary = adaptive_summary,
  adaptive_conditional_summary = adaptive_conditional_summary,
  W_freq = W_freq,
  null_calibration_table = null_calibration_table
)

save(
  master_result,
  file = file.path(outdir, "shannon_null_diagnostic_master.RData")
)

# -----------------------------------------------------------------------------
# 19. Consola
# -----------------------------------------------------------------------------
cat("\n================ RESUMEN =================\n")
cat("Tabla de calibración fija:\n")
print(fixed_summary)

cat("\nResumen adaptive global:\n")
print(adaptive_summary)

cat("\nFrecuencia W_map (toda la imagen e interior):\n")
print(W_freq)

cat("\nResumen adaptive condicionado por W (interior):\n")
print(adaptive_conditional_summary)

cat("\nArchivos guardados en:\n")
cat(normalizePath(outdir), "\n")
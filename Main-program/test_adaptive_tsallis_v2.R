# =============================================================================
# test_adaptive_tsallis_FIXED_ETA.R
#
# Tsallis + adaptive windows (Park) with FIXED eta (no auto-calibration).
# Adjust sigma_n according to the data scale:
#   - amplitude   -> sigma_n = 0.523 / sqrt(L)
#   - intensity   -> sigma_n = 1 / sqrt(L)
#
# Outputs:
#   - difference_values: matrix of the statistic (estimated - theoretical)
#   - L_map: map of the selected local window size (odd side length)
#   - .Rdata with everything needed for later plotting/analysis
# =============================================================================

# --- Libraries and parallel plan ----------------------------------------------
req_pkgs <- c("future.apply", "future")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(future.apply)
  library(future)
})

# Portable execution plan (multisession if >1 core)
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
    plan(sequential); parallel <- FALSE
  }
}, silent = TRUE)

set.seed(1234567890, kind = "Mersenne-Twister")
message(sprintf("Plan: %s con %d worker(s).",
                if (parallel) "multisession" else "sequential",
                future::nbrOfWorkers()))

# --- Parameters ----------------------------------------------------------------
input_type <- "envi"   # "envi" or "sim"

# ENVI (if input_type == "envi")
img_path   <- "../../../Data/SAR/L12_envi_munich_size_1024/Intensity_HH.img"
hdr_path   <- "../../../Data/SAR/L12_envi_munich_size_1024/Intensity_HH.hdr"

# Simulated (if input_type == "sim")
sim_rdata <- "./Data/L9_simulated_image_500.Rdata"
sim_var   <- "Z"   # object name inside the RData (if known)

# --- Data scale and SAR/entropy parameters ------------------------------------
# IMPORTANT: choose the correct scale:
#   "intensity"  -> if your image is intensity/power (HH "Intensity")
#   "amplitude"  -> if your image is amplitude (magnitude)
data_scale <- "intensity"   # "intensity" or "amplitude"

L       <- 12          # looks
lambda  <- 0.85
B       <- 100         # bootstraps for Tsallis
W_min   <- 5           # minimum odd side
W_max   <- 11          # maximum odd side
eta     <- 3          # <-- FIXED ETA 

# sigma_n according to data scale
sigma_n <- if (identical(tolower(data_scale), "amplitude")) {
  0.523 / sqrt(L)
} else if (identical(tolower(data_scale), "intensity")) {
  1 / sqrt(L)
} else {
  stop("data_scale debe ser 'intensity' o 'amplitude'.")
}

cat(sprintf("Escala de datos: %s  |  L=%d  |  sigma_n=%.6f  |  eta fijo=%.3f\n",
            data_scale, L, sigma_n, eta))

# --- Window checks -------------------------------------------------------------
if (W_min %% 2 == 0 || W_max %% 2 == 0) stop("W_min y W_max deben ser impares.")
if (W_min < 3) stop("W_min debe ser >= 3.")
if (W_min > W_max) stop("W_min no puede ser > W_max.")

# --- External functions --------------------------------------------------------
source("./Code/read_ENVI_images.R")
source("./Code/gamma_sar_sample.R")
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")

# --- Load image (ENVI or simulated) -------------------------------------------
load_simulated_matrix <- function(rdata_path, var_name = NULL) {
  stopifnot(file.exists(rdata_path))
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  
  pick_object <- function(obj) {
    if (is.data.frame(obj)) obj <- as.matrix(obj)
    if (is.array(obj) && length(dim(obj)) > 2) obj <- obj[,,1, drop = FALSE]
    if (!is.matrix(obj)) obj <- tryCatch(as.matrix(obj), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    storage.mode(obj) <- "double"
    obj
  }
  
  if (!is.null(var_name)) {
    stopifnot(is.character(var_name), length(var_name) == 1L)
    if (!exists(var_name, envir = env, inherits = FALSE)) {
      stop(sprintf("Objeto '%s' no encontrado en '%s'.", var_name, rdata_path))
    }
    obj <- get(var_name, envir = env, inherits = FALSE)
    obj <- pick_object(obj)
    if (is.null(obj)) stop(sprintf("Objeto '%s' no es una matriz 2D numérica.", var_name))
    return(obj)
  }
  
  # Auto-detect the first usable 2D matrix
  nms <- ls(envir = env, all.names = TRUE)
  for (nm in nms) {
    cand <- get(nm, envir = env, inherits = FALSE)
    obj  <- pick_object(cand)
    if (!is.null(obj) && is.matrix(obj) && all(dim(obj) >= 3)) {
      return(obj)
    }
  }
  stop(sprintf("No se encontró una matriz 2D numérica en '%s'. Ajusta 'sim_var'.", rdata_path))
}

if (identical(tolower(input_type), "envi")) {
  stopifnot(file.exists(img_path), file.exists(hdr_path))
  x <- myread.ENVI(img_path, hdr_path)
  input_label <- "envi"
} else if (identical(tolower(input_type), "sim")) {
  x <- load_simulated_matrix(sim_rdata, sim_var)
  input_label <- "sim"
} else {
  stop("input_type debe ser 'envi' o 'sim'.")
}

# --- Shape checks --------------------------------------------------------------
if (!is.matrix(x) || any(dim(x) < 3)) {
  stop("La imagen 'x' debe ser matriz 2D numérica con al menos 3x3 píxeles.")
}
rows <- nrow(x); cols <- ncol(x)
cat(sprintf("Imagen '%s' cargada: %d x %d\n", input_label, rows, cols))

# --- Utilities -----------------------------------------------------------------
clamp <- function(v, n) pmin(pmax(v, 1), n)

get_window <- function(i, j, N) {
  r0 <- (i - N):(i + N); c0 <- (j - N):(j + N)
  r  <- clamp(r0, rows);  c  <- clamp(c0, cols)
  x[r, c, drop = FALSE]
}

border_stats <- function(win) {
  Lw <- nrow(win)
  is_edge <- (row(win) == 1) | (row(win) == Lw) | (col(win) == 1) | (col(win) == Lw)
  vals <- win[is_edge]
  # minimal robustness
  mu_b <- mean(vals)
  var_b_raw <- mean(vals^2) - mu_b^2
  var_b <- if (var_b_raw < 0) 0 else var_b_raw
  list(mu = mu_b, var = var_b)
}

# Tsallis statistic (bootstrap) minus theoretical value (Gamma-Amplitude/Intensity)
stat_tsallis_window <- function(mat_win) {
  z <- as.vector(mat_win)
  mu_hat <- mean(z)
  est <- bootstrap_tsallis_entropy_optimized(
    z, B = B, lambda = lambda, parallel = FALSE
  )
  theo <- (1 - exp((1 - lambda) * log(mu_hat) +
                     (lambda - 1) * log(L) +
                     lgamma(lambda * (L - 1) + 1) -
                     lambda * lgamma(L) -
                     (lambda * (L - 1) + 1) * log(lambda))) / (lambda - 1)
  est - theo
}

# --- Row-wise adaptive sweep (Park) --------------------------------------------
# Note: following Park, we use N_{i,j+1} starting from N_{i,j}.
# Initialize each row with N_prev = (W_min - 1)/2 (as in prior scripts).
adaptive_row <- function(i) {
  N_prev   <- (W_min - 1) / 2
  row_stat <- numeric(cols)
  row_L    <- integer(cols)
  
  eps_mu <- .Machine$double.eps
  
  for (j in seq_len(cols)) {
    # Threshold is evaluated with the CURRENT size (L_prev)
    win_big <- get_window(i, j, N_prev)
    bs      <- border_stats(win_big)
    mu_b    <- bs$mu
    sd_b    <- sqrt(bs$var)
    Cij     <- if (mu_b > eps_mu) sd_b / mu_b else Inf
    
    L_prev  <- 2 * N_prev + 1
    Tij <- eta * (1 + sqrt((1 + 2 * sigma_n^2) / (8 * (L_prev - 1)))) * sigma_n
    
    if (Cij <= Tij) {
      N_curr <- min(N_prev + 1, (W_max - 1) / 2)
    } else {
      N_curr <- max(N_prev - 1, (W_min - 1) / 2)
    }
    L_curr      <- 2 * N_curr + 1
    row_L[j]    <- L_curr
    
    win_full    <- get_window(i, j, N_curr)
    row_stat[j] <- stat_tsallis_window(win_full)
    
    N_prev <- N_curr
  }
  list(stat = row_stat, Ls = row_L)
}

# --- Execution -----------------------------------------------------------------
start_time <- Sys.time()
res_list <- if (parallel) {
  future_lapply(seq_len(rows), adaptive_row, future.seed = TRUE)
} else {
  lapply(seq_len(rows), adaptive_row)
}
end_time <- Sys.time()

# --- Reconstruct outputs -------------------------------------------------------
difference_values <- t(sapply(res_list, `[[`, "stat"))
L_map             <- t(sapply(res_list, `[[`, "Ls"))

# --- Save ----------------------------------------------------------------------
outdir  <- "./Data"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile <- sprintf(
  "%s/%s_munich1024adaptive_tsallis_FIXED_ETA_%s_L%d_eta_%.2f_%dx%d_B%d_lambda_%s.Rdata",
  outdir,
  if (identical(input_label, "envi")) "envi" else "sim",
  tolower(data_scale), L, eta, W_min, W_max, B, gsub("\\.", "", format(lambda))
)
save(difference_values, x, L_map, data_scale, L, eta, sigma_n, W_min, W_max, B, lambda,
     file = outfile)

# --- Timing and summary --------------------------------------------------------
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat(sprintf("Finished in %.1f sec\n", elapsed))
cat(sprintf("Results saved to '%s'\n", outfile))
h <- floor(elapsed / 3600); m <- floor((elapsed %% 3600) / 60); s <- round(elapsed %% 60, 2)
cat(sprintf("Total time: %02d:%02d:%05.2f (hh:mm:ss)\n", h, m, s))

cat("Rango de W_ij: ", paste(range(L_map), collapse = " "), "\n")
cat("Frecuencia de tamaños de ventana:\n"); print(table(L_map))

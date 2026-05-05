# =============================================================================
# run_shannon_adaptive_dublin.R
#
# Shannon statistic with ADAPTIVE windows for Dublin SAR image.
# Saves:
#   - result$diff_map
#   - result$W_map
#   - metadata
#   - exec_time in mm:ss
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
req_pkgs <- c("future.apply", "future")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(future.apply)
  library(future)
})

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

# -----------------------------------------------------------------------------
# 1. External functions
# -----------------------------------------------------------------------------
source("./Code/read_ENVI_images.R")
source("./Code/gamma_sar_sample.R")
source("./Code/al_omari_1_estimator.R")

# -----------------------------------------------------------------------------
# 2. Parameters
# -----------------------------------------------------------------------------
input_type <- "envi"

scene_name <- "dublin"
image_id   <- "L16_dublin_1600_nuevo"

img_path <- "./Data-SAR/L16_dublin_1600_nuevo/Intensity_HH.img"
hdr_path <- "./Data-SAR/L16_dublin_1600_nuevo/Intensity_HH.hdr"

data_scale    <- "intensity"
L             <- 16
use_bootstrap <- TRUE
B             <- 50

W_min <- 5
W_max <- 11
eta   <- 2

outdir <- "./Data/results/shannon1"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 3. Helper: format execution time as mm:ss
# -----------------------------------------------------------------------------
format_elapsed_mmss <- function(elapsed_sec) {
  total_sec <- round(elapsed_sec)
  mm <- total_sec %/% 60
  ss <- total_sec %% 60
  sprintf("%02d:%02d", mm, ss)
}

# -----------------------------------------------------------------------------
# 4. Derived parameters
# -----------------------------------------------------------------------------
sigma_n <- if (tolower(data_scale) == "amplitude") {
  0.523 / sqrt(L)
} else if (tolower(data_scale) == "intensity") {
  1 / sqrt(L)
} else {
  stop("data_scale must be 'intensity' or 'amplitude'.")
}

if (W_min %% 2 == 0 || W_max %% 2 == 0) stop("W_min and W_max must be odd.")
if (W_min < 3) stop("W_min must be >= 3.")
if (W_min > W_max) stop("W_min cannot be > W_max.")

# -----------------------------------------------------------------------------
# 5. Shannon helpers
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

# -----------------------------------------------------------------------------
# 6. Load image
# -----------------------------------------------------------------------------
stopifnot(file.exists(img_path), file.exists(hdr_path))
x <- myread.ENVI(img_path, hdr_path)

if (!is.matrix(x) || any(dim(x) < 3)) {
  stop("Input image must be a numeric 2D matrix with at least 3x3 pixels.")
}

rows <- nrow(x)
cols <- ncol(x)

cat(sprintf("Image loaded: %s | %d x %d\n", image_id, rows, cols))
cat(sprintf("Mode: ADAPTIVE | W_min = %d | W_max = %d | eta = %d | bootstrap = %s\n",
            W_min, W_max, eta, ifelse(use_bootstrap, "TRUE", "FALSE")))

# -----------------------------------------------------------------------------
# 7. Utilities
# -----------------------------------------------------------------------------
clamp <- function(v, n) pmin(pmax(v, 1), n)

get_window <- function(i, j, N) {
  r0 <- (i - N):(i + N)
  c0 <- (j - N):(j + N)
  r  <- clamp(r0, rows)
  c  <- clamp(c0, cols)
  x[r, c, drop = FALSE]
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

stat_shannon_window <- function(mat_win) {
  z <- as.vector(mat_win)
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
# 8. Row-wise adaptive sweep
# -----------------------------------------------------------------------------
adaptive_row <- function(i) {
  N_prev   <- (W_min - 1) / 2
  row_stat <- numeric(cols)
  row_W    <- integer(cols)
  
  eps_mu <- .Machine$double.eps
  
  for (j in seq_len(cols)) {
    win_big <- get_window(i, j, N_prev)
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
    
    W_curr      <- 2 * N_curr + 1
    row_W[j]    <- W_curr
    win_full    <- get_window(i, j, N_curr)
    row_stat[j] <- stat_shannon_window(win_full)
    
    N_prev <- N_curr
  }
  
  list(stat = row_stat, Ws = row_W)
}

# -----------------------------------------------------------------------------
# 9. Run
# -----------------------------------------------------------------------------
start_time <- Sys.time()

res_list <- if (parallel) {
  future_lapply(seq_len(rows), adaptive_row, future.seed = TRUE)
} else {
  lapply(seq_len(rows), adaptive_row)
}

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
elapsed_mmss <- format_elapsed_mmss(elapsed)

diff_map <- t(sapply(res_list, `[[`, "stat"))
W_map    <- t(sapply(res_list, `[[`, "Ws"))

# -----------------------------------------------------------------------------
# 10. Save
# -----------------------------------------------------------------------------
result <- list(
  diff_map      = diff_map,
  image         = x,
  W_map         = W_map,
  entropy       = "shannon",
  bootstrap     = use_bootstrap,
  B             = if (use_bootstrap) B else 0L,
  window_mode   = "adaptive",
  W_min         = W_min,
  W_max         = W_max,
  crop_margin   = (W_max - 1) / 2,
  eta           = eta,
  sigma_n       = sigma_n,
  L             = L,
  data_scale    = data_scale,
  scene         = scene_name,
  image_id      = image_id,
  img_path      = img_path,
  hdr_path      = hdr_path,
  exec_time     = elapsed_mmss
)

outfile <- sprintf(
  "%s/%s_shannon_adaptive_%s_eta%d_w%d-%d_L%d%s.RData",
  outdir,
  scene_name,
  ifelse(use_bootstrap, "boot", "noboot"),
  eta,
  W_min, W_max,
  L,
  if (use_bootstrap) sprintf("_B%d", B) else ""
)

save(result, file = outfile)

cat(sprintf("Finished in %s (mm:ss)\n", elapsed_mmss))
cat(sprintf("Saved to: %s\n", outfile))
cat("W_map frequencies:\n")
print(table(W_map, useNA = "ifany"))
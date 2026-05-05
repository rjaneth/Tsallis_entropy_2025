# =============================================================================
# run_shannon_fixed_dublin.R
#
# Shannon statistic with FIXED sliding window for Dublin SAR image.
# Saves:
#   - result$diff_map
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
scene_name   <- "dublin"
image_id     <- "L16_dublin_1600_nuevo"

img_path     <- "./Data-SAR/L16_dublin_1600_nuevo/Intensity_HH.img"
hdr_path     <- "./Data-SAR/L16_dublin_1600_nuevo/Intensity_HH.hdr"

data_scale   <- "intensity"
L            <- 16
window_size  <- 7

use_bootstrap <- TRUE
B             <- 80

outdir <- "./Data/results/shannonf_80"
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
# 4. Shannon helpers
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
# 5. Load image
# -----------------------------------------------------------------------------
stopifnot(file.exists(img_path), file.exists(hdr_path))
x <- myread.ENVI(file = img_path, headerfile = hdr_path)

if (!is.matrix(x) || any(dim(x) < 3)) {
  stop("Input image must be a numeric 2D matrix with at least 3x3 pixels.")
}

rows <- nrow(x)
cols <- ncol(x)

if (window_size %% 2 == 0 || window_size < 3) {
  stop("window_size must be odd and >= 3.")
}

out_rows <- rows - window_size + 1
out_cols <- cols - window_size + 1
diff_map <- matrix(NA_real_, nrow = out_rows, ncol = out_cols)

cat(sprintf("Image loaded: %s | %d x %d\n", image_id, rows, cols))
cat(sprintf("Mode: FIXED | window_size = %d | bootstrap = %s\n",
            window_size, ifelse(use_bootstrap, "TRUE", "FALSE")))

# -----------------------------------------------------------------------------
# 6. Statistic function
# -----------------------------------------------------------------------------
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
# 7. Sliding-window computation
# -----------------------------------------------------------------------------
start_time <- Sys.time()

idx <- expand.grid(i = seq_len(out_rows), j = seq_len(out_cols))

compute_one <- function(k) {
  i <- idx$i[k]
  j <- idx$j[k]
  win <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
  stat_shannon_window(win)
}

res <- if (parallel) {
  future_sapply(seq_len(nrow(idx)), compute_one, future.seed = TRUE)
} else {
  sapply(seq_len(nrow(idx)), compute_one)
}

diff_map[] <- matrix(res, nrow = out_rows, byrow = FALSE)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
elapsed_mmss <- format_elapsed_mmss(elapsed)

# -----------------------------------------------------------------------------
# 8. Save
# -----------------------------------------------------------------------------
result <- list(
  diff_map      = diff_map,
  image         = x,
  entropy       = "shannon",
  bootstrap     = use_bootstrap,
  B             = if (use_bootstrap) B else 0L,
  window_mode   = "fixed",
  window_size   = window_size,
  crop_margin   = (window_size - 1) / 2,
  L             = L,
  data_scale    = data_scale,
  scene         = scene_name,
  image_id      = image_id,
  img_path      = img_path,
  hdr_path      = hdr_path,
  exec_time     = elapsed_mmss
)

outfile <- sprintf(
  "%s/%s_shannon_fixed_%s_w%d_L%d%s.RData",
  outdir,
  scene_name,
  ifelse(use_bootstrap, "boot", "noboot"),
  window_size,
  L,
  if (use_bootstrap) sprintf("_B%d", B) else ""
)

save(result, file = outfile)

cat(sprintf("Finished in %s (mm:ss)\n", elapsed_mmss))
cat(sprintf("Saved to: %s\n", outfile))
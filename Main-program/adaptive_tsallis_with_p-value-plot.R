# =============================================================================
# adaptive_tsallis_with_Lmap.R
# Adaptive Tsallis with window-size adaptation — row-parallel + L_map
# Supports ENVI (.img + .hdr) and simulated RData inputs
# Now includes: p-values + plots (imagematrix_visualizer if available)
# =============================================================================

# --- Libraries and parallel plan (portable and auto-adjusting) ----------------
req_pkgs <- c("future.apply", "future")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(future.apply)
  library(future)
})

# --- Worker selection (priority: options -> env -> auto max(1, cores-1)) ------
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
    plan(multisession, workers = workers)   # add gc=TRUE if you want more frequent GC
  } else {
    plan(sequential)
    parallel <- FALSE
  }
}, silent = TRUE)

set.seed(1234567890, kind = "Mersenne-Twister")
message(sprintf("Plan: %s with %d worker(s).",
                if (parallel) "multisession" else "sequential",
                future::nbrOfWorkers()))

# --- Path bootstrap (sibling folders: Code/, Data/, Data-SAR/) ----------------
get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (rstudioapi::isAvailable()) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
      if (nzchar(p)) return(dirname(normalizePath(p)))
    }
  }
  if (requireNamespace("knitr", quietly = TRUE)) {
    p <- tryCatch(knitr::current_input(), error = function(e) NULL)
    if (!is.null(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}

SCRIPT_DIR   <- get_script_dir()
PROJECT_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = FALSE)
CODE_DIR     <- file.path(PROJECT_DIR, "./Code")
DATA_DIR     <- file.path(PROJECT_DIR, "./Data")
DATA_SAR_DIR <- file.path(PROJECT_DIR, "./Data-SAR")

cat("SCRIPT_DIR  :", SCRIPT_DIR,   "\n")
cat("PROJECT_DIR :", PROJECT_DIR,  "\n")
cat("CODE_DIR    :", CODE_DIR,     "\n")
cat("DATA_DIR    :", DATA_DIR,     "\n")
cat("DATA_SAR_DIR:", DATA_SAR_DIR, "\n")

# --- External helper functions ------------------------------------------------
source(file.path(CODE_DIR, "read_ENVI_images.R"))
source(file.path(CODE_DIR, "gamma_sar_sample.R"))
source(file.path(CODE_DIR, "tsallis_estimator_optimized.R"))
source(file.path(CODE_DIR, "bootstrap_tsallis_entropy_optimized.R"))

# imagematrix_visualizer.R is optional; if missing we will fall back to base plots
viz_available <- file.exists(file.path(CODE_DIR, "imagematrix_visualizer_v2.R"))
if (viz_available) {
  source(file.path(CODE_DIR, "imagematrix_visualizer_v2.R"))
}

# --- Parameters ---------------------------------------------------------------
# Input:
#   - "envi": load SAR image from ENVI (.img + .hdr)
#   - "sim" : load simulated image from an .RData file
input_type <- "sim"   # change to "envi" to use ENVI data

# ENVI paths (used when input_type == "envi")
img_path <- file.path(DATA_SAR_DIR, "L16_envi_dublin_size_1100", "Intensity_HH.img")
hdr_path <- file.path(DATA_SAR_DIR, "L16_envi_dublin_size_1100", "Intensity_HH.hdr")

# Simulated RData (used when input_type == "sim")
# sim_var is optional: the object name inside the .RData. If NULL, auto-detects.
sim_rdata <- file.path(DATA_SAR_DIR, "L5_simulated_image_500.Rdata")
sim_var   <- "Z" #NULL   # e.g., "Z" if you know it; otherwise leave NULL

# Tsallis / SAR parameters
L       <- 5            # number of looks
lambda  <- 0.85         # Tsallis order
B       <- 100            # bootstrap replicates
W_min   <- 5            # minimum window size (odd)
W_max   <- 11           # maximum window size (odd)
eta     <- 5           # smoothing parameter for growth/shrink
sigma_n <- 0.523 / sqrt(L)

# p-value parameters / output
p_thr <- 0.05

# --- Load image (ENVI or simulated RData) -------------------------------------
load_simulated_matrix <- function(rdata_path, var_name = NULL) {
  stopifnot(file.exists(rdata_path))
  before <- ls(envir = .GlobalEnv, all.names = TRUE)
  load(rdata_path, envir = .GlobalEnv)
  
  if (!is.null(var_name)) {
    if (!exists(var_name, envir = .GlobalEnv)) {
      stop(sprintf("Object '%s' not found in '%s'.", var_name, rdata_path))
    }
    obj <- get(var_name, envir = .GlobalEnv)
  } else {
    after <- ls(envir = .GlobalEnv, all.names = TRUE)
    new_objs <- setdiff(after, before)
    if (length(new_objs) == 0L) stop(sprintf("No new objects after loading '%s'.", rdata_path))
    pick_obj <- NULL
    for (nm in new_objs) {
      cand <- get(nm, envir = .GlobalEnv)
      if (is.matrix(cand) || is.data.frame(cand) || (is.array(cand) && length(dim(cand)) >= 2)) {
        pick_obj <- cand; break
      }
    }
    if (is.null(pick_obj)) {
      stop(sprintf("No 2D matrix/array-like object found in '%s'. Consider setting sim_var.", rdata_path))
    }
    obj <- pick_obj
  }
  
  if (is.data.frame(obj)) obj <- as.matrix(obj)
  if (is.array(obj) && length(dim(obj)) > 2) obj <- obj[ , , 1, drop = FALSE]  # first band if 3D
  if (!is.matrix(obj)) obj <- as.matrix(obj)
  storage.mode(obj) <- "double"
  obj
}

if (identical(tolower(input_type), "envi")) {
  if (!file.exists(img_path)) stop("ENVI .img not found: ", img_path)
  if (!file.exists(hdr_path)) stop("ENVI .hdr not found: ", hdr_path)
  x <- myread.ENVI(img_path, hdr_path)
  input_label <- "envi"
} else if (identical(tolower(input_type), "sim")) {
  if (!file.exists(sim_rdata)) stop("Simulated .RData not found: ", sim_rdata)
  x <- load_simulated_matrix(sim_rdata, sim_var)
  input_label <- "sim"
} else {
  stop("input_type must be either 'envi' or 'sim'.")
}

# --- Basic shape checks --------------------------------------------------------
if (!is.matrix(x) || any(dim(x) < 3)) {
  stop("Loaded image 'x' must be a 2D numeric matrix with at least 3x3 pixels.")
}
rows <- nrow(x); cols <- ncol(x)
cat(sprintf("Loaded '%s' image with dimensions: %d x %d\n", input_label, rows, cols))

# --- Utilities: clamp, window extraction, border stats, Tsallis test -----------
clamp <- function(v, n) pmin(pmax(v, 1), n)

get_window <- function(i, j, N) {
  r0 <- (i - N):(i + N)
  c0 <- (j - N):(j + N)
  r  <- clamp(r0, rows)
  c  <- clamp(c0, cols)
  x[r, c, drop = FALSE]
}

border_stats <- function(win) {
  Lw      <- nrow(win)
  is_edge <- (row(win) == 1) | (row(win) == Lw) | (col(win) == 1) | (col(win) == Lw)
  vals    <- win[is_edge]
  mu_b    <- mean(vals)
  var_b_r <- mean(vals^2) - mu_b^2
  var_b   <- if (var_b_r < 0) 0 else var_b_r
  list(mu = mu_b, var = var_b)
}

stat_tsallis_window <- function(mat_win) {
  z      <- as.vector(mat_win)
  mu_hat <- mean(z)
  
  est <- bootstrap_tsallis_entropy_optimized(
    z, B = B, lambda = lambda, parallel = FALSE
  )
  
  theo <- (1 - exp((1 - lambda) * log(mu_hat) +
                     (lambda - 1) * log(L) +
                     lgamma(lambda * (L - 1) + 1) -
                     lambda * lgamma(L) -
                     (lambda * (L - 1) + 1) * log(lambda))) /
    (lambda - 1)
  
  est - theo
}

# --- Adaptive pass for a given row i (returns stat and window size L_ij) -------
adaptive_row <- function(i) {
  N_prev   <- (W_min - 1) / 2
  row_stat <- numeric(cols)
  row_L    <- integer(cols)
  
  for (j in seq_len(cols)) {
    win_big <- get_window(i, j, N_prev)
    bs      <- border_stats(win_big)
    Cij     <- sqrt(bs$var) / bs$mu
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

# --- Parallel or sequential execution -----------------------------------------
start_time <- Sys.time()
if (parallel) {
  res_list <- future_lapply(seq_len(rows), adaptive_row, future.seed = TRUE)
} else {
  res_list <- lapply(seq_len(rows), adaptive_row)
}
end_time <- Sys.time()

# --- Reassemble output matrices ------------------------------------------------
difference_values <- t(sapply(res_list, `[[`, "stat"))
L_map             <- t(sapply(res_list, `[[`, "Ls"))

# --- P-values (two-sided normal approximation on standardized differences) -----
# Vectorized version equivalent to your loop; mean()=0 under H0 is implied,
# so we standardize by the sample SD of difference_values.
mean_diff <- mean(difference_values, na.rm = TRUE)
sd_diff   <- sd(difference_values,   na.rm = TRUE)
if (!is.finite(sd_diff) || sd_diff <= 0) {
  stop("Non-positive or non-finite SD for difference_values; cannot compute p-values.")
}

# eps = (diff - 0)/sd ; two-sided p-value = 2*pnorm(-abs(eps))
eps      <- (difference_values - 0) / sd_diff
p_values <- 2 * pnorm(-abs(eps))

# --- Save results (RData) ------------------------------------------------------
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
outfile <- sprintf("%s/%s_adaptive_tsallis_with_Lmap_B%d_eta_%d_%dx%d_%s.Rdata",
                   DATA_DIR,
                   if (identical(input_label, "envi")) "dublin1100" else "simulated",
                   B, eta, W_min, W_max, gsub("\\.", "", format(lambda)))
save(difference_values, p_values, x, L_map, file = outfile)

# --- Timing report and quick checks -------------------------------------------
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat(sprintf("Finished in %.1f sec\n", elapsed))
cat(sprintf("Results saved to '%s'\n", outfile))

h <- floor(elapsed / 3600)
m <- floor((elapsed %% 3600) / 60)
s <- round(elapsed %% 60, 2)
cat(sprintf("Total time: %02d:%02d:%05.2f (hh:mm:ss)\n", h, m, s))

cat("Range of L_ij: ", paste(range(L_map), collapse = " "), "\n")
cat("Window-size frequencies:\n"); print(table(L_map))

# --- Plots: L_map and p-values -------------------------------------------------
#if (viz_available) {
# Using imagematrix_visualizer helpers (if present in Code/)
# Binary mask at alpha = p_thr
mask_title <- sprintf("H (p < %.2f) — %s", p_thr, input_label)
imagematrixPNG(imagematrix(p_values < p_thr),
               name = file.path(DATA_DIR, sprintf("H_mask_p%.2f_%s_%dx%d.png",
                                                  p_thr, input_label, W_min, W_max)))
# Color p-values (viridis-H by default in your helper)
imagematrix_colorPNG(
  imagematrix_color(p_values),
  name           = file.path(DATA_DIR, sprintf("pvalues_color_%s_%dx%d.png",
                                               input_label, W_min, W_max)),
  palette_option = "viridis-H",
  legend_width_px = 540,
  scale_factor = 2,
  direction = -1
)


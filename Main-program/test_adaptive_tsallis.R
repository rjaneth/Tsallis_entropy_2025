# =============================================================================
# test_adaptive_tsallis.R
#
# This script computes the Tsallis entropy–based test statistic on SAR or 
# simulated images using an adaptive windowing strategy. For each pixel, the 
# local window size is adjusted according to scene homogeneity, and the 
# test statistic is calculated. The outputs are:
#   - difference_values: matrix of test statistics
#   - L_map: map of locally selected window sizes
# The results are saved as an .Rdata file for later analysis and visualization.
# =============================================================================

# --- Libraries and parallel plan ------------------------------------------------
req_pkgs <- c("future.apply", "future")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(future.apply)
  library(future)
})

# --- Worker selection (portable) -----------------------------------------------
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
message(sprintf("Plan: %s with %d worker(s).",
                if (parallel) "multisession" else "sequential",
                future::nbrOfWorkers()))

# --- Parameters ----------------------------------------------------------------
input_type <- "sim"   # "envi" or "sim"

# ENVI paths (used when input_type == "envi")
img_path <- "./Data-SAR/L16_envi_dublin_size_1100/Intensity_HH.img"
hdr_path <- "./Data-SAR/L16_envi_dublin_size_1100/Intensity_HH.hdr"

# Simulated RData (used when input_type == "sim")
sim_rdata <- "./Data-SAR/L9_simulated_image_500.Rdata"
sim_var   <- "Z"    # e.g., "Z" if you know the object name; otherwise NULL

# Tsallis / SAR parameters
L       <- 9
lambda  <- 0.85
B       <- 100
W_min   <- 5
W_max   <- 11
eta     <- 5
sigma_n <- 0.523 / sqrt(L)

# --- External functions ---------------------------------------------------------
source("./Code/read_ENVI_images.R")
source("./Code/gamma_sar_sample.R")
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")

# --- Safe loader: load .RData into a local environment -------------------------
load_simulated_matrix <- function(rdata_path, var_name = NULL) {
  stopifnot(file.exists(rdata_path))
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  
  pick_object <- function(obj) {
    # Coerce data.frame/array to matrix when sensible
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
      stop(sprintf("Object '%s' not found in '%s'.", var_name, rdata_path))
    }
    obj <- get(var_name, envir = env, inherits = FALSE)
    obj <- pick_object(obj)
    if (is.null(obj)) stop(sprintf("Object '%s' is not a 2D numeric matrix.", var_name))
    return(obj)
  }
  
  # Auto-detect first 2D matrix-like object
  nms <- ls(envir = env, all.names = TRUE)
  for (nm in nms) {
    cand <- get(nm, envir = env, inherits = FALSE)
    obj  <- pick_object(cand)
    if (!is.null(obj) && is.matrix(obj) && all(dim(obj) >= 3)) {
      return(obj)
    }
  }
  stop(sprintf("No suitable 2D numeric matrix found in '%s'. Consider setting sim_var.", rdata_path))
}

# --- Load image (ENVI or simulated) --------------------------------------------
if (identical(tolower(input_type), "envi")) {
  stopifnot(file.exists(img_path), file.exists(hdr_path))
  x <- myread.ENVI(img_path, hdr_path)
  input_label <- "envi"
} else if (identical(tolower(input_type), "sim")) {
  x <- load_simulated_matrix(sim_rdata, sim_var)
  input_label <- "sim"
} else {
  stop("input_type must be either 'envi' or 'sim'.")
}

# --- Basic shape checks ---------------------------------------------------------
if (!is.matrix(x) || any(dim(x) < 3)) {
  stop("Loaded image 'x' must be a 2D numeric matrix with at least 3x3 pixels.")
}
rows <- nrow(x); cols <- ncol(x)
cat(sprintf("Loaded '%s' image with dimensions: %d x %d\n", input_label, rows, cols))

# --- Utilities: clamp, window extraction, border stats, Tsallis test -----------
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
  mu_b <- mean(vals)
  var_b_raw <- mean(vals^2) - mu_b^2
  var_b <- if (var_b_raw < 0) 0 else var_b_raw
  list(mu = mu_b, var = var_b)
}

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

# --- Adaptive pass for row i (returns stat & per-pixel window side L_ij) -------
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

# --- Execution -----------------------------------------------------------------
start_time <- Sys.time()
res_list <- if (parallel) {
  future_lapply(seq_len(rows), adaptive_row, future.seed = TRUE)
} else {
  lapply(seq_len(rows), adaptive_row)
}
end_time <- Sys.time()

# --- Reassemble outputs ---------------------------------------------------------
difference_values <- t(sapply(res_list, `[[`, "stat"))
L_map             <- t(sapply(res_list, `[[`, "Ls"))

# --- Save results ---------------------------------------------------------------
outdir  <- "./Data"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile <- sprintf("%s/%s_adaptive_tsallis_with_Lmap_B%d_eta_%d_%dx%d_%s.Rdata",
                   outdir,
                   if (identical(input_label, "envi")) "dublin1100" else "simulated",
                   B, eta, W_min, W_max, gsub("\\.", "", format(lambda)))
save(difference_values, x, L_map, file = outfile)

# --- Timing --------------------------------------------------------------------
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat(sprintf("Finished in %.1f sec\n", elapsed))
cat(sprintf("Results saved to '%s'\n", outfile))
h <- floor(elapsed / 3600); m <- floor((elapsed %% 3600) / 60); s <- round(elapsed %% 60, 2)
cat(sprintf("Total time: %02d:%02d:%05.2f (hh:mm:ss)\n", h, m, s))

cat("Range of L_ij: ", paste(range(L_map), collapse = " "), "\n")
cat("Window-size frequencies:\n"); print(table(L_map))

# Optional:
# on.exit(try(future::plan(sequential), silent = TRUE), add = TRUE)

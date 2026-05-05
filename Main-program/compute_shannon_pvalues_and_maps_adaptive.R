# =============================================================================
# compute_shannon_pvalues_and_maps_adaptive.R
#
# Adaptive Shannon p-values:
#   z_ij = (S_ij - mu0(W_ij)) / sd_obs
#   p_ij = 2 * pnorm(-abs(z_ij))
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Parameters
# -----------------------------------------------------------------------------
result_file <- "./Data/results/shannon1/dublin_shannon_adaptive_boot_eta2_w5-11_L16_B50.RData"
alpha       <- 0.05
calibration_file <- "./Data/null_diagnostic_shannon/shannon_null_fixed_summary.csv"

figdir <- "./Data/figures/shannon1"
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
if (!dir.exists("./Data/pvalues/shannon1")) dir.create("./Data/pvalues/shannon1", recursive = TRUE)
if (!dir.exists("./Data/decisions/shannon1")) dir.create("./Data/decisions/shannon1", recursive = TRUE)

img_binary_name  <- file.path(figdir, "dublin_shannon_adaptive_boot_eta3_alpha005.png")
img_pvalues_name <- file.path(figdir, "dublin_shannon_adaptive_boot_eta3_pvalues.png")

# -----------------------------------------------------------------------------
# 1. Utilities
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a)) ||
      (is.character(a) && length(a) == 1L && !nzchar(a))) b else a
}

load_result_object <- function(filepath) {
  env <- new.env(parent = emptyenv())
  load(filepath, envir = env)
  
  if (exists("result", envir = env, inherits = FALSE)) {
    return(get("result", envir = env))
  }
  
  stop("No 'result' object found in the file.")
}

load_calibration_table <- function(filepath) {
  cal <- read.csv(filepath, stringsAsFactors = FALSE)
  req <- c("W", "mu0_hat", "sigma0_hat")
  miss <- setdiff(req, names(cal))
  if (length(miss) > 0) {
    stop(sprintf("Missing columns in calibration file: %s", paste(miss, collapse = ", ")))
  }
  
  cal$W <- as.integer(cal$W)
  cal$mu0_hat <- as.numeric(cal$mu0_hat)
  cal$sigma0_hat <- as.numeric(cal$sigma0_hat)
  
  if ("method" %in% names(cal)) {
    cal <- cal[cal$method == "fixed", ]
  }
  
  cal[, req]
}

# -----------------------------------------------------------------------------
# 2. Load inputs
# -----------------------------------------------------------------------------
res <- load_result_object(result_file)
calibration_tbl <- load_calibration_table(calibration_file)

if (!identical(res$window_mode, "adaptive")) {
  stop("This script is only for adaptive results.")
}
if (is.null(res$W_map)) {
  stop("W_map was not found in result.")
}

diff_map <- res$diff_map
W_map    <- res$W_map

# -----------------------------------------------------------------------------
# 3. Compute adaptive p-values
# -----------------------------------------------------------------------------
finite_vals <- diff_map[is.finite(diff_map)]
sd_obs <- sd(finite_vals, na.rm = TRUE)
if (!is.finite(sd_obs) || sd_obs <= 0) stop("sd_obs is not valid.")

z_map <- matrix(NA_real_, nrow = nrow(diff_map), ncol = ncol(diff_map))
p_values <- matrix(NA_real_, nrow = nrow(diff_map), ncol = ncol(diff_map))

for (w in sort(unique(as.vector(W_map)))) {
  idx <- (W_map == w)
  mu0 <- calibration_tbl$mu0_hat[calibration_tbl$W == w]
  if (length(mu0) == 0) stop(sprintf("No mu0 calibration found for W = %d", w))
  
  vals <- diff_map[idx]
  good <- is.finite(vals)
  
  z_tmp <- rep(NA_real_, length(vals))
  p_tmp <- rep(NA_real_, length(vals))
  
  z_tmp[good] <- (vals[good] - mu0) / sd_obs
  p_tmp[good] <- 2 * pnorm(-abs(z_tmp[good]))
  
  z_map[idx]    <- z_tmp
  p_values[idx] <- p_tmp
}

decisions <- ifelse(!is.na(p_values) & p_values < alpha, 1L, 0L)

# -----------------------------------------------------------------------------
# 4. Save outputs
# -----------------------------------------------------------------------------
pval_result <- list(
  p_values        = p_values,
  z_map           = z_map,
  diff_map        = diff_map,
  W_map           = W_map,
  calibration_tbl = calibration_tbl,
  alpha           = alpha,
  entropy         = "shannon",
  bootstrap       = res$bootstrap,
  B               = res$B,
  window_mode     = "adaptive",
  W_min           = res$W_min,
  W_max           = res$W_max,
  crop_margin     = res$crop_margin,
  eta             = res$eta,
  L               = res$L,
  scene           = res$scene,
  image_id        = res$image_id,
  pvalue_rule     = "adaptive_mu0_sdobs",
  sd_obs          = sd_obs
)

decision_result <- list(
  decisions     = decisions,
  W_map         = W_map,
  alpha         = alpha,
  entropy       = "shannon",
  bootstrap     = res$bootstrap,
  B             = res$B,
  window_mode   = "adaptive",
  W_min         = res$W_min,
  W_max         = res$W_max,
  crop_margin   = res$crop_margin,
  eta           = res$eta,
  L             = res$L,
  scene         = res$scene,
  image_id      = res$image_id
)

tag_boot <- ifelse(isTRUE(res$bootstrap), "boot", "noboot")

pval_file <- sprintf(
  "./Data/pvalues/shannon1/%s_shannon_adaptive_%s_eta%d_pvalues_v2.RData",
  res$scene %||% "scene",
  tag_boot,
  res$eta
)

decision_file <- sprintf(
  "./Data/decisions/shannon1/%s_shannon_adaptive_%s_eta%d_alpha_%s_decisions_v2.RData",
  res$scene %||% "scene",
  tag_boot,
  res$eta,
  gsub("\\.", "", format(alpha))
)

save(pval_result, file = pval_file)
save(decision_result, file = decision_file)

# -----------------------------------------------------------------------------
# 5. Maps
# -----------------------------------------------------------------------------
source("imagematrix_visualizer_v2.R")

imagematrixPNG(
  imagematrix(decisions),
  name = img_binary_name
)

imagematrix_colorPNG(
  imagematrix_color(p_values),
  name = img_pvalues_name,
  palette_option = "viridis-H",
  legend_width_px = 850,
  scale_factor = 2,
  direction = -1
)

cat(sprintf("Saved p-values to: %s\n", pval_file))
cat(sprintf("Saved decisions to: %s\n", decision_file))
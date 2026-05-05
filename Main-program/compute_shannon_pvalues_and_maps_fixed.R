# =============================================================================
# compute_shannon_pvalues_and_maps_fixed.R
#
# Fixed Shannon p-values using the original practical rule:
#   z_ij = S_ij / sd_diff
#   p_ij = 2 * pnorm(-abs(z_ij))
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Parameters
# -----------------------------------------------------------------------------
result_file <- "./Data/results/shannon/dublin_shannon_fixed_boot_w7_L16_B100.RData"
alpha       <- 0.05

figdir <- "./Data/figures/shannon"
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
if (!dir.exists("./Data/pvalues/shannon")) dir.create("./Data/pvalues/shannon", recursive = TRUE)
if (!dir.exists("./Data/decisions/shannon")) dir.create("./Data/decisions/shannon", recursive = TRUE)

img_binary_name  <- file.path(figdir, "dublin_shannon_fixed_boot_w7_alpha005.png")
img_pvalues_name <- file.path(figdir, "dublin_shannon_fixed_boot_w7_pvalues.png")

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

calculate_p_values_matrix <- function(data_matrix, sigma) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  p_values_matrix <- matrix(NA_real_, nrow = rows, ncol = cols)
  
  for (i in seq_len(rows)) {
    for (j in seq_len(cols)) {
      val <- data_matrix[i, j]
      
      if (is.finite(val)) {
        epsilon <- val / sigma
        p_values_matrix[i, j] <- 2 * pnorm(-abs(epsilon))
      } else {
        p_values_matrix[i, j] <- NA_real_
      }
    }
  }
  
  p_values_matrix
}

# -----------------------------------------------------------------------------
# 2. Load result
# -----------------------------------------------------------------------------
res <- load_result_object(result_file)

if (!identical(res$window_mode, "fixed")) {
  stop("This script is only for fixed results.")
}

diff_map <- res$diff_map
if (is.null(diff_map) || !is.matrix(diff_map)) {
  stop("The file does not contain a valid diff_map.")
}

finite_vals <- diff_map[is.finite(diff_map)]

if (length(finite_vals) < 2L) {
  stop("Not enough finite values in diff_map to compute p-values.")
}

mean_diff <- mean(finite_vals)
sd_diff   <- sd(finite_vals)

if (!is.finite(sd_diff) || sd_diff <= 0) {
  stop("sd_diff is not valid.")
}

# -----------------------------------------------------------------------------
# 3. Compute p-values and decisions
# -----------------------------------------------------------------------------
diff_map_clean <- diff_map
diff_map_clean[!is.finite(diff_map_clean)] <- NA_real_

p_values  <- calculate_p_values_matrix(diff_map_clean, sigma = sd_diff)
decisions <- ifelse(!is.na(p_values) & p_values < alpha, 1L, 0L)

# -----------------------------------------------------------------------------
# 4. Save outputs
# -----------------------------------------------------------------------------
pval_result <- list(
  p_values     = p_values,
  diff_map     = diff_map_clean,
  mean_diff    = mean_diff,
  sd_diff      = sd_diff,
  alpha        = alpha,
  entropy      = "shannon",
  bootstrap    = res$bootstrap,
  B            = res$B,
  window_mode  = "fixed",
  window_size  = res$window_size,
  crop_margin  = res$crop_margin,
  L            = res$L,
  scene        = res$scene,
  image_id     = res$image_id,
  pvalue_rule  = "fixed_sdobs_original"
)

decision_result <- list(
  decisions     = decisions,
  alpha         = alpha,
  entropy       = "shannon",
  bootstrap     = res$bootstrap,
  B             = res$B,
  window_mode   = "fixed",
  window_size   = res$window_size,
  crop_margin   = res$crop_margin,
  L             = res$L,
  scene         = res$scene,
  image_id      = res$image_id
)

tag_boot <- ifelse(isTRUE(res$bootstrap), "boot", "noboot")

pval_file <- sprintf(
  "./Data/pvalues/shannon/%s_shannon_fixed_%s_w%d_pvalues.RData",
  res$scene %||% "scene",
  tag_boot,
  res$window_size
)

decision_file <- sprintf(
  "./Data/decisions/shannon/%s_shannon_fixed_%s_w%d_alpha_%s_decisions.RData",
  res$scene %||% "scene",
  tag_boot,
  res$window_size,
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
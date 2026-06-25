# =============================================================================
# compute_tsallis_null_adaptive_pvalue_diagnostic.R
#
# Final adaptive p-value diagnostic under H0 for the Tsallis test.
#
# Uses the output of:
#   run_tsallis_null_diagnostic_sim.R
#
# Produces:
#   1) histogram of final adaptive p-values under H0
#   2) empirical rejection rates at 1%, 5%, and 10%
#   3) CSV and RData outputs
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
req_pkgs <- c("ggplot2")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Input and output paths
# -----------------------------------------------------------------------------
master_file <- "./Data/tsallis_null_diagnostic_master.RData"

outdir <- "./Data/pvalue_diagnostic"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

figdir <- "./Figures-R1"
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. Load master result
# -----------------------------------------------------------------------------
if (!file.exists(master_file)) {
  stop("Master file not found. Run run_tsallis_null_diagnostic_sim.R first.")
}

load(master_file)

if (!exists("master_result")) {
  stop("Object 'master_result' was not found in the loaded file.")
}

cat("Loaded master_result\n")
cat("L =", master_result$L, "\n")
cat("eta =", master_result$eta, "\n")
cat("lambda =", master_result$lambda, "\n")
cat("B =", master_result$B, "\n")
cat("W_min =", master_result$W_min, "\n")
cat("W_max =", master_result$W_max, "\n")

# -----------------------------------------------------------------------------
# 3. Extract adaptive statistic map and W map
# -----------------------------------------------------------------------------
adaptive_result <- master_result$adaptive_result

S_map <- adaptive_result$diff_map
W_map <- adaptive_result$W_map
margin <- adaptive_result$crop_margin

if (!all(dim(S_map) == dim(W_map))) {
  stop("S_map and W_map must have the same dimensions.")
}

# -----------------------------------------------------------------------------
# 4. Helper: extract interior region
# -----------------------------------------------------------------------------
extract_interior_matrix <- function(mat, margin) {
  if (margin <= 0) return(mat)
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  mat[(margin + 1):(nr - margin),
      (margin + 1):(nc - margin),
      drop = FALSE]
}

S_int <- extract_interior_matrix(S_map, margin)
W_int <- extract_interior_matrix(W_map, margin)

cat("Interior dimensions:", paste(dim(S_int), collapse = " x "), "\n")
cat("Selected W in interior:\n")
print(table(W_int))

# -----------------------------------------------------------------------------
# 5. Window-dependent means under H0
# -----------------------------------------------------------------------------
cal_tbl <- master_result$fixed_summary[, c("W", "mu0_hat", "sigma0_hat")]
cal_tbl$W <- as.integer(cal_tbl$W)

cat("Null summary table:\n")
print(cal_tbl)

# -----------------------------------------------------------------------------
# 6. Observed scale sigma_obs
# -----------------------------------------------------------------------------
S_vec <- as.numeric(S_int)
S_vec <- S_vec[is.finite(S_vec)]

sigma_obs <- sd(S_vec)

if (!is.finite(sigma_obs) || sigma_obs <= 0) {
  stop("sigma_obs is not valid.")
}

cat(sprintf("sigma_obs = %.8f\n", sigma_obs))

# -----------------------------------------------------------------------------
# 7. Compute final adaptive z-values and p-values under H0
# -----------------------------------------------------------------------------
z_map <- matrix(NA_real_, nrow = nrow(S_int), ncol = ncol(S_int))
p_map <- matrix(NA_real_, nrow = nrow(S_int), ncol = ncol(S_int))

observed_W <- sort(unique(as.vector(W_int)))
observed_W <- observed_W[is.finite(observed_W)]

for (w in observed_W) {
  mu_w <- cal_tbl$mu0_hat[cal_tbl$W == w]
  
  if (length(mu_w) != 1 || !is.finite(mu_w)) {
    stop(sprintf("No valid mu0_hat found for W = %s.", w))
  }
  
  idx <- (W_int == w) & is.finite(S_int)
  
  z_map[idx] <- (S_int[idx] - mu_w) / sigma_obs
  p_map[idx] <- 2 * pnorm(-abs(z_map[idx]))
}

p_values <- as.numeric(p_map)
p_values <- p_values[is.finite(p_values)]

cat("Number of valid p-values:", length(p_values), "\n")
cat("Summary of p-values:\n")
print(summary(p_values))

# -----------------------------------------------------------------------------
# 8. Empirical rejection rates
# -----------------------------------------------------------------------------
alpha_levels <- c(0.01, 0.05, 0.10)

rejection_table <- data.frame(
  alpha = alpha_levels,
  empirical_rejection_rate = sapply(alpha_levels, function(a) mean(p_values < a)),
  n_valid_pixels = length(p_values)
)

print(rejection_table)

write.csv(
  rejection_table,
  file = file.path(outdir, "tsallis_null_adaptive_pvalue_rejection_rates.csv"),
  row.names = FALSE
)

# LaTeX-friendly version
rejection_table_latex <- data.frame(
  alpha = sprintf("%.2f", rejection_table$alpha),
  empirical_rejection_rate = sprintf("%.4f", rejection_table$empirical_rejection_rate),
  n_valid_pixels = rejection_table$n_valid_pixels
)

write.csv(
  rejection_table_latex,
  file = file.path(outdir, "tsallis_null_adaptive_pvalue_rejection_rates_latex.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. Histogram against Uniform(0,1)
# -----------------------------------------------------------------------------
hist_df <- data.frame(p_value = p_values)

p_hist <- ggplot(hist_df, aes(x = p_value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "grey35"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Adaptive p-value under H0",
    y = "Density"
  ) +
  theme_minimal(base_size = 13)

breaks_p <- seq(0, 1, length.out = 31)

p_hist <- ggplot(hist_df, aes(x = p_value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = breaks_p,
    closed = "right",
    fill = "grey85",
    color = "grey35"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Adaptive p-value under H0",
    y = "Density"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "tsallis_null_adaptive_pvalue_histogram1.pdf"),
  plot = p_hist,
  width = 6.5,
  height = 4.5
)

ggsave(
  filename = file.path(outdir, "tsallis_null_adaptive_pvalue_histogram.png"),
  plot = p_hist,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

# Copy also to Figures-R1 for the manuscript
ggsave(
  filename = file.path(figdir, "tsallis_null_adaptive_pvalue_histogram.pdf"),
  plot = p_hist,
  width = 6.5,
  height = 4.5
)

# -----------------------------------------------------------------------------
# 10. Save full diagnostic object
# -----------------------------------------------------------------------------
pvalue_diagnostic <- list(
  p_values = p_values,
  p_map = p_map,
  z_map = z_map,
  sigma_obs = sigma_obs,
  rejection_table = rejection_table,
  calibration_table = cal_tbl,
  L = master_result$L,
  eta = master_result$eta,
  lambda = master_result$lambda,
  B = master_result$B,
  W_min = master_result$W_min,
  W_max = master_result$W_max,
  margin = margin,
  W_frequency_interior = table(W_int)
)

save(
  pvalue_diagnostic,
  file = file.path(outdir, "tsallis_null_adaptive_pvalue_diagnostic.RData")
)

cat("\nSaved outputs in:\n")
cat(normalizePath(outdir), "\n")
cat("\nFigure also saved in:\n")
cat(normalizePath(figdir), "\n")
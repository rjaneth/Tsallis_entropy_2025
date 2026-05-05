# =============================================================================
# compute_shannon_threshold_sweep_and_auc.R
#
# Generic script for Shannon fixed/adaptive results.
# It computes:
#   - threshold sweep (0.01 to 0.20)
#   - summary metrics at 0.05
#   - ROC curve and AUC
#   - roc_df for later combined plots
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Parameters
# -----------------------------------------------------------------------------
sar_path    <- "./Data-SAR/L16_dublin_1600_nuevo/Intensity_HH.img"
pvalues_rda <- "./Data/pvalues/shannon1/dublin_shannon_adaptive_boot_eta2_pvalues_v2.RData"
rois_path   <- "./Data-paper/rois_dublin1600.gpkg"

outdir <- "./Data/metrics/shannon1"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

thresholds <- seq(0.01, 0.20, by = 0.01)

# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------
req_pkgs <- c("terra", "sf", "caret", "pROC", "ggplot2")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(caret)
  library(pROC)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
safe_scalar <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || is.na(x)) NA_real_ else x
}

load_single_object <- function(filepath) {
  env <- new.env(parent = emptyenv())
  load(filepath, envir = env)
  nms <- ls(envir = env)
  if ("pval_result" %in% nms) return(get("pval_result", envir = env))
  get(nms[1], envir = env)
}

build_raster_from_matrix <- function(mat, sar, crop_margin = 0) {
  rotated_matrix <- t(apply(mat, 2, rev))
  
  nr <- nrow(rotated_matrix)
  nc <- ncol(rotated_matrix)
  resx <- res(sar)[1]
  resy <- res(sar)[2]
  
  ext_new <- ext(
    xmin(sar) + crop_margin * resx,
    xmax(sar) - crop_margin * resx,
    ymin(sar) + crop_margin * resy,
    ymax(sar) - crop_margin * resy
  )
  
  r_out <- rast(nrows = nr, ncols = nc, ext = ext_new, crs = crs(sar))
  values(r_out) <- as.vector(rotated_matrix)
  flip(r_out, direction = "vertical")
}

compute_metrics_at_threshold <- function(y_true, p_values, thr) {
  y_pred <- ifelse(p_values < thr, 1L, 0L)
  
  confR <- confusionMatrix(
    as.factor(y_pred),
    as.factor(y_true),
    positive = "1"
  )
  
  Recall      <- safe_scalar(confR$byClass["Sensitivity"])
  Specificity <- safe_scalar(confR$byClass["Specificity"])
  Precision   <- safe_scalar(confR$byClass["Pos Pred Value"])
  F1          <- safe_scalar(confR$byClass["F1"])
  Kappa       <- safe_scalar(confR$overall["Kappa"])
  OA          <- safe_scalar(confR$overall["Accuracy"])
  Pd          <- Recall
  Pfa         <- if (is.na(Specificity)) NA_real_ else 1 - Specificity
  
  data.frame(
    Threshold = thr,
    Precision = Precision,
    Recall    = Recall,
    Pd        = Pd,
    Pfa       = Pfa,
    F1        = F1,
    Kappa     = Kappa,
    OA        = OA
  )
}

# -----------------------------------------------------------------------------
# 3. Load data
# -----------------------------------------------------------------------------
sar  <- rast(sar_path)
pobj <- load_single_object(pvalues_rda)

p_values    <- pobj$p_values
crop_margin <- if (!is.null(pobj$crop_margin)) pobj$crop_margin else 0
scene       <- if (!is.null(pobj$scene)) pobj$scene else "scene"
entropy     <- if (!is.null(pobj$entropy)) pobj$entropy else "shannon"
bootstrap   <- if (!is.null(pobj$bootstrap)) pobj$bootstrap else NA
window_mode <- if (!is.null(pobj$window_mode)) pobj$window_mode else NA
window_size <- if (!is.null(pobj$window_size)) pobj$window_size else NA
eta         <- if (!is.null(pobj$eta)) pobj$eta else NA
pvalue_rule <- if (!is.null(pobj$pvalue_rule)) pobj$pvalue_rule else NA_character_

# -----------------------------------------------------------------------------
# 4. Rasterize and get ROI vectors
# -----------------------------------------------------------------------------
r_p <- build_raster_from_matrix(p_values, sar = sar, crop_margin = crop_margin)

rois <- st_read(rois_path, quiet = TRUE)
rois <- st_transform(rois, crs(sar))

gt <- rasterize(rois, r_p, field = "class")

idx <- !is.na(values(gt)) & !is.na(values(r_p))
y_true <- values(gt)[idx]
p_vec  <- values(r_p)[idx]

y_score <- 1 - p_vec

# -----------------------------------------------------------------------------
# 5. Threshold sweep
# -----------------------------------------------------------------------------
sweep_list <- lapply(thresholds, function(thr) {
  compute_metrics_at_threshold(y_true = y_true, p_values = p_vec, thr = thr)
})

sweep_df <- do.call(rbind, sweep_list)

# -----------------------------------------------------------------------------
# 6. ROC / AUC
# -----------------------------------------------------------------------------
roc_obj <- roc(
  response  = y_true,
  predictor = y_score,
  quiet = TRUE,
  direction = "<"
)

auc_value <- as.numeric(auc(roc_obj))

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities,
  Threshold = roc_obj$thresholds,
  Scene = tools::toTitleCase(scene),
  Entropy = tools::toTitleCase(entropy),
  Bootstrap = ifelse(isTRUE(bootstrap), "Yes", "No"),
  Window = ifelse(identical(window_mode, "fixed"), "Fixed", "Adaptive"),
  WindowSize = ifelse(identical(window_mode, "fixed"), window_size, NA),
  Eta = ifelse(identical(window_mode, "adaptive"), eta, NA),
  PValueRule = pvalue_rule,
  stringsAsFactors = FALSE
)

summary_df <- data.frame(
  Scene      = tools::toTitleCase(scene),
  Entropy    = tools::toTitleCase(entropy),
  Bootstrap  = ifelse(isTRUE(bootstrap), "Yes", "No"),
  Window     = ifelse(identical(window_mode, "fixed"), "Fixed", "Adaptive"),
  WindowSize = ifelse(identical(window_mode, "fixed"), window_size, NA),
  Eta        = ifelse(identical(window_mode, "adaptive"), eta, NA),
  PValueRule = pvalue_rule,
  AUC        = auc_value,
  stringsAsFactors = FALSE
)

row_005 <- sweep_df[abs(sweep_df$Threshold - 0.05) < 1e-12, ]
summary_df$F1_005         <- row_005$F1
summary_df$Kappa_005      <- row_005$Kappa
summary_df$OA_005         <- row_005$OA
summary_df$Recall_005     <- row_005$Recall
summary_df$Precision_005  <- row_005$Precision
summary_df$Pd_005         <- row_005$Pd
summary_df$Pfa_005        <- row_005$Pfa

# -----------------------------------------------------------------------------
# 7. Method label
# -----------------------------------------------------------------------------
method_label <- paste0(
  "Shannon_",
  if (identical(window_mode, "fixed")) paste0("fixed_w", window_size) else paste0("adaptive_eta", eta),
  "_",
  if (isTRUE(bootstrap)) "boot" else "noboot"
)

# -----------------------------------------------------------------------------
# 8. Plots
# -----------------------------------------------------------------------------
p_sweep <- ggplot(sweep_df, aes(x = Threshold)) +
  geom_line(aes(y = F1, color = "F1"), linewidth = 1) +
  geom_line(aes(y = Kappa, color = "Kappa"), linewidth = 1) +
  geom_line(aes(y = OA, color = "OA"), linewidth = 1) +
  labs(
    title = paste("Threshold sweep -", tools::toTitleCase(scene), "-", tools::toTitleCase(entropy)),
    y = "Metric value",
    color = "Metric"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(outdir, paste0(tolower(scene), "_", tolower(entropy), "_", method_label, "_threshold_sweep.png")),
  plot = p_sweep, width = 8, height = 5, dpi = 300
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.4) +
  annotate(
    "text",
    x = 0.72, y = 0.08,
    label = paste0("AUC = ", sprintf("%.4f", auc_value)),
    size = 5
  ) +
  labs(
    title = paste("ROC curve -", tools::toTitleCase(scene), "-", tools::toTitleCase(entropy)),
    x = "Probability of false alarm",
    y = "Probability of detection"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(outdir, paste0(tolower(scene), "_", tolower(entropy), "_", method_label, "_roc.png")),
  plot = p_roc, width = 6.5, height = 5.5, dpi = 300
)

# -----------------------------------------------------------------------------
# 9. Save
# -----------------------------------------------------------------------------
outfile <- file.path(
  outdir,
  paste0(
    tolower(scene), "_", tolower(entropy), "_", method_label,
    "_threshold_sweep_auc_roc.RData"
  )
)

save(
  sweep_df,
  summary_df,
  roc_obj,
  roc_df,
  file = outfile
)

write.csv(sweep_df,  sub("\\.RData$", "_sweep.csv", outfile), row.names = FALSE)
write.csv(summary_df, sub("\\.RData$", "_summary.csv", outfile), row.names = FALSE)
write.csv(roc_df,     sub("\\.RData$", "_roc.csv", outfile), row.names = FALSE)

print(summary_df)
cat(sprintf("\nSaved to: %s\n", outfile))
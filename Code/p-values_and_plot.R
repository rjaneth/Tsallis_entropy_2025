# =============================================================================
# p-values_and_plot.R
# Load difference_values from .RData, compute p-values, and plot/save outputs.
# Run this file from Main_program/.
# =============================================================================

# --- Install-if-missing 
ensure <- function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# ensure("viridisLite")  # only if your viewer needs it

# --- Paths 
get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m  <- grep("^--file=", ca)
  if (length(m) > 0) return(dirname(normalizePath(sub("^--file=", "", ca[m]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  if (requireNamespace("knitr", quietly = TRUE)) {
    p <- tryCatch(knitr::current_input(), error = function(e) NULL)
    if (!is.null(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}

SCRIPT_DIR <- normalizePath(get_script_dir(), winslash = "/", mustWork = FALSE)

# Code/ or Codes/
possible_code_dirs <- c(file.path(SCRIPT_DIR, "Code"), file.path(SCRIPT_DIR, "Codes"))
CODE_DIR <- possible_code_dirs[dir.exists(possible_code_dirs)][1]
if (is.na(CODE_DIR)) stop("Missing 'Code' (or 'Codes') folder inside Main_program/.")

DATA_DIR <- file.path(SCRIPT_DIR, "Data")
if (!dir.exists(DATA_DIR)) stop("Missing 'Data' folder inside Main_program/.")

cat("SCRIPT_DIR :", SCRIPT_DIR, "\n")
cat("CODE_DIR   :", CODE_DIR,    "\n")
cat("DATA_DIR   :", DATA_DIR,    "\n")

# --- Choose the .RData --------------------------------------------------------
# Option A: set it explicitly (comment out to use Option B)
 rdata_file <- file.path(DATA_DIR, "munich1024_adaptive_tsallis_with_Lmap_B100_eta_4_5x11_09nuevo.Rdata")

# Option B: auto-pick the most recent .Rdata in Data/
# if (!exists("rdata_file")) {
#   rds <- list.files(DATA_DIR, pattern = "\\.Rdata$", full.names = TRUE)
#   if (length(rds) == 0) stop("No .Rdata files found in Data/.")
#   info <- file.info(rds)
#   rdata_file <- rds[order(info$mtime, decreasing = TRUE)][1]
#   message("Auto-selected: ", basename(rdata_file))
# }

p_thr <- 0.05  # threshold for binary mask

# --- Load .RData into a clean env 
if (!file.exists(rdata_file)) stop("File not found: ", rdata_file)
env <- new.env(parent = emptyenv())
obj_names <- load(rdata_file, envir = env)
if (!"difference_values" %in% obj_names) {
  stop(sprintf("No 'difference_values' in '%s'. Available: %s",
               basename(rdata_file), paste(obj_names, collapse = ", ")))
}
difference_values <- get("difference_values", envir = env)

# --- Compute p-values (two-sided, z-score) 
sd_diff <- sd(difference_values, na.rm = TRUE)
if (!is.finite(sd_diff) || sd_diff <= 0) stop("Invalid SD for difference_values: ", sd_diff)
p_values <- 2 * pnorm(-abs(difference_values / sd_diff))

# --- Try custom visualizer; fallback to base R plots --------------------------
viz_file <- file.path(CODE_DIR, "imagematrix_visualizer_v2.R")
viz_file <- file.path(CODE_DIR, "iimagematrix_visualizer_blues_v1.R")
if (file.exists(viz_file)) {
  source(viz_file)  # must define: imagematrix(), imagematrixPNG(), imagematrix_color(), imagematrix_colorPNG()
  
  # Binary mask (p < p_thr)
  out_mask <- file.path(DATA_DIR, sprintf("H_mask-munich13%.2f.png", p_thr))
  imagematrixPNG(imagematrix(p_values < p_thr), name = out_mask)
  
  # Color map of p-values (discrete blues; custom breaks)
  out_color <- file.path(DATA_DIR, "p-values_color-munich13.png")
  
  
  # custom_breaks <- c(0, 0.05, 0.2, 0.4, 0.6, 0.8, 1)  # edit if you want
  # imagematrix_colorPNG(
  #   imagematrix_color(p_values),
  #   name               = out_color,
  #   significance_level = 0.05,
  #   discrete           = TRUE,
  #   sig_blues          = TRUE,
  #   breaks             = custom_breaks,
  #   legend_width_px    = 540,
  #   scale_factor       = 2
  # )
  
  #alpha <- 0.05  # significance level
  # 
  # #custom_breaks  <- c(0, alpha, 0.2, 0.4, 0.6, 0.8, 0.9, 1)
  # #custom_colors  <- c("black", "purple", "red", "orange", "lightgreen", "lightskyblue", "navy")
  
  # custom_breaks  <- c(0, alpha, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  # custom_colors  <- c("black","purple","darkgreen","orange","red","yellow","coral4",
  #                     "deeppink","green","lightskyblue","navy")#"coral4",
  # 
  # # PNG
  # imagematrix_colorPNG(
  #   imagematrix_color(p_values),
  #   name              = out_color,
  #   discrete          = TRUE,
  #   breaks            = custom_breaks,
  #   colors            = custom_colors,
  #   significance_level= alpha,
  #   legend_width_px   = 540,
  #   scale_factor      = 2
  # )
  
  imagematrix_colorPNG(
    imagematrix_color(p_values),
    name            = out_color,
    palette_option  = "viridis-H",
    legend_width_px = 540,
    scale_factor    = 2,
    #discrete = TRUE,
    direction = -1
  )
  
  cat("Saved:", out_mask,  "\n")
  cat("Saved:", out_color, "\n")
  
} else {
  message("Visualizer not found: ", basename(viz_file), " in ", CODE_DIR, " — using base plots.")
  
  # Binary plot
  out_mask <- file.path(DATA_DIR, sprintf("pvalues_binary_p%.2f.png", p_thr))
  png(out_mask, width = 1200, height = 900, res = 150)
  image(p_values < p_thr, col = c("white", "black"), axes = FALSE,
        main = sprintf("H (p < %.2f)", p_thr))
  dev.off()
  
  # Grayscale map
  out_gray <- file.path(DATA_DIR, "pvalues_gray.png")
  png(out_gray, width = 1200, height = 900, res = 150)
  image(p_values, col = gray.colors(256), axes = FALSE, main = "p-values (gray)")
  dev.off()
  
  cat("Saved:", out_mask, "\n")
  cat("Saved:", out_gray, "\n")
}

# --- Save p-values for reuse --------------------------------------------------
saveRDS(p_values, file = file.path(DATA_DIR, "p_values.rds"))
message("P-value outputs saved in ", DATA_DIR)

# =============================================================================
# Optional: p-value binning utilities (for tables/QA)
# =============================================================================

# Default 11 intervals (12 cut points): [0,1]; tweak as needed
.pval_breaks_default <- c(0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.60, 0.80, 0.90, 1.00)

.label_breaks <- function(breaks, right = FALSE, digits = 3) {
  stopifnot(is.numeric(breaks), length(breaks) >= 2)
  fmt <- function(v) sub("0+$", "", sub("\\.$", "", format(round(v, digits), nsmall = digits)))
  L <- head(breaks, -1); U <- tail(breaks, -1)
  if (right) paste0("(", fmt(L), ", ", fmt(U), "]") else paste0("[", fmt(L), ", ", fmt(U), ")")
}

bin_pvals <- function(mat,
                      breaks = .pval_breaks_default,
                      right = FALSE,            # [a,b) by default
                      include.lowest = TRUE,
                      na.value = NA_integer_) {
  stopifnot(is.matrix(mat) || is.array(mat))
  stopifnot(is.numeric(breaks), all(is.finite(breaks)),
            isTRUE(all.equal(min(breaks), 0)), isTRUE(all.equal(max(breaks), 1)))
  if (!is.numeric(mat)) stop("mat must be numeric in [0,1].")
  
  v  <- as.vector(mat)
  f  <- cut(v, breaks = breaks, right = right, include.lowest = include.lowest)
  idx <- as.integer(f)
  idx[is.na(v)] <- na.value
  dim(idx) <- dim(mat)
  
  labs <- .label_breaks(breaks, right = right)
  counts <- as.integer(table(f))
  names(counts) <- labs
  total <- sum(counts)
  proportions <- counts / if (total > 0) total else 1
  
  out_tab <- data.frame(
    interval    = factor(labs, levels = labs),
    count       = counts,
    proportion  = proportions,
    stringsAsFactors = FALSE
  )
  
  list(idx = idx, labels = labs, counts = counts, proportions = proportions, table = out_tab)
}

# Example table 
 res <- bin_pvals(p_values)
 print(res$table)

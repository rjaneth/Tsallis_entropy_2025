##############################################################################
# 0.  Paths
##############################################################################
sar_path      <- "./Data-SAR/L12_envi_munich_size_1024/Intensity_HH.img"
decisions_rda <-  "./Data/decisions_tsallis_munich_eta1_alpha005.Rdata"
rois_path     <- "./Data_paper/rois_munich.gpkg"

##############################################################################
# 1.  Load SAR
##############################################################################
library(terra)
sar <- rast(sar_path)

##############################################################################
# 2.  Load matrix and create raster with rotation + vertical flip
##############################################################################
load(decisions_rda)  # object: decisions_tsallis_munich_eta1

# Rotate 90° clockwise
rotated_matrix <- t(apply(decisions_tsallis_munich_eta1, 2, rev))

# Dimensions and resolution
nr <- nrow(rotated_matrix)
nc <- ncol(rotated_matrix)
resx <- res(sar)[1]
resy <- res(sar)[2]

# Adjusted extent (crop 3 pixels per side)
ext_new <- ext(xmin(sar) + 3*resx,
               xmax(sar) - 3*resx,
               ymin(sar) + 3*resy,
               ymax(sar) - 3*resy)

# Create aligned raster
r_dec <- rast(nrows = nr, ncols = nc,
              ext   = ext_new,
              crs   = crs(sar))

# Assign rotated values
values(r_dec) <- as.vector(rotated_matrix)

# Vertical flip
r_dec <- flip(r_dec, direction = "vertical")

# Save (optional, useful for inspection in QGIS)
writeRaster(r_dec,
            filename = "./Data/decisions_tsallis_munich_eta1_n.tif",
            datatype = "INT1U",
            overwrite = TRUE)

##############################################################################
# 3.  Read ROIs and rasterize them
##############################################################################
library(sf)
rois <- st_read(rois_path)
rois <- st_transform(rois, crs(sar))   # ensure CRS compatibility

gt <- rasterize(rois, r_dec, field = "class")

##############################################################################
# 4.  Performance metrics
##############################################################################
library(caret)

idx     <- !is.na(values(gt))
y_true  <- values(gt)[idx]
y_pred  <- values(r_dec)[idx]

confR <- confusionMatrix(as.factor(y_pred),
                         as.factor(y_true),
                         positive = "1")

# Extract metrics
Pd   <- confR$byClass["Sensitivity"]
Pfa  <- 1 - confR$byClass["Specificity"]
F1   <- confR$byClass["F1"]
κ    <- confR$overall["Kappa"]

# ---- NEW METRICS -----------------------------------------------
Precision <- confR$byClass["Pos Pred Value"]      # PPV
OA        <- confR$overall["Accuracy"]            # Overall accuracy
# (opt) AUC will be computed in the ROC script, not here
# ---------------------------------------------------------------

# Example: MUNICH - TSALLIS
metrics_tsallis_munich_eta1_f <- data.frame(
  Scene  = "Munich",
  Method = "Tsallis_eta1",
  Pd = Pd,
  Pfa = Pfa,
  F1     = F1,
  Precision = Precision,
  Kappa  = κ,
  OA     = OA
)
save(metrics_tsallis_munich_eta1_f, file = "./Data/metrics_tsallis_munich_eta1_f.RData")

# results <- data.frame(
#   Metric = c("Pd","Pfa","Precision","F1","Kappa","OA"),
#   Value  = c(Pd, Pfa, Precision, F1, κ, OA)
# )
print(metrics_tsallis_munich_eta1_f, row.names = FALSE)

# # Save them in case you want to merge tables later
# saveRDS(results, "./Data/metrics_renyi_munich.rds")   # or metrics_shannon_…

##############################################################################
# 5.  Display results
##############################################################################
results <- data.frame(Metric = c("F1","Kappa"),
                      Value  = c(F1, κ))
print(results, row.names = FALSE)
print(confR$table)

# --- cargar función de lectura ENVI (la que ya usas) ---
source("../../../Code/R/Programs/read_ENVI_images.R")



# x <- myread.ENVI(
#   file       = "../../../Data/SAR/dublin_water_100/Intensity_HH.img",
#   headerfile = "../../../Data/SAR/dublin_water_100/Intensity_HH.hdr"
# )

x <- myread.ENVI(
  file       = "../../../Data/SAR/water_100_dublin/Intensity_HH.img",
  headerfile = "../../../Data/SAR/water_100_dublin/Intensity_HH.hdr"
)


v <- as.vector(x)
v <- v[!is.na(v) & v > 0]       

# --- ENL clásico 
mu  <- mean(v)
sig <- sd(v)
cv  <- sig / mu
ENL_classic <- 1 / (cv^2)

# --- ENL robusto (mediana & MAD) ---
med <- median(v)
mad <- mad(v, constant = 1)      # MAD sin factor 1.4826
cv_r <- mad / med
ENL_robust <- 1 / (cv_r^2)

cat("ENL clásico  :", round(ENL_classic, 2), "\n")
cat("ENL robusto  :", round(ENL_robust, 2), "\n")



# --- 1. CV basado en MnAD ---------------------------------------------------
MnAD   <- mean(abs(v - median(v)))
cv_MnAD <- MnAD / median(v)
ENL_MnAD <- 1 / (cv_MnAD^2)



# --- salida -----------------------------------------------------------------
cat(sprintf("ENL  (MnAD/median) : %.2f\n", ENL_MnAD))




# library(raster)      # para focal()
# 
# r <- raster(x)       # crea objeto raster
# w <- 20              # tamaño de ventana
# FUN_ENL <- function(z) {
#   m  <- mean(z, na.rm = TRUE)
#   s  <- sd(z,   na.rm = TRUE)
#   cv <- s / m
#   return(1 / (cv^2))
# }
# enl_map <- focal(r, w = matrix(1, w, w), fun = FUN_ENL, pad = TRUE, na.rm = TRUE)
# # enl_map es un raster con el ENL local
# plot(enl_map, main = "Local ENL (20×20)")

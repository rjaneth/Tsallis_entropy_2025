# imagematrix_visualizer_v2.R
# Unified template for plotting and exporting grayscale and color image matrices in R
# Author: Alejandro C. Frery & Janeth Alpala, extended and documented with enhancements
### 30 April 2025


## imagematrix class definition
## $Header: /database/repository/rimage/R/Attic/imagematrix.R,v 1.1.2.5 2004/03/17 06:35:18 tomo Exp $
## Copyright (c) 2003 Nikon Systems Inc.
## For complete license terms see file LICENSE


# Supported features:
# - Grayscale image matrix visualization (original base implementation)
# - Color image matrix visualization with colorbar (new implementation) 
# - Automatic detection and inversion for binary maps (0 and 1) for consistent black/white display
# - Export functions for PNG and PDF (grayscale and color)

# Usage examples:
# imagematrix(equalize(x))                            # grayscale image
# imagematrix_color(p_values)                         # color image with colorbar
# imagematrix(p_values < 0.05)                        # binary grayscale (black = TRUE)
# imagematrix_color(p_values < 0.05)                  # binary color image (black = TRUE)
# imagematrixPNG(imagematrix(p_values), "img.png")    # export grayscale
#imagematrix_colorPNG(imagematrix_color(p_values), name = "name.png",palette_option = "viridis-H",  ...)
# imagematrix_colorPDF(imagematrix_color(p_values), ...)  # export with colorbar
# previewImagematrix(imagematrix_color(...))          # preview in RStudio plot panel

# -------------------------------------------------------
# REQUIRED PACKAGES
# -------------------------------------------------------
if (!require(fields)) install.packages("fields", dependencies = TRUE)
if (!require(viridis)) install.packages("viridis", dependencies = TRUE)
if (!require(RColorBrewer)) install.packages("RColorBrewer", dependencies = TRUE)

library(fields)
library(viridis)
library(RColorBrewer)

# -------------------------------------------------------
# GRAYSCALE IMAGE FUNCTIONS 
# -------------------------------------------------------

# Creates a grayscale or RGB imagematrix object
imagematrix <- function(mat, type = NULL, ncol = dim(mat)[1], nrow = dim(mat)[2],
                        noclipping = FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
  if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
  if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
  if (type != "rgb" && type != "grey") stop("Type is incorrect.")
  if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")
  
  imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
  
  if (length(imgdim) == 3 && type == "grey") {
    mat <- rgb2grey(mat)
  }
  if (!noclipping && ((min(mat) < 0) || (max(mat) > 1))) {
    warning("Pixel values were automatically clipped because of range over.")
    mat <- clipping(mat)
  }
  mat <- array(mat, dim = imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix", class(mat))
  mat
}

# Prints basic info about imagematrix
print.imagematrix <- function(x, ...) {
  x.dim <- dim(x)
  cat("size: ", x.dim[1], "x", x.dim[2], "\n")
  cat("type: ", attr(x, "type"), "\n")
}

# Plots grayscale or RGB imagematrix
plot.imagematrix <- function(x, ...) {
  colvec <- switch(attr(x, "type"),
                   grey = {
                     u <- sort(unique(c(x)))
                     if (length(u) == 2 && all(u %in% c(0, 1))) {
                       grey(1 - x)  # Invert binary: 1 → black (TRUE), 0 → white (FALSE)
                     } else {
                       grey(x)
                     }
                   },
                   rgb = rgb(x[,,1], x[,,2], x[,,3]))
  if (is.null(colvec)) stop("image matrix is broken.")
  colors <- unique(colvec)
  colmat <- array(match(colvec, colors), dim = dim(x)[1:2])
  image(x = 0:(dim(colmat)[2]), y = 0:(dim(colmat)[1]),
        z = t(colmat[nrow(colmat):1, ]), col = colors,
        xlab = "", ylab = "", axes = FALSE, asp = 1, ...)
}

# Returns image type
imageType <- function(x) {
  attr(x, "type")
}

# Converts RGB array to grayscale
rgb2grey <- function(img, coefs = c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
  if (length(dim(img)) < 3) stop("image matrix isn't rgb image.")
  imagematrix(coefs[1] * img[,,1] + coefs[2] * img[,,2] + coefs[3] * img[,,3],
              type = "grey")
}

# Clipping values between 0 and 1
clipping <- function(img, low = 0, high = 1) {
  img[img < low] <- low
  img[img > high] <- high
  img
}

# Normalize to [0,1]
normalize <- function(img) {
  (img - min(img)) / (max(img) - min(img))
}

# Save grayscale imagematrix to PNG
imagematrixPNG <- function(x, name){
  dimensions <- dim(x)
  zero4 <- rep(0, 4)
  png(file = name, width = dimensions[2], height = dimensions[1])
  par(mar = zero4, oma = zero4, omi = zero4)
  plot.imagematrix(x)
  dev.off()
}

# Save grayscale imagematrix to EPS
imagematrixEPS <- function(x, name){
  dimensions <- dim(x)
  zero4 <- rep(0, 4)
  postscript(file = name, width = dimensions[1], height = dimensions[2], paper = "special")
  par(mar = zero4, oma = zero4, omi = zero4)
  plot(x)
  dev.off()
}

# Histogram equalization (1 band)
equalize <- function(imagem) {
  imagemeq <- ecdf(imagem)(imagem)
  dim(imagemeq) <- dim(imagem)
  return(imagemeq)
}

# Histogram equalization (3 bands independently)
equalize_indep <- function(imagem) {
  req <- ecdf(imagem[,,1])(imagem[,,1])
  geq <- ecdf(imagem[,,2])(imagem[,,2])
  beq <- ecdf(imagem[,,3])(imagem[,,3])
  imagematrix(array(c(req, geq, beq), dim = dim(imagem)))
}

# Normalize each RGB band independently
normalize_indep <- function(imagem) {
  rlin <- normalize(imagem[,,1])
  glin <- normalize(imagem[,,2])
  blin <- normalize(imagem[,,3])
  imagematrix(array(c(rlin, glin, blin), dim = dim(imagem)))
}

# -------------------------------------------------------
# COLOR IMAGE FUNCTIONS WITH COLORBAR
# -------------------------------------------------------

imagematrix_color <- function(mat, type = NULL, ncol = dim(mat)[1], nrow = dim(mat)[2],
                              noclipping = FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
  if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
  if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
  if (type != "rgb" && type != "grey") stop("Type is incorrect.")
  if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")
  
  imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
  if (length(imgdim) == 3 && type == "grey") {
    mat <- rgb2grey(mat)
  }
  if (!noclipping && ((min(mat, na.rm=TRUE) < 0) || (1 < max(mat, na.rm=TRUE)))) {
    warning("Pixel values were automatically clipped because of range over.")
    mat <- pmax(0, pmin(1, mat))
  }
  mat <- array(mat, dim = imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix_color", class(mat))
  mat
}

# [plot.imagematrix_color, previewImagematrix, imagematrix_colorPNG, imagematrix_colorPDF] already include binary detection, so no change needed here.

# Función para obtener el tipo de imagen
imageType <- function(x) {
  attr(x, "type")
}

# Conversión de RGB a escala de grises
rgb2grey <- function(img, coefs = c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
  if (length(dim(img)) < 3) stop("image matrix isn't rgb image.")
  imagematrix_color(drop(img[,,1] * coefs[1] + img[,,2] * coefs[2] + img[,,3] * coefs[3]), type = "grey")
}


.get_palette <- function(ncolors, palette_option, direction = 1) {
  if (grepl("viridis-", palette_option)) {
    opt <- strsplit(palette_option, "-")[[1]][2]
    viridis(ncolors, option = opt, direction = direction)
  } else {
    pal <- brewer.pal(min(ncolors, 9), palette_option)
    if (direction == -1) pal <- rev(pal)
    pal
  }
}


plot.imagematrix_color <- function(x, 
                                   significance_level = 0.05, 
                                   ncolors           = 100,
                                   palette_option    = "viridis-H",
                                   legend_shrink     = 0.9,
                                   direction         = 1,
                                   ...) {
  vals <- unique(as.vector(x))
  if (length(vals) == 2 && all(sort(vals) %in% c(0, 1))) x <- 1 - x
  
  palette_colors <- .get_palette(ncolors, palette_option, direction)
  
  zlim <- c(0, 1)
  base_breaks <- seq(0, 1, by = 0.2)
  all_breaks  <- sort(unique(c(base_breaks, significance_level)))
  all_labels  <- as.character(all_breaks)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(7, 1))  
  par(mar = c(0, 0, 0, 0))
  image(x = 1:ncol(x), y = 1:nrow(x),
        z = t(x[nrow(x):1, , drop=FALSE]),
        col = palette_colors, zlim = zlim,
        axes = FALSE, xlab = "", ylab = "", asp = 1, ...)
  
  par(mar = c(1, 0, 1, 0))  
  image.plot(zlim = zlim, legend.only = TRUE,
             col = palette_colors, horizontal = FALSE,
             axis.args = list(at = all_breaks, labels = all_labels, cex.axis = 2.4))
}





# Apply same fix to previewImagematrix

previewImagematrix <- function(x,
                               significance_level = 0.05,
                               ncolors           = 100,
                               palette_option    = "viridis-H",
                               direction         = 1, ...) {
  vals <- unique(as.vector(x))
  if (length(vals) == 2 && all(sort(vals) %in% c(0, 1))) x <- 1 - x
  
  palette_colors <- .get_palette(ncolors, palette_option, direction)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(7, 1))
  par(mar = c(0, 0, 0, 6))
  image(1:ncol(x), 1:nrow(x),
        t(x[nrow(x):1, , drop = FALSE]),
        col = palette_colors, zlim = c(0,1),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  
  par(mar = c(3, 3, 3, 1))
  ticks <- sort(unique(c(seq(0,1,0.2), significance_level)))
  image.plot(zlim = c(0,1), legend.only = TRUE, col = palette_colors,
             horizontal = FALSE,
             axis.args = list(at = ticks, labels = as.character(ticks), cex.axis = 1.2))
}



previewImagematrixPanel <- function(x,
                                    significance_level = 0.05,
                                    ncolors            = 100,
                                    palette_option     = "viridis-H",
                                    legend_frac        = 0.9,
                                    main               = NULL,
                                    direction          = 1) {
  if (identical(sort(unique(as.vector(x))), c(0,1))) x <- 1 - x
  pal   <- .get_palette(ncolors, palette_option, direction)
  ticks <- sort(unique(c(seq(0,1,.2), significance_level)))
  
  whole <- par("plt"); w <- whole[2]-whole[1]
  img_r <- c(whole[1], whole[1]+w*(1-legend_frac), whole[3], whole[4])
  leg_r <- c(whole[1]+w*(1-legend_frac), whole[2], whole[3], whole[4])
  
  par(plt = img_r, new = FALSE, mar = c(0,0,0,0))
  image(1:ncol(x), 1:nrow(x), t(x[nrow(x):1,, drop=FALSE]),
        col = pal, zlim = c(0,1), axes = FALSE, xlab = "", ylab = "", useRaster = TRUE, asp = 1)
  if (!is.null(main)) mtext(main, side = 3, line = 0.8, cex = 1)
  
  par(plt = leg_r, new = TRUE, mar = c(0,0,0,0))
  fields::image.plot(zlim = c(0,1), legend.only = TRUE, col = pal,
                     horizontal = FALSE,
                     smallplot = c(0.92, 0.95, 0.15, 0.95),
                     axis.args = list(at = ticks, labels = format(ticks), cex.axis = 0.8))
  par(new = FALSE, plt = whole)
}





imagematrix_colorPNG <- function(x, name,
                                 significance_level = 0.05,
                                 ncolors           = 100,
                                 palette_option    = "viridis-H",
                                 legend_width_px   = 100,
                                 scale_factor      = 1,
                                 direction         = 1,
                                 ...) {
  vals <- unique(as.vector(x))
  if (length(vals) == 2 && all(sort(vals) %in% c(0, 1))) x <- 1 - x
  
  img_width  <- ncol(x) * scale_factor
  img_height <- nrow(x) * scale_factor
  total_width <- img_width + legend_width_px
  
  png(file = name, width = total_width, height = img_height, units = "px", res = 200)
  layout(matrix(c(1,2), nrow = 1), widths = c(img_width, legend_width_px))
  par(oma = c(0,0,0,0), mar = c(0,0,0,0))
  
  palette_colors <- .get_palette(ncolors, palette_option, direction)
  
  image(1:ncol(x), 1:nrow(x), t(x[nrow(x):1,, drop=FALSE]),
        col = palette_colors, zlim = c(0,1),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE, ...)
  
  par(oma = c(0,0,0,0), mar = c(0,0,0,0))
  ticks  <- sort(unique(c(seq(0,1,0.2), significance_level)))
  labels <- ifelse(abs(ticks - significance_level) < 1e-8, "0.05", as.character(ticks))
  image.plot(zlim = c(0,1), legend.only = TRUE, col = palette_colors,
             horizontal = FALSE,
             smallplot = c(0.2, 0.53,  0.02, 0.98),
             axis.args = list(at = ticks, labels = labels, cex.axis = 3.1))
  dev.off()
}

  #   #mexico
#   image.plot(zlim = c(0, 1),
#              legend.only = TRUE,
#              col = palette_colors,
#              horizontal = FALSE,
#              smallplot = c(0.3, 0.5,  0.02, 0.98),#smallplot = c(0.30, 0.58, 0.1, 0.9),
#              axis.args = list(
#                at     = ticks,
#                labels = labels,
#                #at     = sort(unique(c(seq(0, 1, by = 0.2), significance_level))),
#                #labels = as.character(sort(unique(c(seq(0, 1, by = 0.2), significance_level)))),
#                cex.axis = 1.5 # 3.3 dublin 3.3, 3.1 munich, 2.1 mexico
#              ))
#   dev.off()
# }

# #london
#   image.plot(zlim = c(0, 1),
#              legend.only = TRUE,
#              col = palette_colors,
#              horizontal = FALSE,
#              smallplot = c(0.2, 0.54,  0.02, 0.98),#smallplot = c(0.30, 0.58, 0.1, 0.9),
#              axis.args = list(
#                at     = ticks,
#                labels = labels,
#                #at     = sort(unique(c(seq(0, 1, by = 0.2), significance_level))),
#                #labels = as.character(sort(unique(c(seq(0, 1, by = 0.2), significance_level)))),
#                cex.axis = 6 # london
#              ))
#   dev.off()
# }
#   
#   #sf
#   image.plot(zlim = c(0, 1),
#              legend.only = TRUE,
#              col = palette_colors,
#              horizontal = FALSE,
#              smallplot = c(0.2, 0.54,  0.02, 0.98),#smallplot = c(0.30, 0.58, 0.1, 0.9),
#              axis.args = list(
#                at     = ticks,
#                labels = labels,
#                #at     = sort(unique(c(seq(0, 1, by = 0.2), significance_level))),
#                #labels = as.character(sort(unique(c(seq(0, 1, by = 0.2), significance_level)))),
#                cex.axis = 6.3 # london
#              ))
#   dev.off()
# }


imagematrix_colorPDF <- function(x, name,
                                 significance_level = 0.05,
                                 ncolors           = 100,
                                 palette_option    = "viridis-H",
                                 legend_width_in   = 2.8,
                                 scale_factor      = 25,
                                 direction         = 1,
                                 ...) {
  vals <- unique(as.vector(x))
  if (length(vals) == 2 && all(sort(vals) %in% c(0, 1))) x <- 1 - x
  
  img_width_in  <- ncol(x) / scale_factor
  img_height_in <- nrow(x) / scale_factor
  total_width_in <- img_width_in + legend_width_in
  
  pdf(file = name, width = total_width_in, height = img_height_in, useDingbats = FALSE)
  layout(matrix(c(1,2), nrow = 1), widths = c(img_width_in, legend_width_in))
  par(oma = c(0,0,0,0), mar = c(0,0,0,0))
  
  palette_colors <- .get_palette(ncolors, palette_option, direction)
  
  image(1:ncol(x), 1:nrow(x), t(x[nrow(x):1,, drop=FALSE]),
        col = palette_colors, zlim = c(0,1),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE, ...)
  
  par(oma = c(0,0,0,0), mar = c(0,0,0,0))
  ticks <- sort(unique(c(seq(0,1,0.2), significance_level)))
  image.plot(zlim = c(0,1), legend.only = TRUE, col = palette_colors,
             horizontal = FALSE,
             smallplot = c(0.25, 0.45, 0.1, 0.9),
             axis.args = list(at = ticks, labels = as.character(ticks), cex.axis = 2.9))
  dev.off()
}

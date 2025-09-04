# -------------------------------------------------------
# Si no los tienes instalados:
# install.packages("fields")
# install.packages("viridis")

library(fields)
library(viridis)

# -------------------------------------------------------
# 1) Definición de la clase 'imagematrix' y auxiliares

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
  
  if (!noclipping && ((min(mat, na.rm=TRUE) < 0) || (1 < max(mat, na.rm=TRUE)))) {
    warning("Pixel values were automatically clipped because of range over.")
    mat <- clipping(mat)
  }
  
  mat <- array(mat, dim = imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix", class(mat))
  mat
}

# Retorna el tipo de la imagematrix ("rgb" o "grey")
imageType <- function(x) {
  attr(x, "type")
}

rgb2grey <- function(img, coefs=c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
  if (length(dim(img)) < 3) stop("image matrix isn't rgb image.")
  imagematrix(coefs[1]*img[,,1] + coefs[2]*img[,,2] + coefs[3]*img[,,3],
              type="grey")
}

clipping <- function(img, low=0, high=1) {
  img[img < low]  <- low
  img[img > high] <- high
  img
}

normalize <- function(img) {
  (img - min(img)) / (max(img) - min(img))
}

# -------------------------------------------------------
# 2) plot.imagematrix: sin márgenes, con barra a la derecha,
#    y posibilidad de "encoger" la barra con legend_shrink

plot.imagematrix <- function(x, 
                             significance_level = 0.05, 
                             ncolors           = 200,
                             palette_colors    = viridis(ncolors, option = "D"), 
                             legend_shrink     = 0.9,
                             ...) {
  
  # Rango de valores
  zlim <- c(0, 1)
  
  # Creamos los breaks y labels que incluyan significance_level
  base_breaks <- seq(0, 1, by = 0.2)
  all_breaks  <- sort(unique(c(base_breaks, significance_level)))
  all_labels  <- as.character(all_breaks)
  
  # Layout a 2 columnas: la imagen (ancho 5) y la barra (ancho 1)
  layout(matrix(c(1, 2), nrow = 1), widths = c(7, 1)) 
  # ↑ Si ves que sigue muy pegada, prueba widths=c(6,1) u otro valor mayor.
  
  # Panel 1: imagen sin márgenes
  par(mar = c(0, 0, 0, 0))
  image(x = 1:ncol(x),
        y = 1:nrow(x),
        z = t(x[nrow(x):1, , drop=FALSE]),
        col    = palette_colors,
        zlim   = zlim,
        axes   = FALSE,  
        xlab   = "", 
        ylab   = "",
        asp    = 1,      
        ...)
  
  # Panel 2: barra de color, ajustando márgenes
  par(mar = c(1, 0, 1, 0))  
  # ↑ Ajusta según la separación que necesites
  
  image.plot(
    zlim         = zlim,
    legend.only  = TRUE,
    col          = palette_colors,
    legend.width = 5,
    horizontal   = FALSE,
    axis.args    = list(
      at     = all_breaks,
      labels = all_labels,
      cex.axis = 2.4  # ↑ Tamaño de las etiquetas numéricas de la barra
    ),
    legend.shrink = legend_shrink
  )
}


# -------------------------------------------------------
# 3) Funciones para GUARDAR la imagematrix en PNG/EPS
#    usando exactamente el plot.imagematrix anterior.

### Added by Alejandro C. Frery
### 24 April 2014
imagematrixPNG <- function(x, name,
                           significance_level = 0.05,
                           ncolors           = 200,
                           palette_colors    = viridis(ncolors, option = "D"),
                           legend_shrink     = 0.9,
                           extra_width       = 200,  # espacio extra p/ la barra
                           ...) {
  
  # Dimensiones base: ancho = ncol(x), alto = nrow(x)
  w <- ncol(x)
  h <- nrow(x)
  
  # Sumamos 'extra_width' para dar espacio a la barra de color
  png(file = name, width = w + extra_width, height = h)
  
  # Llamamos a plot.imagematrix con los mismos parámetros
  plot(x, 
       significance_level = significance_level,
       ncolors            = ncolors,
       palette_colors     = palette_colors,
       legend_shrink      = legend_shrink,
       ...)
  
  dev.off()
}

imagematrixEPS <- function(x, name,
                           significance_level = 0.05,
                           ncolors           = 200,
                           palette_colors    = viridis(ncolors, option = "D"),
                           legend_shrink     = 0.9,
                           width_inch        = 7,
                           height_inch       = 5,
                           ...) {
  
  postscript(file = name, 
             width = width_inch, 
             height = height_inch,
             paper = "special", 
             horizontal = FALSE)
  
  plot(x, 
       significance_level = significance_level,
       ncolors            = ncolors,
       palette_colors     = palette_colors,
       legend_shrink      = legend_shrink,
       ...)
  
  dev.off()
}

# (Si quieres también un PDF, puedes crear imagematrixPDF análogo.)
# -------------------------------------------------------
# 4) Funciones de ecualización (por si las usas):
equalize <- function(imagem) {
  imagemeq <- ecdf(imagem)(imagem)
  dim(imagemeq) <- dim(imagem)
  return(imagemeq)
}
equalize_indep <- function(imagem) {
  req <- ecdf(imagem[,,1])(imagem[,,1])
  geq <- ecdf(imagem[,,2])(imagem[,,2])
  beq <- ecdf(imagem[,,3])(imagem[,,3])
  
  imagematrix(
    array(c(req, geq, beq), dim=dim(imagem)),
    type = imageType(imagem)
  )
}
normalize_indep <- function(imagem) {
  rlin <- normalize(imagem[,,1])
  glin <- normalize(imagem[,,2])
  blin <- normalize(imagem[,,3])
  
  imagematrix(
    array(c(rlin, glin, blin), dim=dim(imagem)),
    type = imageType(imagem)
  )
}

# -------------------------------------------------------
# EJEMPLO DE USO (descomenta para probar):
#
# # 1) Generamos matriz de p-valores
# p_vals <- matrix(runif(600*800, 0, 1), 600, 800)
#
# # 2) Convertimos a imagematrix
# im_pvals <- imagematrix(p_vals, type="grey")
#
# # 3) Graficamos en pantalla
# plot(im_pvals, significance_level=0.05, legend_shrink=0.7)
#
# # 4) Guardamos a PNG (espacio extra 300 px p/ la barra)
# imagematrixPNG(im_pvals, name="example.png",
#                significance_level=0.01,
#                legend_shrink=0.5,
#                extra_width=300)
#
# # 5) Ecualizar, y luego guardar
# eq_im <- imagematrix(equalize(p_vals), type="grey")
# imagematrixPNG(eq_im, name="example_eq.png")

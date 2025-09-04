# plot SAR images

source("read_ENVI_images.R")

x <- myread.ENVI(file='../Data-SAR/L16_envi_dublin_size_1100/Intensity_HH.img', 
                 headerfile='../Data-SAR/L16_envi_dublin_size_1100/Intensity_HH.hdr')#dublin

source("imagematrix_visualizer_v2.R")

#plot(imagematrix(p_values))
#plot(imagematrix(equalize(difference_values)))

plot(imagematrix(equalize(x))) # visualize the SAR image

plot(imagematrix(p_values< 0.05)) # binary map


# save images in grayscale

# imagematrixPNG(imagematrix(equalize(x)), name = "image1.png")
# imagematrixPNG(imagematrix(equalize(difference_values)), name = "test.png")
#imagematrixPNG(imagematrix(p_values), name="H_pvalue.png")
#imagematrixPNG(imagematrix(p_values<0.05), name="H_005_.png")


# save images with color viridis-H gradient

imagematrix_colorPNG(
  imagematrix_color(p_values),
  name = "p-values_image.png",
  palette_option = "viridis-H",
  legend_width_px = 540,
  scale_factor = 2,
  direction = -1
)




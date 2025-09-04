# plot abrir imagen sar

source("read_ENVI_images.R")

x <- myread.ENVI(file='./SAR/envi_dublin_1100_HH/Intensity_HH.img', 
                 headerfile='./SAR/envi_dublin_1100_HH/Intensity_HH.hdr')#dublin

source("../imagematrix.R")
#hist(x)

#plot(imagematrix(difference_values))

plot(imagematrix(equalize(x)))

# 
#imagematrixPNG(imagematrix(equalize(x)), name = "dublin.png")

#imagematrixPNG(imagematrix(p_values>0.05), name="H_005_p_values.png")

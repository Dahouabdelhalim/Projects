##################################################################
# 
# R functions used to create artificial noise in the images
# in the "Noisy OCR Dataset"
#
# Thomas Hegghammer
#
# June 2021 
#
# This file accompanies the "Noisy OCR Dataset" Zenodo repository.
# 
##################################################################

library(glue)
library(magick)
library(shape)
library(showtext)
library(yarrr)
library(dplyr)

## Notes: 

# a) The "Watermark", "Scribbles", and "Ink stains" functions 
# require the image dimensions to be set

img_width <- #<TBC>
img_height <- #<TBC>

# b) The "Scribbles" function requires the proprietary font "Shopping Script" 
# (https://www.dafont.com/shopping-script.font)


## 1. BLUR

blur <- function(source_ims, dest_folder) {
  for (i in source_ims){
    image_read(i) %>%
      image_blur(5, 4) %>%
      image_write(glue("{dest_folder}/{basename(i)}"))
  }
}

## 2. WEAK INK

weaken <- function(source_ims, dest_folder) {
  for (i in source_ims){
    image_read(i) %>%
      image_oilpaint() %>%
      image_write(glue("{dest_folder}/{basename(i)}"))
  }
}

## 3. SALT & PEPPER

snp <- function(source_ims, dest_folder) {
  for (i in source_ims){
    image_read(i) %>%
      image_noise(noisetype = "poisson") %>%
      image_write(glue("{dest_folder}/{basename(i)}"))
  }
}

## 4. WATERMARK

watermark <- function(source_ims, dest_folder) {
  transp_grey <- transparent(orig.col = "gray80", trans.val = 0.4, maxColorValue = 255)
  for (i in source_ims){
    img <- image_read(i)
    tiff(file=glue("{dest_folder}/{basename(i)}"), width=img_width, height=img_height, units="px", res=300)
    plot(img)
    text(1200, 2000, "Watermark", cex=20, srt=50, col=transp_grey)
    dev.off()
  }
}

## 5. SCRIBBLES

scribble <- function(source_ims, dest_folder) {
  font_add(family = "Shopping Script", regular = "~/.fonts/ShoppingScript-Regular.otf")
  showtext_auto()
  for (i in source_ims){
    img <- image_read(i)
    tiff(file=glue("{dest_folder}/{basename(i)}"), width=img_width, height=img_height, units="px", res=300)
    plot(img)
    text(1500, 2000, "fascinating", family = "Shopping Script", cex = 10, srt = 15, col = "gray30")
    text(600, 1000, "__", family = "Shopping Script", cex = 12, col = "gray10", srt = 4)
    text(300, 700, "_ _", family = "Shopping Script", cex = 12, col = "gray20", srt = -4)
    text(1600, 800, "NB!", family = "Shopping Script", cex = 10, srt = -10, col = "gray30")
    text(1650, 1400, "V", family = "Shopping Script", cex =8, srt = 10, col = "gray35")
    text(1700, 500, "V", family = "Shopping Script", cex =12, srt = 30, col = "gray25")
    text(400, 1500, "mmm", family = "Shopping Script", cex=10)
    text(1100, 500, "mnnm", family = "Shopping Script", cex=10, col = "gray35")
    text(1000, 1200, "mmmmmmm", family = "Shopping Script", cex=10, col = "gray25")
    text(50, 1300, "|", family = "Shopping Script", cex=10)
    text(100, 1100, "Z", family = "Shopping Script", cex=10, col = "gray25")
    text(400, 1800, "O", family = "Shopping Script", cex=12, col = "gray35")
    text(1200, 700, "O", family = "Shopping Script", cex=10, col = "gray25")
    text(450, 400, "X", family = "Shopping Script", cex=10, col = "gray15")
    dev.off()
  }
}

## 6. INK STAINS

ink <- function(source_ims, dest_folder) {
  for (i in source_ims){
    tiff(file="ink.tiff", width=img_width, height=img_height, units="px", res=300)
    par(bg=NA)
    emptyplot()
    filledshape(matrix(nc = 4, nr = 4, runif(8)), col = shadepalette(50, "black", "grey20"))
    dev.off()
    main <- image_read(i)
    inset <- image_read("ink.tiff")
    final <- image_composite(main, inset, operator = "atop")
    image_write(final, glue("{dest_folder}/{basename(i)}"))
  }
}

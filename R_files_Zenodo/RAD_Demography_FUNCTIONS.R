
# Plot theme --------------------------------------------------------------

theme_claudia <- function(addlegend = FALSE)
{
  theme_bw()+
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          panel.spacing = unit(0.8, "lines"),
          axis.title.x = element_text(size = 10, colour = "black", vjust = 0),
          axis.title.y = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.5),
          axis.line = element_blank(),
          legend.title = element_blank(),
          legend.text = if(addlegend == TRUE) element_text(size = 9) else element_blank(),
          legend.background = if(addlegend == TRUE) element_rect(fill = "transparent") else element_blank(),
          legend.key.size = if(addlegend == TRUE) unit(0.7, "line") else element_blank())
}



# Plot colours ------------------------------------------------------------

yellow <- "#ffd32a"
dark_blue <- "#2980b9"
dark_green <- "#0A6B28"
purple <- "#833471"
fire_orange <- "#EE5A24"
sky_blue <- "#87CEEB"
grey <- "#57595D"
light_green <- "#2ecc71"
pink <- "#ffc0cb"
mid_green <- "#008000"
navy <- "#123456"
pale <- "#eee8aa"
dark_red <- "#8b0000"
bright_orange <- "#ff8c00"


# Relatedness GCTA --------------------------------------------------------
# read the GRM output from GCTA relatedness 

sum_i=function(i) {
  return(sum(1:i))
}

GRMreader <- function(filenm, flag=1)
{
  xDt.bin <- paste(filenm, ".grm.bin", sep="")
  xDt.nfl <- paste(filenm, ".grm.N.bin", sep="")
  xDt.gid <- paste(filenm, ".grm.id", sep="")
  
  xDt.id <- read.table(xDt.gid)
  xDt.n <- dim(xDt.id)[1]
  xDt.grm <- readBin(file(xDt.bin, "rb"), n=xDt.n*(xDt.n+1)/2, what=numeric(0), size=4)
  
  sn <- sapply(1:xDt.n, sum_i)
  off <- xDt.grm[-sn]
  diag <- xDt.grm[sn]
  if(flag==1) return(list(diag=diag, off=off, n=xDt.n))
  else {
    xDt.mat <- matrix(data=NA, nrow=xDt.n, ncol=xDt.n)
    xDt.mat[upper.tri(xDt.mat)] <- off
    xDt.mat <- t(xDt.mat)
    xDt.mat[upper.tri(xDt.mat)] <- off
    diag(xDt.mat) <- diag
    return(list(mat=xDt.mat, n=xDt.n))
  }
}


# TreeMix -----------------------------------------------------------------

# you need to download the file "plotting_funcs.R" from the TreeMix package release.
source("plotting_funcs.R")
# note. There is no data editing for main TreeMix plots so skip straight to "DO" for these plots. 


# Ruminant functions
# Defining functions that I'll be using throughout the data processing, analysis, and visualization

# loading packages
#############
 library(colorblindcheck)
 library(Rvcg)
 library(vioplot)
 library(geomorph)
 library(Morpho)
 library(stringr)
 library(ape)
 library(phytools)
 library(colorRamps)
 library(MASS)
 library(spatial)
 library(phangorn)
 library(evolqg)
 library(geiger)
 library(nlme)
 library(adephylo)
 library(mvMORPH)
 library(sp)
 library(ggtern)
 library(mixtools)
 library(colorBlindness)
 library(inlmisc)
 library(pracma)
#############


#TPS function code...
###########
# clean TPS code

  # target.lm - the warped landmark configuration
  # reference.lm - the reference landmark configuration
  # scale =  1 - how amplified the landmarks should be relative to the grid
      # this scales the entire landmark configuration, to modify the x and y aspects separately use 'buffer.x' and 'buffer.y'
  # n.grid.col = 25 - the number of horizontal grid vertices
  # n.grid.row = NULL - the number of vertical grid vertices, automaticaly calculated if left NULL
  # palette = NULL - the palette to color the cells
  # shade.method = c("area", "shape.proc", "shape.natural")
    # if both 'palette' and 'shade.method' are empty -> no shade
    # if a palette is defined but no method is specified, the 'area' method is used where cells are colored according to how different in area they are to the unwarped grid
  # shade.transparency = 1 - transparency of the shade, where 1 is not transparent and 0 is invisible
  # abs = FALSE - whether or not you want larger/smaller cells colored similarly -> only works with 'area' shade method
  # grid.color = 'black' -  the color of the grid
  # grid.width = 1 - width of the grid lines
  # links - index of landmarks to be connected with lines
  # link.color = 'black' - color of the links
  # link.width = 2 - width of the links
  # link.transparency = 1 - transparency of the links, where 1 is not transparent and 0 is invisible
  # lm.size = 1 - size of the landmarks
  # lm.color = 'black' - color of the landmarks
  # lm.transparency = 1 - transparency of the landmarks
  # plot.lm = TRUE - logical whether or not to plot the landmarks

  # smooth = FALSE - if true, a grid is plotted over an oversampled heatmap... -> DO THIS!
  # add landmark groups...
  # make it a movie....

  # option to 'snap' the grid to the bounding box of the landmarks OR to the actual landmark configuration
      # probably by doing the bounding box first, and then pruning the grid cells outside of the convex hull of the landmarks
      #
  # option to color landmarks according to variance from reference
  # should make somthing to denote if you should also plot the reference LM...
  #
  # add arrows
    # option to override the 'joinline' thing -> automatically a convex hull of the LM configuration
  # contour plot option
  # Claude's shade method, Zelditch's shade method
  # add 'trap' for when the user tries to define a ton of grid cells, and provide option for them to go back and change the number of grid cells...

  # make it possible to add a bunch of landmark configurations at once, and only plot one, but with the shade colored according to the total range of change in area of cells..
    # aka... the scale is the same across all the shape models

tps <- function(target.lm,
                reference.lm,
                A=NULL,
                shade.scale=FALSE,
                scale.lm = 1, # what to scale the whole landmark configuration by
                mag = 1,
                add = FALSE,
                at = c(0,0), # where to add the TPS grids, if anywhere
                snap = TRUE,
                #snap.buffer = c(0,0,0,0) , # how many rows or columns to add onto the snap as a buffer... left, right, down, up (xmin, xmax, ymin, ymax)
                shade = FALSE,
                shade.trans = 0.7,
                palette = NULL,
                n.grid.col = 24,
                n.grid.row = NULL,
                grid.aes = list(col = "black", lwd = 1, trans = 1),
                links = NULL,
                link.aes = list(col = "black", lwd = 3, trans = 1),
                plot.lm = TRUE,
                lm.aes = list(col = "black", cex = 1, trans = 1),
                plot.ref.links = FALSE,
                ref.link.aes = list(col = "gray", lwd = 3, trans = 0.85),
                plot.ref.lm = FALSE,
                ref.lm.aes = list(col = "gray", cex = 1, trans = 0.5),
                plot.arrows = FALSE,
                arrow.aes = list(col = "black", lwd = 1, trans = 1, code = 2, length = 0.1),
             #   smooth = FALSE,
                buffer.x = 0.3,
                buffer.y = 0.3){

  # helper functions
#####
# from 'Morphometrics with R' by Claude
# computes deformation of grid
  tps2d<-function(M, matr, matt){p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
  P<-matrix(NA, p, p)
  for (i in 1:p)
     {for (j in 1:p){
         r2<-sum((matr[i,]-matr[j,])^2)
         P[i,j]<- r2*log(r2)}}
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
  m2<-rbind(matt, matrix(0, 3, 2))
  coefx<-solve(L)%*%m2[,1]
  coefy<-solve(L)%*%m2[,2]
  fx<-function(matr, M, coef)
     {Xn<-numeric(q)
      for (i in 1:q)
           {Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
           Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
      Xn}
  matg<-matrix(NA, q, 2)
  matg[,1]<-fx(matr, M, coefx)
  matg[,2]<-fx(matr, M, coefy)
  matg}

  # my own internal function to make an array of the cells
  # cell -> array == celery
celery <- function(X){
  for(i in 1:n.row){
  rows <- array(dim=c(4,2,n.col))
  for(j in 1:n.col){
  rows[1,,j] <- X[1+(j-1)+((i-1)*n),]
  rows[2,,j] <- X[1+(j-1)+n+((i-1)*n),]
  rows[3,,j] <- X[1+j+n+((i-1)*n),]
  rows[4,,j] <- X[1+j+((i-1)*n),]
  }
if(i==1){cell.array<-rows}else{
cell.array <- bindArr(cell.array,rows,along=3)}
  }
  return(cell.array)
}
#####


# scaling the landmarks
  reference.lm <- reference.lm*scale.lm

    # rotating the target landmark configuration onto the reference
target.lm <- rotonto(x = reference.lm, y = target.lm, scale = TRUE)$yrot

  # incorporating the magnitude of the shape deformation
  # by modifying the 'target.lm' object
  if(mag != 1){
    target.lm <- reference.lm+((target.lm-reference.lm)*mag)
  }


  #
  if(!add){at <- c(0,0)}
  n <- n.grid.col
  s <- 1/1#zoom
  matr <- translate.lm(reference.lm, at)
  matt <- translate.lm(target.lm, at)
  M2 <- matt
  xm<-min(matt[,1])*s # x min
  ym<-min(matt[,2])*s # x max
  xM<-max(matt[,1])*s # y min
  yM<-max(matt[,2])*s # y max
  rX<-xM-xm; rY<-yM-ym # ranges


  pad.x <- buffer.x*rX
  pad.y <- buffer.y*rY

  a<-seq(xm-pad.x, xM+pad.x, length=n) # x values of grid points

  # if the number of rows isn't defined, this calculates it
  # if it is, this assigns n.grid.row to object 'm'
  if(is.null(n.grid.row)){m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))}else{(m <- n.grid.row)}


  b<-seq(ym-pad.y, yM+pad.y,length=m) # y values of grid points
  M<-as.matrix(expand.grid(a,b))

  ngrid<-tps2d(M,matr,matt)


  n.col <- (n-1)
  n.row <- length(b)-1
  n.cell <- n.col*n.row


# warped array of cells
  # storing the grid-points that define each cell in an array
cell.array <- celery(X=ngrid)

# reference cells
cell.array.ref <- celery(X=M)

# calculating the individual area of a cell in the reference (unwarped) configuration
# for use later, but it feels right to compute it here
area.reference <- abs(polyarea(cell.array.ref[,1,1],cell.array.ref[,2,1]))




# if 'snap' is true, this section trims the grid to only fit cells with landmarks inside
  # probably a much faster way of doing this
if(snap){
xmr <- min(matr[,1])
ymr <- min(matr[,2])
xMr <- max(matr[,1])
yMr <- max(matr[,2])

grid.id <- expand.grid(1:n.col,1:n.row) # making row and columns id's of each grid cell

xm.idr <- which(matr[,1] == xmr)
ym.idr <- which(matr[,2] == ymr)
xM.idr <- which(matr[,1] == xMr)
yM.idr <- which(matr[,2] == yMr)

dists.xmr <- c()
dists.ymr <- c()
dists.xMr<- c()
dists.yMr <- c()
for(i in 1:n.cell){
  dists.xmr[[i]] <- euc.dist(matr[xm.idr,],centroid1(cell.array.ref[3:4,,i]))
  dists.ymr[[i]] <- euc.dist(matr[ym.idr,],centroid1(cell.array.ref[2:3,,i]))
  dists.xMr[[i]] <- euc.dist(matr[xM.idr,],centroid1(cell.array.ref[1:2,,i]))
  dists.yMr[[i]] <- euc.dist(matr[yM.idr,],centroid1(cell.array.ref[c(1,4),,i]))
  }
which.xmr <- which.min(unlist(dists.xmr))
which.ymr <- which.min(unlist(dists.ymr))
which.xMr <- which.min(unlist(dists.xMr))
which.yMr <- which.min(unlist(dists.yMr))

col.min.idr <- grid.id[which.xmr,1]#-snap.buffer[1]
col.max.idr <- grid.id[which.xMr,1]#+snap.buffer[2]
row.min.idr <- grid.id[which.ymr,2]#-snap.buffer[3]
row.max.idr <- grid.id[which.yMr,2]#+snap.buffer[4]

min.cell <- cell.array.ref[1,,which(grid.id[,1]==col.min.idr & grid.id[,2]==row.min.idr)]
max.cell <- cell.array.ref[3,,which(grid.id[,1]==col.max.idr & grid.id[,2]==row.max.idr)]

grid.to.drop <- which(M[,1] < min.cell[1] | M[,1] > max.cell[1] | M[,2] < min.cell[2] | M[,2] > max.cell[2])
cells.to.drop <- which(grid.id[,1] < col.min.idr | grid.id[,1] > col.max.idr | grid.id[,2] < row.min.idr | grid.id[,2] > row.max.idr)

cell.array <- cell.array[,,-cells.to.drop]
n.cell <- n.cell-length(cells.to.drop)

ngrid <- ngrid[-grid.to.drop,]
n <- length(col.min.idr:col.max.idr)+1
m <- length(row.min.idr:row.max.idr)+1
}







# if there's an array:
  # need to consolidate this into a couple of clean functions.
if(!is.null(A)){
  if(is.matrix(A)){arrayspecs(A=A,p=dim(target.lm)[1],k=dim(target.lm)[2])} # if a matrix is provided, this switches it to an array

# landmark variation
target.array <- A
all.diff <- c()

for(i in 1:dim(target.array)[3])
  {

  target.lm1 <- rotonto(x = reference.lm, y = target.array[,,i], scale = TRUE)$yrot

  # incorporating the magnitude of the shape deformation
  # by modifying the 'target.lm' object
  if(mag != 1){
    target.lm1 <- reference.lm+((target.lm1-reference.lm)*mag)
  }

  #
  if(!add){at <- c(0,0)}
  matt1 <- target.lm1*scale.lm + at

  ngrid1<-tps2d(M,matr,matt1)


  # warped array of cells
  # storing the gridpoints that define each cell in an array
cell.array1 <- celery(X=ngrid1)

cell.array1 <- cell.array1[,,-cells.to.drop]

diff1 <- c()
for(i in 1:n.cell){
  area.target1 <- abs(polyarea(cell.array1[,1,i],cell.array1[,2,i]))
  diff1[[i]] <- area.target1-area.reference
  }
diff1 <- as.vector(unlist(diff1))


if(i == 1){all.diff <- diff1}
if(i > 1){all.diff <- c(all.diff, diff1)}
}

cs <- range(unlist(all.diff)) # the color scale, taking into account all the variation in the landmark data provided

}
if(is.null(A)){}




# setting the color palette according to area of the cells
if(!is.null(palette) | shade){

  # setting the pallette if none is given:
  if(is.null(palette)){shade.cols <- colorRampPalette(Blue2DarkRed18Steps)(n.cell)}else{shade.cols <- palette(n.cell)}
  idcolor <- addTrans(shade.cols, 255*shade.trans)

diff <- numeric(n.cell)
for(i in 1:n.cell){
  area.target <- abs(pracma::polyarea(cell.array[,1,i],cell.array[,2,i]))
    diff[i] <- (area.target-area.reference)
  }
#if(ABS){diff <- abs(rescale.numeric(diff,to=c(-1,1)))} # if you want small cells to be the same color as larger cells


if(shade.scale){nm <- max(abs(cs))*1.01 ; brk <- seq(-nm,nm,length.out = n.cell)}else{
  nm <- max(abs(range(diff)))*1.01 ; brk <- seq(-nm,nm,length.out = n.cell)
}


dit2<-cut(diff,breaks = length(brk))
}


# plotting:

# if you aren't adding to an existing plot, this just plots a lil canvas
if(!add){
plot(NA, asp=1,axes=F,xlab="",ylab="",xlim = range(ngrid[,1]), ylim = range(ngrid[,2]))
  }

# move this to be below everything else in the plotting section !
if(!is.null(palette) | shade){

# plotting the shade
for(i in 1:n.cell){
  if(shade.scale){j <- i+1} else{j <- i}
  polygon(x = cell.array[,,i], col = idcolor[as.numeric(dit2[i])], border = NA)
}

}


# plotting the grids
for (i in 1:m){lines(ngrid[(1:n)+(i-1)*n,],lwd=grid.aes$lwd, col=addTrans(grid.aes$col, grid.aes$trans*255))}
for (i in 1:n){lines(ngrid[(1:m)*n-i+1,],lwd=grid.aes$lwd, col=addTrans(grid.aes$col, grid.aes$trans*255))}


# plotting the landmarks. reference landmarks under the target landmarks
if(plot.ref.lm){points(matr,pch = 19,cex=ref.lm.aes$cex, col = addTrans(ref.lm.aes$col, ref.lm.aes$trans*255))}
if(plot.lm){points(M2,pch = 19,cex=lm.aes$cex, col = addTrans(lm.aes$col, lm.aes$trans*255))}


if(plot.ref.links){
for (i in 1:nrow(links)) {
segments(matr[links[i, 1], 1], matr[links[i, 1], 2], matr[links[i, 2], 1], matr[links[i, 2], 2],
         lwd = ref.link.aes$lwd, col = addTrans(ref.link.aes$col, ref.link.aes$trans*255))}
  }
if(!is.null(links)){
for (i in 1:nrow(links)) {
segments(M2[links[i, 1], 1], M2[links[i, 1], 2], M2[links[i, 2], 1], M2[links[i, 2], 2],
         lwd = link.aes$lwd, col = addTrans(link.aes$col, link.aes$trans*255))}
  }

if(plot.arrows){
  for(i in 1:nrow(M2)){
    arrows(x0 = matr[i,1], y0 = matr[i,2], x1 = M2[i,1], y1 = M2[i,2], code = arrow.aes$code,
           col = addTrans(arrow.aes$col, arrow.aes$trans*255), lwd = arrow.aes$lwd, length = arrow.aes$length)}
  }

}


###########

########

  links.dor<-matrix(nrow=14, ncol=2)
  links.dor[,1] <- c(4,1,1,2,14,13,1,11,10,7,7,3,4,5) ; links.dor[,2] <- c(15,15,2,3,15,14,13,12,11,8,9,8,5,5)
  links.lat<-matrix(nrow=18, ncol=2)
  links.lat[,1] <- c(6,1,19,18,1,1,2,12,3,4,8,6,9,16,15,1,10,10) ; links.lat[,2] <- c(20,20,20,19,18,2,3,13,5,5,9,8,21,17,16,15,11,14)
  lat.drop <- c(-4,-5,-24,-23) ; dor.drop <- c(-4,-5,-6,-7,-9,-23,-24,-25,-14,-15)

# ruminant.TPS function
# automatically plots ruminants, making things easier for me (us)
# r can be either the scientific name (with a underscore in the middle)
  # or a landmark configuration
  #####
ruminant.TPS <- function(r, ref = NULL, both = FALSE, just.dorsal = FALSE, title = FALSE, shade.scale = NULL){
  par(bg=NULL)

  links.dor<-matrix(nrow=14, ncol=2)
  links.dor[,1] <- c(4,1,1,2,14,13,1,11,10,7,7,3,4,5) ; links.dor[,2] <- c(15,15,2,3,15,14,13,12,11,8,9,8,5,5)
  links.lat<-matrix(nrow=18, ncol=2)
  links.lat[,1] <- c(6,1,19,18,1,1,2,12,3,4,8,6,9,16,15,1,10,10) ; links.lat[,2] <- c(20,20,20,19,18,2,3,13,5,5,9,8,21,17,16,15,11,14)
  lat.drop <- c(-4,-5,-24,-23) ; dor.drop <- c(-4,-5,-6,-7,-9,-23,-24,-25,-14,-15)

   just.ventral <- just.dorsal # i'm super lazy
  if(is.character(r) && length(r) == 1){
    A <- R[lat.drop,1:2,r]
    A2 <- R[dor.drop,c(1,3),r]}else{
      A <- r[lat.drop,1:2]
      A2 <- r[dor.drop,c(1,3)]}
  if(is.null(ref)){
  B <- mshape(R)[lat.drop,1:2]; B2 <- mshape(R)[dor.drop,c(1,3)]}else
 {B <- ref[lat.drop,1:2]; B2 <- ref[dor.drop,c(1,3)]}

  A2[,2] <- -A2[,2]
  B2[,2] <- -B2[,2]

if(both){par(mfrow=c(2,1))}else{par(mfrow=c(1,1))}

  if(!just.ventral){
  tps(target.lm = A,
      reference.lm = B,
      links=links.lat,
    n.grid.col = 50,
    n.grid.row = NULL,
    grid.aes = list(col="gray", lwd=1, trans=1),
    mag=1,
    shade = TRUE,
    palette = col.pal,
    shade.trans = 0.5,
    plot.ref.links = TRUE
    )
    if(is.character(r) && length(r) == 1 && title){title(main = r)}

# plots lines to show face length
# if(!both & !just.ventral){
# lines(x = c(A[8,1],A[1,1]), y = c(-0.16,-0.16), col = family.colors[1], lwd = 10)
# lines(x = c(A[1,1],A[3,1]), y = c(-0.16,-0.16), col = family.colors[2], lwd = 10)
# middle <- mean(c(A[8,1],A[3,1]))
# lines(x = c(middle,middle), y = c(-0.17,-0.15), col = 'black', lwd = 7)
# #lines(x = c(A[1,1],A[1,1]), y = c(-0.17,-0.15), col = 'gray', lwd = 7)
# }}

if(both | just.ventral){
  tps(target.lm = A2,
      reference.lm = B2,
      links=links.dor,
    n.grid.col = 50,
    n.grid.row = NULL,
    grid.aes = list(col="gray", lwd=1, trans=1),
    mag=1,
    shade = TRUE,
    palette = col.pal,
    shade.trans = 0.5,
    plot.ref.links = TRUE
    )
}

# if(just.ventral){
#   if(is.character(r) && length(r) == 1 && title){title(main = r)}
# lines(x = c(A2[4,1],A2[1,1]), y = c(-0.16,-0.16), col = family.colors[1], lwd = 10)
# lines(x = c(A2[1,1],A2[3,1]), y = c(-0.16,-0.16), col = family.colors[2], lwd = 10)
# middle <- mean(c(A2[4,1],A2[3,1]))
# lines(x = c(middle,middle), y = c(-0.17,-0.15), col = 'black', lwd = 7)
#lines(x = c(A[1,1],A[1,1]), y = c(-0.17,-0.15), col = 'gray', lwd = 7)
}}


#####


# other functions
#############
# takes a vector and scales and returns a color palette accordingly
make.color.vec <- function(x, n, pallette, rev = F){
  x1 <- cut(as.numeric(x),breaks=n)
if(!rev){idcolor<-pallette(n)}else{id.color<-rev(pallette(n))}
x.cols <- 1:length(x1)
for(i in 1:length(x.cols)){x.cols[[i]] <- idcolor[as.numeric(x1[i])]}
names(x.cols) <- names(x)
return(x.cols)
}

centroid1 <- function(lm){apply(lm, 2, mean)}

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

# automatically calculates face length from the landmark configurations
face.length <- function(x){
  return(euc.dist(x[8,],x[18,])/ euc.dist(x[8,],x[3,]))
} # same index as Bibi & Tyler 2022
nasal.retraction <- function(x){
  return(euc.dist(x[8,],x[22,])/ euc.dist(x[8,],x[3,]))
}


# https://www.r-bloggers.com/2011/11/outersect-the-opposite-of-rs-intersect-function/
outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}


# calculating relative size of the face module
face.size <- function(x){
return(cSize(mshape(x)[face,]))
} # I don't use this


# plots a color palette
plot.cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)), col=a, axes=FALSE , xlab="", ylab="")

  # from the source code of the 'rescale' function in the scales package
# for some reason just using 'rescale' gives me errors sometimes so I've just been using this
rescale.numeric <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE), ...) {
  require(scales)
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  }
  (x - from[1]) / diff(from) * diff(to) + to[1]
  }


  # addTrans adds transparency to a vector of colors
#  Retrieved from http://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


# scale a vector to unit length
# https://stackoverflow.com/questions/27455948/in-r-whats-the-simplest-way-to-scale-a-vector-to-a-unit-vector
scalar1 <- function(x) {x / sqrt(sum(x^2))}
# could also use the 'Normalize' function from the 'evolqg' package

# calculates angle between multivariate vectors, assuming direction of the axes is arbitrary (like in PCs)
# returns angle in degrees
vector.angle <- function(x,y){
      x <- as.vector(x)/sqrt(sum(x^2))
    y <- as.vector(y)/sqrt(sum(y^2))
    rho <- acos(x%*%y)
    rho <- (rho * 180) / (pi)
    if(rho>90){rho <- (180-rho)}
    return(rho)
    }
# could also just use morpho functions
# ie angle.calc

# from: https://stackoverflow.com/questions/32370485/convert-radians-to-degree-degree-to-radians
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}


# calculates morphological integration va eigenvalue dispersion
# default: deviation of the relative eigenvalue dispersion given the sample size
eigen.var <- function(lm, sd = FALSE, rel = TRUE, sample = TRUE, cor = FALSE){
  lm2 <- two.d.array(lm)
  if(sample){n <- dim(lm)[3]}else{n <- NULL}
  v <- var(lm2)
  if(cor){v<-cov2cor(v)}
  CalcEigenVar(matrix = v, sd = sd, rel = rel, sample = n)
}




# returns n equally spaced colours from along the color wheel. from stack overflow
colourWheel <- function(n = n, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# vectorizes a landark configuration
flatten <- function(x){
return( two.d.array(A=bindArr(x, x, along = 3))[1,] )
}

unflatten <- function(x,  k = 3){
 return(  arrayspecs(A = rbind(x, x), k = k, p = length(x)/k)[,,1] )
}


# add a tree to a morphospace
# modified version of an internal RRPP function "addtree"
# automatically calculates ancestors using fastAnc()
addtree2 <- function(x, tree, col = "black", alpha = 1, lwd = 1){
  # if(anc == "fast"){
  #   anc.x <- cbind(fastAnc(tree = tree, x = x[,1]), fastAnc(tree = tree, x = x[,2]))}
  # if(anc == "ML"){
  #   anc.x <- # ... anc.ML() ...
  # }
  # if(anc == "Bayes"){
  #   anc.x <- # ... anc.Bayes() ...
  # }
  #
  anc.x <- cbind(fastAnc(tree = tree, x = x[,1]), fastAnc(tree = tree, x = x[,2]))
trans <- alpha*255
pts <- x[,1:2]
anc <- anc.x[,1:2]
ind <- match(tree$tip.label, rownames(pts))
pts <- pts[ind, ]
z <- rbind(pts,anc)
# coordinates of each tip and node in the tree
edges <- as.matrix(tree$edge)
# edges is two column matrix where the first row is the index of the ancestor node
  for (i in 1:NROW(edges)) {
          pts <- z[edges[i, ], ]
          points(pts, type = "l", col = addTrans(color = col,trans = trans), lwd = lwd, lty = 1)}
}


# a depriciated (?) Morpho function
restoreShapes <- function(scores,PC,mshape)
  {
    dims <- dim(mshape)
    PC <- as.matrix(PC)

    if (!is.matrix(scores) && ncol(PC) == 1)
        if (length(scores) > 1)
            scores <- as.matrix(scores)
    if (!is.matrix(scores)){
        if (length(scores) != ncol(PC))
            stop("scores must be of the same length as ncol(PC)")
        predPC <- PC%*%scores
        modell <- mshape+matrix(predPC,dims[1],dims[2])
        return(modell)
    } else {
          n <- nrow(scores)
          outarr <- array(0,dim=c(dims,n))
          for (i in 1:n) {
              outarr[,,i] <- restoreShapes(scores[i,],PC,mshape)
          }
          if (!is.null(rownames(scores)))
              dimnames(outarr)[[3]] <- rownames(scores)
          return(outarr)
    }
}



# calculates the euclidean norm of a vector
norm_vec <- function(x) sqrt(sum(x^2))



# calculates the projected variance of a vector
# ie, how much variance of a P matrix would a given vector explain
# scales the vector to unit length and accounts for the total variance in the P matrix
projected.variance <- function(P, z, scale = TRUE){
  if(scale){z <- (z / sqrt(sum(z^2)))} # scale to unit length
  trace <- sum(diag(P))
  Vproj <- t(z)%*%P%*%z/trace
  return(Vproj)
}



# translate.lm translates a matrix or all the lm configurations in array to a common position
# 'coords' is a 3D coordinate array
# 'centroid' is the xyz coordinate you want to translate to
translate.lm <- function(coords,centroid = c(0,0,0)){
  # 'trans1' is from 'Morphometrics with R' by Julien Claude, translates lm configuration to centroid
  trans1 <- function(M) {M-matrix(centroid1(M),nrow(M),ncol(M),byrow=T)}
  k <- dim(coords)[[2]]

  if(length(dim(coords) == 2)){
  A<-matrix(NA,nrow=nrow(coords),ncol=ncol(coords))
  coords <- trans1(coords)
	A[,1]<-coords[,1]+centroid[1]
	A[,2]<-coords[,2]+centroid[2]
	if(k==3){A[,3]<-coords[,3]+centroid[3]}
	if(!is.null(dimnames(coords))){dimnames(A)<-dimnames(coords)}
	return(A)}

  else{
  A<-array(NA,dim=dim(coords))
	for(i in 1:dim(coords)[[3]]){
	  coords[,,i] <- trans1(coords[,,i])}
  A[,1,]<-coords[,1,]+centroid[1]
	A[,2,]<-coords[,2,]+centroid[2]
	if(k==3){A[,3,]<-coords[,3,]+centroid[3]}
	dimnames(A)<-dimnames(coords)

	return(A)}
}




# an old function I would use to calculate grids along a morphospace, there's a quicker way to do this though (expand.grid() ... )
  pca.grid <- function(pca, PC1.segments = 7, PC2.segments = 5, pc1.min.buffer = 0, pc1.max.buffer = 0, pc2.min.buffer = 0, pc2.max.buffer = 0){
  mshape <- mshape(pca$A)
  p = dim(pca$A)[1]
  k = dim(pca$A)[2]
  n = dim(pca$A)[3]
  eigenvectors <- arrayspecs(t(pca$rotation), p = p, k = k)
  pc1.min <- min(pca$x[,1])*(1+pc1.min.buffer)
  pc1.max <- max(pca$x[,1])*(1+pc1.max.buffer)
  pc2.min <- min(pca$x[,2])*(1+pc2.min.buffer)
  pc2.max <- max(pca$x[,2])*(1+pc2.max.buffer)
  pc1.scores <- seq(pc1.min,pc1.max,(pc1.max-pc1.min)/PC1.segments)
  pc2.scores <- seq(pc2.min,pc2.max,(pc2.max-pc2.min)/PC2.segments)
  total.shapes <- PC1.segments*PC2.segments
  shape.grid <- array(data = NA, dim = c(p, k, total.shapes))
  score.grid <- matrix(nrow=total.shapes,ncol=2)
  for(i in 1:PC2.segments){
    for(j in 1:PC1.segments){
      shape.grid[,,(j+PC1.segments*(i-1))] <- mshape + eigenvectors[,,1]*pc1.scores[j] + eigenvectors[,,2]*pc2.scores[i]
      score.grid[(j+PC1.segments*(i-1)),] <- c(pc1.scores[j],pc2.scores[i])
      }
  }
  return(list(shape.grid = shape.grid,score.grid = score.grid))
}


#############



  # idk why this isn't automaitcally loading
  # but this is from the RRPP package:
  scaleCov <- function(Cov, scale. = 1, exponent = 1,
                     scale.diagonal = FALSE,
                     scale.only.diagonal = FALSE) {

  dims <- dim(Cov)
  if(length(dims) != 2)
    stop("Cov must be a matrix", call. = FALSE)
  if(dims[1] != dims[2])
    stop("Cov must be a square matrix", call. = FALSE)
  if(length(scale.) != 1)
    stop("The scale. argument must be a single numeric value.",
         call. = FALSE)
  if(length(exponent) != 1)
    stop("The exponent argument must be a single numeric value.",
         call. = FALSE)
  if(!is.numeric(scale.) || !is.numeric(exponent))
    stop("Scaling parameters must be numeric.",
         call. = FALSE)

  D <- diag(diag(Cov))
  C <- Cov - D
  if(scale.only.diagonal) {
    D <- scale. * D^(exponent)
  } else {
    C <- scale. * C^(exponent)
    if(scale.diagonal) D <- scale. * D^(exponent)
  }

  C + D
}

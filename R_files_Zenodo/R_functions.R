calculate_curvature <- function(x, y=NULL, z=NULL){

	z_zero <- TRUE
	if(is.null(y) && ncol(x) > 1){
		if(ncol(x) == 2){
			z <- rep(0, length(x))
		}else{
			z <- x[, 3]
			z_zero <- FALSE
		}
		y <- x[, 2]
		x <- x[, 1]
	}

	xyz <- cbind(x,y,z)

	fit_x <- loess(y ~ t, data=data.frame(t=1:length(x), y=x))
	fit_y <- loess(y ~ t, data=data.frame(t=1:length(x), y=y))
	if(!z_zero) fit_z <- loess(y ~ t, data=data.frame(t=1:length(x), y=z))
	
	# Get smoothed points at higher sampling
	idx_s <- seq(1, length(x), by=1)
	x_s <- predict(fit_x, idx_s)
	y_s <- predict(fit_y, idx_s)
	if(z_zero){
		z_s <- rep(0, length(idx_s))
	}else{
		z_s <- predict(fit_z, idx_s)
	}

	diff_x <- diff(x_s)
	dx <- colMeans(rbind(c(NA, diff_x), c(diff_x, NA)), na.rm=TRUE)

	diff_dx <- diff(dx)
	ddx <- colMeans(rbind(c(NA, diff_dx), c(diff_dx, NA)), na.rm=TRUE)
	
	diff_y <- diff(y_s)
	dy <- colMeans(rbind(c(NA, diff_y), c(diff_y, NA)), na.rm=TRUE)

	diff_dy <- diff(dy)
	ddy <- colMeans(rbind(c(NA, diff_dy), c(diff_dy, NA)), na.rm=TRUE)
	
	diff_z <- diff(z_s)
	dz <- colMeans(rbind(c(NA, diff_z), c(diff_z, NA)), na.rm=TRUE)

	diff_dz <- diff(dz)
	ddz <- colMeans(rbind(c(NA, diff_dz), c(diff_dz, NA)), na.rm=TRUE)

	k <- sqrt((ddz*dy - ddy*dz)^2 + (ddx*dz - ddz*dx)^2 + (ddy*dx - ddx*dy)^2) / ((dx^2 + dy^2 + dz^2)^(3/2))

	if(FALSE){

		layout(cbind(1:3))
		par(mar=c(2,2,2,2))
		plot(xyz[, 1])
		points(predict(fit_x, idx_s), type='l', col='blue')
		plot(xyz[, 2])
		points(predict(fit_y, idx_s), type='l', col='blue')
		plot(xyz[, 3])
		points(predict(fit_z, idx_s), type='l', col='blue')
	}

	return(k)
}

draw_checkerboard <- function(nx, ny, square.size, max.circle.size, num.circles, circle.lwd, min.circle.size=30, filename, margin.x = c(round(square.size/2), round(square.size/2)), margin.y = c(round(square.size/2), round(square.size/2)), ...){

	lines.lwd <- 4

	# REQUIRES GRID PACKAGE

	# GET X AND Y COORDINATES FOR CIRCLES
	xy_circle <- c(max.circle.size/2, max.circle.size/2)

	# GET X AND Y COORDINATES FOR SQUARES
	x <- c(0, 0, square.size, square.size) + max.circle.size + margin.x[1]
	y <- c(0, square.size, square.size, 0)
	
	# GET FILE EXTENSION FROM FILENAME
	image_type <- tolower(substr(filename, nchar(filename)-attributes(regexpr(pattern='[A-Za-z]+$', text=filename))$match.length+1, nchar(filename)))
	
	# CHANGE JPG TO JPEG TO MATCH FUNCTION CALL
	if(image_type == 'jpg') image_type <- 'jpeg'
	
	# GET WIDTH AND HEIGHT
	width <- sum(margin.x) + (nx+1)*square.size + max.circle.size + margin.x[1]
	height <- sum(margin.y) + (ny+1)*square.size
	
	# CALL CORRESPONDING IMAGE FUNCTION TO START IMAGE WRITING
	do.call(image_type, list(filename=filename, width=width, height=height, ...))

	# WRITE CHECKERBOARD SQUARES TO IMAGE
	grid.newpage()
	for(i in seq(0, nx, by=2)){
		for(j in seq(0, ny, by=1)){

			# IF nx IS EVEN, i IS LAST NUMBER AND j IS ODD CONTINUE - THIS IS AN EXTRA ROW THAT SHOULD BE SKIPPED
			if(i == nx && j %% 2 == 1) next
			
			# DRAW BLACK SQUARE
			grid.polygon(x=margin.x[1]+x+square.size*i + square.size*(j %% 2), y=margin.y[1]+y+square.size*j, gp=gpar(fill=1, lwd=0), default.units="native")
		}
	}
	
	# DRAW CONCENTRIC CIRCLES
	for(r in seq(min.circle.size/2, max.circle.size/2, length=num.circles)){
		grid.circle(x=margin.x[1]+xy_circle, y=margin.y[1]+xy_circle, r=r, gp=gpar(fill=0, lwd=circle.lwd), default.units="native")
	}
	#for(i in seq(rmin.circles, 1-num.circles*rstep.circles, length=num.circles)){
	#	grid.circle(x=margin.x[1]+xy_circle, y=margin.y[1]+xy_circle, r=(max.circle.size/2)*i, gp=gpar(fill=0, lwd=circle.lwd), default.units="native")
	#}
	
	# ADD CROSSED LINES IN CORNERS
	# TOP LEFT
	#grid.segments(margin.x[1]/2, margin.y[1]/4, margin.x[1]/2, margin.y[1]/4+margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")
	#grid.segments(margin.x[1]/4, margin.y[1]/2, margin.x[1]/4+margin.x[1]/2, margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")

	# BOTTOM RIGHT
	#grid.segments(width-margin.x[1]/2, height-margin.y[1]/4, width-margin.x[1]/2, height-margin.y[1]/4-margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")
	#grid.segments(width-margin.x[1]/4, height-margin.y[1]/2, width-margin.x[1]/4-margin.x[1]/2, height-margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")

	# TOP RIGHT
	#grid.segments(width-margin.x[1]/2, margin.y[1]/4, width-margin.x[1]/2, margin.y[1]/4+margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")
	#grid.segments(width-margin.x[1]/4, margin.y[1]/2, width-margin.x[1]/4-margin.x[1]/2, margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")

	# BOTTOM LEFT
	#grid.segments(margin.x[1]/2, height-margin.y[1]/4, margin.x[1]/2, height-margin.y[1]/4-margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")
	#grid.segments(margin.x[1]/4, height-margin.y[1]/2, margin.x[1]/4+margin.x[1]/2, height-margin.y[1]/2, gp=gpar(lwd=lines.lwd), default.units="native")

	# CLOSE IMAGE CONNECTION
	dev.off();
}

draw_fin_kinematics <- function(xyz_arr, file, window.title, path.connect, animate.duration){

	svg.new(file=file, window.title=window.title, animate.duration=animate.duration)

	# Get landmarks
	lm_arr <- xyz_arr[!dimnames(xyz_arr)[[1]] %in% c('body_centroid', 't_vec', 'a_vec', 'd_vec'), , ]
	lm_arr <- lm_arr[!grepl('_pre_align$', dimnames(lm_arr)[[1]]), , ]

	# Draw landmarks
	svg.points(lm_arr, cex=1)

	for(i in 1:length(path.connect)){
	
		if(length(path.connect[[i]]) == 1){
			paths_connect <- dimnames(lm_arr)[[1]][grepl(path.connect[[i]], dimnames(lm_arr)[[1]])]
		}else{
			paths_connect <- path.connect[[i]][path.connect[[i]] %in% dimnames(lm_arr)[[1]]]
		}
		
		if(length(paths_connect) < 2) next
		
		path_idx <- rep(NA, length(paths_connect))
		for(j in 1:length(path_idx)) path_idx[j] <- which(paths_connect[j] == dimnames(lm_arr)[[1]])
		
		svg.pathsC(path_idx)
	}
	
	# Draw coordinate system
	if('body_centroid' %in% dimnames(xyz_arr)[[1]] && sum(!is.na(xyz_arr['body_centroid', 1, ])) > 0){
		svg.points(t(xyz_arr['body_centroid', , !is.na(xyz_arr['body_centroid', 1, ])]))
		svg.arrows(xyz_arr[c('body_centroid', 'a_vec'), , !is.na(xyz_arr['a_vec', 1, ])], col='red', len=1.5)
		svg.arrows(xyz_arr[c('body_centroid', 'd_vec'), , !is.na(xyz_arr['d_vec', 1, ])], col='green', len=1.5)
		svg.arrows(xyz_arr[c('body_centroid', 't_vec'), , !is.na(xyz_arr['t_vec', 1, ])], col='blue', len=1.5)
	}
	
	# Draw midline plane
	#ranges <- apply(xyz_arr, 2, 'range', na.rm=TRUE)

	if(sum(grepl('_pre_align$', dimnames(xyz_arr)[[1]])) > 0){
		
		lms_pre <- xyz_arr[grepl('_pre_align$', dimnames(xyz_arr)[[1]]), , ]

		svg.points(lms_pre, cex=1, col='red')
	}
	
	# Draw frame
	svg.frame(xyz_arr, z.index=-1)
	
	svg.close()
}

path.connect.curves <- list(
	c('FR[0-9]+_base_R'),
	c('FR[0-9]+_tip_R'),
	c('FR[0-9]+_base_L'),
	c('FR[0-9]+_tip_L')
)

for(i in 1:14) path.connect.curves[[length(path.connect.curves)+1]] <- paste0('FR', formatC(i, width=2, format="d", flag="0"), '_L_[0-9]+')
for(i in 1:14) path.connect.curves[[length(path.connect.curves)+1]] <- paste0('FR', formatC(i, width=2, format="d", flag="0"), '_R_[0-9]+')

for(i in 1:14){
	path.connect.curves[[length(path.connect.curves)+1]] <- c(
			paste0('FR', formatC(i, width=2, format="d", flag="0"), '_base_L'), 
			paste0('FR', formatC(i, width=2, format="d", flag="0"), '_tip_L')
		)
	path.connect.curves[[length(path.connect.curves)+1]] <- c(
			paste0('FR', formatC(i, width=2, format="d", flag="0"), '_base_R'), 
			paste0('FR', formatC(i, width=2, format="d", flag="0"), '_tip_R')
		)
}

path.connect.curves[[length(path.connect.curves)+1]] <- c('dent_midline_inf', 'mid_operc_ser')
path.connect.curves[[length(path.connect.curves)+1]] <- c('mid_operc_ser', 'anal_fin_ant')
path.connect.curves[[length(path.connect.curves)+1]] <- c('anal_fin_ant', 'anal_fin_pos')
path.connect.curves[[length(path.connect.curves)+1]] <- c('pec_fin_ant_L', 'pec_fin_ant_R')
#path.connect.curves[[length(path.connect.curves)+1]] <- c('body_centroid', 't_vec')
#path.connect.curves[[length(path.connect.curves)+1]] <- c('body_centroid', 'a_vec')
#path.connect.curves[[length(path.connect.curves)+1]] <- c('body_centroid', 'd_vec')

path.connect <- list(
	c('FR01_[0-9]+'),
	c('FR11_[0-9]+'),
	c('FR[0-9]+_base'),
	c('FR[0-9]+_tip'),
	c('FR01_base', 'FR01_tip'),
	c('FR02_base', 'FR02_tip'),
	c('FR03_base', 'FR03_tip'),
	c('FR04_base', 'FR04_tip'),
	c('FR05_base', 'FR05_tip'),
	c('FR06_base', 'FR06_tip'),
	c('FR07_base', 'FR07_tip'),
	c('FR08_base', 'FR08_tip'),
	c('FR09_base', 'FR09_tip'),
	c('FR10_base', 'FR10_tip'),
	c('FR11_base', 'FR11_tip'),
	c('FR12_base', 'FR12_tip'),
	c('FR13_base', 'FR13_tip'),
	c('FR14_base', 'FR14_tip'),
	c('FR15_base', 'FR15_tip')
)

path.connect.edges <- list(
	c('FR01_[0-9]+', 'FR02_tip'),
	c('FR11_[0-9]+', 'FR10_tip'),
	c('FR01_00001', 'FR[0-9]+_base', 'FR11_00001'),
	c('FR[0-9]+_tip'),
	c('FR01_base', 'FR01_tip'),
	c('FR02_base', 'FR02_tip'),
	c('FR03_base', 'FR03_tip'),
	c('FR04_base', 'FR04_tip'),
	c('FR05_base', 'FR05_tip'),
	c('FR06_base', 'FR06_tip'),
	c('FR07_base', 'FR07_tip'),
	c('FR08_base', 'FR08_tip'),
	c('FR09_base', 'FR09_tip'),
	c('FR10_base', 'FR10_tip'),
	c('FR11_base', 'FR11_tip'),
	c('FR12_base', 'FR12_tip'),
	c('FR13_base', 'FR13_tip'),
	c('FR14_base', 'FR14_tip'),
	c('FR15_base', 'FR15_tip')
)

set_colors <- function(scheme = c('bluegreen', 'redblack')){

	if(scheme[1] == 'bluegreen'){
		col_control <- c(0,169,0) / 255
		col_control_light <- c(0,200,0) / 255
		col_transect <- c(3,101,255) / 255
		col_transect_light <- c(50,150,255) / 255
	}else{
		col_control <- rep(70,3) / 255
		col_control_light <- rep(120,3) / 255
		col_transect <- c(190,0,0) / 255
		col_transect_light <- c(230,0,0) / 255
	}

	col_control_rgb <- rgb(col_control[1], col_control[2], col_control[3])
	col_control_light_rgb <- rgb(col_control_light[1], col_control_light[2], col_control_light[3])
	col_transect_rgb <- rgb(col_transect[1], col_transect[2], col_transect[3])
	col_transect_light_rgb <- rgb(col_transect_light[1], col_transect_light[2], col_transect_light[3])

	list(
		'col_c_rgb'=col_control_rgb,
		'col_c_light_rgb'=col_control_light_rgb,
		'col_c_rgba'=rgb(col_control_light[1], col_control_light[2], col_control_light[3], 0.3),
		'col_t_rgb'=col_transect_rgb,
		'col_t_light_rgb'=col_transect_light_rgb,
		'col_t_rgba'=rgb(col_transect_light[1], col_transect_light[2], col_transect_light[3], 0.3)
	)
}

std <- function(x, ...) sd(x, ...)/sqrt(sum(!is.na(x)))
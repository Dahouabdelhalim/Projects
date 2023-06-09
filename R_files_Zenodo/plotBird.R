# Function to plot curved lines (graph object) on a 2D image made of polygons
# x <- sumchanges.2
# x = adjacency matrix (counts of transitions to/from different plumage regions)
# direc = whether to use anterior, posterio directions
# col = color
# TODO: add other options for user to input images, polygons, patch names, etc.

plotBird <- function(x, col = 'black', direc = c('both', 'anterior', 'posterior'),
		lwd = 1, maxweight = NULL, pch = 16, pt.col = 'gray', pt.bg = 'gray', cex = 1) {

	direc <- match.arg(direc)

	require(sna)
	require(ggplot2)
	require(Hmisc)
	
	adjacencyMatrix <- x

	codes <- read_excel("~/Dropbox/king_plumage/data/kingfisher_plumage_mechanisms.xlsx", sheet=4)

	# get coordinates of patches
	raw <- readLines("data/alcedo patches.svg")
	nms <- as.character(str_match_all(raw, '<polygon.*?id\\\\=\\\\"(.*?)\\\\"')[[1]][, 2])
	pts <- str_match_all(raw, '<polygon.*?points\\\\=\\\\"(.*?)\\\\"')[[1]][, 2]
	pts <- sapply(seq_along(pts), function(x) {as.numeric(strsplit(pts[[x]], " ")[[1]])})
	coords <- lapply(seq_along(pts), function(i) {
		ss <- pts[[i]]
		x <- ss[seq(1, length(ss),by=2)]
		y <- ss[seq(2, length(ss),by=2)]
		x <- c(x, x[1])
		y <- c(y, y[1])
		cbind(x, y)
	})
	names(coords) <- nms

	xy <- t(sapply(coords, apply, 2, mean))
	xy <- xy[rownames(xy)!="beak", ]

	xy <- xy[match(rownames(xy), codes$illustrator_code), ]

	# rotate points around front of head points
	pt1 <- xy['fronthead', ]
	pt2 <- xy['tail', ]
	dx <- pt1[1] - pt2[1]
	dy <- pt1[2] - pt2[2]
	theta <- atan(dy / dx)
	rotmat <- matrix(c(cos(theta), - sin(theta), sin(theta), cos(theta)), nrow = 2)
	center <- pt2
	xy2 <- t(rotmat %*% apply(xy, 1, '-', center) + center) # centered and rotated coordinates
	colnames(xy2) <- c('x', 'y')

	layoutCoordinates <- xy

	adjacencyList <- reshape2::melt(adjacencyMatrix)  # Convert to list of ties only
	adjacencyList <- adjacencyList[adjacencyList$value > 0, ] # remove zeros along diagonal

# ggplot(data.frame(layoutCoordinates), aes(x, y)) + geom_point() + geom_text(aes(label=rownames(layoutCoordinates))) + coord_fixed()

	# Function to generate paths between each connected node
	edgeMaker <- function(whichRow, len = 100, curved = TRUE){
		fromC <- layoutCoordinates[adjacencyList[whichRow, 1], ]  # Origin
		toC <- layoutCoordinates[adjacencyList[whichRow, 2], ]  # Terminus
		fromC_rot <- xy2[adjacencyList[whichRow, 1], ]  # Origin
		toC_rot <- xy2[adjacencyList[whichRow, 2], ]  # Terminus
		direc <- ifelse(toC_rot['x'] - fromC_rot['x'] > 0, 'posterior', 'anterior')
		wt <- adjacencyList[whichRow, 3]  # Weight
		# Add curve:
		graphCenter <- colMeans(layoutCoordinates)  # Center of the overall graph
		bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
		distance1 <- sum((graphCenter - bezierMid)^2)
		if (distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)) {
			bezierMid <- c(toC[1], fromC[2])
	    }  # To select the best Bezier midpoint
		bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
		if (curved == FALSE) {
			bezierMid <- (fromC + toC) / 2
		}  # Remove the curve
		edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
	                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
				evaluation = len))  # Bezier path coordinates
		edge$weight <- wt
		edge$Sequence <- 1:len  # For size and colour weighting in plot
		edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")
		edge$direc <- direc
		return(edge)
	}

	# Generate a (curved) edge path for each pair of connected nodes
	allEdges <- lapply(1:nrow(adjacencyList), edgeMaker, len = 50, curved = TRUE)
	allEdges <- do.call(rbind, allEdges)  # a fine-grained path ^, with bend ^

	# create plot structure
	if (direc != 'both') {
		zp1 <- ggplot(allEdges[allEdges$direc==direc, ])
	} else {
		zp1 <- ggplot(allEdges)
	}

	# add bird cartoon
	df <- plyr::ldply(coords, .id="patch")
	df$id <- match(df$patch, codes$illustrator_code)
	df$col[df$patch=="beak"] <- NA
	zp1 <- zp1 + geom_polygon(data = df, aes(x=x, y=y, group=patch, fill=col), col="gray", lwd = 0.25)

	# Add nodes
	N <- length(levels(df$patch))-1 # number of unique patches - removing the beak

	if (length(pt.col)==1) {
		pt.col <- rep(pt.col, N)
	}

	# add paths

	if (is.null(maxweight)) {
		maxweight <- max(allEdges$weight)
	}

	zp1 <- zp1 + geom_path(aes(x = x, y = y, group = Group,  # Edges with gradient
					# lwd = (weight/max(weight))^lwd,
					lwd = (weight/maxweight)^lwd,
					alpha = (weight/maxweight)^2),
					colour = col) # colour = Sequence, size = -Sequence))  # and taper
	
	# zp1 <- zp1 + scale_colour_gradient(low = gray(0), high = gray(9/10), guide = "none")
	# zp1 <- zp1 + scale_colour_gradient(low = cols[1], high = cols[2], guide = "none")
	zp1 <- zp1 + scale_size(range = c(1/10, 1), guide = "none")  # Customize taper

	zp1 <- zp1 +
		# geom_polygon(data = df, aes(x=x, y=y, group=patch, fill=col), col="gray", lwd = 0.25) +
		scale_fill_identity() +
		coord_fixed() +
		scale_x_reverse() +
		scale_y_reverse() +
		theme_void() +
		theme(legend.position = "none")

	zp1 <- zp1 + geom_point(data = data.frame(layoutCoordinates),
	                        aes(x = x, y = y), size = cex, pch = pch, fill = pt.bg,
							colour = pt.col)  # Customize gradient v

	zp1

}

library(rpanel)
library(tkrplot)
library(tcltk)

# packages needed to be installed and loaded
# rpanel, tkrplot, tcltk

debugflag <<- TRUE

# for stochastic, make it an exponential distribution e^beta*fi / normalization
# beta is a parameter along with death rate
# look up cran documentation
# contour plotting for tirangles or hex? low priority
# fitness should get from cell i,j on other pop
# read snowdrift journal of theoretical biology to doublecheck parameters and try to replicate game dynamics
# read speciation paper
# tag and preference and tolerance

# every period, highest fitness neighbor wins
# simultaneous or sequential, sequential
# every tick, random neighborhood chosen to update, 
# sweep as well

# prioritize death rates
# frequency dependence would be each cell gets taken over by max fitness
# public goods check paper because death is differently structured
# make CA work within a package style system


drawGrid <- function(panel)
{
	dbg('drawGrid')

    par(mar=rep(0,4))
    plot.new()
    rect(0, 0, 1, 1)
	#print("draw")
    #print(ca$curgrid)
	
	for (n in 1:length(ca$populations))
	{
		cg = ca$populations[[n]]$curgrid
		og = ca$populations[[n]]$oldgrid
		
		for (i in 1:ncol(cg))
		{
			for (j in 1:nrow(cg))
			{
				if (cg[i, j] == 0)
					color <- "#FFFFFF"
				else
					color <- ca$populations[[n]]$popColors[cg[i, j]]
				
				# {
					# if (og[i, j] == 1)
					# {
						# color <- "#FF8080"
					# }
					# else if (og[i, j] == 2)
					# {
						# color <- "#FFD280"
					# }	
					# else if (og[i, j] == 3)
					# {
						# color <- "#8080FF"
					# }
				# }
				# if (cg[i, j] == 1)
				# {
					# color <- "#FF0000"
				# }
				# else if (cg[i, j] == 2)
				# {
					# color <- "#FFA500"
				# }
				# else if (cg[i, j] == 3)
				# {
					# color <- "#0000FF"
				# }
				#n = paste("cell_", as.character(i), "_", as.character(j), sep="")
                if (ca$bwOnly)
                    color = ca$bwColors[cg[i, j] + 1]
				if (n == 1)
                {
					rect((i - 1) / nrow(cg), (j - 1) / ncol(cg), i / nrow(cg), j / ncol(cg), col=color)
					if (ca$viewFitness)
						text((i - 0.5) / nrow(cg), (j - 0.5) / ncol(cg), ca$populations[[n]]$fitnessMatrix[i, j])
                }
				else if (n == 2)
				{
					if (ca$panel$drawStyle == ca$drawStyles[1])
						rect((i - 0.75) / nrow(cg), (j - 0.75) / ncol(cg), (i - 0.25) / nrow(cg), (j - 0.25) / ncol(cg), col=color)
					else if (ca$panel$drawStyle == ca$drawStyles[2])
						rect((i - 0.5) / nrow(cg), (j - 0.5) / ncol(cg), (i - 0) / nrow(cg), (j - 0) / ncol(cg), col=color)
				}

			}
		}
	}
    rect(0, 0, 1, 1)
	text(0.25, -0.02, paste(ca$gameType, ", generation ", ca$tick))
	dbg('end drawGrid')
    panel
}

caStop <- function(panel)
{
	ca$running <<- FALSE
	rp.tkrreplot(ca$panel, tkrp)
	panel
}

caStart <- function(panel)
{
	ca$running <<- TRUE
	panel
}

makeFile <- function()
{
	settingsline = list()
	settingsline[[1]] = paste("game=", ca$gameType, ",", sep="")
	settingsline[[2]] = paste("size=", ca$size$x, "x", ca$size$y, ",", sep="")
	settingsline[[3]] = paste("surface=", ca$world, ",", sep="")
	settingsline[[4]] = "TODO: full game parameters go here!"

	firstline = list()
	firstline[[1]] = "generation"
	for (n in 1:length(ca$populations))
	{
		for (i in 1:ca$populations[[n]]$players)
		{
			firstline[[length(firstline)+1]] = paste("shares_pop", n, "_strat", i, collapse="", sep="")
		}
	}
	for (i in 1:length(settingsline))
	{
		writeLines(settingsline[[i]])
	}
	#s = paste(settingsline, collapse=",", sep="")
	#writeLines(s)
	
	s = paste(firstline, collapse=",", sep="")
	#print(s)
	writeLines(s)
	
	#s = lapply(firstline, paste, collapse="", sep=",")
	#write(s, file="", append=FALSE, ncolumns=1)
	#write(firstline, file="", ncolumns=length(firstline), append=FALSE, sep=",")
}

writeLine <- function()
{
	caline = list()
	caline[[1]] = ca$tick
	for (n in 1:length(ca$populations))
	{
		for (i in 1:ca$populations[[n]]$players)
		{
			caline[[length(caline)+1]] = ca$cellCounts[[n]][i]
		}
	}
	s = paste(caline, collapse=",", sep="")
	#print(s)
	writeLines(s)
	
	#s = lapply(caline, paste, collapse="", sep=",")
	#write(s, file="", append=FALSE, ncolumns=1)
	#write(caline, file="", ncolumns=length(caline), append=TRUE, sep=",")
}

countCells <- function()
{
    counts = list()
    for (n in 1:length(ca$populations))
	{
        counts[[n]] = numeric(length=ca$populations[[n]]$players)
		
		for (i in 1:ncol(ca$populations[[n]]$curgrid))
		{
			for (j in 1:nrow(ca$populations[[n]]$curgrid))
			{
				# first update cell counts for this generation
				counts[[n]][ca$populations[[n]]$curgrid[i, j]] = counts[[n]][ca$populations[[n]]$curgrid[i, j]] + 1
			}
		}
	}
    return(counts)
}

updateFitnessMatrix <- function()
{
    for (n in 1:length(ca$populations))
	{
		ca$populations[[n]]$fitnessMatrix <<- matrix(ncol=ncol(ca$populations[[n]]$curgrid), nrow=nrow(ca$populations[[n]]$curgrid), byrow=TRUE)

		for (i in 1:ncol(ca$populations[[n]]$curgrid))
		{
			for (j in 1:nrow(ca$populations[[n]]$curgrid))
			{
				# calculate fitness for each individual cell
				if (ca$populations[[n]]$curgrid[i, j] == 0)
					ca$populations[[n]]$fitnessMatrix[i, j] <<- 0
				else
				{
					neighborCoords = getNeighborCoords(i, j, ca$populations[[n]])
					neighborCounts = countNeighbors(neighborCoords, ca$populations[[n]])
                    
					ca$populations[[n]]$fitnessMatrix[i, j] <<- ca$populations[[n]]$fitnessFunction(ca$populations[[n]], ca$populations[[n]]$curgrid[i, j], neighborCounts) 
				}
			}
		}
	}
}

caUpdate <- function(panel)
{
	dbg('caUpdate')
	ca$tick <<- ca$tick + 1
	ca$cellCounts <<- countCells()
	
	globalfitness = list()
	for (n in 1:length(ca$populations))
	{	
		globalfitness[[n]] = numeric(length=ca$populations[[n]]$players)
		for (i in 1:ca$populations[[n]]$players)
		{
			globalfitness[[n]][i] = ca$populations[[n]]$fitnessFunction(ca$populations[[n]], i, ca$cellCounts) / (ca$size$x * ca$size$y) * 4
		}
	}
    
	updateFitnessMatrix()
    
	for (n in 1:length(ca$populations))
	{
		dbg('caUpdate 1')
		
		newgrid = matrix(ncol=ncol(ca$populations[[n]]$curgrid), nrow=nrow(ca$populations[[n]]$curgrid), byrow=TRUE)
		
		dbg('caUpdate 2')
		for (i in 1:ncol(ca$populations[[n]]$curgrid))
		{
			for (j in 1:nrow(ca$populations[[n]]$curgrid))
			{
				newgrid[i, j] = updateCell(i, j, ca$populations[[n]], globalfitness)
			}
		}
		dbg('caUpdate 3')
		for (i in 1:ncol(ca$populations[[n]]$curgrid))
		{
			for (j in 1:nrow(ca$populations[[n]]$curgrid))
			{
				ca$populations[[n]]$oldgrid[i, j] <<- ca$populations[[n]]$curgrid[i, j]
				ca$populations[[n]]$curgrid[i, j] <<- newgrid[i, j]
			}
		}
		#ca$oldgrid <<- ca$curgrid
		#ca$curgrid <<- newgrid
	}
	
    updateFitnessMatrix()
    writeLine()
	
	dbg('end caUpdate')
    rp.tkrreplot(ca$panel, tkrp)
    panel
}

highestFitness <- function(i, j, pop, globalfitness)
{
	roll = runif(1, 0.0, 100.0)
	if (roll < ca$globalEffect)
	{
		cell = pop$curgrid[i, j]
		fittest = 1
		for (i in 1:pop$players)
		{
			if (globalfitness[[pop$index]][i] > globalfitness[[pop$index]][fittest])
				fittest = i
		}
		return(fittest)
	}
	else
	{
		neighbors = getNeighborCoords(i, j, pop)
		fittest = 1
		for (i in 1:length(neighbors))
		{
			if (pop$fitnessMatrix[neighbors[[i]]$x, neighbors[[i]]$y] > pop$fitnessMatrix[neighbors[[fittest]]$x, neighbors[[fittest]]$y])
				fittest = i
		}
		return(pop$curgrid[neighbors[[fittest]]$x, neighbors[[fittest]]$y])
	}
	
	# fitnesses = c(pop$fitnessMatrix[neighbors[[fittest]]$x, neighbors[[fittest]]$y])
	# for (i in 1:length(neighbors))
    # {
        # if (i != fittest && pop$fitnessMatrix[neighbors[[i]]$x, neighbors[[i]]$y] == pop$fitnessMatrix[neighbors[[fittest]]$x, neighbors[[fittest]]$y])
			# fitnesses[length(fitnesses) + 1] = pop$fitnessMatrix[neighbors[[i]]$x, neighbors[[i]]$y]
    # }
	# return(fitnesses[pickProbability(fitnesses)])
}

exponentialFitness <- function(i, j, pop, globalfitness)
{
	roll = runif(1, 0.0, 100.0)
	if (roll < ca$globalEffect)
	{
		fitnesses = numeric(length=pop$players)
		for (i in 1:pop$players)
		{
			fitnesses[i] = exp(ca$b * globalfitness[[pop$index]][i])
		}
		return(pickProbability(fitnesses))
	}
	else
	{
		neighbors = getNeighborCoords(i, j, pop)
		fitnesses = numeric(length=pop$players)
		for (i in 1:length(neighbors))
		{
			fitnesses[pop$curgrid[neighbors[[i]]$x, neighbors[[i]]$y]] = fitnesses[pop$curgrid[neighbors[[i]]$x, neighbors[[i]]$y]] + exp(ca$b * pop$fitnessMatrix[neighbors[[i]]$x, neighbors[[i]]$y])
		}
		return(pickProbability(fitnesses))
	}
}

pickProbability <- function(somelist)
{
    total = sum(somelist)
    roll = runif(1, 0, total)

    for (i in 1:length(somelist))
    {
        if (roll < somelist[i])
        {
            return(i)
        }
        else
            roll = roll - somelist[i]
    }
    #ERROR
	return(0)
}

updateCell <- function(i, j, pop, globalfitness)
{
	if (pop$curgrid[i, j] > 0 && runif(1, 0, 1) < ca$deathrate)
		pop$curgrid[i, j] = 0 #death
	if (pop$curgrid[i, j] == 0 && runif(1, 0, 1) < ca$birthrate)
	{
		#print(ca$panel$fillType)
		if (ca$panel$fillType == ca$fillTypes[1])
		{
			return(highestFitness(i, j, pop, globalfitness))
		}
		else if (ca$panel$fillType == ca$fillTypes[2])
		{
			return(exponentialFitness(i, j, pop, globalfitness))
		}
	}
	else
		return(pop$curgrid[i, j])
}

calculatePublicGoods <- function(pop, cell, neighborCounts) #specialized public goods game function for calculating fitness
{
	cooperators = neighborCounts[[pop$index]][1]
	defectors = neighborCounts[[pop$index]][2]
	loners = neighborCounts[[pop$index]][3]
	if (cell == 1) # cooperator
	{
		return(ca$publicGoodsCost * ca$publicGoodsR * cooperators / 4 - ca$publicGoodsCost)
	}
	if (cell == 2) # defector
	{
		return(ca$publicGoodsCost * ca$publicGoodsR * cooperators / 4)
	}
	if (cell == 3) #loner
	{
		return(ca$publicGoodsLoner)
	}
	print("ERROR in calculatePublicGoods")
	mPrint("cell", cell)
	return(-99999) #ERROR
}

calculateFitnessSimple <- function(pop, cell, neighborCounts) # get fitness value for a non-zero square
{
	fitness = 0
	for (i in 1:length(neighborCounts[[pop$index]]))
	{
		#print(cell)
		#print(i)
		fitness = fitness + pop$payoffMatrix[cell, i] * neighborCounts[[pop$index]][i]	
	}
	return(fitness)
}

calculateBuyersSellers <- function(pop, cell, neighborCounts)
{
	fitness = 0
	for (i in 1:length(neighborCounts[[1 - pop$index + 2]])) #switch 1 to 2 and 2 to 1
	{
		fitness = fitness + pop$payoffMatrix[cell, i] * neighborCounts[[1 - pop$index + 2]][i]	
	}
	return(fitness)
}

countNeighbors <- function(neighbors, pop) # takes in output from getNeighborCoords
{
	pops = list()
	for (n in 1:length(ca$populations))
	{
	
		pops[[n]] = numeric(length=ca$populations[[n]]$players)
		#mPrint("neighbors", neighbors)
		for (i in 0:length(neighbors))
		{
			if (i == 0)
			{
				# this is a hack to get a pop to consider i,j in other pops its neighbors
				if (n != pop$index)
				{
					coords = getWrappedCoords(neighbors[[1]]$x, neighbors[[1]]$y + 1, pop)
					value = ca$populations[[n]]$curgrid[coords$x, coords$y]

					if (value > 0)
						pops[[n]][value] = pops[[n]][value] + 1
				}
			}
			else
			{
				value = ca$populations[[n]]$curgrid[neighbors[[i]]$x, neighbors[[i]]$y]

				if (value > 0)
					pops[[n]][value] = pops[[n]][value] + 1
			}
		}
	}
	return(pops)
}

getNeighborCoords <- function(i, j, pop) # point i,j to get neighbors from
{
    # 8-neighbors
    #print(ca$panel)
	#do 4-way always
	neighbors <- list(getWrappedCoords(i, j - 1, pop), getWrappedCoords(i - 1, j, pop), getWrappedCoords(i, j + 1, pop), getWrappedCoords(i + 1, j, pop))
	
	if (ca$panel$neighborhood == "8-way") # add these extra in 8-way case
	{
		neighbors[[5]] <- getWrappedCoords(i - 1, j - 1, pop)
		neighbors[[6]] <- getWrappedCoords(i - 1, j + 1, pop)
		neighbors[[7]] <- getWrappedCoords(i + 1, j + 1, pop)
		neighbors[[8]] <- getWrappedCoords(i + 1, j - 1, pop)
	}
	#print(neighbors)
	return(neighbors)
}

getWrappedCoords <- function(i, j, pop) # wrap i,j to a torus space
{
	p = list()
	p$x <- i
	p$y <- j
    if (i == 0)
		p$x = ncol(pop$curgrid)
	else if (i == ncol(pop$curgrid) + 1)
		p$x = 1
	if (j == 0)
		p$y = nrow(pop$curgrid)
	else if (j == nrow(pop$curgrid) + 1)
		p$y = 1
	return(p)
}

updateNeighbors <- function(panel)
{
	ca$panel$neighborhood <<- panel$neighborhood
	#print(ca$panel)
	panel
}

mPrint <- function(name, something)
{
	print(name)
	print(something)
}

updateGameType <- function()
{
	#caStop(ca$panel) TODO MOVE THIS SOMEWHERE!
	
	#print(ca$panel)
	if (ca$gameType == ca$gameMenuLabels[2])
	{
		makeCa(list(makePopulation(1, ca$size$x, ca$size$y, 3, calculateFitnessSimple, 
			matrix(
			c(0.5, 0, 1,
				1, 0.5, 0,
				0, 1, 0.5),
			nrow=3, ncol=3, byrow=TRUE), c("#FF0000", "#FFA500", "#0000FF"))))
            
        # for (i in 1:ca$size$x)
        # {
            # for (j in 1:ca$size$y)
            # {
                # if (i < 8)
                # {
                    # ca$populations[[1]]$curgrid[i, j] <<- 1
                # }
                # else if (i < 18)
                # {
                    # ca$populations[[1]]$curgrid[i, j] <<- 2
                # }
                # else
                # {
                    # ca$populations[[1]]$curgrid[i, j] <<- 3
                # }
            # }
        # }
	}
	else if (ca$gameType == ca$gameMenuLabels[3]) #pd
	{
		makeCa(list(makePopulation(1, ca$size$x, ca$size$y, 2, calculateFitnessSimple, 
			matrix(
			c(ca$bval - ca$cval, -ca$cval, ca$bval, 0),
			nrow=2, ncol=2, byrow=TRUE), c("#00FF00", "#FF0000"))))
	}
	else if (ca$gameType == ca$gameMenuLabels[4]) #snowdrift
	{
		makeCa(list(makePopulation(1, ca$size$x, ca$size$y, 2, calculateFitnessSimple, 
			matrix(
			c(ca$bval - ca$cval / 2, ca$bval - ca$cval, ca$bval, 0),
			nrow=2, ncol=2, byrow=TRUE), c("#00FF00", "#FF0000"))))
	}
	else if (ca$gameType == ca$gameMenuLabels[5])
	{
        if (ca$pubGoodsMode == 2)
        {
            makeCa(list(makePopulation(1, ca$size$x, ca$size$y, 2, calculatePublicGoods, NULL, c("#00FF00", "#FF0000"))))
        }
        else if (ca$pubGoodsMode == 3)
        {
            makeCa(list(makePopulation(1, ca$size$x, ca$size$y, 3, calculatePublicGoods, NULL, c("#00FF00", "#FF0000", "#FFA500"))))
        }
	}
	else if (ca$gameType == ca$gameMenuLabels[6])
	{
		print(ca$buyerSellerFitness)
		makeCa(list(
			makePopulation(1, ca$size$x, ca$size$y, 2, calculateBuyersSellers, #buyers
				matrix(c(ca$buyerSellerFitness[1], ca$buyerSellerFitness[3], ca$buyerSellerFitness[2], ca$buyerSellerFitness[4]), nrow=2, ncol=2, byrow=TRUE), c("#FFA500", "#0000FF")),
			makePopulation(2, ca$size$x, ca$size$y, 2, calculateBuyersSellers, #sellers
				matrix(c(ca$buyerSellerFitness[5], ca$buyerSellerFitness[7], ca$buyerSellerFitness[6], ca$buyerSellerFitness[8]), nrow=2, ncol=2, byrow=TRUE), c("#FF0000", "#00FF00"))))
	}
    else if (ca$gameType == ca$gameMenuLabels[7])
    {
        makeCa(list(makePopulation(1, ca$size$x, ca$size$y, ca$nstrats, calculateFitnessSimple, 
            ca$customPayoff, c("#FF0000", "#FFA500", "#0000FF", "#00FF00"))))
    }
    if (ca$initialized)
        rp.tkrreplot(ca$panel, tkrp)
}

updateFillType <- function(panel)
{
	ca$panel$fillType <<- panel$fillType
	#print(ca$panel)
	panel
}

updateGlobalEffect <- function(panel)
{
	ca$globalEffect <<- panel$globalEffect
	panel
}

updateDrawStyle <- function(panel)
{
	ca$panel$drawStyle <<- panel$drawStyle
	rp.tkrreplot(ca$panel, tkrp)
	panel
}

refresh <- function(panel)
{
    #update all prefs
    panel
}

timerUpdateNormal <- function(panel)
{
	if (ca$running)
	{
		caUpdate(ca$panel)
	}
	panel
}

returnTrue <- function(panel)
{
	TRUE
}

makeCa <- function(populationsList)
{
	ca$tick <<- 0
	ca$populations <<- populationsList
	makeFile()
}

makePopulation <- function(index, width, height, players, fitnessFunction, payoffMatrix, popColors)
{
	p = list()
	p$popColors = popColors
	p$index = index
	p$width = width
	p$height = height
	p$players = players
	p$curgrid = makeGrid(width, height, players)
	p$oldgrid = makeGrid(width, height, 0)
	p$fitnessFunction = fitnessFunction
	p$payoffMatrix = payoffMatrix
    #print(p)
	return(p)
}

dbg <- function(s)
{
	if (debugflag == TRUE)
		print(s)
}

makeGrid <- function(width, height, players) # grid constructor, initializes random grid
{
    g = matrix(ncol=width,nrow=height,byrow=TRUE)
    for (i in 1:width)
    {
        for (j in 1:height)
        {
			if (players == 0)
				g[i, j] = 0
			else
				g[i, j] <- sample(0:players, 1)
        }
    }
    return(g)
}

makeCustomChooseScreen <- function()
{
	customchoose = rp.control(title="Choose CA type")
	poplbls = c("1", "2")
	stratlbls = c("2", "3", "4")
	rp.radiogroup(customchoose, pops, poplbls, labels=poplbls, action=customOpt, title="Number of populations", initval=poplbls[1])
	rp.radiogroup(customchoose, strats, stratlbls, labels=stratlbls, action=customOpt, title="Number of strategies per population", initval=stratlbls[1])
	rp.button(customchoose, action = customDone, title = "Continue", quitbutton=TRUE)
}

customOpt <- function(panel)
{
	panel
}

customDone <- function(panel)
{
    pops = as.integer(panel$pops)
    if (pops > 1)
        ca$iter <<- 4
	makeMatrixSetupPanel(pops, as.integer(panel$strats))
	panel
}

menuClick <- function(panel)
{
	if (panel$menubtn == ca$gameMenuLabels[7]) #custom
	{
		ca$gameType <<- panel$menubtn
		makeCustomChooseScreen()
	}
	else if (panel$menubtn == ca$gameMenuLabels[2]) #rps
	{
		ca$gameType <<- panel$menubtn
		makeSetupPanel()
	}
	else if (panel$menubtn == ca$gameMenuLabels[3] || panel$menubtn == ca$gameMenuLabels[4]) #pd/sd
	{
		ca$gameType <<- panel$menubtn
		makeBCSetupPanel("bc")
	}
    else if (panel$menubtn == ca$gameMenuLabels[5]) #public goods
    {
		ca$gameType <<- panel$menubtn
		publicGoodsSetupPanel()
    }
	else if (panel$menubtn == ca$gameMenuLabels[6]) #buyers sellers
    {
		ca$gameType <<- panel$menubtn
		buyersSellersSetupPanel()
    }
	else if (panel$menubtn == ca$viewMenuLabels[2]) #toggle fitness view
	{
		ca$viewFitness <<- !(ca$viewFitness)
		rp.tkrreplot(ca$panel, tkrp)
	}
    else if (panel$menubtn == ca$viewMenuLabels[3]) # toggle greyscale
    {
        ca$bwOnly <<- !(ca$bwOnly)
        rp.tkrreplot(ca$panel, tkrp)
    }
}

makeSetupPanel <- function()
{
	setup = rp.control(title="Choose game options")
	
	rp.catitle = rp.text(setup, "CA Dimensions", name="catitle")
	rp.textentry(setup, size, title="dimensions", labels=c("width", "height"), keydown=TRUE, width=5, initval=c(ca$size$x, ca$size$y))
	rp.radiogroup(setup, worldType, ca$worldTypes, labels=ca$worldTypes, title = "Surface type", initval=ca$worldTypes[1])
	
	rp.text(setup, "CA parameters")
	rp.radiogroup(setup, fillType, ca$fillTypes, labels=ca$fillTypes, action = updateFillType, title = "Update type", initval=ca$fillTypes[1])
	rp.textentry(setup, betaRate, labels="beta (learning rate)", width = 5, keydown=TRUE, initval = ca$b)
	rp.textentry(setup, birthRate, labels="birth rate", width = 5, keydown=TRUE, initval=ca$birthrate)
	rp.textentry(setup, deathRate, labels="death rate", width = 5, keydown=TRUE, initval=ca$deathrate)
	
	rp.radiogroup(setup, neighborhood, ca$neighborhoodTypes, labels=ca$neighborhoodTypes, action = updateNeighbors, title = "Neighborhood size", initval=ca$neighborhoodTypes[1])
	
	rp.slider(setup, globalEffect, 0, 100, action=updateGlobalEffect, resolution = 5, showvalue = TRUE, title = "Global fitness", initval=ca$globalEffect)
	
	rp.textentry(setup, outputFile, labels="Output file", width=25, keydown=TRUE, initval = "-", name="filenameText")
	#rp.button(setup, action=saveDialog, title="Browse...", pos="right")
	
	rp.button(setup, action=setupComplete, title="Create", quitbutton=TRUE)
}

saveDialog <- function(panel)
{
	panel$outputFile <- tclvalue(tkgetSaveFile(initialfile = "foo.csv", filetypes = "{{CSV Files} {.csv}}"))
	print(panel)
	print(panel$catitle)
	panel$catitle$text = "it worked!"
	panel
}

setupOpt <- function(panel)
{
	#panel$size = panel$size
	panel
}

setupComplete <- function(panel)
{
	rp.control.put(panel$panelname, panel)
	#print(panel)
	ca$neighborhoodStyle <<- panel$neighborhood
	ca$fillType <<- panel$fillType
	ca$globalEffect <<- panel$globalEffect
	ca$size$x <<- as.integer(panel$size[[1]])
	ca$size$y <<- as.integer(panel$size[[2]])
	ca$world <<- panel$worldType
	#print(ca$size)
	updateGameType()
	ca$birthrate <<- as.numeric(panel$birthRate)
	ca$deathrate <<- as.numeric(panel$deathRate)
	ca$b <<- as.numeric(panel$betaRate)
	
	panel
}

makeBCSetupPanel <- function(style)
{
	setup = rp.control(title="Choose game parameters", style=style)
	#setup$style = style

	if (style == "bc")
		rp.textentry(setup, vals, labels=c("b  ", "c  "), title="game parameters", pos="left", width=5, keydown=TRUE, initval=c(ca$bval, ca$cval))
		
	rp.button(setup, action=bcDone, title="Continue", quitbutton=TRUE)
}

bcDone <- function(panel)
{
	rp.control.put(panel$panelname, panel)
	if (panel$style == "bc")
	{
		ca$bval <<- as.numeric(panel$vals[[1]])
		ca$cval <<- as.numeric(panel$vals[[2]])
	}
	makeSetupPanel()
	panel
}

publicGoodsSetupPanel <- function()
{
	setup = rp.control(title="Choose game parameters")

	lbls = c("2", "3")
	rp.radiogroup(setup, pubgoodsmode, lbls, labels=lbls, action=customOpt, title="Number of populations", initval=lbls[1])
	rp.button(setup, action=publicGoodsDone, title="Continue", quitbutton=TRUE)
}

publicGoodsSetup2 <- function()
{
	setup2 = rp.control(title="Choose game parameters")

	if (ca$pubGoodsMode == 2)
	{
		rp.textentry(setup2, vals, labels=c("R  ", "Cost  "), title="game parameters", pos="left", width=5, keydown=TRUE, initval=c(ca$publicGoodsR, ca$publicGoodsCost))
	}
	else
	{
		rp.textentry(setup2, vals, labels=c("R  ", "Cost  ", "LonerPayout  "), title="game parameters", pos="left", width=5, keydown=TRUE, initval=c(ca$publicGoodsR, ca$publicGoodsCost, ca$publicGoodsLoner))
	}

	rp.button(setup2, action=publicGoodsDone2, title="Continue", quitbutton=TRUE)
}

publicGoodsDone <- function(panel)
{
	rp.control.put(panel$panelname, panel)
    ca$pubGoodsMode <<- as.integer(panel$pubgoodsmode)
	publicGoodsSetup2()
	panel
}

publicGoodsDone2 <- function(panel)
{
	rp.control.put(panel$panelname, panel)
    
	if (ca$pubGoodsMode == 2)
	{
		ca$publicGoodsR <<- as.numeric(panel$vals[[1]])
		ca$publicGoodsCost <<- as.numeric(panel$vals[[2]])
	}
	else
	{
		ca$publicGoodsR <<- as.numeric(panel$vals[[1]])
		ca$publicGoodsCost <<- as.numeric(panel$vals[[2]])
		ca$publicGoodsLoner <<- as.numeric(panel$vals[[3]])
	}
	
	makeSetupPanel()
	panel
}

buyersSellersSetupPanel <- function()
{
	setup = rp.control(title="Choose payoff matrix")

	rp.text(setup, "Payoff Matrix for buyers and sellers")
		
	rp.text(setup, "Buyers vs Sellers")
	rp.text(setup, "         s_cheat    s_honest")
	rp.textentry(setup, payoffs1, labels=c("b_test  ", "b_lazy  "), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c(ca$buyerSellerFitness[1], ca$buyerSellerFitness[3]))
	rp.textentry(setup, payoffs2, labels=c("", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c(ca$buyerSellerFitness[2], ca$buyerSellerFitness[4]))
	
	rp.text(setup, "Sellers vs Buyers")
	rp.text(setup, "                    b_test     b_lazy")
	rp.textentry(setup, payoffs3, labels=c("s_cheat  ", "s_honest  "), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c(ca$buyerSellerFitness[5], ca$buyerSellerFitness[7]))
	rp.textentry(setup, payoffs4, labels=c("", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c(ca$buyerSellerFitness[6], ca$buyerSellerFitness[8]))
	
	rp.button(setup, action = buyersSellersSetupDone, title = "Continue", pos="bottom", quitbutton=TRUE)
}

buyersSellersSetupDone <- function(panel)
{   
	ca$buyerSellerFitness[1] <<- as.numeric(panel$payoffs1[[1]])
	ca$buyerSellerFitness[2] <<- as.numeric(panel$payoffs2[[1]])
	ca$buyerSellerFitness[3] <<- as.numeric(panel$payoffs1[[2]])
	ca$buyerSellerFitness[4] <<- as.numeric(panel$payoffs2[[2]])
	
	ca$buyerSellerFitness[5] <<- as.numeric(panel$payoffs3[[1]])
	ca$buyerSellerFitness[6] <<- as.numeric(panel$payoffs4[[1]])
	ca$buyerSellerFitness[7] <<- as.numeric(panel$payoffs3[[2]])
	ca$buyerSellerFitness[8] <<- as.numeric(panel$payoffs4[[2]])
		
	makeSetupPanel()
	panel
}

makeMatrixSetupPanel <- function(npops, nstrats)
{
	setup = rp.control(title="Choose payoff matrix", nstrats=nstrats, npops=npops)

	rp.text(setup, "Payoff Matrix")
	if (npops > 1)
		rp.text(setup, paste("Fitness matrix for pop ", floor(iter / 2), " vs ", iter %% 2))
	
	# causing crashes???
	#rp.grid(setup, pos=list(row=1, column=1), background="navy", name="p1")
	#rp.grid(setup, pos=list(row=2, column=1), background="green", name="p2")
	
	if (nstrats == 2)
	{
		rp.text(setup, "      1              2")
		rp.textentry(setup, payoffs1, labels=c("1  ", "2  "), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0"))
		rp.textentry(setup, payoffs2, labels=c("", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0"))
	}
	else if (nstrats == 3)
	{
		rp.text(setup, "      1              2              3")
		rp.textentry(setup, payoffs1, labels=c("1  ", "2  ", "3  "), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0"))
		rp.textentry(setup, payoffs2, labels=c("", "", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0"))
		rp.textentry(setup, payoffs3, labels=c("", "", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0"))
	}
	else if (nstrats == 4)
	{
		rp.text(setup, "      1              2              3              4")
		rp.textentry(setup, payoffs1, labels=c("1  ", "2  ", "3  ", "4  "), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0", "0"))
		rp.textentry(setup, payoffs2, labels=c("", "", "", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0", "0"))
		rp.textentry(setup, payoffs3, labels=c("", "", "", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0", "0"))
		rp.textentry(setup, payoffs4, labels=c("", "", "", ""), title="payoff matrix", pos="left", width=5, keydown=TRUE, initval=c("0", "0", "0", "0"))
	}
	rp.button(setup, action = setupDone, title = "Continue", quitbutton=TRUE)
}

setupOpt <- function(panel)
{
	panel
}

setupDone <- function(panel)
{
	#print('panelpayoffs1')
	#print(panel$payoffs1)
	#print('firstelem')
	#print(panel$payoffs1[1])
	#print('deref')
	#print(panel$payoffs1[[1]])
    
    if (panel$nstrats == 2)
	{
        ca$customPayoff <<- matrix(c(panel$payoffs1[[1]], panel$payoffs1[[2]], panel$payoffs2[[1]], panel$payoffs2[[2]]), nrow=2, ncol=2, byrow=TRUE)
	}
	else if (panel$nstrats == 3)
	{
        ca$customPayoff <<- matrix(c(panel$payoffs1[[1]], panel$payoffs1[[2]], panel$payoffs1[[3]],
            panel$payoffs2[[1]], panel$payoffs2[[2]], panel$payoffs2[[3]],
            panel$payoffs3[[1]], panel$payoffs3[[2]], panel$payoffs3[[3]]), nrow=3, ncol=3, byrow=TRUE)
	}
	else if (panel$nstrats == 4)
	{
		ca$customPayoff <<- matrix(c(panel$payoffs1[[1]], panel$payoffs1[[2]], panel$payoffs1[[3]], panel$payoffs1[[4]],
            panel$payoffs2[[1]], panel$payoffs2[[2]], panel$payoffs2[[3]], panel$payoffs2[[4]],
            panel$payoffs3[[1]], panel$payoffs3[[2]], panel$payoffs3[[3]], panel$payoffs3[[4]],
            panel$payoffs4[[1]], panel$payoffs4[[2]], panel$payoffs4[[3]], panel$payoffs4[[4]]), nrow=4, ncol=4, byrow=TRUE)
	}
    ca$nstrats <<- panel$nstrats
    
    if (panel$npops > 1 && ca$iter > 1)
    {
        ca$iter <<- ca$iter - 1
        makeMatrixSetupPanel(panel$nstrats, panel$npops)
    }
    else
    {
        makeSetupPanel()
    }
	panel
}

init <- function()
{
	ca <<- list()
    ca$initialized <<- FALSE
	ca$running <<- FALSE
	ca$tick <<- 0
	ca$panel <<- rp.control(title="Cellular automata simulation", size=c(500, 572, 100, 100))
	ca$gameMenuLabels <<- list("Create", "Rock Paper Scissors", "Prisoners Dilemma", "Snowdrift", "Public Goods", "Buyers Sellers", "Custom")
	ca$viewMenuLabels <<- list("View", "Cell fitness", "Toggle greyscale")
	ca$menuLabels <<- list(ca$gameMenuLabels, ca$viewMenuLabels)
	ca$gameType <<- ca$gameMenuLabels[[2]]
	ca$panel$menubtn <<- ""
	rp.menu(ca$panel, menubtn, labels=ca$menuLabels, action=menuClick)
	rp.grid(ca$panel, pos=list(row=1, column=1, width=500, height=500), name="g1")
	#rp.grid(ca$panel, pos=list(row=1, column=2, width=200, height=600), background="green", name="g2")
	
	ca$viewFitness <<- FALSE
	ca$bwOnly <<- FALSE
    
	ca$bval <<- 1
	ca$cval <<- 1
    ca$pubGoodsMode <<- 2
	ca$publicGoodsR <<- 2
	ca$publicGoodsCost <<- 6
	ca$publicGoodsLoner <<- 3
	
	ca$buyerSellerFitness <<- c(1, -1, -4, 2, -8, 0, 2, 1)
	
	ca$neighborhoodTypes <<- c("4-way", "8-way")
	ca$worldTypes <<- c("Torus", "Flat")
	ca$world <<- ca$worldTypes[1]
	#rp.radiogroup(ca$panel, neighborhood, ca$neighborhoodTypes, labels=ca$neighborhoodTypes, action = updateNeighbors, title = "Neighborhood", parentname="g2", initval=ca$neighborhoodTypes[1])
	ca$panel$neighborhood <<- ca$neighborhoodTypes[1]
	
	ca$fillTypes <<- c("deterministic", "stochastic")
	#rp.radiogroup(ca$panel, fillType, ca$fillTypes, labels=ca$fillTypes, action = updateFillType, title = "Update type", parentname="g2", initval=ca$fillTypes[1])
	ca$panel$fillType <<- ca$fillTypes[1]
    
    #ca$globaleffectTypes <<- c("0%", "20%")
	#rp.radiogroup(ca$panel, globalEffect, ca$globaleffectTypes, labels=ca$globaleffectTypes, action = updateGlobalEffect, title = "Global effects", parentname="g2", initval=ca$globaleffectTypes[1])
	#ca$panel$globalEffect <<- ca$globaleffectTypes[1]
	
	#rp.slider(ca$panel, globalEffect, 0, 100, action=updateGlobalEffect, parentname="g2", resolution = 5, showvalue = TRUE, title = "Global effects", initval=0)
	ca$globalEffect <<- 0
	
	ca$drawStyles <<- c("Centered", "Corner")
	#rp.radiogroup(ca$panel, drawStyle, ca$drawStyles, labels=ca$drawStyles, action = updateDrawStyle, title = "Draw style", parentname="g2", initval=ca$drawStyles[1])
	ca$panel$drawStyle <<- ca$drawStyles[1]
	
	
	#rp.listbox(ca$panel, game, ca$gameTypeLabels, labels=ca$gameTypeLabels, rows=length(ca$gameTypeLabels), initval=ca$gameTypeLabels[1], action=updateGameType, parentname="g2", sleep = 0.01)
	#ca$panel$game <<- ca$gameTypeLabels[1]
	
	rp.grid(ca$panel, pos=list(row=2, column=1, width=500, height=1), background="black", name="sep")
	rp.grid(ca$panel, pos=list(row=3, column=1, width=500, height=54), name="btnrow")
	
	rp.text(ca$panel, " ", pos="right", parentname="btnrow")
	rp.button(ca$panel, action = caStop, title = "Stop", pos="right", parentname="btnrow")
	rp.button(ca$panel, action = caStart, title = "Start", pos="right", parentname="btnrow")
	rp.button(ca$panel, action = caUpdate, title = "Update", pos="right", parentname="btnrow")
	#rp.text(ca$panel, "Rock Paper Scissors - generation #", name="statustext", pos="left", parentname="btnrow")
	rp.timer(ca$panel, 2000, timerUpdateNormal, returnTrue)
	
	ca$size <<- list(x=25, y=25)
	updateGameType()
	updateNeighbors(ca$panel)
	ca$deathrate <<- 0.125
	ca$birthrate <<- 1
	ca$b <<- 0.5
    ca$bwColors <<- c("#FFFFFF", "#000000", "#AAAAAA", "#555555")
    ca$iter <<- 0
	
	rp.tkrplot(ca$panel, tkrp, drawGrid, pos=list(row=1, column=1, width=500, height=500), hscale=1.39, vscale=1.39, parentname="g1")
	ca$initialized <<- TRUE
}

init()


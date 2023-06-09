# our developed model is executed by running the following function "simcode".
# in our study, parameters held constant were specified within the function,
# but the varied parameters were specified and supplied to the function like: simcode(cover=0.1,...,). 

library(raster)
simcode <- function(
  # landscape stuff               
  cover, # habitat amount (or proportion: HP)
  frag,  # habitat fragmentation (FRAG)
  # background (low-quality) matrix stuff
  matrix.move.prob, # movement probability (permeability) of low-quality matrix
  matrix.move.survival, # movement survival of low-quality matrix
  # high quality matrix stuff  
  hqm_cover, # proportion of high-quality matrix
  hqm_frag,  # fragmentation of high-quality matrix
  improv.growth.rate, # productivity of high-quality matrix
  improv.move.prob,   # permeability of high-quality matrix
  improv.move.survival.unsettle     # when high-quality (improved) matrix works as 'dispersal passage' and 'population sink', this variable dictates movement survival and partial settlement into improved matrix, respectively
){
  # set-up parameters
  nyear              = 2000         # number of simulation years
  size               = 50           # length of one side (in cell number)
  # demographic & movement stuff of habitat
  growth.rate        = 0.1          # reproductive rate of habitat
  mortal.rate        = 0.05         # annual mortality of habitat
  emigr.rate         = 0.1          # emigration rate (whether individual depart from the habitat cell) 
  habita.move.prob   = 1            # movement probability (permeability) of habitat cell
  # demographic & movement stuff of high quality matrix
  improv.mortal.rate = 0.1			# annual mortality of high-quality matrix when it supports population sink
  improv.emigr.rate  = 0.5			# emigration rate of high-quality matrix when it supports population sink
  # converting survival to mortality
  matrix.move.mortality        <- 1 - matrix.move.survival
  improv.move.mortality.settle <- 1 - improv.move.survival.unsettle
  
  # other parameters of focal landscape
  ncol  <- nrow <- size             # landscape size
  nhabi <- cover*ncol*nrow          # number of habitat cell
  nhqm  <- hqm_cover*ncol*nrow      # number of high-quality matrix cell
  
  # landscape generation (for habitat)
  # cover (habitat) type code
  hcode <- c(0,10,100) # background low-quality matrix, high-quality matrix, and (original) habitat
  
  # landscape: one additional row and column are added for a practical reason (will be deleted finally)
  land <- matrix(0,ncol+1,nrow+1)
  
  while(sum(land)<nhabi*hcode[3]){ # original habitat will continue to be generated until its total area is met
    # select random cell by referring their coordinates
    # coordinates are designed to be between 1 to ncol or nrow ('zero' is not included)
    # note: since as.integer make non-integral values truncated towards zero, '+1' is required: check -> table(as.integer(runif(1000,min=1,max=(10+1))))
    x.land <- as.integer(runif(1,min=1,max=(ncol+1))) 
    y.land <- as.integer(runif(1,min=1,max=(nrow+1)))
    if(land[x.land,y.land]==hcode[1]){ # if the selected cell is already designated as habitat, do nothing (= else{})
      if(runif(1) < frag){ # habitat configuration: following Fahrig's methodology
        land[x.land,y.land] <- hcode[3]
      } else{
        # examine surrounding cells: 'zero' row and column numbers can be specified here (but row/column+1 cannot be specified)
        surround <- land[c(x.land-1,x.land,x.land+1),c(y.land-1,y.land,y.land+1)]
        if(sum(surround) >= hcode[3]){
          land[x.land,y.land] <- hcode[3]
        } else{}
      }
      # print(sum(land)/hcode[3])
    } else{}
  }
  
  land <- land[-c(nrow+1),-c(ncol+1)]
  land <- raster(land)
  # plot(land,las=1) # check the status visually if needed
  
  # write & read
  writeRaster(land,paste0("land_",cover,"_",frag,"_.asc"),overwrite=TRUE)
  land <- raster(paste0("land_",cover,"_",frag,"_.asc"))
  # plot(land,las=1,main=paste0("HABITAT = ",cover,", FRAG = ",frag))
  land <- as.matrix(land)
  
  #############################
  # landscape generation (for high-quality [improved] matrix)
  # adding one row and column with 'zero' values for the practical reason
  land <- rbind(land,rep(0,dim(land)[1]))
  land <- cbind(land,rep(0,dim(land)[1]))
  
  while(sum(land)<(nhabi*hcode[3]+nhqm*hcode[2])){ # high quality matrix will continue to be generated until its total area is met
    x.land <- as.integer(runif(1,min=1,max=(ncol+1))) 
    y.land <- as.integer(runif(1,min=1,max=(nrow+1)))
    if(land[x.land,y.land]==hcode[1]){
      if(runif(1) < hqm_frag){
        land[x.land,y.land] <- hcode[2]
      } else{
        surround <- land[c(x.land-1,x.land,x.land+1),c(y.land-1,y.land,y.land+1)]
        if(sum(surround) >= hcode[2]){ # the case when original habitat (hcode[3] = 100) exists nearby is included
          land[x.land,y.land] <- hcode[2]
        } else{}
      }
      # print((sum(land)-nhabi*hcode[3])/hcode[2])
    } else{}
  }
  
  land <- land[-c(nrow+1),-c(ncol+1)]
  table(land) # check whether the number of every habitat type is equal to the intended number
  land <- raster(land)
  #  plot(land,las=1,main=paste0("HP=",cover,", FRAG=",frag,", HQM=",hqm_cover,", HQMF=",hqm_frag))
  # write (& read)
  writeRaster(land,paste0("land_",cover,"_",frag,"_",hqm_cover,"_",hqm_frag,"_.asc"), overwrite=TRUE)
  ############################# end of landscape generation  
  
  #############################
  # dispersal & population dynamics
  # given parameters
  ini.pop     <- 10     # initial population size for every habitat cell
  Ceiling     <- 10     # ceiling of local population size (per cell)
  # reading landscape
  land <- raster(paste0("land_",cover,"_",frag,"_",hqm_cover,"_",hqm_frag,"_.asc"))
  # plot(land,las=1,main=paste0("HABITAT = ",cover,", FRAG = ",frag,", HQM = ",hqm_cover,", HQMFRAG = ",hqm_frag))
  land <- as.matrix(land)
  
  #############################
  # dispersal
  # input var.
  matrix.move.resi <- 1/matrix.move.prob
  habita.move.resi <- 1/habita.move.prob
  improv.move.resi <- 1/improv.move.prob
  inj.curr  <- 100 # 100 ampere is injected to every habitat cell to assess possible movement destination

  # coordinates of original habitat & high-quality matrix
  table(land)
  x.coord    <- which(land==hcode[3],arr.ind=TRUE)[,"row"]
  y.coord    <- which(land==hcode[3],arr.ind=TRUE)[,"col"]
  hab.coord  <- rbind(x.coord,y.coord)
  im.x.coord <- which(land==hcode[2],arr.ind=TRUE)[,"row"]
  im.y.coord <- which(land==hcode[2],arr.ind=TRUE)[,"col"]
  im.coord   <- rbind(im.x.coord,im.y.coord)
  hab.im.x.coord <- which((land==hcode[3]|land==hcode[2]),arr.ind=TRUE)[,"row"]
  hab.im.y.coord <- which((land==hcode[3]|land==hcode[2]),arr.ind=TRUE)[,"col"]
  hab.im.coord <- rbind(hab.im.x.coord,hab.im.y.coord)
  im.mat.x.coord <- which((land==hcode[1]|land==hcode[2]),arr.ind=TRUE)[,"row"]
  im.mat.y.coord <- which((land==hcode[1]|land==hcode[2]),arr.ind=TRUE)[,"col"]
  im.mat.coord <- rbind(im.mat.x.coord,im.mat.y.coord)
  
  # movement resistance
  resist <- matrix(matrix.move.resi,nrow=nrow,ncol=ncol)
  resist[land==hcode[3]] <- habita.move.resi # original habitat
  resist[land==hcode[2]] <- improv.move.resi # high-quality matrix
  Resist <- resist
  resist <- raster(resist,xmn=1,ymn=1,xmx=ncol,ymx=nrow)
  #  plot(resist,las=1,main="Resistance")

  # ground resistor is used to designate possible final destination (habitat cells other than source cell), movement mortality and partial settlement
  # a common ground resistor grid is created; specific ground resistor for every iteration of the Circuitscape would be newly created by modifying the common one (in the for loop)
  common.ground <- matrix(NA,nrow,ncol)
  common.ground[Resist==habita.move.resi]            <- 0                # assign zero movement mortality into source habitat, which means direct ground connection (R = 0) (i.e., all individuals reaching the cells settle into the cells and reproduce there)
  # specify the values of exiting resistance by looking at surrounding movement resistance (since exiting flow depends not only on exiting resistance of the exiting cells but also movement resistance of the surrounding cells)
  # since the values of Resist[,51] and Resist[51,] (i.e., at margin) cannot be assessed, resistance matrix is expanded (but Resist[0,0] can be assessed)
  ExpandResist <- cbind(Resist,rep(NA,dim(Resist)[2]))
  ExpandResist <- rbind(ExpandResist,rep(NA,dim(ExpandResist)[2]))
  for(i in 1:dim(im.mat.coord)[2]){
    if(land[im.mat.coord[1,i],im.mat.coord[2,i]]==hcode[1]){ # movement mortality depends on whether the cell is the background low-quality matrix or improved high-quality matrix
      move.mortality <- matrix.move.mortality
    } else{
      move.mortality <- improv.move.mortality.settle
    }
    common.ground[im.mat.coord[1,i],im.mat.coord[2,i]] <- ((1-move.mortality)/move.mortality)/(
      # four cardinals
      sum(
          1/ExpandResist[im.mat.coord[1,i]-1,im.mat.coord[2,i]],
          1/ExpandResist[im.mat.coord[1,i],  im.mat.coord[2,i]-1],
          1/ExpandResist[im.mat.coord[1,i]+1,im.mat.coord[2,i]],
          1/ExpandResist[im.mat.coord[1,i],  im.mat.coord[2,i]+1],na.rm=T
      ) + 
        # then, four diagonals
        (1/sqrt(2))*(
          sum(
              1/ExpandResist[im.mat.coord[1,i]-1,im.mat.coord[2,i]-1],
              1/ExpandResist[im.mat.coord[1,i]+1,im.mat.coord[2,i]-1],
              1/ExpandResist[im.mat.coord[1,i]+1,im.mat.coord[2,i]+1],
              1/ExpandResist[im.mat.coord[1,i]-1,im.mat.coord[2,i]+1],na.rm=T
          )
        )
    )
  }
  # plot(raster(common.ground))

  if(improv.growth.rate==0){  # when high-quality matrix does not allow the reproduction, only dispersals from original habitat cells are considered below
    hab.im.coord <- hab.coord
  } else{}

  # from making source & ground rasters to running Circuitscape, and retrieve dispersal success
  # scenario 10 produces large numbers of output files, which can be larger than the upper limit of the number of files in the single folder
  # therefore, the series of the following procedure is to be divided into five parts to lower the number of files
  # first, divide the reproductive cells into five groups
  group  <- seq(1,dim(hab.im.coord)[2])
  group1 <- seq(1,c(length(group)/5))    # denominator 5 should work in typical cases since the cell size is 50 
  group2 <- seq(tail(group1,1)+1,c(length(group)/5*2))
  group3 <- seq(tail(group2,1)+1,c(length(group)/5*3))
  group4 <- seq(tail(group3,1)+1,c(length(group)/5*4))
  group5 <- seq(tail(group4,1)+1,c(length(group)/5*5))
  all(group==c(group1,group2,group3,group4,group5)) # check the division
  
  nsite <- dim(hab.im.coord)[2] # number of possible reproduction cells (including source habitat and population sink)
  # dispersal success (number of immigrants)
  disp.succ <- matrix(NA,nsite,nsite) # row = source (departure or injected) cell; column = recipient cell
  
  # since the same procedure would be repeated five times, the procedure is to be wrapped as a function
  RasterCircuitDisSucc <- function(GROUP){
    for(i in GROUP){ # when high-quality matrix does allow the reproduction, this includes high-quality matrix since it can support potential emigrants
      # current is injected into every reproductive cell to assess possible movement from it (injected cell is called 'source')
      source <- matrix(NA,nrow,ncol)
      source[hab.im.coord[1,i],hab.im.coord[2,i]] <- inj.curr
      source <- raster(source,xmn=1,ymn=1,xmx=ncol,ymx=nrow)
      # plot(source,las=1)
      # common ground is customized to consider that injected cell should have 'NA' ground resistance
      ground <- common.ground
      ground[hab.im.coord[1,i],hab.im.coord[2,i]] <- NA  # injected cell from which dispersers come (can comprise source & population sink) should have 'NA' resistance
      ground <- raster(ground,xmn=1,ymn=1,xmx=ncol,ymx=nrow)
      # plot(ground,las=1)
      # set projection (Circuitscape requires rasters to be projected)
      projection(resist) <- projection(source) <- projection(ground)  <- "+init=epsg:3100"
      # save raster files
      writeRaster(resist,"resistance.asc",overwrite=TRUE)
      writeRaster(source,paste0("source_",i,"_.asc"),overwrite=TRUE)
      writeRaster(ground,paste0("ground_",i,"_.asc"),overwrite=TRUE)
    } # produce the bunch of raster files for each group
    
    ##################
    # run Circuitscape
    for(i in GROUP){ # when high-quality matrix does allow the reproduction, this includes high-quality matrix since it can support potential emigrants
      # settings and options
      cs_ini <- c("[circuitscape options]","data_type = raster","scenario = advanced","write_cur_maps = 1","write_volt_maps = 1",
                  paste(c("ground_file =","source_file =","habitat_file =","output_file ="),
                        paste(getwd(),c(
                          paste0("ground_",i,"_.asc"),
                          paste0("source_",i,"_.asc"),
                          "resistance.asc",
                          paste0("flow_",i,".out")
                        ),sep="/")))
      writeLines(cs_ini,paste0("circuitini_",i,".ini"))
    }
    # Circuitscape is executed: specific commands depend on the computation environment
    # for Windows environment, something like as follows:
    cs_exe <- 'C:/"Program Files"/Circuitscape/cs_run.exe'
    for(i in GROUP){
      wd <- getwd()
      cs_run <- paste(cs_exe,paste(getwd(),paste0("circuitini_",i,".ini"),sep="/"))
      system(cs_run)
    }
    # end of Circuitscape
 
    # retrieve dispersal success
    # note: dispersal success of source habitat (with 'zero' ground resistance) and population sink (with 'nonzero' ground resistance) is assessed by current and voltage map, respectively 
    for(i in GROUP){
      cur <- raster(paste0("flow_",i,"_curmap.asc"))  # current will always take values >=0
      # plot(cur,las=1,main=paste0("Current (from ",i,")"))
      # text(cur)
      cur <- as.matrix(cur)
      vol <- raster(paste0("flow_",i,"_voltmap.asc")) # voltage will always take values >=0
      # plot(vol,las=1,main=paste0("Voltage (from ",i,")"))
      # text(vol)
      vol <- as.matrix(vol)
      grnd <- raster(paste0("ground_",i,"_.asc"))
      # plot(grnd,las=1,main=paste0("Ground (from ",i,")")) # ground will take values of NA (injected cell), 0 (possible recipient cell [source habitat]) or >=0 (low- or high-quality matrix)
      # text(grnd)
      grnd <- as.matrix(grnd)
      exitcur <- vol/grnd     # exiting current from ground resistor with 'non-zero' resistance
      # exitcur <- raster(exitcur) # NA (exiting resistance = NA), NaN (exiting resistance = 0 - original habitat)
      # plot(exitcur,las=1,main=paste0("Exiting current (from ",i,")"))
      # text(exitcur)
      
      for(j in 1:nsite){
        # if target cell is original habitat (hcode[3]) and high-quality matrix (hcode[2]), retrieve the values from entering current ('current' output file) and exiting current from ground resistor (voltage/ground), respectively 
        if(land[hab.im.coord[1,j],hab.im.coord[2,j]]==hcode[3]){
          disp.succ[i,j] <<- cur[hab.im.coord[1,j],hab.im.coord[2,j]]
        } else{
          if(land[hab.im.coord[1,j],hab.im.coord[2,j]]==hcode[2]){ # when HQ matrix does not and does allow reproduction, exiting current represents movement mortality and proportion of settling individuals, respectively
            disp.succ[i,j] <<- exitcur[hab.im.coord[1,j],hab.im.coord[2,j]]
          } else{}
        } 
      }
      disp.succ[i,i] <<- cur[hab.im.coord[1,i],hab.im.coord[2,i]] # diagonal element should be injected current (exiting current would yield 'NA' since ground target cell has 'NA' resistance)
    }  
    
    # here, source, ground, and flow files are to be deleted to reduce the number of output files
    # list-up all files
    files <- list.files()
    # list-up files with '_.asc' (ground & source asc files), 'curmap.asc', 'voltmap.asc', '.ini' or '.out'
    gs.files     <- grep("\\\\_.asc$",       files)
    cur.files    <- grep("\\\\curmap.asc$",  files)
    volt.files   <- grep("\\\\voltmap.asc$", files)
    ini.files    <- grep("\\\\.ini$",        files)
    out.files    <- grep("\\\\.out$",        files)
    # remove these files
    file.remove(files[gs.files])
    file.remove(files[cur.files])
    file.remove(files[volt.files])
    file.remove(files[ini.files])
    file.remove(files[out.files])
  }
  
  # repeat the above procedure
  RasterCircuitDisSucc(group1)
  RasterCircuitDisSucc(group2)
  RasterCircuitDisSucc(group3)
  RasterCircuitDisSucc(group4)
  RasterCircuitDisSucc(group5)

  # # making plots to check individual movements
  #  pdf(file="flow.pdf")
  #  plot(raster(land),las=1,main=paste0("HABITAT = ",cover,", FRAG = ",frag,", HQM = ",hqm_cover,", HQMFRAG = ",hqm_frag))
  #  plot(resist,las=1,main="Resistance")
  #  for(i in 1:nsite){
  #    cur <- raster(paste0("flow_",i,"_curmap.asc"));  cur <- as.matrix(cur); flow <- cur
  #    vol <- raster(paste0("flow_",i,"_voltmap.asc")); vol <- as.matrix(vol)
  #    grnd <- raster(paste0("ground_",i,"_.asc"));     grnd <- as.matrix(grnd)
  #    exitcur <- vol/grnd; exitcur <- as.matrix(exitcur)
  #    flow[Resist==matrix.move.resi] <- exitcur[Resist==matrix.move.resi]
  #    flow[Resist==improv.move.resi] <- exitcur[Resist==improv.move.resi]
  #    flow[hab.im.coord[1,i],hab.im.coord[2,i]] <- cur[hab.im.coord[1,i],hab.im.coord[2,i]] # because HQ matrix can yield NA value
  #    Nflow <- sum(flow,na.rm=T)
  #    plot(raster(flow),las=1,main=paste0("Current flow (from ",i,": sum = ",round(Nflow,0),")"))
  #    if(improv.growth.rate==0){ # when HQ matrix does not allow reproduction, settlement only occurs in original habitat
  #      settle <- flow; settle[Resist!=habita.move.resi] <- NA # excluding low & high quality matrix
  #      text(raster(settle))
  #      text(raster(exitcur),col="red")
  #    } else{ # when HQ matrix does allow reproduction, settlement occurs in original habitat as well as HQ matrix
  #      settle <- flow; settle[Resist==matrix.move.resi] <- NA # excluding low quality matrix
  #      text(raster(settle))
  #      exitcur[Resist==improv.move.resi] <- NA
  #      text(raster(exitcur),col="red")
  #    }
  #  }
  #  dev.off()
  
  # mathematical matrix dictating dispersal success from each OR habitat & HQ matrix cell to the other OR habitat & HQ matrix cells
  disp.mort <- 1 - (apply(disp.succ,1,sum)-inj.curr)/inj.curr
  disp.mort <- ifelse(disp.mort<0,0,disp.mort) # when improv.move.mortality.settle = 0, minus mortality (but ~0) can occur and the multinomial distribution fails. minus mortality is replaced by zero
  disp.succ <- disp.succ/inj.curr; diag(disp.succ) <- disp.mort
  disp.succ <- ifelse(disp.succ<0,0,disp.succ) # sink to sink can also yield minus mortality
  
  #############################
  # population dynamics
  # habitat type of each original habitat or high quality matrix cells
  hab.type <- rep(NA,nsite) 
  for(i in 1:nsite){
    hab.type[i] <- land[hab.im.coord[1,i],hab.im.coord[2,i]]
  }
  
  # number of individuals in each OR habitat & HQ matrix cell
  N <- matrix(NA,nsite,nyear) 
  N[hab.type==hcode[3],1] <- ini.pop # original habitat has individuals
  N[hab.type==hcode[2],1] <- ini.pop # high-quality matrix also supports individuals
  
  # organize demographic parameters
  hab.im.growth.rate <- rep(NA,nsite)
  hab.im.growth.rate[hab.type==hcode[3]] <- growth.rate        
  hab.im.growth.rate[hab.type==hcode[2]] <- improv.growth.rate
  hab.im.emigr.rate <- rep(NA,nsite)
  hab.im.emigr.rate[hab.type==hcode[3]] <- emigr.rate
  hab.im.emigr.rate[hab.type==hcode[2]] <- improv.emigr.rate
  hab.im.mortal.rate <- rep(NA,nsite)
  hab.im.mortal.rate[hab.type==hcode[3]] <- mortal.rate
  hab.im.mortal.rate[hab.type==hcode[2]] <- improv.mortal.rate
  
  for(i in 1:(nyear-1)){
    
    # reproduction
    expect_grow <- N[,i]*hab.im.growth.rate
    recruit <- rpois(nsite,expect_grow)
    Ncombine <- N[,i] + recruit
    
    # natural death (winter/annual mortality)
    Nsurvive <- rbinom(n=nsite,size=Ncombine,prob=c(1-hab.im.mortal.rate))
    
    # dispersal
    Nceiling <- Nsurvive
    Nceiling[Nceiling>Ceiling] <- Ceiling  # ceiling of each habitat cell
    Nsurplus <- Nsurvive - Nceiling        # density-dependent/active dispersal fraction (surplus)
    
    Ndisp <- rbinom(n=nsite,size=Nceiling,hab.im.emigr.rate)  # intrinsic/density-independent/passive dispersal
    Nremain <- Nceiling - Ndisp
    Ndisp <- Ndisp + Nsurplus  # total dispersal fraction
    
    # dispersal
    NsuccD <- matrix(NA,nsite,nsite)
    for(j in 1:nsite){
      NsuccD[j,] <- rmultinom(n=1,size=Ndisp[j],prob=disp.succ[j,])
    }
    NdeadD <- diag(NsuccD)
    diag(NsuccD) <- 0
    N[,i+1] <- Nremain + apply(NsuccD,2,sum)
    
    # overabundance mortality
    N[N[,i+1]>Ceiling,i+1] <- Ceiling
  }
  
  # pdf(file="dynamics.pdf")
  # plot(N[1,]~seq(1,nyear),ylim=c(0,max(N)),
  #      type="l",las=1,xlab="Year",ylab="N. of individuals in each patch",
  #      main="Population dynamics")
  # for(i in 2:nsite){
  #   lines(N[i,]~seq(1,nyear),col=i)
  # }
  # SumN <- apply(N,2,sum); SumN <- (max(N)/max(SumN))*SumN
  # lines(SumN~seq(1,nyear),lwd=3,col="grey")
  # legend("topright",lwd=c(1,3),lty=1,col=c("black","grey"),legend=c("Each","Total"))
  # dev.off()
  
  # removing raster and Circuitscape files to save the memory
  # list-up all files
  files <- list.files()
  # list-up files with '.asc'
  asc.files <- grep("\\\\.asc$", files)
  # remove these files
  file.remove(files[asc.files])
  
  # saving parameters
  # background parameters
  write.table(size,                          "size.csv",                          row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(cover,                         "cover.csv",                         row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(frag,                          "frag.csv",                          row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(hqm_cover,                     "hqm_cover.csv",                     row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(hqm_frag,                      "hqm_frag.csv",                      row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(matrix.move.prob,              "matrix.move.prob.csv",              row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(habita.move.prob,              "habita.move.prob.csv",              row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(improv.move.prob,              "improv.move.prob.csv",              row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(matrix.move.survival,          "matrix.move.survival.csv",          row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(improv.move.survival.unsettle, "improv.move.survival.unsettle.csv", row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(growth.rate,                   "growth.rate.csv",                   row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(emigr.rate,                    "emigr.rate.csv",                    row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(mortal.rate,                   "mortal.rate.csv",                   row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(improv.growth.rate,            "improv.growth.rate.csv",            row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(improv.emigr.rate,             "improv.emigr.rate.csv",             row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(improv.mortal.rate,            "improv.mortal.rate.csv",            row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  write.table(l,                             "repl.csv",                          row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  
  # outputs
  SumN <- apply(N,2,sum) # total population size
  write.table(t(SumN),             "SumN.csv",              row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
  
}

## Analysis on Individual behavior (1D experiment)
## The analysis is a post-hoc test arising from an observation of the data

## Setting
{
  ## packages
  library(data.table)
  library(lme4)
  library(car)
  library(ggplot2)
  library(ggridges)
  library("survminer")
  library(survival)
  
  ## function
  # Standard Error
  se  <-  function(x){
    y  <-  x[!is.na(x)]  #  remove  the  missing  values
    sqrt(var(as.vector(y))/length(y))
  }
  
  
  ## Data
  {
    Space <- "Location of Data in IndividualLocations.zip"
    Species <- c(rep("Hetero", 6), rep("Paraneo",6), rep("Reticuli",6))
    Colony <- rep(rep(c("A","B"), each=3), 3)
    Rep <- rep(c("1","2","3"),6)
    
    setwd(Space)
  }
}

## Data preparation
{
  DataScale <- function(Place, Scaling, Tunnel){
    Data <- data.frame(as.matrix(fread(Place, header=T)))
    na.omit(Data)
    Data <- Data[Data$position < 18001,]
    Data[,2:length(Data[1,])] <- Data[,2:length(Data[1,])]/Scaling
    Data[, seq(2,20,2)] <- Data[, seq(2,20,2)] - Tunnel/Scaling
    Data[, seq(3,21,2)] <- Data[, seq(3,21,2)] - max(Data[, seq(3,21,2)])
    return(Data)
  }
  Raw <- "Rawdata\\\\"
  Place <- paste(Space, Species, "\\\\", Colony, "-", Rep, "\\\\", sep="")
  #Scaling information to adjust pixel to mm
  #Hetero 191-1 (1208.2-283.8)/50 = 18.488
  #Hetero 191-3 (1233.1-294.2)/50 = 18.778
  #Hetero 191-4 (1194.6-280.6)/50 = 18.28
  #Hetero 198-2 (1192.2-254.9)/50 = 18.746, 305.5
  #Hetero 198-3 (1231.5-309.5)/50 = 18.44, 309
  #Hetero 198-4 (1223.5-302.3)/50 = 18.424, 302
  #Paraneo 203-2 (1158.4-263.7)/50 = 17.894, 350
  #Paraneo 203-3 (1218.6-374.5)/50 = 16.882, 374
  #Paraneo 203-5 (1164.8-340)/50 = 16.496, 340
  #Reticuli 206-2 (1229.1-293.4)/50 = 18.714
  #Reticuli 206-4 (1241.1-334.4)/50 = 18.134
  #Reticuli 206-5 (1237.1-332.0)/50 = 18.102
  Hscale <- c(18.488, 18.778, 18.28, 18.746, 18.44, 18.424)
  Pscale <- c(16.34, 16.25, 17.22, 17.894, 16.882, 16.496)
  Rscale <- c(19.89, 18.61, 18.61, 18.714, 18.134, 18.102)
  ScaleValue <- c(Hscale, Pscale, Rscale)
  
  #Position of the entrance of the tunnel
  Hlimit <- c(283.8, 294.2, 280.6, 305.5, 309, 302)
  Plimit <- c(336, 321, 336, 350, 374, 340)
  Rlimit <- c(278.5, 283, 283, 293.4, 334.4, 332.0)
  TunnelLimit <- c(Hlimit, Plimit, Rlimit)
  
  #Read data
  AllData <- NULL
  for(i in 1:length(ScaleValue)){
    f.namesplace <- list.files(Place[i], pattern=".csv",full.names=T)
    for(j in 1:length(f.namesplace)){
      d <- DataScale(f.namesplace[j], ScaleValue[i], TunnelLimit[i])
      sec = rep(d[,1]/30, 10)
      x = as.matrix(d[,seq(2,20,2)])[1:6010]
      y = as.matrix(d[,seq(3,21,2)])[1:6010]
      dis = c(NA, sqrt( diff(x)^2 + diff(y)^2 ))
      dis[sec==0] <- NA
      d <- data.frame(species = rep(Species[i], 6010),
                      colony = rep(Colony[i], 6010),
                      rep = rep(Rep[i], 6010),
                      hour = rep(j-1, 6010),
                      ind = rep(0:9, each=601),
                      sec=sec, x=x, y=y, dis=dis)
      AllData <- rbind(AllData, d)
    }
  }
  
  Hetero = AllData[AllData$species=="Hetero",]
  Reticuli = AllData[AllData$species=="Reticuli",]
  Paraneo = AllData[AllData$species=="Paraneo",]
}

## Analysis
{
  All <- rbind(Hetero,Paraneo,Reticuli)
  videoname <- paste(All[,1], All[,2], All[,3], All[,4], sep="-")
  videolist <- unique(videoname)
  dtunnel <- "data of tunnel length for 1D experiment"
  
  #### data.frame for analysis
  {
    d_species <- data.frame(species = c("Hetero", "Paraneo", "Reticuli"))
    d_group <- data.frame(
      species = c(rep("Hetero",6), rep("Paraneo",6), rep("Reticuli",6)),
      colony = rep( rep(c("A","B"), each=3), 3), rep = rep( c(1,2,3), 6))
    d_group <- cbind(d_group, scolony = paste(d_group[,1], d_group[,2], sep="-"),
                     srep = paste(d_group[,1], d_group[,2], d_group[,3], sep="-"))
    d_ind <- data.frame( 
      species = c(rep("Hetero",60), rep("Paraneo",60), rep("Reticuli",60)),
      colony = rep( rep(c("A","B"), each=30), 3),
      rep = rep( rep( rep(c(1,2,3),each=10), 2), 3), ind = rep( rep(0:9,6), 3))
    d_ind <- cbind(d_ind, scolony = paste(d_ind$species, d_ind$colony, sep="-"),
                   srep = paste(d_ind$species, d_ind$colony, d_ind$rep, sep="-"))
    d_video <- data.frame(
      video = videolist, dtunnel[,2:7]
    )
    d_video <- data.frame(d_video, group=paste(d_video$species, d_video$colony, d_video$rep, sep="-"))
    d_video_ind <- data.frame(
      video = rep(videolist, each=10), ind=0:9
    )
  }
  
  ### analysis for 1st individual
  # definition of front dig
  # get into dig zone (1mm from the front edge of the tunnel)
  # then back more than 2mm or 3mm from the dig zone
  Res = Res2 <- NULL
  ResCount<-NULL
  Texist <- NULL
  for(j in 1:length(videolist)){
    # data of each video
    df <- All[videoname==videolist[j],]
    begin <- dtunnel[dtunnel[,1]==videolist[j],6]
    end <- dtunnel[dtunnel[,1]==videolist[j],7]
    
    # order of individuals
    xorder <- tapply(df$x, df[,c(6,5)], sum)
    texist <- xorder>0
    xorder[1,][order(xorder[1,])] <- 10:1
    for(i in 2:length(xorder[,1])){
      if(sum(xorder[i,]>0)<1){
        xorder[i,] = xorder[i-1,]
      } else {
        xorder[i,][order(xorder[i,])] <- 10:1
      }
    }
    
    # individuals within tunnels
    Texist <- c(Texist, mean(apply(texist,1,sum)))
    
    # detect the change of the 1st individual
    count <- 0
    for(i in 2:length(xorder[,1])){
      if( (0:9)[xorder[i,]==1] != (0:9)[xorder[i-1,]==1]){
        count = count +1
      }
    }
    if(sum(tapply(df$x, df[,c(6,5)], sum)[1,]>0)<1){count <- count -1}
    ResCount <- c(ResCount, count)
    
    xorder <- xorder[1:6010]
    df <- cbind(df,xorder)
    
    # detect if 1st is at tunnel face
    if(begin>1.5){
      dig <- df$xorder==1 & df$x > begin-1.5
    } else {
      dig <- df$xorder==1 & df$x > 0
    }
    if(df[1,1]=="Paraneo"){ thresh = 3} else { thresh = 2}
    
    # detect if 2nd is interacting with 1st
    if(df[1,1]=="Paraneo"){ size = 6.5} else if(df[1,1]=="Hetero") { size = 4}  else {size=4.5}
    second <- (df$xorder==2) & (df$x > begin-1.5-size)
    
    # check the back behavior of 1st from the tunnel face
    lmax = lmin = rep(0, 6010)
    for(i in 1:6010){
      if(df[i,]$sec==0){ # initialization
        afterdig = 0; # detect if first dig or not
        digcount = 0; # detect if in dig or not
        trcount = 0;  # detect if in transport or not
        mem = 0;      # memorize the position of min or max
      }
      if(dig[i]){
        # just after geting in dig zone?
        if(afterdig==0 || trcount==1){ 
          if(trcount==1){
            lmin[mem] = 1;
          }
          afterdig=1; digcount=1; mem=i; trcount=0
        }
        if(df[mem,]$x < df[i,]$x){
          mem = i
        }
      } else {
        # exclude first phase
        if(afterdig==1){
          # just after getting into transpott?
          if(digcount==1){
            if(df[mem,]$x-df[i,]$x > thresh || df[i,]$x<0){
              lmax[mem] = 1;
              trcount=1; mem=i; digcount=0
            }
          }
          if(trcount == 1){
            if(df[mem,]$x > df[i,]$x){
              mem = i
            }
            if(df[i,]$sec==600){
              if(df[mem,]$x>0){
                lmin[mem] = 1;
              } else {
                lmin[mem] = 1;
              }
            }
          }
        }
      }
    }
    
    # summarize all results
    res1 <- cbind(df,lmax,lmin)[lmax>0,]
    res2 <- cbind(df,lmax,lmin)[lmin>0,]
    res1 <- res1[res2$lmin!=2,]; res2 <- res2[res2$lmin!=2,]
    res2$x[res2$x<0] = 0
    res <- data.frame(res1[,1:6],
                      backdis = res1$x-res2$x,
                      spenttime = res2$sec-res1$sec,
                      tunnellength = rep(end,length(res1[,1])),
                      tunnellength2 = rep(begin,length(res1[,1])))
    res
    Res <- rbind(Res,res)
    Res2 <- rbind(Res2, df[second,])
  }
  
  # summarize all resutls
  d_video <- data.frame(d_video, Texist)
  
  tgroup <- paste(Res$species, Res$colony, Res$rep, sep="-")
  Res <- cbind(Res, tgroup)
  d_back <- Res 
  
  bintunnellength <- floor(d_back$tunnellength/10)*10
  bintunnellength2 <- floor(d_back$tunnellength2/10)*10
  d_back <- data.frame(d_back, bintunnellength=bintunnellength, bintunnellength2=bintunnellength2)
  
  BackCount=as.vector(table(tgroup))
  d_group <- cbind(d_group, BackCount)
  d_species <- cbind(d_species, BackCountMean = tapply(d_group$BackCount, d_group$species, mean), 
                     BackCountMeanSE = tapply(d_group$BackCount, d_group$species, se))
  
  d_video <- data.frame(d_video, FrontChange = ResCount)
  FrontChange <- tapply(d_video$FrontChange, d_video[,4:2], sum)[1:18]
  d_group <- cbind(d_group, FrontChange)
  d_species <- data.frame(d_species, FrontChangeMean = tapply(d_group$FrontChange, d_group$species, mean),
                          FrontChangeSE = tapply(d_group$FrontChange, d_group$species, se))
  
  
  ### calcurate moved distance
  d_ind <- data.frame(d_ind,
                      # all distance
                      sumdis = c( tapply(Hetero$dis, Hetero[,c(5,3,2)], sum, na.rm=T)[1:60],
                                  tapply(Paraneo$dis, Paraneo[,c(5,3,2)], sum, na.rm=T)[1:60],
                                  tapply(Reticuli$dis, Reticuli[,c(5,3,2)], sum, na.rm=T)[1:60]),
                      
                      # moved distance inside tunnel
                      tunneldis = c( tapply(Hetero[Hetero$x>0,]$dis, Hetero[Hetero$x>0,c(5,3,2)], sum, na.rm=T)[1:60],
                                     tapply(Paraneo[Paraneo$x>0,]$dis, Paraneo[Paraneo$x>0,c(5,3,2)], sum, na.rm=T)[1:60],
                                     tapply(Reticuli[Reticuli$x>0,]$dis, Reticuli[Reticuli$x>0,c(5,3,2)], sum, na.rm=T)[1:60]
                      ),
                      
                      # x-axis movement inside tunnel
                      tunnelXdis = c( tapply((c(NA, abs(diff(Hetero$x)))*(Hetero$dis-Hetero$dis+1))[Hetero$x>0], Hetero[Hetero$x>0,c(5,3,2)], sum, na.rm=T)[1:60],
                                      tapply((c(NA, abs(diff(Paraneo$x)))*(Paraneo$dis-Paraneo$dis+1))[Paraneo$x>0], Paraneo[Paraneo$x>0,c(5,3,2)], sum, na.rm=T)[1:60],
                                      tapply((c(NA, abs(diff(Reticuli$x)))*(Reticuli$dis-Reticuli$dis+1))[Reticuli$x>0], Reticuli[Reticuli$x>0,c(5,3,2)], sum, na.rm=T)[1:60]
                      ),
                      
                      # time spent inside tunnel
                      tunneltime = c( tapply(Hetero$x>0, Hetero[,c(5,3,2)], sum, na.rm=T)[1:60],
                                      tapply(Paraneo$x>0, Paraneo[,c(5,3,2)], sum, na.rm=T)[1:60],
                                      tapply(Reticuli$x>0, Reticuli[,c(5,3,2)], sum, na.rm=T)[1:60]
                      )
  )
  d_ind$tunneldis[is.na(d_ind$tunneldis)] <- 0
  d_ind$tunnelXdis[is.na(d_ind$tunnelXdis)] <- 0
  
  d_group <- cbind(d_group,
                   sumdis = tapply(d_ind$sumdis, d_ind[,3:1], sum)[1:18],
                   tunneldis = tapply(d_ind$tunneldis, d_ind[,3:1], sum)[1:18],
                   tunnelXdis = tapply(d_ind$tunnelXdis, d_ind[,3:1], sum)[1:18])
  
  ## moved distance (total)
  d_species <- cbind(d_species, 
                     group_dis_mean = tapply(d_group$sumdis, d_group$species, mean),
                     group_dis_se = tapply(d_group$sumdis, d_group$species, se))
  
  ## moved distance within tunnel (x)
  d_species <- cbind(d_species, 
                     group_txdis_mean = tapply(d_group$tunnelXdis, d_group$species, mean),
                     group_txdis_se = tapply(d_group$tunnelXdis, d_group$species, se))
  
  ## how long each group takes to dig 50 mm
  digtime <- tapply(d_video$hour, d_video$group, max)
  d_group <- data.frame(d_group, digtime)
  
  
  # max-speed for each individual
  d_ind <- data.frame(d_ind, maxspeed = c( 
    tapply(Hetero$dis, Hetero[,c(5,3,2)], max, na.rm=T)[1:60],
    tapply(Paraneo$dis, Paraneo[,c(5,3,2)], max, na.rm=T)[1:60],
    tapply(Reticuli$dis, Reticuli[,c(5,3,2)], max, na.rm=T)[1:60]
  ))
  d_species <- cbind(d_species, 
                     maxspeed_mean = tapply(d_ind$maxspeed, d_ind$species, mean),
                     maxspeed_se = tapply(d_ind$maxspeed, d_ind$species, se))
  
  
  ### spending time at the front
  Res <- NULL
  for(j in 1:length(videolist)){
    # data of each video
    df <- All[videoname==videolist[j],]
    # order of individuals
    xorder <- tapply(df$x, df[,c(6,5)], sum)
    xorder[1,][order(xorder[1,])] <- 10:1
    for(i in 2:length(xorder[,1])){
      if(sum(xorder[i,]>0)<1){
        xorder[i,] = xorder[i-1,]
      } else {
        xorder[i,][order(xorder[i,])] <- 10:1
      }
    }
    dd <- (xorder==1)&(tapply(df$x, df[,c(6,5)], sum)>0)
    ind <- NULL
    fronttime <- NULL
    count <- 1
    for(i in 2:length(dd[,1])){
      if(sum(dd[i,])<1){
        if(sum(dd[i-1,])>0){
          ind <- c(ind, (0:9)[dd[i-1,]])
          fronttime <- c(fronttime, count)
          count <- 1
        }
      } else {
        if(sum(dd[i-1,])>0){
          if( (0:9)[dd[i,]] == (0:9)[dd[i-1,]] ){
            count <- count + 1
          } else {
            ind <- c(ind, (0:9)[dd[i-1,]])
            fronttime <- c(fronttime, count)
            count <- 1
          }
        }
      }
    }
    if(sum(dd[length(dd[,1]),])>0){
      ind <- c(ind, (0:9)[dd[length(dd[,1]),]])
      fronttime <- c(fronttime, count)
    }
    res <- data.frame(
      species = rep(df[1,1],length(ind)),
      colony = rep(df[1,2],length(ind)),
      rep =  rep(df[1,3],length(ind)),
      hour =  rep(df[1,4],length(ind)),
      ind, fronttime)
    Res <- rbind(Res, res)
  }
  ggplot(Res, aes(species, fronttime)) + geom_boxplot()
  
  cens <- 0
  for(i in 3:length(Res[,1])-1){
    if(Res[i,]$hour != Res[i-1,]$hour || Res[i,]$hour != Res[i+1,]$hour){
      cens <- c(cens,0)
    } else {cens <- c(cens,1)}
  }
  cens <- c(cens, 0)
  
  Res <- cbind(Res,cens)
  d_front <- Res
  
  # binning the tunnellength for some figures
  bin = floor( d_video$beginlength/10 ) * 10
  d_video <- data.frame(d_video, bin)
  
  # summarize results and output
  d_species <- transform(d_species, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  d_group <- transform(d_group, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  d_ind <- transform(d_ind, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  d_video <- transform(d_video, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  d_species <- transform(d_species, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  d_video_ind <- transform(d_video_ind, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  d_back <- transform(d_back, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  
  save(d_species, file="d_species.Rdata")
  save(d_group, file="d_group.Rdata")
  save(d_ind, file="d_ind.Rdata")
  save(d_video, file="d_video.Rdata")
  save(d_video_ind, file="d_video_ind.Rdata")
  save(d_back, file="d_back.Rdata")
  save(d_front, file="d_front.Rdata")
}

## Plots
{
  ## plot all trajectories
  {
    ggplot(data=Paraneo,aes(x=x,y=y,colour=as.factor(ind))) + 
      geom_path() + xlim(-25, 55) + ylim(-25,0) + ggtitle("Paraneo") + geom_vline(xintercept = 0, color = "black", size=0.5) + 
      facet_grid(hour~colony + rep) + coord_fixed() + theme_bw()
    
    p <- ggplot(data=Reticuli,aes(x=x,y=y,colour=as.factor(ind))) + geom_vline(xintercept = 0, color = "black", size=0.5) +
      geom_path() + xlim(-20, 55) + ylim(-20,0) + 
      facet_grid(hour~colony + rep) + coord_fixed() + theme_bw()
    
    p <- ggplot(data=Hetero,aes(x=x,y=y,colour=as.factor(ind))) + geom_vline(xintercept = 0, color = "black", size=0.5)
    p <- p + geom_path() + xlim(-20, 55) + ylim(-20,0)
    p + facet_grid(hour~colony + rep) + coord_fixed()+ theme_bw()
    
    ### plot all movement in x-axis
    p <- ggplot(data=Paraneo,aes(x=x,y=sec,colour=as.factor(ind)))
    p <- p + geom_path() + xlim(-25, 55) + ylim(600,0) + ggtitle("Paraneo")
    p <- p + facet_grid(hour~colony + rep)
    p + theme_bw() + geom_vline(xintercept = 0, color = "black", size=0.5)
    
    p <- ggplot(data=Reticuli,aes(x=x,y=sec,colour=as.factor(ind)))
    p <- p + geom_path() + xlim(-20, 55) + ylim(600,0) + ggtitle("Reticuli")
    p <- p + facet_grid(hour~colony + rep)
    p + theme_bw() + geom_vline(xintercept = 0, color = "black", size=0.5)
    
    p <- ggplot(data=Hetero,aes(x=x,y=sec,colour=as.factor(ind)))
    p <- p + geom_path() + xlim(-20, 55) + ylim(600,0) + ggtitle("Hetero")
    p <- p + facet_grid(hour~colony + rep)
    p + theme_bw() + geom_vline(xintercept = 0, color = "black", size=0.5)
  }
  
  ## plot separately video
  {
    for(j in 1:length(videolist)){
      # data of each video
      df <- All[videoname==videolist[j],]
      begin <- dtunnel[dtunnel[,1]==videolist[j],6]
      end <- dtunnel[dtunnel[,1]==videolist[j],7]
      dfm <- data.frame(ind = 0:9, x=tapply(df$x, df$ind, mean), y=tapply(df$y, df$ind, mean))
      
      pt <- ggplot(data=df,aes(x=x,y=y,colour=as.factor(ind))) + 
        geom_vline(xintercept = 0, color = "black", size=0.5) +
        geom_path() + xlim(-20, 60) + ylim(-20,0) + ggtitle(videolist[j]) + 
        coord_fixed() + xlab("x (mm)") + ylab("y (mm)") + guides(color=FALSE) + theme_bw() +
        theme(aspect=1/4, axis.title.x = element_blank())
      
      px <- ggplot(data=df,aes(x=x,y=sec,colour=as.factor(ind))) +
        geom_path() + xlim(-20, 60) + ylim(650,-50) + theme(legend.position = "bottom") +
        geom_vline(xintercept = 0, color = "black", size=0.5) + 
        xlab("x (mm)") + ylab("time (sec)") + theme_bw() + theme(aspect=1/4, axis.title.x = element_blank()) +
        guides(color=guide_legend(title="ind")) + 
        annotate("segment",x=begin,xend=begin,y=-50,yend=0,colour=1, size=0.1,arrow=arrow(length=unit(0.2, "cm"))) + 
        annotate("segment",x=end,xend=end,y=650,yend=600,colour=1, size=0.1,arrow=arrow(length=unit(0.2, "cm"))) +
        theme(legend.position = 'none')
      
      px2 <- ggplot(data=df, aes(x=ind, y=x, color=as.factor(ind))) + guides(color=FALSE) + 
        geom_point(alpha=0.1) + geom_path(alpha=0.4) + coord_flip() + 
        ylab("x (mm)") + xlab("individual") + theme_bw() + theme(aspect=1/4) +
        ylim(-20,60) + geom_hline(yintercept = 0, color = "black", size=0.5) +
        geom_point(data = dfm, size = 5, alpha=0.1, shape=21, color="black", aes(fill=as.factor(ind))) + 
        geom_point(data = dfm, size = 3, shape=48:57, col=1)+ theme(legend.position = 'none')
      
      (pt/px/px2) 
    }
  }
  
  ## plot separately group
  grouplist <- unique(paste(All$species, All$colony, All$rep, sep="-"))
  for( i in 1:18){
    showdata <- paste(All$species, All$colony, All$rep, sep="-") == grouplist[i]
    ggplot(data=All[showdata,],aes(x=x,y=y,colour=as.factor(ind))) + 
      geom_path(size=0.1) + xlim(-25, 55) + 
      ggtitle(grouplist[i]) + geom_vline(xintercept = 0, color = "black", size=0.5) + 
      scale_y_continuous(breaks=c(0,-10,-20), limits=c(-25,5)) +
      facet_wrap(~hour, nrow=7, ncol=3, labeller = label_both) + coord_fixed() + 
      theme_bw() + theme(strip.background = element_rect(colour="#00000000", fill="#00000000"), legend.position = 'none')
    ggsave(paste(grouplist[i],".pdf",sep=""), width=6)
  }
  for( i in 1:18){
    showdata <- paste(All$species, All$colony, All$rep, sep="-") == grouplist[i]
    ggplot(data=All[showdata,],aes(x=x,y=sec,colour=as.factor(ind))) +
      geom_path(size=0.1) + xlim(-20, 60) + ylim(620,-20) + 
      facet_wrap(~hour, nrow=7, ncol=3, labeller = label_both) +
      theme_bw() + theme(strip.background = element_rect(colour="#00000000", fill="#00000000"),
                         legend.position = 'none', aspect=3/8) +
      geom_vline(xintercept = 0, color = "black", size=0.5) + 
      xlab("x (mm)") + ylab("time (sec)")
    
    ggsave(paste("x_",grouplist[i],".pdf",sep=""), width=6)
  }
  
  ## back distance
  {
    ggplot(d_back[d_back$species=="Hetero",], aes(x=backdis, y=as.factor(bintunnellength2))) + 
      geom_density_ridges(fill=3, stat = "binline", binwidth=1, alpha=0.2, draw_baseline = T) +
      theme_ridges() + theme(aspect.ratio = 0.75, legend.position = 'none') +
      ylab("Tunnel length(mm)") + scale_y_discrete(breaks=c("0","10","20","30","40"), labels=c("0-10", "10-20", "20-30", "30-40", "40-50")) +
      ggtitle(bquote(bolditalic(.("H. aureus")))) + xlab("Back distance (mm)")
    ggsave("BackDisHetero.pdf", width = 3, height = 3)
    ggplot(d_back[d_back$species=="Paraneo",], aes(x=backdis, y=as.factor(bintunnellength2))) + 
      geom_density_ridges(fill=2, stat = "binline", binwidth=1, alpha=0.2, draw_baseline = T) +
      theme_ridges() + theme(aspect.ratio = 0.75, legend.position = 'none') +
      ylab("Tunnel length(mm)") + scale_y_discrete(breaks=c("0","10","20","30","40"), labels=c("0-10", "10-20", "20-30", "30-40", "40-50")) +
      ggtitle(bquote(bolditalic(.("P. simplicicornis")))) + xlab("Back distance (mm)")
    ggsave("BackDisParaneo.pdf", width = 3, height = 3)
    ggplot(d_back[d_back$species=="Reticuli",], aes(x=backdis, y=as.factor(bintunnellength2))) + 
      geom_density_ridges(fill=4, stat = "binline", binwidth=1, alpha=0.2, draw_baseline = T) +
      theme_ridges() + theme(aspect.ratio = 0.75, legend.position = 'none') +
      ylab("Tunnel length(mm)") + scale_y_discrete(breaks=c("0","10","20","30","40"), labels=c("0-10", "10-20", "20-30", "30-40", "40-50")) +
      ggtitle(bquote(bolditalic(.("R. tibialis")))) + xlab("Back distance (mm)")
    ggsave("BackDisReticuli.pdf", width = 3, height = 3)
  }
  
  ## moved distance within tunnel (x)
  ggplot(d_species, aes(species, group_txdis_mean*6)) + 
    geom_bar(stat="identity", color=1, fill=c(2,4,3), alpha=0.4) +
    geom_errorbar(aes(ymin=group_txdis_mean*6-group_txdis_se*6, ymax=group_txdis_mean*6+group_txdis_se*6), width=0.1) +
    theme_bw() + ylab("Group moved dis (mm)") + theme(aspect=1, legend.position = 'none') + 
    ggtitle("Moved distance in tunnel") +
    geom_jitter(data=d_group, aes(species, tunnelXdis*6, shape=colony), width=0.075)
  ggsave("TunnelMovedDisX.pdf", width=3, height = 3)
  
  ## dig time
  df<-survfit(Surv(digtime)~species, type = "kaplan-meier", data=d_group)
  ggsurvplot(fit = df, data = d_group,
             pval = F, pval.method = F,
             risk.table = F, conf.int = FALSE,
             ncensor.plot = FALSE, size = 1.5, linetype = 1:3,
             legend.title = "Species", legend.lab=c("Hetero","Paraneo","Reticuli"), 
             ylab="Digging probability", xlab="Time (hours)",
             ggtheme = theme_bw()  + theme(aspect.ratio = 0.75)) 
  ggsave("DigTime.pdf", width=4, height = 3)
  
  ### max moving speed
  ggplot(d_species, aes(species, maxspeed_mean, fill=species)) + geom_bar(stat="identity", color=1, alpha=0.4) +
    theme_bw() + ylab("Maximum speed (mm/sec)") + theme(aspect=1) + 
    ggtitle("Maximum moving speed") + ylim(0,10)+
    #geom_jitter(data=d_ind, aes(x=species, y=maxspeed, color=species), width=0.1, alpha=0.6) + 
    geom_errorbar(aes(ymin=maxspeed_mean-maxspeed_se, ymax=maxspeed_mean+maxspeed_se), width=0.1) 
  ggsave("MaximumSpeed.pdf", width=4, height = 4)
  
  ## Number of visiting tunnel face and Switching the 1st row
  {
    d_group=data.frame(d_group, tunnel_length=c(50,47.6,rep(50,16)))
    ggplot(d_group, aes(BackCount/tunnel_length*6, FrontChange/tunnel_length*6, col=species)) +
      geom_point() +
      theme_bw() + theme(aspect.ratio = 1) +xlim(0,30) + ylim(0,30)+
      xlab("The number of visiting the tunnel face (per 1 mm dig)")+
      ylab("The number of switching the 1st-row individual (per 1 mm dig)")
    ggsave("VisitTunnelFace.pdf", width=4, height = 4)
    
    r<-lm(BackCount/tunnel_length*6 ~ species, data=d_group)
    Anova(r)
    #Anova Table (Type II tests)
    #Response: BackCount/tunnel_length * 6
    #Sum Sq Df F value    Pr(>F)    
    #species   552.41  2  24.892 1.717e-05 ***
    #  Residuals 166.44 15    
    
    t.test(BackCount/tunnel_length*6 ~ colony, data=d_group[d_group$species=="Paraneo",])
    t.test(BackCount/tunnel_length*6 ~ colony, data=d_group[d_group$species=="Hetero",])
    t.test(BackCount/tunnel_length*6 ~ colony, data=d_group[d_group$species=="Reticuli",])
    
    
    r<-lm(FrontChange/tunnel_length*6 ~ species, data=d_group)
    Anova(r)
    #Analysis of Deviance Table (Type II Wald chisquare tests)
    #Response: FrontChange/tunnel_length * 6
    #Sum Sq Df F value    Pr(>F)    
    #species   687.04  2  24.307 1.968e-05 ***
    #  Residuals 211.99 15                      
    
    t.test(FrontChange/tunnel_length*6 ~ colony, data=d_group[d_group$species=="Paraneo",])
    t.test(FrontChange/tunnel_length*6 ~ colony, data=d_group[d_group$species=="Hetero",])
    t.test(FrontChange/tunnel_length*6 ~ colony, data=d_group[d_group$species=="Reticuli",])
  }
}

## Individual behavior (Fig. 3B)
{
  ## packages
  library(ggplot2)
  library(binom)
  
  ## data
  df <- data.frame(species=c("Hetero", "Paraneo", "Reticuli"), jaw=c(160,17,159), kick=c(0, 59, 0), 
                   compress=c(151,0,134), trans=c(160,76,159), nontrans=c(4,2,18))
  df <- transform(df, species= factor(species, levels = c("Paraneo", "Reticuli", "Hetero")))
  
  ## statistics
  x <- matrix(c(df$jaw, df$kick), 2, 3, byrow=T)
  fisher.test(x)
  # p-value < 2.2e-16
  
  x <- matrix(c(df$compress, df$trans-df$compress), 2, 3, byrow=T)
  fisher.test(x)
  # p-value < 2.2e-16
  
  ## plots
  CIs <- binom.confint(x=kick, n=jaw+kick, methods="wilson")
  df.kick <- data.frame(df,CIs)
  qplot(x=species, y=kick/(jaw+kick), ymin=lower, ymax=upper, data=df.kick, ylim=c(0,1), geom = "pointrange")  + 
    geom_bar(aes(species, kick/(jaw+kick)), stat="identity", col="black", alpha=0.4, fill=c(2,4,3)) +
    theme_bw() + theme(aspect=1)+ylab("Use of kicking behavior")
  #ggsave("200411_KickingBehavior.pdf", width = 3, height = 3)
  
  CIs <- binom.confint(x=df$compress, n=df$trans, methods="wilson")
  df.compress <- data.frame(df,CIs)
  qplot(x=species, y=compress/trans, ymin=lower, ymax=upper, data=df.compress, ylim=c(0,1), geom = "pointrange")  + 
    geom_bar(aes(species, compress/trans), stat="identity", col="black", alpha=0.4, fill=c(2,4,3)) +
    theme_bw() + theme(aspect=1)+ylab("Use of compressing behavior")
  #ggsave("200411_CompressingBehavior.pdf", width = 3, height = 3)
}

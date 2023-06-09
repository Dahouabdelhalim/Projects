################################################################################
#Fig S2 Hourly production and consumption profiles for all scenarios in the Online Appendix
#similar to Fig 8 but with 13 sample days instead of 3
################################################################################

#directory <- c("outputs_theta01","outputs_theta05","outputs_theta2")
directory <- c("outputs_theta01")
for (dirfold0 in directory) {
  dirfold <- paste("./outputshpc/",dirfold0,"/",sep="")
  print(getwd())
  thetaval <- as.numeric(substr(dirfold0, 14,16))
  thetaval=ifelse(thetaval!=2,thetaval/10,thetaval)
  print(paste("working with theta:",thetaval))
  
  #to produce all graphs with all combination, uncomment below
  # evsq <- c("half", "full","2016")
  # loadsq <- c("2007", "2045")
  # categorysq <- c("rps_100","fossil", "free")
  # scenariosq <- seq(1,3)
  # costsq <- c("current", "future")
  # pricingsq <- c("flat","dynamic")
  
  evsq <- c("half")
  loadsq <- c("2045")
  categorysq <- c("rps_100","fossil", "free")
  scenariosq <- c(2)
  costsq <- c("future")
  pricingsq <- c("flat","dynamic")
  
  scalem <- 1000000
  
  for (cat in categorysq) {
    for (ev in evsq) {
      for(flatdyn in pricingsq) {
        for (load in loadsq) {
          for (cost in costsq) {
            for (scen in scenariosq){
              tryCatch({ 
                if (load=="2007" & ev!="half"){
                  print("Skip")
                } else {
                  
                  print(paste(thetaval,load,ev,cat,cost,flatdyn,scen))
                  
                  #gather final iteration in housely energy sources file
                  filename_e <- paste(dirfold,'/energy_sources_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                  dfen <- read.csv(filename_e, header=TRUE, sep=",")
                  colnames(dfen)[colnames(dfen) == "timepoint_label"] <- "timepoint"
                  dfen <- transform(dfen, hour=as.numeric(format(
                    strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%H")))
                  dfen <- transform(dfen,month =as.numeric(format(
                    strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%m")))
                  dfen <- transform(dfen,day =as.numeric(format(
                    strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%d")))
                  dfen <- transform(dfen, net.final.mc=final.mc.energy+final.mc.energy.down-final.mc.energy.up)
                  dfen <- transform(dfen, net.final.price=offered.price.energy+offered.price.energy.down-offered.price.energy.up)
                  dfen <- transform(dfen, net.final.q=final.q.energy) 
                  dfen <- transform(dfen, batnet=Battery-StorageNetCharge)
                  dfen <- transform(dfen, DischargeBattery=ifelse(batnet>=0, batnet,0))
                  dfen <- transform(dfen, ChargeBattery=ifelse(batnet<0, batnet*-1,0))
                  dfen <- transform(dfen, evexp=dfen$ChargeEVs*dfen$net.final.price)
                  
                  #dfen <- dfen[(dfen$month==1&dfen$day==12)|(dfen$month==5&dfen$day==2)|(dfen$month==11&dfen$day==22),]
                  cat <- ifelse(cat=="rps","rps_100",cat) 
                  #making plots
                  ###############################################################################################
                  #make ggplot hourly production and consumption
                  #sum energy if 0 then exclude
                  
                  #1. production(LSFO	Biodiesel	Diesel	Coal	Pellet-Biomass	LNG	MSW	SUN	WND	
                  #curtail_MSW	curtail_SUN	curtail_WND)
                  
                  colnames(dfen)[colnames(dfen) == "ZoneTotalCentralDispatch"] <- "LZ_NetDispatch"
                  ymaxl <- roundUp(max(dfen$LZ_NetDispatch))
                  ifelse(ymaxl<4000, ymax <- 2500, ymax <- as.numeric(ymaxl))
                  #ymax <- 3000
                  ifelse(cat=="rps",renew <- "100% Clean", ifelse(cat=="free", renew <- "Unconstrained", renew <- cat))
                  ifelse(ev=="2016",evl <- "0.5%", ifelse(ev=="half", evl <- "50%", evl <- "100%"))
                  nam <- simpleCap(paste(renew,"-",flatdyn,"Pricing -", cost,"Cost Assumption -",evl,"EV -",load,"Load - Aggregate Elasticity",thetaval))
                  
                  ifelse(cost=="future"&load=="2045"&ev=="half",base <-"base", base <-"")
                  
                  pdf(paste("plot/FigS2_theta=",thetaval,"load",load,"ev",ev,cat,cost,flatdyn,"scen",scen,".pdf",sep=""), height = 3, width = 12)
                  
                  par(mfrow=(c(2,1)))
                  names(dfen)[names(dfen) == 'WND'] <- 'Wind'
                  names(dfen)[names(dfen) == 'SUN'] <- 'SolarPV'
                  names(dfen)[names(dfen) == 'MSW'] <- 'H.Power'
                  names(dfen)[names(dfen) == 'DispatchFuelCellMW'] <- 'FuelCell'
                  dfp <- dfen[,c("month","day","hour","curtail_WND", "curtail_SUN","curtail_MSW","FuelCell","DischargeBattery",
                                 "SolarPV","Wind","Pellet.Biomass","Biodiesel","Diesel",
                                 "LSFO","LNG", "H.Power","Coal")]
                  dfp_in <- dfp[,colSums(dfp) > 0]
                  dfp <- transform(dfp, monthhour=seq(1,nrow(dfp)))
                  dfpm <- melt(dfp, id=c("month","day","hour","monthhour"))
                  names(dfpm)[names(dfpm) == 'variable'] <- 'source'
                  
                  dfsel<- subset(dfpm, source==colnames(dfp_in)[3:ncol(dfp_in)][1])
                  for (i in 2:(ncol(dfp_in)-2)){
                    dfsel1<- subset(dfpm, source==colnames(dfp_in)[3:ncol(dfp_in)][i])
                    dfsel <- rbind(dfsel, dfsel1)
                  }
                  
                  
                  dfpm <- dfsel
                  dfpm <- transform(dfpm, value = round(value,2))
                  dfpm <- transform(dfpm, label=ifelse(hour==12,paste(day,"/",month,"/2045",sep=""),""))
                  p1 <- ggplot(dfpm[order(dfpm$source),], aes(monthhour, value, fill = source, order = (source))) +
                    geom_area(position = 'stack')+
                    geom_text(aes(label = label, y = -200), size = 3, colour="#330033")+
                    ylim(-200, ymax)+
                    scale_fill_manual(name="Source",
                                      values=c(LSFO="#CCCCCC", Biodiesel="#99CC00", Diesel="#999999", Coal="#000000",
                                               Pellet.Biomass="#003300", LNG="#666666", H.Power="#CC9900",
                                               Wind="#3399FF", SolarPV="#FFFF00", DischargeBattery="#990000",FuelCell="#DDA0DD",
                                               curtail_MSW="#FFCC66",curtail_SUN="#FFFF99",curtail_WND="#B0C4DE"))+
                    geom_vline(xintercept=24, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=48, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=72, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=96, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=120, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=144, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=168, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=192, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=216, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=240, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=264, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=288, col = "gray60", linetype = "longdash") +
                    labs( x="", y= "MW",title="Hourly Power Production (MW)")+
                    theme_classic()+
                    theme(axis.text.x=element_blank(), axis.ticks=element_blank(), legend.text = element_text(size=7),
                          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),legend.key.size = unit(0.4, 'cm'),
                          legend.position="bottom",legend.margin = margin(6, 6, 6, 6))+
                    guides(fill=guide_legend(nrow=1,byrow=TRUE))
                  p1
                  #3. price (base_price offered price energy final price energy)
                  dfp <- dfen[,c("month","day","hour","base_price",
                                 #"offered.price.energy", 
                                 "final.price.energy")]
                  dfp <- transform(dfp, monthhour=seq(1,nrow(dfp)))
                  dfpm <- melt(dfp, id=c("month","day","hour","monthhour"))
                  dfpm <- transform(dfpm, labelbase=ifelse(hour==5 & month==2 & variable=="base_price","Base price",""))
                  #dfpm <- transform(dfpm, labeloffp=ifelse(hour==12 & month==11 & variable=="offered.price.energy","Offered price",""))
                  dfpm <- transform(dfpm, labelofff=ifelse(hour==1 & month==11 & variable=="final.price.energy","Final price",""))
                  p2 <- ggplot(dfpm, aes(monthhour, value, color=variable)) + 
                    geom_line(linetype = 'twodash', size=0.65) + 
                    xlab("") +
                    ylab("") + 
                    ylim(-20, 600) +
                    geom_text(aes(y=value+15,label=labelbase), size=3, colour="#330033")+
                    #geom_text(aes(y=value+20, label=labeloffp), size=3, colour="#FF3300")+
                    geom_text(aes(y=value+30, label=labelofff), size=3, colour="#FF9999")+
                    scale_color_manual(name="Price",
                                       values=c(base_price="#330033",
                                                #offered.price.energy="#FF3300",
                                                final.price.energy="#FF9999"))+
                    theme(axis.text.x=element_blank(), legend.position="bottom", 
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(fill = NA))
                  grid.newpage()
                  # extract gtable
                  g1 <- ggplot_gtable(ggplot_build(p1))
                  g2 <- ggplot_gtable(ggplot_build(p2))
                  
                  # overlap the panel of 2nd plot on that of 1st plot
                  pp <- c(subset(g1$layout, name == "panel", se = t:r))
                  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                                       pp$l, pp$b, pp$l)
                  
                  # axis tweaks
                  ia <- which(g2$layout$name == "axis-l")
                  ga <- g2$grobs[[ia]]
                  ax <- ga$children[[2]]
                  ax$widths <- rev(ax$widths)
                  ax$grobs <- rev(ax$grobs)
                  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(1, "cm")
                  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
                  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
                  
                  # draw it
                  grid.draw(g)
                  #dev.off()
                  #2. consumption (DispatchFuelCellMW	DischargeBattery	DRUnservedLoad	FlexibleDemand	
                  #DumpPower	ChargeEVs	RunElectrolyzerMW	LiquifyHydrogenMW	ChargeBattery)
                  #subset for consumption 
                  #pdf(paste("HourlyCons",cat,cost,flatdyn,load,ev,"scen",scen,"iter",maxi,"theta=",thetaval,".pdf",sep=""), height = 4, width = 12)
                  names(dfen)[names(dfen) == 'LiquifyHydrogenMW'] <- 'LiquifyHidrogren'
                  names(dfen)[names(dfen) == 'RunElectrolyzerMW'] <- 'ProduceHydrogen'
                  dfp <- dfen[,c("month","day","hour", "ProduceHydrogen","ChargeEVs", 
                                 "LiquifyHidrogren", "ChargeBattery", "FlexibleDemand")]
                  
                  dfp_in <- dfp[,colSums(dfp) > 0]
                  dfp <- transform(dfp, monthhour=seq(1,nrow(dfp)))
                  dfpm <- melt(dfp, id=c("month","day","hour","monthhour"))
                  names(dfpm)[names(dfpm) == 'variable'] <- 'source'
                  
                  #subset for base load and bid q
                  dfpbid <- dfen[,c("month","hour","base_load"
                                    #,"bid.q.energy"
                  )]
                  dfpbid <- transform(dfpbid, monthhour=seq(1,nrow(dfpbid)))
                  dfpmb <- melt(dfpbid, id=c("month","hour","monthhour"))
                  dfpmb <- transform(dfpmb, labelbase=ifelse(hour==1 & month==12 & variable=="base_load","Base load",""))
                  dfpmb <- transform(dfpmb, labeloffp=ifelse(hour==13 & month==12 & variable=="bid.q.energy","Last bid",""))
                  
                  #delete category = 0
                  dfsel<- subset(dfpm, source==colnames(dfp_in)[3:ncol(dfp_in)][1])
                  for (i in 2:(ncol(dfp_in)-2)){
                    dfsel1<- subset(dfpm, source==colnames(dfp_in)[3:ncol(dfp_in)][i])
                    dfsel <- rbind(dfsel, dfsel1)
                  }
                  
                  dfpm <- dfsel
                  dfpm <- transform(dfpm, label=ifelse(hour==12,paste(day,"/",month,"/2045",sep=""),""))
                  p3 <- ggplot(dfpm[order(dfpm$source),], aes(monthhour, value, fill = source)) +
                    geom_area(position = 'stack')+
                    xlab("") +
                    ylab("MW") + 
                    geom_text(aes(label = label, y = -200), size = 3, colour="#330033")+
                    ylim(-200, ymax)+
                    geom_vline(xintercept=24, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=48, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=72, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=96, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=120, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=144, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=168, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=192, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=216, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=240, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=264, col = "gray60", linetype = "longdash") +
                    geom_vline(xintercept=288, col = "gray60", linetype = "longdash") +
                    scale_fill_manual(name="Usage",
                                      values=c(FlexibleDemand="#CCCCCC",
                                               ChargeBattery="#990000",
                                               LiquifyHidrogren="#FF9933",
                                               ChargeEVs="#666600",
                                               ProduceHydrogen="#0033FF"))+
                    ggtitle("Hourly Power Consumption (MW)")+ theme_classic()+
                    theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
                          plot.title = element_text(hjust = 0.5), legend.position="bottom",legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=8))+
                    guides(fill=guide_legend(nrow=1,byrow=TRUE))
                  
                  p4 <- ggplot(dfpmb, aes(monthhour, value, color=variable))+
                    geom_line(linetype = 'dotted', size=0.65) + 
                    labs(x="", y="") +
                    ylim(-200, ymax)+
                    geom_text(aes(y=value-100, label=labelbase), size=3, colour="#330033")+
                    geom_text(aes(y=value-100, label=labeloffp), size=3, colour="#FF3300")+
                    scale_color_manual(name="Usage",
                                       values=c(base_load="#330033",
                                                bid.q.energy="#FF3300"))+
                    theme(axis.text.x=element_blank(), legend.position="bottom",
                          axis.text.y=element_text(colour="white"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(fill = NA))
                  grid.newpage()
                  
                  # extract gtable
                  g1 <- ggplot_gtable(ggplot_build(p3))
                  g2 <- ggplot_gtable(ggplot_build(p4))
                  
                  # overlap the panel of 2nd plot on that of 1st plot
                  pp <- c(subset(g1$layout, name == "panel", se = t:r))
                  g4 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                                        pp$l, pp$b, pp$l)
                  
                  # axis tweaks
                  ia <- which(g2$layout$name == "axis-l")
                  ga <- g2$grobs[[ia]]
                  ax <- ga$children[[2]]
                  ax$widths <- rev(ax$widths)
                  ax$grobs <- rev(ax$grobs)
                  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(1, "cm")
                  g1 <- gtable_add_cols(g4, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
                  g1 <- gtable_add_grob(g4, ax, pp$t, length(g$widths) - 1, pp$b)
                  
                  # draw it
                  grid.draw(g1)
                  dev.off()
                  
                }}, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
            }
          }
        }
      }
    }
  }
}
dev.off()

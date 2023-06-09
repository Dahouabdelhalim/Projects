#This code produces:
#Table 1: Assumptions about flexible demand and demand-side reserves
#Figure 1: Demand flexibility scenarios by hour and month
#Table 2: Generators Cost Assumption
#Table 3: Fuel Cost Assumption
#Figure 3: Average output and potential capacity of renewable energy sources on Oahu

#######################################################################
#Figure 3: Average output and potential capacity of renewable energy sources on Oahu
#######################################################################
flexshares <- read.csv(file="inputs/flexshares.csv", header=TRUE, sep=",")
flexshares$hour <- seq(1,24)
flexshares <- melt(flexshares, id=c("hour"))
flexshares$variable <- gsub("X", "", flexshares$variable)
names(flexshares)[names(flexshares) == 'variable'] <- 'month'
names(flexshares)[names(flexshares) == 'value'] <- 'flex'
flexshares$inflex <- 1-flexshares$flex
flexshares$hour <- gsub(24,0, flexshares$hour)

setwd("./inputs")
#capacity factor graph 
#x-axis MW available
#y-axis project capacity factor
filename <- 'generation_projects_info.csv'
capacitydf <- read.csv(filename, header=TRUE, sep=",")
unique(capacitydf$gen_tech)
capacitydf <- capacitydf[which(capacitydf$gen_tech=="CentralTrackingPV"|
                                 capacitydf$gen_tech=="FlatDistPV"|
                                 capacitydf$gen_tech=="OffshoreWind"|
                                 capacitydf$gen_tech=="OnshoreWind"|
                                 capacitydf$gen_tech=="SlopedDistPV"),c("gen_tech","GENERATION_PROJECT","gen_capacity_limit_mw")] 
filename <- 'variable_capacity_factors.csv'
capacityfdf <- read.csv(filename, header=TRUE, sep=",")
setwd('..')

idlist <- list(capacityfdf$GENERATION_PROJECT)
capacityfdf <- aggregate(capacityfdf[,c("gen_max_capacity_factor")], idlist, mean)
names(capacityfdf)[names(capacityfdf) == 'Group.1'] <- 'GENERATION_PROJECT'
names(capacityfdf)[names(capacityfdf) == 'x'] <- 'cf'
all<- merge(capacitydf, capacityfdf, by=c("GENERATION_PROJECT"))
idlist <- list(all$gen_tech)
allmin <- aggregate(all[,c("cf")], idlist, max)
allmin$gen_capacity_limit_mw<-0
names(allmin)[names(allmin) == 'Group.1'] <- 'gen_tech'
names(allmin)[names(allmin) == 'x'] <- 'cf'
all <- rbind.fill(allmin,all)

#plot
scale <- 1
myColors <- brewer.pal(5,"Set1")
all <- all[order(all[,c("gen_tech")],all[,c("cf")],decreasing = TRUE),]
all$Capacity<-ave(all$gen_capacity_limit_mw,all$gen_tech,FUN=cumsum)
all$gen_tech<-as.factor(all$gen_tech)
names(myColors) <- levels(all$gen_tech)
colScale <- scale_colour_manual(name = "Technology",values = myColors)
all$Technology <- all$gen_tech

p1 <- ggplot(all, aes(Capacity, cf,color=Technology, shape=Technology)) +
  geom_point()+geom_line()+
  theme_classic()+ylab("Project Capacity Factor")+xlab("MW Available")+
  theme(legend.position = c(0.7, 0.80),
        legend.title = element_text(colour = "black", size = 10, face = "bold"),
        legend.text = element_text(colour = "black", size = 10), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))
p1 + colScale
p1 <- p1 + scale_x_continuous(expand = c(0, 0))
print(p1)
ggsave("plot/Fig3_RenewablePotentialCapacity.pdf",height=5.1,width=10)
dev.off()

########################################################################
#Table 2: Generators Cost Assumption
########################################################################
setwd("./inputs")
#table cost assumptions 
#gen_storage_energy_overnight_cost=battery_cost_per_mwh_cycled

filename <- 'gen_build_costs.csv'
costdf <- read.csv(filename, header=TRUE, sep=",")
costdf$gen_storage_energy_overnight_cost <- sub(".", "0", costdf$gen_storage_energy_overnight_cost)
costdf$gen_storage_energy_overnight_cost <- as.numeric(costdf$gen_storage_energy_overnight_cost)
idlist <- list(costdf$GENERATION_PROJECT,costdf$build_year)
costdf <- aggregate(costdf[,c("gen_overnight_cost","gen_storage_energy_overnight_cost","gen_fixed_om")], idlist, mean)
names(costdf)[names(costdf) == 'Group.1'] <- 'GENERATION_PROJECT'
names(costdf)[names(costdf) == 'Group.2'] <- 'Year'

filename <- 'generation_projects_info.csv'
capacitydf <- read.csv(filename, header=TRUE, sep=",")
capacitydf <- capacitydf[,c("GENERATION_PROJECT","gen_tech",
                            "gen_scheduled_outage_rate" ,"gen_max_age",
                            "gen_variable_om")]
costdf<- merge(costdf, capacitydf, by=c("GENERATION_PROJECT"))
costdf <- transform(costdf, gen_tech=ifelse(gen_tech=="Waiau_9","Waiau",gen_tech))
costdf <- transform(costdf, gen_tech=ifelse(gen_tech=="Waiau_10","Waiau",gen_tech))
idlist <- list(costdf$gen_tech,costdf$Year)
costdf <- aggregate(costdf[,c("gen_overnight_cost",
                              "gen_fixed_om","gen_variable_om","gen_max_age",
                              "gen_scheduled_outage_rate"
                              )], idlist, mean)
scalem <- 1000
costdf$gen_overnight_cost <- round(costdf$gen_overnight_cost/scalem,2)
costdf$gen_fixed_om <- round(costdf$gen_fixed_om/scalem,2)
costdf$gen_scheduled_outage_rate <- round(costdf$gen_scheduled_outage_rate,2)
costdf$gen_variable_om <- round(costdf$gen_variable_om,2)
names(costdf)[names(costdf) == 'Group.1'] <- 'gen_tech'
names(costdf)[names(costdf) == 'Group.2'] <- 'year_build'
dat <- costdf
dat[dat==0] <- ""
dat$gen_max_age <- as.character(dat$gen_max_age)
dat <- dat[order(dat$gen_tech,dat$year_build),]
dat$tech <- ""

lis <- c("Airport_DSG","Waiau","CC_152","CIP_CT","IC_Barge","IC_MCBH",
         "IC_Schofield", "Kalaeloa_CC1", "Kalaeloa_CC2","Kalaeloa_CC3")
for (z in lis) {
nam <- "Fossil-fueled"
dat <- transform(dat, tech=ifelse(gen_tech==z,nam,tech))
}  
lis <- c("Battery_Bulk","Battery_Conting","Battery_Reg","DistBattery")
for (z in lis) {
  nam <- "Storage"
  dat <- transform(dat, tech=ifelse(gen_tech==z,nam,tech))
}  
lis <- c("OffshoreWind","OnshoreWind")
for (z in lis) {
  nam <- "Wind"
  dat <- transform(dat, tech=ifelse(gen_tech==z,nam,tech))
} 
lis <- c("FlatDistPV","SlopedDistPV","CentralTrackingPV" )
for (z in lis) {
  nam <- "Solar PV"
  dat <- transform(dat, tech=ifelse(gen_tech==z,nam,tech))
} 
dat <- transform(dat, tech=ifelse(gen_tech=="H-Power","Biomass",tech))
dat <- dat[,c(8,1,2,3,4,5,6,7)]
dat <- dat[order(dat$tech),]

filename <- 'hydrogen.csv'
hyddf <- read.csv(filename, header=TRUE, sep=",")
df <- t(hyddf)
hyd1 <- data.frame("Storage","HydrogenElectrolyzer",2045,df[1,1],df[2,1],df[3,1],df[4,1],"")
colnames(hyd1) <- c("Technology","Generators","Year Build","Capital Cost ($/KW)"
                   ,"Fixed Cost (Mio$/Year)","Variable Cost ($/MWh)","Age (Years)","Outage Rate")
hyd2 <- data.frame("Storage","HydrogenFuelCell",2045,df[6,1],df[7,1],df[9,1],df[8,1],"")
colnames(hyd2) <- c("Technology","Generators","Year Build","Capital Cost ($/KW)"
                    ,"Fixed Cost (Mio$/Year)","Variable Cost ($/MWh)","Age (Years)","Outage Rate")
hyd3 <- data.frame("Storage","HydrogenLiquifier",2045,df[11,1],df[12,1],df[14,1],df[13,1],"")
colnames(hyd3) <- c("Technology","Generators","Year Build","Capital Cost ($/KW)"
                    ,"Fixed Cost (Mio$/Year)","Variable Cost ($/MWh)","Age (Years)","Outage Rate")
df <- rbind(hyd1,hyd2)
df <- rbind(df,hyd3)
df$'Capital Cost ($/KW)' <- round(df$'Capital Cost ($/KW)'/scalem,2)
df$'Fixed Cost (Mio$/Year)' <- round(df$'Fixed Cost (Mio$/Year)'/scalem,2)
df$'Variable Cost ($/MWh)' <- round(df$'Variable Cost ($/MWh)',2)

colnames(dat) <- c("Technology", "Generators","Year Build","Capital Cost ($/KW)"
                   ,"Fixed Cost (Mio$/Year)","Variable Cost ($/MWh)","Age (Years)"
                   ,"Outage Rate")
dat <- rbind(df,dat)
dat <- dat[order(dat$Technology),]
latextab <- xtable(dat,caption=paste("Generators Cost Assumption"),sep="")
print(getwd())
setwd('..')
print(latextab,
      type="latex", 
      include.rownames=FALSE, 
      file=paste("tables/Tab2_GeneratorCosts.tex"))

########################################################################
#Table 3: Fuel Cost Assumption
########################################################################
setwd("./inputs")
filename <- 'fuel_supply_curves.csv'
fueldf <- read.csv(filename, header=TRUE, sep=",")
idlist <- list(fueldf$fuel,fueldf$period)
fueldf <- aggregate(fueldf[,c("unit_cost")], idlist, mean)
colnames(fueldf) <- c("Fuel","Year","Unit cost ($/MMBtu)")
latextab <- xtable(fueldf,caption=paste("Fuel Cost Assumption"),sep="")
setwd('..')
print(latextab,
      type="latex", 
      include.rownames=FALSE, 
      file=paste("tables/Tab3_FuelCost.tex"))

########################################################################
#Table 1: Assumptions about flexible demand and demand-side reserves
#Figure 1: Demand flexibility scenarios by hour and month
########################################################################
#loads graph
setwd("./inputs")
filename <- 'timepoints.csv'
timelabeldf <- read.csv(filename, header=TRUE, sep=",")
filename <- 'loads.csv'
loadsdf <- read.csv(filename, header=TRUE, sep=",")
names(loadsdf)[names(loadsdf) == 'TIMEPOINT'] <- 'timepoint_id'
timelabeldf <- merge(timelabeldf, loadsdf, by=c("timepoint_id"))

filename <- 'timeseries.csv'
sampleweightdf <- read.csv(filename, header=TRUE, sep=",")
sampleweightdf$samplew <- round(sampleweightdf$ts_scale_to_period/sum(sampleweightdf$ts_scale_to_period),2)

names(sampleweightdf)[names(sampleweightdf) == 'TIMESERIES'] <- 'timeseries'
timelabeldf <- merge(timelabeldf, sampleweightdf, by=c("timeseries"))
timelabeldf$hour <- hour(timelabeldf$timestamp)
timelabeldf$month <- month(timelabeldf$timestamp)
timelabeldf<- merge(timelabeldf, flexshares, by=c("month","hour"))

setwd('..')

set.scenario <- function(scen) {
  ifelse(scen==1, {#optimistic
    shareflex <<- c(0.67, 0.05, 0.28) #sigma level - highly flex, midflex, inflex
    shareother <<- c(0.15, 0.05, 0.8)},
    ifelse(scen==2, 
           {#2nd scenario=baseline/moderate
             shareflex <<- c(0.33, 0.05, 0.62)
             shareother <<- c(0.075, 0.05, 0.875)},
           ifelse(scen==3,
                  {#3rd scenario=pessimistic
                    shareflex <<- c(0.15, 0.05, 0.8)
                    shareother <<- c(0, 0.05, 0.95)},
                  print("Scenario is not available"))))
}
################################################################################
#creating the Table 1
scen=1
set.scenario(scen)
timelabeldf <- timelabeldf[order(timelabeldf$timestamp,timelabeldf$month,timelabeldf$hour),]
x <- timelabeldf[,c("timestamp","zone_demand_mw","ts_scale_to_period","samplew","flex","inflex")]

x$id<-seq(1,nrow(x))
x$HighFlexible <- ((shareflex[1]*x$flex)+(shareother[1]*x$inflex))*x$zone_demand_mw
x$MidFlexible <- ((shareflex[2]*x$flex)+(shareother[2]*x$inflex))*x$zone_demand_mw
x$InFlexible <- ((shareflex[3]*x$flex)+(shareother[3]*x$inflex))*x$zone_demand_mw
x <- x[,c("timestamp","zone_demand_mw","id","ts_scale_to_period","samplew","HighFlexible","MidFlexible","InFlexible")]

df1 <- x %>%
  summarize(across(c(HighFlexible,MidFlexible,InFlexible), 
                   list(w.mean = ~weighted.mean(., w = ts_scale_to_period), 
                        sd = ~sd(., na.rm = TRUE)), .groups = 'drop'))
df1$scen <- 1

scen=2
set.scenario(scen)
x <- timelabeldf[,c("timestamp","zone_demand_mw","ts_scale_to_period","samplew","flex","inflex")]

x$id<-seq(1,nrow(x))
x$HighFlexible <- ((shareflex[1]*x$flex)+(shareother[1]*x$inflex))*x$zone_demand_mw
x$MidFlexible <- ((shareflex[2]*x$flex)+(shareother[2]*x$inflex))*x$zone_demand_mw
x$InFlexible <- ((shareflex[3]*x$flex)+(shareother[3]*x$inflex))*x$zone_demand_mw
x <- x[,c("timestamp","zone_demand_mw","id","ts_scale_to_period","samplew","HighFlexible","MidFlexible","InFlexible")]

df2 <- x %>%
  summarize(across(c(HighFlexible,MidFlexible,InFlexible), 
                   list(w.mean = ~weighted.mean(., w = ts_scale_to_period), 
                        sd = ~sd(., na.rm = TRUE)), .groups = 'drop'))
df2$scen <- 2

scen=3
set.scenario(scen)
x <- timelabeldf[,c("timestamp","zone_demand_mw","ts_scale_to_period","samplew","flex","inflex")]

x$id<-seq(1,nrow(x))
x$HighFlexible <- ((shareflex[1]*x$flex)+(shareother[1]*x$inflex))*x$zone_demand_mw
x$MidFlexible <- ((shareflex[2]*x$flex)+(shareother[2]*x$inflex))*x$zone_demand_mw
x$InFlexible <- ((shareflex[3]*x$flex)+(shareother[3]*x$inflex))*x$zone_demand_mw
x <- x[,c("timestamp","zone_demand_mw","id","ts_scale_to_period","samplew","HighFlexible","MidFlexible","InFlexible")]

df3 <- x %>%
  summarize(across(c(HighFlexible,MidFlexible,InFlexible), 
                   list(w.mean = ~weighted.mean(., w = ts_scale_to_period), 
                        sd = ~sd(., na.rm = TRUE)), .groups = 'drop'))
df3$scen <- 3
df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- df %>% mutate(across(1:6, round, 0))
colnames(df) <- c("high(mean)","high(sd)","mid(mean)","mid(sd)","inflex(mean)","inflex(sd)")
write.csv(df, file = paste("tables/Tab1_shareflexload.csv", sep = ""), row.names = FALSE)

################################################################################
#Creating Figure 1
################################################################################
#for each scenario
scen=1
set.scenario(scen)
timelabeldf <- timelabeldf[order(timelabeldf$timestamp,timelabeldf$month,timelabeldf$hour),]
x <- timelabeldf[,c("timestamp","zone_demand_mw","samplew","flex","inflex")]

x$id<-seq(1,nrow(x))
x$HighFlexible <- ((shareflex[1]*x$flex)+(shareother[1]*x$inflex))*x$zone_demand_mw
x$MidFlexible <- ((shareflex[2]*x$flex)+(shareother[2]*x$inflex))*x$zone_demand_mw
x$InFlexible <- ((shareflex[3]*x$flex)+(shareother[3]*x$inflex))*x$zone_demand_mw
x <- x[,c("timestamp","zone_demand_mw","id","samplew","HighFlexible","MidFlexible","InFlexible")]
all <- melt(x, id=c("timestamp","zone_demand_mw","id","samplew"))

all$hour <- hour(all$timestamp)
all$day <- day(all$timestamp)
all$month <- month(all$timestamp)
all$year <- year(all$timestamp)
all <- transform(all, labeldate=ifelse(hour==12,paste(day,"/",month,"/",year,sep=""),""))
all <- transform(all, labelw=ifelse(hour==12&id==13,paste("Weight=",samplew,sep=""),
                                    ifelse(hour==12,paste(samplew,sep=""),"")))
all$Flexibility <- all$variable

#pdf(paste("plot/Fig1_baseload_","scen",scen,".pdf",sep=""), height = 3, width = 12)
p1 <- ggplot(all, aes(id, value, fill = Flexibility)) + 
  geom_area(position = 'stack')+xlab("Sample Days") +
  ylab("MW") + theme_classic()+
  geom_text(aes(label = labeldate, y = -200), size = 3, colour="#330033")+
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
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position="bottom")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
print(p1)
ggsave(paste("plot/Fig1_baseload_","scen",scen,".pdf",sep=""),height=3,width=12)
dev.off()
############################################################################################################
scen=2
set.scenario(scen)
timelabeldf <- timelabeldf[order(timelabeldf$timestamp,timelabeldf$month,timelabeldf$hour),]
x <- timelabeldf[,c("timestamp","zone_demand_mw","samplew","flex","inflex")]

x$id<-seq(1,nrow(x))
x$HighFlexible <- ((shareflex[1]*x$flex)+(shareother[1]*x$inflex))*x$zone_demand_mw
x$MidFlexible <- ((shareflex[2]*x$flex)+(shareother[2]*x$inflex))*x$zone_demand_mw
x$InFlexible <- ((shareflex[3]*x$flex)+(shareother[3]*x$inflex))*x$zone_demand_mw
x <- x[,c("timestamp","zone_demand_mw","id","samplew","HighFlexible","MidFlexible","InFlexible")]
all <- melt(x, id=c("timestamp","zone_demand_mw","id","samplew"))

all$hour <- hour(all$timestamp)
all$day <- day(all$timestamp)
all$month <- month(all$timestamp)
all$year <- year(all$timestamp)
all <- transform(all, labeldate=ifelse(hour==12,paste(day,"/",month,"/",year,sep=""),""))
all <- transform(all, labelw=ifelse(hour==12&id==13,paste("Weight=",samplew,sep=""),
                                    ifelse(hour==12,paste(samplew,sep=""),"")))
all$Flexibility <- all$variable

#pdf(paste("plot/Fig1_baseload_","scen",scen,".pdf",sep=""), height = 3, width = 12)
p2 <- ggplot(all, aes(id, value, fill = Flexibility)) + 
  geom_area(position = 'stack')+xlab("Sample Days") +
  ylab("MW") + theme_classic()+
  geom_text(aes(label = labeldate, y = -200), size = 3, colour="#330033")+
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
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position="bottom")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
p2
ggsave(paste("plot/Fig1_baseload_","scen",scen,".pdf",sep=""),height=3,width=12)
#dev.off()

############################################################################################################
scen=3
set.scenario(scen)
timelabeldf <- timelabeldf[order(timelabeldf$timestamp,timelabeldf$month,timelabeldf$hour),]
x <- timelabeldf[,c("timestamp","zone_demand_mw","samplew","flex","inflex")]

x$id<-seq(1,nrow(x))
x$HighFlexible <- ((shareflex[1]*x$flex)+(shareother[1]*x$inflex))*x$zone_demand_mw
x$MidFlexible <- ((shareflex[2]*x$flex)+(shareother[2]*x$inflex))*x$zone_demand_mw
x$InFlexible <- ((shareflex[3]*x$flex)+(shareother[3]*x$inflex))*x$zone_demand_mw
x <- x[,c("timestamp","zone_demand_mw","id","samplew","HighFlexible","MidFlexible","InFlexible")]
all <- melt(x, id=c("timestamp","zone_demand_mw","id","samplew"))

  all$hour <- hour(all$timestamp)
  all$day <- day(all$timestamp)
  all$month <- month(all$timestamp)
  all$year <- year(all$timestamp)
  all <- transform(all, labeldate=ifelse(hour==12,paste(day,"/",month,"/",year,sep=""),""))
  all <- transform(all, labelw=ifelse(hour==12&id==13,paste("Weight=",samplew,sep=""),
                                      ifelse(hour==12,paste(samplew,sep=""),"")))
  all$Flexibility <- all$variable
  
  p3 <- ggplot(all, aes(id, value, fill = Flexibility)) + 
    geom_area(position = 'stack')+xlab("Sample Days") +
    ylab("MW") + theme_classic()+
    geom_text(aes(label = labeldate, y = -200), size = 3, colour="#330033")+
    geom_text(aes(label = labelw, y = 100), size = 3, colour="#330033")+
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
    theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
          plot.title = element_text(hjust = 0.5), legend.position="bottom")+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  print(p3)
  ggsave(paste("plot/Fig1_baseload_","scen",scen,".pdf",sep=""),height=3,width=12)
  
dev.off()

## Calculating the energy gain of a humpback whale foraging on different prey types##
setwd("C:/Users/emchenoweth/Desktop/New folder/Step 1 Propoprtion Captured")

##packages##
library("ggplot2", lib.loc="~/R/win-library/3.5")
library("plyr", lib.loc="~/R/win-library/3.5")
library("doBy", lib.loc="~/R/win-library/3.5")
library("directlabels", lib.loc="~/R/win-library/3.5")

##basic functions##
assignVec <- Vectorize("assign",c("x","value")) # from StackOverflow user  Joran: https://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line-in-r/13353585
deg2rad <- function(deg) {(deg * pi) / (180)} # from StackOverflow user Iris: https://stackoverflow.com/questions/32370485/r-convert-radians-to-degree-degree-to-radians
se <- function(x, na.rm = TRUE) {
  sd(x,na.rm = )/sqrt(sum(!is.na(x)))}

##original functions##
Create.var<- function(Lbody) {
  #this function creates a dataframe of relevant humpback whale measurements 
  
  #variables are defined in Potevin et al 2012. 
  
  #Lbody: the total body length of a humpback whale (m). 
  #Ljaw: length of humback whale's lowerjaw (m).
  #Whead: width of humpback whale's head (m).3
  #Mc: body mass
  #Abody: area of the humpback whale's body as it moves though the water
  #L0: length of ventral groove blubber
  #Xjd: humpback whale jaw disarticulation factor
  #Hsync: synchronization factor
  #Dwater: density of seawater
  
  #Allometric relationships from Table 3 in this paper

  Rname<-c("Lbody","Ljaw", "Whead", "Mc", "Abody", "L0", "Xjd","Hsync", "Dwater")
  value<- c(Lbody,
            0.134*Lbody^1.21,
            0.155*Lbody^1.04,
            0.7683*Lbody^4.17,
            0.00674*Lbody^2.75,
            0.376*Lbody^1.19,
            1.03,
            46.69,
            1025)
    assign("var",as.data.frame(cbind(Rname, value)))
    var$value<-as.numeric(as.character(var$value))
    var$Rname<-as.character(var$Rname)
    assign("var",var, envir = .GlobalEnv)}

Create.df<- function(var,Prey, Vc_0){
  #This function yields humpback whale speed, acceleration, distance and engulfment volume at each 0.01 sec increment of engulfment.  
  #Code is based off of Potvin et al 2012 in PLOS one.
  
  #var: dataframe of whale morpholoigcal parameters allometrically scaled
  #Prey: catagorical variable options are 2: "krill" or "fish" and are used to determine
  #Vc_0: initial speed of the whale
  
  #internal functions
  
  mouthopening<-function(){
    for (t in seq(0.01, Gapmax.t, dt)) {
      #THETAgape.t<- Hsync1*Vc_0/(sin(deg2rad(THETAgape))*Xjd*Ljaw)*t #gape angle at time t if opening EQ6
      THETAgape.t<-(THETAgape/Topen)*t
      .THETAgape.t<-(THETAgape/Topen)
      ifelse(t<=0.66*Topen, PHE<-1.6, PHE<-0 )
      Ac.t<-0.5*pi*Whead/2*(Xjd*Ljaw*sin(deg2rad(THETAgape.t))) # EQ11 mouth area in meters cubed (squared?)
      .Mw.t<-Dwater*Ac.t*(Vc.t-Vw.t)*PHE
      Mw.t<-Mw.t+.Mw.t*dt
      Fbc.t<-Kopen*4*Ac.t*Mw.t/(pi*Whead*Tengulf^2)# EQ16 Force applied to engulfed mass by buccal cavity walls
      Fed.t<-Fbc.t
      Fsd.t<-(Cd_open*Ac.t+Cd_body*Abody)*0.5*Dwater*Vc.t^2 #EQ20
      Fextthrust<- 4000*(Lbody-0.77)/(14.40-0.77)*(Lbody/14.40)^0.82 #Thrust plus bouyancy
      
      
      ac.t<-(Fextthrust-Fed.t-Fsd.t)/Mc
      Fww<-0
      aw.t<- (Fbc.t-Fww-Vw.t*.Mw.t)/Mw.t
      Vw.t<-Vw.t+aw.t*dt
      Vc.t<- Vc.t+ac.t*dt
      
      df[n,1]<-t # time
      df[n,2]<-Ac.t #Area of mouth at time = t
      df[n,3]<-Mw.t #Mass of engulfed water at time = t
      df[n,4]<-.Mw.t # filling rate of the mouth at time t
      df[n,5]<-Vc.t # speed of the whale at tiem t
      df[n,6]<-Vw.t # speed of the water at time t
      df[n,7]<-ac.t #acceleration of whale at time t
      df[n,8]<-aw.t # acceleration of water at time t
      df[n,9]<-THETAgape.t # angle of jaw at time t
      df[n,10]<-Fed.t
      df[n,11]<-Fsd.t
      df[n,12]<-Fextthrust
      df[n,13]<-.THETAgape.t
      n<-n+1
    }
    assign("df",df, envir = .GlobalEnv)
    assign("n",n, envir = .GlobalEnv)
    assign("Vc.t",Vc.t, envir = .GlobalEnv)
    assign("Vw.t",Vw.t, envir = .GlobalEnv)
    assign("ac.t",ac.t, envir = .GlobalEnv)
    assign("aw.t",aw.t, envir = .GlobalEnv)
    assign("Ac.t",Ac.t, envir = .GlobalEnv)
    assign("Mw.t",Mw.t, envir = .GlobalEnv)
    assign(".Mw.t",.Mw.t, envir = .GlobalEnv)
    assign("THETAgape.t",THETAgape.t, envir = .GlobalEnv)
    assign(".THETAgape.t",THETAgape.t, envir = .GlobalEnv)
    assign("Fed.t",Fed.t, envir = .GlobalEnv)
    assign("Fsd.t",Fsd.t, envir = .GlobalEnv)
    assign("Fextthrust",Fextthrust, envir = .GlobalEnv)
    
  }
  mouthopen<-function(){
    #*note Potvin et al 2012 does not have the mouth open phase.
    #This was added due to infomration from Cade et al 2016.
    
    for (t in seq(Gapmax.t+dt,Gapmax.t+GapemaxDur.mean , dt)) {
      THETAgape.t<-THETAgape #gape angle at time t if opening EQ6
      PHE<-0 
      Ac.t<-0.5*pi*Whead/2*(Xjd*Ljaw*sin(deg2rad(THETAgape.t))) # EQ11 mouth area in meters cubed (squared?)
      .Mw.t<-Dwater*Ac.t*(Vc.t-Vw.t)*PHE
      Mw.t<-Mw.t+.Mw.t*dt
      Fbc.t<-Kopen*4*Ac.t*Mw.t/(pi*Whead*Tengulf^2)# EQ16 Force applied to engulfed mass by buccal cavity walls
      Fed.t<-Fbc.t
      Fsd.t<-(Cd_open*Ac.t+Cd_body*Abody)*0.5*Dwater*Vc.t^2 #EQ20
      Fextthrust<- 4000*(Lbody-0.77)/(14.40-0.77)*(Lbody/14.40)^0.82 #Thrust plus bouyancy
      
      ac.t<-(Fextthrust-Fed.t-Fsd.t)/Mc
      Fww<-0
      aw.t<- (Fbc.t-Fww-Vw.t*.Mw.t)/Mw.t
      Vw.t<-Vw.t+aw.t*dt
      Vc.t<- Vc.t+ac.t*dt
      
      df[n,1]<-t # time
      df[n,2]<-Ac.t #Area of mouth at time = t
      df[n,3]<-Mw.t #Mass of engulfed water at time = t
      df[n,4]<-.Mw.t # filling rate of the mouth at time t
      df[n,5]<-Vc.t # speed of the whale at tiem t
      df[n,6]<-Vw.t # speed of the water at time t
      df[n,7]<-ac.t #acceleration of whale at time t
      df[n,8]<-aw.t # acceleration of water at time t
      df[n,9]<-THETAgape.t # angle of jaw at time t
      df[n,10]<-Fed.t
      df[n,11]<-Fsd.t
      df[n,12]<-Fextthrust
      n<-n+1
    }
    assign("df",df,envir = .GlobalEnv)
    assign("n",n, envir = .GlobalEnv)
    assign("Vc.t",Vc.t, envir = .GlobalEnv)
    assign("Vw.t",Vw.t, envir = .GlobalEnv)
    assign("ac.t",ac.t, envir = .GlobalEnv)
    assign("aw.t",aw.t, envir = .GlobalEnv)
    assign("Fed.t",Fed.t, envir = .GlobalEnv)
    assign("Fsd.t",Fsd.t, envir = .GlobalEnv)
    assign("Fextthrust",Fextthrust, envir = .GlobalEnv)
    assign("Mw.t",Mw.t, envir = .GlobalEnv)
    assign(".Mw.t",.Mw.t, envir = .GlobalEnv)
    assign("Ac.t",Ac.t, envir = .GlobalEnv)
    assign("THETAgape.t",THETAgape.t, envir = .GlobalEnv)
  }
  emptydf<-function(){
    assign("Kclose",Kopen/(Kopen__Kclose),envir = .GlobalEnv)
    assign("Vc.t",Vc_0,envir = .GlobalEnv)
    t<-dt
    assign("Vw.t",0,envir = .GlobalEnv)
    assign("Fww",0,envir = .GlobalEnv)#  "Ocean to engulfed mass" drag
    assign("Mw.t",0,envir = .GlobalEnv)
    assign("Ac.t",0,envir = .GlobalEnv)
    THETAgape.t<-0
    df <- data.frame(t=numeric(1),
                     Ac.t=numeric(1), 
                     Mw.t=numeric(1), 
                     .Mw.t=numeric(1),
                     Vc.t=numeric(1),
                     Vw.t=numeric(1),
                     ac.t=numeric(1),
                     aw.t=numeric(1),
                     THETAgape.t= numeric(1),
                     Fed.t= numeric(1),
                     Fsd.t= numeric(1),
                     Fextthrust.t = numeric(1))
    df[1,1]<-0
    df[1,2]<-Ac.t
    df[1,3]<-Mw.t
    df[1,4]<-0
    df[1,5]<-Vc.t
    df[1,6]<-Vw.t
    df[1,7]<-0
    df[1,8]<-0
    df[1,9]<-THETAgape.t
    assign("df",df, envir = .GlobalEnv)
    assign("n",2, envir = .GlobalEnv)
    assign("t",dt, envir = .GlobalEnv)
    assign("THETAgape.t",THETAgape.t, envir = .GlobalEnv)
    #Hsync<-(Xjd/Topen)*(Ljaw/Vc_0)*(sin(deg2rad(THETAgape))*THETAgape)
  }
  
  #import values
  var$values<-as.numeric(as.character(var$values))
  var$Rname<-as.character(var$Rname)
  assignVec(var$Rname,var$values, envir = .GlobalEnv)  
  
  #Calculate engulf duration (equ 8)and total engulfment capacity (equn 9 and 10)
  Tengulf<- (1/(THETAgape*sin(deg2rad(THETAgape)))*Ljaw/(L0-Ljaw)*Ljaw/Vc_0*sin(deg2rad(THETAgape))*THETAgape)*2  #EQ 8 with substitution for Hsync
  Mw_postTMJ<-Dwater*pi/3*(L0-Ljaw)*(Xjd*Ljaw*sin(deg2rad(THETAgape)))*(0.5*Whead) #Eq 9 engulfed mass post TMJ
  Mw_antTMJ<-Dwater*pi/3*(Ljaw)*(Xjd*Ljaw*sin(deg2rad(THETAgape)))*(0.5*Whead) #Eq 10
  Mw<-Mw_postTMJ+Mw_antTMJ
  
  if (Prey == "krill"){
    #from Cade et al 2016
    Gapmax.t<-Tengulf*0.35
    GapemaxDur.mean<-Tengulf*0.25
  }
  if (Prey == "fish"){
    #from Cade et al 2016
    Gapmax.t<-Tengulf*0.25
    GapemaxDur.mean<-Tengulf*0.44
  }
  
  Topen<-Gapmax.t

  assign("Tclose",GapemaxDur.mean+Gapmax.t, envir = .GlobalEnv)

  Kopen<-10.7
  emptydf()
  mouthopening()
  mouthopen()
  
  
  while (Vc.t-Vw.t > 0.0001|Vc.t-Vw.t < -0.0001 ){
    
    #this while loop solves for Kopen as shown in Potvin et al. 2012 
    Kopen<- Kopen+(Vc.t-Vw.t)*0.1
    emptydf()
    mouthopening()
    mouthopen()
    print(Vc.t-Vw.t)
  }
  
  #Mouthclosing stage
  for (t in seq(Gapmax.t+GapemaxDur.mean+dt,Tengulf, dt)) {
   THETAgape.t <-(-THETAgape)/(Tengulf-Tclose)*(t-Gapmax.t-GapemaxDur.mean-dt)+THETAgape
    .THETAgape.t<-(-THETAgape)/(Tengulf-Tclose)

    PHE<-0
    Ac.t<-0.5*pi*Whead/2*(Xjd*Ljaw*sin(deg2rad(THETAgape.t))) # EQ11 mouth area in meters cubed
    
    .Mw_antTMJ.t<-(-Mw_antTMJ/THETAgape)*(.THETAgape.t) #EQ12b
    .Mw.t<-.Mw_antTMJ.t
    Mw.t<-Mw.t+.Mw.t*dt
    Fbc.t<-Kclose*4*Ac.t*Mw.t/(pi*Whead*Topen)# EQ17 Force applied to engulfed mass by buccal cavity walls
    Fed.t<-Fbc.t
    Kam<-0.2
    Fsd.t<-(Cd_close*Ac.t+Cd_body*Abody)*(0.5*Dwater*Vc.t^2)+ 0.2*Dwater*(pi/12)*L0*((Whead/2)+Ljaw)^2*ac.t #EQ21

    ac.t<-(Fextthrust-Fed.t-Fsd.t)/Mc
    Fww<-Fed.t*(1+(Mw.t/Mc))-Vw.t*.Mw_antTMJ.t #EQ27
    aw.t<- (Fbc.t-Fww-Vw.t*.Mw.t)/Mw.t # 
    ac.t<-aw.t
    Fextthrust.t<-ac.t*Mc+Fsd.t+Fed.t
    Vw.t<-Vw.t+aw.t*dt
    Vc.t<- Vc.t+ac.t*dt
    
    df[n,1]<-t
    df[n,2]<-Ac.t
    df[n,3]<-Mw.t
    df[n,4]<-.Mw.t
    df[n,5]<-ifelse(Vc.t>0,Vc.t,0)
    df[n,6]<-ifelse(Vw.t>0,Vw.t,0)
    df[n,7]<-ac.t
    df[n,8]<-aw.t
    df[n,9]<-THETAgape.t
    df[n,10]<-Fed.t
    df[n,11]<-Fsd.t
    df[n,12]<-Fextthrust.t
    n<-n+1
  }
  df$d.t<-df$Vc.t*dt #distance traveld by the jaw hinge at time step t
  df$D.t<-NA #total distance traveld by jaw hinge at time step t
  df$Dcapture.area.t<-NA #total distance traveled by the capture area at time t
  for(i in 1:nrow(df)){
    df[i,]$D.t<-sum(df[1:i,]$d.t)
  }
  df$vol.t<-df$d.t*df$Ac.t
  df$Dlj.t<- df$D.t+Xjd*Ljaw*cos(df$THETAgape.t*3.14159/180)# total horizontal distance traveld by the lower jaw at time step t
  df$dlj.t<-NA
  for(i in 2:nrow(df)){
    df[i,]$dlj.t<-df[i,]$Dlj.t-df[i-1,]$Dlj.t
  }
  df$Vlj.t<-df$dlj.t/dt
  df[1,]$Vlj.t<-df[1,]$Vc.t
  df$Vcapture.area.t<-df$Vc.t
  df[df$t> Tclose,]$Vcapture.area.t<-df[df$t> Tclose,]$Vlj.t
  df$dcapture.area.t<-df$Vcapture.area.t*dt
  for(i in 1:nrow(df)){
    df[i,]$Dcapture.area.t<-sum(df[1:i,]$dcapture.area.t)
  }
  assign("df", df, envir = .GlobalEnv)
}



Capture.prop.stochastic<-function (Vp, Ddetect, Vp.sd, Ddetect.sd, df, lungetype, Tclose, var,grid = 15){
  #Vp = velocity of the prey.  units = m/s. Choose a burst rather than milling speed
  #Ddetect = Distance at which prey react to a predator (begin to swim away). units = m
  # Vp.sd = unbiased population standard deviation of Vp. units = m/s
  #Ddetect.sd = unbiased population standard deviation of Ddetect. units = m
  #df = output of the Create df. function.  A dataframe containing many hydrodynamic characteristics at each 0.01 sec time step within engulfment
  #lungetype = Catagorical factor variable.  Options are "subsurface, shallow, or group surface." THis allows for modeling different escape paths of prey 
  #Tclose = The time at which the whale's mouth begins to close.  Also an output of the Create.df
  #var = a dataframe with allometrically-scaled morphological parameters for engulfment
  #grid = the x by x dimentions of a grid used to model different possible positions of prey within the whale's mouth
  
  if (lungetype == "subsurface"){
    #df$capture.prob<-NA
    df$capture.prob2<-NA
    df$escape.dist<-NA
    df$escape.dist.sd1<-NA
    df$escape.dist.sd2<-NA
    for (r in 2:nrow(df)){
      
      print(r)
      r1<-Whead/2
      r2<-Xjd*Ljaw*sin(deg2rad(df[r,]$THETAgape.t))
      x<- seq (from = -r1, to = r1, length.out = grid)
      y<- seq (0, r2, length.out = grid)
      points<-expand.grid(x,y)
      #plot(points$Var1,points$Var2)
      points<-points[points$Var1^2/(r1^2)+ points$Var2^2/r2^2 <=1, ]
      points$dist2<-NA
      for (n in 1:nrow(points)){
        elipse.dist<-function(x){
          sqrt(
            (points[n,]$Var1-x)^2+
                 (points[n,]$Var2- sqrt(r2^2-((x^2*r2^2)/r1^2)))^2
                  )
        }
        points[n,]$dist2<-optimize(elipse.dist,c(-r1,r1))$objective
      }
            #ggplot()+geom_contour(aes(x = Var1, y = Var2, z = dist),data = points)+theme_classic()
      #hist(points$dist)
      df[r,]$escape.dist<- Vp*Ddetect/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t)
      df[r,]$escape.dist.sd1<- sqrt((Ddetect.sd^2*Vp.sd^2+Vp.sd^2*Ddetect^2+Ddetect.sd^2*Vp^2)*(1/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t))^2)
      #df[r,]$escape.dist.sd2<-sqrt((Vp.sd^2*Ddetect^2+Ddetect.sd^2*Vp^2+2*Cov*Ddetect*Vp)*(1/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t))^2)
      points$prob2<-pnorm(points$dist2, df[r,]$escape.dist, df[r,]$escape.dist.sd1)
      
      #plot 
      #p1<-ggplot(aes(x = Var1, y = -Var2),data = points)+geom_contour(aes(z = prob,colour = ..level..),size = 1,data = points)+scale_colour_gradient(limits=c(0, 1))+
       # theme_classic()+theme(aspect.ratio=max(points$Var2)/max(points$Var1*2),
                            #  panel.border = element_blank(),
                             # panel.grid.major = element_blank(),
                              #panel.grid.minor = element_blank(),
                              #axis.line.x = element_line(color="black", size = 0.5),
                              #axis.line.y = element_line(color="black", size = 0.5),
                              #panel.background = element_blank(),text = element_text(size=12))
       #direct.label(p1, list("bottom.pieces", colour='black'))
       #ggsave("StochasticMod2.tiff", path = "C:/Users/Ellen/Desktop/Analysis/Chapter3/", width =3.5, height = 4.5, unit = "in")
       
      #p2<-ggplot(aes(x = Var1, y = -Var2),data = points)+geom_point(aes(color = prob),size = 2)+
     #   theme_classic()+theme(aspect.ratio=max(points$Var2)/max(points$Var1*2),
     # panel.border = element_blank(),
     # panel.grid.major = element_blank(),
     # panel.grid.minor = element_blank(),
     # axis.line.x = element_line(color="black", size = 0.5),
     # axis.line.y = element_line(color="black", size = 0.5),
     # panel.background = element_blank(),text = element_text(size=12))+scale_colour_gradient(limits=c(0, 1))
      #ggsave("StochasticMod.tiff",p2, path = "C:/Users/Ellen/Desktop/Analysis/Chapter3/", width =3.5, height = 4.5, unit = "in")
      
      df[r,]$capture.prob2<-mean(points$prob2,na.rm = TRUE)##}
    }
    
    df2<-df[df$Dcapture.area.t< max(df$Dcapture.area.t)-df[nrow(df),]$escape.dist,]
    assign("prop_captured",sum(df2$Ac.t*df2$Vcapture.area.t*dt*df2$capture.prob2, na.rm = TRUE)/sum(df$Ac.t*df$Vcapture.area.t*dt, na.rm = TRUE), envir = .GlobalEnv)
    }
  if (lungetype == "shallow"){
    df$capture.prob2<-NA
    df$escape.dist<-NA
    df$escape.dist.sd1<-NA
    df$escape.dist.sd2<-NA
    for (r in 2:nrow(df)){
      print(r)
      r1<-Whead/2
      r2<-Xjd*Ljaw*sin(deg2rad(df[r,]$THETAgape.t))
      x<- seq (from = -r1, to = r1, length.out = grid)
      y<- seq (0, r2, length.out = grid)
      points<-expand.grid(x,y)
      #plot(points$Var1,points$Var2)
      points<-points[points$Var1^2/(r1^2)+ points$Var2^2/r2^2 <=1, ]
      points$dist2<-NA
      
      for (n in 1:nrow(points)){
        elipse.dist<-function(x){
          sqrt(
            (points[n,]$Var1-abs(x))^2+
              (points[n,]$Var2- sqrt(r2^2-((abs(x)^2*r2^2)/r1^2)))^2
          )
        }
        points[n,]$dist2<-optimize(elipse.dist,c(-r1,r1))$objective
      }
      #ggplot()+geom_contour(aes(x = Var1, y = Var2, z = dist),data = points)
      #hist(points$dist)
      df[r,]$escape.dist<- Vp*Ddetect/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t)
      df[r,]$escape.dist.sd1<- sqrt((Ddetect.sd^2*Vp.sd^2+Vp.sd^2*Ddetect^2+Ddetect.sd^2*Vp^2)*(1/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t))^2)
      #df[r,]$escape.dist.sd2<-sqrt((Vp.sd^2*Ddetect^2+Ddetect.sd^2*Vp^2+2*Cov*Ddetect*Vp)*(1/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t))^2)
      #points$prob<-pnorm(points$dist, df[r,]$escape.dist, df[r,]$escape.dist.sd1)
      points$prob2<-pnorm(points$dist2, df[r,]$escape.dist, df[r,]$escape.dist.sd1)
      
      #ggplot()+geom_contour(aes(x = Var1, y = Var2, z = dist),size = 1,data = points)+geom_point(aes(x = Var1, y = Var2),data = points)      
      df[r,]$capture.prob2<-mean(points$prob2,na.rm = TRUE)
      #p1<-ggplot(aes(x = -Var2, y = -Var1),data = points)+geom_contour(aes(z = prob,colour = ..level..),size = 1,data = points)+scale_colour_gradient(limits=c(0, 1))+
      #theme_classic()+theme(aspect.ratio=max(points$Var1*2)/max(points$Var2),
      #panel.border = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #axis.line.x = element_line(color="black", size = 0.5),
      #axis.line.y = element_line(color="black", size = 0.5),
      #panel.background = element_blank(),text = element_text(size=12))
#direct.label(p1, list("bottom.pieces", colour='black'))
#ggsave("StochasticMod3.tiff", path = "C:/Users/Ellen/Desktop/Analysis/Chapter3/", width =3.5, height = 4.5, unit = "in")

    #ggplot(aes(x =-Var2, y =  -Var1),data = points)+geom_point(aes(color = prob),size = 2)+theme_classic()+scale_colour_gradient(limits=c(0, 1))+theme(aspect.ratio=max(points$Var1*2)/max(points$Var2),
                                                                                                                 #panel.border = element_blank(),
                                                                                                                # panel.grid.major = element_blank(),
                                                                                                                # panel.grid.minor = element_blank(),
                                                                                                                # axis.line.x = element_line(color="black", size = 0.5),
                                                                                                                # axis.line.y = element_line(color="black", size = 0.5),
                                                                                                                # panel.background = element_blank(),text = element_text(size=12))                                                                                                          
#ggsave("StochasticMod4.tiff", path = "C:/Users/Ellen/Desktop/Analysis/Chapter3/", width =3.5, height = 4.5, unit = "in")
      }
    
    #df2<-df[df$Dcapture.area.t< max(df$Dcapture.area.t)-df[nrow(df),]$escape.dist,]
    df2<-df
    assign("prop_captured",sum(df2$Ac.t*df2$Vcapture.area.t*dt*df2$capture.prob, na.rm = TRUE)/sum(df$Ac.t*df$Vcapture.area.t*dt, na.rm = TRUE), envir = .GlobalEnv)
  }
  if(lungetype == "group lunge"){
    df$capture.prob<-1
    assign("prop_captured",1,envir = .GlobalEnv)
  }
  if(lungetype == "group subsurface"){
    df$escape.dist<-NA
    for (r in 2:nrow(df)){
    df[r,]$escape.dist<- Vp*Ddetect/mean(df[df$Dcapture.area.t< df[r,]$Dcapture.area.t & df$Dcapture.area.t > df[r,]$Dcapture.area.t - Ddetect,]$Vcapture.area.t)
    }
    df$capture.prob<-1
    df2<-df[df$Dcapture.area.t< max(df$Dcapture.area.t)-df[nrow(df),]$escape.dist,]
    assign("prop_captured",sum(df2$Ac.t*df2$Vcapture.area.t*dt*df2$capture.prob, na.rm = TRUE)/sum(df$Ac.t*df$Vcapture.area.t*dt, na.rm = TRUE), envir = .GlobalEnv)
  }
  assign("df",df,envir = .GlobalEnv)
}

  
##Read inputs##

var<-read.csv("WhaleParameters.csv")

EnergyParameters<-read.csv("PreyParameters.csv")
EnergyParameters$Prey<-as.character(EnergyParameters$Prey)

EnergyParameters$cov<-0.5*sqrt(EnergyParameters$SD.Prey.speed^2*EnergyParameters$SD.Reaction.distance^2)
EnergyParameters[5,]$Whale.initial.speed..subsurface.<-3.1
  
max<-nrow(EnergyParameters)

#without whale behavior
EnergyParameters$CaptureProp_NoWhaleBehave<-NA
for (a in 1:10){
  print(a)
    Create.df(var = var,Prey = EnergyParameters[a,]$Prey, Vc_0 = mean(EnergyParameters$Whale.initial.speed..subsurface., na.rm = TRUE))
    Capture.prop.stochastic(Vp = EnergyParameters[a,]$Prey.speed,
                            Vp.sd = EnergyParameters[a,]$SD.Prey.speed,
                            Ddetect = EnergyParameters[a,]$Mean.Reaction.distance,
                            Ddetect.sd = EnergyParameters[a,]$SD.Reaction.distance,
                            df = df,
                            lungetype="subsurface",
                            Tclose=Tclose,
                            var = var)
    EnergyParameters[a,]$CaptureProp_NoWhaleBehave<-prop_captured
  }

#with whale behavior
for (a in 1:10){
  print(a)
  if(!is.na(EnergyParameters[a,]$Whale.initial.speed..subsurface.)){
  Create.df(var = var,Prey = EnergyParameters[a,]$Prey, Vc_0 = EnergyParameters[a,]$Whale.initial.speed..subsurface.)
  Capture.prop.stochastic(Vp = EnergyParameters[a,]$Prey.speed,
                          Vp.sd = EnergyParameters[a,]$SD.Prey.speed,
                          Ddetect = EnergyParameters[a,]$Mean.Reaction.distance,
                          Ddetect.sd = EnergyParameters[a,]$SD.Reaction.distance,
                          df = df,
                          lungetype="subsurface",
                          Tclose=Tclose,
                          var = var)
  EnergyParameters[a,]$CaptureProp_subsurface<-prop_captured
  }
  if(!is.na(EnergyParameters[a,]$Whale.initial.speed..surface.)){
  Create.df(var = var,Prey = EnergyParameters[a,]$Prey, Vc_0 = EnergyParameters[a,]$Whale.initial.speed..surface.)
  Capture.prop.stochastic(Vp = EnergyParameters[a,]$Prey.speed,
                             Vp.sd = EnergyParameters[a,]$SD.Prey.speed,
                             Ddetect = EnergyParameters[a,]$Mean.Reaction.distance,
                             Ddetect.sd = EnergyParameters[a,]$SD.Reaction.distance,
                             df = df,
                             lungetype="shallow",
                             Tclose=Tclose,
                             var = var)
  EnergyParameters[a,]$CaptureProp_surface<-prop_captured
  }
}
EnergyParameters[2,]$CaptureProp_surface<-1

  Create.df(var = var,Prey = EnergyParameters[2,]$Prey, Vc_0 = EnergyParameters[2,]$Whale.initial.speed..subsurface.)
  Capture.prop.stochastic(Vp = EnergyParameters[2,]$Prey.speed,
                          Vp.sd = EnergyParameters[2,]$SD.Prey.speed,
                          Ddetect = EnergyParameters[2,]$Mean.Reaction.distance,
                          Ddetect.sd = EnergyParameters[2,]$SD.Reaction.distance,
                          df = df,
                          lungetype="group subsurface",
                          Tclose=Tclose,
                          var = var)
  EnergyParameters[2,]$CaptureProp_subsurface<-prop_captured

EnergyParameters[1:4,]$CaptureProportion_WA<-EnergyParameters[1:4,]$CaptureProp_surface*EnergyParameters[1:4,]$Proportion_surface+EnergyParameters[1:4,]$CaptureProp_subsurface*(1-EnergyParameters[1:4,]$Proportion_surface)
EnergyParameters[EnergyParameters$Proportion_surface ==1,]$CaptureProportion_WA<-EnergyParameters[EnergyParameters$Proportion_surface ==1,]$CaptureProp_surface
EnergyParameters[EnergyParameters$Proportion_surface ==0,]$CaptureProportion_WA<-EnergyParameters[EnergyParameters$Proportion_surface ==0,]$CaptureProp_subsurface

#duplicate rows where shallow and subsurface will be dealt with separately
EnergyParameters$depth<-NA
EnergyParameters[is.na(EnergyParameters$CaptureProp_subsurface),]$depth<-"shallow"
EnergyParameters[is.na(EnergyParameters$CaptureProp_surface),]$depth<-"subsurface"
EnergyParameters[!is.na(EnergyParameters$CaptureProp_surface)& !is.na(EnergyParameters$CaptureProp_subsurface),]$depth<-"shallow"
duplicate<-EnergyParameters[!is.na(EnergyParameters$CaptureProp_surface)& !is.na(EnergyParameters$CaptureProp_subsurface),]
duplicate$depth<-"subsurface"
EnergyParameters<-rbind(EnergyParameters,duplicate)
EnergyParameters$Prey.Group<-paste0(EnergyParameters$Prey.Group,"_",EnergyParameters$depth)
EnergyParameters$Data<-as.character(EnergyParameters$Data)
EnergyParameters[EnergyParameters$Data =="none",]$Data<-"CS2014Coho"

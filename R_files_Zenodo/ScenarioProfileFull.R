# Scenario setting -----------------------------------------------------

# Initialization for Scenario 
TT <- 78
ScenarioFlag1 <- substring(Scenario,1,3)
ScenarioFlag2 <- substring(Scenario,5,7)
ScenarioFlag3 <- substring(Scenario,9,11) 
ScenarioFlag4 <- substring(Scenario,13,15) 
ScenarioFlag5 <- substring(Scenario,17,19) 
ScenarioFlag7 <- substring(Scenario,29,30)
ScenarioFlag8 <- substring(Scenario,22,23) 


load("./Input/Labor_Cons_v6.RData")
load("./Input/Demand_Cons_v6.RData")

Scename <- substring(ScenarioFlag4,2,3)
Scename <- paste("./Input/202107/Scenario_v",Scename,".RData",sep = "") #Set as S30 structure is S25
load(Scename)

if(ScenarioFlag8 == "30")
{
   D_TES <- D_TES*(1-0.1);  D_LCS <- D_LCS*(1-0.1);  D_REA <- D_REA*(1-0.1)
}else if(ScenarioFlag8 == "50")
{
   D_TES <- D_TES*(1+0.1);  D_LCS <- D_LCS*(1+0.1);  D_REA <- D_REA*(1+0.1)
}else{

}

unsub <- c(91,96:112,120:127)#viechle prodution,electricity,96:112,transport,120:127,[91,96:112,120:127]

tr <- 35

Time_ED <- unlist(read_excel('./Data/parameter_data_v3.xlsx',sheet = 'Labor-Demand_p',"BI4:BI53",col_names = T))

if(ScenarioFlag3 == "D06")
{
  Time_s <- Time_ED+1+0
}else if(ScenarioFlag3 == "D08")
{
  Time_s <- Time_ED+1+10
}else if(ScenarioFlag3 == "D10")
{
  Time_s <- Time_ED+1+15
}else if(ScenarioFlag3 == "D12")
{
  Time_s <- Time_ED+1+6
}else
{
  Time_s <- Time_ED
}
# Scenario for supply (Labor Constrains) ----------------------------------

# Labor Definition -------------------------------------------------------------

SectLbSk_Label = c("SectLbSkL","SectLbSkM","SectLbSkH")

SectLbSkL = c(1:34,55,90,96:112,132:138)#0.1
SectLbSkM = c(56:85,127:130,135,139:158)#0.5
SectLbSkH = c(35:54,86:89,91:95,113:126,131,159:NN)#1

SectLbSk = list("SectLbSkL" = SectLbSkL,"SectLbSkM" = SectLbSkM,"SectLbSkH" = SectLbSkH)

RateLbSk = c("SectLbSkL" = 0.1,"SectLbSkM" = 0.5,"SectLbSkH" = 1)

Labor_Cons = array(1,dim = c(NN,RR,TT))

Labor_TT = cbind(Labor_r[,1:tr],matrix(0,nrow = RR, ncol = (TT-tr)))

temp  <- matrix(0,nrow <- RR, ncol <- 1)

for(j in 1:RR)
{
  for(t in tr:TT)
  {
    if(Labor_TT[j,t] == 0)
    {
      temp[j,] <- t
      break
    }
  }
}

temp[31,] <- 6

if(ScenarioFlag6 == "N")
{
  temp <- matrix(Time_s,ncol = 1,byrow = T)
  temp[31,] <- 6
}else{
  temp <- temp
}
# Labor Shock -------------------------------------------------------------
for(Sect in SectLbSk_Label)
{
  for(j in 1:RR)
  {
    for(i in SectLbSk[[Sect]])
    {
      tt <- temp[j,1]
      
      Labor_Cons[i,j,1:(tt-1)] <- 1+RateLbSk[Sect]*Labor_TT[j,1:(tt-1)]/100
    }
  }
}

# Labor recover -----------------------------------------------------------

SectLbRc_Label = c("SectLbRcL","SectLbRcM","SectLbRcH")

SectLbRcL = c(1:34,55,90,96:112,132:138)#0.1
SectLbRcM = c(56:85,127:130,135,139:158)#0.5
SectLbRcH = c(35:54,86:89,91:95,113:126,131,159:NN)#1

SectLbRc = list("SectLbRcL" = SectLbRcL,"SectLbRcM" = SectLbRcM,"SectLbRcH" = SectLbRcH)

if(ScenarioFlag2 == "LR2")
{
  RateLbRc_LR1 = c("SectLbRcL" = 2,"SectLbRcM" = 2,"SectLbRcH" = 2)
}else if(ScenarioFlag2 == "LR4")
{
  RateLbRc_LR1 = c("SectLbRcL" = 4,"SectLbRcM" = 4,"SectLbRcH" = 4)
}else if(ScenarioFlag2 == "LR8")
{
  RateLbRc_LR1 = c("SectLbRcL" = 8,"SectLbRcM" = 8,"SectLbRcH" = 8)
}else
{
  RateLbRc_LR1 = c("SectLbRcL" = 6,"SectLbRcM" = 6,"SectLbRcH" = 6)
}

TimeLag_LT1 = c("SectLbRcL" = 0,"SectLbRcM" = 0,"SectLbRcH" = 0)
TimeLag_LT2 = c("SectLbRcL" = 0,"SectLbRcM" = 0,"SectLbRcH" = 0)
TimeLag_LT3 = c("SectLbRcL" = 0,"SectLbRcM" = 0,"SectLbRcH" = 0)

SectLbRc_s = list("SectLbRc_TES" = L_TES,"SectLbRc_LCS" = L_LCS,
                  "SectLbRc_LDS" = L_LDS,"SectLbRc_REA" = L_REA)

RateLbRc_s = c("SectLbRc_TES" = as.numeric(RateLbRc_LR1[1]),
               "SectLbRc_LCS" = as.numeric(RateLbRc_LR1[1]),
               "SectLbRc_DES" = as.numeric(RateLbRc_LR1[1]),
               "SectLbRc_LDS" = as.numeric(RateLbRc_LR1[1]),
               "SectLbRc_REA" = as.numeric(RateLbRc_LR1[1]))#

TimeLagLb_s = c("SectLbRc_TES" = 0,"SectLbRc_LCS" = 0,"SectLbRc_DES" = 0,
                "SectLbRc_LDS" = 0,"SectLbRc_REA" = 0)

if(ScenarioFlag5 =="T12")
{
  TopLb_s = c("SectLbRc_TES" = 1.2,
              "SectLbRc_LCS" = 1.2,
              "SectLbRc_DES" = 1.2,
              "SectLbRc_LDS" = 1.2,
              "SectLbRc_REA" = 1.2)
  
}else if(ScenarioFlag5 =="T10")
{
  TopLb_s = c("SectLbRc_TES" = 1.01,
              "SectLbRc_LCS" = 1.01,
              "SectLbRc_DES" = 1.01,
              "SectLbRc_LDS" = 1.01,
              "SectLbRc_REA" = 1.01)
}else if(ScenarioFlag5 =="T11")
{
  TopLb_s = c("SectLbRc_TES" = 1.1,
              "SectLbRc_LCS" = 1.1,
              "SectLbRc_DES" = 1.1,
              "SectLbRc_LDS" = 1.1,
              "SectLbRc_REA" = 1.1)
}else if(ScenarioFlag5 =="T15")
{
  TopLb_s = c("SectLbRc_TES" = 1.5,
              "SectLbRc_LCS" = 1.5,
              "SectLbRc_DES" = 1.5,
              "SectLbRc_LDS" = 1.5,
              "SectLbRc_REA" = 1.5)
}else
{
  TopLb_s = c("SectLbRc_TES" = 2,
              "SectLbRc_LCS" = 2,
              "SectLbRc_DES" = 2,
              "SectLbRc_LDS" = 2,
              "SectLbRc_REA" = 2)
}

for(Sect in SectLbRc_Label)
{
  for(i in SectLbRc[[Sect]])
  {
    for(j in 1:RR)
    {
      tt <- temp[j,1]
      
      temp_a <- Labor_Cons[i,j,(tt-1)]#
      
      Lb <- Labor_Cons[i,j,(tt-1)]#
      
      Tl <- TimeLag_LT1[Sect] #LT
      
      for(t in tt:TT)
      {
        if((t <= (tt+Tl-1))&&(Lb<1.01))#
        {
          Labor_Cons[i,j,t] <- Lb
          
        }else if((t > (tt+Tl-1))&&(Lb<1.01))#
        {
          Lb <- temp_a+RateLbRc_LR1[Sect]*(t-tt-Tl)/100#LR
          
          Labor_Cons[i,j,t] <- Lb
          
        }else
        {
          Labor_Cons[i,j,t] <- 1.01
        }
      }
    }
  }
}

if(ScenarioFlag1 != "BAU")
{
  Sect_s <- paste("SectLbRc_",ScenarioFlag1,sep = "")
  
  for(j in c(1:RR))#1:RR
  {
    for(i in SectLbRc_s[[Sect_s]][[j]])
    {
      tt <- temp[j,1]
      
      temp_a <- Labor_Cons[i,j,(tt-1)]#
      
      Lb <- Labor_Cons[i,j,(tt-1)]#
      
      Tl <- TimeLagLb_s[Sect_s] #LT
      
      for(t in tt:TT)
      {
        if((t <= (tt+Tl-1))&&(Lb<TopLb_s[Sect_s]))#
        {
          Labor_Cons[i,j,t] <- Lb
          
        }else if((t > (tt+Tl-1))&&(Lb<TopLb_s[Sect_s]))#
        {
          Lb <- temp_a+RateLbRc_s[Sect_s]*(t-tt-Tl)/100#LR
          
          Labor_Cons[i,j,t] <- Lb
          
        }else
        {
          Lb <-TopLb_s[Sect_s]
          
          Labor_Cons[i,j,t] <- Lb
        }
      }
    }
  }
}

# Scenario for Demand (Demand decrease/increase) ----------------------------------

# Demand Definition --------------------------------------------------------------

Demand_Cons =  array(0,dim =  c(NNRR,RRFF,TT))

Demand_Shk =  array(0,dim =  c(NN,RRFF,TT))

Demand_TT =  array(0,dim =  c(NN,RR,TT))

IOF_TT0 =  array(rep(IOF_0,TT),dim =  c(NNRR,RRFF,TT))

IOF_TT =  IOF_TT0

# Demand Shock -------------------------------------------------------------

SectDmSk_Label =  c("SectDmSk_RR","SectDmSk_TS","SectDmSk_RS","SectDmSk_MD","SectDmSk_OT","SectDmSk_PK")

SectDmSk_RR =  c(115:119) # retail_and_recreation_20:34,43:44,56:95,113:119,127:135,138
SectDmSk_PK =  c(126,160) # parks
SectDmSk_TS =  c(120:125) # transit
SectDmSk_RS =  c(96:112) # residential(electricity demand +)
SectDmSk_MD =  c(90,138) # health

SectDmSk_AL =  c(1:NN)
SectDmSk_OT =  SectDmSk_AL[-c(SectDmSk_RR,SectDmSk_PK,SectDmSk_TS,SectDmSk_RS,SectDmSk_MD)]

# ... Demand real/policy timetable input

Demand_r_MD =  Demand_r_RS

Demand_r_OT =  (Demand_r_RR+Demand_r_TS)/2

Demand_r_OT[31,] = Demand_r_OT_CN

SectDmSk =  list("SectDmSk_RR" =  SectDmSk_RR,"SectDmSk_PK" =  SectDmSk_PK,"SectDmSk_TS" =  SectDmSk_TS,
                 "SectDmSk_RS" =  SectDmSk_RS,"SectDmSk_MD" =  SectDmSk_MD,"SectDmSk_OT" =  SectDmSk_OT)

Demand_r =  list("SectDmSk_RR" =  Demand_r_RR,"SectDmSk_PK" =  Demand_r_PK,"SectDmSk_TS" =  Demand_r_TS,
                 "SectDmSk_RS" =  Demand_r_RS,"SectDmSk_MD" =  Demand_r_MD,"SectDmSk_OT" =  Demand_r_OT)

# ... Demand timetable & Shock

for(Sect in SectDmSk_Label)
{
  Demand_rp <- cbind(Demand_r[[Sect]][,1:tr],matrix(0,nrow <- RR, ncol <- (TT-tr)))
  
  for(i in SectDmSk[[Sect]])
  {
    
    Demand_rp <- data.matrix(Demand_rp)
    
    Demand_TT[i,,] <- Demand_rp
    
    temp <- matrix(0, nrow <- RR, ncol <- 1)
    
    for(j in 1:RR)
    {
      for(t in tr:TT)
      {
        if(Demand_TT[i,j,t] == 0)
        {
          temp[j,] <- t
          break
        }
      }
      
      temp[31,] <- 11
      
      if(ScenarioFlag6 == "N")
      {
        temp <- matrix(Time_s,ncol = 1,byrow = T)
      }else{
        temp <- temp
      }
      
      tt <- temp[j,1]
      
      for(f in (2:3))
      {
        a <- (j-1)*FF+f
        
        Demand_Shk[i,a,1:(tt-1)] <- Demand_TT[i,j,1:(tt-1)]
      }
    }
  }
}


# Demand Recover -----------------------------------------------------------
Time_ED <- unlist(read_excel('./Data/parameter_data_v3.xlsx',sheet = 'Labor-Demand_p',"BI4:BI53",col_names = T))

if(ScenarioFlag3 == "D06")
{
  Time_s <- Time_ED+1+0
}else if(ScenarioFlag3 == "D08")
{
  Time_s <- Time_ED+1+10
}else if(ScenarioFlag3 == "D10")
{
  Time_s <- Time_ED+1+15
}else if(ScenarioFlag3 == "D12")
{
  Time_s <- Time_ED+1+6
}else
{
  Time_s <- Time_ED
}

SectDmRcC_Label = c("SectDmRcC1","SectDmRcC2","SectDmRcC3")
SectDmRcI_Label = c("SectDmRcI1","SectDmRcI2","SectDmRcI3")

SectDmRcC1 = c(SectDmSk_RR,SectDmSk_PK,SectDmSk_TS)
SectDmRcC2 = c(SectDmSk_RS,SectDmSk_MD)
SectDmRcC3 = c(SectDmSk_OT)
SectDmRcI1 = c(SectDmSk_RR,SectDmSk_PK,SectDmSk_TS)
SectDmRcI2 = c(SectDmSk_RS,SectDmSk_MD)
SectDmRcI3 = c(SectDmSk_OT)

SectDmRcC = list("SectDmRcC1" = SectDmRcC1,"SectDmRcC2" = SectDmRcC2,"SectDmRcC3" = SectDmRcC3)
SectDmRcI = list("SectDmRcI1" = SectDmRcC1,"SectDmRcI2" = SectDmRcI2,"SectDmRcI3" = SectDmRcI3)
RateDmRcC = c("SectDmRcC1" = 2,"SectDmRcC2" = 2,"SectDmRcC3" = 2)
RateDmRcI = c("SectDmRcI1" = 2,"SectDmRcI2" = 2,"SectDmRcI3" = 2)

TimeLagC = c("SectDmRcC1" = 0,"SectDmRcC2" = 0,"SectDmRcC3" = 0)
TimeLagI = c("SectDmRcI1" = 0,"SectDmRcI2" = 0,"SectDmRcI3" = 0)

num1 <- seq(1,97,2)
num2 <- seq(2,98,2)

TopDmC_TES = as.matrix(D_TES[,num1])
TopDmC_LCS = as.matrix(D_LCS[,num1])
TopDmC_LDS = as.matrix(D_LDS[,num1])
TopDmC_REA = as.matrix(D_REA[,num1])

TopDmC_s = list("TopDmC_TES" = TopDmC_TES,"TopDmC_LCS" = TopDmC_LCS,
                "TopDmC_LDS" = TopDmC_LDS,"TopDmC_REA" = TopDmC_REA)

TopDmI_TES = as.matrix(D_TES[,num2])
TopDmI_LCS = as.matrix(D_LCS[,num2])
TopDmI_LDS = as.matrix(D_LDS[,num2])
TopDmI_REA = as.matrix(D_REA[,num2])


TopDmI_s = list("TopDmI_TES" = TopDmI_TES,"TopDmI_LCS" = TopDmI_LCS,
                "TopDmI_LDS" = TopDmI_LDS,"TopDmI_REA" = TopDmI_REA)

for(f in 2:3)
{
  if(f == 2)
  {
    for(Sect in SectDmRcC_Label)
    {
      for(i in SectDmRcC[[Sect]])
      {
        for(j in 1:RR)
        {
          a = (j-1)*FF+f
          
          tt <- temp[j,1]#
          
          temp_a = Demand_TT[i,j,(tt-1)]
          
          Dm <- Demand_TT[i,j,(tt-1)]
          
          Tl <- TimeLagC[Sect]
          
          ta <- tt+Tl-1
          
          tb <- tt+Tl-1+(5-temp_a)/RateDmRcC[Sect]
          
          if(temp_a >= 0)
          {
            for(t in tt:TT)
            {
              if(Dm > 0)
              {
                Dm <- 0
                
                Demand_Shk[i,a,t] <- Dm
              }else
              {
                Dm <- 0
                
                Demand_Shk[i,a,t] <- Dm
              }
            }
          }else
          {
            for(t in tt:TT)
            {
              if(t <= ta)
              {
                Demand_Shk[i,a,t] <- Dm
              }else if((t > ta) && (t <= tb))
              {
                if(Dm < 0)
                {
                  Dm <- Dm + RateDmRcC[Sect]
                  
                  Demand_Shk[i,a,t] <- Dm
                  
                }else
                {
                  Dm <- 0
                  
                  Demand_Shk[i,a,t] <- Dm
                }
              }else
              {
                if(Dm >= 0)
                {
                  Dm <- 0
                  
                  Demand_Shk[i,a,t] <- Dm
                }else{break}
              }
            }
          }
        }
      }
    }
  }
  if(f == 3)
  {
    for(Sect in SectDmRcI_Label)
    {
      for(i in SectDmRcI[[Sect]])
      {
        for(j in 1:RR)
        {
          a = (j-1)*FF+f
          
          tt <- temp[j,1]
          
          temp_a = Demand_TT[i,j,(tt-1)]
          
          Dm <- Demand_TT[i,j,(tt-1)]
          
          Tl <- TimeLagI[Sect]
          
          ta <- tt+Tl-1
          
          tb <- tt+Tl-1+(5-temp_a)/RateDmRcI[Sect]
          
          if(temp_a >= 0)
          {
            for(t in tt:TT)
            {
              if(Dm > 0)
              {
                Dm <- 0
                
                Demand_Shk[i,a,t] <- Dm
              }else
              {
                Dm <- 0
                
                Demand_Shk[i,a,t] <- Dm
              }
            }
          }else
          {
            for(t in tt:TT)
            {
              if(t <= ta)
              {
                Demand_Shk[i,a,t] <- Dm
              }else if((t > ta) && (t <= tb))
              {
                if(Dm < 0)
                {
                  Dm <- Dm + RateDmRcI[Sect]
                  
                  Demand_Shk[i,a,t] <- Dm
                  
                }else
                {
                  Dm <- 0
                  
                  Demand_Shk[i,a,t] <- Dm
                }
              }else
              {
                if(Dm >= 0)
                {
                  Dm <- 0
                  
                  Demand_Shk[i,a,t] <- Dm
                }else{break}
              }
            }
          }
        }
      }
    }
  }
}

for(i in 1:NN)
{
  for(j in 1:RR)
  {
    Demand_Cons[(j-1)*NN+i,,] <- Demand_Shk[i,,]
  }
}

IOF_TT = IOF_TT0*(1+Demand_Cons/100)
IOF_TT_BAU = IOF_TT


if(ScenarioFlag1 != "BAU")
{
  
  for(j in c(1:RR))#1:RR
  {
    for(f in 2:3)
    {
      if(f == 2)
      { 
        nameS_s <- paste("TopDmC_",ScenarioFlag1,sep = "")
        
        a = (j-1)*FF+f
        
        for(i in c(1:NN))
        {
          c <- (j-1)*NN+i
          
          Top <- TopDmC_s[[nameS_s]][i,j]/(TT-Time_s[j]+1)
          
          if(Top != 0)
          {
            IOF_TT[c,a,Time_s[j]:TT] <- IOF_TT[c,a,Time_s[j]:TT]+Top
          }else
          {
            IOF_TT[c,a,Time_s[j]:TT] <- IOF_TT[c,a,Time_s[j]:TT]
          }
          
          if(j == 31)
          {
            cat("IOF_Rate_C"," ",i,sum(IOF_TT[c,a,Time_s[j]:TT]-IOF_TT_BAU[c,a,Time_s[j]:TT])/sum(IOF_0[c,a]*26),"\\n")
          }
        }
      }
      
      if(f == 3)
      { 
        nameS_s <- paste("TopDmI_",ScenarioFlag1,sep = "")
        
        a = (j-1)*FF+f
        
        for(i in c(1:NN))
        {
          c <- (j-1)*NN+i
          
          Top <- TopDmI_s[[nameS_s]][i,j]/(TT-Time_s[j]+1)
          
          if(Top != 0)
          {
            IOF_TT[c,a,Time_s[j]:TT] <- IOF_TT[c,a,Time_s[j]:TT]+Top
          }else
          {
            IOF_TT[c,a,Time_s[j]:TT] <- IOF_TT[c,a,Time_s[j]:TT]
          }
          
          if(j == 31)
          {
            cat("IOF_Rate_I"," ",i,sum(IOF_TT[c,a,Time_s[j]:TT]-IOF_TT_BAU[c,a,Time_s[j]:TT])/sum(IOF_0[c,a]*26),"\\n")
          }
        }
      }
    }
  }
}


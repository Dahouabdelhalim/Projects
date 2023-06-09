library("multcomp")

## GEN1
GEN1<-read.csv('[pathname]/GEN1colonies.csv',h=T)
GEN1$EXP<-factor(GEN1$Treatment) # 3 levels

  # GEN1 specialization
  gen1_m1<-glm(Specialization~EXP,data=GEN1)
  shapiro.test(gen1_m1$resid)
  drop1(gen1_m1,test="Chisq")
  summary(glht(gen1_m1, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))

  # GEN1 behavioral variation
  gen1_m2<-glm(std_RMSD^0.5~EXP,data=GEN1)
  shapiro.test(gen1_m2$resid) 
  drop1(gen1_m2,test="Chisq")
  summary(glht(gen1_m2, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))
  
  
## GEN2 
GEN2<-read.csv('[pathname]/GEN2colonies.csv',h=T)
GEN2$EXP<-factor(GEN2$Treatment) # 3 levels

  # GEN2 specialization
  gen2_m1<-glm(Specialization~EXP,data=GEN2)
  shapiro.test(gen2_m1$resid)
  drop1(gen2_m1,test="Chisq")
  summary(glht(gen2_m1, linfct=mcp(EXP="Tukey")),test=adjusted("holm")) 
 
  # GEN2 behavioral variation
  gen2_m2<-glm(std_RMSD~EXP,data=GEN2)
  shapiro.test(gen2_m2$resid) 
  drop1(gen2_m2,test="Chisq")
  summary(glht(gen2_m2, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))

  
## IW 
IW<-read.csv('[pathname]/IWcolonies.csv',h=T)
IW$EXP<-factor(IW$Treatment) # 3 levels
  
  # IW specialization
  IW_m1<-glm(Specialization~EXP,data=IW)
  shapiro.test(IW_m1$resid)
  drop1(IW_m1,test="Chisq")
  summary(glht(IW_m1, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))

  # IW behavioral variation
  IW_m2<-glm(std_RMSD~EXP,data=IW)
  shapiro.test(IW_m2$resid) 
  drop1(IW_m2,test="Chisq")
  summary(glht(IW_m2, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))

  
## AGE 
AGE<-read.csv('[pathname]/AGEcolonies.csv',h=T)
AGE$EXP<-factor(AGE$Treatment) # 3 levels

  # AGE specialization
  AGE_m1<-glm(Specialization~EXP,data=AGE)
  shapiro.test(AGE_m1$resid)
  drop1(AGE_m1,test="Chisq")
  summary(glht(AGE_m1, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))
  
  # AGE behavioral variation
  AGE_m2<-glm(std_RMSD~EXP,data=AGE)
  shapiro.test(AGE_m2$resid) 
  drop1(AGE_m2,test="Chisq")
  summary(glht(AGE_m2, linfct=mcp(EXP="Tukey")),test=adjusted("holm"))
 
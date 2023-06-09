# RUN AFTER 4_TableS3-S7.R

t=coeftest(md_pf, vcov = vcovHC(md_pf, "HC1"))

var1_pf <- data.frame(Variable = "Advice Doctor",
                      Coefficient = t[2,1],
                      SE = t[2,2], 
                      p = t[2,4])
var2_pf <- data.frame(Variable = "Advice Scientist",
                      Coefficient = t[3,1],
                      SE = t[3,2], 
                      p = t[3,4])
var3_pf <- data.frame(Variable = "Advice Anti-Vaccine Propagist",
                      Coefficient = t[4,1],
                      SE = t[4,2], 
                      p = t[4,4])
var4_pf <- data.frame(Variable = "Advice Politicians",
                      Coefficient = t[5,1],
                      SE = t[5,2], 
                      p = t[5,4])
var5_pf <- data.frame(Variable = "Advice Family",
                      Coefficient = t[6,1],
                      SE = t[6,2], 
                      p = t[6,4])
var6_pf <- data.frame(Variable = "Advice Friends",
                      Coefficient = t[7,1],
                      SE = t[7,2], 
                      p = t[7,4])
var7_pf <- data.frame(Variable = "Advice Journalists",
                      Coefficient = t[8,1],
                      SE = t[8,2], 
                      p = t[8,4])
var8_pf <- data.frame(Variable = "Advice Celebrities",
                      Coefficient = t[9,1],
                      SE = t[9,2], 
                      p = t[9,4])
var9_pf <- data.frame(Variable = "Age",
                      Coefficient = t[10,1],
                      SE = t[10,2], 
                      p = t[10,4])
var10_pf <- data.frame(Variable = "Male",
                       Coefficient = t[11,1],
                       SE = t[11,2], 
                       p = t[11,4])
var11_pf <- data.frame(Variable = "Schooling",
                       Coefficient = t[12,1],
                       SE = t[12,2], 
                       p = t[12,4])
var12_pf <- data.frame(Variable = "City Category by Size",
                       Coefficient = t[13,1],
                       SE = t[13,2], 
                       p = t[13,4])
var13_pf <- data.frame(Variable = "Smoking",
                       Coefficient = t[14,1],
                       SE = t[14,2], 
                       p = t[14,4])
var14_pf <- data.frame(Variable = "Chronic Illness",
                       Coefficient = t[15,1],
                       SE = t[15,2], 
                       p = t[15,4])
var15_pf <- data.frame(Variable = "COVID Previously",
                       Coefficient = t[16,1],
                       SE = t[16,2], 
                       p = t[16,4])
var16_pf <- data.frame(Variable = "Serious COVID Previously",
                       Coefficient = t[17,1],
                       SE = t[17,2], 
                       p = t[17,4])

t=coeftest(md_mo, vcov = vcovHC(md_mo, "HC1"))

var1_mo <- data.frame(Variable = "Advice Doctor",
                      Coefficient = t[2,1],
                      SE = t[2,2], 
                      p = t[2,4])
var2_mo <- data.frame(Variable = "Advice Scientist",
                      Coefficient = t[3,1],
                      SE = t[3,2], 
                      p = t[3,4])
var3_mo <- data.frame(Variable = "Advice Anti-Vaccine Propagist",
                      Coefficient = t[4,1],
                      SE = t[4,2], 
                      p = t[4,4])
var4_mo <- data.frame(Variable = "Advice Politicians",
                      Coefficient = t[5,1],
                      SE = t[5,2], 
                      p = t[5,4])
var5_mo <- data.frame(Variable = "Advice Family",
                      Coefficient = t[6,1],
                      SE = t[6,2], 
                      p = t[6,4])
var6_mo <- data.frame(Variable = "Advice Friends",
                      Coefficient = t[7,1],
                      SE = t[7,2], 
                      p = t[7,4])
var7_mo <- data.frame(Variable = "Advice Journalists",
                      Coefficient = t[8,1],
                      SE = t[8,2], 
                      p = t[8,4])
var8_mo <- data.frame(Variable = "Advice Celebrities",
                      Coefficient = t[9,1],
                      SE = t[9,2], 
                      p = t[9,4])
var9_mo <- data.frame(Variable = "Age",
                      Coefficient = t[10,1],
                      SE = t[10,2], 
                      p = t[10,4])
var10_mo <- data.frame(Variable = "Male",
                       Coefficient = t[11,1],
                       SE = t[11,2], 
                       p = t[11,4])
var11_mo <- data.frame(Variable = "Schooling",
                       Coefficient = t[12,1],
                       SE = t[12,2], 
                       p = t[12,4])
var12_mo <- data.frame(Variable = "City Category by Size",
                       Coefficient = t[13,1],
                       SE = t[13,2], 
                       p = t[13,4])
var13_mo <- data.frame(Variable = "Smoking",
                       Coefficient = t[14,1],
                       SE = t[14,2], 
                       p = t[14,4])
var14_mo <- data.frame(Variable = "Chronic Illness",
                       Coefficient = t[15,1],
                       SE = t[15,2], 
                       p = t[15,4])
var15_mo <- data.frame(Variable = "COVID Previously",
                       Coefficient = t[16,1],
                       SE = t[16,2], 
                       p = t[16,4])
var16_mo <- data.frame(Variable = "Serious COVID Previously",
                       Coefficient = t[17,1],
                       SE = t[17,2], 
                       p = t[17,4])

t=coeftest(md_as, vcov = vcovHC(md_as, "HC1"))

var1_as <- data.frame(Variable = "Advice Doctor",
                      Coefficient = t[2,1],
                      SE = t[2,2], 
                      p = t[2,4])
var2_as <- data.frame(Variable = "Advice Scientist",
                      Coefficient = t[3,1],
                      SE = t[3,2], 
                      p = t[3,4])
var3_as <- data.frame(Variable = "Advice Anti-Vaccine Propagist",
                      Coefficient = t[4,1],
                      SE = t[4,2], 
                      p = t[4,4])
var4_as <- data.frame(Variable = "Advice Politicians",
                      Coefficient = t[5,1],
                      SE = t[5,2], 
                      p = t[5,4])
var5_as <- data.frame(Variable = "Advice Family",
                      Coefficient = t[6,1],
                      SE = t[6,2], 
                      p = t[6,4])
var6_as <- data.frame(Variable = "Advice Friends",
                      Coefficient = t[7,1],
                      SE = t[7,2], 
                      p = t[7,4])
var7_as <- data.frame(Variable = "Advice Journalists",
                      Coefficient = t[8,1],
                      SE = t[8,2], 
                      p = t[8,4])
var8_as <- data.frame(Variable = "Advice Celebrities",
                      Coefficient = t[9,1],
                      SE = t[9,2], 
                      p = t[9,4])
var9_as <- data.frame(Variable = "Age",
                      Coefficient = t[10,1],
                      SE = t[10,2], 
                      p = t[10,4])
var10_as <- data.frame(Variable = "Male",
                       Coefficient = t[11,1],
                       SE = t[11,2], 
                       p = t[11,4])
var11_as <- data.frame(Variable = "Schooling",
                       Coefficient = t[12,1],
                       SE = t[12,2], 
                       p = t[12,4])
var12_as <- data.frame(Variable = "City Category by Size",
                       Coefficient = t[13,1],
                       SE = t[13,2], 
                       p = t[13,4])
var13_as <- data.frame(Variable = "Smoking",
                       Coefficient = t[14,1],
                       SE = t[14,2], 
                       p = t[14,4])
var14_as <- data.frame(Variable = "Chronic Illness",
                       Coefficient = t[15,1],
                       SE = t[15,2], 
                       p = t[15,4])
var15_as <- data.frame(Variable = "COVID Previously",
                       Coefficient = t[16,1],
                       SE = t[16,2], 
                       p = t[16,4])
var16_as <- data.frame(Variable = "Serious COVID Previously",
                       Coefficient = t[17,1],
                       SE = t[17,2], 
                       p = t[17,4])

t=coeftest(md_sp, vcov = vcovHC(md_sp, "HC1"))

var1_sp <- data.frame(Variable = "Advice Doctor",
                      Coefficient = t[2,1],
                      SE = t[2,2], 
                      p = t[2,4])
var2_sp <- data.frame(Variable = "Advice Scientist",
                      Coefficient = t[3,1],
                      SE = t[3,2], 
                      p = t[3,4])
var3_sp <- data.frame(Variable = "Advice Anti-Vaccine Propagist",
                      Coefficient = t[4,1],
                      SE = t[4,2], 
                      p = t[4,4])
var4_sp <- data.frame(Variable = "Advice Politicians",
                      Coefficient = t[5,1],
                      SE = t[5,2], 
                      p = t[5,4])
var5_sp <- data.frame(Variable = "Advice Family",
                      Coefficient = t[6,1],
                      SE = t[6,2], 
                      p = t[6,4])
var6_sp <- data.frame(Variable = "Advice Friends",
                      Coefficient = t[7,1],
                      SE = t[7,2], 
                      p = t[7,4])
var7_sp <- data.frame(Variable = "Advice Journalists",
                      Coefficient = t[8,1],
                      SE = t[8,2], 
                      p = t[8,4])
var8_sp <- data.frame(Variable = "Advice Celebrities",
                      Coefficient = t[9,1],
                      SE = t[9,2], 
                      p = t[9,4])
var9_sp <- data.frame(Variable = "Age",
                      Coefficient = t[10,1],
                      SE = t[10,2], 
                      p = t[10,4])
var10_sp <- data.frame(Variable = "Male",
                       Coefficient = t[11,1],
                       SE = t[11,2], 
                       p = t[11,4])
var11_sp <- data.frame(Variable = "Schooling",
                       Coefficient = t[12,1],
                       SE = t[12,2], 
                       p = t[12,4])
var12_sp <- data.frame(Variable = "City Category by Size",
                       Coefficient = t[13,1],
                       SE = t[13,2], 
                       p = t[13,4])
var13_sp <- data.frame(Variable = "Smoking",
                       Coefficient = t[14,1],
                       SE = t[14,2], 
                       p = t[14,4])
var14_sp <- data.frame(Variable = "Chronic Illness",
                       Coefficient = t[15,1],
                       SE = t[15,2], 
                       p = t[15,4])
var15_sp <- data.frame(Variable = "COVID Previously",
                       Coefficient = t[16,1],
                       SE = t[16,2], 
                       p = t[16,4])
var16_sp <- data.frame(Variable = "Serious COVID Previously",
                       Coefficient = t[17,1],
                       SE = t[17,2], 
                       p = t[17,4])

t=coeftest(md_si, vcov = vcovHC(md_si, "HC1"))

var1_si <- data.frame(Variable = "Advice Doctor",
                      Coefficient = t[2,1],
                      SE = t[2,2], 
                      p = t[2,4])
var2_si <- data.frame(Variable = "Advice Scientist",
                      Coefficient = t[3,1],
                      SE = t[3,2], 
                      p = t[3,4])
var3_si <- data.frame(Variable = "Advice Anti-Vaccine Propagist",
                      Coefficient = t[4,1],
                      SE = t[4,2], 
                      p = t[4,4])
var4_si <- data.frame(Variable = "Advice Politicians",
                      Coefficient = t[5,1],
                      SE = t[5,2], 
                      p = t[5,4])
var5_si <- data.frame(Variable = "Advice Family",
                      Coefficient = t[6,1],
                      SE = t[6,2], 
                      p = t[6,4])
var6_si <- data.frame(Variable = "Advice Friends",
                      Coefficient = t[7,1],
                      SE = t[7,2], 
                      p = t[7,4])
var7_si <- data.frame(Variable = "Advice Journalists",
                      Coefficient = t[8,1],
                      SE = t[8,2], 
                      p = t[8,4])
var8_si <- data.frame(Variable = "Advice Celebrities",
                      Coefficient = t[9,1],
                      SE = t[9,2], 
                      p = t[9,4])
var9_si <- data.frame(Variable = "Age",
                      Coefficient = t[10,1],
                      SE = t[10,2], 
                      p = t[10,4])
var10_si <- data.frame(Variable = "Male",
                       Coefficient = t[11,1],
                       SE = t[11,2], 
                       p = t[11,4])
var11_si <- data.frame(Variable = "Schooling",
                       Coefficient = t[12,1],
                       SE = t[12,2], 
                       p = t[12,4])
var12_si <- data.frame(Variable = "City Category by Size",
                       Coefficient = t[13,1],
                       SE = t[13,2], 
                       p = t[13,4])
var13_si <- data.frame(Variable = "Smoking",
                       Coefficient = t[14,1],
                       SE = t[14,2], 
                       p = t[14,4])
var14_si <- data.frame(Variable = "Chronic Illness",
                       Coefficient = t[15,1],
                       SE = t[15,2], 
                       p = t[15,4])
var15_si <- data.frame(Variable = "COVID Previously",
                       Coefficient = t[16,1],
                       SE = t[16,2], 
                       p = t[16,4])
var16_si <- data.frame(Variable = "Serious COVID Previously",
                       Coefficient = t[17,1],
                       SE = t[17,2], 
                       p = t[17,4])


allModelFrame <- data.frame(rbind(
  var1_pf, var2_pf, var3_pf, var4_pf, var5_pf, var6_pf, var7_pf, var8_pf, var9_pf, var10_pf, var11_pf, var12_pf, var13_pf, var14_pf, var15_pf, var16_pf,
  var1_mo, var2_mo, var3_mo, var4_mo, var5_mo, var6_mo, var7_mo, var8_mo, var9_mo, var10_mo, var11_mo, var12_mo, var13_mo, var14_mo, var15_mo, var16_mo,
  var1_as, var2_as, var3_as, var4_as, var5_as, var6_as, var7_as, var8_as, var9_as, var10_as, var11_as, var12_as, var13_as, var14_as, var15_as, var16_as,
  var1_sp, var2_sp, var3_sp, var4_sp, var5_sp, var6_sp, var7_sp, var8_sp, var9_sp, var10_sp, var11_sp, var12_sp, var13_sp, var14_sp, var15_sp, var16_sp,
  var1_si, var2_si, var3_si, var4_si, var5_si, var6_si, var7_si, var8_si, var9_si, var10_si, var11_si, var12_si, var13_si, var14_si, var15_si, var16_si))


allModelFrame$col = "lightblue"
allModelFrame$col[allModelFrame$p<0.05] = "darkblue"

allModelFrame$fcol <- factor(allModelFrame$col)


# allModelFrame$ypos=c(2.5,6.5,10.5,14.5,18.5,22.5,
#                      3.5,7.5,11.5,15.5,19.5,23.5)

allModelFrame$ypos=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5,
                     80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5,
                     80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5,
                     80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5,
                     80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5)


allModelFrame$shape=c(16,16,16,16,16,16,16,16,16,16,16, 16, 16, 16, 16, 16,
                      16,16,16,16,16,16,16,16,16,16,16, 16, 16, 16, 16, 16,
                      16,16,16,16,16,16,16,16,16,16,16, 16, 16, 16, 16, 16,
                      16,16,16,16,16,16,16,16,16,16,16, 16, 16, 16, 16, 16,
                      16,16,16,16,16,16,16,16,16,16,16, 16, 16, 16, 16, 16)

# Point estimates figure

library(ggplot2)
library(wooldridge)
library(dplyr)
library(tibble)
library(lemon)
library(ggstance)


png("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot1.png", width=800, height=1200)
ggplot(allModelFrame[1:16,]) + 
  geom_vline(xintercept = 0, linetype="solid", size=4, colour="gray80") +

  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[1:16], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[1:16]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[1:16]) + 
  
  ggtitle("Pfizer") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))

dev.off()


png("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot2.png", width=800, height=1200)
ggplot(allModelFrame[17:32,]) + 
  geom_vline(xintercept = 0, linetype="solid", size=4, colour="gray80") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[17:32], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[17:32]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[17:32]) + 
  
  ggtitle("Moderna") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))

dev.off()

png("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot3.png", width=800, height=1200)
ggplot(allModelFrame[33:48,]) + 
  geom_vline(xintercept = 0, linetype="solid", size=4, colour="gray80") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[33:48], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[33:48]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[33:48]) + 
  
  ggtitle("AstraZeneca") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()

png("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot4.png", width=800, height=1200)
ggplot(allModelFrame[49:64,]) + 
  geom_vline(xintercept = 0, linetype="solid", size=4, colour="gray80") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[49:64], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[49:64]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[49:64]) + 
  
  ggtitle("Sputnik") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()

png("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot5.png", width=800, height=1200)
ggplot(allModelFrame[65:80,]) + 
  geom_vline(xintercept = 0, linetype="solid", size=4, colour="gray80") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[65:80], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[65:80]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[65:80]) + 
  
  ggtitle("Sinopharm") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))

dev.off()


# Save plots in EPS

cairo_ps("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot1.eps", width=20, height=30)
ggplot(allModelFrame[1:16,]) + 
  geom_vline(xintercept = 0, linetype="dashed", size=4, colour="gray50") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[1:16], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[1:16]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[1:16]) + 
  
  ggtitle("Pfizer") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()


cairo_ps("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot2.eps", width=20, height=30)
ggplot(allModelFrame[17:32,]) + 
  geom_vline(xintercept = 0, linetype="dashed", size=4, colour="gray50") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[17:32], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[17:32]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[17:32]) + 
  
  ggtitle("Moderna") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()

cairo_ps("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot3.eps", width=20, height=30)
ggplot(allModelFrame[33:48,]) + 
  geom_vline(xintercept = 0, linetype="dashed", size=4, colour="gray50") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[33:48], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[33:48]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[33:48]) + 
  
  ggtitle("AstraZeneca") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()

cairo_ps("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot4.eps", width=20, height=30)
ggplot(allModelFrame[49:64,]) + 
  geom_vline(xintercept = 0, linetype="dashed", size=4, colour="gray50") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[49:64], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[49:64]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[49:64]) + 
  
  ggtitle("Sputnik") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()

cairo_ps("c:/Users/bleng/Documents/1 munka/0 járvány/vaccine/preferences/coeff_plot5.eps", width=20, height=30)
ggplot(allModelFrame[65:80,]) + 
  geom_vline(xintercept = 0, linetype="dashed", size=4, colour="gray50") +
  geom_point(aes(x = Coefficient, y = ypos), shape = allModelFrame$shape[65:80], 
             show.legend = FALSE, size=12, colour=allModelFrame$fcol[65:80]) +
  geom_linerangeh(aes(xmin = Coefficient - 1.96*SE, y = ypos, 
                      xmax = Coefficient + 1.96*SE), show.legend = FALSE,
                  lwd = 3,  colour=allModelFrame$fcol[65:80]) + 
  
  ggtitle("Sinopharm") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=50)) + 
  theme(axis.text.x =element_text(size=20)) +
  theme(axis.text.y =element_text(size=30)) +
  theme(axis.title.x =element_text(size=20)) +
  labs(x="Coefficients of Rejection", y="") +
  scale_y_continuous(labels=c("Advice Doctor",
                              "Advice Scientist",
                              "Advice Anti-Vaccine Propagators",
                              "Advice Politicians",
                              "Advice Family",
                              "Advice Friends",
                              "Advice Journalists",
                              "Advice Celebrities",
                              "Age",
                              "Female",
                              "University",
                              "City Category by Size",
                              "Smoking",
                              "Chronic Illness",
                              "Covid Previously",
                              "Serious Covid Previously"), 
                     breaks=c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5))
dev.off()

# Spring code (with parameters determined for spring season)


# Annotated code to determine the 
# cost of byssal thread production from 
# the relationship between thread production and growth.

# The linear model is set up as follows:
# B1 = -thread_num
# B2 = a'*mass^d*time
# B3 = -b*mass^e*time
# y = growth in Joules
# y = B1 * h + B2 * food + B3

# Set up workspace
rm(list = ls())
library(multcompView)
library(RColorBrewer)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(lsmeans)
library(car)
library(multcomp)

# Define functions
Data_plot <- function(df, thread_count, growth_WT_len_cubed_g, treatment, new.pal) {
  gg <- ggplot(df, aes(x=thread_count, 
                       y=growth_WT_len_cubed_g*1000, 
                       color = treatment, 
                       fill = treatment
  )) + 
    geom_point() +
    geom_smooth(method = "lm") +
    scale_x_continuous(limits = c(0,600), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-5,105), expand = c(0, 0)) +
    xlab("Thread production (#/day)") +
    ylab("Tissue growth (mg DW)") 
  gg <- gg +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      axis.line = element_line(color = "black", 
                               size = 0.5, linetype = "solid")
    )+ scale_color_manual(values=new.pal) + scale_fill_manual(values=new.pal)
  gg
}

mypalette<-brewer.pal(8,"YlGnBu")
mypalette2<-brewer.pal(8,"Greys")
new.pal <- c(mypalette[c(4,6)], mypalette2[8])
palette(c(mypalette[c(4,6)], mypalette2[8]))

# Parameter values ####
# Parameters are as described here as in the paper. Means plus or minus SE are marked as high and low values.

# conversion factor ####
gWW_per_gDW <- 3.98 # From Summer 2015 collection, shape coeff estimation


# exponents ####
d <- 0.69 # intake exponent
d_high <- 0.70 # mean + SE
d_low <- 0.68 # mean - SE

e <- 1


# constants
b <- 0.158 # J/(day*mgDW) #Now in dry weight, b differs in spring and autumn
b_low <- 0.158-0.043 # J/(day*mgDW) #Now in dry weight
b_high <- 0.158+0.043 # J/(day*mgDW) #Now in dry weight

ED <- 21.6 # J/mg soft tissue
ED_low <- 21.6 - 1.6
ED_high <- 21.6 + 1.6

Wopt_measured_gDW <- .72 #gDW 
Wopt_measured_gDW_high <- .72+.06
Wopt_measured_gDW_low <- .72-.06  

Wopt_measured_mgDW <- Wopt_measured_gDW * 1000
a_fromWopt <- (b*e)/(( Wopt_measured_mgDW)^(d-e)*d) # backwards calculation of a
(a_prime <- a_fromWopt) 


# import data ####
season="Spring"
setwd("~/BE/BE_2019_06_26/Datasets")
# df_all <- read.csv(file="RawData_AutSpr 2019ks4.csv", stringsAsFactors = FALSE)
df_all <- read.csv(file="Byssus_allocation_experiment_data.csv", stringsAsFactors = FALSE)
df_all$treatment <- as.factor(df_all$treatment)
df_all$season <- as.factor(df_all$season)
df_all$treatment <- ordered(df_all$treatment, levels = c("never", "weekly","daily"))

df <- df_all[df_all$season==season,]
df.never <- df[df$treatment=="never",]
cage <- as.factor(df[df$season==season,]$cage)
thread_num <- df[df$season==season,]$thread_count
time <- as.Date(df$date_final, "%d-%b-%y")-as.Date(df$date_init,"%d-%b-%y")

mass_mg_dw <- df$initial_dry_wt_calc_buoy_wt_g*1000 #gDW
growth_tissue <- df$growth_WT_len_cubed_g*1000 #mgDW
growth_tissue_J <- growth_tissue*(ED) #Converted Joules (J/mgDW =ED) J per month
# qqp(growth_tissue_J)
# qqp(thread_num)

# Plot
par(mfrow = c(1,1))
Data_plot(df, thread_count, growth_WT_len_cubed_g, treatment, new.pal)

# Step 1 of model ####
y <- lm(growth_tissue_J~thread_num)
(new <- summary(y))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 847.3827   118.7129   7.138 1.05e-08 ***
#   thread_num   -1.0183     0.3781  -2.693   0.0102 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 443.7 on 41 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.1503,	Adjusted R-squared:  0.1296 
# F-statistic: 7.252 on 1 and 41 DF,  p-value: 0.01021
(cost_per_thread <- -1*new$coefficients[2])
# 1.018307 SE = 0.3781


# Build model ####
Th <- -thread_num
B1 <- Th
Intake <- a_prime*((mass_mg_dw)^d)*29 #good, this is already in mass_mg_dw... hope this is g not mg
B2 <- Intake
Metab <- -(b)*((mass_mg_dw)^e)*29 #*(1-baseline_byssus_multiplier) # b is now the cost per mgWW
B3 <- Metab
reproduction <- 0

dat_lm <- data.frame(
  growth_tissue_J = growth_tissue_J, #mgDW converted to J
  B1 = B1, #byssus
  B2 = B2, #intake
  B3 = B3, #cost
  cage = cage
)
str(dat_lm)
dat_lm <- dat_lm[!is.na(dat_lm$growth_tissue_J),]
dat_lm <- dat_lm[!is.na(dat_lm$B1),]

# 2nd of two steps: ####
# y <- lm(growth_tissue_J ~ B2 + offset(B3+01.13*B1) + 0, data = dat_lm)
y <- lm(growth_tissue_J ~ B2 + offset(B3+cost_per_thread*B1) + 0, data = dat_lm)
(sum_y<-summary(y))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# B2  1.42957    0.09968   14.34   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 516.7 on 42 degrees of freedom
# Multiple R-squared:  0.559,	Adjusted R-squared:  0.5485 
# F-statistic: 53.24 on 1 and 42 DF,  p-value: 5.494e-09
f_intake <- sum_y$coefficients[1]
obs_growth <-growth_tissue #mg_DW
pred_growth_J <- (cost_per_thread*B1+f_intake*B2+B3) #J converted into mgDW
pred_growth_mg_DW <- pred_growth_J / ED #Converts from J to mgDW
pred_growth_2 <- pred_growth_mg_DW
mass_dw <- (mass_mg_dw)
plot(f_intake*B2~mass_mg_dw, 
     ylim = c(0,1500), 
     xlim = c(0,260), 
     col = "blue", type = "b",
     ylab = "J over 29 days") #Intake
#points(-B3~mass_mg_dw, type = "b") # Non-byssus metabolic costs
points((f_intake*B2+B3)~mass_mg_dw, col = "light blue", type = "b") #Intake minus non-byssus metabolic costs
points((f_intake*B2+B3+B1)~mass_mg_dw, col = "purple", type = "p") #Intake minus all costs
points((growth_tissue_J)~mass_mg_dw, col = "black", pch = 20, type = "p") #Intake minus all costs

# y <- lm(growth_tissue_J ~ B2 + B3 + offset(cost_per_thread*B1) + 0, data = dat_lm)
# (sum_y<-summary(y))
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # B2   4.5155     0.7779   5.805 8.20e-07 ***
# #   B3   8.1590     1.7937   4.549 4.73e-05 ***
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 443.8 on 41 degrees of freedom
# # Multiple R-squared:  0.6651,	Adjusted R-squared:  0.6488 
# # F-statistic: 40.72 on 2 and 41 DF,  p-value: 1.818e-10
# 
# This demonstrates that metabolic cost is likely 8 times 
# what we think it is given respiration from the lab experiment. 
# f_intake <- sum_y$coefficients[2]
# maint_mult <- sum_y$coefficients[1]
# obs_growth <-growth_tissue #mg_DW
# pred_growth_J <- (cost_per_thread*B1+f_intake*B2+maint_mult * B3) #J converted into mgDW
# pred_growth_mg_DW <- pred_growth_J / ED #Converts from J to mgDW
# pred_growth_2 <- pred_growth_mg_DW
# mass_dw <- (mass_mg_dw)
# plot(f_intake*B2~mass_mg_dw, 
#      ylim = c(0,9000), 
#      xlim = c(0,260), 
#      col = "blue", type = "b",
#      ylab = "J over 29 days") #Intake
# #points(-B3~mass_mg_dw, type = "b") # Non-byssus metabolic costs
# points((f_intake*B2 + maint_mult * B3)~mass_mg_dw, col = "light blue", type = "b") #Intake minus non-byssus metabolic costs
# points((f_intake*B2 + maint_mult * B3 + B1)~mass_mg_dw, col = "purple", type = "p") #Intake minus all costs
# points((growth_tissue_J)~mass_mg_dw, col = "black", pch = 20, type = "p") #Intake minus all costs





#Plot observed growth and residuals as a function of thread production ####
par(mfrow = c(3,2), mar = c(5,4, 0,0) + 0.1)
plot(obs_growth~thread_num, col=df$treatment, pch = 20, 
     xlab = "Thread production (#/mussel)",
     ylab = "Observed growth (mg DW)",
     xlim = c(0,600),
     ylim = c(0,90))
m <- lm(obs_growth~thread_num)
abline(m)
summary(m)
# m <- lm(obs_growth~thread_num*df$treatment)



plot(pred_growth_2~thread_num, col=df$treatment, pch = 20, 
     xlab = "Thread production (#/mussel)",
     ylab = "Model growth prediction (mg DW)",
     xlim = c(0,600),
     ylim = c(0,90))
m <- lm(pred_growth_2~thread_num)
abline(m)
summary(m)

# x <- seq(min(thread_num),max(thread_num), by = 1)
# y1 <- m$coefficients[1] + x * m$coefficients[2]
# lines(y1~x)
plot(growth_tissue~mass_dw, col=df$treatment, pch = 20,
     xlab = "Initial tissue mass (mg DW)",
     ylab = "Observed growth (mg DW)",
     xlim = c(50,250),
     ylim = c(0,90))
m <- lm(growth_tissue~mass_dw)
abline(m, lty = 2)
summary(m)


plot(pred_growth_2~mass_dw, col=df$treatment, pch = 20,
     xlab = "Tissue mass (mg DW)",
     ylab = "Model growth prediction (mg DW)",
     xlim = c(50,250),
     ylim = c(0,90))
m <- lm(pred_growth_2~mass_dw)
abline(m, lty = 2)
summary(m)

plot((pred_growth_2-obs_growth)~mass_dw, #In mg DW
     #ylim = c(-0,.5),
     #xlim = c(-0,.5),
     ylab = "Residuals (mg DW)",
     xlab = "Initial tissue mass (mg DW)",
     col=df$treatment,
     pch = 20
)
x <- seq(from = 0, to=400, by=.1)
y0 <- rep(0, len = length(x))
lines(x,y0, lty = 2)


plot((pred_growth_2-obs_growth)~obs_growth, #In mg DW
     ylim = c(-80,40),
     xlim = c(0,105),
     ylab = "Residuals (mg DW)",
     xlab = "Observed growth (mg DW)",
     col=df$treatment,
     pch = 20
)
x <- seq(from = 0, to=400, by=.1)
y0 <- rep(0, len = length(x))
lines(x,y0, lty = 2)
#End plot ####


# Energy budget calculations ####

# Daily energy budget (J/day)
(intake <- sum_y$coefficients[1]*a_prime*((mass_mg_dw)^d)) #good, this is already in mass_mg_dw... hope this is g not mg
(intake <- (sum_y$coefficients[1]*B2)/29)
(cost <- (b)*((mass_mg_dw)^e)) #*(1-baseline_byssus_multiplier)
(cost <- -B3/29)
#cost <- sum_y$coefficients[3]*(b)*((mass_mg_dw)^e) #*(1-baseline_byssus_multiplier)
# 

# byssus_baseline <- b*(mass_mg_dw^e)*baseline_byssus_multiplier
byssus_induced <- (thread_num*cost_per_thread)/29
reproduction <- 0
(model.predG_J_daily <- intake-cost-byssus_induced-reproduction) #predicts mussel growth in J
pred_growth_J /29
#model.predG_mg_DW <- model.predG_J/ED  #converts to gDW

gonad_proportion <- df$final_gonad_wt_dry_g/df$final_total_wt_dry_g

proportions.df <- data.frame(
  treat=df$treat,
  thread_num = thread_num/29,
  mass = (mass_mg_dw),
  cost_gonad=cost*gonad_proportion, # Cost of gonad
  cost_somatic = cost*(1-gonad_proportion), # Cost of somatic tissue
  byssus = byssus_induced, # Cost of byssus
  growth = model.predG_J_daily, # Surplus to growth
  intake
)


setwd("~/BE/BE_2019_06_26/Datasets")
write.csv(file = "Spring_proportions.df.csv", proportions.df)


library(dplyr)
library(tidyr)
library(ggplot2)
library(sjPlot)
library(sjlabelled)
library(sjmisc)




proportions.df <- proportions.df[proportions.df$thread_num != 102,]
expand(proportions.df,treat)
new.df <- gather(proportions.df,prop,measurement,cost_gonad:intake)
y <- aggregate(measurement ~ treat + prop, data = new.df, mean)

#setwd("~/BE/BE_2019_06_26/stats")
mod1<-lm(measurement~thread_num*prop, data = new.df)
summary(mod1)
Anova(mod1)
theme_set(theme_sjplot())
# plot_model(mod1)
# write.csv(file = "stats_temp.csv",mod1$coefficients)
# write.csv(file = "stats_temp_ANOVA.csv",Anova(mod1))


setwd("~/BE/BE_2019_06_26/figs_201908")
#pdf(file = "Daily_Energy_Budget_Spring_2_step.pdf", width = 5.7, height = 6.2)
gg <- ggplot(new.df, aes(x=thread_num,
                         y=measurement,
                         color = prop,
                         fill = prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0,20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,70), expand = c(0, 0)) +
  xlab("Rate of thread production (threads per day)") +
  # scale_color_manual(values=c('#FD6000','#000000', '#000000','#FFB601', '#00CC3B' )) +
  # scale_fill_manual(values=c('#FD6000','#000000', '#000000','#FFB601', '#00CC3B' )) +
  ylab("Rate (J / day)")
gg <- gg +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.line = element_line(color = "black",
                             size = 0.5, linetype = "solid")
  )
gg
#dev.off()

proportions.df$intake_prop <- proportions.df$byssus/(proportions.df$intake)

(a <- aggregate(intake_prop ~ treat, data = proportions.df, FUN =mean))
b <- aggregate(intake_prop ~ treat, data = proportions.df, FUN =sd)
c <- aggregate(intake_prop ~ treat, data = proportions.df, FUN = length)
(SE <- b[,2]/sqrt(c[,2]-1))

proportions.df$cost_prop <- proportions.df$byssus/(proportions.df$cost_somatic+proportions.df$byssus)

(a <- aggregate(cost_prop ~ treat, data = proportions.df, FUN =mean))
b <- aggregate(cost_prop ~ treat, data = proportions.df, FUN =sd)
c <- aggregate(cost_prop ~ treat, data = proportions.df, FUN = length)
(SE <- b[,2]/sqrt(c[,2]-1))


setwd("~/BE/BE_2019_06_26/Datasets")
write.csv(file = "Spring_proportions.df.csv", proportions.df)

# Fit curve to proportion of somatic cost ####
output <- nls(cost_prop~Vmax*(1-exp(-thread_num/tau)), start=list(tau=10,Vmax=100), data = proportions.df)
AIC(output)
output_MM <- nls(cost_prop~(Vmax*thread_num)/(Km+thread_num), start=list(Km=10,Vmax=80), data = proportions.df)
AIC(output_MM)
summary(output)

head(proportions.df)
gg <- ggplot(proportions.df, aes(x=thread_num, 
                                 y=cost_prop*100
                                 # color = treat, 
                                 # fill = treat
)) + 
  geom_point(aes(color = treat), shape = 16) +
  geom_smooth(method="nls", 
              formula=y~Vmax*(1-exp(-x/tau)), # this is an nls argument
              method.args = list(start=c(tau=10,Vmax=100)), # this too
              se=FALSE, color = "gray")+
  scale_x_continuous(limits = c(0,20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0,110), expand = c(0, 0)) +
  xlab("Thread production (#/day)") +
  ylab("Percent of somatic cost (%)") 
gg <- gg +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.line = element_line(color = "black", 
                             size = 0.5, linetype = "solid")
  ) + scale_color_manual(values=new.pal) + scale_fill_manual(values=new.pal)
gg

# Save image ####
#setwd("~/BE/BE_2019_06_26/figs_201908")
#pdf(file = "Spring_percent.pdf", width = 5, height = 4.5)
#gg
#dev.off()


#-------- Nesting check using lme
# I can nest in this way, which is pretty cool. 
# require(lme4)
# require(nlme)
# require(lmerTest)

# y1 <- lme(growth_tissue_J ~ B1 + B2 + offset(B3) + 0, random =~ 1|cage, data = dat_lm)
# #y2 <- lme(growth_tissue_J ~ B1 + B2 + offset(B3) + 0, random =~ B2|cage, data = dat_lm)
# y3 <- gls(growth_tissue_J ~ B1 + B2 + offset(B3) + 0, method = "REML", data = dat_lm)
# anova(y1, y3) # The model that does not include the cage is better. 
# # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# # y1     1  4 -53.66790 -46.81361 30.83395                        
# # y3     2  3 -50.18675 -45.04604 28.09338 1 vs 2 5.481149  0.0192

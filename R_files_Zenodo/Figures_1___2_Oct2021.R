
#############################################################################################
#Plots for African wild dog feeding priority paper
#Figures 1 & 2

#install packages
install.packages("Rtools")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("ggforce")
install.packages("wesanderson")
install.packages("colorRamps")
install.packages("viridis")
install.packages("cowplot")
install.packages("ggsci")
install.packages("LaplacesDemon")

#load packages
library(ggplot2)
library(RColorBrewer)
library(ggforce)
library(wesanderson)
library(colorRamps) 
library(viridis)
library(cowplot)
library(ggsci)
library(LaplacesDemon)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Set theme
theme_set(theme_minimal())

col <- "black"
colBG <- "white"
font <- "Helvetica"
fontSize <- 20

pointsize <- 3 # define size of points in the plot
lineSize <- 3 # define line size
alphaRibbon <- 0.3


################################
# Figure 1: Model fit from model averaged GLMM outputs 


df.estimates.fig1 <- data.frame(term = c("intercept","CarcassCondition","RelativePFQ","Noon_kill", "CarcassCondition:Noonkill_Corrected", "RelativePFQ:CarcassCondition"),
                                estimate = c(3.4381764,-0.0120271, 0.4926998, -0.1224332, 0.0009608, -0.0002976),
                                stringsAsFactors = FALSE)
#3.119518
#-0.009847
#0.52847
#-0.000507


df.fig1 <- as.data.frame(expand.grid(SuccessLikelihood.mean = as.numeric(NA), #mean estimate of success likelihood
                                     RelativePFQ = c(-6, -4, -2, 0, 2, 4),
                                     CarcassCondition = rev(seq(0,100,0.1)),
                                     Noon_Kill = as.numeric(3), #set as median
                                     interaction1 = as.numeric(NA),
                                     interaction2 = as.numeric(NA)
                            
))

# we calculate the interaction
df.fig1$interaction <- df.fig1$RelativePFQ*df.fig1$CarcassCondition

# we create an additional column and turn RelativePFQ into a factor
df.fig1$RelativePFQ.fac <- factor(x = df.fig1$RelativePFQ, labels = c("-6","-4", "-2", "0", "2","4"))

head(df.fig1)

# now, let's estimate the success likelihoods given the values of covariates/interaction in df.fig1 
# and the coefficient estimates in df.estimates.fig1:


# we calculate by row:
for(i in 1:nrow(df.fig1)){
  #mean estimate of success likelihood --> we use estimates (i.e. coefficient means) given in df.estimates.fig1
  p <- df.estimates.fig1$estimate[which(df.estimates.fig1$term=="intercept")] + 
    df.estimates.fig1$estimate[which(df.estimates.fig1$term=="CarcassCondition")]*df.fig1$CarcassCondition[i] + 
    df.estimates.fig1$estimate[which(df.estimates.fig1$term=="RelativePFQ")]*df.fig1$RelativePFQ[i] + 
    df.estimates.fig1$estimate[which(df.estimates.fig1$term=="Noon_kill")]*3 + 
    df.estimates.fig1$estimate[which(df.estimates.fig1$term=="RelativePFQ:CarcassCondition")]*df.fig1$interaction[i]
  
  # response is binomial, so we have to back transform using the inverse logit link function:
  df.fig1$SuccessLikelihood.mean[i] <- invlogit(p)
}  




#we create the plot (Fig 1):
fig1 <- ggplot(data = df.fig1, aes(x = CarcassCondition))+
  
  # we add line of mean estimates:
  geom_line(aes(y = SuccessLikelihood.mean, color = RelativePFQ.fac), size = 1)+
  
  
  # we customize axes:
  scale_y_continuous(name = "Likelihood of success\\n(approacher joins carcass)", limits = c(0,1), breaks = seq(0,1,0.1))+
  scale_x_reverse(name = "Carcass condition (% remaining)", limits = c(100,0), breaks = seq(0,100,10))+
  
  # We use colorblind-friendly colors from the 'ggsci' package for fill colors:
  scale_fill_jco()+
  # ... and for line colors:
  scale_color_jco()+
  
  # customize colors, font size, font, positions, etc. of axes, text, legends, ...
  theme(axis.text = element_text(colour = col, size = fontSize, family = font),
        axis.title.x = element_text(colour = col, size = fontSize, family = font, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = col, size = fontSize, family = font, margin = margin(0,15,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = col),
        axis.ticks = element_line(colour = col),
        legend.position = c(0.65, 0.1), 
        legend.justification = c(0,0), 
        legend.box.margin=margin(c(0,20,0,20)),
        legend.text = element_text(colour = col, size = fontSize, family = font),
        legend.title = element_blank(),
        legend.key.width = unit(35,"pt"),
        panel.border = element_blank(),
        plot.title = element_text(colour = col, size = fontSize, family = font, hjust = -0.068, margin = margin(0,0,10,0)),
        panel.background = element_blank(),
        plot.background = element_blank(),
  ) 



summary(df.fig1)
df.fig1$RelativePFQ<-as.factor(df.fig1$RelativePFQ.fac)


# create plot:
fig1R<-ggplot((data = df.fig1), aes(CarcassCondition, SuccessLikelihood.mean, color = RelativePFQ)) +
  geom_line(aes(size = RelativePFQ)) +
  scale_size_manual(values = c("-6" = 4, "-4" = 2, "-2" = 1,  "0" = 1, "2" = 1, "4" =2)) +
  scale_colour_manual(values = c("-6" = "red", "-4" = "red", "-2" = "red",  "0" = "black", "2" = "blue", "4" ="blue"))+
  
  scale_y_continuous(name = "Likelihood of success\\n(approacher joins carcass)", limits = c(0,1), breaks = seq(0,1,0.1))+
  scale_x_reverse(name = "Carcass condition (% remaining)", limits = c(100,0), breaks = seq(0,100,10))+
  # customize colors, font size, font, positions, etc. of axes, text, legends, ...
  theme(axis.text = element_text(colour = col, size = fontSize, family = font),
        axis.title.x = element_text(colour = col, size = fontSize, family = font, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = col, size = fontSize, family = font, margin = margin(0,15,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = col),
        axis.ticks = element_line(colour = col),
        legend.position = c(0.65, 0.05), 
        legend.justification = c(0,0), 
        legend.box.margin=margin(c(0,20,0,20)),
        legend.text = element_text(colour = col, size = fontSize, family = font),
        legend.title = element_text(colour = col, size = fontSize, family = font),
        legend.key.width = unit(35,"pt"),
        panel.border = element_blank(),
        plot.title = element_text(colour = col, size = fontSize, family = font, hjust = -0.068, margin = margin(0,0,10,0)),
        panel.background = element_blank(),
        plot.background = element_blank(),
  ) 

# view plot:
print(fig1R)



# export figure as PDF:
pdf(file = "Fig1R4.pdf", width = 10, height = 7, family = "Helvetica", fonts = "Helvetica")
fig1R
dev.off()



################################
# Figure 2: Proportion of approaches to carcasses:

df.fig2 <- data.frame(PropApproaches.mean = c(0.078544061, 	0.352380952,	0.094562648, 	0.03358209, 0.467741935,	0.135135135, 0.036842105, 0.571428571,	0.292857143), 
                      ConsumptionStage = c(rep(1,3),rep(2,3),rep(3,3)),
                      RelativePFQ = rep(c(-1:1),3))

df.fig2
# we create an additional column and turn RelativePFQ into a factor
df.fig2$RelativePFQ.fac <- factor(x = df.fig2$RelativePFQ, labels = c("-1 (approacher further down queue)", "0 (approacher same position in queue as incumbent)", "+1 (approacher higher in queue)"))

# the same for consumption stage
df.fig2$ConsumptionStage.fac <- factor(x = df.fig2$ConsumptionStage, labels = c("Late","Mid","Early"))


# we create the plot (Fig 2):
df.fig2$ConsumptionStage.fac <- factor(df.fig2$ConsumptionStage.fac, levels = c("Early", "Mid", "Late"))
fig2 <- ggplot(data = df.fig2, aes(y = PropApproaches.mean, x = ConsumptionStage.fac, fill = RelativePFQ.fac))+
  geom_bar(stat="identity", position=position_dodge(), color = "white") +
  
  # we customize axes:
  scale_y_continuous(name = "Proportion of all approaches to carcass", limits = c(0,0.7), breaks = seq(0,1,0.1))+
  xlab("Carcass Consumption Stage")+
  
  # We use colorblind-friendly colors from the 'ggsci' package for fill colors:
  scale_fill_jco()+
  
  # customize colors, font size, font, positions, etc. of axes, text, legends, ...
  theme(axis.text = element_text(colour = col, size = fontSize, family = font),
        axis.title.x = element_text(colour = col, size = fontSize, family = font, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = col, size = fontSize, family = font, margin = margin(0,15,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = col),
        axis.ticks = element_line(colour = col),
        legend.position = c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.margin=margin(c(0,20,0,20)),
        legend.text = element_text(colour = col, size = fontSize, family = font),
        legend.title = element_blank(),
        legend.key.width = unit(35,"pt"),
        plot.title = element_text(colour = col, size = fontSize, family = font, hjust = -0.068, margin = margin(0,0,10,0)),
        panel.border = element_blank(),
        plot.background = element_blank(),
  ) 

# view plot:
print(fig2)

# export figure as PDF:
pdf(file = "Fig2.pdf", width = 15, height = 10, family = "Helvetica", fonts = "Helvetica")
fig2
dev.off()

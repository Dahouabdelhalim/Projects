#Ryan Stephens; Finalized December 24, 2020 
rm(list=ls()) # clears workspace

setwd("~/Isotopic_Routing_Small_Mammals/TDF_Estimates_Comparisons")


#################################################################################################
#Get data ready
#################################################################################################
TDF<- read.csv("Bartlett_estimated_TDF.csv",header=T)
head(TDF)

library(dplyr)
library(forcats)
#re-order and re-name levels
TDF<-mutate(TDF, Method = fct_relevel(Method,"SIDER", "SIDER_update", "Meta_analysis", "Field"))#change level order
TDF<-mutate(TDF, Method = fct_recode(Method,"Meta-analysis" = "Meta_analysis", "SIDER+Meta" = "SIDER_update"))#rename levels
TDF<-mutate(TDF, Species = fct_relevel(Species,"M. gapperi", "N. insignis", "P. maniculatus","B. brevicauda"))#rename levels
TDF<-mutate(TDF, Diet =ifelse(Species == "M. gapperi", "Herbivore", #text for diet class
                              ifelse(Species == "B. brevicauda", "Carnivore", "Omnivore")))

head(TDF)
TDF<-TDF%>%mutate(Line=ifelse(Isotope == "d13C", 1, 3))#commonly used TDF values

#split d13C and d15N data sets to make it easier to graph
d15N_TDF<-filter(TDF, Isotope == "d15N")
d13CN_TDF<-filter(TDF, Isotope == "d13C")
#################################################################################################




#################################################################################################
#Graph TDF-d13C
#################################################################################################
library(ggplot2)
library(lemon)
library(cowplot)
d13C_TDF_plot<-ggplot(d13CN_TDF, aes(x=Method, y=TDF_Mean,fill=Method))+ 
#points and bars
  geom_errorbar(aes(ymin=TDF_Mean-TDF_SD, ymax=TDF_Mean+TDF_SD),color="black",width=0)+#error bars
  geom_point(size=3,shape=21)+#points
  scale_fill_manual(values =c("mediumorchid3","darkmagenta", "chartreuse4","orange"))+#custom colors
#facets and background
  coord_capped_cart(bottom='both', left = 'both')+#use package "lemon" to keep tick marks
  facet_grid(. ~ Species, switch = "y")+#facets
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+
  theme_cowplot(12)+#make background transparent for mammal images
#line and text
  geom_hline(aes(yintercept=Line), color='grey50',linetype = "dotted")+#line for commonly used TDFs
  geom_text(aes(x = 2.5, y = 7.8,label = Diet),color="black", size=2.8)+#text for diet class
#strip text
  theme(strip.background = element_blank())+#remove background
  theme(strip.placement = "outside")+#put label outside of y axis 
  theme(strip.text.y = element_text(size = 10))+#size of strip text
  theme(strip.text.x = element_text(size = 8, face="italic", color = "gray65"))+#size of strip text
  theme(panel.spacing.x=unit(0.75, "lines"))+#spacing to match mixing space
#legend
  theme(legend.position = "none")+
#Axes portion
  scale_y_continuous(expand = c(0, 0),limits = c(-.5, 11),breaks=c(0,2,4,6,8))+#custom range
  ylab(expression("TDF " * delta^13 * "C"))+
  theme(axis.title.x= element_blank(),#remove x axis title
        axis.text.x  = element_blank())+#remove x axis text
  theme(axis.title.y = element_text(colour="black", size=10),#adjusts size and color of y axis title
        axis.text.y  = element_text(angle=0, vjust=.5, size=8, colour="black"))#adjusts angles, size, etc of labels for y axis

d13C_TDF_plot#print plot
#################################################################################################




#################################################################################################
#Graph TDF-d15N
#################################################################################################
d15N_TDF_plot<- ggplot(d15N_TDF, aes(x=Method, y=TDF_Mean,fill=Method))+ 
#points and bars
  geom_errorbar(aes(ymin=TDF_Mean-TDF_SD, ymax=TDF_Mean+TDF_SD),color="black",width=0)+#error bars
  geom_point(size=3,shape=21)+#points
  scale_fill_manual(values =c("mediumorchid3","darkmagenta", "chartreuse4","orange"))+#custom colors
#facets and background
  coord_capped_cart(bottom='both', left = 'both')+#use package "lemon" to keep tick marks
  facet_grid(.~Species,switch = "y")+#,labeller = label_parsed
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+
  theme_cowplot(12)+#make background transparent (not needed for this plot but keeps axis lines consistent with d13C plot)
#line
  geom_hline(aes(yintercept=Line), color='grey50',linetype = "dotted")+
#strip text
  theme(strip.background = element_blank())+#remove background
  theme(strip.placement = "outside")+#put label outside of y axis 
  theme(strip.text.x = element_blank())+#remove strip text
  theme(panel.spacing.x=unit(0.75, "lines"))+#spacing to match mixing space
#legend
  theme(legend.position = "none")+
#Axes portion
  scale_y_continuous(expand = c(0, 0),limits = c(-.5, 11),breaks=c(0,2,4,6,8))+#custom range
  ylab(expression("TDF " * delta^15 * "N"))+
  theme(axis.title.x= element_blank(),#adjusts size and color of x axis title
        axis.text.x  = element_blank())+
  theme(axis.title.y = element_text(colour="black", size=10),#adjusts size and color of y axis title
        axis.text.y  = element_text(angle=0, vjust=.5, size=8, colour="black"))#adjusts angles, size, etc of labels for y axis

d15N_TDF_plot#print plot
#################################################################################################




#################################################################################################
#Mixing space test plot
#################################################################################################
Mixing_space<- read.csv("Mixing_space_summary.csv",header=T)#Code to produce these values starts on line 190 of this script
head(Mixing_space)

#re-order and re-name levels
Mixing_space<-mutate(Mixing_space, Method = fct_relevel(Method,"SIDER", "SIDER_update", "Meta_analysis", "Field"))#re-order levels
Mixing_space<-mutate(Mixing_space, Method = fct_recode(Method,"Meta-analysis" = "Meta_analysis", "SIDER+Meta" = "SIDER_update"))#re-name levels
Mixing_space<-mutate(Mixing_space, Species = fct_relevel(Species,"MYGA", "NAIN", "PEMA","BLBR"))#re-order levels


Mixing_space_plot<-ggplot(data=Mixing_space, aes(x=Method, y=Proportion, fill = Method, color = Method)) +
#bars
  geom_bar(stat="identity", width = .7, size = .7, alpha = .5)+#bars
  scale_fill_manual(values =c("mediumorchid3","darkmagenta", "chartreuse4","orange"))+#custom color for fill
  scale_color_manual(values =c("mediumorchid3","darkmagenta", "chartreuse4","orange"))+#custom color for line
#facets and background
  coord_capped_cart(bottom='both', left = 'both')+#use package "lemon" to keep tick marks
  facet_grid(.~Species,switch = "y")+
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+
  theme_cowplot(12)+#make background transparent (not needed for this plot but keeps axis lines consistent with d13C plot)
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1.01),breaks=c(0, .2, .4, .6, .8, 1))+
  theme(panel.spacing.x=unit(.75, "lines"))+#spacing
#strip text
  theme(strip.background = element_blank())+#remove background
  theme(strip.text.x = element_blank())+#remove strip text
#legend
  theme(legend.position = "none")+
#Axes portion
  labs(x="Method",y="Samples within 95%\\nof mixing region") +#parse text onto two lines
  theme(axis.title.x= element_text(colour="black", size=10, hjust=.5),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=50,vjust=1,  hjust = 1,size=8,colour="black"))+
  theme(axis.title.y = element_text(colour="black", size=10),#adjusts size and color of y axis title
        axis.text.y  = element_text(angle=0, vjust=.5, size=8, colour="black"))#adjusts angles, size, etc of labels for y axis

Mixing_space_plot#print plot
#################################################################################################




#################################################################################################
#Combine plots
#################################################################################################
library(cowplot)
library(magick)

ggdraw() +
#mammal illustrations (available from Ryan Stephens upon request if used for non-profit endeavors)
  draw_image("~/Isotopic_Routing_Small_Mammals/TDF_Estimates_Comparisons/Southern_Redbacked_Vole_White.jpg",  x = -.28, y = 0.425, scale = .26) +
  draw_image("~/Isotopic_Routing_Small_Mammals/TDF_Estimates_Comparisons/Woodland_Jumping_Mouse_White.jpg",  x = -.04, y = 0.425, scale = .25) +
  draw_image("~/Isotopic_Routing_Small_Mammals/TDF_Estimates_Comparisons/Northern_short-tailed_shrew_white.jpg",  x = .40, y = 0.425, scale = .27) +
  draw_image("~/Isotopic_Routing_Small_Mammals/TDF_Estimates_Comparisons/Woodland_Deer_Mouse_White.jpg",  x = .18, y = 0.425, scale = .24) +
#plots
  draw_plot(Mixing_space_plot, x = 0.000, y = 0.00, width = 1.00, height = .38)+
  draw_plot(d15N_TDF_plot    , x = 0.035, y = 0.38, width = .965, height = .33)+
  draw_plot(d13C_TDF_plot    , x = 0.035, y = 0.64, width = .965, height = .36)+
#letters
  draw_label("(a)",x = .045, y = .900, size = 10, hjust = 0, color="black", fontface = "bold")+#TDF d13C
  draw_label("(b)",x = .045, y = .635, size = 10, hjust = 0, color="black", fontface = "bold")+#TDF d15N  
  draw_label("(c)",x = .045, y = .390, size = 10, hjust = 0, color="black", fontface = "bold")#Mixing space test

ggsave("Figure_6_TDF_Comparison.tiff",
       plot = last_plot(), width =4.1, height = 6, units = "in",
       dpi = 600) 
#################################################################################################













#################################################################################################
#Mixing space polygons
#################################################################################################
#Ryan Stephens; Finalized December 22, 2020
#Code modified from Mixing Polygon Simulation (2 isotopes); v1.2, 29.6.2015, R version 3.1.2; J.A.Smith et. al. 2013 
#For updates and instructions, visit http://www.famer.unsw.edu.au/downloads.html

rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/TDF_Estimates_Comparisons")
library(dplyr)
library(sp)
library(splancs)
#################################################################################################





#################################################################################################
#Food source information - this is the same for each small mammal species
#################################################################################################
Food_sources <- read.table("Bartlett_diet_items.csv",header=T,sep=",")#load isotope data for food sources
head(Food_sources)

Food_sources_means_sd<-Food_sources %>%
  group_by (Type) %>%
  summarise (d13C_Av = mean(d13C),#average d13C values
             d13C_Sd = sd(d13C),#sd d13C values
             d15N_Av = mean(d15N),#average d15N values
             d15N_Sd = sd(d15N)) %>%#sd d15N values
  as.data.frame()#make into a dataframe
sources<-select(Food_sources_means_sd, -Type)#remove food source name for analysis
#################################################################################################




#################################################################################################
#MYGA - Mixing Polygon Simulations
#################################################################################################

#Hair samples from live trapping- used for all MYGA simulations
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)#Isotope data
MYGA<-filter(Hair, Species=="MYGA")#select species
MYGA_hair<-select(MYGA, d13C, d15N)#select d13C and d15N columns for analysis
MYGA_hair<-na.omit(MYGA_hair)


############################################
#MYGA SIDER
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
MYGA_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER" & Isotope == "d13C")
MYGA_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_SIDER <- as.data.frame(Type)#data frame for the 5 food sources
TDF_MYGA_SIDER$d13C_Av<-MYGA_SIDER_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_MYGA_SIDER$d13C_Sd<-MYGA_SIDER_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_MYGA_SIDER$d15N_Av<-MYGA_SIDER_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_MYGA_SIDER$d15N_Sd<-MYGA_SIDER_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_MYGA_SIDER<-select(TDF_MYGA_SIDER, -Type)#remove names of food sources


#Mixing polygon simulation for MYGA with SIDER TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(MYGA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_MYGA_SIDER),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_MYGA_SIDER[j,1], sd=TDF_MYGA_SIDER[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_MYGA_SIDER[j,3], sd=TDF_MYGA_SIDER[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(MYGA_hair[,1], MYGA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_MYGA_SIDER
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(MYGA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="MYGA_SIDER.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
MYGA_SIDER_data<-as.data.frame(p)
MYGA_SIDER_data.df <- tibble::rownames_to_column(MYGA_SIDER_data, "Iteration")
head(MYGA_SIDER_data.df)
library(tidyr)
MYGA_SIDER_long <- gather(MYGA_SIDER_data.df, Sample, In_out, V1:ncol(MYGA_SIDER_data.df), factor_key=TRUE)#wide to long format
MYGA_SIDER_long$Species<-"MYGA"
MYGA_SIDER_long$Method<-"SIDER"

head(MYGA_SIDER_long)
###########################################


############################################
#MYGA SIDER_update
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
MYGA_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER_update" & Isotope == "d13C")
MYGA_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER_update" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_SIDER_update <- as.data.frame(Type)#data frame for the 5 food sources
TDF_MYGA_SIDER_update$d13C_Av<-MYGA_SIDER_update_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_MYGA_SIDER_update$d13C_Sd<-MYGA_SIDER_update_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_MYGA_SIDER_update$d15N_Av<-MYGA_SIDER_update_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_MYGA_SIDER_update$d15N_Sd<-MYGA_SIDER_update_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_MYGA_SIDER_update<-select(TDF_MYGA_SIDER_update, -Type)#remove names of food sources


#Mixing polygon simulation for MYGA with SIDER_update TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(MYGA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_MYGA_SIDER_update),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_MYGA_SIDER_update[j,1], sd=TDF_MYGA_SIDER_update[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_MYGA_SIDER_update[j,3], sd=TDF_MYGA_SIDER_update[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(MYGA_hair[,1], MYGA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_MYGA_SIDER_update
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(MYGA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="MYGA_SIDER_update.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
MYGA_SIDER_update_data<-as.data.frame(p)
MYGA_SIDER_update_data.df <- tibble::rownames_to_column(MYGA_SIDER_update_data, "Iteration")
head(MYGA_SIDER_update_data.df)
library(tidyr)
MYGA_SIDER_update_long <- gather(MYGA_SIDER_update_data.df, Sample, In_out, V1:ncol(MYGA_SIDER_update_data.df), factor_key=TRUE)#wide to long format
MYGA_SIDER_update_long$Species<-"MYGA"
MYGA_SIDER_update_long$Method<-"SIDER_update"

head(MYGA_SIDER_update_long)
###########################################


############################################
#MYGA Meta_analysis
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
MYGA_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Meta_analysis" & Isotope == "d13C")
MYGA_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Meta_analysis" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_Meta_analysis <- as.data.frame(Type)#data frame for the 5 food sources
TDF_MYGA_Meta_analysis$d13C_Av<-MYGA_Meta_analysis_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_MYGA_Meta_analysis$d13C_Sd<-MYGA_Meta_analysis_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_MYGA_Meta_analysis$d15N_Av<-MYGA_Meta_analysis_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_MYGA_Meta_analysis$d15N_Sd<-MYGA_Meta_analysis_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_MYGA_Meta_analysis<-select(TDF_MYGA_Meta_analysis, -Type)#remove names of food sources


#Mixing polygon simulation for MYGA with Meta_analysis TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(MYGA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_MYGA_Meta_analysis),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_MYGA_Meta_analysis[j,1], sd=TDF_MYGA_Meta_analysis[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_MYGA_Meta_analysis[j,3], sd=TDF_MYGA_Meta_analysis[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(MYGA_hair[,1], MYGA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_MYGA_Meta_analysis
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(MYGA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="MYGA_Meta_analysis.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
MYGA_Meta_analysis_data<-as.data.frame(p)
MYGA_Meta_analysis_data.df <- tibble::rownames_to_column(MYGA_Meta_analysis_data, "Iteration")
head(MYGA_Meta_analysis_data.df)
library(tidyr)
MYGA_Meta_analysis_long <- gather(MYGA_Meta_analysis_data.df, Sample, In_out, V1:ncol(MYGA_Meta_analysis_data.df), factor_key=TRUE)#wide to long format
MYGA_Meta_analysis_long$Species<-"MYGA"
MYGA_Meta_analysis_long$Method<-"Meta_analysis"

head(MYGA_Meta_analysis_long)
###########################################


############################################
#MYGA Field
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
MYGA_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Field" & Isotope == "d13C")
MYGA_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Field" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_Field <- as.data.frame(Type)#data frame for the 5 food sources
TDF_MYGA_Field$d13C_Av<-MYGA_Field_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_MYGA_Field$d13C_Sd<-MYGA_Field_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_MYGA_Field$d15N_Av<-MYGA_Field_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_MYGA_Field$d15N_Sd<-MYGA_Field_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_MYGA_Field<-select(TDF_MYGA_Field, -Type)#remove names of food sources


#Mixing polygon simulation for MYGA with Field TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(MYGA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_MYGA_Field),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_MYGA_Field[j,1], sd=TDF_MYGA_Field[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_MYGA_Field[j,3], sd=TDF_MYGA_Field[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(MYGA_hair[,1], MYGA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_MYGA_Field
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(MYGA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="MYGA_Field.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
MYGA_Field_data<-as.data.frame(p)
MYGA_Field_data.df <- tibble::rownames_to_column(MYGA_Field_data, "Iteration")
head(MYGA_Field_data.df)
library(tidyr)
MYGA_Field_long <- gather(MYGA_Field_data.df, Sample, In_out, V1:ncol(MYGA_Field_data.df), factor_key=TRUE)#wide to long format
MYGA_Field_long$Species<-"MYGA"
MYGA_Field_long$Method<-"Field"

head(MYGA_Field_long)
###########################################

#MYGA model output for all four TDFs
MYGA_model_output<-rbind(MYGA_SIDER_long, MYGA_SIDER_update_long, MYGA_Meta_analysis_long, MYGA_Field_long)
#################################################################################################




#################################################################################################
#NAIN - Mixing Polygon Simulations
#################################################################################################

#Hair samples from live trapping- used for all NAIN simulations
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)#discrimination factors
NAIN<-filter(Hair, Species=="NAIN")#select species
NAIN_hair<-select(NAIN, d13C, d15N)#select d13C and d15N columns for analysis
NAIN_hair<-na.omit(NAIN_hair)


############################################
#NAIN SIDER
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
NAIN_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER" & Isotope == "d13C")
NAIN_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_NAIN_SIDER <- as.data.frame(Type)#data frame for the 5 food sources
TDF_NAIN_SIDER$d13C_Av<-NAIN_SIDER_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_NAIN_SIDER$d13C_Sd<-NAIN_SIDER_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_NAIN_SIDER$d15N_Av<-NAIN_SIDER_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_NAIN_SIDER$d15N_Sd<-NAIN_SIDER_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_NAIN_SIDER<-select(TDF_NAIN_SIDER, -Type)#remove names of food sources


#Mixing polygon simulation for NAIN with SIDER TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(NAIN_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_NAIN_SIDER),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_NAIN_SIDER[j,1], sd=TDF_NAIN_SIDER[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_NAIN_SIDER[j,3], sd=TDF_NAIN_SIDER[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(NAIN_hair[,1], NAIN_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_NAIN_SIDER
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(NAIN_hair, pch=19, cex=1.3)
dev.copy2pdf(file="NAIN_SIDER.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
NAIN_SIDER_data<-as.data.frame(p)
NAIN_SIDER_data.df <- tibble::rownames_to_column(NAIN_SIDER_data, "Iteration")
head(NAIN_SIDER_data.df)
library(tidyr)
NAIN_SIDER_long <- gather(NAIN_SIDER_data.df, Sample, In_out, V1:ncol(NAIN_SIDER_data.df), factor_key=TRUE)#wide to long format
NAIN_SIDER_long$Species<-"NAIN"
NAIN_SIDER_long$Method<-"SIDER"

head(NAIN_SIDER_long)
###########################################


############################################
#NAIN SIDER_update
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
NAIN_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER_update" & Isotope == "d13C")
NAIN_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER_update" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_NAIN_SIDER_update <- as.data.frame(Type)#data frame for the 5 food sources
TDF_NAIN_SIDER_update$d13C_Av<-NAIN_SIDER_update_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_NAIN_SIDER_update$d13C_Sd<-NAIN_SIDER_update_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_NAIN_SIDER_update$d15N_Av<-NAIN_SIDER_update_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_NAIN_SIDER_update$d15N_Sd<-NAIN_SIDER_update_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_NAIN_SIDER_update<-select(TDF_NAIN_SIDER_update, -Type)#remove names of food sources


#Mixing polygon simulation for NAIN with SIDER_update TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(NAIN_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_NAIN_SIDER_update),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_NAIN_SIDER_update[j,1], sd=TDF_NAIN_SIDER_update[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_NAIN_SIDER_update[j,3], sd=TDF_NAIN_SIDER_update[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(NAIN_hair[,1], NAIN_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_NAIN_SIDER_update
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(NAIN_hair, pch=19, cex=1.3)
dev.copy2pdf(file="NAIN_SIDER_update.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
NAIN_SIDER_update_data<-as.data.frame(p)
NAIN_SIDER_update_data.df <- tibble::rownames_to_column(NAIN_SIDER_update_data, "Iteration")
head(NAIN_SIDER_update_data.df)
library(tidyr)
NAIN_SIDER_update_long <- gather(NAIN_SIDER_update_data.df, Sample, In_out, V1:ncol(NAIN_SIDER_update_data.df), factor_key=TRUE)#wide to long format
NAIN_SIDER_update_long$Species<-"NAIN"
NAIN_SIDER_update_long$Method<-"SIDER_update"

head(NAIN_SIDER_update_long)
###########################################


############################################
#NAIN Meta_analysis
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
NAIN_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Meta_analysis" & Isotope == "d13C")
NAIN_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Meta_analysis" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_NAIN_Meta_analysis <- as.data.frame(Type)#data frame for the 5 food sources
TDF_NAIN_Meta_analysis$d13C_Av<-NAIN_Meta_analysis_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_NAIN_Meta_analysis$d13C_Sd<-NAIN_Meta_analysis_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_NAIN_Meta_analysis$d15N_Av<-NAIN_Meta_analysis_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_NAIN_Meta_analysis$d15N_Sd<-NAIN_Meta_analysis_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_NAIN_Meta_analysis<-select(TDF_NAIN_Meta_analysis, -Type)#remove names of food sources


#Mixing polygon simulation for NAIN with Meta_analysis TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(NAIN_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_NAIN_Meta_analysis),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_NAIN_Meta_analysis[j,1], sd=TDF_NAIN_Meta_analysis[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_NAIN_Meta_analysis[j,3], sd=TDF_NAIN_Meta_analysis[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(NAIN_hair[,1], NAIN_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_NAIN_Meta_analysis
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(NAIN_hair, pch=19, cex=1.3)
dev.copy2pdf(file="NAIN_Meta_analysis.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
NAIN_Meta_analysis_data<-as.data.frame(p)
NAIN_Meta_analysis_data.df <- tibble::rownames_to_column(NAIN_Meta_analysis_data, "Iteration")
head(NAIN_Meta_analysis_data.df)
library(tidyr)
NAIN_Meta_analysis_long <- gather(NAIN_Meta_analysis_data.df, Sample, In_out, V1:ncol(NAIN_Meta_analysis_data.df), factor_key=TRUE)#wide to long format
NAIN_Meta_analysis_long$Species<-"NAIN"
NAIN_Meta_analysis_long$Method<-"Meta_analysis"

head(NAIN_Meta_analysis_long)
###########################################


############################################
#NAIN Field
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
NAIN_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Field" & Isotope == "d13C")
NAIN_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Field" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_NAIN_Field <- as.data.frame(Type)#data frame for the 5 food sources
TDF_NAIN_Field$d13C_Av<-NAIN_Field_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_NAIN_Field$d13C_Sd<-NAIN_Field_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_NAIN_Field$d15N_Av<-NAIN_Field_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_NAIN_Field$d15N_Sd<-NAIN_Field_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_NAIN_Field<-select(TDF_NAIN_Field, -Type)#remove names of food sources


#Mixing polygon simulation for NAIN with Field TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(NAIN_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_NAIN_Field),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_NAIN_Field[j,1], sd=TDF_NAIN_Field[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_NAIN_Field[j,3], sd=TDF_NAIN_Field[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(NAIN_hair[,1], NAIN_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_NAIN_Field
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(NAIN_hair, pch=19, cex=1.3)
dev.copy2pdf(file="NAIN_Field.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
NAIN_Field_data<-as.data.frame(p)
NAIN_Field_data.df <- tibble::rownames_to_column(NAIN_Field_data, "Iteration")
head(NAIN_Field_data.df)
library(tidyr)
NAIN_Field_long <- gather(NAIN_Field_data.df, Sample, In_out, V1:ncol(NAIN_Field_data.df), factor_key=TRUE)#wide to long format
NAIN_Field_long$Species<-"NAIN"
NAIN_Field_long$Method<-"Field"

head(NAIN_Field_long)
###########################################

#NAIN model output for all four TDFs
NAIN_model_output<-rbind(NAIN_SIDER_long, NAIN_SIDER_update_long, NAIN_Meta_analysis_long, NAIN_Field_long)
#################################################################################################


#################################################################################################
#PEMA - Mixing Polygon Simulations
#################################################################################################
#Hair samples from live trapping- used for all PEMA simulations
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)#discrimination factors
PEMA<-filter(Hair, Species=="PEMA")#select species
PEMA_hair<-select(PEMA, d13C, d15N)#select d13C and d15N columns for analysis
PEMA_hair<-na.omit(PEMA_hair)


############################################
#PEMA SIDER
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
PEMA_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER" & Isotope == "d13C")
PEMA_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_PEMA_SIDER <- as.data.frame(Type)#data frame for the 5 food sources
TDF_PEMA_SIDER$d13C_Av<-PEMA_SIDER_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_PEMA_SIDER$d13C_Sd<-PEMA_SIDER_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_PEMA_SIDER$d15N_Av<-PEMA_SIDER_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_PEMA_SIDER$d15N_Sd<-PEMA_SIDER_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_PEMA_SIDER<-select(TDF_PEMA_SIDER, -Type)#remove names of food sources


#Mixing polygon simulation for PEMA with SIDER TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(PEMA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_PEMA_SIDER),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_PEMA_SIDER[j,1], sd=TDF_PEMA_SIDER[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_PEMA_SIDER[j,3], sd=TDF_PEMA_SIDER[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(PEMA_hair[,1], PEMA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_PEMA_SIDER
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(PEMA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="PEMA_SIDER.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
PEMA_SIDER_data<-as.data.frame(p)
PEMA_SIDER_data.df <- tibble::rownames_to_column(PEMA_SIDER_data, "Iteration")
head(PEMA_SIDER_data.df)
library(tidyr)
PEMA_SIDER_long <- gather(PEMA_SIDER_data.df, Sample, In_out, V1:ncol(PEMA_SIDER_data.df), factor_key=TRUE)#wide to long format
PEMA_SIDER_long$Species<-"PEMA"
PEMA_SIDER_long$Method<-"SIDER"

head(PEMA_SIDER_long)
###########################################


############################################
#PEMA SIDER_update
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
PEMA_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER_update" & Isotope == "d13C")
PEMA_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER_update" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_PEMA_SIDER_update <- as.data.frame(Type)#data frame for the 5 food sources
TDF_PEMA_SIDER_update$d13C_Av<-PEMA_SIDER_update_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_PEMA_SIDER_update$d13C_Sd<-PEMA_SIDER_update_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_PEMA_SIDER_update$d15N_Av<-PEMA_SIDER_update_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_PEMA_SIDER_update$d15N_Sd<-PEMA_SIDER_update_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_PEMA_SIDER_update<-select(TDF_PEMA_SIDER_update, -Type)#remove names of food sources


#Mixing polygon simulation for PEMA with SIDER_update TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(PEMA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_PEMA_SIDER_update),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_PEMA_SIDER_update[j,1], sd=TDF_PEMA_SIDER_update[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_PEMA_SIDER_update[j,3], sd=TDF_PEMA_SIDER_update[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(PEMA_hair[,1], PEMA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_PEMA_SIDER_update
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(PEMA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="PEMA_SIDER_update.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
PEMA_SIDER_update_data<-as.data.frame(p)
PEMA_SIDER_update_data.df <- tibble::rownames_to_column(PEMA_SIDER_update_data, "Iteration")
head(PEMA_SIDER_update_data.df)
library(tidyr)
PEMA_SIDER_update_long <- gather(PEMA_SIDER_update_data.df, Sample, In_out, V1:ncol(PEMA_SIDER_update_data.df), factor_key=TRUE)#wide to long format
PEMA_SIDER_update_long$Species<-"PEMA"
PEMA_SIDER_update_long$Method<-"SIDER_update"

head(PEMA_SIDER_update_long)
###########################################


############################################
#PEMA Meta_analysis
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
PEMA_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Meta_analysis" & Isotope == "d13C")
PEMA_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Meta_analysis" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_PEMA_Meta_analysis <- as.data.frame(Type)#data frame for the 5 food sources
TDF_PEMA_Meta_analysis$d13C_Av<-PEMA_Meta_analysis_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_PEMA_Meta_analysis$d13C_Sd<-PEMA_Meta_analysis_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_PEMA_Meta_analysis$d15N_Av<-PEMA_Meta_analysis_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_PEMA_Meta_analysis$d15N_Sd<-PEMA_Meta_analysis_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_PEMA_Meta_analysis<-select(TDF_PEMA_Meta_analysis, -Type)#remove names of food sources


#Mixing polygon simulation for PEMA with Meta_analysis TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(PEMA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_PEMA_Meta_analysis),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_PEMA_Meta_analysis[j,1], sd=TDF_PEMA_Meta_analysis[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_PEMA_Meta_analysis[j,3], sd=TDF_PEMA_Meta_analysis[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(PEMA_hair[,1], PEMA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_PEMA_Meta_analysis
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(PEMA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="PEMA_Meta_analysis.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
PEMA_Meta_analysis_data<-as.data.frame(p)
PEMA_Meta_analysis_data.df <- tibble::rownames_to_column(PEMA_Meta_analysis_data, "Iteration")
head(PEMA_Meta_analysis_data.df)
library(tidyr)
PEMA_Meta_analysis_long <- gather(PEMA_Meta_analysis_data.df, Sample, In_out, V1:ncol(PEMA_Meta_analysis_data.df), factor_key=TRUE)#wide to long format
PEMA_Meta_analysis_long$Species<-"PEMA"
PEMA_Meta_analysis_long$Method<-"Meta_analysis"

head(PEMA_Meta_analysis_long)
###########################################


############################################
#PEMA Field
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
PEMA_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Field" & Isotope == "d13C")
PEMA_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Field" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_PEMA_Field <- as.data.frame(Type)#data frame for the 5 food sources
TDF_PEMA_Field$d13C_Av<-PEMA_Field_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_PEMA_Field$d13C_Sd<-PEMA_Field_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_PEMA_Field$d15N_Av<-PEMA_Field_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_PEMA_Field$d15N_Sd<-PEMA_Field_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_PEMA_Field<-select(TDF_PEMA_Field, -Type)#remove names of food sources


#Mixing polygon simulation for PEMA with Field TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(PEMA_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_PEMA_Field),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_PEMA_Field[j,1], sd=TDF_PEMA_Field[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_PEMA_Field[j,3], sd=TDF_PEMA_Field[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(PEMA_hair[,1], PEMA_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_PEMA_Field
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(PEMA_hair, pch=19, cex=1.3)
dev.copy2pdf(file="PEMA_Field.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
PEMA_Field_data<-as.data.frame(p)
PEMA_Field_data.df <- tibble::rownames_to_column(PEMA_Field_data, "Iteration")
head(PEMA_Field_data.df)
library(tidyr)
PEMA_Field_long <- gather(PEMA_Field_data.df, Sample, In_out, V1:ncol(PEMA_Field_data.df), factor_key=TRUE)#wide to long format
PEMA_Field_long$Species<-"PEMA"
PEMA_Field_long$Method<-"Field"

head(PEMA_Field_long)
###########################################

#PEMA model output for all four TDFs
PEMA_model_output<-rbind(PEMA_SIDER_long, PEMA_SIDER_update_long, PEMA_Meta_analysis_long, PEMA_Field_long)
#################################################################################################




#################################################################################################
#BLBR - Mixing Polygon Simulations
#################################################################################################

#Hair samples from live trapping- used for all BLBR simulations
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)#discrimination factors
BLBR<-filter(Hair, Species=="BLBR")#select species
BLBR_hair<-select(BLBR, d13C, d15N)#select d13C and d15N columns for analysis
BLBR_hair<-na.omit(BLBR_hair)


############################################
#BLBR SIDER
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
BLBR_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER" & Isotope == "d13C")
BLBR_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_SIDER <- as.data.frame(Type)#data frame for the 5 food sources
TDF_BLBR_SIDER$d13C_Av<-BLBR_SIDER_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_BLBR_SIDER$d13C_Sd<-BLBR_SIDER_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_BLBR_SIDER$d15N_Av<-BLBR_SIDER_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_BLBR_SIDER$d15N_Sd<-BLBR_SIDER_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_BLBR_SIDER<-select(TDF_BLBR_SIDER, -Type)#remove names of food sources


#Mixing polygon simulation for BLBR with SIDER TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(BLBR_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_BLBR_SIDER),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_BLBR_SIDER[j,1], sd=TDF_BLBR_SIDER[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_BLBR_SIDER[j,3], sd=TDF_BLBR_SIDER[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(BLBR_hair[,1], BLBR_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_BLBR_SIDER
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(BLBR_hair, pch=19, cex=1.3)
dev.copy2pdf(file="BLBR_SIDER.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
BLBR_SIDER_data<-as.data.frame(p)
BLBR_SIDER_data.df <- tibble::rownames_to_column(BLBR_SIDER_data, "Iteration")
head(BLBR_SIDER_data.df)
library(tidyr)
BLBR_SIDER_long <- gather(BLBR_SIDER_data.df, Sample, In_out, V1:ncol(BLBR_SIDER_data.df), factor_key=TRUE)#wide to long format
BLBR_SIDER_long$Species<-"BLBR"
BLBR_SIDER_long$Method<-"SIDER"

head(BLBR_SIDER_long)
###########################################


############################################
#BLBR SIDER_update
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
BLBR_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER_update" & Isotope == "d13C")
BLBR_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER_update" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_SIDER_update <- as.data.frame(Type)#data frame for the 5 food sources
TDF_BLBR_SIDER_update$d13C_Av<-BLBR_SIDER_update_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_BLBR_SIDER_update$d13C_Sd<-BLBR_SIDER_update_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_BLBR_SIDER_update$d15N_Av<-BLBR_SIDER_update_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_BLBR_SIDER_update$d15N_Sd<-BLBR_SIDER_update_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_BLBR_SIDER_update<-select(TDF_BLBR_SIDER_update, -Type)#remove names of food sources


#Mixing polygon simulation for BLBR with SIDER_update TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(BLBR_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_BLBR_SIDER_update),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_BLBR_SIDER_update[j,1], sd=TDF_BLBR_SIDER_update[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_BLBR_SIDER_update[j,3], sd=TDF_BLBR_SIDER_update[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(BLBR_hair[,1], BLBR_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_BLBR_SIDER_update
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(BLBR_hair, pch=19, cex=1.3)
dev.copy2pdf(file="BLBR_SIDER_update.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
BLBR_SIDER_update_data<-as.data.frame(p)
BLBR_SIDER_update_data.df <- tibble::rownames_to_column(BLBR_SIDER_update_data, "Iteration")
head(BLBR_SIDER_update_data.df)
library(tidyr)
BLBR_SIDER_update_long <- gather(BLBR_SIDER_update_data.df, Sample, In_out, V1:ncol(BLBR_SIDER_update_data.df), factor_key=TRUE)#wide to long format
BLBR_SIDER_update_long$Species<-"BLBR"
BLBR_SIDER_update_long$Method<-"SIDER_update"

head(BLBR_SIDER_update_long)
###########################################


############################################
#BLBR Meta_analysis
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
BLBR_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Meta_analysis" & Isotope == "d13C")
BLBR_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Meta_analysis" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_Meta_analysis <- as.data.frame(Type)#data frame for the 5 food sources
TDF_BLBR_Meta_analysis$d13C_Av<-BLBR_Meta_analysis_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_BLBR_Meta_analysis$d13C_Sd<-BLBR_Meta_analysis_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_BLBR_Meta_analysis$d15N_Av<-BLBR_Meta_analysis_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_BLBR_Meta_analysis$d15N_Sd<-BLBR_Meta_analysis_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_BLBR_Meta_analysis<-select(TDF_BLBR_Meta_analysis, -Type)#remove names of food sources


#Mixing polygon simulation for BLBR with Meta_analysis TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(BLBR_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_BLBR_Meta_analysis),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_BLBR_Meta_analysis[j,1], sd=TDF_BLBR_Meta_analysis[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_BLBR_Meta_analysis[j,3], sd=TDF_BLBR_Meta_analysis[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(BLBR_hair[,1], BLBR_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_BLBR_Meta_analysis
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(BLBR_hair, pch=19, cex=1.3)
dev.copy2pdf(file="BLBR_Meta_analysis.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
BLBR_Meta_analysis_data<-as.data.frame(p)
BLBR_Meta_analysis_data.df <- tibble::rownames_to_column(BLBR_Meta_analysis_data, "Iteration")
head(BLBR_Meta_analysis_data.df)
library(tidyr)
BLBR_Meta_analysis_long <- gather(BLBR_Meta_analysis_data.df, Sample, In_out, V1:ncol(BLBR_Meta_analysis_data.df), factor_key=TRUE)#wide to long format
BLBR_Meta_analysis_long$Species<-"BLBR"
BLBR_Meta_analysis_long$Method<-"Meta_analysis"

head(BLBR_Meta_analysis_long)
###########################################


############################################
#BLBR Field
###########################################
#TDF
Est_DFT<- read.csv("Bartlett_estimated_TDF.csv",header=T)#discrimination factors
#Select mammal species and method of estimated TDF
BLBR_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Field" & Isotope == "d13C")
BLBR_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Field" & Isotope == "d15N")
#dataframe of TDFs for each food source
Type<- c("AM Fungi", "Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_Field <- as.data.frame(Type)#data frame for the 5 food sources
TDF_BLBR_Field$d13C_Av<-BLBR_Field_d13C_TDF$TDF_Mean#fill column with TDF for d13C
TDF_BLBR_Field$d13C_Sd<-BLBR_Field_d13C_TDF$TDF_SD#fill column with sd of TDF for d13C			
TDF_BLBR_Field$d15N_Av<-BLBR_Field_d15N_TDF$TDF_Mean#fill column with TDF for d15N			
TDF_BLBR_Field$d15N_Sd<-BLBR_Field_d15N_TDF$TDF_SD#fill column with sd of TDF for d15N			
TDF_BLBR_Field<-select(TDF_BLBR_Field, -Type)#remove names of food sources


#Mixing polygon simulation for BLBR with Field TDF
its <- 1500  #specify the number of iterations ("its")
min_C <- -40  #specify the dimensions and resolution for the mixing region figure
max_C <- -12    #choose values outside the 95% mixing region
min_N <- -10  
max_N <- 21
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(BLBR_hair))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TDF_BLBR_Field),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TDF_BLBR_Field[j,1], sd=TDF_BLBR_Field[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TDF_BLBR_Field[j,3], sd=TDF_BLBR_Field[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(BLBR_hair[,1], BLBR_hair[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
  if (i %% 10 ==0) cat(paste("iteration", i, "\\n")) 
}

#Check variance
Iterations<-data.frame(Iteration=unlist(Par_values[,ncol(Par_values)-1]))#put iterations into dataframe
Variance<-data.frame(Variance=unlist(Par_values[,ncol(Par_values)]))#put variance into dataframe
Iterations_variance<-cbind(Iterations, Variance)#bind iterations with variance
head(Iterations_variance)
library(ggplot2)
ggplot(Iterations_variance, aes(x=Iteration, y=Variance))+#plot iterations versus variance
  geom_line(size=1, color = "blue") 

#visual check
mix_reg <- mix_reg/its
mix_reg[mix_reg==0] <- NA
mix_regt <- t(mix_reg[ncol(mix_reg):1,])
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TDF_BLBR_Field
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(BLBR_hair, pch=19, cex=1.3)
dev.copy2pdf(file="BLBR_Field.pdf")

#proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
BLBR_Field_data<-as.data.frame(p)
BLBR_Field_data.df <- tibble::rownames_to_column(BLBR_Field_data, "Iteration")
head(BLBR_Field_data.df)
library(tidyr)
BLBR_Field_long <- gather(BLBR_Field_data.df, Sample, In_out, V1:ncol(BLBR_Field_data.df), factor_key=TRUE)#wide to long format
BLBR_Field_long$Species<-"BLBR"
BLBR_Field_long$Method<-"Field"

head(BLBR_Field_long)
###########################################

#BLBR model output for all four TDFs
BLBR_model_output<-rbind(BLBR_SIDER_long, BLBR_SIDER_update_long, BLBR_Meta_analysis_long, BLBR_Field_long)
#################################################################################################




#################################################################################################
#Compile mixing space data for all mammal species
#################################################################################################
Mixing_space_data<-rbind(MYGA_model_output, NAIN_model_output, PEMA_model_output, BLBR_model_output) 

head(Mixing_space_data)

Mixing_space_summary<-Mixing_space_data%>%
  group_by(Species, Method, Sample)%>%
  summarise(sum_In = sum(In_out), n = n())%>%#Total number within 95% and total number of model iterations
  mutate(Percent = sum_In/n)%>%#percent of models runs within 95%
  mutate(within_95 = ifelse(Percent < 0.05, 0, 1))#if out less than 5% of the times than 0, if greater than 1

Mixing_space_95<-Mixing_space_summary%>%
  group_by(Species, Method)%>%
  summarise(Total =n(), within_95_sum = sum(within_95))%>%#total number of hair samples and number that are within 95%
  mutate(Proportion = within_95_sum/Total)#proportion of samples within 95%

write.csv(Mixing_space_95, "Mixing_space_summary.csv", row.names = F)
#################################################################################################

















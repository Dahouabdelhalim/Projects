#script for Figure 3 (relationship between brain volume and whole body SB)

library(ggplot2)
library(extrafont)
library(ape)
library(caper)
library(gridExtra)
library(lattice)
library(grid) 

#read in the tree and the data
parrottree <- read.nexus("consensus_parrot_tree.nex")
parrot_data<-read.csv("parrot_comp_data.csv", header=TRUE)
#removing NS values from the dataset
parrot_data_narm<-parrot_data[!is.na(parrot_data$Stand_cage) & !is.na(parrot_data$BSB) & !is.na(parrot_data$Brain_vol)
                              & !is.na(parrot_data$Body_mass),] 
attach(parrot_data_narm)

#make a comparative object
parrot_comp_data <- comparative.data(phy=parrottree, data= parrot_data_narm, names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
parrot_comp_data$dropped

#running the whole body SB and brain volume final hypothesis testing model
m1 <-pgls(BSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, parrot_comp_data, lambda ='ML')

#because the model has >1 predictor variable, we will use the model to predict whole body SB based on those
#variables, and then plot brain volume against these predicted values
predicted<-predict(m1)
#and adding the predicted values of whole body SB into the dataset
parrot_data_narm$predicted_BSB=predicted

#running the oral SB and brain volume final hypothesis testing model
m2 <-pgls(OSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, parrot_comp_data, lambda ='ML')

#because the model has >1 predictor variable, we will use the model to predict oral SB based on those
#variables, and then plot food handling against these predicted values
predicted2<-predict(m2)

#and adding the predicted values of oral SB in
parrot_data_narm$predicted_OSB=predicted2

#plotting (predicted) whole body SB v brain volume
graph1<-ggplot(parrot_data_narm, aes(x=log(Brain_vol), y=predicted_BSB))+
  theme_bw()+
  theme(text = element_text(family="Calibri", size=30)) +#text size
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"), axis.line = element_line(color = "black"),
        plot.background = element_rect(fill="white")) + 
        geom_point(shape=19, cex=5)+
  scale_x_continuous(limits=c(0,3.5 ),breaks=seq(0, 3.5, 1))+#set the limits for the x axis
  scale_y_continuous(limits=c(-0.02,0.4),breaks=seq(0, 0.4, 0.1))+#set the limits for the y axis
  geom_smooth(method=lm, color="black", alpha=0.2) +
  labs(y="Predicted whole body SB prevalence", x="")

#plotting (predicted) oral SB v brain volume
graph2<-ggplot(parrot_data_narm, aes(x=log(Brain_vol), y=predicted_OSB))+
  theme_bw()+
  theme(text = element_text(family="Calibri", size=30)) +#text size
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"), axis.line = element_line(color = "black"),
        plot.background = element_rect(fill="white")) + 
  geom_point(shape=19, cex=5)+
  scale_x_continuous(limits=c(0,3.5 ),breaks=seq(0, 3.5, 1))+#set the limits for the x axis
  scale_y_continuous(limits=c(-0.02,0.4),breaks=seq(0, 0.4, 0.1), position = "right")+#set the limits for the y axis
  geom_smooth(method=lm, color="black", alpha=0.2) +
  labs(y="Predicted oral SB prevalence", x="")

#arrange them side-by-side into one pane 
graphs<-grid.arrange(graph1, graph2, nrow=1)

#save image as a file
png(filename="BSB_and_OSB4graph.png", type="cairo", units="in",
    width=15, height=8, pointsize = 12, res = 150)
plot(graphs)
dev.off()



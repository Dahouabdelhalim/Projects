#script for Figure 2 (relationship between food handling and FDB)
library(ggplot2)
library(extrafont)
library(ape)
library(caper)
#read in tree then data
parrottree <- read.nexus("consensus_parrot_tree.nex")
parrot_data<-read.csv("parrot_comp_data.csv", header=TRUE)

#removing NA values from the dataset
parrot_data_narm<-parrot_data[!is.na(parrot_data$Prop_adult) & !is.na(parrot_data$Prop_female),] 
#make comparative object
parrot_comp_data <- comparative.data(phy=parrottree, data= parrot_data_narm, names.col=Species_name, vcv=TRUE, 
                                     vcv.dim=3, na.omit=F)

attach(parrot_data_narm)

#running the FDB and food handling final hypothesis testing model
m1 <-pgls(FDB~ Prop_adult+ Prop_female+Stand_cage+sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m1)
#because the model has >1 predictor variable, we will use the model to predict FDB based on those
#variables, and then plot food handling against these predicted values
predicted<-predict(m1)

#and adding the predicted values of FDB in
parrot_data_narm$predicted_FDB=predicted

#plotting (predicted) FDB v food handling
graph<-ggplot(parrot_data_narm, aes(x=sqrt(Food_handling), y=predicted_FDB))+
  theme_bw()+
  theme(text = element_text(family="Calibri", size=30)) +#text size
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"), axis.line = element_line(color = "black"),
        plot.background = element_rect(fill="white")) + 
        geom_point(shape=19, cex=5)+
geom_smooth(method=lm, color="black", alpha=0.2) +
  labs(y="Predicted FDB prevalence", x="Percentage diet needing extensive handling (square-root transformed)")

#save image as a file
png(filename="FDB_graph3.png", type="cairo", units="in",
    width=15, height=8, pointsize = 12, res = 150)
plot(graph)
dev.off()



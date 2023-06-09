library(ggtern)
#Import measurements data
Hindlimb <- read.csv("Hindlimb_data_ternaryplot_FINAL.csv")
View(Hindlimb)
pdf(file="Ternary_labels.pdf",width=20,height=20)
#Create a data frame including the femur, tibia and metatarsal data
df = data.frame(femur = Hindlimb$femur,
                tibia = Hindlimb$tibia,
                metatarsal = Hindlimb$Longest.metatarsal,
                Group = Hindlimb$Major.clade.1)
#Create base plot with points
base = ggtern(data=df,aes(femur,tibia,metatarsal))
base + geom_point()
#Plot points in different colours + scale size by femur 
base + 
  geom_point(aes(fill=femur,shape=Group),color='black',size=5) +
  scale_fill_gradientn(colours=c("blue","purple","red"), trans="log") +
  scale_shape_manual(values=c(21,22,23,24,25,20)) +
  #Define limits of scale
  tern_limits(0.7,0.7,0.5)+
  theme_showarrows()+
  #Plot individual groups on separate ternary diagrams
  facet_wrap(~Group)
dev.off()


#Additional options
geom_text(label=Hindlimb$taxon, size=1) #include labels 

#plot femur against tibia
qplot(
  x = femur,
  y = tibia,
  data = df,
  color = df$Group # color by factor color 
)
ggsave(file="femurtibia_incTeleocrater.pdf") 



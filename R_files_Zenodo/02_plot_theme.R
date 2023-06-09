# plot theme
#set the ggplot theme
OYSTERcolor<-rgb(0,94,157,maxColorValue = 255) #rgb colour from OYSTER logo letters https://www.euromarinenetwork.eu/sites/default/files/header%20images/OYSTER_col_logo_spaced_v2.png

OYSTERpalette<- c(rgb(232,89,42, maxColorValue = 255), #red ball logo
                  rgb(140,198,62, maxColorValue = 255), #yellow ball logo
                  rgb(242,189,18, maxColorValue = 255)) #green ball logo

OYSTERpalette4<- c(rgb(232,89,42, maxColorValue = 255), #red ball logo
                   rgb(140,198,62, maxColorValue = 255), #yellow ball logo
                   rgb(242,189,18, maxColorValue = 255), #green ball logo
                   rgb(0,94,157,maxColorValue = 255)) #rgb colour from OYSTER logo letters)

OYSTERtheme<-theme_minimal()+
  theme(text = element_text(colour=OYSTERcolor, family = "sans", size = 10), # arial 10 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1.3)),
        plot.title=element_text(size=rel(1.3), lineheight=.9,
                                face="bold", colour=OYSTERcolor),
        axis.text = element_text(size = rel(1.30)),
        legend.title.align=0.5)


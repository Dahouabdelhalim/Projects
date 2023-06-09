# mtexti <- function(text, side, off = 0.25,
#                    srt = if(side == 2) 90  else
#                      if(side == 4) 270 else 0, ...) {
#   # dimensions of plotting region in user units
#   usr <- par('usr')
#   # dimensions of plotting region in inches
#   pin <- par('pin')
#   # user units per inch
#   upi <- c(usr[2]-usr[1],
#            usr[4]-usr[3]) / pin
#   # default x and y positions
#   xpos <- (usr[1] + usr[2])/2
#   ypos <- (usr[3] + usr[4])/2
#   if(1 == side)
#     ypos <- usr[3] - upi[2] * off
#   if(2 == side)
#     xpos <- usr[1] - upi[1] * off
#   if(3 == side)
#     ypos <- usr[4] + upi[2] * off
#   if(4 == side)
#     xpos <- usr[2] + upi[1] * off
#   text(x=xpos, y=ypos, text, xpd=NA, srt=srt, ...)
# }




blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    legend.position="none",
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


t_col <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
  
}

theme_figs<- theme_bw()+
  theme(panel.border= element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(colour="black", size=12),
        axis.title = element_text(colour="black", size=13),
        axis.line = element_line(color="black", size = 0.5)) 


mytheme_figure<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(face="bold", colour="black", size=11),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5),
        panel.grid = element_blank())



mytheme<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(face="bold", colour="black", size=11),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

mytheme_big<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=20),
        axis.text.x = element_text(face="bold", colour="black", size=20),
        axis.title.y = element_text(face="bold", colour="black", size=20),
        axis.title.x = element_text(face="bold", colour="black", size=20),
        axis.line.y = element_line(color="black", size = 1),
        axis.line.x = element_line(color="black", size = 1),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


mytheme_90<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(colour="black", size=5, angle = 90, vjust=0.5, hjust=1),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

mytheme_45<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(colour="black", size=5, angle = 45, vjust=1, hjust=1),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))



mytheme_45_medium<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(colour="black", size=7, angle = 48, vjust=1.02, hjust=1),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


mytheme_45_big<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(colour="black", size=15, angle = 45, vjust=1, hjust=1),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


mytheme_45_leg<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=10),
    axis.text.x = element_text(colour="black", size=5, angle = 45, vjust=1, hjust=1),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


mytheme_45_leg_medium<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(colour="black", size=5),
    axis.text.x = element_text(colour="black", face = "bold", size=8, angle = 45, vjust=1, hjust=1),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


mytheme_leg_90<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=10),
    axis.text.x = element_text(colour="black", size=5, angle = 90, vjust=0.5, hjust=1),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


mytheme_small<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=5),
        axis.text.x = element_text(face="bold", colour="black", size=5),
        axis.title.y = element_text(face="bold", colour="black", size=5),
        axis.title.x = element_text(face="bold", colour="black", size=5),
        axis.line.y = element_line(color="black", size = 0.3),
        axis.line.x = element_line(color="black", size = 0.3),
        plot.title = element_text(lineheight=.8, face="bold", size = 5, hjust = 0.5))


mytheme_leg<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=10),
    axis.text.x = element_text(face="bold", colour="black", size=11),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

mytheme_leg_smally<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=8),
    axis.text.x = element_text(face="bold", colour="black", size=11),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

mytheme_leg_big<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=20),
    axis.text.x = element_text(face="bold", colour="black", size=20),
    axis.title.y = element_text(face="bold", colour="black", size=20),
    axis.title.x = element_text(face="bold", colour="black", size=20),
    axis.line.y = element_line(color="black", size = 1),
    axis.line.x = element_line(color="black", size = 1),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15))



mytheme_55_big<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=20),
        axis.text.x = element_text(colour="black", size=8, angle = 55, vjust=1, hjust=1),
        axis.title.y = element_text(face="bold", colour="black", size=20),
        axis.title.x = element_text(face="bold", colour="black", size=1),
        axis.line.y = element_line(color="black", size = 1),
        axis.line.x = element_line(color="black", size = 1),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


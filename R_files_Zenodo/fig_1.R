require(ggplot2)
require(ggforce)

p1 <- ggplot()+
  ggtitle("Top view")+xlab("Length (cm)")+ylab("Width (cm)")+
  geom_rect(mapping=aes(xmin=0, xmax=100, ymin=0, ymax=50, fill="Tank"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 50, y0 = 25, r = 18.3/2, fill = "Speaker"), color="black")+
  geom_circle(aes(x0 = 38, y0 = 30, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 38, y0 = 20, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 62, y0 = 30, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 62, y0 = 20, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 45, y0 = 13, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 55, y0 = 13, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 45, y0 = 37, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_circle(aes(x0 = 55, y0 = 37, r = 7/2, fill = "Bottle"), color="black", alpha=0.5)+
  geom_rect(mapping=aes(xmin=37.7, xmax=38.3, ymin=11, ymax=39, fill="Stick"), color="black")+
  geom_rect(mapping=aes(xmin=61.7, xmax=62.3, ymin=11, ymax=39, fill="Stick"), color="black")+
  geom_rect(mapping=aes(xmin=36, xmax=64, ymin=36.7, ymax=37.3, fill="Stick"), color="black")+
  geom_rect(mapping=aes(xmin=36, xmax=64, ymin=12.7, ymax=13.3, fill="Stick"), color="black")+
  # wires left bottom
  geom_line(aes(x = c(35.5,40.5), y = c(17.5,22.5), color = "Wire"))+
  geom_line(aes(x = c(40.5,35.5), y = c(17.5,22.5), color = "Wire"))+
  # Wires left top
  geom_line(aes(x = c(35.5,40.5), y = c(27.5,32.5), color = "Wire"))+
  geom_line(aes(x = c(40.5,35.5), y = c(27.5,32.5), color = "Wire"))+
  # Wires right top
  geom_line(aes(x = c(59.5,64.5), y = c(27.5,32.5), color = "Wire"))+
  geom_line(aes(x = c(64.5,59.5), y = c(27.5,32.5), color = "Wire"))+
  # Wires right bottom
  geom_line(aes(x = c(59.5,64.5), y = c(17.5,22.5), color = "Wire"))+
  geom_line(aes(x = c(64.5,59.5), y = c(17.5,22.5), color = "Wire"))+
  # Wires middle right bottom
  geom_line(aes(x = c(52.5,57.5), y = c(10.5,15.5), color = "Wire"))+
  geom_line(aes(x = c(57.5,52.5), y = c(10.5,15.5), color = "Wire"))+
  # Wires middle left bottom
  geom_line(aes(x = c(42.5,47.5), y = c(10.5,15.5), color = "Wire"))+
  geom_line(aes(x = c(47.5,42.5), y = c(10.5,15.5), color = "Wire"))+
  # Wires middle left top
  geom_line(aes(x = c(42.5,47.5), y = c(34.5,39.5), color = "Wire"))+
  geom_line(aes(x = c(47.5,42.5), y = c(34.5,39.5), color = "Wire"))+
  # Wires middle right top
  geom_line(aes(x = c(52.5,57.5), y = c(34.5,39.5), color = "Wire"))+
  geom_line(aes(x = c(57.5,52.5), y = c(34.5,39.5), color = "Wire"))+
  scale_fill_manual(breaks = c("Bottle", "Speaker", "Stick", "Tank"), 
                    values=c("White", "#00b8ff", "yellow", "Lightblue"), name = "Item")+
  scale_color_manual(breaks = c("Wire"), 
                     values=c("#734326"), name = "")+
  scale_x_continuous(breaks=seq(0,100,10), limits = c(0,100))+
  scale_y_continuous(breaks=seq(0,50,10), limits = c(0,55))+
  theme_bw()+coord_fixed(ratio = 1)

p1
#ggsave('Set-up_top.png', p1, png(), units = "in", width = 8.3/2, height = 2.77)


p2 <- ggplot()+
  ggtitle("Side view")+xlab("Length (cm)")+ylab("Height (cm)")+
  geom_rect(mapping=aes(xmin=0, xmax=100, ymin=0, ymax=50), fill = NA, color="black", alpha=0.5)+
  geom_rect(mapping=aes(xmin=0.15, xmax=99.85, ymin=0.12, ymax=42, fill = "Tank"), color=NA, alpha=0.5)+
  geom_rect(mapping=aes(xmin=34.5, xmax=41.4, ymin=30, ymax=44, fill="Bottle"), color="black", alpha=0.75)+
  geom_rect(mapping=aes(xmin=41.6, xmax=48.5, ymin=30, ymax=44, fill="Bottle"), color="black", alpha=0.75)+
  geom_rect(mapping=aes(xmin=51.5, xmax=58.4, ymin=30, ymax=44, fill="Bottle"), color="black", alpha=0.75)+
  geom_rect(mapping=aes(xmin=58.6, xmax=65.5, ymin=30, ymax=44, fill="Bottle"), color="black", alpha=0.75)+
  # Mussel + valve monitor
  geom_circle(aes(x0 = 38, y0 = 32, r = 1.5, fill = "Mussel"), color = NA)+
  geom_path(aes(x = c(37.7,38,38.3), y = c(33.6,31.8,33.6)), size = 1.1, color = "white")+
  geom_circle(aes(x0 = 45, y0 = 32, r = 1.5, fill = "Mussel"), color = NA)+
  geom_path(aes(x = c(44.7,45,45.3), y = c(33.6,31.8,33.6)), size = 1.1, color = "white")+
  #geom_circle(aes(x0 = 55, y0 = 32, r = 1.5, fill = "Mussel"), color = NA)+
  #geom_path(aes(x = c(54.7,55,55.3), y = c(33.6,31.8,33.6)), size = 1.1, color = "white")+
  geom_circle(aes(x0 = 62, y0 = 32, r = 1.5, fill = "Mussel"), color = NA)+
  geom_path(aes(x = c(61.7,62,62.3), y = c(33.6,31.8,33.6)), size = 1.1, color = "white")+
  # Cables
  geom_path(aes(x = c(36.8,38.5,38.2,-10), 
                y = c(33.0,52.6,53.2,54), color = "Cable"))+
  geom_path(aes(x = c(39.2,38.5,38.2,-10), 
                y = c(33.0,52.6,53.2,55), color = "Cable"))+
  geom_path(aes(x = c(43.8,45.0,45.0,-10), 
                y = c(33.0,53.0,53.2,54), color = "Cable"))+
  geom_path(aes(x = c(46.2,45.0,45.0,-10), 
                y = c(33.0,53.0,53.2,55), color = "Cable"))+
  #geom_path(aes(x = c(53.8,55.0,55.0,-10), 
   #             y = c(33.0,53.0,53.2,54), color = "Cable"))+
  #geom_path(aes(x = c(56.2,55.0,55.0,-10), 
   #             y = c(33.0,53.0,53.2,55), color = "Cable"))+ 
  geom_path(aes(x = c(60.8,62.5,62.2,-10), 
                y = c(33.0,52.6,53.2,54.3), color = "Cable"))+
  geom_path(aes(x = c(63.2,62.5,62.2,-10), 
                y = c(33.0,52.6,53.2,54.7), color = "Cable"))+
  geom_path(aes(x = c(36.2,36.8), y = c(32,33.3), color = "VG Monitor"), size = 1.2)+
  geom_path(aes(x = c(39.2,39.8), y = c(33.3,32), color = "VG Monitor"), size = 1.2)+
  geom_path(aes(x = c(43.2,43.8), y = c(32,33.3), color = "VG Monitor"), size = 1.2)+
  geom_path(aes(x = c(46.2,46.8), y = c(33.3,32), color = "VG Monitor"), size = 1.2)+
  #geom_path(aes(x = c(53.2,53.8), y = c(32,33.3), color = "VG Monitor"), size = 1.2)+
  #geom_path(aes(x = c(56.2,56.8), y = c(33.3,32), color = "VG Monitor"), size = 1.2)+
  geom_path(aes(x = c(60.3,60.8), y = c(32,33.3), color = "VG Monitor"), size = 1.2)+
  geom_path(aes(x = c(63.2,63.8), y = c(33.3,32), color = "VG Monitor"), size = 1.2)+
  # Speaker
  geom_rect(mapping=aes(xmin=42.85, xmax=57.15, ymin=2, ymax=6, fill="Speaker"), color="black", alpha=1)+
  geom_rect(mapping=aes(xmin=40.85, xmax=59.15, ymin=4.5, ymax=5, fill="Speaker"), color="black", alpha=1)+
  geom_rect(mapping=aes(xmin=42.95, xmax=57.05, ymin=2.5, ymax=5.7, fill="Speaker"), color=NA, alpha=1)+
  # Speaker mesh
  geom_line(aes(x = c(42.1,42.1), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(57.9,57.9), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(42.1,57.9), y = c(4,4), color = "Mesh"))+
  geom_line(aes(x = c(42.1,43.3), y = c(1.008,0), color = "Mesh"))+
  geom_line(aes(x = c(42.3,47.3), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(46.3,51.3), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(50.3,55.3), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(54.3,57.9), y = c(4.2,1.176), color = "Mesh"))+
  geom_line(aes(x = c(45.6,42.1), y = c(4.2,1.26), color = "Mesh"))+
  geom_line(aes(x = c(49.6,44.6), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(53.6,48.6), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(57.6,52.6), y = c(4.2,0), color = "Mesh"))+
  geom_line(aes(x = c(57.9,56.6), y = c(1.092,0), color = "Mesh"))+
  # Wires to ceiling + sticks
  geom_line(aes(x = c(38,38), y = c(53.4,65), color = "Wire"))+
  geom_line(aes(x = c(62,62), y = c(53.4,65), color = "Wire"))+
  geom_rect(mapping=aes(xmin=36, xmax=64, ymin=53, ymax=53.6, fill="Stick"), color="black", alpha=1)+
  geom_circle(aes(x0 = 38, y0 = 52.6, r = 0.4, fill = "Stick"), color="black")+
  geom_circle(aes(x0 = 62, y0 = 52.6, r = 0.4, fill = "Stick"), color="black")+
  # Wires 1st
  geom_line(aes(x = c(34.5,37.7), y = c(44,52.6), color = "Wire"))+
  geom_line(aes(x = c(38.3,41.4), y = c(52.6,44), color = "Wire"))+
  # Wires 2nd
  geom_line(aes(x = c(41.6,44.9), y = c(44,53.6), color = "Wire"))+
  geom_line(aes(x = c(45.1,48.5), y = c(53.6,44), color = "Wire"))+
  # Wires 3th
  geom_line(aes(x = c(51.5,54.9), y = c(44,53.6), color = "Wire"))+
  geom_line(aes(x = c(55.1,58.5), y = c(53.6,44), color = "Wire"))+
  # Wires 4th
  geom_line(aes(x = c(58.6,61.7), y = c(44,52.6), color = "Wire"))+
  geom_line(aes(x = c(62.3,65.6), y = c(52.6,44), color = "Wire"))+
  geom_line(aes(x = c(1,99), y = c(42,42), color = "Water level"), size = 1, linetype = "dashed")+
  scale_fill_manual(breaks = c("Bottle", "Mussel", "Speaker", "Stick", "Tank"), 
                    values=c("White", "#707070", "#00b8ff", "yellow", "Lightblue"), name = "Item")+
  scale_color_manual(breaks = c("Cable", "VG Monitor", "Water level", "Wire", "Mesh"), 
                    values=c("Grey", "Darkgreen", "pink", "#088bcc", "#734326"), name = NULL)+ # "#D98594", "#088bcc", "#734326", "Darkgreen"
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(breaks=seq(0,50,10))+
  guides(fill = guide_legend(order = 1),col = guide_legend(order = 2))+coord_fixed(ratio = 1, ylim=c(0, 55), xlim=c(0,100))+
  theme_bw()+theme(legend.margin = margin(-.2,0,-.3,0, unit="cm"))
p2

require(gridExtra)
# Get legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(p2)

# Side and top view next to each other
overview1 <- grid.arrange(p1+theme(legend.position="none"), p2+theme(legend.position="none"), nrow = 2, top = "Experimental arena")
overview <- grid.arrange(overview1, mylegend, ncol = 2, widths = c(0.8,0.2))
overview2 <- grid.arrange(ggplot()+theme_minimal(), overview1, mylegend, ggplot()+theme_minimal(), ncol = 4, widths = c(.13, 0.6, 0.14, .13))
ggsave('Figures/Fig_1_HQ.jpg', overview2, png(), units = "in", width = 8.3, height = 6.7, dpi = 600)


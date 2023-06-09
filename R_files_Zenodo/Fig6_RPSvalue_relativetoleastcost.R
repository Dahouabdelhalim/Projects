
#######################################################################
#Figure 6: The social cost of clean power relative to a fossil future with at pricing.
#needs: demandsummaryrps.csv, basesum*.csv, energysourcesCO2*_rps.csv
#######################################################################
baselineexp <- 1300000000
print(getwd())
  recent <- read.csv(paste("data/demandsummaryrps.csv"), header=TRUE, sep=",")
  basesum <- read.csv(paste("data/basesum0.1.csv"), header=TRUE, sep=",")
  basesum <- basesum[,c("theta","ev","scen","load","totexpb","co2")]
  basesum5 <- read.csv(paste("data/basesum0.5.csv"), header=TRUE, sep=",")
  basesum5 <- basesum5[,c("theta","ev","scen","load","totexpb","co2")]
  basesum2 <- read.csv(paste("data/basesum2.csv"), header=TRUE, sep=",")
  basesum2 <- basesum2[,c("theta","ev","scen","load","totexpb","co2")]
  basesum <- rbind(basesum, basesum2)
  basesum <- rbind(basesum, basesum5)
  pricedf <- read.csv(paste("data/energysourcesCO20.1_rps.csv"), header=TRUE, sep=",")

#Multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#######################################################################################
{
#EV 2016 and Full EV - theta 0.1
sumtab <- recent
colnames(sumtab)[colnames(sumtab) == "renewable_share_2045"] <- "Renewable_Share"
colnames(sumtab)[colnames(sumtab) == "total_cost"] <- "tots"
jointh <- sumtab[which(sumtab$load==2045&sumtab$ev!="half"&sumtab$theta==0.1&sumtab$scen!=2),]
jointh <- merge(jointh, basesum,by=c("scen","load","theta","ev"), all.x=TRUE )
jointh <- transform(jointh,tots=round(tots/totexpb*100,0))

jointh_flat <- jointh[which(jointh$pricing=="flat" & jointh$scen==3),]
jointh_dyn <- jointh[which(jointh$pricing=="dynamic"),]
jointh <- rbind(jointh_flat,jointh_dyn)
jointh <- jointh[,c("pricing","cost","load","theta","ev","scen","Renewable_Share","tots")]

jointh <- transform(jointh, scen=ifelse(jointh$scen==1,"Optimistic", ifelse(jointh$scen==2,"Moderate","Pessimistic")))
jointh <- transform(jointh, pricing=ifelse(jointh$pricing=="dynamic","Dynamic","Flat"))
jointh <- transform(jointh, cost=ifelse(jointh$cost=="current","2016 Cost","2045 Cost"))
jointh<-transform(jointh, ev=ifelse(ev==2016,"EV = 0.5%",ifelse(ev=="full","EV = 100%","EV = 50%")))

jointhf <- jointh
p1 <- ggplot(jointhf[which(jointhf$pricing=="Dynamic"),])+
  geom_line(aes(Renewable_Share, tots, color=factor(scen)))+
  geom_point(aes(Renewable_Share, tots, color=scen), alpha=0.75, shape=2, size=2, show.legend=F)+
  geom_point(data=jointhf[which(jointhf$pricing=="Flat"&jointhf$scen=="Pessimistic"),],
             aes(Renewable_Share, tots), alpha=0.75, shape=1, size=2, show.legend=F)+
  geom_line(data=jointhf[which(jointhf$pricing=="Flat"&jointhf$scen=="Pessimistic"),],aes(Renewable_Share, tots, color="black"))+
  facet_grid(ev~cost, scale = "fixed")+
  #geom_hline(yintercept=70, col = "black", linetype = "dashed") +
  ylim(min(50), max(100))+
  scale_y_continuous(breaks = seq(50, 100, by = 10))+
  xlim(min(0), max(1))+
  scale_colour_manual("", labels=c("Flat pricing","RTP - Optimistic","RTP - Pessimistic"),
                      values = c("black",brewer.pal(6, "Paired")[3:4],"black",brewer.pal(6, "Paired")[3:4]))+
  xlab("Clean Power shares") +
  ylab("Total cost (in % of baseline expenditure)") + theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  theme(legend.position = c(0.62, 0.80),
        legend.title = element_text(colour = "black", size = 10, face = "bold"),
        legend.text = element_text(colour = "black", size = 10), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
p1
#####################################################################################

#greekslegend=list(bquote(theta==.(0.1)~" - Pessimistic"),bquote(theta==.(2)~" - Pessimistic"))
greekslegend=list(bquote(theta==.(2)),bquote(theta==.(thetaval)))

#https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
#\\u03B8 -> theta
sumtab <- recent
colnames(sumtab)[colnames(sumtab) == "renewable_share_2045"] <- "Renewable_Share"
colnames(sumtab)[colnames(sumtab) == "total_cost"] <- "tots"
jointh <- sumtab[which(sumtab$load==2045&sumtab$ev=="half"&
                         (sumtab$theta==0.1|sumtab$theta==2)&sumtab$scen!=2),]
jointh <- merge(jointh, basesum,by=c("scen","load","theta","ev"), all.x=TRUE )
jointh <- transform(jointh,tots=round(tots/totexpb*100,0))

jointh_flat <- jointh[which(jointh$pricing=="flat" & jointh$scen==3),]
jointh_dyn <- jointh[which(jointh$pricing=="dynamic"),]
jointh <- rbind(jointh_flat,jointh_dyn)
jointh <- jointh[,c("pricing","cost","load","theta","ev","scen","Renewable_Share","tots")]

jointh <- transform(jointh, scen=ifelse(jointh$scen==1,"Optimistic", ifelse(jointh$scen==2,"Moderate","Pessimistic")))
jointh <- transform(jointh, pricing=ifelse(jointh$pricing=="dynamic","Dynamic","Flat"))
jointh <- transform(jointh, cost=ifelse(jointh$cost=="current","2016 Cost","2045 Cost"))
jointh <- transform(jointh, ev=ifelse(jointh$ev=="full","Full",ifelse(jointh$ev=="half","Half",2016)))

jointhf <- jointh
p2 <- ggplot(jointhf[which(jointhf$pricing=="Dynamic"),])+
  geom_line(aes(Renewable_Share, tots, color=factor(scen)))+
  geom_point(aes(Renewable_Share, tots, color=scen), alpha=0.75, shape=2, size=2, show.legend=F)+
  geom_point(data=jointhf[which(jointhf$pricing=="Flat"&jointhf$scen=="Pessimistic"),],
             aes(Renewable_Share, tots), alpha=0.75, shape=1, size=2, show.legend=F)+
  geom_line(data=jointhf[which(jointhf$pricing=="Flat"&jointhf$scen=="Pessimistic"),],aes(Renewable_Share, tots, color="black"))+
  facet_grid(theta~cost, scale = "fixed",labeller = label_bquote(theta == .(theta)))+
  #geom_hline(yintercept=50, col = "black", linetype = "dashed") +
  ylim(min(-25), max(100))+
  scale_y_continuous(breaks = seq(-25, 100, by = 25))+
  xlim(min(0), max(1))+
  scale_colour_manual("", labels=c("Flat pricing","RTP - Optimistic","RTP - Pessimistic"),
                      values = c("black",brewer.pal(6, "Paired")[3:4],"black",brewer.pal(6, "Paired")[3:4]))+
  xlab("Clean Power shares") +
  ylab("Total cost (in % of baseline expenditure)") + theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  theme(legend.position = c(0.62, 0.80),
        legend.title = element_text(colour = "black", size = 10, face = "bold"),
        legend.text = element_text(colour = "black", size = 10), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
p2
pdf("plot/Fig6_costofrenewable.pdf", height = 10, width = 8)
multiplot(p1,p2)
dev.off()
}


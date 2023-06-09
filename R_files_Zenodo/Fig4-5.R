#######################################################################
#Figure 4: Surplus gain from real time pricing under different policy, cost and demand flexibility scenarios.
#Figure 5: Cost of 100 percent renewable energy system under different policy, cost and demand flexibility scenarios.
#######################################################################

#Figure 4: Surplus gain from real time pricing

joinall <- graphdata[which(graphdata$theta==thetaval),]
maxvallabel <- max(as.numeric(joinall$dyntot))+5
joinall <- joinall[which(joinall$pricing=="dynamic"),c("category","scen","cost","load","ev","dyntot")]
joinall <- transform(joinall, category=ifelse(joinall$category=="free","unconstrained",paste(joinall$category)))
joinall <- transform(joinall, category=ifelse(joinall$category=="rps_100","100% Clean",paste(joinall$category)))
joinall <- transform(joinall, category=ifelse(joinall$category=="fossil","Fossil",paste(joinall$category)))
joinall <- unique(joinall, fromLast=TRUE)
idlist <- list(joinall$category, joinall$scen, joinall$cost,joinall$load,joinall$ev)
colnames(joinall) <- c("category","scen","cost","load","ev","dyntot")
joinall <- transform(joinall, ev=ifelse(joinall$ev=="2016","0.5%",ifelse(joinall$ev=="half","50%","100%")))
joinall <-  as.data.frame(sapply(joinall, str_to_title))

draw <- as.data.table(joinall)
###############################################################################
#plot benefit dynamic pricing
###############################################################################
draw <- melt(draw, id=c("category","cost","scen","load","ev", value.name="dyntot"))
draw <- transform(draw, idlabel=(paste(draw$category, draw$scen, draw$load, sep = "")))
draw <- dcast(draw,category + cost +  scen + load + idlabel ~ ev, value.var="dyntot")  
colnames(draw) <- c("category","cost","scen","load","idlabel","low","full","half")

#get value from full ev
#get value from half ev
#get value from 0.5ev
draw <- transform(draw,mid=as.numeric(draw$half))
draw <- transform(draw,opt=as.numeric(draw$low))
draw <- transform(draw,pes=as.numeric(draw$full))
draw <- transform(draw, evlabel=ifelse(draw$category=="Unconstrained" & draw$cost=="Current" &
                                         draw$scen==2, "100% EV",
                                       ifelse(draw$category=="100% Clean" & draw$cost=="Current" &
                                                draw$scen==2, "50% EV",
                                              ifelse(draw$category=="Fossil" & draw$cost=="Current" &
                                                       draw$scen==2, "0.5% EV",""))))
draw <- transform(draw, evy=ifelse(draw$category=="Unconstrained" & draw$cost=="Current" &
                                     draw$scen==2, draw$pes,
                                   ifelse(draw$category=="100% Clean" & draw$cost=="Current" &
                                            draw$scen==2, draw$mid,
                                          ifelse(draw$category=="Fossil" & draw$cost=="Current" &
                                                   draw$scen==2, draw$opt,""))))
draw <- transform(draw, evy=as.numeric(as.character(draw$evy)))                                       
draw <- transform(draw, loadlabel=ifelse(draw$load=="2007" & draw$category=="100% Clean" & draw$cost=="Current" &
                                           draw$scen==1, "2007 load",""))
draw <- transform(draw, levelone=factor(1))
draw <- transform(draw, Demand_flexibility=ifelse(draw$scen==1,"1-Optimistic", ifelse(draw$scen==2,"2-Moderate","3-Pessimistic")))
draw <- transform(draw, cost=ifelse(draw$cost=="Current","2016 Cost","2045 Cost"))
draw <- transform(draw, categoryord=ifelse(draw$category=="100% Clean",2, ifelse(draw$category=="Fossil",1,3)))
draw <- transform(draw, category = reorder(category, categoryord))

draw2045 <- draw[which(draw$load==2045),c("category","cost","scen","idlabel","low","full","half","mid","opt","pes","evlabel","evy","loadlabel","levelone","Demand_flexibility","categoryord")]
draw2007 <- draw[which(draw$load==2007),c("category","cost","scen","idlabel","low","full","half","mid","opt","pes","evlabel","evy","loadlabel","levelone","Demand_flexibility","categoryord")]

ggplot(data = draw2045, aes(category, mid, fill=Demand_flexibility, linetype=loadlabel, shape=levelone), shape=18) + 
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(cost~ ., scale = "free_y", space = "fixed") +
  geom_point(data = draw2007, shape=18, aes(category, mid), size=4, color="red", alpha=0.7, position=position_dodge(0.9),show.legend=F) +
  geom_errorbar(aes(ymin=pes, ymax=opt), width=.2, position=position_dodge(.9), show.legend=F ) +
  #geom_text(data = draw2045,aes(category, evy, label = evlabel)) +
  ylim(min(0), max(maxvallabel))+
  #geom_hline(yintercept=10, col = "gray80", linetype = "longdash") +
  #geom_hline(yintercept=20, col = "gray80", linetype = "longdash") +
  #geom_hline(yintercept=30, col = "gray80", linetype = "longdash") +
  #geom_hline(yintercept=40, col = "gray80", linetype = "longdash") +
  scale_size_area() + 
  xlab("") +
  ylab("Percentage points of baseline expenditure") + theme_classic()+
  scale_fill_manual("",values=c(brewer.pal(3,"BuGn")),labels=c("Optimistic","Moderate","Pessimistic"),
                    guide=guide_legend(override.aes = list(color=c(brewer.pal(3,"BuGn")))))+
  geom_hline(mapping=aes(yintercept=c(0), line="95% CI"),fill="white",color="white", alpha=1,show.legend=T)+
  geom_point(mapping=aes(1,1, shape="90% CI"),fill="white", color="black", shape=18, alpha=0,show.legend=T)+
  scale_shape_manual("", values=1,labels=c("EV level (0.5%-100%)"),guide=guide_legend(override.aes = list(linetype="solid", size=1, color="black", alpha=1)))+
  theme(legend.direction="horizontal",
        legend.title = element_text(colour = "black", size = 0),
        legend.text = element_text(colour = "black", size = 10),
        legend.position = "bottom", axis.ticks=element_blank(), 
        plot.title = element_text(hjust = -10))+
  geom_hline(yintercept=0, col = "gray80")+
  scale_linetype_manual("", values=c(1,2),labels=c("2007 Load",""),
                        guide=guide_legend(override.aes = list(linetype="blank", fill="white",size=4, shape=18,color=c("red","white"), alpha=1)))
ggsave(paste("plot/Fig",num1,"_theta=",thetaval,".pdf",sep=""),height=6.1,width=10)
#dev.off()

###############################################################################
#Figure 5: Cost of 100 percent renewable energy system
###############################################################################
joinall <- graphdata[which(graphdata$theta==thetaval),]
rpsbase <- joinall[which(joinall$category=="rps_100"),c("category","scen","cost","pricing","load","ev","tots")]
colnames(rpsbase) <- c("category","scen","cost","pricing","load","ev","totsb")
joinall <- merge(joinall, rpsbase[,c("scen","cost","pricing","load","ev","totsb")], by=c("scen","cost","pricing","load","ev"), all.x=TRUE)
joinall <- transform(joinall, rpsperc = round((tots-totsb),2))
maxvallabel <- max(as.numeric(joinall$rpsperc),na.rm=TRUE)+5
minvallabel <- min(as.numeric(joinall$rpsperc),na.rm=TRUE)-5

joinall <- joinall[,c("category","scen","cost","pricing","load","ev","rpsperc")] 
joinall <- transform(joinall, category=ifelse(joinall$category=="free","unconstrained",paste(joinall$category)))
joinall <- transform(joinall, category=ifelse(joinall$category=="rps_100","100% Clean",paste(joinall$category)))
joinall <- transform(joinall, category=ifelse(joinall$category=="fossil","Fossil",paste(joinall$category)))
idlist <- list(joinall$category, joinall$scen, joinall$cost,joinall$pricing,joinall$load,joinall$ev)
joinall <- aggregate(joinall$rpsperc, idlist, mean)
colnames(joinall) <- c("category","scen","cost","pricing","load","ev","rpsperc")
joinall <- unique(joinall, fromLast=TRUE)
joinall <- transform(joinall, ev=ifelse(joinall$ev=="2016","0.5%",ifelse(joinall$ev=="half","50%","100%")))
joinall <-  as.data.frame(sapply(joinall, str_to_title))

draw <- joinall[joinall$category!="100% Clean" ,]

draw <- melt(data = as.data.table(draw[,c("category","cost","scen","pricing","load","ev","rpsperc")]), 
             id=c("category","cost","scen","pricing","load","ev"))
draw <- transform(draw, idlabel=(paste(draw$category, draw$scen, draw$pricing, draw$load, sep = "")))
draw <- dcast(draw,category + cost +  scen + pricing + load + variable + idlabel ~ ev, value.var="value")[,-6] 
colnames(draw) <- c("category","cost","scen","pricing","load","idlabel","opt","pes","mid")
draw <- transform(draw, Demand_flexibility=ifelse(draw$scen==1,"1-Optimistic", ifelse(draw$scen==2,"2-Moderate","3-Pessimistic")))
draw <- transform(draw, cost=ifelse(draw$cost=="Current","2016 Cost","2045 Cost"))
draw <- transform(draw, pricing=ifelse(draw$pricing=="Dynamic","Real Time",paste(draw$pricing)))
draw <- transform(draw, pricing=paste(draw$pricing, "pricing"))
#draw <- draw[complete.cases(draw),]
draw <- transform(draw,mid=as.numeric(draw$mid))
draw <- transform(draw,opt=as.numeric(draw$opt))
draw <- transform(draw,pes=as.numeric(draw$pes))
draw <- transform(draw, categoryord=ifelse(draw$category=="100% Clean",2, ifelse(draw$category=="Fossil",1,3)))
draw <- transform(draw, category = reorder(category, categoryord))
draw <- transform(draw, loadlabel=ifelse(draw$load=="2007" & draw$category=="100% Clean" & draw$cost=="Current" &
                                           draw$scen==1, "2007 load",""))
draw <- transform(draw, evlabel=factor(1))

#basecase=future cost

draw2045 <- draw[which(draw$load==2045&draw$category!="100% Clean"),]
draw2007 <- draw[which(draw$load==2007&draw$category!="100% Clean"),]

ggplot(data = draw2045,aes(category, mid, fill=Demand_flexibility, linetype=loadlabel, shape=evlabel), shape=18) + 
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(cost~ pricing, scale = "free_y", space = "fixed") +
  geom_errorbar(aes(ymin=pes, ymax=opt), width=.2, position=position_dodge(.9)) +
  geom_point(data = draw2007, shape=18, aes(category, mid), size=3, color="red", alpha=0.7, position=position_dodge(0.9),show.legend=F) +
  ylim(min(minvallabel), max(maxvallabel))+
  #geom_hline(yintercept=100, col = "gray80", linetype = "longdash") +
  #geom_hline(yintercept=50, col = "gray80", linetype = "longdash") +
  #geom_hline(yintercept=-50, col = "gray80", linetype = "longdash") +
  #geom_hline(yintercept=-100, col = "gray80", linetype = "longdash") +
  scale_size_area() + 
  xlab("") +
  ylab("Percentage points of baseline expenditure") + theme_classic()+
  scale_fill_manual("",values=c(brewer.pal(3,"BuGn")),labels=c("Optimistic","Moderate","Pessimistic"),
                    guide=guide_legend(override.aes = list(color=c(brewer.pal(3,"BuGn")))))+
  geom_hline(mapping=aes(yintercept=c(0), line="95% CI"),fill="white",color="white", alpha=1,show.legend=T)+
  geom_point(mapping=aes(1,1, shape="90% CI"),fill="white", color="black", shape=18, alpha=0,show.legend=T)+
  scale_shape_manual("", values=1,labels=c("EV level (0.5%-100%)"),guide=guide_legend(override.aes = list(linetype="solid", size=1, color="black", alpha=1)))+
  scale_linetype_manual("", values=1,labels=c("Load Profile=2007"),guide=guide_legend(override.aes = list(linetype="blank", fill="white",size=4, shape=18,color=c("red"), alpha=1))) +
  theme(legend.direction="horizontal",
        legend.title = element_text(colour = "black", size = 0),
        legend.text = element_text(colour = "black", size = 10),
        legend.position = "bottom", axis.ticks=element_blank(), 
        plot.title = element_text(hjust = -10))+
  geom_hline(yintercept=0, col = "gray80")
ggsave(paste("plot/Fig",num2,"_theta=",thetaval,".pdf",sep=""),height=6.1,width=10)

########END########################################################################


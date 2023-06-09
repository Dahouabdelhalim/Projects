library(plyr)
library(data.table)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggExtra)
library(gtable)

scaleFUN <- function(x) sprintf("%.1f", x)
setwd("~/Dropbox/L/RScripts/Supplemental_endoderm_plots")
source("bjp_functions.R")


# PCA ---------------------------------------------------------------------

bjp<-
theme(
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  axis.text.y =  element_text(size = 10,face = "bold",color = "black"),
  axis.text.x =  element_text(size = 10,face = "bold",color = "black"),
  strip.text.x = element_text(size = 14,face = "bold"),
	axis.title.y = element_text(size = 11,face = "bold"),
  axis.title.x = element_text(size = 11,face = "bold"),
	#strip.text.x = element_blank(),
 	#strip.text.y = element_blank(),
	panel.spacing = unit(1, "lines"),
  strip.text.y = element_text(size = 14,face = "bold"),
	panel.grid.major = element_line(color = "#E5E5E5",size=0.43),
	panel.grid.minor= element_line(color = "#FAFAFA",size=1.07),
  strip.background = element_blank(),
	legend.position="none")


	

a<- fread("Data_fig3A.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
b<- fread("Data_fig3B.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
c<- fread("Data_figS3C.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)



meta<- fread("meta.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)

a$fig="A"
b$fig="B"
c$fig="C"


rbind(a,b,c)->combined
combined %>% select(fig,pc1=Perc_exp_PC1,pc2=Perc_exp_PC2) %>% unique()->varexp

a<-combined%>%filter(fig=="A")%>%
ggplot(., aes(x=PC1, y=PC2, color=as.factor(Day), shape=as.factor(Species)),show.legend = F)+
  geom_point(aes(),size=3) + theme_bw()+bjp+
	xlab(paste("PC1 (",varexp$pc1[1]," % ","of variance) ",sep=""))+
	ylab(paste("PC2 (",varexp$pc2[1]," % ","of variance) ",sep=""))


b<-combined%>%filter(fig=="B")%>%
ggplot(., aes(x=PC1, y=PC2, color=as.factor(Day), shape=as.factor(Species)),show.legend = F)+
  geom_point(aes(),size=3) + theme_bw()+bjp+
	xlab(paste("PC1 (",varexp$pc1[2]," % ","of variance) ",sep=""))+
	ylab(paste("PC2 (",varexp$pc2[2]," % ","of variance) ",sep=""))

c<-combined%>%filter(fig=="C")%>%
ggplot(., aes(x=PC1, y=PC2, color=as.factor(Day), shape=as.factor(Species)),show.legend = F)+
  geom_point(aes(),size=3) + theme_bw()+bjp+
	xlab(paste("PC1 (",varexp$pc1[3]," % ","of variance) ",sep=""))+
	ylab(paste("PC2 (",varexp$pc2[3]," % ","of variance) ",sep=""))


grid.arrange(a,b,c, ncol=3)

#Save as 8.36 x 2.57

 
# Figure 12A ----------------------------------------------------------------
a<- fread("Data_figS12A.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
b<- fread("Data_figS12B.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)

a$sp="Chimp"
b$sp="Human"

rbind(a,b)->c


#chimps
c %>% filter(sp =="Chimp")->d
d$day<- factor(d$day,levels = c("Days 0 to 1","Days 2 to 3","Days 1 to 2","Days 1 to 3"))
d %>% select(line,day,sp)%>%unique()->lines





d$day<- factor(d$day,levels = c("Days 0 to 1","Days 2 to 3","Days 1 to 2","Days 1 to 3"))
ggplot(NULL,aes())+ 
geom_histogram(data=d,aes(x = pval, ..density..), fill="white", colour="black", breaks=seq(0, 1, by = 0.02)) +
	theme_bw() + scale_x_continuous(limits = c(0,1))+facet_wrap(~day,scales="free_y",ncol=2)+bjp+
	theme(plot.title = element_text(face = "bold")) +  labs(x = "Unadjusted p-values") +
	geom_hline(data=lines,aes(yintercept = line), size=1.5, colour = "#E77642")+ylab("Density")+
	scale_y_continuous(labels=scaleFUN)+
	 annotate("text", label='bold(bolditalic(hat(pi)[0][", chimp"])=="0.73")', parse=TRUE, x=0.80, y=(1.25), size = 4, colour = "#E77642") +bjp
  


#Humans
c %>% filter(sp =="Human")->d
d$day<- factor(d$day,levels = c("Days 0 to 1","Days 2 to 3","Days 1 to 2","Days 1 to 3"))
d %>% select(line,day,sp)%>%unique()->lines

d$day<- factor(d$day,levels = c("Days 0 to 1","Days 2 to 3","Days 1 to 2","Days 1 to 3"))
ggplot(NULL,aes())+ 
geom_histogram(data=d,aes(x = pval, ..density..), fill="white", colour="black", breaks=seq(0, 1, by = 0.02)) +
	theme_bw() + scale_x_continuous(limits = c(0,1))+facet_wrap(~day,scales="free_y",ncol=2)+bjp+
	theme(plot.title = element_text(face = "bold")) +  labs(x = "Unadjusted p-values") +
	geom_hline(data=lines,aes(yintercept = line), size=1.5, colour = "#00A4F4")+ylab("Density")+
	scale_y_continuous(labels=scaleFUN)+
	 annotate("text", label='bold(bolditalic(hat(pi)[0][", human"])=="0.23")', parse=TRUE, x=0.80, y=(2), size = 4, colour = "#00A4F4") +bjp
  

#Long format
d$day<- factor(d$day,levels = c("Days 0 to 1","Days 2 to 3","Days 1 to 2","Days 1 to 3"))
ggplot(NULL,aes())+ 
geom_histogram(data=d,aes(x = pval, ..density..), fill="white", colour="black", breaks=seq(0, 1, by = 0.02)) +
	theme_bw() + scale_x_continuous(limits = c(0,1))+facet_wrap(~day,scales="free_y",ncol=4)+bjp+
	theme(plot.title = element_text(face = "bold")) +  labs(x = "Unadjusted p-values") +
	geom_hline(data=lines,aes(yintercept = line), size=1.5, colour = "#00A4F4")+ylab("Density")+
	scale_y_continuous(labels=scaleFUN)+
	 annotate("text", label='bold(bolditalic(hat(pi)[0][", human"])=="0.23")', parse=TRUE, x=0.80, y=(2), size = 4, colour = "#00A4F4") +bjp
  


# SuppFigure 13 ----------------------------------------------------------------
#Theme for plots:
bjp<-
theme(
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  axis.text.y =  element_text(size = 10,face = "bold",color = "black"),
  axis.text.x =  element_text(size = 10,face = "bold",color = "black"),
  strip.text.x = element_text(size = 14,face = "bold"),
	axis.title.y = element_text(size = 14,face = "bold"),
  axis.title.x = element_text(size = 14,face = "bold"),
	#strip.text.x = element_blank(),
 	#strip.text.y = element_blank(),
	panel.spacing = unit(0.8, "lines"),
  strip.text.y = element_text(size = 14,face = "bold"),
	panel.grid.major = element_line(color = "#E5E5E5",size=0.43),
	panel.grid.minor= element_line(color = "#FAFAFA",size=1.07),
  strip.background = element_blank())




a<- fread("Data_figS13A.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
b<- fread("Data_figS13B.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)

names(a)<-c("Species","pval","day","line")
names(b)<-c("Species","pval","day","line")

#Specify order to plot
a$day<- factor(a$day,levels = c("Days 1 to 2","Days 1 to 3","Days 2 to 3"))
b$day<- factor(b$day,levels = c("Days 1 to 2","Days 1 to 3","Days 2 to 3"))

#Label species
a$sp2="Chimpanzee"
b$sp2="Human"



#Plot humans and chimpanzee
rbind(a,b)->comb
comb$Species=revalue(comb$Species,c("(Chimpanzee|Significant in Humans)"="csh","(Human|Significant in Chimps)"="hsc"))

#When creating a named vector for custom palettes, the vector is the color codes, the names are the assignments of the colors
pal1 = as.factor(c("Chimpanzee","Human","csh","hsc"))
pal2 <- c("#E77642","#00A4F4", "#FFFFFF","#FFFFFF") 
setNames(pal2,pal1)->pal


#Separate out species from the conditional overlaps
comb %>% filter(Species =="Human" | Species =="Chimpanzee")->chs
comb %>% filter(Species !="Human" & Species !="Chimpanzee")->nchs

ab<-ggplot(NULL, aes()) +
	geom_density(data=chs,aes(x = pval,fill=Species),show.legend = F, alpha =1,linetype="22",size=0.4) +
  geom_density(data=nchs,aes(x = pval,fill=Species),show.legend = F, alpha = 0.8,linetype=1,size=0.63) +
	geom_density(data=chs,aes(x = pval),show.legend = F, alpha =1,linetype="22",size=0.4) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + theme_bw()+
	labs(x = "Unadjusted p-values") + labs(y = "Density") + scale_y_continuous(labels=scaleFUN) +
 facet_grid(sp2~day,scales = "free")+bjp+
	scale_fill_manual(values = pal)

#Plot C-D
c<- fread("Data_figS13C.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
d<- fread("Data_figS13D.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)

names(c)<-c("Species","pval","day","line")
names(d)<-c("Species","pval","day","line")


#Specify order to plot
c$day<- factor(c$day,levels = c("Days 0 to 1","Days 1 to 3","Days 2 to 3"))
d$day<- factor(d$day,levels = c("Days 0 to 1","Days 1 to 3","Days 2 to 3"))

#Label species
c$sp2="Chimpanzee"
d$sp2="Human"


#Plot humans and chimpanzee
rbind(c,d)->comb
comb$Species=revalue(comb$Species,c("(Chimpanzee|Significant in Humans)"="csh","(Human|Significant in Chimps)"="hsc"))
#When creating a named vector for custom palettes, the vector is the color codes, the names are the assignments of the colors
pal1 = as.factor(c("Chimpanzee","Human","csh","hsc","Shared"))
pal2 <- c("#E77642","#00A4F4", "#FFFFFF","#FFFFFF","#000000") 
setNames(pal2,pal1)->pal


#Separate out species from the conditional overlaps
comb %>% filter(Species =="Human" | Species =="Chimpanzee")->chs
comb %>% filter(Species !="Human" & Species !="Chimpanzee")->nchs


#Make data for lines and value labeling:
comb %>% select(line,day,Species,sp2)%>%unique()->lines
lines$Species=revalue(lines$Species,c("csh"="Shared","hsc"="Shared"))
# Add coordinates for text labels:
lines[which(lines['Species'] =="Chimpanzee" | lines['Species'] =="Human"),"add"] <- 0.1
lines[which(lines['Species'] =="Shared"),"add"] <- (-0.1)

# Make text and final coordinates:
lines%>% filter(!is.na(line))%>%mutate(txt=paste('bold(bolditalic(hat(pi)[0][", ',Species,'"])=="',line,'")',sep=""),yx=1.6+add)->lines


cd<-ggplot(NULL, aes()) +
	geom_density(data=chs,aes(x = pval,fill=Species,color=Species),show.legend = F, alpha =1) +
  geom_density(data=nchs,aes(x = pval,fill=Species),show.legend = F, alpha = 0.8,linetype=1,size=0.63) +
	geom_density(data=chs,aes(x = pval),show.legend = F, alpha =1,linetype="22",size=0.4) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + theme_bw()+
	geom_hline(data=lines,aes(yintercept = line),color="#666666", size = 1.9,show.legend = F)+ 
	geom_hline(data=lines,aes(yintercept = line,color=Species), size = 1.5,show.legend = F) + 
	labs(x = "Unadjusted p-values") + labs(y = "Density") + scale_y_continuous(labels=scaleFUN) +
 facet_grid(sp2~day,scales = "free")+bjp+
	scale_fill_manual(values = pal)+
	scale_color_manual(values = pal)+
	geom_text(data=lines, aes(x=0.6, y=yx,label = txt,color=Species),parse=T,show.legend=F,size=3.5)


grid.arrange(ab,cd, ncol=1)






# SuppFigure 14 ----------------------------------------------------------------
#Theme for plots:
bjp<-
theme(
  panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
  axis.text.y =  element_text(size = 8,face = "bold",color = "black"),
  axis.text.x =  element_text(size = 8,face = "bold",color = "black"),
  strip.text.x = element_text(size = 12,face = "bold"),
	axis.title.y = element_text(size = 14,face = "bold"),
  axis.title.x = element_text(size = 14,face = "bold"),
	panel.spacing = unit(0.65, "lines"),
  strip.text.y = element_text(size = 12,face = "bold"),
	panel.grid.major = element_line(color = "#E5E5E5",size=0.43),
	panel.grid.minor= element_line(color = "#FAFAFA",size=1.07),
  strip.background = element_blank())




e<- fread("Data_figS14E_with_days12.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
f<- fread("Data_figS14F.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
g<- fread("Data_figS14G.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
h<- fread("Data_figS14H.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
i<- fread("Data_figS14I.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
j<- fread("Data_figS14J.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
k<- fread("Data_figS14K.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)
l<- fread("Data_figS14L.csv", sep=",",stringsAsFactors = FALSE,header=T,data.table=FALSE)

names(e)<-c("Species","pval","day","line")
names(f)<-c("Species","pval","day","line")
names(g)<-c("Species","pval","day","line")
names(h)<-c("Species","pval","day","line")
names(i)<-c("Species","pval","day","line")
names(j)<-c("Species","pval","day","line")
names(k)<-c("Species","pval","day","line")
names(l)<-c("Species","pval","day","line")

#Specify order to plot
e$day<- factor(e$day,levels = c("Days 0 to 1","Days 1 to 2","Days 1 to 3","Days 2 to 3"))
f$day<- factor(f$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))
g$day<- factor(g$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))
h$day<- factor(h$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))
i$day<- factor(i$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))
j$day<- factor(j$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))
k$day<- factor(k$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))
l$day<- factor(l$day,levels = c("Days 0 to 1","Days 1 to 2","Days 2 to 3","Days 1 to 3"))

#Label species
e$sp2="Chimpanzee"
f$sp2="Chimpanzee"
i$sp2="Chimpanzee"
j$sp2="Chimpanzee"
g$sp2="Human"
h$sp2="Human"
k$sp2="Human"
l$sp2="Human"

e$f="e"
f$f="f"
i$f="i"
j$f="j"
g$f="g"
h$f="h"
k$f="k"
l$f="l"

#e %>% filter(day !="Days 1 to 2") ->e_f

#Rest:
#Plot humans and chimpanzee
rbind(e,f)->comb1
comb1$Species=revalue(comb1$Species,c("Shared (Chimpanzee|Significant in Humans)"="csh","Shared (Human|Significant in Chimps)"="hsc"))
#When creating a named vector for custom palettes, the vector is the color codes, the names are the assignments of the colors
pal1 = as.factor(c("Chimpanzee","Human","csh","hsc","Shared"))
pal2 <- c("#E77642","#00A4F4", "#FFFFFF","#FFFFFF","#000000") 
setNames(pal2,pal1)->pal

comb1->comb
#Separate out species from the conditional overlaps
comb %>% filter(Species =="Human" | Species =="Chimpanzee")->chs
comb %>% filter(Species !="Human" & Species !="Chimpanzee")->nchs
#Make data for lines and value labeling:
comb %>% select(line,day,Species,f)%>%unique()->lines
lines$Species=revalue(lines$Species,c("csh"="Shared","hsc"="Shared"))
# Add coordinates for text labels:
lines[which(lines['Species'] =="Chimpanzee" | lines['Species'] =="Human"),"add"] <- 0.1
lines[which(lines['Species'] =="Shared"),"add"] <- (-0.1)

# Make text and final coordinates:
lines%>% filter(!is.na(line))%>%mutate(txt=paste('bold(bolditalic(hat(pi)[0][", ',Species,'"])=="',line,'")',sep=""),yx=1.6+add)->lines


ef<-ggplot(NULL, aes()) +
	geom_density(data=chs,aes(x = pval,fill=Species,color=Species),show.legend = F, alpha =1) +
  geom_density(data=nchs,aes(x = pval,fill=Species),show.legend = F, alpha = 0.8,linetype=1,size=0.63) +
	geom_density(data=chs,aes(x = pval),show.legend = F, alpha =1,linetype="22",size=0.4) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + theme_bw()+
	geom_hline(data=lines,aes(yintercept = line),color="#666666", size = 1.9,show.legend = F)+ 
	geom_hline(data=lines,aes(yintercept = line,color=Species), size = 1.5,show.legend = F) + 
	labs(x = "Unadjusted p-values") + labs(y = "Density") + scale_y_continuous(labels=scaleFUN) +
 facet_grid(f~day,scales = "free")+bjp+
	scale_fill_manual(values = pal)+
	scale_color_manual(values = pal)+
	geom_text(data=lines, aes(x=0.6, y=yx,label = txt,color=Species),parse=T,show.legend=F,size=2.5)


rbind(g,h)->comb1
comb1$Species=revalue(comb1$Species,c("Shared (Chimpanzee|Significant in Humans)"="csh","Shared (Human|Significant in Chimps)"="hsc"))
#When creating a named vector for custom palettes, the vector is the color codes, the names are the assignments of the colors
comb1->comb
#Separate out species from the conditional overlaps
comb %>% filter(Species =="Human" | Species =="Chimpanzee")->chs
comb %>% filter(Species !="Human" & Species !="Chimpanzee")->nchs
#Make data for lines and value labeling:
comb %>% select(line,day,Species,f)%>%unique()->lines
lines$Species=revalue(lines$Species,c("csh"="Shared","hsc"="Shared"))
# Add coordinates for text labels:
lines[which(lines['Species'] =="Chimpanzee" | lines['Species'] =="Human"),"add"] <- 0.1
lines[which(lines['Species'] =="Shared"),"add"] <- (-0.1)

# Make text and final coordinates:
lines%>% filter(!is.na(line))%>%mutate(txt=paste('bold(bolditalic(hat(pi)[0][", ',Species,'"])=="',line,'")',sep=""),yx=1.6+add)->lines


gh<-ggplot(NULL, aes()) +
	geom_density(data=chs,aes(x = pval,fill=Species,color=Species),show.legend = F, alpha =1) +
  geom_density(data=nchs,aes(x = pval,fill=Species),show.legend = F, alpha = 0.8,linetype=1,size=0.63) +
	geom_density(data=chs,aes(x = pval),show.legend = F, alpha =1,linetype="22",size=0.4) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + theme_bw()+
	geom_hline(data=lines,aes(yintercept = line),color="#666666", size = 1.9,show.legend = F)+ 
	geom_hline(data=lines,aes(yintercept = line,color=Species), size = 1.5,show.legend = F) + 
	labs(x = "Unadjusted p-values") + labs(y = "Density") + scale_y_continuous(labels=scaleFUN) +
 facet_grid(f~day,scales = "free")+bjp+
	scale_fill_manual(values = pal)+
	scale_color_manual(values = pal)+
	geom_text(data=lines, aes(x=0.6, y=yx,label = txt,color=Species),parse=T,show.legend=F,size=2.5)

rbind(i,j)->comb1
comb1$Species=revalue(comb1$Species,c("Shared (Chimpanzee|Significant in Humans)"="csh","Shared (Human|Significant in Chimps)"="hsc"))
#When creating a named vector for custom palettes, the vector is the color codes, the names are the assignments of the colors

comb1->comb
#Separate out species from the conditional overlaps
comb %>% filter(Species =="Human" | Species =="Chimpanzee")->chs
comb %>% filter(Species !="Human" & Species !="Chimpanzee")->nchs
#Make data for lines and value labeling:
comb %>% select(line,day,Species,f)%>%unique()->lines
lines$Species=revalue(lines$Species,c("csh"="Shared","hsc"="Shared"))
# Add coordinates for text labels:
lines[which(lines['Species'] =="Chimpanzee" | lines['Species'] =="Human"),"add"] <- 0.1
lines[which(lines['Species'] =="Shared"),"add"] <- (-0.1)

# Make text and final coordinates:
lines%>% filter(!is.na(line))%>%mutate(txt=paste('bold(bolditalic(hat(pi)[0][", ',Species,'"])=="',line,'")',sep=""),yx=1.6+add)->lines


ij<-ggplot(NULL, aes()) +
	geom_density(data=chs,aes(x = pval,fill=Species,color=Species),show.legend = F, alpha =1) +
  geom_density(data=nchs,aes(x = pval,fill=Species),show.legend = F, alpha = 0.8,linetype=1,size=0.63) +
	geom_density(data=chs,aes(x = pval),show.legend = F, alpha =1,linetype="22",size=0.4) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + theme_bw()+
	geom_hline(data=lines,aes(yintercept = line),color="#666666", size = 1.9,show.legend = F)+ 
	geom_hline(data=lines,aes(yintercept = line,color=Species), size = 1.5,show.legend = F) + 
	labs(x = "Unadjusted p-values") + labs(y = "Density") + scale_y_continuous(labels=scaleFUN) +
 facet_grid(f~day,scales = "free")+bjp+
	scale_fill_manual(values = pal)+
	scale_color_manual(values = pal)+
	geom_text(data=lines, aes(x=0.6, y=yx,label = txt,color=Species),parse=T,show.legend=F,size=2.5)

rbind(k,l)->comb1
comb1$Species=revalue(comb1$Species,c("Shared (Chimpanzee|Significant in Humans)"="csh","Shared (Human|Significant in Chimps)"="hsc"))
#When creating a named vector for custom palettes, the vector is the color codes, the names are the assignments of the colors

comb1->comb
#Separate out species from the conditional overlaps
comb %>% filter(Species =="Human" | Species =="Chimpanzee")->chs
comb %>% filter(Species !="Human" & Species !="Chimpanzee")->nchs
#Make data for lines and value labeling:
comb %>% select(line,day,Species,f)%>%unique()->lines
lines$Species=revalue(lines$Species,c("csh"="Shared","hsc"="Shared"))
# Add coordinates for text labels:
lines[which(lines['Species'] =="Chimpanzee" | lines['Species'] =="Human"),"add"] <- 0.1
lines[which(lines['Species'] =="Shared"),"add"] <- (-0.1)

# Make text and final coordinates:
lines%>% filter(!is.na(line))%>%mutate(txt=paste('bold(bolditalic(hat(pi)[0][", ',Species,'"])=="',line,'")',sep=""),yx=1.6+add)->lines


kl<-ggplot(NULL, aes()) +
	geom_density(data=chs,aes(x = pval,fill=Species,color=Species),show.legend = F, alpha =1) +
  geom_density(data=nchs,aes(x = pval,fill=Species),show.legend = F, alpha = 0.8,linetype=1,size=0.63) +
	geom_density(data=chs,aes(x = pval),show.legend = F, alpha =1,linetype="22",size=0.4) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + theme_bw()+
	geom_hline(data=lines,aes(yintercept = line),color="#666666", size = 1.9,show.legend = F)+ 
	geom_hline(data=lines,aes(yintercept = line,color=Species), size = 1.5,show.legend = F) + 
	labs(x = "Unadjusted p-values") + labs(y = "Density") + scale_y_continuous(labels=scaleFUN) +
 facet_grid(f~day,scales = "free")+bjp+
	scale_fill_manual(values = pal)+
	scale_color_manual(values = pal)+
	geom_text(data=lines, aes(x=0.6, y=yx,label = txt,color=Species),parse=T,show.legend=F,size=2.5)


grid.arrange(ef,gh,ij,kl, ncol=1)
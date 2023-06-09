###### Parameter ######
chrmax<-25
snum<-(chrmax+1)*(chrmax+2)/2-1

###### Library ######

library(diversitree)
library(Matrix)
library(colorRamps)
library(ggplot2)
library(stats4)
library(dplyr)
library(extrafont)

###### Font preparation ######

subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")

###### Function ######

chr_arm_vec<-function(svec){
  chr_arm_mat<-c()  
  for(i in 1:length(svec)){
    chr_arm_mat<-rbind(chr_arm_mat,chr_arm_cal(svec[i]))
  }
  return(chr_arm_mat)
}

coloring_prob<-function(ps){
  ordered_asr<-sort(ps,decreasing=T)
  range50<-sum(cumsum(ordered_asr)<0.50)+1
  range75<-sum(cumsum(ordered_asr)<0.75)+1
  range90<-sum(cumsum(ordered_asr)<0.90)+1
  pc<-rep("100%",length(ps))
  pc[ps>=ordered_asr[range90]]<-rep("90%")
  pc[ps>=ordered_asr[range75]]<-rep("75%")
  pc[ps>=ordered_asr[range50]]<-rep("50%")
  return(pc)
}

###### Data preparation for karyograph of ASR ######
### ASR information
load("st_nodes_bras_M4_YK2021.Robj")

pdata1<-data.frame(chr_arm_vec(1:snum))
colnames(pdata1)<-c("py","px")
pdata1$ps<-st.nodes.bras[,1]
pdata1$pc<-coloring_prob(st.nodes.bras[,1])
pdata1$pc<-factor(pdata1$pc,level=c("50%","75%","90%","100%"))
pdata1$Taxon<-rep("The 1st node (MRCA)")

pdata2<-data.frame(chr_arm_vec(1:snum))
colnames(pdata2)<-c("py","px")
pdata2$ps<-st.nodes.bras[,2]
pdata2$pc<-coloring_prob(st.nodes.bras[,2])
pdata2$pc<-factor(pdata2$pc,level=c("50%","75%","90%","100%"))
pdata2$Taxon<-rep("The 2nd node")

pdata<-rbind(pdata1,pdata2)

###### Figure karyograph of ASR ######

postscript("Fig_karyograph_plot_EC_ASR_1st_2nd_bras_CM25_w_Poly_XXXXXX.eps", height = 3.8, width = 5.1,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)

xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata)+geom_point(aes(x=px,y=py,size=ps,col=pc)) +
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  scale_size(range=c(0,1.5),limits=c(0,0.2)) +
  facet_grid(Taxon~.,scales="free",space="free") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        legend.position="right",
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        strip.background = element_rect(colour="black"))
dev.off()

###### 95% range of chromosome and arm number ######

range95cal<-function(probs){
  cumprob<-0
  range95<-c()
  chrindex<-1:length(probs)
  while(cumprob<0.95){
    index<-which(probs==max(probs))
    range95<-c(range95,chrindex[index])
    cumprob<-cumprob+probs[index]
    chrindex<-chrindex[-index]
    probs<-probs[-index]
  }
  lower<-min(range95)
  upper<-max(range95)
  return(c(lower,upper))
}

pdata1 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata1 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata2 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata2 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal

###### Data Preparation for tree w/ polyploidy ######

## Phylogeny

phy.bras<-read.tree(file="Brassica_tree_YK2021.tre")
kdata<-read.table(file="Brassicaceae_karyotype_data_YK2021.txt",header=T,sep="\\t")
kdata$state<-scal(kdata$Chromosome_haploid,kdata$Arm_haploid)
unknownsp<-which(is.na(match(phy.bras$tip.label,kdata$Species)))
phy.bras<-drop.tip(phy.bras,tip=unknownsp)
phy.bras$tip.state<-kdata$state[match(phy.bras$tip.label,kdata$Species)]
names(phy.bras$tip.state)<-phy.bras$tip.label

## Remove the tip with chr no out of range in PCM

snum<-(chrmax+1)*(chrmax+2)/2-1
phy.bras<-keep.tip(phy.bras,which(phy.bras$tip.state<=snum))
phy.bras$tip.state<-phy.bras$tip.state[phy.bras$tip.state<=snum]

## Matrix for change to chromosome or arm number

stch.mat<-matrix(rep(0),nrow=snum,ncol=(chrmax*2))
star.mat<-matrix(rep(0),nrow=snum,ncol=(chrmax*2))
stch.vec<-c()
star.vec<-c()
rnum<-1
for(d in 1:chrmax){
  for(j in d:(2*d)){
    stch.vec<-c(stch.vec,d)
    star.vec<-c(star.vec,j)
    stch.mat[rnum,d]<-1
    star.mat[rnum,j]<-1
    rnum<-rnum+1
  }
}


max(chr_arm_vec(phy.bras$tip.state)[,1])
max(chr_arm_vec(phy.bras$tip.state)[,2])

## Inverstigation of chrmosome number 

ch.node.mat<-t(st.nodes.bras) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:(chrmax*2))
istate<-nodestate[phy.bras$edge[,1]-Ntip(phy.bras)]
fstate<-ifelse(phy.bras$edge[,2]-Ntip(phy.bras)>0,
               nodestate[abs(phy.bras$edge[,2]-Ntip(phy.bras))],
               stch.vec[phy.bras$tip.state[phy.bras$edge[,2]]])
max(nodestate)

## Different thresholds for polyploidization

(th1<-which(fstate/istate>1.8))
(th1<-which(fstate/istate>1.7))
(th1<-which(fstate/istate>1.6))
(th2<-which(fstate/istate>1.5))
(th3<-which(fstate/istate>1.4))
(th4<-which(fstate/istate>1.3))
(th5<-which(fstate/istate>1.2))
(th5<-which(fstate/istate>1.1))
cbind(istate,fstate)[th1,]
cbind(istate,fstate)[th2,]
cbind(istate,fstate)[th3,]

## Inverstigation of arm number 

ar.node.mat<-t(st.nodes.bras) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:(chrmax*2))
istate<-nodestate[phy.bras$edge[,1]-Ntip(phy.bras)]
fstate<-ifelse(phy.bras$edge[,2]-Ntip(phy.bras)>0,
               nodestate[abs(phy.bras$edge[,2]-Ntip(phy.bras))],
               star.vec[phy.bras$tip.state[phy.bras$edge[,2]]])
max(nodestate)

## Different thresholds for polyploidization

(th1<-which(fstate/istate>1.8))
(th1<-which(fstate/istate>1.7))
(th1<-which(fstate/istate>1.6))
(th2<-which(fstate/istate>1.5))
(th3<-which(fstate/istate>1.4))
(th4<-which(fstate/istate>1.3))
(th5<-which(fstate/istate>1.2))
(th5<-which(fstate/istate>1.1))
cbind(istate,fstate)[th1,]
cbind(istate,fstate)[th2,]
cbind(istate,fstate)[th3,]

thX<-c(15,17,26,31,63,73,75)

###### Main figure trees with ASR (Fig6) ######

chrmaxfig<-20
colvec<-matlab.like((chrmaxfig*2-1)*10+1)

postscript("Fig_bras_tree_with_ancestral_states_w_Poly_wo_spname_XXXXXX.eps", height = 3.5, width = 5.8,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ew<-0.2
panel_label_par<-c(-42,42,-42,42,2.5)

par(oma=c(1.5,2,2,0))
layout(matrix(c(1,2,3,4),2,2,byrow=T),
       widths=c(3.65,3.65),heights=c(3.65,0.2))

ch.node.mat<-t(st.nodes.bras) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:(chrmax*2))
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[stch.vec[phy.bras$tip.state]*10-9]

plot(phy.bras,"fan",no.margin=T,edge.width=1,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=1.2)
nodelabels(pch=16,col=nodecol,frame="n",cex=1.2)
edgelabels(pch=17,edge=thX,col="black",frame="n",cex=1.2)
text(panel_label_par[1],panel_label_par[2],labels="A",cex=panel_label_par[5])
mtext("Chromosome number evolution",side=3,cex=1)
mtext("Brassicaceae (M4)",side=2,cex=1)

ar.node.mat<-t(st.nodes.bras) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:(chrmax*2))
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.bras$tip.state]*10-9]

plot(phy.bras,"fan",no.margin=T,edge.width=1,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=1.2)
nodelabels(pch=16,col=nodecol,frame="n",cex=1.2)
edgelabels(pch=17,edge=thX,col="black",frame="n",cex=1.2)
text(panel_label_par[1],panel_label_par[2],labels="B",cex=panel_label_par[5])
mtext("Arm number evolution",side=3,cex=1)

indexvec<-matlab.like(chrmaxfig*2)
barplot(rep(2,length(indexvec)),
        axes = F, 
        border=NA,
        space = 0,
        col = indexvec)
axis(1,lwd=1,cex.axis=1,mgp=c(0,0.5,0))

dev.off()

###### Supplimentary figure trees with ASR (S16 Fig) ######

chrmaxfig<-20
colvec<-matlab.like((chrmaxfig*2-1)*10+1)

postscript("Fig_bras_tree_with_ancestral_states_spname_two_datasets_XXXXXX.eps", height = 7.4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ew<-0.2
panel_label_par<-c(-50,50,-50,50,3)
psize<-60
par(oma=c(1.5,2,2,0))
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T),
       widths=c(3.65,3.65),heights=c(3.65,3.65,0.2))

load("Dataset/st_nodes_bras_polyp_CM25_unisample_210208.Robj")
phy.bras<-read.tree(file="Dataset/Brassica_tree.tre")
kdata<-read.table(file="Dataset/Brassicaceae_karyotype_data_2021.txt",header=T,sep="\\t")
kdata$state<-scal(kdata$Chromosome_haploid,kdata$Arm_haploid)
unknownsp<-which(is.na(match(phy.bras$tip.label,kdata$Species)))
phy.bras<-drop.tip(phy.bras,tip=unknownsp)
phy.bras$tip.state<-kdata$state[match(phy.bras$tip.label,kdata$Species)]
names(phy.bras$tip.state)<-phy.bras$tip.label
snum<-(chrmax+1)*(chrmax+2)/2-1
phy.bras<-keep.tip(phy.bras,which(phy.bras$tip.state<=snum))
phy.bras$tip.state<-phy.bras$tip.state[phy.bras$tip.state<=snum]
thX<-c(15,17,26,31,63,73,75)

ch.node.mat<-t(st.nodes.bras) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:(chrmax*2))
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[stch.vec[phy.bras$tip.state]*10-9]

plot(phy.bras,"fan",no.margin=T,edge.width=1,show.tip.label=T,cex=0.2,label.offset=3,x.lim=c(-psize,psize),y.lim=c(-psize,psize)) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=1.2)
nodelabels(pch=16,col=nodecol,frame="n",cex=1.2)
edgelabels(pch=17,edge=thX,col="black",frame="n",cex=1.2)
text(panel_label_par[1],panel_label_par[2],labels="A",cex=panel_label_par[5])
mtext("Chromosome number evolution",side=3,cex=1)
mtext("Brassicaceae (All species, M4)",side=2,cex=1)

ar.node.mat<-t(st.nodes.bras) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:(chrmax*2))
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.bras$tip.state]*10-9]

plot(phy.bras,"fan",no.margin=T,edge.width=1,show.tip.label=T,cex=0.2,label.offset=3,x.lim=c(-psize,psize),y.lim=c(-psize,psize)) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=1.2)
nodelabels(pch=16,col=nodecol,frame="n",cex=1.2)
edgelabels(pch=17,edge=thX,col="black",frame="n",cex=1.2)
text(panel_label_par[1],panel_label_par[2],labels="B",cex=panel_label_par[5])
mtext("Arm number evolution",side=3,cex=1)

## with removal of polyploid species

load("st_nodes_bras_removal_poly_species_M4_YK2021.Robj")
phy.bras<-read.tree(file="Brassica_tree_YK2021.tre")
kdata<-read.table(file="Brassicaceae_karyotype_data_YK2021.txt",header=T,sep="\\t")
kdata<-kdata[kdata$Polyploid=="no",]
kdata$state<-scal(kdata$Chromosome_haploid,kdata$Arm_haploid)
unknownsp<-which(is.na(match(phy.bras$tip.label,kdata$Species)))
phy.bras<-drop.tip(phy.bras,tip=unknownsp)
phy.bras$tip.state<-kdata$state[match(phy.bras$tip.label,kdata$Species)]
names(phy.bras$tip.state)<-phy.bras$tip.label
phy.bras<-keep.tip(phy.bras,which(phy.bras$tip.state<=snum))
phy.bras$tip.state<-phy.bras$tip.state[phy.bras$tip.state<=snum]
thX<-c(44,55)

ch.node.mat<-t(st.nodes.bras) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:(chrmax*2))
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[stch.vec[phy.bras$tip.state]*10-9]

plot(phy.bras,"fan",no.margin=T,edge.width=1,show.tip.label=T,cex=0.2,label.offset=3,x.lim=c(-psize,psize),y.lim=c(-psize,psize)) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=1.2)
nodelabels(pch=16,col=nodecol,frame="n",cex=1.2)
edgelabels(pch=17,edge=thX,col="black",frame="n",cex=1.2)
text(panel_label_par[3],panel_label_par[4],labels="C",cex=panel_label_par[5])
mtext("Brassicaceae w/o reported polyploidy (M4)",side=2,cex=1)

ar.node.mat<-t(st.nodes.bras) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:(chrmax*2))
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.bras$tip.state]*10-9]

plot(phy.bras,"fan",no.margin=T,edge.width=1,show.tip.label=T,cex=0.2,label.offset=3,x.lim=c(-psize,psize),y.lim=c(-psize,psize)) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=1.2)
nodelabels(pch=16,col=nodecol,frame="n",cex=1.2)
edgelabels(pch=17,edge=thX,col="black",frame="n",cex=1.2)
text(panel_label_par[3],panel_label_par[4],labels="D",cex=panel_label_par[5])

indexvec<-matlab.like(chrmaxfig*2)
barplot(rep(2,length(indexvec)),
        axes = F, 
        border=NA,
        space = 0,
        col = indexvec)
axis(1,lwd=1,cex.axis=1,mgp=c(0,0.5,0))

dev.off()


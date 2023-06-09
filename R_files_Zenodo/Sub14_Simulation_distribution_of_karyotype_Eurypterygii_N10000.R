###### Library ######
library(ggpubr)
library(ggplot2)
library(diversitree)
library(Matrix)
library(stats4)


###### Function, q_generation_2020 ######
### making Q-matrix from the 4 parameters of karytype evolution
### chrmax: maximum chromosome number

q_generation_2020<-function(k1,k2,k3,k4,chrmax){
  stanum<-(chrmax+1)*(chrmax+2)/2-1
  q<-matrix(rep(0),nrow=stanum,ncol=stanum)
  
  #        k2
  #        ^
  #        |
  # k4 <---+---> k3
  #        |
  #        v
  #        k1
  #
  
  #transition rates from states with limited directions of transition
  #in the beginning are defined one by one.
  
  q[1,2]<- k3
  q[1,1]<- -k3
  q[2,1]<- k4
  q[2,3]<- k2
  q[2,2]<- -k2 -k4
  q[3,2]<- k1 ##2020
  q[3,4]<- 2*k3
  q[3,3]<- -k1-2*k3 ##2020
  q[4,3]<- k4
  q[4,5]<- k3
  q[4,6]<- k2
  q[4,4]<- -k4-k3-k2
  q[5,4]<- 2*k4
  q[5,7]<- 2*k2
  q[5,5]<- -2*k4-2*k2
  
  for(d in 3:(chrmax-1)){
    #d, chromosome number
    
    #d*(d+1)/2 state with no metacentric (no fission, no M-A transition)
    q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1/2 ##2020
    q[d*(d+1)/2,d*(d+1)/2]<- (-d*k3 -(d*(d-1)*k1/2)) ##2020
    
    for(j in (d+1):(2*d-2)){
      #j, arm number
      #loop for the states with 4 directions
      
      #k,state number
      #up, upper state
      #down, lower state
      #anum, the number of acrocentric
      #mnum, the number of metacentric
      
      k<-d*(d+1)/2+j-d
      up<-(d+1)*(d+2)/2+j-d-1
      down<-d*(d-1)/2+j-d+1
      anum<-2*d-j
      mnum<-j-d
      q[k,k+1]<-anum*k3
      q[k,k-1]<-mnum*k4
      q[k,down]<-anum*(anum-1)*k1/2 ##2020
      q[k,up]<-mnum*k2
      q[k,k]<- -q[k,k+1]-q[k,k-1]-q[k,down]-q[k,up]
    }
    
    #states with 1 acrocentric (no fusion)
    k<-d*(d+1)/2+2*d-1-d
    up<-(d+1)*(d+2)/2+2*d-1-d-1
    q[k,k+1]<-k3
    q[k,k-1]<-(d-1)*k4
    q[k,up]<-(d-1)*k2
    q[k,k]<- -q[k,k+1]-q[k,k-1]-q[k,up]
    
    #states with no acrocentric (no fusion, no A-M transition)
    k<-d*(d+1)/2+2*d-d
    up<-(d+1)*(d+2)/2+2*d-d-1
    q[k,k-1]<-d*k4
    q[k,up]<-d*k2
    q[k,k]<- -q[k,k-1]-q[k,up]
    
  }
  
  d<-chrmax
  
  #state in the chrmax (no fission)
  #no metacentric case (no fission, no M-A transition)
  q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1/2 ##2020
  q[d*(d+1)/2,d*(d+1)/2]<- -d*k3-d*(d-1)*k1/2 ##2020
  
  for(j in (d+1):(2*d-2)){
    #3 directions
    k<-d*(d+1)/2+j-d
    down<-d*(d-1)/2+j-d+1
    anum<-2*d-j
    mnum<-j-d
    q[k,k+1]<-anum*k3
    q[k,k-1]<-mnum*k4
    q[k,down]<-anum*(anum-1)*k1/2 ##2020
    q[k,k]<- -q[k,k+1]-q[k,k-1]-q[k,down]
  }
  
  #one acrocentric case (no fission, no fusion)
  k<-d*(d+1)/2+2*d-1-d
  q[k,k+1]<-k3
  q[k,k-1]<-(d-1)*k4
  q[k,k]<- -q[k,k+1]-q[k,k-1]
  
  #no acrocentric case (no fission, no fusion, no A-M transition)
  k<-d*(d+1)/2+2*d-d
  q[k,k-1]<-d*k4
  q[k,k]<- -q[k,k-1]
  
  return(q)
}



###### Tree preparation ######
##!!Need to run RS17 to define the functions first.

load("phy_o_neot_YK2021.Robj")

chrmax<-35
snum<-(chrmax+1)*(chrmax+2)/2-1

tipnum<-Ntip(phy.o.neot)
edgeindex<-which(phy.o.neot$edge[,2]<=tipnum)
meanedge<-mean(node.depth.edgelength(phy.o.neot)[1:tipnum])
excesslen<-node.depth.edgelength(phy.o.neot)[1:tipnum]-meanedge
phy.o.neot$edge.length[edgeindex]<-phy.o.neot$edge.length[edgeindex]-excesslen

###### Data preparation for MuSSE Eurypterygii ######

load(file="fit_musse_neot_M0_YK2021.Robj")
coef_res<-c(coef(fit.musse.neot)[c(6,5,3,4,1,2)]) #2020 definition
names(coef_res)<-c("k1","k2","k3","k4","lambda","mu")
qmatrix<-q_generation_2020(coef_res["k1"],coef_res["k2"],coef_res["k3"],coef_res["k4"],chrmax)
res_counter<-rep(0,snum)
res_300num<-c()
for(trial in 1:10000){
  sim_traits<-sim.character(phy.sim,qmatrix,x0=353,model="mkn")
  sim_table<-table(sim_traits)
  res_300num<-c(res_300num,sum(sim_traits==300))
  res_counter[as.numeric(names(sim_table))]<-res_counter[as.numeric(names(sim_table))]+as.vector(sim_table)
  if(trial %% 100 ==0){
    print(trial)
  }
}

save(res_counter,file="res_counter_CM35_neot_from353_sim_10000_XXXXXX.Robj")
hist_300num<-table(res_300num)
save(hist_300num,file="hist_300num_CM35_neot_from353_sim_10000_XXXXXX.Robj")

###### Histogram of no of (24, 24) (S11A Fig) ######

hdata<-data.frame(cbind(names(hist_300num),as.vector(hist_300num)))
hdata$X2<-as.numeric(hdata$X2)/10000
hdata$X1<-as.integer(hdata$X1)
xlabel<-expression(paste("% karyotype ",bold("(24, 24)")))
g1<-ggplot(hdata)+geom_bar(aes(x=X1/815*100,y=X2),stat="identity")+
  scale_y_continuous(limit=c(0,0.25),expand=c(0,0))+
  xlim(c(-0.2,3))+
  labs(x=xlabel,y="Probability")+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        title=element_text(size=10))

sum(hdata$X2*hdata$X1) #2.29
length(which(phy.o.neot$tip.state==300))

###### Plot of mean distribution (S11B Fig) ######

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

pdata<-data.frame(chr_arm_vec(1:snum))
colnames(pdata)<-c("py","px")
pdata$ps<-res_counter/815/10000
pdata$pc<-coloring_prob(pdata$ps)
pdata$pc<-factor(pdata$pc,level=c("50%","75%","90%","100%"))
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
g2<-ggplot(pdata)+geom_point(aes(x=px,y=py,size=ps,col=pc)) +
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Mean Proportion",color="Top") +
  scale_size(range=c(0,1.5),limits=c(0,0.05)) +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        axis.title.y =element_text(margin=unit(c(0,0,0,7),"pt")),
        legend.position="right",
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"))

###### Combine figures ######

postscript("Fig_simulation_distribution_karyotype_N10000_XXXXXX.eps", height = 5.1, width = 5.1,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggarrange(g1,g2,ncol=1,labels=c("A","B"),font.label=list(size=20,familty="Arial",face="plain"),hjust=-0.1,vjust=1.3)
dev.off()

max(res_counter/815/10000)
chr_arm_vec(which(res_counter/815/10000==max(res_counter/815/10000)))
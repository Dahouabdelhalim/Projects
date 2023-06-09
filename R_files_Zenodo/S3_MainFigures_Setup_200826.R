####NEW FOSSILIZE FUNCTION ----

#setwd("C:/Users/jos23/Google Drive/PERSONAL_NETWORK-FOSSILIZATION/Manuscript/NetFos_Revision_200713/SupplementaryData")

source("S3_MainFigures_ComparisonFunctions.R")
source("S3_MainFigures_FossilizationFunctions.R")

theme_set(theme_bw())
library(patchwork)

####Packages----
library(reshape2)




####ASSIGNING PROBABILITIES ----

#
SpAtt_names<-as.data.frame(rbind(c("BenthicPelagic","1","Pelagic","Life habit"),
                                 c("BenthicPelagic","2","Mixed","Life habit"),
                                 c("BenthicPelagic","3","Benthic","Life habit"),
                                 c("HardSoft","1","Soft-bodied","Body type"),
                                 c("HardSoft","2","Mixed","Body type"),
                                 c("HardSoft","3","Hard-bodied","Body type")))
colnames(SpAtt_names)<-c("type","level","NewLevel","NewSpAtt")

NetAtt_names<-as.data.frame(rbind(c("trophlev","Mean Tr. Lev."),
                                  c("swtl","Mean TL"),
                                  c("oindex","Mean Omn. Ind."),
                                  c("SOI","SOI"),
                                  c("diam","Diameter"),
                                  c("C","Connectance"),
                                  c("degree","Degree"),
                                  c("cc","Clustering"),
                                  c("avpath","CPL"),
                                  c("S","S")
))
colnames(NetAtt_names)<-c("OldNetAtt","NewNetAtt")

NodeAtt_names<-as.data.frame(rbind(c("SWTL","TL"),
                                   c("OI","OI"),
                                   #c("deg_in","Generality"),
                                   c("deg_in_norm","Generality"),
                                   #c("deg_out","Vulnerability"),
                                   c("deg_out_norm","Vulnerability"),
                                   c("btw","Btw.")
                                   
))
colnames(NodeAtt_names)<-c("OldNodeAtt","NewNodeAtt")


## Assign global parameters
p.eql = 0.5
p.hard = 0.75
p.med = 0.5
p.soft = 0.25
runs = 10

att_list<-c(
  "p.hs"
)

####IMPORT ALL WEBS ----

#Pull in all webs

  Sanak_fw<-read.csv("S1_sanak_nearshore_marine.csv", row.names = 1, stringsAsFactors = FALSE) %>% fixnames.fw
  Sanak_nodes <- read.csv("S1_sanak_nearshore_marine_meta.csv", stringsAsFactors = FALSE) %>% fixnames.nodes
  Sanak_comb<-list(FW=Sanak_fw,NODE=Sanak_nodes)
  
  Weddell_fw<-read.csv("S1_Weddell_sea.csv", row.names = 1, stringsAsFactors = FALSE) %>% fixnames.fw
  Weddell_nodes <- read.csv("S1_Weddell_sea_meta.csv", stringsAsFactors = FALSE) %>% fixnames.nodes
  Weddell_comb<-list(FW=Weddell_fw,NODE=Weddell_nodes)
  
  LR_fw<-read.csv("S1_Little_Rock.csv", row.names = 1, stringsAsFactors = FALSE) %>% fixnames.fw
  LR_nodes <- read.csv("S1_Little_Rock_meta.csv", stringsAsFactors = FALSE) %>% fixnames.nodes
  LR_comb<-list(FW=LR_fw,NODE=LR_nodes)
  
  Caiman_fw<-read.csv("S1_coralreef_caiman_caymanislands.csv", row.names = 1, stringsAsFactors = FALSE)  %>% fixnames.fw
  Caiman_nodes <- read.csv("S1_coralreef_caiman_meta_updated.csv", stringsAsFactors = FALSE) %>% mutate(Species=Nodename) %>% fixnames.nodes
  Caiman_nodes<-Caiman_nodes %>% filter(Species %in% rownames(Caiman_fw)) %>% mutate(Species=Nodename) %>% dplyr::select(Species,class) %>% rename("Class"="class")
  Caiman_nodes<-cbind(X=1:nrow(Caiman_nodes),Caiman_nodes)
  Caiman_comb<-list(FW=Caiman_fw,NODE=Caiman_nodes)
  
  Chilean_fw<-read.csv("S1_chilean.csv", row.names = 1, stringsAsFactors = FALSE)  %>% fixnames.fw
  Chilean_nodes <- read.csv("S1_chilean_meta.csv", stringsAsFactors = FALSE)%>% fixnames.nodes
  Chilean_comb<-list(FW=Chilean_fw,NODE=Chilean_nodes)
  
  Messel_fw<-read.csv("S1_Messel_lake.csv", row.names = 1, stringsAsFactors = FALSE)  %>% fixnames.fw
  Messel_nodes <- read.csv("S1_Messel_lake_meta.csv", stringsAsFactors = FALSE)%>% fixnames.nodes
  Messel_comb<-list(FW=Messel_fw,NODE=Messel_nodes)
  
  Burgess_fw<-read.csv("S1_Burgess.csv", row.names = 1, stringsAsFactors = FALSE)  %>% fixnames.fw
  Burgess_nodes <- read.csv("S1_Burgess_meta.csv", stringsAsFactors = FALSE)%>% fixnames.nodes
  Burgess_comb<-list(FW=Burgess_fw,NODE=Burgess_nodes)
  
  Chengjiang_fw<-read.csv("S1_Chengjiang.csv", row.names = 1, stringsAsFactors = FALSE)  %>% fixnames.fw
  Chengjiang_nodes <- read.csv("S1_Chengjiang_meta.csv", stringsAsFactors = FALSE)%>% fixnames.nodes
  Chengjiang_comb<-list(FW=Chengjiang_fw,NODE=Chengjiang_nodes)


####SET FW AND NODE NAMES TO BE THE SAME
Burgess_comb<-set.sp.names(Burgess_comb)
Caiman_comb<-set.sp.names(Caiman_comb)
Chengjiang_comb<-set.sp.names(Chengjiang_comb)
Chilean_comb<-set.sp.names(Chilean_comb)
LR_comb<-set.sp.names(LR_comb)
Messel_comb<-set.sp.names(Messel_comb)
Sanak_comb<-set.sp.names(Sanak_comb)
Weddell_comb<-set.sp.names(Weddell_comb)


####ADD THE OTHER WEBS
allwebs_comb<-list(
  Burgess=Burgess_comb
  ,Caiman=Caiman_comb
  ,Chengjiang=Chengjiang_comb
  ,Chilean=Chilean_comb
  ,LR=LR_comb
  ,Messel=Messel_comb
  ,Sanak=Sanak_comb
  ,Weddell=Weddell_comb
)


modern_comb<-list(
  Caiman=Caiman_comb
  ,Chilean=Chilean_comb
  ,LR=LR_comb
  ,Sanak=Sanak_comb
  ,Weddell=Weddell_comb
)

fossil_comb<-list(
  Burgess=Burgess_comb
  ,Chengjiang=Chengjiang_comb
  ,Messel=Messel_comb
)

web.order<-c(
  "Burgess"
  ,"Chengjiang"
  ,"Messel"
  ,"LR"
  ,"Caiman"
  ,"Chilean"
  ,"Sanak"
  ,"Weddell"
)

web.agesize.order<-c(
  "Chengjiang"
  ,"Burgess"
  ,"Messel"
  ,"Chilean"
  ,"LR"
  ,"Caiman"
  ,"Weddell"
  ,"Sanak"
)

size.order<-c(
  "Chengjiang"
  ,"Chilean"
  ,"Messel"
  ,"Burgess"
  ,"LR"
  ,"Caiman"
  ,"Weddell"
  ,"Sanak"
)

net.var.order<-c()


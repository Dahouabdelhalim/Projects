setwd("~/EMBER_GcGc_Data_Nov_2019/Final_Clinical_Data_and_Peak_Table_April2020")
dd<-read.table(file="Pat_ID_Group.csv",header=TRUE,sep=",")
library(randomizr)
Z <- complete_ra(N = nrow(dd))
table(Z)
dd$Grp_Aloc<-Z

blocks<-paste(dd$Grp_Aloc,dd$Diagnosis)
Z <- block_ra(blocks = blocks, prob_each =  c(0.4693141, 1-0.4693141))
dd$alloc_ind<-as.numeric(Z==0)
for(i in 1:nrow(dd)){
if(dd$alloc_ind[i]==0){dd$Data_Set[i]<-"Non-Template"}
else if(dd$alloc_ind[i]==1){dd$Data_Set[i]<-"Template"}
}
table(dd$Data_Set)
table(dd$Data_Set,dd$Final_adjudicated_diagnosis)
table(dd$Final_adjudicated_diagnosis)/277
prop.table(table(dd$Data_Set,dd$Final_adjudicated_diagnosis),margin=1)
dd_template<-subset(dd,dd$Data_Set=="Template")
dd_nontemplate<-subset(dd,dd$Data_Set=="Non-Template")
#intersect should be 0
intersect(as.character(dd_template$Study.number),as.character(dd_nontemplate$Study.number))
table(dd_template$Final_adjudicated_diagnosis)
table(dd_nontemplate$Final_adjudicated_diagnosis)

write.table(dd_template,file="Subjects_for_Template_n130.csv",row.names=FALSE,sep=",")


#02 basic general changes
D<-D%>%
  mutate(gender=ifelse(Q9==0, "Female", ifelse(Q9==1, "Male", "Other")))%>%
  mutate(stage=as.factor(ifelse(Q6==0, "PhD student",
                      ifelse(Q6==1, "Postdoctoral researcher","Research assistant/Researcher/other ECR"))))
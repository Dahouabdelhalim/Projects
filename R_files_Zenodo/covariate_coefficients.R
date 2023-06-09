######
setwd("~/Dropbox/Things we generated during Shuxi's visit/R code/Codes_for_Covariate_Coefficients_show")
library(ggplot2)
rm(list=ls())
f=1;adv_id=8
for (m in c(1:5)){
source("varcodebook_0923.R")
load(file=paste(mediation.index[m],"covariate_info.RData",sep="_"))
theta_M_mean=apply(theta_M,2,mean)[-1]/X_M_sd
theta_M_up=apply(theta_M,2,FUN=function(x){quantile(x,0.975)})[-1]/X_M_sd
theta_M_down=apply(theta_M,2,FUN=function(x){quantile(x,0.025)})[-1]/X_M_sd


library(ggplot2)
df=data.frame(name=model_x_name,
              mean=theta_M_mean,
              min=theta_M_down,
              max=theta_M_up)

#df$name=factor(df$name,levels=df$name[order(df$mean)])
df$name=factor(df$name,levels=df$name)
p <- ggplot(df,aes(name,mean,shape=name))

#Added horizontal line at y=0, error bars to points and points with size two
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=min, ymax=max), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=rep(1,12))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Covariate") + ylab("Size of effect on mediator variable")
p <- p+ coord_flip()

#pdf(file=paste("../../Figures/Coefficients_on_Covariates/",mediation.index[m],"_covariate_coefficient.pdf",sep=""))
ggsave(p,file=paste("../../Figures/Coefficients_on_Covariates/",mediation.index[m],"_covariate_coefficient.pdf",sep=""))
#dev.off()

write.csv(df,file=paste("../../Tables/covariate_coefficient/",mediation.index[m],"_covariate_coefficient.csv",sep=""))

if(m==1){
theta_Y_mean=apply(theta_Y,2,mean)[-c(1,ncol(theta_Y))]/X_Y_sd
theta_Y_up=apply(theta_Y,2,FUN=function(x){quantile(x,0.975)})[-c(1,ncol(theta_Y))]/X_Y_sd
theta_Y_down=apply(theta_Y,2,FUN=function(x){quantile(x,0.025)})[-c(1,ncol(theta_Y))]/X_Y_sd


library(ggplot2)
df=data.frame(name=model_y_name,
              mean=theta_Y_mean,
              min=theta_Y_down,
              max=theta_Y_up)

#df$name=factor(df$name,levels=df$name[order(df$mean)])
df$name=factor(df$name,levels=df$name)
p <- ggplot(df,aes(name,mean,shape=name))

#Added horizontal line at y=0, error bars to points and points with size two
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=min, ymax=max), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=rep(1,12))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Covariate") + ylab("Size of effect on outcome variable (fGCs)")
p <- p+ coord_flip()
 
#pdf(file="../../Figures/Coefficients_on_Covariates/gc_covariate_coefficient.pdf",onefile=FALSE)
ggsave(p, file="../../Figures/Coefficients_on_Covariates/gc_covariate_coefficient.pdf")
#dev.off()

write.csv(df,file="../../Tables/covariate_coefficient/gc_covariate_coefficient.csv")
}
}
# # # 
# rm(list=ls())
# f=1;adv_id=8;
# for (m in c(1:5)){
# source("varcodebook_0923.R")
# load(paste(mediation.index[m],"8_1_Collect_Details.RData",sep="_"))
# theta_M=M_MCMC_Result$theta
# theta_Y=Y_MCMC_Result$theta
# 
# X_M = model.matrix(as.formula(paste(
#   "~", paste(covariate.index.dsi, collapse = "+")
# )), Complete.data)
# model_x_name=colnames(X_M)[-1]
# X_M=X_M[,-1]
# X_M=scale(X_M)
# X_M_mean=attr(X_M,"scaled:center")
# X_M_sd=attr(X_M,"scaled:scale")
# 
# X_Y = model.matrix(as.formula(paste(
#   "~", paste(covariate.index, collapse = "+")
# )), Complete.data)
# model_y_name=colnames(X_Y)[-1]
# X_Y=X_Y[,-1]
# #Centralize
# X_Y=scale(X_Y)
# X_Y_mean=attr(X_Y,"scaled:center")
# X_Y_sd=attr(X_Y,"scaled:scale")
# 
# 
# save(theta_M, theta_Y, X_M_sd, X_Y_sd, model_x_name, model_y_name,
#      file=paste(mediation.index[m],"covariate_info.RData",sep="_"))
# }
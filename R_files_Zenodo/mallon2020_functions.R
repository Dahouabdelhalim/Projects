# to subset a vector of columns
subset_cols <- function(df, col.list){
  idx <- match(col.list, names(df))
  if(is.data.table(df)==T){
    df <- df[,na.omit(idx), with = F]
  }else{
    df <- df[,na.omit(idx)]
  }
  return(df)
}

# to plot an ellipse (dependent function)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# to plot data and color based on group, rename group, and plot ellipse
nmds_ellipse_plot<-function(group1, group2 = NA, MDS, group.names = NA, 
                            xlim = NA, ylim = NA){
  
  if(is.na(group2[1])){
    
    NMDS = data.frame(NMDS1 = MDS$points[,1], 
                      NMDS2 = MDS$points[,2],
                      group = as.factor(group1))
  }else{
    NMDS = data.frame(NMDS1 = MDS$points[,1], 
                      NMDS2 = MDS$points[,2],
                      group1 = as.factor(group1),
                      group2 = as.factor(group2),
                      group = interaction(group1, group2))
    NMDS$group<-droplevels(NMDS$group)
  }
  
  
  ord<-plot.ellipse.invisible(MDS,NMDS$group)
  
  df_ell <- data.frame()
  for(g in levels(NMDS$group)){
    df_ell <- rbind(df_ell, 
                    cbind(as.data.frame(
                      with(NMDS[NMDS$group==g,],
                           veganCovEllipse(ord[[g]]$cov,
                                           ord[[g]]$center,
                                           ord[[g]]$scale))), 
                      group=g))
  }
  
  
  p<-ggplot() + 
    geom_point(data = NMDS, aes(NMDS1, NMDS2), color = "grey", size = 0.75)
  p<-p +
    geom_point(data = subset(NMDS, !is.na(group)), 
               aes(NMDS1, NMDS2, color = group), size = 0.85) +
    geom_path(data=df_ell, 
              aes(x=NMDS1, y=NMDS2,colour=group), 
              size=0.85, linetype=1)+
    theme_classic()+ theme(legend.title=element_blank())
  
  if(!is.na(group.names[1])){
    grps<-levels(NMDS$group)
    p<- p + scale_colour_discrete(breaks=c(grps),
                                  labels=c(group.names)) 
  }
  
  if(!(is.na(ylim[1]))){
    p<- p+ xlim(xlim)+ylim(ylim)
  }
  
  print(p)
  
}

# dependent function for ellipse plotting
plot.ellipse.invisible <- function(MDS, MDSgroup){
  ff <- tempfile()
  png(filename=ff)
  plot.new()
  ord<-ordiellipse(MDS, MDSgroup, display = "sites", draw = "none", 
                   kind = "se", conf = 0.90, label = F)
  dev.off()
  unlink(ff)
  return(ord)
}

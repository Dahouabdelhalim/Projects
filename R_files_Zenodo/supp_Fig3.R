# load packages
library(rethinking)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)


###################################################
# import data
d <- read.csv("survival_dataset.csv")

# load fitted models
load("survival_models.RDATA")
###################################################

# model names that you need to loop through
models <- c("mgg", "mgr", "msg", "msr", "mf")

# Fix age to 10 and 30 (standardized)
age_seq <- c(-0.3898737, 2.3499895)

#  Fix rank: the middle of the hierarchy, standardized
rank <- 0.2409084

# fix proportion of year the mother was alive for each age, standardize
# for 10-year-old fix 1 (mother alive entire year), for 30-year-old fix 0 (standardized)
mother <- c(1.0667025, -0.9843545)

# fix number to average for each age group: 0 for both, since it is centered
daughter <- 0

# fix median group size, standardized
group_size <- 0.06674918

# create an empty list to store sociality measures for each behavior
soc_seq_true <- list()

# compute predictions:
for( i in models){
  
  # extract posterior samples
  post <- extract.samples( eval( parse( text=paste0( i ) ) ) )
  
  # make empty lists t store the results
  results_lambda_median <- list()
  results_lambda_HPDI <- list()
  
  # create sociality sequence from SURV true_rates
  mean_true_rates <- sapply(1:dim(d)[1], function(x) mean(post$true_rate[1:2000, x ] ) )
  soc_seq_true[[ which(i==models) ]] <- seq( min(mean_true_rates), max(mean_true_rates), length.out = 20 )

  # loop through sociality measures and compute predictions
  for( a in 1:length( age_seq ) ){
    
    mu_link <- function(x) exp( post$base_rate + 
                                  post$b_sociality * x + 
                                  post$b_rank * rank + 
                                  post$b_age * age_seq[ a ] + # either 10 or 30
                                  post$b_mother * mother[ a ] + # either alive or dead
                                  post$b_daughters * daughter +
                                  post$b_group_size * group_size )
                                  
    
    mu <- sapply( soc_seq_true[[ which(i==models) ]], mu_link ) 
    
    # number of days
    mu_median <- apply( mu , 2 , median )
    mu_HPDI <- apply( mu , 2 , HPDI , prob=0.95 )
    
    # daily death rate
    lambda_day <- 1/mu
    lambda <- 1 - (1 - lambda_day)^365 # approx. annual probability of dying

    # approx annual probability of dying
    lambda_median <- apply( lambda , 2 , median )
    lambda_HPDI <- apply( lambda , 2 , HPDI , prob=0.95 )
    
     # save results
    results_lambda_median[[ a ]] <- lambda_median 
    results_lambda_HPDI[[ a ]] <- lambda_HPDI
  
    
  }
  
  # store the results outside of the loop
  assign( paste0(i, "_med"), results_lambda_median )
  assign( paste0(i, "_HPDI"), results_lambda_HPDI )
  
}

# the data is ready for plotting

# Unstandardize sociality sequence
g <- read.csv("groom_individual_post_estimates.csv", stringsAsFactors = F) 
s <- read.csv("support_individual_post_estimates.csv", stringsAsFactors = F) 
p <- read.csv("forage_individual_post_estimates.csv", stringsAsFactors = F) 


soc_nat_scale <- list()
soc_nat_scale[[1]] <- mean( g$post_mean_giver[which(g$sample=="female")] ) + soc_seq_true[[1]] * sd( g$post_mean_giver[which(g$sample=="female")] )
soc_nat_scale[[2]] <- mean( g$post_mean_receiver[which(g$sample=="female")] ) + soc_seq_true[[2]] * sd( g$post_mean_receiver[which(g$sample=="female")] )
soc_nat_scale[[3]] <- mean( s$post_mean_giver[which(s$sample=="female")] ) + soc_seq_true[[3]] * sd( s$post_mean_giver[which(s$sample=="female")] )
soc_nat_scale[[4]] <- mean( s$post_mean_receiver[which(s$sample=="female")] ) + soc_seq_true[[4]] * sd( s$post_mean_receiver[which(s$sample=="female")]  )
soc_nat_scale[[5]] <- mean( p$post_mean[which(p$sample=="female")] ) + soc_seq_true[[5]] * sd( p$post_mean[which(p$sample=="female")] )


###########################################################################################

# adjust plot parameters for prediction plot
shade_alpha = 0.2 # transparency of HPDI 
y2 <-1 # y axis limit
line_width <- 2 # line width
color1 <- "antiqueqhite2"
color2 <- "cadetblue3"
color3 <- "gold"
shade <- "burlywood"
# the colors might slightly differ from the original colors in the paper

# axes titles
axes_title_size <- 16 # 
panel_label_size <- 8

# change soc_seq to natural scale (if necessary, otherwise skip this line)
soc_seq <- soc_nat_scale

# PLOT
png("Figure3.png",width=750, height=850 )

# make a legend
legend_data <- data.frame(sociality=rep(c("10-year-old female","30-year-old female"),1),
                          prob = c(0.35,0.5) )

p0 <- ggplot(legend_data, aes(x=sociality,y=prob, fill=sociality) )+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c(color3,color2), name="")+ #"sociality level\\npercentile"
  theme(legend.text=element_text(size=20),)+ # pdf=20, png=24
  theme_bw() +
  theme(legend.title = element_blank())+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = "right")

legend <- cowplot::get_legend(p0)



# plot for grooming giving predictions
p1 <- ggplot(NULL, aes(x=soc_seq[[1]], y=mgg_med[[1]] ) )+
  geom_line(color=color3, lwd=line_width)+
  geom_line( aes(x=soc_seq[[1]], y= mgg_med[[2]]), color=color2, lwd=line_width ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=mgg_HPDI[[1]][1,], ymax=mgg_HPDI[[1]][2,]), alpha=shade_alpha)+
  geom_ribbon( aes(ymin=mgg_HPDI[[2]][1,], ymax=mgg_HPDI[[2]][2,]), alpha=shade_alpha)+
  ylim(0, y2) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="Grooming giving posterior mean", breaks=c(-4:4))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=11) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))+
  annotate("text",  x=Inf, y =Inf, label = "A", hjust=1.5,vjust=1.5,size=panel_label_size)

# plot for grooming receive predictions
p2 <- ggplot(NULL, aes(x=soc_seq[[2]], y=mgr_med[[1]] ) )+
  geom_line(color=color3, lwd=line_width)+
  geom_line( aes(x=soc_seq[[2]], y= mgr_med[[2]]), color=color2, lwd=line_width ) +
  theme_bw() +
  scale_x_continuous(name ="Grooming receiving posterior mean", breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=mgr_HPDI[[1]][1,], ymax=mgr_HPDI[[1]][2,]), alpha=shade_alpha)+
  geom_ribbon( aes(ymin=mgr_HPDI[[2]][1,], ymax=mgr_HPDI[[2]][2,]), alpha=shade_alpha)+
  ylim(0, y2) +
  ylab("")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=11) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))+
  annotate("text",  x=Inf, y =Inf, label = "B", hjust=1.5,vjust=1.5,size=panel_label_size)


# plot for support giving predictions
p3 <- ggplot(NULL, aes(x=soc_seq[[3]], y=msg_med[[1]] ) )+
  geom_line(color=color3, lwd=line_width)+
  geom_line( aes(x=soc_seq[[3]], y= msg_med[[2]]), color=color2, lwd=line_width ) +
  theme_bw() +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=msg_HPDI[[1]][1,], ymax=msg_HPDI[[1]][2,]), alpha=shade_alpha)+
  geom_ribbon( aes(ymin=msg_HPDI[[2]][1,], ymax=msg_HPDI[[2]][2,]), alpha=shade_alpha)+
  ylim(0, y2) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="Support giving posterior mean", breaks=c(-2:3))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=11) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))+
  annotate("text",  x=Inf, y =Inf, label = "C", hjust=1.5,vjust=1.5,size=panel_label_size)

# plot for support receiving predictions
p4 <- ggplot(NULL, aes(x=soc_seq[[4]], y=msr_med[[1]] ) )+
  geom_line(color=color3, lwd=line_width)+
  geom_line( aes(x=soc_seq[[4]], y= msr_med[[2]]), color=color2, lwd=line_width ) +
  theme_bw() +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=msr_HPDI[[1]][1,], ymax=msr_HPDI[[1]][2,]), alpha=shade_alpha)+
  geom_ribbon( aes(ymin=msr_HPDI[[2]][1,], ymax=msr_HPDI[[2]][2,]), alpha=shade_alpha)+
  ylim(0, y2) +
  ylab("")+
  scale_x_continuous(name ="Support receiving posterior mean", breaks=c(-1,-0.5,0,0.5,1))+
  #ggtitle("Groom receive")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=11) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))+
  annotate("text",  x=Inf, y =Inf, label = "D", hjust=1.5,vjust=1.5,size=panel_label_size)

# plot for groom giving predictions
p5 <- ggplot(NULL, aes(x=soc_seq[[5]], y=mf_med[[1]] ) )+
  geom_line(color=color3, lwd=line_width)+
  geom_line( aes(x=soc_seq[[5]], y= mf_med[[2]]), color=color2, lwd=line_width ) +
  theme_bw() +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=mf_HPDI[[1]][1,], ymax=mf_HPDI[[1]][2,]), alpha=shade_alpha)+
  geom_ribbon( aes(ymin=mf_HPDI[[2]][1,], ymax=mf_HPDI[[2]][2,]), alpha=shade_alpha)+
  ylim(0, y2) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="Foraging posterior mean", breaks=c(-3:3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=11) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))+
  annotate("text",  x=Inf, y =Inf, label = "E", hjust=1.5,vjust=1.5,size=panel_label_size)


# arrange the plots
grid.arrange( p1, p2, p3, p4, p5, legend, nrow=3)

dev.off()






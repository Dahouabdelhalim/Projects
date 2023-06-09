

library(ape)
library(deeptime)
library(ggplot2)
library(ggthemes)

# Phylogenies
bifI <- read.tree("bivalve_family_tree_bifI_dates.tre")
bifII <- read.tree("bivalve_family_tree_bifII_dates.tre")
budI <- read.tree("bivalve_family_tree_budI_dates.tre")
budII <- read.tree("bivalve_family_tree_budII_dates.tre")

drops <- c("Scaphopoda",  "Gastropoda","Gastropoda",    
"Cephalopoda" , "Cephalopoda"  ,"Chitonidae"  ,"Monoplacophora")

bifI <- drop.tip(bifI, drops)
bifII <- drop.tip(bifII, drops)
budI <- drop.tip(budI, drops)
budII <- drop.tip(budII, drops)

cal3.res <- readRDS("ttree.RDS")
family.dat <- read.csv("family_age_placement_data.csv")
extinct <- family.dat$family[family.dat$source=="FOSSIL"]

cal3.res <- lapply(cal3.res, drop.tip, extinct)

# Coordinates for the LTT plot
cc <- ltt.plot.coords(bifI)
cc <- as.data.frame(cc)
cc$time <- abs(cc$time)
cc$`Diversification type` <- "Bifurcation I"

cc1 <- ltt.plot.coords(bifII)
cc1 <- as.data.frame(cc1)
cc1$time <- abs(cc1$time)
cc1$`Diversification type` <- "Bifurcation II"

cc2 <- ltt.plot.coords(budI)
cc2 <- as.data.frame(cc2)
cc2$time <- abs(cc2$time)
cc2$`Diversification type` <- "Budding I"

cc3 <- ltt.plot.coords(budII)
cc3 <- as.data.frame(cc3)
cc3$time <- abs(cc3$time)
cc3$`Diversification type` <- "Budding II"

cc.all <- rbind(cc, cc1, cc2, cc3)

cc.store <- cc.all

for(i in 1:length(cal3.res)){
  bind <- ltt.plot.coords(cal3.res[[i]])
  bind <- as.data.frame(bind)
  bind$time <- abs(bind$time)
  bind$`Diversification type` <- paste0("cal3_", i)
  cc.all <- rbind(cc.all, bind)
}


max.age <- max(cc.all$time)


########

geo.ltt.plot <- ggplot(cc.all, aes(x = time, 
                                   y = N, 
                                   linetype=`Diversification type`,
                                   group = `Diversification type`)) +
  geom_vline(xintercept = c(251.902,66),
             linetype = "dashed",
             alpha = 0.5) +   
  geom_line(aes(color=`Diversification type`,
                linetype = `Diversification type`), 
            size=1.5) +
  labs(x = "Time before present (Ma)",
       y = "Number of lineages") +
  theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=18),
        legend.position = c(.2,.8),
        legend.key.width=unit(2,"cm"),
        legend.box.background = element_rect(colour = "black")) +
  xlim(c(501, 0)) +
  scale_color_colorblind() +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed"))
  


geo.ltt.plot.new <- ggplot(cc.all, aes(x = time, 
                                   y = N, 
                                   linetype=`Diversification type`,
                                   group = `Diversification type`)) +
  geom_vline(xintercept = c(251.902,66),
             linetype = "dashed",
             alpha = 0.5) +   
  geom_line(aes(color=`Diversification type`,
                linetype = `Diversification type`), 
            size=1.5) +
  labs(x = "Time before present (Ma)",
       y = "Number of lineages") +
  theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=18),
        legend.position = "none",
        legend.key.width=unit(2,"cm"),
        legend.box.background = element_rect(colour = "black")) +
  xlim(c(501, 0)) +
  scale_color_manual(values = c("royalblue3","royalblue3","firebrick2","firebrick2", rep("gray", 100))) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed", rep("dotted", 100)))

  

# gggeo_scale(geo.ltt.plot)


png("geo_ltt.png", width=10,height = 9, units="in", res=500)
gggeo_scale(geo.ltt.plot.new)
dev.off()

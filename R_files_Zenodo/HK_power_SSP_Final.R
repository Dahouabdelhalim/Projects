#load packages
library(SSP)
library(tidyr)
library(ggplot2)
library(dplyr)
library(visreg)

#NOTE: 2 runs (2x response variables) at ~1hr per run

#####data management####
#load data
O2_primer<- read.csv("ReadyForPrimer_O2Consumption.csv")
feeding<- read.csv("ReadyForPrimer_Feeding.csv")

#convert all to positive
O2_primer$O2Consumption = O2_primer$O2Consumption + 1
feeding$FeedingRate = feeding$FeedingRate+ 1

######create O2 dataset#####
O2_primer$sector <- 1 
names(O2_primer)[names(O2_primer)=="tank.n"] <- "site"

O2_primer<- O2_primer %>% 
  mutate(sector = case_when(treatment == "Ambient" & species == "Chlorostoma" ~ "Ambient Chlorostoma", 
                            treatment == "Ambient" & species == "Lunella" ~ "Ambient Lunella",
                            treatment == "Elevated" & species == "Chlorostoma" ~ "Elevated Chlorostoma",
                            treatment == "Elevated" & species == "Lunella" ~ "Elevated Lunella"))

O2_power <- select(O2_primer, sector, site, O2Consumption)
O2_power$sector<- as.factor(O2_power$sector)
O2_power$extracol <- 1
str(O2_power)

######O2 power analysis#####

sectors <- levels(O2_power$sector)
sectors

#Defining arguments for simulation and sampling
N = 30
sites = 30
cases = 20
n = 20
m = 12
k = 10

#Lists to store results
sum.l <- opt.l <- qua.l <- vector(mode = "list", length = 3)

#Loop SSP at each sector
for (i in 1:length(sectors)){
  dat <- O2_power[O2_power$sector==sectors[i],2:length(O2_power)]
  #parameters for simulation
  par <- assempar(data = dat, type = "cover", Sest.method = "chao")
  # Simulation of data
  sim <- simdata(Par = par, cases = cases, N = N, sites = sites)
  # Quality of simulated data
  qua <- datquality(data = dat, dat.sim = sim, Par = par, transformation = "none", method = "euclidean")
  qua$sector <- rep(sectors[i], nrow(qua))
  qua.l[[i]] <- qua
  # Sampling and estimation of multse for each data set
  samp <- sampsd(dat.sim = sim, Par = par, transformation = "none",
                 method = "euclidean", n = n, m = m, k = k)
  # average of multse for each potential sampling design
  sum <- summary_ssp(samp, multi.site = TRUE)
  #Optimal sample sizes
  opt <- ioptimum(sum)
  opt <- as.data.frame(opt)
  opt$sv <- c("sites", "samples")
  opt <- pivot_longer(opt, cols = c("c1", "c2", "c3"), names_to = "cut", values_to = "effort")
  opt$sector <- rep(sectors[i], nrow(opt))
  opt.l[[i]] <- opt
  #arrangement to plot
  sum$sector <- rep(sectors[i], nrow(sum))
  sum.l[[i]] <- sum
}

beepr::beep() # makes a sound when done

#combine summary into a data frame
sum.df <- do.call(rbind.data.frame, sum.l)
sum.df$sector <- factor(sum.df$sector, levels = c("Ambient Chlorostoma", "Ambient Lunella", "Elevated Chlorostoma", "Elevated Lunella"))
#combine optimal sample sizes into a data frame
opt.df <- do.call(rbind.data.frame, opt.l)
opt.df$sector <- factor(opt.df$sector, levels = c("Ambient Chlorostoma", "Ambient Lunella", "Elevated Chlorostoma", "Elevated Lunella"))
#Combine quality features into a data frame
qua.df <- do.call(rbind.data.frame, qua.l)
# Generation of plot
my_breaks <- function(x) {
  y <- seq(min(x), max(x), 1)
  y <- round(y,0)
}

#Definition of values for shade areas
shade.opt <- opt.df %>%
  group_by(sector, sv) %>%
  filter(cut != "c1") %>%
  summarise(xmin = min(effort), xmax = max(effort))
shade.sub <- opt.df %>%
  group_by(sector, sv) %>%
  filter(cut != "c3") %>%
  summarise(xmin = min(effort), xmax = max(effort))

#####O2 plot####

#saving analyses results 
sum.df.o2 <- sum.df
write.csv(sum.df.o2, file = "sum.df.o2.csv")

shade.opt.o2 <- shade.opt
write.csv(shade.opt.o2, file = "shade.opt.o2 .csv")

shade.sub.o2 <- shade.sub
write.csv(shade.sub.o2, file = "shade.sub.o2 .csv")

#plot
sum.df.o2<- subset(sum.df.o2, sv == "sites") 
shade.opt.o2<- subset(shade.opt.o2, sv == "sites") 
shade.sub.o2<- subset(shade.sub.o2, sv == "sites") 

 fig.o2.power <- ggplot(sum.df.o2, aes(x=samples, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
  facet_wrap(~sector, scales = "free",
             ncol = 2,
             strip.position = "right")+
  theme_bw(base_size=16) +
  ylab ("Multivariate pseudo SE")+
  xlab("Sampling effort")+
  scale_y_continuous(breaks=seq(0.0, max(sum.df.o2$upper), 0.004))+
  scale_x_continuous(breaks = my_breaks)+
  theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2),
        strip.text.y = element_text(face = "italic"))+
  #geom_rect(data = shade.opt, aes_(x = NULL,y = NULL,
                                   #xmin=~xmin, xmax=~xmax+0.05, ymin=-Inf, ymax=Inf), alpha=0.5, fill="grey10")+
  #geom_rect(data = shade.sub, aes_(x = NULL,y = NULL,
                                   #xmin=~xmin, xmax=~xmax, ymin=-Inf, ymax=Inf), alpha=0.5, fill="grey50")+
  geom_text(aes_(x = ~samples, y = ~upper + 0.001, label = ~cum), na.rm = TRUE) +
  #geom_vline(data = shade.opt, aes(xintercept =xmax), linetype = 'dashed', color = '#80cdc1', size = 1)+ #max
  geom_vline(data = shade.opt.o2, aes(xintercept =xmin), color = '#80cdc1', size = 1, alpha = 0.8)+ #optimal
  geom_vline(data = shade.sub.o2, aes(xintercept =xmin), color = '#dfc27d', size = 1, linetype = "dashed", alpha = 0.8) #sub-optimal

fig.o2.power 

ggsave("fig.o2.power.jpeg",
       width = 15, height = 10)
   


######create feeding dataset#####
feeding$sector <- 1 
names(feeding)[names(feeding)=="tank.n"] <- "site"

feeding<- feeding %>% 
  mutate(sector = case_when(treatment == "Ambient" & species == "Chlorostoma" ~ "Ambient Chlorostoma", 
                            treatment == "Ambient" & species == "Lunella" ~ "Ambient Lunella",
                            treatment == "Elevated" & species == "Chlorostoma" ~ "Elevated Chlorostoma",
                            treatment == "Elevated" & species == "Lunella" ~ "Elevated Lunella"))

feeding <- select(feeding, sector, site, FeedingRate)
feeding$sector<- as.factor(feeding$sector)
feeding$extracol <- 1
str(feeding)

######Feeding power analysis#####

sectors <- levels(feeding$sector)
sectors

#Defining arguments for simulation and sampling
N = 30
sites = 30
cases = 20
n = 20
m = 12
k = 10

#Lists to store results
sum.l <- opt.l <- qua.l <- vector(mode = "list", length = 3)

#Loop SSP at each sector
for (i in 1:length(sectors)){
  dat <- feeding[feeding$sector==sectors[i],2:length(feeding)]
  #parameters for simulation
  par <- assempar(data = dat, type = "cover", Sest.method = "chao")
  # Simulation of data
  sim <- simdata(Par = par, cases = cases, N = N, sites = sites)
  # Quality of simulated data
  qua <- datquality(data = dat, dat.sim = sim, Par = par, transformation = "none", method = "euclidean")
  qua$sector <- rep(sectors[i], nrow(qua))
  qua.l[[i]] <- qua
  # Sampling and estimation of multse for each data set
  samp <- sampsd(dat.sim = sim, Par = par, transformation = "none",
                 method = "euclidean", n = n, m = m, k = k)
  # average of multse for each potential sampling design
  sum <- summary_ssp(samp, multi.site = TRUE)
  #Optimal sample sizes
  opt <- ioptimum(sum)
  opt <- as.data.frame(opt)
  opt$sv <- c("sites", "samples")
  opt <- pivot_longer(opt, cols = c("c1", "c2", "c3"), names_to = "cut", values_to = "effort")
  opt$sector <- rep(sectors[i], nrow(opt))
  opt.l[[i]] <- opt
  #arrangement to plot
  sum$sector <- rep(sectors[i], nrow(sum))
  sum.l[[i]] <- sum
}

beepr::beep() # makes a sound when done

#combine summary into a data frame
sum.df <- do.call(rbind.data.frame, sum.l)
sum.df$sector <- factor(sum.df$sector, levels = c("Ambient Chlorostoma", "Ambient Lunella", "Elevated Chlorostoma", "Elevated Lunella"))
#combine optimal sample sizes into a data frame
opt.df <- do.call(rbind.data.frame, opt.l)
opt.df$sector <- factor(opt.df$sector, levels = c("Ambient Chlorostoma", "Ambient Lunella", "Elevated Chlorostoma", "Elevated Lunella"))
#Combine quality features into a data frame
qua.df <- do.call(rbind.data.frame, qua.l)
# Generation of plot
my_breaks <- function(x) {
  y <- seq(min(x), max(x), 1)
  y <- round(y,0)
}

#Definition of values for shade areas
shade.opt <- opt.df %>%
  group_by(sector, sv) %>%
  filter(cut != "c1") %>%
  summarise(xmin = min(effort), xmax = max(effort))
shade.sub <- opt.df %>%
  group_by(sector, sv) %>%
  filter(cut != "c3") %>%
  summarise(xmin = min(effort), xmax = max(effort))

#####Feeding plot#####

#saving analyses results 
sum.df.feeding <- sum.df
write.csv(sum.df.feeding, file= "sum.df.feeding.csv")

shade.opt.feeding <- shade.opt
write.csv(shade.opt.feeding, file = "shade.opt.feeding.csv")

shade.sub.feeding <- shade.sub
write.csv(shade.sub.feeding, file = "shade.sub.feeding.csv")

#plot
sum.df.feeding <- subset(sum.df.feeding, sv == "sites") 
shade.opt.feeding <- subset(shade.opt.feeding, sv == "sites") 
shade.sub.feeding <- subset(shade.sub.feeding, sv == "sites") 

fig.feeding.power <- ggplot(sum.df.feeding, aes(x=samples, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
  facet_wrap(~sector, 
             scales = "free",
             ncol = 2,
             strip.position = "right")+
  theme_bw(base_size=16) +
  ylab ("Multivariate pseudo SE")+
  xlab("Sampling effort")+
  scale_y_continuous(breaks=seq(0.0, max(sum.df.feeding$upper), 0.001))+
  scale_x_continuous(breaks = my_breaks)+
  theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2),
        strip.text.y = element_text(face= "italic") )+
  #geom_rect(data = shade.opt, aes_(x = NULL,y = NULL,
            #xmin=~xmin, xmax=~xmax+0.05, ymin=-Inf, ymax=Inf), alpha=0.5, fill="grey10")+
  #geom_rect(data = shade.sub, aes_(x = NULL,y = NULL,
            #xmin=~xmin, xmax=~xmax, ymin=-Inf, ymax=Inf), alpha=0.5, fill="grey50")+
  geom_text(aes_(x = ~samples, y = ~upper + 0.001, label = ~cum), na.rm = TRUE) +
  #geom_vline(data = shade.opt, aes(xintercept =xmax), linetype = 'dashed', color = '#80cdc1', size = 1)+ #max
  geom_vline(data = shade.opt.feeding, aes(xintercept =xmin), color = '#80cdc1', size = 0.7)+ #optimal
  geom_vline(data = shade.sub.feeding, aes(xintercept =xmin), linetype = 'dashed', color = '#dfc27d', size = 1) #sub-optimal


fig.feeding.power
ggsave("fig.feeding.power.jpeg",
       width = 15, height = 10)

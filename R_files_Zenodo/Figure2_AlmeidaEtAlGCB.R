#Comparison of rarefaction-based models
#IMPORTANT: Set wd() to wd containing output from simulation analyses

#tidyverse for data manipulation and plotting
library(tidyverse)
#Cairo for exporting high quality figures
library(Cairo)

############ Upload results from separate analyses----
#Individual based rarefaction
IBR_richness <- read.csv("IBR_TotalRichness_AlmeidaEtAlGCB.csv")
IBR_extinctions <- read.csv("IBR_TotalExtinction_AlmeidaEtAlGCB.csv")
IBR_Pe <- read.csv("IBR_TotalPe_AlmeidaEtAlGCB.csv")
IBR_CI <- read.csv("IBR_TotalNull95CI_AlmeidaEtAlGCB.csv")

#Spatial rarefaction
SCR_richness <- read.csv("SCR_TotalRichness_AlmeidaEtAlGCB.csv")
SCR_extinctions <- read.csv("SCR_TotalExtinction_AlmeidaEtAlGCB.csv")
SCR_Pe <- read.csv("SCR_TotalPe_AlmeidaEtAlGCB.csv")
SCR_CI <- read.csv("SCR_TotalNull95CI_AlmeidaEtAlGCB.csv")

#Sample based rarefaction
SBR_richness <- read.csv("SBR_TotalRichness_AlmeidaEtAlGCB.csv")
SBR_extinctions <- read.csv("SBR_TotalExtinction_AlmeidaEtAlGCB.csv")
SBR_ENS <- read.csv("SBR_TotalENS_AlmeidaEtAlGCB.csv")
SBR_Pe <- read.csv("SBR_TotalPe_AlmeidaEtAlGCB.csv")
SBR_CI <- read.csv("SBR_TotalNull95CI_AlmeidaEtAlGCB.csv")
############ Calculate Pre-disturbance richness ----
PreS <- filter(IBR_richness, survey == "Pre", treatment == "Treatment")
PreS_mean <- summarise(group_by(PreS, scale), 
                       mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(10)), 
                       upper = mean(richness) + 1.96*(sd(richness)/sqrt(10)))
PreS_mean$group <- "Pre-disturbance"

############ Plot Figure 2A (Richness) ----
#Let's compare richness predictions across rarefaction methods

### IBR richness
#Extract null richness data without colonists
IBR_S <- filter(IBR_richness, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
IBR_S_mean <- summarise(group_by(IBR_S, scale), 
                        mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(10)), 
                        upper = mean(richness) + 1.96*(sd(richness)/sqrt(10)))
IBR_S_mean$group <- "Individual based rarefaction"
IBR_S_CIs <- filter(IBR_CI, site == "Region", colonists == 0, survey  == "Post1")
IBR_S_mean[2,]$lower <- IBR_S_CIs$lower
IBR_S_mean[2,]$upper <- IBR_S_CIs$upper

### SBR richness
#Extract null richness data without colonists
SBR_S <- filter(SBR_richness, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
SBR_S_mean <- summarise(group_by(SBR_S, scale), 
                        mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(10)), 
                        upper = mean(richness) + 1.96*(sd(richness)/sqrt(10)))
SBR_S_mean$group <- "Sampled based rarefaction"
SBR_S_CIs <- filter(SBR_CI, site == "Region", colonists == 0, survey  == "Post1")
SBR_S_mean[2,]$lower <- SBR_S_CIs$lower
SBR_S_mean[2,]$upper <- SBR_S_CIs$upper


### SCR richness
#Extract null richness data without colonists
SCR_S <- filter(SCR_richness, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
SCR_S_mean <- summarise(group_by(SCR_S, scale), 
                        mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(10)), 
                        upper = mean(richness) + 1.96*(sd(richness)/sqrt(10)))
SCR_S_mean$group <- "Spatially-clustered rarefaction"
SCR_S_CIs <- filter(SCR_CI, site == "Region", colonists == 0, survey  == "Post1")
SCR_S_mean[2,]$lower <- SCR_S_CIs$lower
SCR_S_mean[2,]$upper <- SCR_S_CIs$upper

Rarefaction_richness <- rbind(IBR_S_mean, SBR_S_mean, SCR_S_mean)

Fig2A <- ggplot(Rarefaction_richness, aes(x = scale, y = mean, color = group, group = group, shape = group, linetype = group)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width  = 0.05) +
  geom_line(size = 1.1) + 
  geom_point(size = 2.7) +
  ylab("Species richness") + #label y axis
  xlab(element_blank())+ #remove label from x axis
  theme_bw() + #remove greyscale background
  theme(legend.title = element_blank()) + #remove legend title
  theme(panel.grid = element_blank()) + #remove gridlines
  ylim(c(0,60))

#Check out Figure 2A
Fig2A

# Use Cairo package to export high quality figure

# tiff(filename="Fig2A_AlmeidaEtAlGCB.tiff",
#      type="cairo",
#      units = 'mm',
#      width=150,
#      height=100,
#      pointsize=12,
#      res = 600)

# Fig2A

# dev.off()


############ Plot Figure 2B (Extinctions) ----
### IBR extinctions
#Extract null extinctions data without colonists
IBR_extinctions <- filter(IBR_extinctions, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
IBR_extinctions_mean <- summarise(group_by(IBR_extinctions, scale), 
                        mean = mean(extinctions), lower = mean(extinctions) - 1.96*(sd(extinctions)/sqrt(10)), 
                        upper = mean(extinctions) + 1.96*(sd(extinctions)/sqrt(10)))
IBR_extinctions_mean$group <- "Individual based rarefaction"
IBR_extinctions_CIs <- filter(IBR_CI, site == "Region", colonists == 0, survey  == "Post1")
IBR_extinctions_mean[2,]$lower <- PreS_mean$mean[2] - IBR_extinctions_CIs$upper
IBR_extinctions_mean[2,]$upper <- PreS_mean$mean[2] - IBR_extinctions_CIs$lower

### SBR extinctions
#Extract null extinctions data without colonists
SBR_extinctions <- filter(SBR_extinctions, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
SBR_extinctions_mean <- summarise(group_by(SBR_extinctions, scale), 
                        mean = mean(extinctions), lower = mean(extinctions) - 1.96*(sd(extinctions)/sqrt(10)), 
                        upper = mean(extinctions) + 1.96*(sd(extinctions)/sqrt(10)))
SBR_extinctions_mean$group <- "Sampled based rarefaction"
SBR_extinctions_CIs <- filter(SBR_CI, site == "Region", colonists == 0, survey  == "Post1")
SBR_extinctions_mean[2,]$lower <- PreS_mean$mean[2] - SBR_extinctions_CIs$upper
SBR_extinctions_mean[2,]$upper <- PreS_mean$mean[2] - SBR_extinctions_CIs$lower


### SCR extinctions
#Extract null extinctions data without colonists
SCR_extinctions <- filter(SCR_extinctions, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
SCR_extinctions_mean <- summarise(group_by(SCR_extinctions, scale), 
                         mean = mean(extinctions), lower = mean(extinctions) - 1.96*(sd(extinctions)/sqrt(10)), 
                         upper = mean(extinctions) + 1.96*(sd(extinctions)/sqrt(10)))
SCR_extinctions_mean$group <- "Spatially-clustered rarefaction"
SCR_extinctions_CIs <- filter(SCR_CI, site == "Region", colonists == 0, survey  == "Post1")
SCR_extinctions_mean[2,]$lower <- PreS_mean$mean[2] - SCR_extinctions_CIs$upper
SCR_extinctions_mean[2,]$upper <- PreS_mean$mean[2] - SCR_extinctions_CIs$lower

Rarefaction_extinctions <- rbind(IBR_extinctions_mean, SBR_extinctions_mean, SCR_extinctions_mean)

Fig2B <- ggplot(Rarefaction_extinctions, aes(x = scale, y = mean, color = group, group = group, shape = group, linetype = group)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width  = 0.05) +
  geom_line(size = 1.1) + 
  geom_point(size = 2.7) +
  ylab("Species extinctions") + #label y axis
  xlab(element_blank())+ #remove label from x axis
  theme_bw() + #remove greyscale background
  theme(legend.title = element_blank()) + #remove legend title
  theme(panel.grid = element_blank()) + #remove gridlines
  ylim(c(0,30))

#Check out Fig2B
Fig2B

# Use Cairo package to export high quality figure

# tiff(filename="Fig2B_AlmeidaEtAlGCB.tiff",
#      type="cairo",
#      units = 'mm',
#      width=150,
#      height=100,
#      pointsize=12,
#      res = 600)

# Fig2B

#dev.off()

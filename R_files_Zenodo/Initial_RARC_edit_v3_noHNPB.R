## Code written by B Hogan 2022. This script generates cone-catches and color contrasts
# for Richard Childer's butterfly spectra, and exports xyz coordinates for plotting in matlab.

### load packages
library(pavo)
library(readxl)
library(tidyverse)
library(ggplot2)

### deal with sensitivity data

sens <- readxl::read_xlsx('../dat/Pavo jalmenus evagoras sensitivities.xlsx', skip = 1) %>%
  filter(complete.cases(uv)) %>%
  mutate(wl = as.numeric(wl)) %>%
  as.rspec()

# see sensitivities
plot(sens)

# okay, so they are organized to max 1. Best to integrate them to sum to 1
plot(pavo::sensdata(visual = "bluetit")) # as in the bluetit here
colSums(pavo::sensdata(visual = "bluetit"))

# Ah, well, do we first integrate and then cut at 360 or inverse?
# I think we integrate, then abridge sensitivities; otherwise we are effectively changing the 
# shape of the UV curve I think; shifting a lot of weight above the 360nm point.

# do that integration. 
sens[,2:dim(sens)[2]] <- sweep(sens[,2:dim(sens)[2]],2,colSums(sens[,2:dim(sens)[2]]),`/`)

# each column should sum to 1
colSums(sens)

# plot
plot(sens)

# and now, we've found out we're stuck with a 360nm lower limit in spectra, so must cut sensitivities too
sens <- sens %>%
  as.rspec(lim = c(360, 700)) 

plot(sens) # plot result

#### deal with spectral data, ah these are in a long format not great for use in pavo
# and 360nm + 

# read data, tidy a bit, make wide format
spec_data <- read.csv('../dat/7-27-20_new_spectra_smooth_R_ready_compiled_with_robots_standardillumination_only_noHNPB.csv') %>%
  pivot_wider(names_from = wavelength, values_from = Reflectance) %>%
  mutate(Wing = na_if(Wing, "Not_applicable")) %>%
  mutate(Side = na_if(Side, "Not_applicable")) %>%
  mutate(Agg_by = paste0(Individual, '_', Wing)) # for simplicity, make grouping variable for spectral averaging

# get index of start and end of spectral data (column called '360' and '810')
ss <- which(names(spec_data) == '360')
se <- which(names(spec_data) == '810')

# get spectral data separated out from metadata
specs <- data.frame(t(spec_data[,ss:se])) %>% # this is poor code since the column to stary
  rownames_to_column("wl") %>%
  mutate(wl = as.numeric(wl)) %>%
  as.rspec(lim = c(360, 700)) # cut to 700nm max

# zero out negative values
specs <- procspec(specs, fixneg = "zero") 
# aggregate/mean across side measurements
specs <- aggspec(specs, by = spec_data$Agg_by) 

# aggregate the metadata too
dat <- spec_data %>% 
  group_by(Agg_by) %>%
  summarise(Sex = first(Sex), 
            Specimen_or_robot = first(Specimen_or_robot),
            Individual = first(Individual),
            Wing = first(Wing))%>%
  arrange(match(Agg_by, names(specs[2:dim(specs)[2]]))) %>% # arrange in the order of spectral columns
  mutate(spec_name = names(specs[2:dim(specs)[2]])) # double check that the order of metadata corresponds with spectra, spec_name and Agg_by should line up

# apply visual sensitivities to get cone catches
# Assumption: since we are interested in the colors themselves, more than their appearance in a 
# given light environment we use ideal background and illuminant. It might make sense to re-consider this
# assumption. We also are currently doing a strictly color based analysis - basically saying 
# that we are not interested in brightness, which may or may not be true.

vis <- vismodel(specs, visual = sens, illum = 'ideal', bkg = 'ideal', relative = FALSE)
# so this is usml weights for each spectrum, not-normalized to sum to 1 across all cones, which
# comes in a next step. These usml values essentially encode both brightness (absolute usml values, or more
# likely the absolute value of the l channel, since brightness pathways usually assumed to be the longest cone type) 
# -and- color (relative usml). 

###### Get color coordinates, example tetrahedral plot

# get coordinates for colors in a tetrahedral space (and calculate hue, saturation, achieved saturation, and relative cone catches)
col <- colspace(vis) %>% # warnings expected (quantum catch, could not find columns)
  rownames_to_column('spec_name') %>%
  left_join(dat) 

# show a simple tetrahedral plot
plot(col)

# # here's an example, filtering out robots and coloring the dots by sex
# p <- col %>% 
#   filter(Specimen_or_robot == 'Specimen') %>%
#   mutate(Sex = factor(Sex)) %>%
#   mutate(Color = case_when(Sex == "M" ~ "Blue", Sex == "F" ~ "Red"))
# 
# tetraplot(p, labels = T, theta = 0, phi = 10, col = p$Color)

# save out a copy of col to send Richard 
# write.csv(col, file = "../out/Childers_butterfly_colspace.csv", row.names = F)

#### Generate some noise-weighted euclidean distances, these are one potential measure of color distance
# generally, given some assumptions a pair-wise distance above a threshold, typically 1 or 3, are considered
# discriminable.

# Following reference "True UV color vision in a butterfly with two UV opsins", Finkbeiner & Briscoe 2020
# https://www.biorxiv.org/content/10.1101/2020.11.14.382507v1
# "Parameters for the butterfly visual models were as follows: 
# Weber fraction=0.05 (Koshitaka et al. 2008) and relative abundances of photoreceptors, 
# V=0.13, B=0.2, L=1 (male) or UV=0.09, V=0.07, B=0.17, L=1 (females) (McCulloch et al. 2016a)."

# Now, we should reiterate that assumed values for cone densities are problematic, as noted by Gary Bernard,
# such that we're not hanging our coat on the absolute values for discriminability, but using them as
# a cue to the relative similarity between male and female, and robot to real colors

# Hm, note that the authors state that "The red receptor found in both sexes is not included in the
# visual modeling because their relative abundance is unknown." Well, given that and the 
# fact that values are only from H. erata, it may be best to steer clear of this JND approach.

# That said, in https://academic.oup.com/mbe/article/34/9/2271/3827455?login=false, the authors seem to include
# the longest wavelength? Apparently with these same densities. However, there they also have values 
# for G = 1 and L = 1 depending on erato sex (or a typo). Are our sensitivities strictly homologous to 
# these, given diff families and everything? 

# I am not satisfied that this can be reliably done without a lot of extra work and thought.
# A proper lit search might turn up values that the authors are happy with. For now I'll leave this 
# orphaned here, in case that happens. If so, values in jnd_dist below can be swapped in easily enough.

weber_fraction <- 0.05
cone_densities <- c(0.09, 0.07, 0.17, 1) 

# pavo relative densities are formatted such that rarest is set to 1
cone_densities <- cone_densities / min(cone_densities) 
  
# generate weighted euclidean distances
jnd_dist <- coldist(vis, weber = 0.05, n = cone_densities)

#### Instead, why not work with euclidean distances.

euc_dist <- coldist(colspace(vis))

# This results in a unitless measure of euclidean distance between each spectrum in tetrahedral colorspace
# in general, colors that are more similar should have smaller euclidean distances in this space. The JND stuff
# above only introduces modelled non-linearity in similarity, and provides a cut-off for colors that are 
# not considered discriminable. If we only need to say A is closer to B than is C, euclidean distance comes 
# without the headaches and assumptions of the modelling.

# I would say that the X,Y,Z coordinates and hue (theta, phi) and saturation (r.achieved) for each spectrum found in col,
# alongside the euclidean distances should be the data that Richard requires for any required statistical analysis. 
# X, Y, and Z positions or hue (note these vals are circular) and saturation (note 0-1)
# might be swapped into his multivariate analyses (swapping for PCA axes). Hue and sat, or XYZ both
# define the position of the spectrum in terms of color (omitting brightness). I lean toward XYZ here I think for
# numerical simplicity, and ease of relation to a tetrahedral figure.

# In addition some summary of distances between
# and within different groups would be valid. In order to summarize across replicates, I think it is perfectly
# legitamate to average XYZ coordinates, as long as a visual check shows that the averaging is reasonable (i.e. that 
# the average is not some spurious color)

#### Get distribution of euclidean distances between model spectra and all specimen spectra

# list of model spec_name 
models <- col %>% 
  filter(Specimen_or_robot == 'Robot') %>% 
  select(spec_name)

# get the relevant values out (where one patch is a model, and one a specimen, but not both)
a <- euc_dist %>% 
  filter(xor(patch1 %in% unlist(models), patch2 %in% unlist(models)))

# just double-check our selection
any(a$patch1 %in% unlist(models)) # make sure that no models are listed under patch1, looking for false here
all(a$patch2 %in% unlist(models)) # make sure that ALL patch2 correspond to models, looking for true here

# get a lookup for the model names
models_names <- col %>% 
  filter(Specimen_or_robot == 'Robot') %>% 
  select(Individual, spec_name) %>%
  rename(patch2 = spec_name)

# make a violinplot of each distance, and order by distance too
# Richard suggested facet by wing, and color by sex



a_dat<-a %>% 
  left_join(models_names) %>%
  left_join(dat, by = c('patch1'='spec_name')) 

a_dat1<-a_dat %>% 
  mutate(Individual = fct_reorder(Individual.x, dS))


violin_1=ggplot(
  data = a_dat1
  , mapping = aes(
    y =dS
    , x=Sex
    , colour=factor(Individual)
  )
)+
  labs(
    x = 'Specimen Sex'
    , y = 'Distance from specimens (dS)'
  )+
  geom_point(position=position_jitterdodge(dodge.width = .9),
             mapping = aes(x=Sex,
                           colour =factor(Individual)
             )              , size=1,alpha=0.25)+
  geom_violin(alpha=.4, draw_quantiles = c(0.25, 0.5, 0.75),
              mapping=aes(x=Sex,fill=factor(Individual),colour=factor(Individual)))+
  #  geom_smooth(method="loess",alpha=.4, se =TRUE)+
  facet_grid(.~.)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+
  scale_fill_manual(values=c("#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+theme_bw()
violin_1_2

violin_1_3=violin_1_2

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=8, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 8, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=8, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=8, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 8, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7

# 
# a %>% 
#   left_join(models_names) %>%
#   left_join(dat, by = c('patch1'='spec_name')) %>%
#   mutate(Individual = fct_reorder(Individual.x, dS)) %>%
#   ggplot(aes(x = Individual, y = dS, color = Sex)) +
#   geom_violin() +
#   geom_jitter(width = 0.15, size = 1) + 
#   stat_summary(fun.data=mean_sdl,
#        geom="pointrange", color="red") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_grid(vars(Wing))

# for me, this plot (of distributions of unweighted euclidean distances between specimen spectra, and model spectra)
# tells the story that HNPB are most similar to the specimens in terms of color. If we did find values we are happy
# with for JND type analysis, we can simply swap in weighted euclidean distances for unweighted ones. Format in any
# way that makes sense to you.

# ggsave('../out/distance_violin.png', plot = last_plot(), units = 'cm', scale = 1.75, width = 10, height = 10)

# what remains then is a nice plot of the locations of specimens in the tetracolorspace. R/pavo has poor plotting 
# functions in my opinion for this, so I think I'll use matlab, on a saved datasheet from here. I think we'll
# want to plot only specimens

plot_col <- col %>%
  filter(Specimen_or_robot == 'Specimen') 

# write.csv(plot_col, '../out/interim_xyz_coords.csv', row.names = F)

#####Make a dot plot showing USML stimulation values of wings versus robot models
usml<-read.csv('../dat/Childers_butterfly_colspace_usml_RARC_edit_noHNPB.csv')
usml_long<- gather(usml, opsin, stimulation, u:l, factor_key=TRUE)
usml_long$Specimen_robot<-factor(usml_long$Specimen_robot
                , levels = c("F_FW", "M_FW", "F_HW","M_HW",
                             "3M_Blue_R374_depol","3M_Blue_R374_pol",
                             "3M_Yellow_R312_depol","3M_Yellow_R312_pol","3M_Red_G280_depol","3M_Red_G280_pol"))

usml_long$Category<-factor(usml_long$Category
                                 , levels = c("F_FW", "M_FW", "F_HW","M_HW","Model"))
violin_1=ggplot(
  data = usml_long
  , mapping = aes(
    y =stimulation
    , x=Category
    , colour=factor(Specimen_robot)
  )
)+
  labs(
    x = 'Category'
    , y = 'opsin stimulation'
  )+
  # geom_point(position=position_jitterdodge(dodge.width = .9),
  #            mapping = aes(x=Category,
  #                          colour =factor(Specimen_robot)
  #            )              , size=1.8,alpha=0.25)  +
  stat_summary(fun.data = "mean_se", position=position_jitterdodge(dodge.width = .8),aes(colour = factor(Specimen_robot)), size = 0.9)+
  stat_summary(fun.data = "mean_se",geom='errorbar', aes(colour = factor(Specimen_robot)), size = 1.5)+
  # geom_violin(alpha=.4, draw_quantiles = c(0.25, 0.5, 0.75),
  #             mapping=aes(x=Sex,fill=factor(Individual),colour=factor(Individual)))+
  #  geom_smooth(method="loess",alpha=.4, se =TRUE)+
  facet_grid(.~opsin)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#fdae61","#2c7bb6","#fdae61","#2c7bb6","#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+
  scale_fill_manual(values=c("#fdae61","#2c7bb6","#fdae61","#2c7bb6","#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+theme_bw()
violin_1_2


violin_1_3=violin_1_2

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 16, face = "bold"))+
  theme(axis.text.x=element_text(colour="Black", size=16, face="plain",angle=90,hjust=0))+
  theme(axis.title=element_text(colour="Black", size=18, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 16, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             # panel.grid.major = element_blank(),
                             # panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
# theme(strip.background=element_blank())
violin_1_7


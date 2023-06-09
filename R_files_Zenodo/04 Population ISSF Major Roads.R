
# Code orignally from:
#' @article{doi:10.1111/1365-2656.13087,
#'   author = {Muff, Stefanie and Signer, Johannes and Fieberg, John},
#'   title = {Accounting for individual-specific variation in habitat-selection studies: Efficient estimation of mixed-effects models using Bayesian or frequentist computation},
#'   journal = {Journal of Animal Ecology},
#'   volume = {89},
#'   number = {1},
#'   pages = {80-92},
#'   keywords = {conditional logistic regression, glmmTMB, integrated nested Laplace approximations (INLA), multinomial regression, random effects, resource-selection functions, step-selection functions},
#'   doi = {10.1111/1365-2656.13087},
#'   url = {https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13087},
#'   eprint = {https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2656.13087},
#'   abstract = {Abstract Popular frameworks for studying habitat selection include resource-selection functions (RSFs) and step-selection functions (SSFs), estimated using logistic and conditional logistic regression, respectively. Both frameworks compare environmental covariates associated with locations animals visit with environmental covariates at a set of locations assumed available to the animals. Conceptually, slopes that vary by individual, that is, random coefficient models, could be used to accommodate inter-individual heterogeneity with either approach. While fitting such models for RSFs is possible with standard software for generalized linear mixed-effects models (GLMMs), straightforward and efficient one-step procedures for fitting SSFs with random coefficients are currently lacking. To close this gap, we take advantage of the fact that the conditional logistic regression model (i.e. the SSF) is likelihood-equivalent to a Poisson model with stratum-specific fixed intercepts. By interpreting the intercepts as a random effect with a large (fixed) variance, inference for random-slope models becomes feasible with standard Bayesian techniques, or with frequentist methods that allow one to fix the variance of a random effect. We compare this approach to other commonly applied alternatives, including models without random slopes and mixed conditional regression models fit using a two-step algorithm. Using data from mountain goats (Oreamnos americanus) and Eurasian otters (Lutra lutra), we illustrate that our models lead to valid and feasible inference. In addition, we conduct a simulation study to compare different estimation approaches for SSFs and to demonstrate the importance of including individual-specific slopes when estimating individual- and population-level habitat-selection parameters. By providing coded examples using integrated nested Laplace approximations (INLA) and Template Model Builder (TMB) for Bayesian and frequentist analysis via the R packages R-INLA and glmmTMB, we hope to make efficient estimation of RSFs and SSFs with random effects accessible to anyone in the field. SSFs with individual-specific coefficients are particularly attractive since they can provide insights into movement and habitat-selection processes at fine-spatial and temporal scales, but these models had previously been very challenging to fit.},
#'   year = {2020}
#' }
#' 

# Edited by Benjamin Michael Marshall, Samantha Nicole Smith and Max Dolton Jones 2021-07-14 for publication of: 
# title = {How do King Cobras move across a major highway? Unintentional wildlife crossing structures may facilitate movement.}
## authors= {# Max Dolton Jones, Benjamin Michael Marshall, Samantha Nicole Smith, Matt Crane, InÃªs Silva, Taksin Artchawakom, Pongthep Suwanwaree, Surachit Waengsothorn, Wolfgang Wüster, Matt Goode, Colin Thomas Strine}
## year = {2021}


library(dplyr)
library(raster)
library(amt)
library(ggplot2)
install.packages("INLA", repos=c(getOption("repos"),
                                 INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

# Create Folders ----------------------------------------------------------

loc.data <- paste0("./DATA/")
loc.new <- paste0("./Figures/")
loc.shape <- paste0("./Shapefiles")

# Read in animal tracking data --------------------------------------------

# Replace file name with the data you are interested in.
# Either adult males or females during the breeding or non-breeding season
snake.data <- read_csv(file = paste0(loc.data, "ISSF_AM_breeding.csv"),
                     locale = locale(tz = "Asia/Bangkok"))

crs.proj <- CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

dat <- snake.data

# Load in raster covariates -----------------------------------------------

## reading inverted tifs from ISSF code
covariates <-  list(paste0(loc.shape, "/Major_Roads_dist_INVERT.tif"))

rp <- stack(x = covariates)

cov.names <- c("dist_major_road")
names(rp) <- cov.names 

# START OF IMPORTED (and edited) CODE from Muff et al --------------------------------------------------

# select only the information we need
dat <- dat %>% 
  select(x = x, y = y,
         t = datetime, id = id)


# separate the individuals out into a list-column of dataframes, each item an animal
dat_all <- dat %>% 
  nest(-id) 
# if you want to add meta data, like sex, do so here. We had a second dataframe
# that included metadata
dat_all <- left_join(dat_all,
                     snake.meta %>% 
                       select(id = deployment.id, sex = animal.sex))
# map operates a lot like an apply or loop. It repeats the function to each item
# in the list. In this case we make all the individual dataframes into track
# objects. Check dat_all to see the list of dataframes now has a second
# dataframe/tibble for the track.
dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, x, y, t, crs = crs.proj)
  }))

# Here the summarize_sampling_rate is repeated on each track object to give you
# an individual level summary.
dat_all %>% 
  mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  select(id, sr) %>%
  unnest(cols = c(sr))

# So dat_all is set up ready for the generation of random steps and the
# extraction of raster covariate values

# Step-Selection function -------------------------------------------------

# To make the random steps we use the map function again to do this on an
# individual basis. We calculate the step lengths so we have data to draw from
# to generate the random locations. Then we drop out the zero step lengths that
# break things, or we can drop the sub-GPS error points because we cannot be
# confident they were moves (ie not a new selection choice from the animal).
# Then the random steps, can be generous here because we aren't using GPS data
# with the high resolution that often entails. We have use more randoms to pick
# up smaller changes and maximise the high res rasters we have. Final step is to
# get the covariate values. After the mutate+map bit is a few lines to compile
# the data into a single v large dataframe for the model.
dat_ssf <- dat_all %>% 
  mutate(stps = map(trk, ~ .x %>%
                      steps() %>% 
                      filter(sl_>0) %>% # removing the non-moves, or under GPS error
                      random_steps(n = 200) %>% 
                      extract_covariates(rp))) %>% 
  select(id, stps) %>%
  unnest(cols = c(stps)) %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"),
    cos_ta = cos(ta_), 
    log_sl = log(sl_))

# So dat_ssf is the correct format for the modelling. We have the id, a heap of
# movement columns, and the critical case_ that describes whether they used a
# location or not.


# Running the INLA model --------------------------------------------------

# We can run the INLA model using the priors and set-up from Muff et al.
# Precision for the priors of slope coefficients
prec.beta.trls <- 1e-4

# "In the model formula for INLA, we set the stratum-specific intercept variance
# to $10^6$ (or rather: the precision to $10^{Ã¢Ë†â€™6}$) by fixing it (`fixed=T`)
# to an initial value. The other precision is given a PC(1,0.05) prior:"

###### For PYBI data
### making formulas for each habitat feature
 
# road 
formula.random.r <- y ~  -1 + 
  dist_major_road + # fixed covariate effect
  dist_major_road:log_sl + # covar iteractions
  dist_major_road:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, dist_major_road, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 




# Last two parts of the model are the gaussian processes to deal with step_id
# (each move), and the id of the animal. So steps_id is the covariate, then in
# the second one is id that is weighted by the raster cov.

# Fit the models -----------------------------------------------------------


# road model 
inla.ssf.r <- inla(formula.random.r, family = "Poisson", data = dat_ssf, 
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls, verbose = TRUE))
)



# Model results -----------------------------------------------------------

# We can see a dataframe of the fixed effects here, so each covar and any interaction terms
inla.ssf.r$summary.fixed



# "Since variances are parameterized and treated as precisions, the summary of
# the respective posterior distributions is given for the precisions:"
inla.ssf.r$summary.hyperpar




inla_mmarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.mmarginal(inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mode of variance"))
  results
}
inla_emarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mean of variance"))
  results
}

# "Posterior mean and mode are obtained as"

inla_emarginal(inla.ssf.r)
inla_mmarginal(inla.ssf.r)


fixed.df.r <- inla.ssf.r$summary.fixed




fixed.df.r <- inla.ssf.r$summary.fixed
fixed.df.r$term <- row.names(fixed.df.r)
names(fixed.df.r) <- c("mean", "sd", "q025", "q50", "q975",
                        "mode", "kld", "term")

write_csv(fixed.df.r, "Pop_ssf_mods.csv")



##### figures for population level ssf
library(readr)

# read in new datafiles that have new Season and Sex columns

#edited data for adult males
popssf.am <- read_csv(file = paste0(loc.data, 
                                       "Pop_ssf_mods_AM_all.csv"))

#edited data for adult females
popssf.af <- read_csv(file = paste0(loc.data, 
                                       "Pop_ssf_mods_AF_all.csv"))

#edited data for all snakes
popssf.all <- read_csv(file = paste0(loc.data, 
                                   "Pop_ssf_mods_all.csv"))

#### Distance_feature plot

colors1 <- c("#00aedb", "#f37735")

popssf.all %>% 
  filter(term %in% grep("log_sl|cos_ta", term, value = TRUE,
                        invert = TRUE)) %>% 
  mutate(term = factor(term, levels = c("dist_major_road"))) %>%
  arrange(term) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5) +
  geom_point(aes(x = term, y = mean, colour = season, group = season),
             size = 3, position = position_dodge(width = 1)) +
  geom_errorbar(aes(x = term, ymin = q025, ymax = q975,
                    colour = season, group = season), size = 1,
                position = position_dodge(width = 1)) +
  facet_wrap(.~sex) +
  scale_shape_manual(values = c(3, 16)) +
  scale_x_discrete(labels = c()
  ) +
  labs(x = "Distance to major roads", y = expression(beta), colour = "Season") +
  scale_colour_manual(values = c("#00aedb",
                                 "#f37735",
                                 labels = c("Breeding", "Non-breeding"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = 4, hjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_line(colour = "grey65", linetype = 1),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = 2),
        legend.text = element_text(lineheight = 1),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 1, face = 2),
        axis.title.x = element_text(hjust = 0.5, face = 2, margin = margin(10,0,0,0)),
        plot.title = element_text(face = 4),
        strip.text.y = element_blank()
  ) +
  guides(shape = guide_none())

ggsave(file = paste0("./Figures/Population season ISSF plot.png"), width = 200, height = 120,
       dpi = 600, units = "mm")
ggsave(file = paste0("./Figures/Population season ISSF plots.pdf"), width = 200, height = 120,
       units = "mm")

# end
#### Code for reproducing analysis presented in: Baker H.K., Bruggeman C.E.F., and Shurin J.B. (2022). Population niche width is driven by within-individual niche expansion and individual specialization in introduced Brook Trout in mountain lakes. Oecologia. ####

library(tidyverse)
library(betareg)
library(bbmle)
library(viridis)
library(patchwork)


######## Read in gut data and clean up ########
gut <- read_csv("gut_deposit.csv") # gut contents

# Remove Blue since it only has 1 fish
gut <- gut %>%
  filter(lake != "Blue")

### Check taxa names
unique(gut$Order)

### Remove digested or unidentifiable rows
unid <- c("Digestive.matter","Unknown.bits","Unknown.full")
gut$Order <- as.character(gut$Order)
gut <- gut %>% 
  filter(!Order %in% unid)

### Set factors as such
cols <- c("id","lake", "Order", "sample.date")
gut[cols] <- lapply(gut[cols],factor)
rm(cols)



### Sum across multiple records for same id-prey combination
gut <- gut %>% 
  group_by(id, lake, sample.date, Order) %>% 
  summarise(count = sum(count))


# make a list of taxa for later use
taxa <- gut$Order
lakes <- as.character(unique(gut$lake))

####### Spread gut data to wide format #######
gut.wide <- gut %>% 
  select(id, lake, sample.date, Order, count) %>% 
  spread(key = Order, value = count) %>%  # spread to wide
  mutate_if(is.numeric, ~replace_na(., 0))

############## Read in fish data and cleanup #############
## Read in
fish_data <- read_csv("fish_deposit.csv") # fish data


## Set factors as such
cols <- c("id","lake", "sample.date")
fish_data[cols] <- lapply(fish_data[cols],factor)
rm(cols)

## Fix lake names
levels(fish_data$lake)
fish_data <- fish_data %>% 
  mutate(lake = fct_recode(lake, 
                           "Gaylor 1" = "Gaylor 1 (Lower)",
                           "Gaylor 1" = "Gaylor1_lower",
                           "Gaylor 2" = "Gaylor 2 (Upper)"))

######## Read in isotope data ###########
iso <- read_csv("isotope_deposit.csv")
iso$id <- factor(iso$id) # set fish ID as factor

#### Merge isotope and fish data ######
fish_data <- left_join(x = fish_data, y = iso, by = "id")
fish_data$id <- factor(fish_data$id)

########### Merge gut and fish data ##################
gut.wide <- left_join(x = gut.wide, y = fish_data, by = c("id", "lake", "sample.date")) #wide
gut.tidy <- left_join(x = gut, y = fish_data, by = c("id", "lake", "sample.date"))

# Change prey columns to integers
gut.wide <- gut.wide %>% mutate_at(vars(Annelidae:Trombidiformes),funs(as.integer))
# Change id to integer
gut.wide$id <- as.integer(as.character(gut.wide$id))
# Change from tibble to dataframe. Very important--will not work if you don't do this.
gut.wide <- gut.wide %>% 
  select(id,lake,sample.date, length.sl.mm, weight.wet.g, delta13C, delta15N,everything()) %>% 
  select(-processed.date) 
gut.wide$lake <- as.character(gut.wide$lake) # make lake character vector
gut.wide <- as.data.frame(gut.wide) 

################ Now for Calculating Specialization and Niche Width Metrics #####################
library(RInSp)
######## Calculate E #############
set.seed(2)
## Function for calculating E (individual specialization metric) 
calc.E <- function(data,lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2:7), 
                    subset.rows = c("lake", lake))
  Emc(a, popd.type = "sum") 
} #end calc.E

calc.NODF <- function(data,lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2:7), 
                    subset.rows = c("lake", lake))
  NODF(a)
} # end calc.NODF 

calc.popd <- function(data, lake) {
  a <- import.RInSp(data, row.names = 1, info.cols = c(2:7), 
                    subset.rows = c("lake", lake))
  pop.diet(a, prop = "sum") 
} # end calc.popd

#### Edit WTdMC to return WIC ####
WTdMC.custom <- function (dataset, pop.diet = "sum", replicates = 999, print.ris = TRUE) 
{
  if (!inherits(dataset, "RInSp")) 
    stop("The input must be an object of class RInSp")
  if (dataset$data.type != "integer") 
    stop("Input data type must be integer.")
  if (pop.diet %in% c("sum", "average") == FALSE) 
    stop("The specified population diet type is wrong.")
  if (pop.diet == "sum") 
    diet.pop <- 0
  else diet.pop <- 1
  replicates <- as.integer(replicates)
  if (replicates <= 1) 
    stop("Wrong number of replicates.")
  if (print.ris == TRUE) 
    cat("\\n If your dataset is big, this can take time. Please be patient. \\n")
  if (!is.double(dataset$resources)) 
    dataset$resources <- matrix(as.double(dataset$resources), 
                                dataset$num.individuals, dataset$num.prey)
  if (!is.integer(replicates)) 
    replicates <- abs(as.integer(replicates))
  Ris <- .Call("CWTdMC", dataset$resources, as.vector(diet.pop), 
               as.vector(replicates), PACKAGE = "RInSp")
  attributes(Ris)$dimnames[[2]] <- c("Zero", "WIC", "BIC", 
                                     "TNW", "WonT")
  cum.distr <- ecdf(Ris[, 5])
  pvalue <- cum.distr(Ris[1, 5])
  checkZero <- (dataset$proportion == 1)
  if (sum(checkZero) > 0) 
    Zeros <- dataset$ind.names[rowSums(checkZero) * c(1:length(dataset$ind.names))]
  else Zeros <- Ris[1, 1]
  Ris2 <- list(WonT = Ris[1, 5], Zeros = Zeros, p.value = cum.distr(Ris[1, 
                                                                        5]), montecarlo = Ris, parameter = 5)
  class(Ris2) <- "RInSp"
  if (print.ris == TRUE) {
    cat("\\n Using Roughgarden's 1979 equations, based on Shannon-Weaver diversity index: ")
    cat("\\n Within-individual component          = ", Ris[1, 
                                                          2])
    cat("\\n Between-individual component         = ", Ris[1, 
                                                          3])
    cat("\\n Total Niche Width for the population = ", Ris[1, 
                                                          4])
    cat("\\n The value of WIC/TNW is: ", Ris[1, 5])
    cat("\\n The p-value is: ", pvalue, "\\n")
    if (Ris[1, 1] > 0) {
      if (as.integer(Ris[1, 1]) == 1) 
        cat("\\n Warning: ", as.integer(Ris[1, 1]), " individual out of your population of ", 
            as.integer(dataset$num.individuals))
      else cat("\\n Warning: ", as.integer(Ris[1, 1]), " individuals out of your population of ", 
               as.integer(dataset$num.individuals))
      cat("\\n have Shannon-Weaver scores equal to zero.  This may exaggerate")
      cat("\\n the apparent degree of individual specialization, because these")
      cat("\\n scores will drag the mean SWi score down and thus reduce WPC.")
      cat("\\n This can be particularly misleading when an individual specializes")
      cat("\\n entirely on one resource, which happens to be the most commonly")
      cat("\\n consumed resource. This individual will have a Shannon-Weaver score")
      cat("\\n of 0, even though its diet proportions are not unlike the proportions")
      cat("\\n of the population as a whole.\\n")
      cat("\\n")
      if (as.integer(Ris[1, 1]) == 1) 
        cat("\\n The name of the individual is: \\n")
      else cat("\\n The names of the individuals are: \\n")
      cat(Zeros, "\\n")
    }
  }
  return(c(Ris[1,],pvalue))
} # End WTdMC.custom 

calc.WICd <- function(data,lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2:7), 
                    subset.rows = c("lake", lake))
  WTdMC.custom(a, pop.diet = "sum") 
} #end calc.WIC.d

##### loop through lakes and generate a dataframe of lake and E
lakes <- unique(gut.wide$lake)

results <- as.data.frame(lakes) %>% 
  add_column(E = NA) %>% 
  add_column(meannullE = NA) %>% 
  add_column(Eadj = NA) %>% 
  add_column(p.value = NA) %>% 
  add_column(NODF = NA) %>% 
  add_column(richness = NA) %>% 
  add_column(Levin.D = NA) %>% 
  add_column(WIC.d = NA) %>% 
  add_column(BIC.d = NA) %>% 
  add_column(TNW.d = NA) %>% 
  add_column(WonT.d = NA) %>% 
  add_column(p.WonT.d = NA)


for(i in 1:length(lakes)){
  a <- calc.E(gut.wide, lakes[i])
  b <- calc.NODF(gut.wide, lakes[i])
  c <- calc.popd(gut.wide, lakes[i])
  d <- calc.WICd(gut.wide, lakes[i])
  results[i,2] <- a$E
  results[i,3] <- a$meannullE
  results[i,4] <- a$Eadj
  results[i,5] <- a$p.value
  results[i,6] <- b$NODF
  results[i,7] <- c$richness
  results[i,8] <- c$D
  results[i,9] <- d[2]
  results[i,10] <- d[3]
  results[i,11] <- d[4]
  results[i,12] <- d[5]
  results[i,13] <- d[6]
} #i

## Calculate median individual taxa richness
count1 <- gut.tidy %>% 
  group_by(id, lake) %>% 
  summarise(n = n())
count2 <- count1 %>% 
  group_by(lake) %>% 
  summarise(median.richness = median(n))
results <- left_join(results, count2, by = c("lakes" = "lake"))


############# Analysis of prey length data ###################
length.tidy <- read_csv("gut_length.csv")
length <- length.tidy %>% 
  group_by(id) %>% 
  mutate(x = seq(1:n())) %>% 
  pivot_wider(names_from = x, 
              values_from = length.mm,
              values_fill = list(length.mm = 0))

length$lake <- as.factor(as.character(length$lake))

# make df
length.df <- as.data.frame(length)

######## Calculate E continuous #############
#### Function for calculating E (individual specialization metric) 
calc.E.c <- function(data, lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2), 
                    subset.rows = c("lake", lakes[i]))
  Eindex(a)
} # end Calc.e.c

calc.NODF.c <- function(data, lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2), 
                    subset.rows = c("lake", lakes[i]))
  NODF(a)
} #end calc.NODF.c

calc.wont <- function(data, lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2), 
                    subset.rows = c("lake", lakes[i]))
  WTcMC(a, weight = "N_items", replicates = 50000)
} #end calc.NODF.c


##### loop through lakes and generate a dataframe of lake and E ######
results1 <- results %>% 
  add_column(E.c = NA) %>% 
  add_column(NODF.c = NA) %>% 
  add_column(WonT = NA) %>% 
  add_column(p.WonT = NA) %>% 
  add_column(TNW = NA)

set.seed(23)
for(i in 1:length(lakes)){
  a <- calc.E.c(length.df, lakes[i])
  b <- calc.NODF.c(length.df, lakes[i])
  c <- calc.wont(length.df, lakes[i])
  results1[i,15]  <-  a$E
  results1[i,16] <- b$NODF
  results1[i,17] <- c$WonT
  results1[i,18] <- c$p.value
} # continuous 

## Variance of length from each lake
length.stats <- length.tidy %>% 
  group_by(lake) %>% 
  summarise(var.length = var(length.mm), 
            min.length = min(length.mm), 
            max.length = max(length.mm)) %>% 
  mutate(range.length = max.length - min.length)

results2 <- left_join(results1, length.stats, by = c("lakes" = "lake")) %>% 
  mutate(TNW = var.length)

### Function for extracting null values from MCMC ###
calc.wont.new <- function(data, lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2), 
                    subset.rows = c("lake", lakes[i]))
  WTcMC(a, weight = "N_items", replicates = 50000)
} #end calc.wont.new

##### 9/30/2021 loop through lakes to add WIC, BIC, NULL values for all ######
results2.new <- results2 %>% 
  add_column(null.WIC = NA) %>% 
  add_column(WIC = NA) %>% 
  add_column(null.BIC = NA) %>% 
  add_column(BIC = NA) %>% 
  add_column(null.WonT = NA) 


set.seed(23)
for(i in 1:length(lakes)){
  a <- calc.wont.new(length.df, lakes[i])
  results2.new[i,24]  <-  mean(a$montecarlo[2:nrow(a$montecarlo),1]) # null wic
  results2.new[i,25] <- a$montecarlo[1,1] #real WIC
  results2.new[i,26] <- mean(a$montecarlo[2:nrow(a$montecarlo),2]) # null BIC
  results2.new[i,27] <- a$montecarlo[1,2] #real BIC
  results2.new[i,28] <- mean(a$montecarlo[2:nrow(a$montecarlo),4]) # null WonT
} # continuous 

results2.new <- results2.new %>% 
  mutate(null.BIC.TNW = 1 - null.WonT)

#### 9/30/2021 loop through lakes to add NULL values for discrete data ######
results2.new2 <- results2.new %>% 
  add_column(null.WIC.d = NA) %>% 
  add_column(null.BIC.d = NA) %>% 
  add_column(null.WonT.d = NA) 

calc.WICd.new <- function(data,lake){
  a <- import.RInSp(data, row.names = 1, info.cols = c(2:7), 
                    subset.rows = c("lake", lakes[i]))
  WTdMC(a, pop.diet = "sum", replicates = 50000) 
} #end calcWICd.nw

set.seed(23)
for(i in 1:length(lakes)){
  a <- calc.WICd.new(gut.wide, lakes[i])
  results2.new2[i,30]  <-  mean(a$montecarlo[2:nrow(a$montecarlo),2]) # null wic
  results2.new2[i,31] <- mean(a$montecarlo[2:nrow(a$montecarlo),3]) # null BIC
  results2.new2[i,32] <- mean(a$montecarlo[2:nrow(a$montecarlo),5]) # null WonT
} # continuous 

results2.new2 <- results2.new2 %>% 
  mutate(null.BIC.TNW.d = 1 - null.WonT)


################################### Isotope Ellipses ####################################
library(SIBER)
set.seed(1)
iso <- gut.wide %>% 
  select(id, lake, delta13C, delta15N) 
lake.id <- tibble(lake = unique(iso$lake), community = factor(seq(1,13)))
iso1 <- left_join(iso, lake.id, by = "lake")
iso2 <- iso1 %>%
  select(delta13C, delta15N, lake, community) %>% 
  rename(iso1 = delta13C, iso2 = delta15N, group = lake, community = community) %>% 
  drop_na()

siber1 <- createSiberObject(iso2)
group.ML <- groupMetricsML(siber1)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber1, parms, priors)
SEA.B <- siberEllipses(ellipses.posterior)
SEA.B.modes <- lapply( #extract modes 
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)
SEA.B.modes1 <- data.frame(SEA.B.modes) %>% 
  pivot_longer(1:13)


siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

iso.results <- data.frame(Lake = unique(iso1$lake),
                          SEA = SEA.B.modes1$value)

#### Read in environmental data ####
env <- read_csv("environmental_data.csv")


##### Combine Results and Data into One Dataset #######
env_results <- left_join(results2.new2, env, by = c("lakes" = "Lake"))
env_results <- left_join(env_results, iso.results, by = c("lakes" = "Lake"))

data <- env_results
data <- data %>% 
  mutate(BIC.TNW = 1 - WonT) %>% 
  mutate(BIC.TNW.d = BIC.d/TNW.d) %>% 
  rename(Lake = lakes)


############################# Models and Figures ########################################
#### How many pops showed specialization in at least one dimension (alpha = 0.05) ####
spec <- data %>% 
  filter(p.value < 0.05 | p.WonT < 0.05) %>% 
  select(Lake, p.value, p.WonT)
print(nrow(spec)/nrow(data))

## Function to standardize (z-score) regression inputs
stand <- function(x){
  a <- rep(NA,length(x))
  m <- mean(x)
  s <- sd(x)
  for(i in 1:length(x)){
    a[i] <- (x[i]-m)/s
  }#i
  a
}#stand

## standardize predictor variables
data.s <- data %>% 
  mutate(Lake.area.ha.s = stand(Lake.area.ha), 
         Elevation.s = stand(Elevation),
         TNW.s = stand(TNW),
         richness.s = stand(richness))

## PNW(size)
summary(lm(data = data.s, TNW ~ Elevation.s + Lake.area.ha.s))

## PNW(taxa)
summary(lm(data = data.s, richness ~ Elevation.s + Lake.area.ha.s))

## PNW(iso)
summary(lm(data = data.s, SEA ~ Elevation.s + Lake.area.ha.s))

## IS(taxa) (Eadj)
summary(t1 <- betareg(Eadj ~ richness.s, link = "logit", data = data.s)) #this one; deltaAIC < 2
summary(t2 <- betareg(Eadj ~ richness.s, link = "loglog", data = data.s))

## IS(size) (BIC/PNW)
summary(t3 <- betareg(BIC.TNW ~ TNW.s, link = "loglog", data = data.s))


summary(t4 <- lm(TNW ~ BIC + WIC, data = data))
### Linear correlations between INW and IS
cor.test(data$WIC, data$BIC, alternative = "two")
cor.test(data$Eadj, data$median.richness, alternative = "two")

################# Figures #################################################################
############ Plot PNW vs PNW #############
cor.test(data$SEA, data$richness, alternative = "two.sided")
pnw1 <- data %>% 
  ggplot(aes(x = SEA, y = richness)) +
  geom_point(size = 4) +
  labs(x = expression(PNW[iso]),
       y = expression(PNW[taxa])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  annotate("text",
           x = 1.66,
           y = 14,
           label = "r = 0.02\\nP = 0.94")

cor.test(data$SEA, data$TNW, alternative = "two.sided")
pnw2 <- data %>% 
  ggplot(aes(x = SEA, y = TNW)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm",
              se = F,
              size = 2,
              color = "black") +
  labs(x = expression(PNW[iso]),
       y = expression(PNW[size])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  annotate("text",
           x = 1.5,
           y = 35,
           label = "r = -0.62\\nP < 0.05")

cor.test(data$richness, data$TNW, alternative = "two.sided")
pnw3 <- data %>% 
  ggplot(aes(x = richness, y = TNW)) +
  geom_point(size = 4) +
  labs(x = expression(PNW[taxa]),
       y = expression(PNW[size])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  annotate("text",
           x = 14,
           y = 35,
           label = "r = -0.39\\nP = 0.19")

pnw4 <- data %>% #TNW.d
  ggplot(aes(x = TNW, y = TNW.d)) +
  geom_point(size = 4) +
  labs(x = expression(PNW[size]),
       y = expression(PNW[taxa])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) 

cor.test(data$Eadj, data$BIC.TNW, alternative = "two.sided")
is.cor <- data %>% 
  ggplot(aes(x = Eadj, y = BIC.TNW)) +
  geom_point(size = 4) +
  labs(x = expression(Ind.~Spec.[taxa]~(E[adj])),
       y = expression(Ind.~Spec.[size]~(BIC/PNW[size]))) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  annotate("text",
           x = 0.7,
           y = 0.1,
           label = "r = 0.16\\nP = 0.61")

pnw_all <- pnw3 + pnw2 + pnw1 +
  plot_annotation(tag_levels = "A")
pnw_and_is_cor <-(pnw3 + pnw2)/(pnw1 + is.cor) +
  plot_annotation(tag_levels = "A")


tnwdat1.2 <- data %>% #update to incorporate null predictions on plot
  select(TNW, BIC, WIC, null.WIC, null.BIC) %>% 
  pivot_longer(cols = c(BIC, WIC, null.WIC, null.BIC), names_to = "type", values_to = "val") 

tnw1_2 <- tnwdat1.2 %>% 
  filter(type %in% c("BIC", "WIC")) %>% 
  ggplot(aes(x = TNW, y = val, color = type)) +
  geom_point(aes(x = TNW, y = val, color = type), 
             size = 4) +
  geom_smooth(method = "lm",
              size = 2, 
              se = F) +
  stat_smooth(data = filter(tnwdat1.2, type == "null.BIC"), 
              geom = "line", method = "lm",  color = "dodgerblue4", linetype= "dashed", se = F,
              alpha = 0.6) +
  scale_color_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  labs(x = expression(PNW[size]),
       y = expression ("Components of PNW"[size]),
       color = "PNW\\nComponent") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.text = element_text(size = rel(0.8)),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  annotate("text", 
           x = 38, 
           y = 14,
           label = "WIC",
           color = "darkgoldenrod") +
  annotate("text", 
           x = 38, 
           y = 21,
           label = "BIC",
           color = "dodgerblue4") +
  annotate("text", 
           x = 38, 
           y = 7,
           label = "Null BIC",
           color = "dodgerblue4",
           alpha = 0.8) 

### Combine tnw3 and tnw4 into one plot that shows both BIC and WIC/TNW ###
# make prediction line from betaregression model 
WICmod2 <- betareg(data = data, WonT ~ TNW, link = "loglog")
BICmod2 <- betareg(data = data, BIC.TNW ~ TNW, link = "loglog")
null.BICmod2 <- betareg(data = data, null.BIC.TNW ~ TNW, link = "loglog")
WICpreds <- data.frame(WIC.TNW = predict(WICmod2, data.frame(TNW = seq(from = 0.1, 
                                                                       to =50,
                                                                       by =  0.01))),
                       TNW = seq(0.1,50, 0.01))
BICpreds <- data.frame(BIC.TNW = predict(BICmod2, data.frame(TNW = seq(from = 0.1, 
                                                                       to =50,
                                                                       by =  0.01))),
                       TNW = seq(0.1,50, 0.01))

null.BICpreds <- data.frame(null.BIC.TNW = 
                              predict(null.BICmod2, data.frame(TNW = seq(from = 0.1, 
                                                                         to =50,
                                                                         by =  0.01))),
                            TNW = seq(0.1,50, 0.01))
prop.data <- data %>%
  select(WonT, BIC.TNW, TNW) %>% 
  pivot_longer(cols = c(WonT, BIC.TNW), 
               names_to = "type", values_to = "prop") %>% 
  mutate(type1 = ifelse(type == "WonT", "WIC", "BIC"))

pr <- summary(betareg(BIC.TNW ~ TNW, link = "loglog", data = data))$pseudo.r.squared

tnw5 <- prop.data %>% 
  filter(type1 %in% c("BIC", "WIC")) %>% 
  ggplot(aes(x = TNW, y = prop, color = type1)) +
  geom_point(aes(x = TNW, y = prop, color = type1), 
             size = 4) +
  scale_color_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  geom_line(data = BICpreds, aes(y = BIC.TNW, x = TNW), 
            col="dodgerblue4",
            size = 2)  +
  stat_smooth(data = null.BICpreds, aes(y = null.BIC.TNW, x = TNW), 
              geom = "line", se = F, col="dodgerblue4",
              linetype = "dashed", alpha = 0.6)  +
  geom_line(data = WICpreds, aes(y = WIC.TNW, x = TNW), 
            col="darkgoldenrod1",
            size = 2) +
  labs(x = expression(PNW[size]),
       y = expression ("Proportion of PNW"[size]),
       color = "PNW\\nComponent") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.text = element_text(size = rel(0.8)),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  annotate("text",
           x = 45, 
           y = 0.82, 
           label = "BIC", 
           color = "dodgerblue4") +
  annotate("text",
           x = 45, 
           y = 0.18, 
           label = "WIC", 
           color = "darkgoldenrod") +
  annotate("text",
           x = 45, 
           y = 0.48, 
           label = "Null BIC", 
           color = "dodgerblue4",
           alpha = 0.6) 


cor.test(data$BIC, data$WIC, alternative = "greater")

tnw6 <- data %>% 
  ggplot(aes(x = WIC, y = BIC)) +
  geom_point(size = 4, color = "black") + 
  geom_smooth(method = "lm",
              se = F,
              size = 2,
              color = "black") +
  geom_smooth(aes(x = WIC, y = null.BIC),
              method = "lm",
              se = F,
              size= 0.8,
              linetype = "dashed",
              color = "darkgray") +
  labs(x = expression(WIC[size]),
       y = expression(BIC[size]))+
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  annotate("text", 
           x = 10, 
           y = 25, 
           label = "r = 0.90\\nP < 0.0001", 
           size = 4, color = "black")


##### TNW for prey taxa #####
## Combine tnw1 and tnw2 into one plot ##
tnwdat1.2.d <- data %>% 
  select(TNW.d, BIC.d, WIC.d, null.WIC.d, null.BIC.d) %>% 
  pivot_longer(cols = c(BIC.d, WIC.d, null.WIC.d, null.BIC.d), 
               names_to = "type", values_to = "val") 

summary(lm(data = data, BIC.d ~ TNW.d))
summary(lm(data = data, WIC.d ~ TNW.d))

tnw1_2.d <- tnwdat1.2.d %>% 
  filter(type %in% c("WIC.d", "BIC.d")) %>% 
  ggplot(aes(x = TNW.d, y = val, color = type)) +
  geom_point(aes(x = TNW.d, y = val, color = type), 
             size = 4) +
  geom_smooth(method = "lm",
              size = 2, 
              se = F) +
  stat_smooth(data = filter(tnwdat1.2.d, type == "null.BIC.d"),
              geom = "line", method = "lm", color = "dodgerblue4", linetype= "dashed", se = F,
              alpha = 0.6) +
  scale_color_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  labs(x = expression(PNW[taxa]),
       y = expression ("Components of PNW"[taxa]),
       color = "PNW\\nComponent") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.text = element_text(size = rel(0.8)),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  annotate("text",
           x = 1.5,
           y = 1.0,
           label = "WIC",
           color = "darkgoldenrod") +
  annotate("text",
           x = 1.5,
           y = 0.5,
           label = "BIC",
           color = "dodgerblue4") +
  annotate("text",
           x = 1.5,
           y = 0.2,
           label = "Null BIC",
           color = "dodgerblue4",
           alpha = 0.8)

### Combine tnw3 and tnw4 into one plot that shows both BIC and WIC/TNW ###
# make prediction line from betaregression model 
data <- data %>% 
  mutate(BIC.TNW.d = BIC.d/TNW.d)

WICmod2.d <- betareg(data = data, WonT.d ~ TNW.d, link = "logit")
summary(WICmod2.d)

WICpreds.d <- data.frame(WIC.TNW.d = predict(WICmod2.d, data.frame(TNW.d = seq(from = 0.01, 
                                                                               to =2,
                                                                               by =  0.01))),
                         TNW.d = seq(0.01,2, 0.01))

BICmod2.d <- betareg(data = data, BIC.TNW.d ~ TNW.d, link = "logit")
summary(BICmod2.d)

BICpreds.d <- data.frame(BIC.TNW.d = predict(BICmod2.d, data.frame(TNW.d = seq(from = 0.01, 
                                                                               to =2,
                                                                               by =  0.01))),
                         TNW.d = seq(0.01, 2, 0.01))

null.BICmod2.d <- betareg(data = data, null.BIC.TNW ~ TNW.d, link = "logit")


null.BICpreds.d <- data.frame(null.BIC.TNW.d = 
                                predict(null.BICmod2.d, data.frame(TNW.d = seq(from = 0.01, 
                                                                               to =2,
                                                                               by =  0.01))),
                              TNW.d = seq(0.01, 2, 0.01))

prop.data.d <- data %>%
  select(WonT.d, BIC.TNW.d, TNW.d) %>% 
  pivot_longer(cols = c(WonT.d, BIC.TNW.d), names_to = "type", values_to = "prop") %>% 
  mutate(type1 = ifelse(type == "WonT.d", "WIC", "BIC"))

pr <- summary(betareg(BIC.TNW.d ~ TNW.d, link = "logit", data = data))$pseudo.r.squared
summary(lm(data = data, BIC.TNW.d ~ TNW.d))


tnw5.d <- prop.data.d %>% 
  filter(type1 %in% c("BIC", "WIC")) %>% 
  ggplot(aes(x = TNW.d, y = prop, color = type1)) +
  geom_point(aes(x = TNW.d, y = prop, color = type1), 
             size = 4) +
  scale_color_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  geom_line(data = BICpreds.d, aes(y = BIC.TNW.d, x = TNW.d), 
            col="dodgerblue4",
            size = 2)  +
  geom_line(data = null.BICpreds.d, aes(y = null.BIC.TNW.d, x = TNW.d),
            col="dodgerblue4",
            alpha = 0.6, linetype = "dashed")  +
  geom_line(data = WICpreds.d, aes(y = WIC.TNW.d, x = TNW.d), 
            col="darkgoldenrod1",
            size = 2) +
  labs(x = expression(PNW[taxa]),
       y = expression ("Proportion of PNW"[taxa]),
       color = "PNW\\nComponent") +
  ylim(0,1) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.text = element_text(size = rel(0.8)),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  annotate("text",
           x = 0.4, 
           y = 0.28, 
           label = "BIC",
           color = "dodgerblue4") +
  annotate("text",
           x = 0.4, 
           y = 0.11, 
           label = "Null BIC",
           color = "dodgerblue4", alpha = 0.8) +
  annotate("text",
           x = 0.4, 
           y = 0.80, 
           label = "WIC",
           color = "darkgoldenrod")


cor.test(data$WIC.d, data$BIC.d)
tnw6.d <- data %>% 
  ggplot(aes(x = WIC.d, y = BIC.d)) +
  geom_point(size = 4, color = "black") + 
  geom_smooth(method = "lm",
              se = F,
              size = 2,
              color = "black") +
  geom_smooth(aes(x = WIC.d, y = null.BIC.d),
              method = "lm",
              se = F,
              size = 0.8,
              linetype = "dashed",
              color = "darkgray") +
  labs(x = expression(WIC[taxa]),
       y = expression(BIC[taxa]))+
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  annotate("text", 
           x = 1.215, 
           y = 0.98, 
           label = "r = 0.60\\nP < 0.05", 
           size = 4, color = "black")

###### NMDS Plot ######
library(vegan)
# Make community matrix
gut.com <- gut.wide %>% 
  select(Annelidae:Trombidiformes)
gut.groups <- gut.wide %>% 
  select(lake)

## Run ordination
set.seed(2)
nmds1 <- metaMDS(as.matrix(gut.com),
                 k = 3,
                 maxit = 10000, #max number of iterations
                 trymax = 5000, # max number of random starts
                 wascores = TRUE # uses wascores function to calculate species scores)
)
mds.out <- as.data.frame(scores(nmds1)) # extract NMDS coordinates for each "site", where site is unique tank-date combination
mds.full <- bind_cols(gut.groups, mds.out) 
# add elevation data 
env1 <- env_results %>% 
  select(lakes, Elevation, Lake.area.ha, TNW.d)
mds.full <- left_join(mds.full, env1, by = c("lake" = "lakes"))

mds.avg <- mds.full %>% 
  group_by(lake) %>% 
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2),
            Elevation = mean(Elevation), area = mean(Lake.area.ha), 
            TNW.d = mean(TNW.d))

mds.spec <- tibble(species = row.names(nmds1$species), NMDS1 = nmds1$species[,1])

spec <- mds.spec %>% 
  ggplot(aes(x = NMDS1, y = 2)) +
  geom_text(aes(label = species), position = position_jitter())

#### Make plot ####
library(RColorBrewer)

base.size <- 22
ord.plot1 <- ggplot() +
  geom_point(data = mds.avg,
             aes(x = NMDS1, y = NMDS2, fill = TNW.d),
             size = 5,
             shape = 21,
             color = "black") +
  geom_point(data = mds.full,
             aes(x = NMDS1, y = NMDS2, fill = TNW.d),
             pch = 21,
             size = 2,
             alpha = 0.3) +
  coord_cartesian(xlim = c(-0.5,1.1),
                  ylim = c(-1.1,1.5)) +
  scale_fill_viridis_c(direction = -1, option = "B", end = 0.8, 
                       name = expression(PNW[taxa])) +
  theme_classic(base_size = base.size) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top") 


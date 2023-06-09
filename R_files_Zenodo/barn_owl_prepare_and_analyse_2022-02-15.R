
## ----libraries--------------------------------------------------------------------------------
# Load: Libraries

library(dplyr)
library(lme4)
library(ggplot2)
library(stringr)
library(lubridate)
library(car)
library(reshape2)
library(cowplot)
library(effects)
library(knitr)
library(plot3D)
source("code/barn_owl_programs.R")

#





## ----data and fix, warning=FALSE, message=FALSE-----------------------------------------------
# Load: Data

# Load data set
dd_f <- read.csv2("data/barn_owl_females_switzerland_20171129.txt", dec = ".")
# These data was obtained from Isabelle Henry Dufresnes on the 29/11-17

# View
dim(dd_f)
names(dd_f)
head(dd_f[, 1:10])

# Load population sizes
pop <- read.csv2("data/barn_owl_popsize_switzerland_20171129.txt")

# View
head(pop)
tail(pop)

#



## ---------------------------------------------------------------------------------------------
# Check: Duplicates

# Duplicated individuals within years
table(duplicated(dd_f[, c("ring", "year")]))
#dd_f$ring[duplicated(dd_f[, c("ring", "year")])]
#View(subset(dd_f, ring == "M031552"))
# Two broods per year for some females

#



## ---------------------------------------------------------------------------------------------

# Rename
dd_f$genotype[dd_f$genotype == "VV"] <- "ww"
dd_f$genotype[dd_f$genotype == "VI"] <- "wr"
dd_f$genotype[dd_f$genotype == "II"] <- "rr"
# Table
table(dd_f$genotype)



## ---------------------------------------------------------------------------------------------
qplot(subset(dd_f, !duplicated(ring) & !is.na(genotype))$genotype) + 
  xlab("MC1R Genotype") +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 420))


## ----allele obs-------------------------------------------------------------------------------
# Create: Annual allele frequencies

# Define the number of each allele, r and w
dd_f$allele_r <- NA
dd_f$allele_r[dd_f$genotype == "ww"] <- 0
dd_f$allele_r[dd_f$genotype == "wr"] <- 1
dd_f$allele_r[dd_f$genotype == "rr"] <- 2
dd_f$allele_w <- NA
dd_f$allele_w[dd_f$genotype == "rr"] <- 0
dd_f$allele_w[dd_f$genotype == "wr"] <- 1
dd_f$allele_w[dd_f$genotype == "ww"] <- 2

# Table
table(dd_f$allele_r)
table(dd_f$allele_w)
table(dd_f$allele_w, dd_f$allele_r)

# Annual allele frequencies (removing duplicated individuals)
dd_allele_obs <- dd_f %>%
  filter(!duplicated(paste(ring,year))) %>%
  group_by(year) %>%
  summarise(n_tot = n(),
            missing = sum(is.na(allele_r)),
            n = n_tot-missing,
            r = sum(allele_r, na.rm = TRUE)/(n*2),
            w = sum(allele_w, na.rm = TRUE)/(n*2))

# To data frame
dd_allele_obs <- as.data.frame(dd_allele_obs)

# Replace NaN
dd_allele_obs[is.nan(dd_allele_obs)] <- NA

# View
head(dd_allele_obs)
#tail(dd_allele_obs)
#str(dd_allele_obs)

# To long format
dd_allele_obs <- melt(dd_allele_obs, id.vars = c("year", "n_tot", "missing", "n"),
                variable.name = "allele", value.name = "freq")

# View
head(dd_allele_obs)
tail(dd_allele_obs)

# Mean allele frequencies
dd_allele_obs %>%
  group_by(allele) %>%
  summarise(mean_allele_freq = mean(freq))

#



## ----geno obs---------------------------------------------------------------------------------

# Create: Annual genotype observations

# Annual frequencies (removing duplicated individuals within years)
dd_geno_obs <- dd_f %>%
    filter(!duplicated(paste(ring,year))) %>%
  group_by(year) %>%
  summarise(n_tot = n(),
            missing = sum(is.na(allele_r)),
            n = n_tot-missing,
            rr_n = sum(genotype == "rr", na.rm = TRUE),
            wr_n = sum(genotype == "wr", na.rm = TRUE),
            ww_n = sum(genotype == "ww", na.rm = TRUE),
            rr = rr_n/n,
            wr = wr_n/n,
            ww = ww_n/n)

# To data frame
dd_geno_obs <- as.data.frame(dd_geno_obs)

# Replace NaN (when there are no observations)
dd_geno_obs[is.nan(dd_geno_obs)] <- NA

# View
#head(dd_geno_obs)
tail(dd_geno_obs)
#str(dd_geno_obs)

# To long format
# Genotype frequencies
dd_geno_obs_long <- melt(select(dd_geno_obs, -rr_n, -wr_n, -ww_n, -n_tot, -missing, -n), id.vars = c("year"), variable.name = "genotype", value.name = "freq")
# Sample size for each genotype
dd_geno_obs_long <- data.frame(dd_geno_obs_long, n = melt(select(dd_geno_obs, -rr, -wr, -ww, -n_tot, -missing, -n), id.vars = c("year"), variable.name = "genotype", value.name = "n")$n)

# View
#head(dd_geno_obs_long)
tail(dd_geno_obs_long)

#


## ----graphs geno allele, fig.width=6, fig.height=4--------------------------------------------
# Explore: Explore allele and genotype data

# Graph allele frequencies
allele_plot <- ggplot(data = subset(dd_allele_obs, allele == "w" & year > 1995),
       mapping = aes(y = freq, x = year)) +
  geom_line(size = 0.25) +
  #geom_text(aes(label = n), nudge_y = c(-0.02,-0.02, 0.02,0.02,0.02, -0.02, 0.02,0.02,0.02, -0.02, 0.02, -0.02, 0.02,0.02,0.02,0.02, -0.02,-0.02, 0.02, -0.02, 0.02), size = 8/(14/5)) +
  geom_point(size = 2.5, shape = 21, stroke = 0.25, fill = "black", colour = "white") +
  scale_x_continuous(breaks = seq(1990, 2016, by = 5), limits = c(1990,2016)) +
  scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1), limits = c(0.7,1)) +
  xlab("Year") +
  ylab("Allele frequency (W-allele)") +
  theme_cowplot(font_size = 10, line_size = 0.25)

allele_plot

ggsave(plot = allele_plot, filename = "output/allele_frequency_year.pdf", height = 80, width = 120, units = "mm")

# Graph genotype frequencies
genotype_plot <- ggplot(data = subset(dd_geno_obs_long, year > 1995),
       mapping = aes(y = freq, x = year, group = genotype)) +
  geom_line(size = 0.25) +
  geom_point(size = 2.5, shape = 21, stroke = 0.25, fill = "black", colour = "white") +
  scale_x_continuous(breaks = seq(1990, 2016, by = 5), limits = c(1990,2016)) +
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  xlab("Year") +
  ylab("Genotype frequency") +
  theme_cowplot(font_size = 10, line_size = 0.25)

genotype_plot

ggsave(plot = genotype_plot, filename = "output/genotype_frequency_year.pdf", height = 80, width = 120, units = "mm")


#


## ---------------------------------------------------------------------------------------------
# Diameter
qplot(factor(age), diabreast, data = dd_f, geom = "boxplot")
qplot(factor(age), diawing, data = dd_f, geom = "boxplot")

# Number
qplot(factor(age), nbbreast, data = dd_f, geom = "boxplot")
qplot(factor(age), nbwing, data = dd_f, geom = "boxplot")

# Colour
qplot(factor(age), colourbreast, data = dd_f, geom = "boxplot")
qplot(factor(age), colourbelly, data = dd_f, geom = "boxplot")
qplot(factor(age), colourflank, data = dd_f, geom = "boxplot")
qplot(factor(age), colourwing, data = dd_f, geom = "boxplot")



## ---------------------------------------------------------------------------------------------
# Count
with(dd_f, table(!is.na(diabreast), age))
with(dd_f, table(!is.na(diawing), age))
with(dd_f, table(!is.na(nbbreast), age))
with(dd_f, table(!is.na(nbwing), age))
with(dd_f, table(!is.na(colourbreast), age))
with(dd_f, table(!is.na(colourwing), age))



## ---------------------------------------------------------------------------------------------
# Diameter
with(dd_f, table(!is.na(diabreast), year)) # Missing years: 1990-1993
with(dd_f, table(!is.na(diabelly), year)) # Missing years: 1990-1993, 2005, 2012
with(dd_f, table(!is.na(diaflank), year)) # Missing years: 1990-1993, 2005, 2012, 2016
with(dd_f, table(!is.na(diawing), year)) # Missing years: 1990-1993, 2005, 2012, 2016
# Number
with(dd_f, table(!is.na(nbbreast), year)) # Missing years: 1990-1993
with(dd_f, table(!is.na(nbbelly), year)) # Missing years: 1990-1993, 2005
with(dd_f, table(!is.na(nbflank), year)) # Missing years: 1990-1993, 2005, 2016
with(dd_f, table(!is.na(nbwing), year)) # Missing years: 1990-1993, 2005, 2016
# Colour
with(dd_f, table(!is.na(colourbreast), year)) # Missing years: 1990-1993
with(dd_f, table(!is.na(colourbelly), year)) # Missing years: 1990-1993
with(dd_f, table(!is.na(colourflank), year)) # Missing years: 1990-1993
with(dd_f, table(!is.na(colourwing), year)) # Missing years: 1990-1993



## ---------------------------------------------------------------------------------------------
# Traits only missing in 1990-1993
dd_f_age <- dd_f %>%
  filter(!is.na(age) & !is.na(colourbreast) & !is.na(nbbreast) & !is.na(diabreast) & age <= 3) %>%
  group_by(ring, age, year) %>%
  summarise(colourbreast = mean(colourbreast),
            colourbelly = mean(colourbelly),
            colourflank = mean(colourflank),
            colourwing = mean(colourwing),
            diabreast = mean(diabreast),
            nbbreast = mean(nbbreast)) %>%
  group_by(ring) %>%
  mutate(n=n()) %>%
  filter(n >= 2)

table(dd_f_age$age)

# Traits missing in 1990-1993 and some later years
dd_f_age2 <- dd_f %>%
  filter(!is.na(age) & !is.na(diawing) & !is.na(nbwing) & age <= 3) %>%
  group_by(ring, age, year) %>%
  summarise(diabelly = mean(diabelly),
            diaflank = mean(diaflank),
            diawing = mean(diawing),
            nbbelly = mean(nbbelly),
            nbflank = mean(nbflank),
            nbwing = mean(nbwing)) %>%
  group_by(ring) %>%
  mutate(n=n()) %>%
  filter(n >= 2)

table(dd_f_age2$age)


## ---------------------------------------------------------------------------------------------
# Colour
summary(mod_cbr <- lmer(colourbreast~factor(age) + (age|ring) + (1|year), data = dd_f_age))
summary(mod_cbe <- lmer(colourbelly~factor(age) + (age|ring) + (1|year), data = dd_f_age))
summary(mod_cf <- lmer(colourflank~factor(age) + (age|ring) + (1|year), data = dd_f_age))
summary(mod_cw <- lmer(colourwing~factor(age) + (age|ring) + (1|year), data = dd_f_age))
drop1(mod_cbr, test = "Chisq")
drop1(mod_cbe, test = "Chisq")
drop1(mod_cf, test = "Chisq")
drop1(mod_cw, test = "Chisq")

# Diameter
summary(mod_dbr <- lmer(diabreast~factor(age) + (age|ring) + (1|year), data = dd_f_age))
summary(mod_dbe <- lmer(diabelly~factor(age) + (age|ring) + (1|year), data = dd_f_age2))
summary(mod_df <- lmer(diaflank~factor(age) + (age|ring) + (1|year), data = dd_f_age2))
summary(mod_dw <- lmer(diawing~factor(age) + (age|ring) + (1|year), data = dd_f_age2))
drop1(mod_dbr, test = "Chisq")
drop1(mod_dbe, test = "Chisq")
drop1(mod_df, test = "Chisq")
drop1(mod_dw, test = "Chisq")

# Number
summary(mod_nbr <- lmer(nbbreast~factor(age) + (age|ring) + (1|year), data = dd_f_age))
summary(mod_nbe <- lmer(nbbelly~factor(age) + (age|ring) + (1|year), data = dd_f_age2))
summary(mod_nf <- lmer(nbflank~factor(age) + (age|ring) + (1|year), data = dd_f_age2))
summary(mod_nw <- lmer(nbwing~factor(age) + (age|ring) + (1|year), data = dd_f_age2))
drop1(mod_nbr, test = "Chisq")
drop1(mod_nbe, test = "Chisq")
drop1(mod_nf, test = "Chisq")
drop1(mod_nw, test = "Chisq")
# No age effects for spot number.


## ---------------------------------------------------------------------------------------------
# Colour
fit_cbr <- predict(mod_cbr, newdata=data.frame(age=1:3), re.form=NA)
fit_cbe <- predict(mod_cbe, newdata=data.frame(age=1:3), re.form=NA)
fit_cf <- predict(mod_cf, newdata=data.frame(age=1:3), re.form=NA)
fit_cw <- predict(mod_cw, newdata=data.frame(age=1:3), re.form=NA)
qplot(1:3, fit_cbr, geom = "line", ylab = "Colour breast")
qplot(1:3, fit_cbe, geom = "line", ylab = "Colour belly")
qplot(1:3, fit_cf, geom = "line", ylab = "Colour flank")
qplot(1:3, fit_cw, geom = "line", ylab = "Colour wing")

# Diameter
fit_dbr <- predict(mod_dbr, newdata=data.frame(age=1:3), re.form=NA)
fit_dbe <- predict(mod_dbe, newdata=data.frame(age=1:3), re.form=NA)
fit_df <- predict(mod_df, newdata=data.frame(age=1:3), re.form=NA)
fit_dw <- predict(mod_dw, newdata=data.frame(age=1:3), re.form=NA)
qplot(1:3, fit_dbr, geom = "line", ylab = "Diameter breast")
qplot(1:3, fit_dbe, geom = "line", ylab = "Diameter belly")
qplot(1:3, fit_df, geom = "line", ylab = "Diameter flank")
qplot(1:3, fit_dw, geom = "line", ylab = "Diameter wing")


## ---------------------------------------------------------------------------------------------
dd_f1 <- dd_f %>%
  filter(!is.na(age)) %>%
  group_by(age) %>%
  mutate(colourbreast1 = round(colourbreast - fit_cbr[ifelse(age <= 3, age, 3)] + fit_cbr[1], 2),
         colourbelly1 = round(colourbelly - fit_cbe[ifelse(age <= 3, age, 3)] + fit_cbe[1], 2),
         colourflank1 = round(colourflank - fit_cf[ifelse(age <= 3, age, 3)] + fit_cf[1], 2),
         colourwing1 = round(colourwing - fit_cw[ifelse(age <= 3, age, 3)] + fit_cw[1], 2),
         diabreast1 = round(diabreast - fit_dbr[ifelse(age <= 3, age, 3)] + fit_dbr[1], 2),
         diabelly1 = round(diabelly - fit_dbe[ifelse(age <= 3, age, 3)] + fit_dbe[1], 2),
         diaflank1 = round(diaflank - fit_df[ifelse(age <= 3, age, 3)] + fit_df[1], 2),
         diawing1 = round(diawing - fit_dw[ifelse(age <= 3, age, 3)] + fit_dw[1], 2),
         nbbreast1 = nbbreast,
         nbbelly1 = nbbelly,
         nbflank1 = nbflank,
         nbwing1 = nbwing)



## ----ind data---------------------------------------------------------------------------------

# Summarise within years for each individual
dd_i <- dd_f1 %>%
  arrange(ring, age) %>%
  group_by(ring) %>%
  mutate(any1 = any(age == 1),
         birthyear = unique(year-age)[1],
         diabreast1 = mean(diabreast1),
         diabelly1 = mean(diabelly1),
         diaflank1 = mean(diaflank1),
         diawing1 = mean(diawing1),
         nbbreast1 = mean(nbbreast1),
         nbbelly1 = mean(nbbelly1),
         nbflank1 = mean(nbflank1),
         nbwing1 = mean(nbwing1),
         cbreast1 = mean(colourbreast1),
         cbelly1 = mean(colourbelly1),
         cflank1 = mean(colourflank1),
         cwing1 = mean(colourwing1)) %>%
  group_by(ring, year) %>%
  summarise(any1 = first(any1),
            sex = first(sex),
            age = first(age),
            birthyear = first(birthyear),
            allele_r = first(allele_r),
            allele_w = first(allele_w),
            genotype = first(genotype),
            diabreast = first(diabreast1),
            diabelly = first(diabelly1),
            diaflank = first(diaflank1),
            diawing = first(diawing1),
            nbbreast = first(nbbreast1),
            nbbelly = first(nbbelly1),
            nbflank = first(nbflank1),
            nbwing = first(nbwing1),
            cbreast = first(cbreast1),
            cbelly = first(cbelly1),
            cflank = first(cflank1),
            cwing = first(cwing1)) %>%
  ungroup()

# View
dim(dd_i)
#head(dd_i)
tail(dd_i)

# Define reference genotype
dd_i$genotype <- factor(dd_i$genotype, levels = c("ww", "wr", "rr"))

# Group genotypes WR and RR (they have similar phenotypes)
dd_i$genotype_g <- dd_i$genotype
dd_i$genotype_g[dd_i$genotype_g == "rr"] <- "wr"

# Add rownumbers
dd_i$obs_no <- 1:dim(dd_i)[1]


## ----rec and surv-----------------------------------------------------------------------------
# Create: Estimate of individual fecundity and survival (pre-breeding census)

# Calculate the number of recruits per brood
dd_f$recruits <- rowSums(dd_f[, c("recruits_f", "recruits_m")], na.rm=TRUE)

# Table
with(dd_f, table(recruits, exclude = NULL))

# Arrange
dd_f <- arrange(dd_f, ring, year)


# Calculate survival

# Define survival and set no NA for last year of data (i.e. don't know surivival until 2017, when 2016 is the last year of data)
dd_f$survival <- 0
dd_f$survival[dd_f$year == max(dd_f$year)] <- NA

# Define last year for each individual
dd_f <- dd_f %>%
  group_by(ring) %>%
  mutate(lastyear = max(year))

# Survival = 1 when year is larger than max(year)
dd_f$survival[dd_f$year < dd_f$lastyear] <- 1

# Table
table(dd_f$survival, exclude = NULL)
table(dd_f$year, dd_f$survival, exclude = NULL)

# Cross fostering
table(dd_f$cross_fostered, exclude = NULL)

# Summarise within years for each individual
dd_a <- dd_f %>%
  group_by(ring, year) %>%
  summarise(recruits = sum(recruits),
            cross_fostered = max(cross_fostered, na.rm = TRUE),
            survival = first(survival))

# View
head(dd_a)

#



## ----rec and surv join------------------------------------------------------------------------
# Join
dd_i <- left_join(dd_i, dd_a, by = c("ring", "year"))

# View
dim(dd_i)
head(dd_i[100:120, c(1:2, (ncol(dd_i)-8):ncol(dd_i))])

# Table to check missing
with(dd_i, table(year, is.na(survival)))
with(dd_i, table(year, is.na(recruits)))

#


## ---------------------------------------------------------------------------------------------
# Year gaps in individual records (count the difference between number of records and the number of years between first and last record)
dd_i %>%
  group_by(ring) %>%
  summarise(gaps = c(1+max(year)-min(year))-n(),
            ind = 1) %>%
  ungroup() %>%
  summarise(gaps_sum = sum(gaps),
            ind_sum = sum(ind))


## ---------------------------------------------------------------------------------------------
with(pop, mean(popsize))
with(pop, se(popsize))
with(pop, sd(popsize))
with(pop, range(popsize))


## ---------------------------------------------------------------------------------------------
# Plot
n_plot <- ggplot(data = pop, mapping = aes(x = year, y = popsize)) +
  geom_line(size = 0.25) +
#  geom_hline(yintercept=100, size = 0.25) +
  geom_point(shape = 21, size = 2.5, stroke = 0.25, fill = "black", colour = "white") +
  scale_x_continuous(breaks = seq(1990, 2016, by = 5)) +
  xlab("Year") +
  ylab("Population size (N)") +
  theme_cowplot(font_size = 10, line_size = 0.25)

n_plot

ggsave(plot= n_plot, filename = "output/popsize_year.pdf", height = 80, width = 120, units = "mm")



## ----density----------------------------------------------------------------------------------
pop$delta_n <- c(diff(pop$popsize), NA)
summary(mod_deltaN <- lm(I(delta_n/popsize)~popsize, data = pop))
drop1(mod_deltaN, test = "F")


## ----density plot-----------------------------------------------------------------------------
# Range
range(pop$popsize)
# Predict delta N from model for plotting (with 95 % CI)
pred_deltaN <- data.frame(popsize = 52:187, predict(mod_deltaN, newdata = data.frame(popsize = 52:187), interval = "confidence"))
pred_deltaN <- rename(pred_deltaN, delta_nn = fit)

deltann_plot <- ggplot(aes(x = popsize, y = delta_nn), data = pred_deltaN) +
  geom_point(aes(y = delta_n/popsize), shape = 21, size = 2.5, stroke =0, fill = "black", colour = "white", data = pop) +
  geom_line() +
  geom_line(aes(x=popsize, y = lwr), linetype = 2) +
  geom_line(aes(x=popsize, y = upr), linetype = 2) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=46, y=1.1, lab="a"), fontface = "bold", size = 10/(14/5)) +
  xlab("Population size (N)") +
  ylab(expression(Delta~N~"/N"~''%+-%''~95~'%'~CI)) +
  theme_cowplot(font_size = 10, line_size = 0.25)

deltann_plot

ggsave(plot= deltann_plot, filename = "output/density_dependence_deltann.pdf", height = 80, width = 90, units = "mm")


## ----growth rate------------------------------------------------------------------------------

# Calculate growth rate from population censuses
pop$lambda <- lead(pop$popsize)/pop$popsize

# Mean lambda
mean(pop$lambda, na.rm = TRUE)

# Mean lambda during increases
mean(pop$lambda[pop$lambda > 1], na.rm = TRUE)

# Mean lambda during decreases
mean(pop$lambda[pop$lambda < 1], na.rm = TRUE)



## ----pop join---------------------------------------------------------------------------------
# Join in populations size to data set

# Individ data set
dd_i <- left_join(dd_i, pop[, c("year", "popsize", "delta_n")], by = "year")

# View
dim(dd_i)
head(dd_i)

# Mean center popsize
dd_i <- dd_i %>%
  mutate(popsize_c = scale(popsize, scale = FALSE))

#


## ----age--------------------------------------------------------------------------------------
# Keep real age classes
dd_i$age_real <- dd_i$age

# Age counts
with(dd_i, table(age_real))

# Age counts in subset of non cross_fostered
with(subset(dd_i, cross_fostered %in% 0), table(age_real))
with(subset(dd_i, cross_fostered %in% 0 & !is.na(genotype)), table(age_real))
with(subset(dd_i, cross_fostered %in% 0 & !is.na(diabreast)), table(age_real))
with(subset(dd_i, cross_fostered %in% 0 & !is.na(diawing)), table(age_real))
with(subset(dd_i, cross_fostered %in% 0 & !is.na(genotype) & !is.na(diawing)), table(age_real))

# Age in other subsets
with(subset(dd_i, !is.na(genotype)), table(age_real))
with(subset(dd_i, !is.na(cbreast)), table(age_real))



## ----age collapse-----------------------------------------------------------------------------

# Collapsed age class
dd_i$age[dd_i$age >= 8] <- 8

# View
with(dd_i, table(age))


## ----age genotype-----------------------------------------------------------------------------

with(dd_i, table(age, genotype))



## ----age duplicated---------------------------------------------------------------------------

dd_i %>%
  group_by(age) %>%
  summarise(n = n(),
            n_unique = n_distinct(ring),
            n_duplicated = n - n_unique)



## ----vital rates------------------------------------------------------------------------------
# Vital rates
# Calculate with recruits divided by 2
vit <- dd_i %>%
  filter(age_real <= max(age)) %>%
  group_by(age) %>%
  summarise(mean_fec = mean(recruits/2, na.rm=T),
            se_fec = se(recruits/2, na.rm=T),
            mean_sur = mean(survival, na.rm=T),
            se_sur = se(survival, na.rm=T),
            n = n())

# To data frame
vit <- as.data.frame(vit)

# Growth rate and reproductive values observed
luv_f_obs <- eigenl2(pm = promat2(pc = vit[, c("age", "mean_fec", "mean_sur")]))
luv_f_obs

# View
knitr::kable(vit, digits = 3)


## ---------------------------------------------------------------------------------------------

# Data frame for estimates
vit_d <- data.frame(age = 1:max(dd_i$age), mean_fec = NA, se_fec = NA, mean_sur = NA, se_sur = NA)

# Add needed variables
dd_i$offset2 <- log(2)
dd_i$weight2 <- 1/2

# Recruits
# Test
drop1(glm(recruits~factor(age) + I(delta_n/popsize), offset = offset2, weights = weight2, data = subset(dd_i, age_real <= max(dd_i$age)), family = "poisson"), test = "Chisq")
# Model
recK_mod <- glm(recruits~factor(age) + I(delta_n/popsize), offset = offset2, weights = weight2, data = subset(dd_i, age_real <= max(dd_i$age)), family = "poisson")
coef(summary(recK_mod))
confint(recK_mod)
# Produce latex table
#xtable::xtable(data.frame(coef(summary(recK_mod))[, 1:2], confint(recK_mod)), digits = 4)
# Predict with delta_n/popsize = 0
vit_d[,2:3] <- data.frame(predict(recK_mod, newdata = data.frame(delta_n = 0, popsize = 40, age = 1:max(dd_i$age), offset2 = log(1), weight2 = 1/1), type = "response", se.fit = TRUE))[, 1:2]

# Survival
# Test
drop1(glm(survival~factor(age) + I(delta_n/popsize), data = subset(dd_i, age_real <= max(dd_i$age)), family = "binomial"), test = "Chisq")
# Model
surK_mod <- glm(survival~factor(age) + I(delta_n/popsize), data = subset(dd_i, age_real <= max(dd_i$age)), family = "binomial")
coef(summary(surK_mod))
confint(surK_mod)
# Produce latex table
#xtable::xtable(data.frame(coef(summary(surK_mod))[, 1:2], confint(surK_mod)), digits = 4)
# Predict with delta_n
vit_d[,4:5] <- data.frame(predict(surK_mod, newdata = data.frame(delta_n = 0, popsize = 40, age = 1:max(dd_i$age)), type = "response", se.fit = TRUE))[, 1:2]

# View
knitr::kable(vit_d, digits = 3)

# Growth rate at K
eigenl2(pm = promat2(pc = vit_d[, c("age", "mean_fec", "mean_sur")]))



## ---------------------------------------------------------------------------------------------
# Range of delta_n/popsize
with(pop, range(delta_n/popsize, na.rm = T))


# Recruitment
summary(recK_mod)

# Predict recruitment from model for plotting (with 95 % CI)
pred_recK <- data.frame(popsize = 1, delta_n = seq(-0.7219251, 1.1011236, length.out = 100), predict(recK_mod, newdata = data.frame(popsize = 1, delta_n = seq(-0.7219251, 1.1011236, length.out = 100), age = 4, offset2 = log(1), weight2 = 1/1), se.fit = TRUE, type = "response"))
pred_recK <- rename(pred_recK, recruits = fit)

# Plot for recruitment
recK_plot <- ggplot(aes(x = delta_n/popsize, y = recruits), data = pred_recK) +
  geom_point(aes(y = recruits/2), shape = 20, size = 2.5, colour = "black", data = dd_i, position = position_jitter(width = 0.025), alpha = 0.1) +
  geom_line() +
  geom_line(aes(y = recruits-1.96*se.fit), linetype = 2) +
  geom_line(aes(y = recruits+1.96*se.fit), linetype = 2) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=-0.75, y=2.6, lab="b"), fontface = "bold", size = 10/(14/5)) +
  xlab(expression(Delta~N~"/N"~"(multiplicative growth rate - 1)")) +
  ylab(expression("Recruits (B/2)"~''%+-%''~95~'%'~CI)) +
  theme_cowplot(font_size = 10, line_size = 0.25)

recK_plot

ggsave(plot= recK_plot, filename = "output/density_dependence_recruits.pdf", height = 80, width = 90, units = "mm")

# Survival
summary(surK_mod)

# Predict survival from model for plotting (with 95 % CI)
pred_surK <- data.frame(popsize = 1, delta_n = seq(-0.7219251, 1.1011236, length.out = 100), predict(surK_mod, newdata = data.frame(popsize = 1, delta_n = seq(-0.7219251, 1.1011236, length.out = 100), age = 4), se.fit = TRUE, type = "response"))
pred_surK <- rename(pred_surK, survival = fit)

# Plot for survival
surK_plot <- ggplot(aes(x = delta_n/popsize, y = survival), data = pred_surK) +
  geom_point(aes(y =survival), shape = 20, size = 2.5, colour = "black", data = dd_i, position = position_jitter(width = 0.02, height = 0.02), alpha = 0.1) +
  geom_line() +
  geom_line(aes(y = survival-1.96*se.fit), linetype = 2) +
  geom_line(aes(y = survival+1.96*se.fit), linetype = 2) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=-0.75, y=1.1, lab="c"), fontface = "bold", size = 10/(14/5)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(Delta~N~"/N"~"(multiplicative growth rate - 1)")) +
  ylab(expression("Survival"~''%+-%''~95~'%'~CI)) +
  theme_cowplot(font_size = 10, line_size = 0.25)

surK_plot

ggsave(plot = surK_plot, filename = "output/density_dependence_survival.pdf", height = 80, width = 90, units = "mm")



## ---------------------------------------------------------------------------------------------
# Align plots
density_plot <- plot_grid(deltann_plot, recK_plot, surK_plot, align = "v", ncol = 1)
# View
density_plot

# Save to file
ggsave(plot = density_plot, filename = "output/density_dependence.pdf", height = 210, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Projection matrix
l_f <- promat2(pc = vit_d[, c("age", "mean_fec", "mean_sur")])

# View
knitr::kable(data.frame(l_f$l))

# Growth rate and reproductive values at K (without scaling)
luv_f_K <- eigenl2(l_f)

# Scaling konstant to get lambda = 1
konstant_f <- eulerlotka(lambda = 1, leslie = l_f$l, census = "pre", konstant = 1)$konstant
# View konstant
konstant_f

# View vital rates
cbind(vit_d[, 1],
      vit_d[, 2:3]*konstant_f,
      vit_d[, 4:5])

# Reproductive values
#eigenl2(pm = l_f$l)
luv_f <- eigenl2(pm = rbind(l_f$l[1, ]*konstant_f, l_f$l[-1,]))

# View
luv_f

# Compare to observed
luv_f_obs

# Compare to K (unscaled)
luv_f_K

# Create data frame with reproductive values for age classes and their transitions
v_age_f <- data.frame(age = luv_f$age, v = luv_f$v, v1 = luv_f$v[1], v_next = lead(luv_f$v, default = luv_f$v[max(luv_f$age)]))

# View
knitr::kable(v_age_f)

# Rename
v_age <- v_age_f


## ---------------------------------------------------------------------------------------------
ggplot(data = data.frame(age = luv_f$age, vu = c(luv_f$v, luv_f$u), type = rep(c("v", "u"), each = length(luv_f$age))), mapping = aes(x = age, y = vu)) +
  geom_line() +
  facet_grid(type~., scale = "free", switch = "y") +
  xlab("Age class") +
  ylab("Rep. values (v) and stable age dist. (u)") +
  theme(strip.placement = "outer",
        strip.background = element_blank())


## ---- cache=TRUE------------------------------------------------------------------------------
# DeltaN
deltaN <- pop$delta_n[-27]
N <- pop$popsize[-27]
quantile(pop$delta_n/pop$popsize, na.rm = T, probs = c(0.025, 0.25, 0.75, 0.975))

# Vital rates data frame
vit_N <- data.frame(age = 1:8, f = NA, s = NA)

# Vital rates
v_list_N <- data.frame(deltaN = deltaN, N = N, year = pop$year[-27], v1 = NA, v2 = NA, v3 = NA, v4 = NA, v5 = NA, v6 = NA, v7 = NA, v8 = NA)

for (i in 1:length(deltaN)){
  
  # Recruits
vit_N$f <- predict(recK_mod,
                      type = "response",
                      newdata = data.frame(delta_n = deltaN[i],
                                           popsize = N[i],
                                           age = 1:max(dd_i$age),
                                           offset2 = log(1),
                                           weight2 = 1/1))
# Survival
vit_N$s <- predict(surK_mod,
                      type = "response",
                      newdata = data.frame(delta_n = deltaN[i],
                                           popsize = N[i],
                                           age = 1:max(dd_i$age)))

# Projection matrix
l_N <- promat2(pc = vit_N)

# Scale by constant (lambda = 1)
l_N <- rbind(l_N$l[1, ]*konstant_f, l_N$l[-1,])

# Reproductive values
v_list_N[i, 4:11] <- eigenl2(pm = l_N)$v
  
}

# To long format
v_list_N <- melt(v_list_N, id.vars= c("deltaN", "N", "year"), variable.name = "v")

# Join in reproductive values at K
v_list_N <- left_join(v_list_N, data.frame(age = 1:8, v = factor(paste0("v", 1:8)), value_k = luv_f$v))
head(v_list_N)

# Test if mean rep.values deviate from rep.values at K
anova(lm(I(value-value_k)~factor(age), data = v_list_N))
with(v_list_N, mean(value-value_k))
with(v_list_N, se(value-value_k))

# Calculate means
v_list_N_mean <- v_list_N %>%
  group_by(v, age) %>%
  summarise(value_mean = mean(value),
            value_se = se(value),
            value_lwr = quantile(value,0.025),
            value_upr = quantile(value,0.975),
            value_k = value_k[1])
v_list_N_mean

# Plot v vs year
ggplot(data = v_list_N, mapping = aes(x=year, y=value, colour = v)) +
  geom_line() +
  scale_colour_manual(values = c("#01665e", "#35978f", "#80cdc1", "#c7eae5", "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a"), name = NULL) +
  geom_hline(yintercept = luv_f$v, colour = c("#01665e", "#35978f", "#80cdc1", "#c7eae5", "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a"), size = 1) +
  ylab("Reproductive value") +
  xlab("Year") +
  theme_cowplot()



## ---------------------------------------------------------------------------------------------

# Plot the difference between v at observed N and v at K
v_plot <- ggplot(data = v_list_N_mean, mapping = aes(x=age, y=value_mean-value_k, group=1)) +
  geom_point(aes(x=age, y=value-value_k),data = v_list_N, alpha = 0.2, size = 2.5) +
  geom_line(size = 0.25) +
  geom_line(aes(y=value_mean-value_k-value_se*1.96), linetype = 2, size = 0.25) +
  geom_line(aes(y=value_mean-value_k+value_se*1.96), linetype = 2, size = 0.25) +
  #geom_ribbon(aes(ymin=value_lwr-value_k, ymax=value_upr-value_k), alpha = 0.1) +
  geom_hline(yintercept=0, size = 0.25) +
  scale_x_continuous(breaks = 1:8)+
  ylab(expression("Difference ("*v[x]^"'"*"-"*v[x]*")"~''%+-%''~95~'%'~CI)) +
  xlab("Age") +
  theme_cowplot(font_size = 10, line_size = 0.25)

v_plot

# To file
ggsave(plot = v_plot, filename = "output/diff_repvalue_plot.pdf", height = 80, width = 85, units = "mm")



## ----individual rep value---------------------------------------------------------------------
# First we define 2W^*
dd_i$ws <- with(dd_i, survival + recruits/2)
dd_i$ws2 <- dd_i$ws*2

# View
table(dd_i$ws)
table(dd_i$ws2)

# Join in reproductive values
dd_i <- left_join(dd_i, v_age, by = c("age"))

# View
head(dd_i[, c("ring", "year", "v", "v1", "v_next")])

# Then we define the individual reproductive value
dd_i$w <- with(dd_i, survival*v_next + recruits*v1/2)

# Divide w by age specific reproductive value
dd_i$lambda <- with(dd_i, w/v)

# As data frame
dd_i <- as.data.frame(dd_i)

# View
head(dd_i[, c("ring", "year", "w", "lambda")])


## ----compare fitness--------------------------------------------------------------------------
# Estimate the age specific means for each measure of fitness
mean_fitness <- dd_i %>%
  filter(age_real <= max(age)) %>%
  group_by(age) %>%
  summarise(mean_lambda = round(mean(lambda, na.rm = T), 3),
            mean_w = round(mean(w, na.rm = T), 3),
            mean_ws = round(mean(ws, na.rm = T), 3))

# Plot
ggplot(data = melt(data.frame(mean_fitness, v_age[,2]), id.vars = c("age")),
       mapping = aes(x = age, y = value, linetype = variable)) +
  geom_line() +
  scale_y_continuous(limits = c(0,1.6))


## ---------------------------------------------------------------------------------------------
# The diameter measurements are all given in mm*10

# Calculate spottiness
dd_i$sbreast <- with(dd_i, 100*(pi*((diabreast/10)/2)^2*nbbreast)/(60*40))
dd_i$sbelly <- with(dd_i, 100*(pi*((diabelly/10)/2)^2*nbbelly)/(60*40))
dd_i$sflank <- with(dd_i, 100*(pi*((diaflank/10)/2)^2*nbflank)/(60*40))
dd_i$swing <- with(dd_i, 100*(pi*((diawing/10)/2)^2*nbwing)/(60*40))

# Center and SD-scale spottiness
dd_i <- dd_i %>%
  mutate(sbreast_c = scalevar(sbreast),
         sbelly_c = scalevar(sbelly),
         sflank_c = scalevar(sflank),
         swing_c = scalevar(swing)) %>%
  as.data.frame()

# Spottiness
dd_i <- dd_i %>%
  rowwise() %>%
  mutate(spottiness = weighted.mean(x = c(sbreast, sbelly, sflank, swing), w = c(1, 1, 2, 2)))

# Center and SD-scale spottiness
dd_i <- dd_i %>%
  mutate(spottiness_c = scalevar(spottiness)) %>%
  as.data.frame()

mean(dd_i$spottiness, na.rm=T)
se(dd_i$spottiness, na.rm=T)
sd(dd_i$spottiness, na.rm=T)



## ---------------------------------------------------------------------------------------------
# Center and SD-scale
dd_i$cbreast_c <- scalevar(dd_i$cbreast)
dd_i$cbelly_c <- scalevar(dd_i$cbelly)
dd_i$cflank_c <- scalevar(dd_i$cflank)
dd_i$cwing_c <- scalevar(dd_i$cwing)

# Mean colour
dd_i <- dd_i %>%
  rowwise() %>%
  mutate(colour = weighted.mean(x = c(cbreast, cbelly, cflank, cwing), w = c(1, 1, 2, 2))) %>%
  ungroup()

# Center and SD-scale colour
dd_i <- dd_i %>%
  mutate(colour_c = scalevar(colour)) %>%
  as.data.frame()

mean(dd_i$colour, na.rm=T)
se(dd_i$colour, na.rm=T)
sd(dd_i$colour, na.rm=T)



## ---------------------------------------------------------------------------------------------
# Correlations spottiness and colour
cor.prob(dd_i[!duplicated(dd_i$ring), c("sbreast", "sbelly", "sflank", "swing", "spottiness", "cbreast", "cbelly", "cflank", "cwing", "colour")], na.rm = TRUE, cor = TRUE, pval = TRUE)
# Sample size for correlations
n_distinct(dd_i[!duplicated(dd_i$ring) & !is.na(dd_i$spottiness) & !is.na(dd_i$colour),])


## ---------------------------------------------------------------------------------------------

# Spots
ggplot(data = subset(dd_i, !duplicated(ring)), mapping = aes(x = spottiness)) +
  geom_histogram()

# Colours
ggplot(data = subset(dd_i, !duplicated(ring)), mapping = aes(x = colour)) +
  geom_histogram()


## ----pheno geno-------------------------------------------------------------------------------
# Mean of all phenotypes for each genotype
mean_p <- dd_i %>%
  filter(!duplicated(ring) & !is.na(genotype)) %>%
  select(., genotype, sbreast, sbelly, sflank, swing, 
         cbreast, cbelly, cflank, cwing) %>%
  group_by(genotype) %>%
  summarise_all(.funs = funs(mean = mean, se = se), na.rm = TRUE)

# To long format
mean_p <- melt(mean_p)

# Add indicator of different group of phenotypes
mean_p$measure <- "se"
mean_p$measure[grepl(pattern = "_mean", x = mean_p$variable)] <- "mean"

# Add indicator of different group of phenotypes
mean_p$type <- "Spottiness"
mean_p$type[grepl(pattern = "^c", x = mean_p$variable)] <- "Colour"

# Remove dia, nb, c, mean and se from variable
mean_p$variable <- gsub(pattern = "^s", replace = "", x = mean_p$variable)
mean_p$variable <- gsub(pattern = "^c", replace = "", x = mean_p$variable)
mean_p$variable <- gsub(pattern = "_mean", replace = "", x = mean_p$variable)
mean_p$variable <- gsub(pattern = "_se", replace = "", x = mean_p$variable)

# Cast
mean_p <- dcast(mean_p, genotype+variable+type~measure, value.var = "value")

# Graph
ggplot(data = mean_p, mapping = aes(x = genotype, y = mean,
                                    ymin = mean+se*1.96, ymax = mean-se*1.96,
                                    group = variable, colour = variable)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.25) +
  geom_line(position = position_dodge(width = 0.3)) +
  facet_grid(type~., scale = "free", switch = "y") +
  xlab("Genotype") +
  ylab("Phenotype") +
  scale_colour_discrete(name = "Trait") +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))


## ----spottiness pheno geno--------------------------------------------------------------------
# Mean of compound phenotypes for each genotype
mean_sp <- dd_i %>%
  filter(!duplicated(ring) & !is.na(genotype)) %>%
  select(., genotype, colour, spottiness) %>%
  group_by(genotype) %>%
  summarise_all(.funs = funs(mean = mean, se = se), na.rm = TRUE)

# To long format
mean_sp <- melt(mean_sp)

# Add indicator of different group of phenotypes
mean_sp$measure <- "se"
mean_sp$measure[grepl(pattern = "_mean", x = mean_sp$variable)] <- "mean"

# Remove mean and se from variable
mean_sp$variable <- gsub(pattern = "_mean", replace = "", x = mean_sp$variable)
mean_sp$variable <- gsub(pattern = "_se", replace = "", x = mean_sp$variable)

# Cast
mean_sp <- dcast(mean_sp, genotype+variable~measure, value.var = "value")

# Rename variable
mean_sp$variable[mean_sp$variable == "spottiness"] <- "Spottiness"
mean_sp$variable[mean_sp$variable == "colour"] <- "Colour"

# Change order
mean_sp$variable <- factor(mean_sp$variable, levels = c("Spottiness", "Colour"))

# Rename genotype
levels(mean_sp$genotype) <- c("WW", "WR", "RR")

# Labels
mean_sp_labels <- data.frame(mean = c(4.5, -2), se = 1, genotype = "WW",
                           variable = c("Spottiness", "Colour"),
                           lab = c("A", "B"))

# Restructure data for plotting of raw data points
dd_i_long <- dd_i %>%
  filter(!duplicated(ring) & !is.na(genotype)) %>%
  select(., genotype, colour, spottiness) %>%
  melt(.) %>%
  mutate(variable=dplyr::recode(variable, spottiness="Spottiness", colour="Colour"),
         genotype=dplyr::recode(genotype, ww="WW", wr="WR", rr="RR"),
         se = 0,
         mean= 0)

# Graph spottiness
spot_plot <- ggplot(data = subset(mean_sp, variable == "Spottiness"),
       mapping = aes(x = genotype, y = mean, 
                     ymin = mean+se*1.96, ymax = mean-se*1.96, group = variable)) +
  geom_point(data=subset(dd_i_long, variable == "Spottiness"&!is.na(value)), aes(y=value, x = genotype), position = position_jitter(width = 0.15), size = 1, alpha=0.1) +
  geom_point(size = 1) +
  geom_errorbar(width = 0.1, size = 0.25) +
  geom_line(size = 0.25) +
  scale_y_continuous(breaks = c(14, 12, 10, 8, 6, 4, 2, 0), labels = c("14 (many)", 12, 10, 8, 6, 4, 2, "0 (few)")) +
  ylab("Spottiness") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        axis.title.x=element_blank(),
        plot.margin = margin(11/2, 11/2, 0, 
            11/2))

# Graph colour
col_plot <- ggplot(data = subset(mean_sp, variable == "Colour"),
       mapping = aes(x = genotype, y = mean, 
                     ymin = mean+se*1.96, ymax = mean-se*1.96, group = variable)) +
  geom_point(data=subset(dd_i_long, variable == "Colour"&!is.na(value)), aes(y=value, x = genotype), position = position_jitter(width = 0.15), size = 1, alpha=0.1) +
  geom_point(size = 1) +
  geom_errorbar(width = 0.1, size = 0.25) +
  geom_line(size = 0.25) +
  scale_y_continuous(breaks = c(-1, -2, -3, -4, -5, -6, -7, -8), labels = c("-1 (red)", -2, -3, -4, -5, -6, -7, "-8 (white)"), limits = c(-8,-1)) +
  xlab("MC1R genotype") +
  ylab("Colour") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        plot.margin = margin(0, 11/2, 11/2, 
            11/2))

# Combine plots
spot_col_plot <- plot_grid(spot_plot, col_plot, align = "v", ncol = 1, rel_heights = c(0.94,1))
spot_col_plot

# To file
ggsave(plot = spot_col_plot, filename = "output/pheno_geno_col_spot.pdf", height = 120, width = 89, units = "mm")
ggsave(plot = spot_col_plot, filename = "output/pheno_geno_col_spot.jpg", height = 120, width = 89, units = "mm", dpi = 300)


## ---------------------------------------------------------------------------------------------
# Overall spottiness
summary(aov(spottiness~genotype, data = subset(dd_i, !duplicated(ring))))
summary.lm(aov(spottiness~genotype, data = subset(dd_i, !duplicated(ring))))
# Multiple comparisons
TukeyHSD(aov(spottiness~genotype, data = subset(dd_i, !duplicated(ring))))
# Grouped genotypes
summary(aov(spottiness~genotype_g, data = subset(dd_i, !duplicated(ring))))
summary.lm(aov(spottiness~genotype_g, data = subset(dd_i, !duplicated(ring))))


## ---------------------------------------------------------------------------------------------
# Overall colouration
summary(aov(colour~genotype, data = subset(dd_i, !duplicated(ring))))
summary.lm(aov(colour~genotype, data = subset(dd_i, !duplicated(ring))))
# Multiple comparisons
TukeyHSD(aov(colour~genotype, data = subset(dd_i, !duplicated(ring))))


## ---------------------------------------------------------------------------------------------
with(subset(dd_i, !duplicated(ring)), table(genotype))
with(subset(dd_i, !duplicated(ring)), table(genotype, !is.na(lambda)))
with(subset(dd_i, !duplicated(ring)), table(genotype, !is.na(colour)))


## ---------------------------------------------------------------------------------------------
with(subset(dd_i, !duplicated(ring)), table(!is.na(genotype), year))


## ----weights----------------------------------------------------------------------------------
# Calculate weight
dd_i$c <- with(dd_i, ws2/lambda)

# Define weights to 1 when w2 and lambda both are 0.
dd_i$c[is.nan(dd_i$c)] <- 1

# View
mean(dd_i$c, na.rm = TRUE)
quantile(dd_i$c, na.rm = TRUE)


## ----melt-------------------------------------------------------------------------------------
# Melt
dd_il <- melt(select(dd_i, ring, age, sex, genotype, genotype_g,
                            sbreast_c, sbelly_c, sflank_c, swing_c,
                            spottiness_c, cbreast_c, cbelly_c,
                            cflank_c, cwing_c, colour_c,
                            ws2, v, c, popsize, year),
              id.vars = c("ring", "age", "sex", "genotype", "genotype_g", "ws2", "c", "v", "popsize", "year"),
              variable.name = "trait")

# View
head(dd_il)

# Number of traits
n_distinct(dd_il$trait)

# Mean center population size
dd_il <- dd_il %>%
  mutate(popsize_c = scale(popsize, scale = F))


## ---------------------------------------------------------------------------------------------

# Define parameters
param2 <- c("(Intercept)", "value", "value2", "popsize", "value:popsize")

# Create dataframe
pheno_trait2 <- data.frame(trait = rep(c("sbreast_c", "sbelly_c", "sflank_c", "swing_c",
                                        "spottiness_c", "cbreast_c", "cbelly_c",
                                        "cflank_c", "cwing_c", "colour_c"),
                                       each = length(param2)),
                          parameter = param2, estimate = NA, std_err = NA,
                          z_value = NA, p_value = NA, star = NA)


# List for models
pheno_trait_mod2 <- vector("list", length(unique(dd_il$trait)))
names(pheno_trait_mod2) <- unique(dd_il$trait)


# Selection on each trait
for (i in unique(dd_il$trait))
{
  
  # Estimate
  pheno_trait_mod2[[i]] <- glmer(ws2 ~ value + I(value^2) + I(popsize/100) +
                                  value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == i),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead"))
  # Extract
  pheno_trait2[pheno_trait2$trait == i, 3:6] <-
    round(summary(pheno_trait_mod2[[i]])$coef, 4)
  # Insert stars
  pheno_trait2[pheno_trait2$trait == i, 7] <-
    pval_stars(pheno_trait2[pheno_trait2$trait == i, 6])
  
}



## ---------------------------------------------------------------------------------------------
# Summary
subset(pheno_trait2, parameter == "value:popsize")
subset(pheno_trait2, trait == "colour_c")

# Colour
mod_f_pheno_colour2 <- glmer(ws2 ~ value + I(value^2) + I(popsize/100) +
                                  value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == "colour_c" & !is.na(trait)),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead"))
summary(mod_f_pheno_colour2)
# Confidence interval
mod_f_pheno_colour2_ci <- confint(mod_f_pheno_colour2)
# Transform from SD to variance
mod_f_pheno_colour2_ci[1,] <- mod_f_pheno_colour2_ci[1,]^2

# Backtransform estimates
# Collect in table
mod_f_pheno_colour2_coef_ci <- data.frame(Estimate = c(VarCorr(mod_f_pheno_colour)$year[1], coef(summary(mod_f_pheno_colour2))[,1]), mod_f_pheno_colour2_ci)
# Backtransform
mod_f_pheno_colour2_coef_ci[3,] <- mod_f_pheno_colour2_coef_ci[3,]/sd(dd_i$colour,na.rm=T)
mod_f_pheno_colour2_coef_ci[4,] <- mod_f_pheno_colour2_coef_ci[4,]/var(dd_i$colour,na.rm=T)
mod_f_pheno_colour2_coef_ci[5,] <- mod_f_pheno_colour2_coef_ci[5,]/100
mod_f_pheno_colour2_coef_ci[6,] <- (mod_f_pheno_colour2_coef_ci[6,]/100)/sd(dd_i$colour,na.rm=T)
# View
mod_f_pheno_colour2_coef_ci

# For print
# xtable::xtable(mod_f_pheno_colour2_coef_ci, digits = 4)


# Spottiness
mod_f_pheno_spottiness2 <- glmer(ws2 ~ value + I(value^2) + I(popsize/100) +
                                  value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == "spottiness_c" &
                                                !is.na(trait)),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "bobyqa"))
summary(mod_f_pheno_spottiness2)
# CI
mod_f_pheno_spottiness2_ci <- confint(mod_f_pheno_spottiness2)
# SD to Var
mod_f_pheno_spottiness2_ci[1,] <- mod_f_pheno_spottiness2_ci[1,]^2

# Backtransform estimates
# Collect in table
mod_f_pheno_spottiness2_coef_ci <- data.frame(Estimate = c(VarCorr(mod_f_pheno_spottiness2)$year[1], coef(summary(mod_f_pheno_spottiness2))[,1]), mod_f_pheno_spottiness2_ci)
# Backtransform
mod_f_pheno_spottiness2_coef_ci[3,] <- mod_f_pheno_spottiness2_coef_ci[3,]/sd(dd_i$spottiness,na.rm=T)
mod_f_pheno_spottiness2_coef_ci[4,] <- mod_f_pheno_spottiness2_coef_ci[4,]/var(dd_i$spottiness,na.rm=T)
mod_f_pheno_spottiness2_coef_ci[5,] <- mod_f_pheno_spottiness2_coef_ci[5,]/100
mod_f_pheno_spottiness2_coef_ci[6,] <- (mod_f_pheno_spottiness2_coef_ci[6,]/100)/sd(dd_i$spottiness,na.rm=T)
# View
mod_f_pheno_spottiness2_coef_ci

# For print
# xtable::xtable(mod_f_pheno_spottiness2_coef_ci, digits = 4)



## ---------------------------------------------------------------------------------------------
# Colour
test_chi_col <- drop1(glmer(ws2 ~ value + I(value^2) + I(popsize/100) +
                                  value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == "colour_c"),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead")),
                      test = "Chisq")

# Spottiness
test_chi_spo <- drop1(glmer(ws2 ~ value + I(value^2) + I(popsize/100) +
                                  value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == "spottiness_c"),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead")),
                      test = "Chisq")


## ---------------------------------------------------------------------------------------------

# Define parameters
param <- c("(Intercept)", "value", "popsize", "value:popsize")

# Create dataframe
pheno_trait <- data.frame(trait = rep(c("sbreast_c", "sbelly_c", "sflank_c", "swing_c",
                                        "spottiness_c", "cbreast_c", "cbelly_c",
                                        "cflank_c", "cwing_c", "colour_c"),
                                      each = length(param)),
                          parameter = param, estimate = NA, std_err = NA,
                          z_value = NA, p_value = NA, star = NA)


# List for models
pheno_trait_mod <- vector("list", length(unique(dd_il$trait)))
names(pheno_trait_mod) <- unique(dd_il$trait)


# Selection on each trait
for (i in unique(dd_il$trait))
{
  
  # Estimate
  pheno_trait_mod[[i]] <- glmer(ws2 ~ value + I(popsize/100) +
                                  value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == i),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead"))
  # Extract
  pheno_trait[pheno_trait$trait == i, 3:6] <-
    round(summary(pheno_trait_mod[[i]])$coef, 4)
  # Insert stars
  pheno_trait[pheno_trait$trait == i, 7] <-
    pval_stars(pheno_trait[pheno_trait$trait == i, 6])
  
}


# Summary
subset(pheno_trait, parameter == "value:popsize")
subset(pheno_trait, trait == "colour_c")



## ---------------------------------------------------------------------------------------------
# Fit model
mod_f_pheno_colour <- glmer(ws2 ~ value + I(popsize/100) + value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == "colour_c"),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead"))

# Summary
summary(mod_f_pheno_colour)
drop1(mod_f_pheno_colour, test = "Chisq")
confint(mod_f_pheno_colour)

# Effects of phenotype and density
plot(allEffects(mod_f_pheno_colour))


# Backtransform estimates
# Collect in table
mod_f_pheno_colour_coef <- fixef(mod_f_pheno_colour)
# Backtransform
mod_f_pheno_colour_coef[2] <- mod_f_pheno_colour_coef[2]/sd(dd_i$colour,na.rm=T)
mod_f_pheno_colour_coef[3] <- mod_f_pheno_colour_coef[3]/100
mod_f_pheno_colour_coef[4] <- (mod_f_pheno_colour_coef[4]/100)/sd(dd_i$colour,na.rm=T)
# View
mod_f_pheno_colour_coef



## ---------------------------------------------------------------------------------------------
# Fit model
mod_f_pheno_spot <- glmer(ws2 ~ value + I(popsize/100) + value:I(popsize/100) + (1|year),
                                data = subset(dd_il, trait == "spottiness_c"),
                                offset = log(c), weights = v*1/c, family = poisson,
                                control = glmerControl(optimizer = "Nelder_Mead"))

# Summary
summary(mod_f_pheno_spot)



## ---------------------------------------------------------------------------------------------
# Phenotype model
summary(mod_f_pheno_colour2)

# Mean growth rate (ind. fitness) (ln mean(W) - 1/2sigma2.e)
# From data
exp(log(mean(subset(dd_i, !is.na(colour)& !is.na(genotype))$lambda,na.rm=T)) - 0.5*0.1544)
# From observed Leslie matrix
exp(log(luv_f_obs$lambda) - 0.5*0.1544)

# Migration rate (-ln mean(W) + 1/2sigma2.e)
# From data
exp(-log(mean(subset(dd_i, !is.na(genotype)&!is.na(colour))$lambda,na.rm=T)) + 0.5*0.1544)
# From Leslie matrices
exp(-log(luv_f_K$lambda) + 0.5*0.1544)


## ---------------------------------------------------------------------------------------------
# Check mean and distribution of population sizes
mean(pop$popsize)
quantile(pop$popsize)

# Data frame with given population sizes and a range of phenotypic values for plotting of mean Malthusian fitness
res_plot <- data.frame(popsize = rep(c(50, 75, 100, 125, 150), each = 100), value = seq(min(dd_i$colour_c, na.rm = T), max(dd_i$colour_c, na.rm = T), length.out = 100), colour = seq(min(dd_i$colour, na.rm = T), max(dd_i$colour, na.rm = T), length.out = 100), m = NA)

# Predict mean Malthusian fitness
res_plot$m <- predict(pheno_trait_mod$colour_c,
                      newdata = res_plot,
                      re.form = NA, type = "link") + log(1.89) - 0.1544/2



## ---------------------------------------------------------------------------------------------
# Mean and range of colour
mean(dd_i$colour, na.rm=T)
range(dd_i$colour_c, na.rm=T)

# Plot
mzn_plot <- ggplot(data = res_plot, mapping = aes(x = colour, y = m)) +
  geom_line(aes(linetype = factor(popsize))) +
  geom_rug(data = subset(dd_i, !is.na(colour_c)&!duplicated(ring)), aes(x = colour, y = NULL), sides = "b", alpha = 0.1, size = 0.2) +
  ylab(expression("mean Malthusian fitness  ["~bar(m)(tilde(z),N)*"]")) +
  xlab(expression("Colour"~(tilde(z)))) +
  #xlab("") +
  scale_x_continuous(breaks = c(-7.5, -6, -4.5, -3, -1.5), labels = c("-7.5 (white)", "-6", "-4.5", "-3", "-1.5 (red)")) +
  scale_linetype(guide = F) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=-3, y=c(0.25, -0.07, -0.37, -0.64, -0.91)+log(1.747273), lab=c("N = 50", "N = 75", "N = 100", "N = 125", "N = 150"), popsize = 100), size = 10/(14/5), hjust = 0) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=-7.5, y=0.24+log(1.747273), lab="a", popsize = 100), fontface = "bold", size = 10/(14/5)) +
  theme_cowplot(font_size = 10, line_size = 0.25)
# View
mzn_plot

ggsave(plot = mzn_plot, filename = "output/pheno_mzn.pdf", height = 70, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Define s(z), gamma(z) and Q(z)
sz <- function(p, z, Vp, sig_e, mig = 1){
  p[1] + p[2]*z - sig_e/2 + log(mig)
  }
gz <- function(p, z, mig = 1){
  -p[3] - p[4]*z
  # Should have the negative of the parameters here as gamma is defined as the strength of density dependence m(z,N)=s(z)-gamma(z)N
}
qz <- function(p, z, Vp, sig_e, mig = 1){
  sz(p=p, z=z, Vp=Vp, sig_e=sig_e, mig=mig)/gz(p=p, z=z, mig=mig)
}


## ---------------------------------------------------------------------------------------------
# Check mean and distribution of population sizes
mean(pop$popsize)
quantile(pop$popsize)

# Data frame with population sizes and phenotypic values for plotting of Q(z)
res_plot_qz <- data.frame(value = seq(min(dd_i$colour_c, na.rm = T), max(dd_i$colour_c, na.rm = T), length.out = 100), colour = seq(min(dd_i$colour, na.rm = T), max(dd_i$colour, na.rm = T), length.out = 100), qz = NA)

# Parameter estimates from model
p_all <- c(fixef(mod_f_pheno_colour)/c(1,1,100,100),
           VarCorr(mod_f_pheno_colour)$year[1])
exp(-mean(sz(p = p_all[1:4], z = subset(dd_i, !is.na(colour))$colour_c, Vp = 0, sig_e = p_all[5], mig = 1)))

# Calculate sz, gz and Qz
res_plot_qz$qz <- qz(p = p_all[1:4], z = res_plot_qz$value, Vp = 0, sig_e = p_all[5], mig = 1.89)
res_plot_qz$sz <- sz(p = p_all[1:4], z = res_plot_qz$value, Vp = 0, sig_e = p_all[5], mig = 1.89)
res_plot_qz$gz <- gz(p = p_all[1:4], z = res_plot_qz$value)



## ---------------------------------------------------------------------------------------------
# 
qz_opt <- data.frame(optimize(f = function(z)qz(p = p_all[1:4], z = z, Vp = 0, sig_e = p_all[5], mig = 1.89), interval = c(min(res_plot_qz$value), max(res_plot_qz$value)), maximum = TRUE, tol = .Machine$double.eps^0.5))
names(qz_opt) <- c("z", "qz")
qz_opt$qmin <- 104
qz_opt$z <- -7.5


## ---------------------------------------------------------------------------------------------
# Plot
qz_plot <- ggplot(data = res_plot_qz, mapping = aes(x = colour, y = qz)) +
  geom_line() +
  geom_segment(mapping = aes(x = z, y = qmin, xend = z, yend = qz), data = qz_opt, linetype = 3) +
  geom_point(mapping = aes(x = z, y = qz), data = qz_opt, shape = 8) +
  geom_rug(data = subset(dd_i, !is.na(colour_c)&!duplicated(ring)), aes(x = colour, y = NULL), sides = "b", alpha = 0.1, size = 0.2) +
  ylab(expression(Q(tilde(z)))) +
  xlab(expression("Colour"~(tilde(z)))) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=-7.5, y=c(148), lab="b", popsize = 100), fontface = "bold", size = 10/(14/5)) +
  scale_y_continuous(limits = c(104,149), breaks = c(105, 115,  125, 135, 145), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-7.5, -6, -4.5, -3, -1.5), labels = c("-7.5 (white)", "-6", "-4.5", "-3", "-1.5 (red)")) +
  theme_cowplot(font_size = 10, line_size = 0.25)
# View
qz_plot

ggsave(plot= qz_plot, filename = "output/pheno_Qz.pdf", height = 70, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Data frame with a range of population sizes for plotting the selection gradient
sel_plot <- data.frame(popsize = seq(52, 187, length.out = 100), s = NA)

# Extract parameter values and create function for the selection gradient
mod_f_pheno_colour_coef

# Estimate selection gradient
sel_plot$s <- 0.179671877-0.001754904*(sel_plot$popsize)

# Plot
sz_plot <- ggplot(data = sel_plot, mapping = aes(x=popsize, y = s)) +
  geom_line() +
  #geom_point(data = subset(sel_plot, abs(s) == min(abs(s))), size = 2.5, shape = 8) +
  geom_rug(data = pop, aes(x = popsize, y = NULL), sides = "b", alpha = 0.3, size = 0.2) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=52, y=0.12, lab="c"), fontface = "bold", size = 10/(14/5)) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=80, y=0.06, lab="Red favoured"), size = 10/(14/5), angle = 321) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=145, y=-0.10, lab="White favoured"), size = 10/(14/5), angle = 321) +
  scale_y_continuous(limits = c(-0.15, 0.12), breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1)) +
  ylab("Selection gradient") +
  xlab("Population size (N)") +
  theme_cowplot(font_size = 10, line_size = 0.25)

sz_plot

ggsave(plot= sz_plot, filename = "output/pheno_sz.pdf", height = 70, width = 85, units = "mm")



## ---------------------------------------------------------------------------------------------
# Align plots
pheno_sel_plot <- plot_grid(mzn_plot, qz_plot, sz_plot, align = "v", ncol = 1)
# View
pheno_sel_plot

# Save to file
ggsave(plot = pheno_sel_plot, filename = "output/pheno_sel.pdf", height = 210, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Check mean and distribution of population sizes
mean(pop$popsize)
quantile(pop$popsize)

# Model
summary(pheno_trait_mod$colour_c)

# Predict values on xy grid
grid.lines <- 25
col.pred_original <- seq(-8, max(dd_i$colour,na.rm=T), length.out = grid.lines)
col.pred <- seq(min(dd_i$colour_c,na.rm=T)-0.345719, max(dd_i$colour_c,na.rm=T), length.out = grid.lines)
pop.pred <- seq(min(dd_i$popsize), max(dd_i$popsize), length.out = grid.lines)
xy <- expand.grid(value = col.pred, popsize = pop.pred)
m.pred <- matrix(predict(pheno_trait_mod$colour_c, newdata = xy, re.form = NA, type = "response"), 
                 nrow = grid.lines, ncol = grid.lines)

# Scatter plot with regression plane
# To file
pdf("output/pheno_lambda_3D.pdf", width = 12/2.54, height = 12/2.54, pointsize = 10)
par(mar = c(0,3,0,2))
scatter3D(x = dd_i$colour[!is.na(dd_i$colour)], y = dd_i$popsize[!is.na(dd_i$colour)], z = dd_i$lambda[!is.na(dd_i$colour)],
          pch = 20, cex = 2, theta = 145, phi = 10, ticktype = "detailed", nticks = 6, xlab = "Colour", ylab = "Population size (N)", zlab = "Individual fitness (lambda)",
          colkey = list(side = 4, plot = TRUE, length = 0.3, width = 0.5, dist = 0.01, shift = 0.1, addlines = F, col.clab = NULL, cex.clab = par("cex.lab"), side.clab = NULL, line.clab = NULL, adj.clab = NULL, font.clab = NULL), 
          alpha = 0.5, bty = "b", xlim = c(-8.2,-0.8), ylim = c(45,194),
          surf = list(x = col.pred_original, y = pop.pred, z = m.pred, facets = NA, col = "grey", lwd = 0.5))
dev.off()



## ---------------------------------------------------------------------------------------------
with(subset(dd_i, !is.na(colour)), mean(lambda, na.rm = T))

subset(dd_i, !is.na(colour)) %>%
  filter(!is.na(genotype)) %>%
  group_by(genotype) %>%
  summarise(mean(lambda, na.rm = T))


## ---------------------------------------------------------------------------------------------
# Names
names(dd_i)

# Fit complete
geno_all <- glmer(ws2 ~ genotype * I(popsize/100) + 
                     (1|year), data = subset(dd_i, !is.na(colour) & year > 1995),
                 offset = log(c), weights = v*1/c, family = poisson,
                 control = glmerControl(optimizer = "Nelder_Mead"))

# Summary model
summary(geno_all)

# CI
geno_all_f_ci <- confint(geno_all)
# SD to var
geno_all_f_ci[1,] <- geno_all_f_ci[1,]^2

# LRT
drop1(geno_all, test = "Chisq")

# For printing
# xtable(cbind(coefficients(summary(geno_all))[,1], geno_all_f_ci[-1,]))


## ---------------------------------------------------------------------------------------------
test_chi_gen <- drop1(geno_all, test = "Chisq")



## ---------------------------------------------------------------------------------------------
# Genotype model
summary(geno_all)
VarCorr(geno_all)$year[1]

# Mean growth rate (ind. fitness) (ln mean(W) - 1/2sigma2.e)
# From data
exp(log(mean(subset(dd_i, !is.na(colour)& !is.na(genotype))$lambda,na.rm=T)) - 0.5*VarCorr(geno_all)$year[1])
# From observed Leslie matrix
exp(log(luv_f_obs$lambda) - 0.5*VarCorr(geno_all)$year[1])

# Migration rate (-ln mean(W) + 1/2sigma2.e)
# From data
exp(-log(mean(subset(dd_i, !is.na(genotype)&!is.na(colour))$lambda,na.rm=T)) + 0.5*VarCorr(geno_all)$year[1])
# From Leslie matrices
exp(-log(luv_f_K$lambda) + 0.5*VarCorr(geno_all)$year[1])



## ---------------------------------------------------------------------------------------------
# Range and mean population size
with(dd_i, range(popsize, na.rm = TRUE))
with(pop, mean(popsize, na.rm = TRUE))

# Prediction data frame
geno_pred <- data.frame(genotype = rep(c("ww", "wr", "rr"), each = length(52:187)),
                        popsize = 52:187, c = 1)

# Predictions for whole population
geno_pred$lambda <- predict(geno_all, newdata = geno_pred, re.form = NA,
                    type = "link") - 0.5*VarCorr(geno_all)$year[1] + log(1.927207)


## ---------------------------------------------------------------------------------------------
geno_mpn_plot <- ggplot(data = geno_pred, mapping = aes(x = popsize, y = lambda)) +
  geom_line(aes(linetype = genotype)) +
  geom_rug(data = pop, aes(x = popsize, y = NULL), sides = "b", alpha = 0.3, size = 0.2) +
  xlab("Population size (N)") +
  ylab(expression("mean Malthusian fitness  ["~bar(m)(g[kl],N)*"]")) +
  scale_linetype_manual(values = c(3, 2, 1), guide = F) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=173, y=c(-0.75, -1.55, -2.7)+log(1.927207), lab=c("WW", "WR", "RR"), genotype = "ww"), size = 10/(14/5), hjust = 0) +
  #geom_vline(xintercept = 113.07) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=52, y=1.85+log(1.927207), lab="a", genotype = "ww"), fontface = "bold", size = 10/(14/5)) +
  theme_cowplot(font_size = 10, line_size = 0.25)
# View
geno_mpn_plot

# Save to file
ggsave(plot = geno_mpn_plot, filename = "output/geno_mpn.pdf", height = 70, width = 85, units = "mm")



## ---------------------------------------------------------------------------------------------
summary(geno_all)
fixef(geno_all)/c(1,1,1,100,100,100)
print(VarCorr(geno_all),comp="Variance")

# K(G_kl) - This is the deterministic carry capasity (K) for each genotype
# For WW
(0.2091237 + log(1.927207))/-(-0.005930454)
# For WR
(0.2091237+0.301741688 + log(1.927207))/-(-0.005930454-0.002738984)
# For RR
(0.2091237+3.107104834 + log(1.927207))/-(-0.005930454-0.030429690)

# Q(G_kl) - This is the long term carry capasity (K) for each genotype
# For WW
(0.2091237 + log(1.927207)-0.5*0.1934)/-(-0.005930454)
# For WR
(0.2091237+0.301741688 + log(1.927207)-0.5*0.1934)/-(-0.005930454-0.002738984)
# For RR
(0.2091237+3.107104834 + log(1.927207)-0.5*0.1934)/-(-0.005930454-0.030429690)

# Mean fitness given an allele frequency pp and assuming HW-equilibrium (at the mean population size ~106)
pp <- seq(0, 1, by = 0.1)
n <- 106
exp((0.2091237-0.005930454*n)*pp^2+(0.2091237+0.301741688)+(-0.005930454-0.002738984*n)*2*pp*(1-pp)+(0.2091237+3.107104834+(-0.005930454-0.030429690)*n)*(1-pp)^2 + log(1.927207) - 0.1934/2)
rm(pp, n)



## ---------------------------------------------------------------------------------------------
# Get parameter estimates from model
summary(geno_all)
par_geno_all <- fixef(geno_all)/c(1,1,1,100,100,100)

# Extract the parameter estimates and create function for the Q-function
qp <- function(p, par = par_geno_all, sig.e = 0.1934, mig = 1.927207){
  b1 <- par[1]
  b2 <- par[2]
  b3 <- par[3]
  # Should have the negative of the parameters for gamma here as gamma is defined as the strength of density dependence
  b4 <- -par[4]
  b5 <- -par[5]
  b6 <- -par[6]
  mi <- log(mig)
  r0 <- (b1+mi+(b1+b3+mi))/2
  g0 <- (b4+(b4+b6))/2
  ar <- b1+mi - r0
  ar <- ar
  ag <- b4 - g0
  ag <- ag
  dr <- (b1+b2)+mi - r0
  dg <- (b4+b5) - g0
  q <- 1-p
  # Q-function
  s_p <- ar*(p-q)+2*p*q*dr+r0-1/2*sig.e
  gamma_p <- ag*(p-q)+2*p*q*dg+g0
  s_p/gamma_p
}

# Q for an allele frequency pp and assuming HW-equilibrium
round(optimize(f = qp, interval = 0:1, maximum=T)$maximum, 3)


## ---------------------------------------------------------------------------------------------
# Data frame with allele frequencies for estimating Q(p)
res_plot_qp <- data.frame(p = seq(0, 1, length.out = 100), qp = NA)

# Estimate Q(p)
res_plot_qp$qp <- qp(par = par_geno_all, p = res_plot_qp$p, sig.e = 0.1934, mig = 1.927207)



## ---------------------------------------------------------------------------------------------
#
qp_opt <- data.frame(optimize(f = qp, interval = c(0,1), maximum = TRUE, tol = .Machine$double.eps^0.5))
names(qp_opt) <- c("p", "qp")
qp_opt$qmin <- 104

# Annual allele frequencies (removing duplicated individuals)
dd_allele_rv <- dd_i %>%
  filter(!is.na(genotype) & year > 1995) %>%
  group_by(year) %>%
  summarise(v_r = sum(allele_r*v, na.rm = TRUE),
            v_w = sum(allele_w*v, na.rm = TRUE),
            v_tot = v_r+v_w,
            freq_r = v_r/v_tot,
            freq_w = v_w/v_tot) %>%
  as.data.frame()
# View
head(dd_allele_rv)
tail(dd_allele_rv)
range(dd_allele_rv$freq_w)


## ---------------------------------------------------------------------------------------------
# Plot (with mean allele frequencies in the period)
qp_plot <- ggplot(data = res_plot_qp, mapping = aes(x = p, y = qp)) +
  geom_line() +
  geom_segment(mapping = aes(x = p, y = qmin, xend = p, yend = qp), data = qp_opt, linetype = 3) +
  geom_point(data = qp_opt, shape = 8) +
  geom_rug(data = dd_allele_rv, aes(x = freq_w, y = NULL), sides = "b", alpha = 0.5, size = 0.2) +
  ylab(expression(Q(tilde(p)))) +
  xlab(expression("Allele frequency"~(tilde(p)))) +
  scale_y_continuous(limits = c(104, 127), breaks = c(105, 115,  125), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.0, 1.0)) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=0.015, y=126, lab="b"), fontface = "bold", size = 10/(14/5)) +
  theme_cowplot(font_size = 10, line_size = 0.25)
# View
qp_plot

ggsave(plot=qp_plot, filename = "output/geno_Qp.pdf", height = 70, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Extract the parameter estimates and create function for the selection gradient
sel_p_fnc <- function(par, p, popsize, mig = 1.927207){
  b1 <- par[1]
  b2 <- par[2]
  b3 <- par[3]
  # Should have the negative of the parameters for gamma here as gamma is defined as the strength of density dependence
  b4 <- -par[4]
  b5 <- -par[5]
  b6 <- -par[6]
  mi <- log(mig)
  r0 <- (b1+mi+(b1+b3+mi))/2
  g0 <- (b4+(b4+b6))/2
  ar <- b1+mi - r0
  ag <- b4 - g0
  dr <- (b1+b2)+mi - r0
  dg <- (b4+b5) - g0
  q <- 1-p
  # Selection gradient
  2*(ar - ag*(popsize/100) + (q-p)*(dr - dg*(popsize/100)))
}
fixef(geno_all)


## ---------------------------------------------------------------------------------------------
# Range of population size
with(dd_i, range(popsize, na.rm = TRUE))

# Data frame with a range of population sizes and allele frequencies for plotting the selection gradient
sel_p_plot <- data.frame(popsize = seq(52, 187, length.out = 1000), p = rep(c(0, 0.25, 0.5, 0.75,1), each = 1000), s = NA)

# Calculate selection gradient
sel_p_plot$s <- sel_p_fnc(par = fixef(geno_all), p = sel_p_plot$p, popsize = sel_p_plot$popsize)

# Find the equilibriums
sel_p_plot0 <- sel_p_plot %>% group_by(p) %>% filter(abs(s) == min(abs(s)))
sel_p_plot0

# As factor
sel_p_plot$p <- factor(sel_p_plot$p)

# Plot
sp_plot <- ggplot(data = sel_p_plot, mapping = aes(x =popsize, y = s)) +
  geom_line(aes(linetype = p)) +
  #geom_point(data = sel_p_plot0, size = 2.5, shape = 8) +
  geom_rug(data = subset(pop, year > 1995), aes(x = popsize, y = NULL), sides = "b", alpha = 0.3, size = 0.2) +
  ylab("Selection gradient") +
  xlab("Population size (N)") +
  labs(linetype=expression("Allele frequency"~(tilde(p)))) +
  scale_linetype(guide = guide_legend(keywidth = 1.6)) +
  scale_x_continuous(limits=c(52,187)) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=85, y=-1.4, lab="Red favoured"), size = 10/(14/5), angle = 39) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=145, y=2.7, lab="White favoured"), size = 10/(14/5), angle = 39) +
  geom_text(mapping = aes(x=x,y=y,label=lab), data = data.frame(x=52, y=4.5, lab="c"), fontface = "bold", size = 10/(14/5)) +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(legend.position=c(.05,.65))
# View
sp_plot

# Save to file
ggsave(plot = sp_plot, filename = "output/geno_s_plot.pdf", height = 70, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Align plots
geno_sel_plot <- plot_grid(geno_mpn_plot, qp_plot, sp_plot, align = "v", ncol = 1)
# View
geno_sel_plot

# Save to file
ggsave(plot=geno_sel_plot, filename = "output/geno_sel.pdf", height = 210, width = 85, units = "mm")


## ---------------------------------------------------------------------------------------------
# Range and mean population size
with(dd_i, range(popsize, na.rm = TRUE))
with(pop, mean(popsize, na.rm = TRUE))

# Prediction data frame
geno_lambda_pred <- data.frame(genotype = rep(c("ww", "wr", "rr"), each = length(52:187)),
                        popsize = 52:187, c = 1)

# Predictions for whole population
geno_lambda_pred$lambda <- predict(geno_all, newdata = geno_pred, re.form = NA,
                    type = "response")

# Change genotype to capital letters
geno_lambda_pred$genotype <- toupper(geno_lambda_pred$genotype)

# Rename genotype
geno_lambda_pred <- rename(geno_lambda_pred, Genotype = genotype)

# Extract data points
dd_i_genotype <- subset(dd_i, !is.na(colour) & !is.na(genotype) & !is.na(lambda))
dd_i_genotype <- select(dd_i_genotype, genotype, lambda, popsize)
dd_i_genotype$genotype <- toupper(dd_i_genotype$genotype)
dd_i_genotype <- rename(dd_i_genotype, Genotype = genotype)



## ---------------------------------------------------------------------------------------------
geno_lambda_plot <- ggplot(data = geno_lambda_pred, mapping = aes(x = popsize, y = lambda)) +
  geom_line() +
  geom_point(data = dd_i_genotype, size = 2, alpha = 0.4) +
  xlab("Population size (N)") +
  ylab(expression("Individual fitness  ("*lambda*")")) +
  scale_linetype_manual(values = c(3, 2, 1), guide = F) +
  scale_shape_manual(values = c(21, 17, 16), guide = F) +
  scale_x_continuous(limits=c(51,188)) +
  scale_y_continuous(sec.axis = sec_axis(~ . + 0, name = "Genotype")) +
  facet_wrap(Genotype~., ncol = 1, strip.position = "right") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = rel(1), vjust = -1),
        axis.text = element_text(size = rel(1)),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y.right = element_text(vjust = -1))
# View
geno_lambda_plot

# Save to file
ggsave(plot = geno_lambda_plot, filename = "output/geno_lambda.pdf", height = 140, width = 85, units = "mm")



## ---------------------------------------------------------------------------------------------
summary(glmer(survival ~ colour_c + I(popsize/100) + colour_c:I(popsize/100) + (1|year),
                                data = dd_i,
                                family = binomial, control = glmerControl(optimizer = "bobyqa")))


## ---------------------------------------------------------------------------------------------
summary(glmer(recruits ~ colour_c + I(popsize/100) + colour_c:I(popsize/100) + factor(age) + (1|year),
                                data = dd_i,
                                family = poisson, control = glmerControl(optimizer = "bobyqa")))


## ---------------------------------------------------------------------------------------------
summary(glmer(survival ~ genotype * I(popsize/100) + factor(age) + 
                     (1|year),data = subset(dd_i, !is.na(colour)),family = binomial,
                 control = glmerControl(optimizer = "bobyqa")))


## ---------------------------------------------------------------------------------------------
summary(glmer(recruits ~ genotype * I(popsize/100) + factor(age) + 
                     (1|year),data = subset(dd_i, !is.na(colour)),family = poisson,
                 control = glmerControl(optimizer = "bobyqa")))



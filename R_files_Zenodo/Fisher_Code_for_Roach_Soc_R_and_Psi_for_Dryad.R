#START####

#Load packages& set working directory
library(dplyr)
library(lubridate)
library(ggplot2)
library(glmmTMB)
library(QGglmm)
library(car)
library(huxtable)
library(asnipe)
library(network)
library(sna)
library(MCMCglmm)

setwd("...")


#Load and manipulate data####

roach_video_raw = read.csv("roach_video_dryad.csv", header=T)

box_greg = roach_video_raw %>% 
  mutate(loc = as.numeric(loc)) %>%
  group_by(date, box) %>%
  summarise(mean_loc = mean(loc, na.rm=T),
            sum_loc = sum(loc, na.rm=T),
            n_obs = n())

roach_trials = read.csv("roach_trials_dryad.csv", header=T)
trial_info = read.csv("roach_trial_info_dryad.csv", header=T)
summary(trial_info) #there is a 2 for end.temp, checked lab book and it is a typo, should be 22
trial_info$end.temp[trial_info$end.temp == 2] = 22
roach_origins = read.csv("roach_origins_dryad.csv", header=T) %>%
  select(id, sex, origin, date_coll = date)

#need mass of ID. Create "week number " variable the use that to join so that individuals get their mass recording for that week

id_mass = roach_trials %>%
  mutate(date = dmy(date),
         week = week(date) ) %>%
  select(week, id = partner, mass_id2 = mass_partner)

roach_greg = roach_trials %>%
  mutate(date2 = dmy(date),
         week = week(date2) ) %>%
  left_join(id_mass, by = c("id", "week")) %>%
  left_join(box_greg, by = c("date", "box")) %>%
  left_join(trial_info, by = "date") %>% 
  left_join(roach_origins, by = "id") %>%
  filter(!(is.na(mean_loc))) %>% #remove trials where they did not record anything due to breach of barrier
  mutate(id = factor(id),
         mean_temp = (start.temp + end.temp) /2) %>%
  group_by(id) %>%
  mutate(mean_greg = mean(mean_loc, na.rm=T)) %>%
  arrange(mean_greg, by_group=T)

with(roach_greg, cor.test(mass_id, mass_id2)) #r = 0.994

#Fitting GLMMs####

#Fit sum_loc as response and offset of number of observations, poisson model
roach_greg_m1 = glmmTMB(sum_loc ~ sex*scale(mass_id2) + 
                                  sex*scale(mass_partner) + 
                                  scale(mean_temp) + 
                          (1 | id) +  (1 | partner) + (1 | date), 
                        offset = log(n_obs),
                        family = "poisson",
                        data = roach_greg)
summary(roach_greg_m1)

vcov(roach_greg_m1)
Anova(roach_greg_m1, type="III") 
#mass of both focal and partner both matter. Heavier indivs have lower scores i.e. more sociable
#sex borderline, males more sociable

hist(residuals(roach_greg_m1)) #good
plot(residuals(roach_greg_m1), fitted(roach_greg_m1)) # no directional correlation

#export model output as table

ranef_m1 = as.data.frame(VarCorr(roach_greg_m1)[[1]])[1,1:3]
hux_m1 = hux(summary(roach_greg_m1)$coef$cond) %>%
  add_columns(c("Intercept", "Sex (Male contrast)", "Body mass of focal", "Body mass of partner",
                "Temperature", "Sex : Body mass of focal", "Sex : Body mass of partner"), after = 0) %>%
  add_rows(c("Variable", "Estimate", "Standard Error","Chi-squared","P value"  ), after = 0)

hux_m1 = hux_m1   %>%
  set_contents(2:8, 4, Anova(roach_greg_m1, type= "III")$Chisq) %>%
  set_contents(2:8, 5, Anova(roach_greg_m1, type= "III")$Pr) %>%
  add_rows(rep(NA, ncol(hux_m1)), after = nrow(hux_m1))  %>%
  add_rows(rep(NA, ncol(hux_m1)), after = nrow(hux_m1))  %>%
  add_rows(rep(NA, ncol(hux_m1)), after = nrow(hux_m1))  %>%
  add_rows(rep(NA, ncol(hux_m1)), after = nrow(hux_m1)) %>%
  set_contents(9:12, 3, c("Random effect", "Individual", "Partner", "Date")) %>%
  set_contents(9:12, 4, c("Variance",ranef_m1[,1], ranef_m1[,2], ranef_m1[,3])) %>%
  set_number_format(3) %>%
  set_bottom_border(8, everywhere) 

quick_docx(hux_m1, file = "Sociability direct & indirect phenotypic effects Table 1.docx")

#compare to version with mean_loc and Gaussian to see how the former is preferable:
roach_greg_m2 = glmmTMB(mean_loc ~ sex*scale(mass_id2) +
                          sex*scale(mass_partner) + 
                          scale(mean_temp) + 
                          (1 | id) +  (1 | partner) +  (1 | date), 
                        data = roach_greg)

hist(residuals(roach_greg_m2)) #not acceptable
plot(residuals(roach_greg_m2), fitted(roach_greg_m2)) #some relationship, positive correlation

#compare residuals
par(mfrow=c(2,2))
hist(residuals(roach_greg_m1), xlab = "Residuals", main = "")
text(x = -7.8, y = 77, labels = "A", cex = 2, font = 2 ) #note precise location will depend on the size of your plotting window
plot(residuals(roach_greg_m1), fitted(roach_greg_m1), xlab = "Residuals", ylab = "Fitted values") 
hist(residuals(roach_greg_m2), xlab = "Residuals", main = "") #not acceptable
text(x = -3.8, y = 48, labels = "B", cex = 2, font = 2 ) #as above
plot(residuals(roach_greg_m2), fitted(roach_greg_m2), xlab = "Residuals", ylab = "Fitted values") 
par(mfrow=c(1,1))

#Plot sex-specific effect of mass on sociability
roach_partner_mass = ggplot(roach_greg, aes(x = mass_partner, 
                                            y = mean_loc,
                                            col = sex)) +
  geom_point() + #note clear separation due to sex
  geom_smooth(method = lm) +
  xlab("Partner body mass (g)") + 
  ylab ("Average distance from barrier") + 
  labs(col = "Sex") +
  theme_classic() + 
  scale_color_manual(values = c("coral", "dodgerblue"), 
                     labels = c("Females", "Males"))

#Repeatability & Repeatable influence####

#use package QGglmm to get variance components on observed scale (values in model summary are latent scale)

AI_var = as.numeric(VarCorr(roach_greg_m1)[[c('cond', 'id')]])
partner_var = as.numeric(VarCorr(roach_greg_m1)[[c('cond', 'partner')]])
date_var = as.numeric(VarCorr(roach_greg_m1)[[c('cond', 'date')]])
p_var = AI_var + partner_var + date_var
intercept = fixef(roach_greg_m1)[[1]]["(Intercept)"]

#for focal ID - Repeatability
QGicc(mu = intercept, var.comp = AI_var, var.p = p_var, model = "Poisson.log")

#for partner ID - Repeatable influence
QGicc(mu = intercept, var.comp = partner_var, var.p = p_var, model = "Poisson.log")


#Remove focal ID
roach_greg_m1b = glmmTMB(sum_loc ~ sex*scale(mass_id2) + 
                           sex*scale(mass_partner) + 
                          scale(mean_temp) + 
                           (1 | partner) + (1 | date), #removed id from model
                        offset = log(n_obs),
                        family = "poisson",
                        data = roach_greg)
AIC(roach_greg_m1, roach_greg_m1b)

AIC(roach_greg_m1b) - AIC(roach_greg_m1)  

#with ID is much better 


#remove partner ID
roach_greg_m1c = glmmTMB(sum_loc ~ sex*scale(mass_id2) +
                           sex*scale(mass_partner) + 
                           scale(mean_temp) + 
                           (1 | id) + (1 | date), #removed partner id from model
                         offset = log(n_obs),
                         family = "poisson",
                         data = roach_greg)
AIC(roach_greg_m1, roach_greg_m1c)

AIC(roach_greg_m1c) - AIC(roach_greg_m1)
#with partner ID is much better 

#try with mass of partner removed
roach_greg_m3 = glmmTMB(sum_loc ~ sex*scale(mass_id2) + 
                          scale(mean_temp) + 
                          (1 | id) +  (1 | partner) + (1 | date), 
                        offset = log(n_obs),
                        family = "poisson",
                        data = roach_greg)
summary(roach_greg_m3)

as.numeric(VarCorr(roach_greg_m1)[[c('cond', 'partner')]])
as.numeric(VarCorr(roach_greg_m3)[[c('cond', 'partner')]]) #an increase as expected, from 2.16 to 2.24


#Psi by sex####

#Get standardised estimates of psi for each sex
#Need to standardised both x and y therefore Poisson error distribution cannot be used
#Instead we use mean_loc with a Gaussian error distribution

roach_greg_f = roach_greg %>% 
  ungroup() %>%
  filter(sex == "f") %>%
  mutate(mean_loc_z = as.numeric(scale(mean_loc)))

roach_greg_psi_f = glmmTMB(mean_loc_z ~ scale(mass_id2) +
                             scale(mass_partner) + 
                           scale(mean_temp) + 
                           (1 | id) + (1 | partner) + (1 | date), 
                         data = roach_greg_f)

summary(roach_greg_psi_f)
fixef(roach_greg_psi_f)[[1]]["scale(mass_partner)"] #0.071


roach_greg_m = roach_greg %>%
  ungroup() %>%
  filter(sex == "m") %>%
  mutate(mean_loc_z = as.numeric(scale(mean_loc)))

roach_greg_psi_m = glmmTMB(mean_loc_z ~ scale(mass_id2) +
                             scale(mass_partner) + 
                             scale(mean_temp) + 
                             (1 | id) + (1 | partner) + (1 | date), 
                           data = roach_greg_m)

summary(roach_greg_psi_m)
fixef(roach_greg_psi_m)[[1]]["scale(mass_partner)"] # -0.129

#Social networks####

roach_assocs_raw = read.csv("roach_assocs_dryad.csv", header=T, stringsAsFactors = F)

roach_assocs = roach_assocs_raw %>% 
  filter(indiv != "??") %>% #remove indiviudals who could not be identified
  mutate(indiv_box = paste0(indiv,"_",box))

roach_assocs_box = roach_assocs %>% group_split(box)

to_gbi_to_net = function(edgelist){
  tryCatch({ get_network(
    get_group_by_individual(data.frame(edgelist[,c("indiv_box", "group","date")]),
                            data_format="individuals"),
    data_format = "GBI",
    association_index = "SRI")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\\n")})
}

box_mats = lapply(roach_assocs_box , to_gbi_to_net)

strength = as.data.frame (unlist(lapply(box_mats, rowSums)) )
colnames(strength)[1] = "stren"

strength$id_box = rownames(strength)

#Comparing sociability in dyadic trials and social networks####

roach_tags = read.csv("roach_tags_dryad.csv", header=T)

soc_comp = roach_tags %>%
  select(id, id_box) %>%
  mutate(id = as.factor(id)) %>%
  left_join(strength, by="id_box") %>% #stren is as expected, higher means more associates/bigger groups
  right_join(roach_greg, by = "id")  #note that sum_loc is inverse, so more sociable have lower values 

prior = list(R = list(V = diag(c(1,0.0001),2,2), nu = 2, fix = 2), #residual variance is a 2x2 (co)variance matrix, 
             #with the variance of the 2nd trait (strength) set to very small, as it does not vary within-individuals
             G = list(G1 = list(V = diag(2), nu = 2, #among-individual variance is a 6x6 (co)variance matrix
                                alpha.mu = rep(0,2),
                                alpha.V = diag(1,2,2))),
             B = list(mu=matrix(c(0,0,1,0),4), V = diag(4)*(10))) #this part of the prior is for fixed effects
diag(prior$B$V)[3:4]<-1e-9 #relationship with offset does not vary for either

roach_greg_mv = MCMCglmm(cbind(sum_loc, stren) ~
                           trait-1 + trait:log(n_obs), 
                         random =~ us(trait):id,
                         rcov =~ idh(trait):units,
                         family = c("poisson","gaussian"),
                         prior = prior,
                         nitt=550000,
                         burnin=50000,
                         thin=100,
                         verbose = F,
                         pr = F,
                         data = soc_comp)
summary(roach_greg_mv)
plot(roach_greg_mv)
posterior.mode(roach_greg_mv$VCV)
posterior.mode(roach_greg_mv$Sol)
cov2cor(matrix(posterior.mode(roach_greg_mv$VCV)[1:4],2,2)) 

AI_corr = roach_greg_mv$VCV[, "traitstren:traitsum_loc.id"] / 
  sqrt(roach_greg_mv$VCV[, "traitsum_loc:traitsum_loc.id"] * 
       roach_greg_mv$VCV[, "traitstren:traitstren.id"])

posterior.mode(AI_corr) 
mean(AI_corr)
HPDinterval(AI_corr, 0.95) 

#Repeat model to test convergence

roach_greg_mv2 = MCMCglmm(cbind(sum_loc, stren) ~
                            trait-1 + trait:log(n_obs), 
                          random =~ us(trait):id,
                          rcov =~ idh(trait):units,
                          family = c("poisson","gaussian"),
                          prior = prior,
                          nitt=550000,
                          burnin=50000,
                          thin=100,
                          verbose = F,
                          pr = F,
                          data = soc_comp)
roach_greg_mv3 = MCMCglmm(cbind(sum_loc, stren) ~
                            trait-1 + trait:log(n_obs), 
                          random =~ us(trait):id,
                          rcov =~ idh(trait):units,
                          family = c("poisson","gaussian"),
                          prior = prior,
                          nitt=550000,
                          burnin=50000,
                          thin=100,
                          verbose = F,
                          pr = F,
                          data = soc_comp)


gelman.diag(mcmc.list(roach_greg_mv$Sol,
                      roach_greg_mv2$Sol,
                      roach_greg_mv3$Sol), 
            confidence=0.95, autoburnin=T)

#MCMCglmm results to table####

#from: https://gkhajduk.github.io/2017-10-25-cleanMCMCglmm/
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  ## pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}


hux_mv = hux(clean.MCMC(roach_greg_mv)) %>% 
  add_columns(c("Post_mode",posterior.mode(roach_greg_mv$Sol),posterior.mode(roach_greg_mv$VCV)  ),
              after = 1) %>%
  add_rows(c("soc-stren AI corr", posterior.mode(AI_corr) ,
             mean(AI_corr),
             HPDinterval(AI_corr, 0.95), "", "","correlation"  ) , after = 9)  %>%
  set_contents(1,everywhere, c("Variable", "Posterior mode", "Posterior mean", "Lower 95% CrI", 
                              "Upper 95% CrI", "Effective sample size", "pMCMC", "Effect")) %>%
  set_contents(2:12,1,c("Pairwise sociability Intercept", "Network strength Intercept", "Pairwise sociability Offset*",
                          "Network strength Offset*", "Pairwise sociability Individual ID", 
                        "Pairwise sociability – Network strength Among-individual covariance",
                        "Pairwise sociability – Network strength Among-individual covariance","Network strength Individual ID",
                            "Pairwise sociability – Network strength Among-individual correlation",
                            "Pairwise sociability Residual", "Network strength Residual*")) %>%
  set_number_format(3)

quick_docx(hux_mv[c(1:6,9,7,10:12),], file = "Sociability direct & indirect phenotypic effects Table S1.docx")

#Plotting networks####

indiv_mean_greg = roach_greg %>%
  group_by(id) %>%
  summarise(mean_greg = max(mean_greg)) %>%
  mutate(id = as.integer(id),
         mean_greg = mean_greg*-1 + max(mean_greg)) %>% #so those with smallest scores now have biggest scores, and all above 0
  left_join(roach_tags, by = "id") %>%
  select(id_box , mean_greg)

l_to_network = function(l_adjacency) {
  tryCatch({ network(l_adjacency,
                     matrix.type="adjacency",
                     directed=F, 
                     ignore.eval=F, 
                     names.eval="weight")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\\n")})
}

box_nets = lapply(box_mats, l_to_network)

for (i in 1:length(unique(roach_assocs$box))) {
  tryCatch({
    box_nets[[i]]  %v% "greg" = indiv_mean_greg$mean_greg[match((box_nets[[i]] %v% "vertex.names"),
                                                                indiv_mean_greg$id_box)] 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\\n")})
}

box_colours = c("coral", "coral4","dodgerblue4", "dodgerblue")

nets_for_plot = list()

for (i in 1:length(unique(roach_assocs$box))) {
  nets_for_plot[[i]] = delete.edges(box_nets[[i]], which(box_nets[[i]] %e% "weight" < 0.12))
}

par(mfrow=c(2,2), mar = c(rep(1,4)))
for (i in 1:length(unique(roach_assocs$box))) {
  gplot(nets_for_plot[[i]], gmode="graph", 
        vertex.col = box_colours[i],
        vertex.cex = (box_nets[[i]] %v% "greg")/2.5, 
        vertex.border = "black",
        edge.col = "darkgrey", 
        edge.lwd = (nets_for_plot[[i]] %e% "weight")*2)
}
par(mfrow=c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))



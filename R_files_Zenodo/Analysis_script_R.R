# Read in data
FC_choices<-read.csv2("FC_choices.csv")
str(FC_choices)


# Change age as a factor and set "yearling" as the reference category
FC_choices$FC_F_age<-factor(FC_choices$FC_F_age)
FC_choices$FC_F_age<-relevel(FC_choices$FC_F_age, ref="yearling")
str(FC_choices)


# Take a subset of observations where great tits had started egg-laying or incubation
FC_choices_GT_nests<-FC_choices[FC_choices$GT_clutch>0,]
str(FC_choices_GT_nests)


#####

# Fitting a model
m.full<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch * FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices_GT_nests)
summary(m.full)

# Diagnostic plots
plot(m.full)

# Derivation of submodels and ranking them on the grounds of AICc
library(MuMIn)
options(na.action = "na.fail")
models.m.full<-dredge(global.model=m.full, beta="none", rank="AICc")

# Derive the 95% confidence set of models
avg.m.full<-model.avg(models.m.full, cumsum(weight) <= 0.95, fit=TRUE)
avg.m.full$msTable

# Remove those models that contain extra parameters compared to the best model 
# (i.e. more complex models) but are worse when measured in AICc, and all models
# that are >6 AICc units from the best model. Finally, average over 
# that set of models.

# First, fit the models to be included in averaging
m12347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch + FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices_GT_nests)
m2347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices_GT_nests)
m1347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch + FC_F_age,
   	family = binomial(link = logit), data = FC_choices_GT_nests)
m347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch,
   	family = binomial(link = logit), data = FC_choices_GT_nests)

# Model comparison and averaging
avg.models1<-model.avg(list(m12347, m2347, m1347, m347), fit=TRUE)
summary(avg.models1)

# Model m12347 is superior to all the other models in the set; evidence
# ratio is 5.46, so there is no need for model averaging and inferences
# can be based solely on model m12347.

summary(m12347)
plot(m12347)


#####

## Analysis of all data

M.full<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch * FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)
summary(M.full)

# Diagnostic plots
plot(M.full)


# Derivation of submodels and ranking them on the grounds of AICc
options(na.action = "na.fail")
models.M.full<-dredge(global.model=M.full, beta="none", rank="AICc")

# Derive the 95% confidence set of models
avg.M.full<-model.avg(models.M.full, cumsum(weight) <= 0.95, fit=TRUE)
avg.M.full$msTable

# Remove those models that contain extra parameters compared to the best model 
# (i.e. more complex models) but are worse when measured in AICc, and all models
# that are >6 AICc units from the best model. Finally, average over 
# that set of models.

# First, fit the models to be included in averaging
M12347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch + FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)
M2347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)
M12<-glm(formula = matching_choice ~ FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)
M1347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch + FC_F_age,
   	family = binomial(link = logit), data = FC_choices)
M123<-glm(formula = matching_choice ~ GT_clutch + FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)
M2<-glm(formula = matching_choice ~ FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)
M1<-glm(formula = matching_choice ~ FC_F_age,
   	family = binomial(link = logit), data = FC_choices)
M347<-glm(formula = matching_choice ~ prop_visible_eggs * GT_clutch,
   	family = binomial(link = logit), data = FC_choices)
M.null<-glm(formula = matching_choice ~ 1,
   	family = binomial(link = logit), data = FC_choices)
M124<-glm(formula = matching_choice ~ prop_visible_eggs + FC_F_age + FC_F_tarsus,
   	family = binomial(link = logit), data = FC_choices)

# Model comparison and averaging
avg.models2<-model.avg(list(M12347, M2347, M12, M1347, M123, M2, M1, M347, M.null, M124), fit=TRUE)
summary(avg.models2)

# Model M12347 is superior to all the other models in the set; evidence
# ratio is 5.2, so there is no need for model averaging and inferences
# can be based solely on model M12347.

summary(M12347)
plot(M12347)





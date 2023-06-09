########### Model 1
library(lme4)
test.data.model1=x### load data file for model 1
test.data.model2=x### load data file for model 2

contr=glmerControl(optCtrl = list(maxfun=10000000), calc.derivs=F, optimizer = "nloptwrap")

# Random slopes Model 1-1 (Both Species - Absolute Ranks and Relationships)
### Please note that the random slopes were selected after running the full model with full random intercept and slope structure (all combinations random intercept - fixed effects) and removing those slopes with below 0.00001 variance
reduced.model.slope.model1.1="(1 + z.rank.partner||focal) +
(1 + z.rank.focal + z.rank.partner + z.average.ddsi.partner + z.max.ddsi.partner||pot.partner) +
(1 + z.rank.focal + z.rank.partner + z.dsi.dyad + z.max.ddsi.partner + z.rank.focal * z.rank.partner||dyad) +
(1|gr.bout.index)"

# Model Specification Model 1-1
model1.1.terms="gr.initiated~z.rank.focal * (z.rank.partner + I(z.rank.partner^2)) * group + z.ddsi.dyad * group + repr.state.partner * group + z.average.ddsi.partner * group + z.max.ddsi.partner * group + prev.aggression * group +
sex.partner * group + sex.focal * sex.partner + offset(log(1/eff.gr.size))"
model1.1=as.formula(paste(model1.1.terms, reduced.model.slope.model1.1, sep="+"))
null1.1.terms="gr.initiated~sex.partner * group + sex.focal * sex.partner + offset(log(1/eff.gr.size))"
null1.1=as.formula(paste(null1.1.terms, reduced.model.slope.model1.1, sep="+"))

full.model1.1=glmer(model1.1, data=test.data.model1, family=binomial, control=contr)
null.model1.1=glmer(null1.1, data=test.data.model1, family=binomial, control=contr)
anova(full.model1.1, null.model1.1) #full null model comparison
drop1.model1.1=drop1(full.model1.1, test="Chisq") ###After, remove non-significant interactions from model and rerun to get final model


# Random slopes Model 1-2 (Both Species - Relative Ranks and Relationships)
### Please note that the random slopes were selected after running the full model with full random intercept and slope structure (all combinations random intercept - fixed effects) and removing those slopes with below 0.00001 variance
reduced.model.slope.model1.2="(1 + z.rel.rank.partner + z.rel.rank.focal * z.rel.rank.partner||focal) +
(1 + z.rel.rank.focal + I(z.rel.rank.partner^2) + z.average.ddsi.partner + z.max.ddsi.partner||pot.partner) +
(1 + z.rel.rank.focal + z.rel.rank.partner + z.rel.ddsi.dyad + z.max.ddsi.partner + z.rel.rank.focal * z.rel.rank.partner||dyad) +
(1 + z.rel.rank.partner + z.rel.ddsi.dyad||gr.bout.index)"

# Model Specification Model 1-2
model1.2.terms="gr.initiated~z.rel.rank.focal * (z.rel.rank.partner + I(z.rel.rank.partner^2)) * group + z.rel.ddsi.dyad * group + repr.state.partner * group + z.average.ddsi.partner * group + z.max.ddsi.partner * group + prev.aggression * group +
sex.partner * group + sex.focal * sex.partner + offset(log(1/eff.gr.size))"
model1.2=as.formula(paste(model1.2.terms, reduced.model.slope.model1.2, sep="+"))
null1.2.terms="gr.initiated~sex.partner * group + sex.focal * sex.partner + offset(log(1/eff.gr.size))"
null1.2=as.formula(paste(null1.2.terms, reduced.model.slope.model1.2, sep="+"))

full.model1.2=glmer(model1.2, data=test.data.model1, family=binomial, control=contr)
null.model1.2=glmer(null1.2, data=test.data.model1, family=binomial, control=contr)
anova(full.model1.2, null.model1.2)#full null model comparison
drop1.model1.2=drop1(full.model1.2, test="Chisq") ###After, remove non-significant interactions from model and rerun to get final model



####### Model 2

# Random slopes Model 2-1 (Only Chimpanzees - Absolute Ranks and Relationships)
### Please note that the random slopes were selected after running the full model with full random intercept and slope structure (all combinations random intercept - fixed effects) and removing those slopes with below 0.00001 variance
reduced.model.slope.model2.1="(1 + z.rank.partner||focal) +
(1 + z.rank.focal + z.rank.partner + z.average.ddsi.partner + z.max.ddsi.partner + I(z.rank.partner^2)||pot.partner) +
(1 + z.rank.focal + z.rank.partner + z.ddsi.dyad + z.max.ddsi.partner||dyad) +
(1|gr.bout.index)"

# Model Specification Model 2-1
model2.1.terms="gr.initiated~z.rank.focal * (z.rank.partner + I(z.rank.partner^2)) * group * sex.focal + z.ddsi.dyad * sex.focal + repr.state.partner * sex.focal + z.average.ddsi.partner * sex.focal + z.max.ddsi.partner * sex.focal + prev.aggression * group +
  sex.focal * sex.partner * group + offset(log(1/eff.gr.size))"
model2.1=as.formula(paste(model2.1.terms, reduced.model.slope.model2.1, sep="+"))
null2.1.terms="gr.initiated~sex.focal * sex.partner * group + offset(log(1/eff.gr.size))"
null2.1=as.formula(paste(null2.1.terms, reduced.model.slope.model2.1, sep="+"))

full.model2.1=glmer(model2.1, data=test.data.model2, family=binomial, control=contr)
null.model2.1=glmer(null2.1, data=test.data.model2, family=binomial, control=contr)
anova(full.model2.1, null.model2.1) #full null model comparison
drop1.model2.1=drop1(full.model2.1, test="Chisq") ###After, remove non-significant interactions from model and rerun to get final model


# Random slopes Model 2-2 (Only Chimpanzees - Relative Ranks and Relationships)
### Please note that the random slopes were selected after running the full model with full random intercept and slope structure (all combinations random intercept - fixed effects) and removing those slopes with below 0.00001 variance
reduced.model.slope.model2.2="(1 + z.rel.rank.partner||focal) +
(1 + z.average.ddsi.partner + z.max.ddsi.partner + I(z.rel.rank.partner^2)||pot.partner) +
(1 + z.rel.rank.focal + z.rel.rank.partner + z.rel.ddsi.dyad + z.rel.rank.focal*z.rel.rank.partner||dyad) +
(1 + z.rel.rank.partner + z.rel.ddsi.dyad||gr.bout.index)"

# Model Specification Model 2-2
model2.2.terms="gr.initiated~z.rel.rank.focal * (z.rel.rank.partner + I(z.rel.rank.partner^2)) * sex.focal + z.rel.ddsi.dyad * sex.focal + repr.state.partner * sex.focal + z.average.ddsi.partner * sex.focal + z.max.ddsi.partner * sex.focal + prev.aggression * group +
sex.focal * sex.partner * group + offset(log(1/eff.gr.size))"
model2.2=as.formula(paste(model2.2.terms, reduced.model.slope.model2.2, sep="+"))
null2.2.terms="gr.initiated~sex.focal * sex.partner * group + offset(log(1/eff.gr.size))"
null2.2=as.formula(paste(null2.2.terms, reduced.model.slope.model2.2, sep="+"))

full.model2.2=glmer(model2.2, data=test.data.model2, family=binomial, control=contr)
null.model2.2=glmer(null2.2, data=test.data.model2, family=binomial, control=contr)
anova(full.model2.2, null.model2.2) #full null model comparison
drop1.model2.2=drop1(full.model2.2, test="Chisq") ###After, remove non-significant interactions from model and rerun to get final model


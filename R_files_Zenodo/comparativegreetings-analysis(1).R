library(lme4); library(optimx)

# read data
xdata <- read.table("comparativegreetings.txt", header = TRUE, sep = "\\t")

# fit initial full model
full <- glmer(vocalized ~ species *isapproached *(Elodiff +affiliationindex +aud_friend +aud_highrank +focalElo)
             +(1|focalid) +(0+isapproached|focalid)
             +(0+affiliationindex|focalid) +(0+Elodiff|focalid)
             +(1|otherid) +(0+isapproached|otherid)
             , data = xdata, family = binomial,
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

# fit the null model
null <- glmer(vocalized ~ species
              +(1|focalid) +(0+isapproached|focalid)
              +(0+affiliationindex|focalid) +(0+Elodiff|focalid)
              +(1|otherid) +(0+isapproached|otherid)
              , data = xdata, family=binomial,
              control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                     optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

# compare full with null model
anova(null, full, test = "Chisq")

# model simplification, removing non-interpretable interaction terms
drop1(full, test="Chisq")

red <- update(full, .~. -species:isapproached:aud_friend)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:isapproached:Elodiff)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:isapproached:aud_highrank)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:isapproached:affiliationindex)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:isapproached:focalElo)
drop1(red, test="Chisq")
length(fixef(red))


# now continue with the two-way
red <- update(red, .~. -isapproached:aud_highrank)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:aud_friend)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -isapproached:aud_friend)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:aud_highrank)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:affiliationindex)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -isapproached:focalElo)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -isapproached:affiliationindex)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -species:focalElo)
drop1(red, test="Chisq")
length(fixef(red))

red <- update(red, .~. -isapproached:Elodiff)
drop1(red, test="Chisq")
length(fixef(red))


# final model renamed to 'res'
res <- red

# model parameters as reported in manuscript
coefficients(summary(full))
coefficients(summary(res))

VarCorr(full)
VarCorr(res)
VarCorr(null)


# test results as reported in manuscript

anova(null, full)

anova(update(res, .~. -focalElo),
      res)

anova(update(res, .~. -species:isapproached),
      res)

anova(update(res, .~. -species:Elodiff),
      res)

anova(update(res, .~. -affiliationindex),
      res)

anova(update(res, .~. -aud_friend),
      res)

anova(update(res, .~. -aud_highrank),
      res)



# analysis for 2011 study

dat2011 <- read.csv("dryaddata_2011.csv", h=T, stringsAsFactors = FALSE)

m.both <- lme(pis.f ~ age + sex + audio_ex + treat_no + stim*dist, random = ~1|uid, data = dat2011, method = "ML")
summary(m.both)
anova.lme(m.both, type = "marginal", adjustSigma = F)


# analysis for 2012 study

dat2012 <- read.csv("dryaddata_2012.csv", h = T, stringsAsFactors = FALSE)

m.both <- lme(pis.f ~ age + sex + audio_ex + treat_no + stim*noise_level, random = ~1|uid, data = dat2012, method = "ML")
summary(m.both)
anova.lme(m.both, type = "marginal", adjustSigma = F)

dat.lnoise <- dat2012[dat2012$noise_level =="low",]
dat.hnoise <- dat2012[dat2012$noise_level =="hi",]

dat.lnoise$stim = as.factor(dat.lnoise$stim)
dat.hnoise$stim = as.factor(dat.hnoise$stim)


# analysis for 2013 study

dat2013 <- read.csv("dryaddata_2013.csv", h = T, stringsAsFactors = FALSE)

m.both <- lme(pis.f ~ age + sex + audio_ex + treat_no + stim + stim*for_ret, random = ~1|uid, data = dat2013, method = "ML")
summary(m.both)
anova.lme(m.both, type = "marginal", adjustSigma = F)
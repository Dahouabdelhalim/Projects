#comparing daily rhythm between 'hard' & 'soft' weevils (using daytime movement)
#'act.day' represent data 'daily rhythm2dayonly.csv'
wilcox.test(act.score~act.type, mu=0, alt="two.sided",conf.level=0.95, paired=FALSE, conf.int=TRUE, exact=TRUE, correct=T, data=act.day)
#comparing daily rhythm between 'hard' & 'soft' weevils (using daytime movement)
#'act2' represent data 'daily rhythm2'
wilcox.test(act.score~act.type, mu=0, alt="two.sided",conf.level=0.95, paired=FALSE, conf.int=TRUE, exact=TRUE, correct=T, data=act2)


#comapring reflectance (hue, brightness and saturation)  between 'hard' & 'soft' weevils
hue.hard<-scan()
569.923
542.647
534.483
574.390
hue.soft<-scan()
556.620
543.803
551.460
536.697
t.test(hue.hard, hue.soft, var.equal=T)

brightness.hard<-scan()
135.335
126.174
129.858
124.521
brightness.soft<-scan()
121.899
123.211
124.056
130.833
t.test(brightness.hard, brightness.soft, var.equal=T)

saturation.hard<-scan()
0.635
0.698
0.680
0.644
saturation.soft<-scan()
0.694
0.665
0.650
0.680
t.test(saturation.hard, saturation.soft, var.equal=T)

#comparing survival rate: Fisher's exact test
#dt represent the data "data.csv"
pmat <- combn(levels(dt$trt), 2)
fisher.post <- vector("list", ncol(pmat))
for(i in 1:ncol(pmat)){
  thisInd <- 
    which(colnames(tabs) == pmat[1,i] | colnames(tabs) == pmat[2,i])
  fisher.post[[i]] <-
    fisher.test(tabs[, thisInd])
}
res <- sapply(fisher.post, function(x){
  list(
    odds.ratio = x$estimate,
    p.value = x$p.value
  )
})
colnames(res) <- apply(pmat, 2, function(x) paste(x, collapse = " vs "))

res <- rbind(res, p.value.adjusted = p.adjust(res["p.value", ], "bonferroni"))

print(res)

#comparing number of bites: Negative binomial regression
#dt represent the data "data.csv"
negps.1way <- glm.nb( bite ~ group + scale(elytra.length, scale = F)+scale(lizard.svl, scale = F), contrasts = list(group = contr.helmert), data = dt)
negps.1way.null<- glm.nb( bite ~ 1, data = dt)
drop1(negps.1way, ~., test="Chisq")

negps.1way.posthoc <- glht(negps.1way, linfct = mcp(group = "Tukey"))
summary(negps.1way.posthoc, test = adjusted("holm"))

#comparing prey handling duration: Survival analysis
#dt represent the data "data.csv"
m.1way <-survreg(Surv(time, event) ~ group + scale(elytra.length) + scale(lizard.svl),data = dt)
summary(m.1way)
Anova(m.1way)

m.1way.posthoc <- glht(m.1way, linfct = mcp(group = "Tukey"))
summary(m.1way.posthoc, test = adjusted("holm"))


#comapring hardness between 'hard' & 'soft' weevils
t.test(hardness ~ factor(treatment), data=pachy.hardness)

#linear regression on 'log10(force at failure)' ~ 'log10(body length)'
#hdt3 represent data 'hardness.bodylength.remove.hard.csv'
formu3 <- lm(log10(hardness) ~ log10(bodylength), data = hdt3)
summary(formu3)
AIC(formu3)
#linear regression on 'force at failure' ~ 'body length'
#hdt3 represent data 'hardness.bodylength.remove.hard.csv'
formu2<-lm(hardness~bodylength, data=hdt3)
summary(formu2)
AIC(formu2)




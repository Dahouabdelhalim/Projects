#Note: even parallelised, models 1.1 and 1.2 take several hours to run, so please contact author before attempting
library(lme4)
test.data.model1.rank=x### load data file for model 1.1
test.data.model1.ddsi=x### load data file for model 1.2

###Model 1.1
## in this model, "partner" refers to the groomer with the higher rank, and "competitor" to the groomer with the lower rank
#create random slopes for random effects of partner and competitor ID with those variables for which enough variation was present
model1.rank.slope="(0 + z.rank.partner|competitor) + (0 + z.rank.partner|pot.intervener) + (0 + z.rank.pot.intervener|competitor) + (0 + z.rank.pot.intervener|partner) + (0 + z.rank.competitor|competitor) + (0 + z.rank.competitor|partner) + (0 + z.rank.competitor|pot.intervener) + (0 + z.dsi.partner.pot.intervener|competitor) + (0 + z.dsi.partner.pot.intervener|partner) + (0 + z.dsi.partner.pot.intervener|pot.intervener) + (0 + z.dsi.competitor.pot.intervener|competitor) + (0 + z.dsi.competitor.pot.intervener|partner) + (0 + z.dsi.competitor.pot.intervener|pot.intervener) + (0 + z.dsi.subject.1.subject.2|competitor) + (0 + z.dsi.subject.1.subject.2|partner) + (0 + z.dsi.subject.1.subject.2|pot.intervener) + (0 + z.rank.partner*z.rank.pot.intervener|competitor) + (0 + z.rank.pot.intervener*z.rank.competitor|competitor) + (0 + z.rank.pot.intervener*z.rank.competitor|partner)"
model1.rank.full.terms<-"intervention.dummy~z.rank.partner*z.rank.pot.intervener*Group + z.rank.competitor*z.rank.pot.intervener*Group + z.dsi.partner.pot.intervener + z.dsi.competitor.pot.intervener + z.dsi.subject.1.subject.2 + sex.pot.intervener*Group + sex.competitor*Group + sex.partner*Group + repr.state.groomers +(1|partner)+(1|competitor)+(1|pot.intervener) + (1|gr.bout.index) + (1|dyad.partner.pot.intervener) + (1|dyad) + (1|dyad.competitor.pot.intervener) + offset(log(time.present))"
model1.rank.full.formula<-as.formula(paste(model1.rank.full.terms, model1.rank.slope, sep="+"))
model1.rank.null.terms<-"intervention.dummy~Group + z.dsi.partner.pot.intervener + z.dsi.competitor.pot.intervener + z.dsi.subject.1.subject.2 + sex.pot.intervener*Group + sex.competitor*Group + sex.partner*Group + repr.state.groomers +(1|partner)+(1|competitor)+(1|pot.intervener) + (1|gr.bout.index) + (1|dyad.partner.pot.intervener) + (1|dyad) + (1|dyad.competitor.pot.intervener) + offset(log(time.present))"
model1.rank.null.formula<-as.formula(paste(model1.rank.null.terms, model1.rank.slope, sep="+"))


contr=glmerControl(optCtrl = list(maxfun=10000000), calc.derivs=F, optimizer = "nloptwrap")
model1.rank.full=glmer(model1.rank.full.formula, data=test.data.model1.rank, family="binomial", control=contr)
model1.rank.null=glmer(model1.rank.null.formula, data=test.data.model1.rank, family="binomial", control=contr)

anova(model1.rank.null, model1.rank.full) #full null model comparison
model1.rank.drop1=drop1(model1.rank.full, test = "Chisq") ###After, remove non-significant interactions from model and rerun to get final model


### Model 1.2
## in this model, "partner" refers to the groomer with the higher DDSI with the bystander, and "competitor" to the groomer with the lower DDSI with the bystander
model1.ddsi.slope="(0 + z.dsi.partner.pot.intervener|competitor) + (0 + z.dsi.partner.pot.intervener|partner) + (0 + z.dsi.partner.pot.intervener|pot.intervener) + (0 + z.dsi.competitor.pot.intervener|competitor) + (0 + z.dsi.competitor.pot.intervener|partner) + (0 + z.dsi.competitor.pot.intervener|pot.intervener) + (0 + z.dsi.subject.1.subject.2|competitor) + (0 + z.dsi.subject.1.subject.2|partner) + (0 + z.dsi.subject.1.subject.2|pot.intervener) + (0 + z.dsi.partner.pot.intervener:z.dsi.competitor.pot.intervener|competitor) + (0 + z.dsi.partner.pot.intervener:z.dsi.competitor.pot.intervener|partner) + (0 + z.dsi.partner.pot.intervener:z.dsi.competitor.pot.intervener|pot.intervener)"
model1.ddsi.full.terms<-"intervention.dummy~z.dsi.partner.pot.intervener*z.dsi.competitor.pot.intervener*Group + z.dsi.subject.1.subject.2*Group + z.rank.pot.intervener + z.rank.partner + z.rank.competitor + sex.pot.intervener*Group + sex.competitor*Group + sex.partner*Group + repr.state.groomers +(1|partner)+(1|competitor)+(1|pot.intervener)+(1|gr.bout.index) + (1|dyad.partner.pot.intervener) + (1|dyad) + (1|dyad.competitor.pot.intervener) + offset(log(time.present))"
model1.ddsi.full.formula<-as.formula(paste(model1.ddsi.full.terms, model1.ddsi.slope, sep="+"))
model1.ddsi.null.terms<-"intervention.dummy~Group + z.rank.pot.intervener + z.rank.partner + z.rank.competitor + sex.pot.intervener*Group + sex.competitor*Group + sex.partner*Group + repr.state.groomers +(1|partner)+(1|competitor)+(1|pot.intervener)+(1|gr.bout.index) + (1|dyad.partner.pot.intervener) + (1|dyad) + (1|dyad.competitor.pot.intervener) + offset(log(time.present))"
model1.ddsi.null.formula<-as.formula(paste(model1.ddsi.null.terms, model1.ddsi.slope, sep="+"))

contr=glmerControl(optCtrl = list(maxfun=10000000), calc.derivs=F, optimizer = "nloptwrap")
model1.ddsi.full=glmer(model1.ddsi.full.formula, data=test.data.model1.ddsi, family="binomial", control=contr)
model1.ddsi.null=glmer(model1.ddsi.null.formula, data=test.data.model1.ddsi, family="binomial", control=contr)

anova(model1.ddsi.null, model1.ddsi.full) #full null model comparison
model1.ddsi.drop1=drop1(model1.ddsi.full, test = "Chisq") ###After, remove non-significant interactions from model and rerun to get final model


### Model 3
test.data.model3=x### load data file for model 3


###Model 3
## in this model, "groomer1" refers to the groomer targeted by the intervener, and "groomer2" to the other groomer
#create random slopes for random effects of target and competitor ID with those variables for which enough variation was present

#Model 3.1
model3.rank.slope="(0 + z.rank.groomer1|groomer2) + (0 + z.rank.groomer1|Intervener) + (0 + z.rank.intervener|groomer1) + (0 + z.rank.intervener|groomer2) + (0 + z.rank.groomer2|groomer1) + (0 + z.rank.groomer2|Intervener) + (0 + z.rank.groomer1:z.rank.intervener|groomer2) + (0 + z.rank.intervener:z.rank.groomer2|groomer1)"
model3.rank.full.terms<-"success~z.rank.groomer1*z.rank.intervener*Group + z.rank.groomer2*z.rank.intervener*Group + sex.groomer1 + sex.groomer2 + sex.Intervener + (1|groomer1)+(1|groomer2)+(1|Intervener)+(1|dyad)+(1|dyad.groomer1.intervener)+(1|dyad.groomer2.intervener)"
model3.rank.full.formula<-as.formula(paste(model3.rank.full.terms, model3.rank.slope, sep="+"))
model3.rank.null.terms<-"success~Group + sex.groomer1 + sex.groomer2 + sex.Intervener + (1|groomer1)+(1|groomer2)+(1|Intervener)+(1|dyad)+(1|dyad.groomer1.intervener)+(1|dyad.groomer2.intervener)"
model3.rank.null.formula<-as.formula(paste(model3.rank.null.terms, model3.rank.slope, sep="+"))

#GLMM

library(lme4)
glmer.control=glmerControl(optCtrl = list(maxfun=10000000), calc.derivs=F, optimizer = "nloptwrap")
model3.rank.full=glmer(model3.rank.full.formula, data=test.data.model3, family="binomial", control=glmer.control)
model3.rank.null=glmer(model3.rank.null.formula, data=test.data.model3, family="binomial", control=glmer.control)

anova(model3.rank.null, model3.rank.full)

#Model 3.2
model3.ddsi.slope="(0 + z.dsi.groomer1.groomer2|groomer1) + (0 + z.dsi.groomer1.groomer2|groomer2) + (0 + z.dsi.groomer1.groomer2|Intervener) + (0 + z.dsi.groomer1.intervener|groomer1) + (0 + z.dsi.groomer1.intervener|groomer2) + (0 + z.dsi.groomer1.intervener|Intervener) + (0 + z.dsi.groomer2.intervener|groomer1) + (0 + z.dsi.groomer2.intervener|groomer2) + (0 + z.dsi.groomer2.intervener|Intervener)"
model3.ddsi.full.terms<-"success~z.dsi.groomer1.groomer2*Group + z.dsi.groomer1.intervener*Group + z.dsi.groomer2.intervener*Group + sex.groomer1 + sex.groomer2 + sex.Intervener + (1|groomer1)+(1|groomer2)+(1|Intervener)+(1|dyad)+(1|dyad.groomer1.intervener)+(1|dyad.groomer2.intervener)"
model3.ddsi.full.formula<-as.formula(paste(model3.ddsi.full.terms, model3.ddsi.slope, sep="+"))
model3.ddsi.null.terms<-"success~Group + sex.groomer1 + sex.groomer2 + sex.Intervener + (1|groomer1)+(1|groomer2)+(1|Intervener)+(1|dyad)+(1|dyad.groomer1.intervener)+(1|dyad.groomer2.intervener)"
model3.ddsi.null.formula<-as.formula(paste(model3.ddsi.null.terms, model3.ddsi.slope, sep="+"))


library(lme4)
glmer.control=glmerControl(optCtrl = list(maxfun=10000000), calc.derivs=F, optimizer = "nloptwrap")
model3.ddsi.full=glmer(model3.ddsi.full.formula, data=test.data.model3, family="binomial", control=glmer.control)
model3.ddsi.null=glmer(model3.ddsi.null.formula, data=test.data.model3, family="binomial", control=glmer.control)

anova(model3.ddsi.null, model3.ddsi.full)




##### Model 2: Repeated Measure Design
###For this model, every intervention is represented by two entries: one for the target, one for the non-target. To test impact of variables on selection,
###we repeatedly randomly select either the target or the non-target to represent the intervention and run the full and null model and the drop1 function

test.data.model2=x ### load dataset for model 2
## in this model, "target" refers to the groomer that was first groomed by the intervener, and "competitor" to the groomer who was not chosen
#create random slopes for random effects of target and competitor ID with those variables for which enough variation was present
model2.slope="(0 + z.rank.target|competitor) + (0 + z.rank.target|Intervener) + (0 + z.dsi.target.intervener|target) + (0 + z.dsi.target.intervener|competitor) + (0 + z.dsi.target.intervener|Intervener)"
model2.full.terms<-"choice~z.rank.target * z.rank.intervener * Group + z.dsi.target.intervener * Group + sex.target + sex.competitor + sex.intervener + z.rank.intervener + (1|target) + (1|Intervener) + (1|competitor) + (1|dyad) + (1|dyad.competitor.intervener) + (1|dyad.target.intervener)"
model2.full.formula<-as.formula(paste(model2.full.terms, model2.slope, sep="+"))
model2.null.terms<-"choice~Group + sex.target + sex.competitor + sex.intervener + z.rank.intervener + (1|target) + (1|Intervener) + (1|competitor) + (1|dyad) + (1|dyad.competitor.intervener) + (1|dyad.target.intervener)"
model2.null.formula<-as.formula(paste(model2.null.terms, model2.slope, sep="+"))

library(lme4)
glmer.control=glmerControl(optCtrl = list(maxfun=10000000), calc.derivs=F, optimizer = "nloptwrap")
full.coefficients=list()
full.null=list()

for(i in 1:1000){
  data.set=c()
  for(j in 1:length(unique(test.data.model2$bout))){
    sel.data=subset(test.data.model2, bout==j)
    sel.row=sel.data[sample(nrow(sel.data), 1), ]
    data.set=rbind(data.set, sel.row)
  }
  full=glmer(model2.full.formula, data=data.set, family="binomial", control=glmer.control)
  null=glmer(model2.null.formula, data=data.set, family="binomial", control=glmer.control)
  full.coefficients[[i]]=as.data.frame(summary(full)$coefficients)
  full.null[[i]]=as.data.frame(anova(null, full, test="Chisq"))
}

mean.null.pvalue=rowMeans(sapply(full.null, function(x) x[,"Pr(>Chisq)"]))
mean.null.df=rowMeans(sapply(full.null, function(x) x[,"Chi Df"]))
mean.null.chisq=rowMeans(sapply(full.null, function(x) x[,"Chisq"]))
factors.null=row.names(full.null[[1]])
result.null=data.frame(factors.null, mean.null.chisq, mean.null.df, mean.null.pvalue)
names(result.null)=c("factor","Chi Sq","df","p-value")

mean.estimate=rowMeans(sapply(full.coefficients, function(x) x$Estimate))
mean.StdErr=rowMeans(sapply(full.coefficients, function(x) x[,"Std. Error"]))
mean.zvalue=rowMeans(sapply(full.coefficients, function(x) x[,"z value"]))
mean.pvalue=rowMeans(sapply(full.coefficients, function(x) x[,"Pr(>|z|)"]))
factors=row.names(full.coefficients[[1]])
result=data.frame(factors, mean.estimate, mean.StdErr, mean.zvalue, mean.pvalue)
names(result)=c("factor","Estimate","StdError","z-value","p-value")

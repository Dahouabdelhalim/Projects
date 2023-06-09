##Read in the Gemsnakes tree
read.tree("Tree_Dated_Point_PL_Astral_topology_130.txt" )->MG
MG<- phangorn::nnls.tree(cophenetic(MG), MG, rooted = TRUE)
drop.tip(MG,MG$tip.label[124:130])->MG2
drop.tip(MG2,MG2$tip.label[1:14])->MG3
MG3->MG4


###run tess
TESSF(MG4,0.96,"Gemsnakes_TESS2", 2)->outMG
outMG$my_model->my_model
###get the posterior
read_Tess_post(outMG$tess.output$`speciation rates`,outMG$tess.output$`extinction rates`,MG4)->posterior



#####EXTINCTION MODELS

###1) temporally autocorelated extinction rates, check rate median to make sure it work for your group

###make sure there is overlap with your estimates change max.rate is best

times <- seq(0, max(my_model$time), length.out = 100)
extinction_rate_samples <- function() {
  sample.basic.models( times = times,
                       model="MRF",
                       max.rate=1)
}


test1<- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="extinction",
  sample.extinction.rates=extinction_rate_samples)

plot(test1)

X1 <- summarize.trends(test1,
                      threshold = 0.02)
plot(X1)


samples_temp_auto_ext <- sample.congruence.class.posterior(posterior, 
                                             num.samples = 100,
                                             rate.type = "extinction",
                                             model = "MRF",
                                             MRF.type = "GMRF",
                                             max.rate = 1)

model_congruence_PP(posterior, samples_temp_auto_ext,thresh=0.02)->out1


####2) linear decreasing


times <- seq(0, max(my_model$time), length.out = 100)
extinction_rate_samples <- function() {
  sample.basic.models( times = times,
                       model="linear",
                       direction="decrease",
                       max.rate=1,
                       fc.mean=2)
}
test2 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="extinction",
  sample.extinction.rates=extinction_rate_samples)

plot(test2)


X2<- summarize.trends(test2,
                       threshold = 0.02)
plot(X2)


samples_lin_decr_ext <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "extinction",
                                                      direction="decrease",
                                                      model = "linear",
                                                      MRF.type = "GMRF",
                                                      max.rate = 1,
                                                      fc.mean=2)

model_congruence_PP(posterior, samples_lin_decr_ext,thresh=0.02)->out2


####3) exponentially decreasing

times <- seq(0, max(my_model$time), length.out = 100)
extinction_rate_samples <- function() {
  sample.basic.models( times = times,
                       model="exponential",
                       direction="decrease",
                       max.rate=1,
                       fc.mean=2)
}
test3 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="extinction",
  sample.extinction.rates=extinction_rate_samples)

plot(test3)
X3<- summarize.trends(test3,
                      threshold = 0.02)
plot(X3)


samples_exp_decr_ext <- sample.congruence.class.posterior(posterior, 
                                             num.samples = 100,
                                             rate.type = "extinction",
                                             direction="decrease",
                                             model = "exponential",
                                             MRF.type = "GMRF",
                                             max.rate = 1,
                                             fc.mean=2)


model_congruence_PP(posterior, samples_exp_decr_ext,thresh=0.02)->out3


####4)Linearly increasing
times <- seq(0, max(my_model$time), length.out = 100)
extinction_rate_samples <- function() {
  sample.basic.models( times = times,
                       direction="increase",
                       model="linear",
                       monotonic=TRUE,
                       MRF.type="GMRF",
                       max.rate=1,
                       fc.mean=2)
}
test4 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="extinction",
  sample.extinction.rates=extinction_rate_samples)


plot(test4)

X4<- summarize.trends(test4,
                      threshold = 0.02)
plot(X4)

samples_lin_incr_ext <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "extinction",
                                                      direction="increase",
                                                      model = "linear",
                                                      monotonic=TRUE,
                                                      MRF.type = "GMRF",
                                                      max.rate = 1,
                                                      fc.mean=2)

model_congruence_PP(posterior, samples_lin_incr_ext,thresh=0.02)->out4

####5)exponentially increasing
times <- seq(0, max(my_model$time), length.out = 100)
extinction_rate_samples <- function() {
  sample.basic.models( times = times,
                       direction="increase",
                       model="exponential",
                       monotonic=TRUE,
                       MRF.type="GMRF",
                       max.rate=1,
                       fc.mean=2)
}
test5 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="extinction",
  sample.extinction.rates=extinction_rate_samples)


plot(test5)

X5<- summarize.trends(test5,
                      threshold = 0.02)
plot(X5)

samples_exp_incr_ext <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "extinction",
                                                      direction="increase",
                                                      model = "exponential",
                                                      monotonic=TRUE,
                                                      MRF.type = "GMRF",
                                                      max.rate = 1,
                                                      fc.mean=2)



model_congruence_PP(posterior, samples_exp_incr_ext,thresh=0.02)->out5




#######SPECIATION MODELS



###6) Temporally autocorelated speciation rates, check rate median to make sure it work for your group

###make sure there is overlap with your estimates change max.rate is best
times <- seq(0, max(my_model$time), length.out = 100)
speciation_rate_samples <- function() {
  sample.basic.models( times = times,
                       rate0 = my_model$lambda(0.0),
                       model = "MRF",
                       MRF.type = "GMRF",
                       max.rate = 1.0)
}


test6 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="speciation",
  sample.speciation.rates = speciation_rate_samples)

plot(test6)
X6 <- summarize.trends(test6,
                       threshold = 0.02)
plot(X6)




samples_temp_auto_sp <- sample.congruence.class.posterior(posterior, 
                                                       num.samples = 100,
                                                       rate.type = "speciation",
                                                       rate0.median = 0.05,
                                                       model = "MRF",
                                                       MRF.type = "GMRF",
                                                       max.rate = 1)

model_congruence_PP(posterior, samples_temp_auto_sp,thresh=0.02)->out6

####7) linearly decreasing speciation

times <- seq(0, max(my_model$time), length.out = 100)
speciation_rate_samples <- function() {
  sample.basic.models( times = times,
                       model="linear",
                       rate0 = my_model$lambda(0.0),
                       direction="decrease",
                       max.rate=1,
                       fc.mean=1)
}
test7 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="speciation",
  sample.speciation.rates=speciation_rate_samples)

plot(test7)

X7<- summarize.trends(test7,
                       threshold = 0.02)
plot(X7)


samples_lin_decr_spec <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "speciation",
                                                      direction="decrease",
                                                      model = "linear",
                                                     noisy=FALSE,
                                                      max.rate = 1,
                                                      fc.mean=2)

model_congruence_PP(posterior, samples_lin_decr_spec,thresh=0.02)->out7

####8) exponentially decreasing speciation rates

times <- seq(0, max(my_model$time), length.out = 100)
speciation_rate_samples <- function() {
  sample.basic.models( times = times,
                       model="exponential",
                       rate0 = my_model$lambda(0.0),
                       direction="decrease",
                       max.rate=1,
                       fc.mean=2)
}
test8 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="speciation",
  sample.speciation.rates=speciation_rate_samples)

plot(test2)
X8 <- summarize.trends(test8,
                      threshold = 0.02)
plot(X8)


samples_exp_decr_speciation<- sample.congruence.class.posterior(posterior, 
                                                                num.samples = 100,
                                                                rate.type = "extinction",
                                                                direction="decrease",
                                                                model = "exponential",
                                                                noisy=FALSE,
                                                                max.rate = 1,
                                                                fc.mean=2)



model_congruence_PP(posterior, samples_exp_decr_speciation,thresh=0.02)->out8

####9)Linearly increasing

times <- seq(0, max(my_model$time), length.out = 100)
speciation_rate_samples <- function() {
  sample.basic.models( times = times,
                       direction="increase",
                       model="linear",
                       rate0 = my_model$lambda(0.0),
                       monotonic=TRUE,
                       MRF.type="GMRF",
                       max.rate=.5,
                       fc.mean=2)
}
test9 <- sample.congruence.class(my_model,num.samples=100,rate.type="speciation",sample.speciation.rates=speciation_rate_samples)

plot(test9)

X9 <- summarize.trends(test9,
                      threshold = 0.02)
plot(X9)

samples_lin_incr <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "speciation",
                                                      direction="increase",
                                                      model = "linear",
                                                      monotonic=TRUE,
                                                      MRF.type = "GMRF",
                                                      max.rate = .5,
                                                      fc.mean=2)


model_congruence_PP(posterior, samples_lin_incr,thresh=0.02)->out9

####10)exponentially increasing speciation
times <- seq(0, max(my_model$time), length.out = 100)
speciation_rate_samples <- function() {
  sample.basic.models( times = times,
                       direction="increase",
                       model="exponential",
                       monotonic=TRUE,
                       rate0 = my_model$lambda(0.0),
                       MRF.type="GMRF",
                       max.rate=.5,
                       fc.mean=2)
}
test10 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="speciation",
  sample.speciation.rates=speciation_rate_samples)


plot(test10)
X10 <- summarize.trends(test10,
                       threshold = 0.02)
plot(X10)

samples_exp_incr <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "speciation",
                                                      direction="increase",
                                                      model = "exponential",
                                                      monotonic=TRUE,
                                                      MRF.type = "GMRF",
                                                      max.rate = .5,
                                                      fc.mean=2)

model_congruence_PP(posterior, samples_exp_incr,thresh=0.02)->out10


####EPISODIC

times <- seq(0, max(my_model$time), length.out = 100)
extinction_rate_samples <- function() {
  sample.basic.models( times = times,
                       direction="increase",
                       model="episodic2",
                       noisy=FALSE,
                       max.rate=1,
                       fc.mean=2)
}
test11 <- sample.congruence.class(
  my_model,
  num.samples=100,
  rate.type="extinction",
  sample.extinction.rates=extinction_rate_samples)

plot(test11)
X11 <- summarize.trends(test11,
                       threshold = 0.02)
plot(X11)

samples_episodic_extinction <- sample.congruence.class.posterior(posterior, 
                                                      num.samples = 100,
                                                      rate.type = "extinction",
                                                      direction="increase",
                                                      monotonic=TRUE,
                                                      model="episodic2",
                                                      MRF.type = "GMRF",
                                                      max.rate = 1,
                                                      fc.mean=2)

model_congruence_PP(posterior, samples_episodic_extinction,thresh=0.02)->out11

namesProb<-c("PP_Model_1", "PP_Model_2", "PP_Model_3", "PP_Model_4", "PP_Model_5", "PP_Model_6", "PP_Model_8", "PP_Model_9", "PP_Model_10", "PP_Model_11")
postProb<-list(out1$plot,out2$plot,out3$plot,out4$plot,out5$plot,out6$plot,out7$plot,out8$plot,out9$plot,out10$plot,out11$plot)

plotter<-function(name, pp) 
{
pdf(name)
print(pp)
dev.off()
}


plotter(namesProb[1], postProb[1])

lapply(c(1:length(namesProb)), function(x) plotter(namesProb[x], postProb[x]))

namesSims<-c("Sims_Model_1","Sims_Model_2","Sims_Model_3","Sims_Model_4","Sims_Model_5","Sims_Model_6","Sims_Model_7","Sims_Model_8","Sims_Model_9","Sims_Model_10","Sims_Model_11")
sims<-list(test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11)

plotter2<-function(name, simu)
{pdf(name)
  print(plot(simu))
  dev.off()
}
medNames<-c("median_model_1","median_model_2","median_model_3","median_model_4","median_model_5","median_model_6","median_model_7","median_model_8","median_model_9","median_model_10","median_model_11")
med<-list(X1, X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)

lapply(c(1:length(namesSims)), function(x) plotter2(namesSims[x], sims[[x]]))

plotter3<-function(name, medi)
{pdf(name)
  print(plot(medi))
  dev.off()
}

lapply(c(1:length(medNames)), function(x) plotter3(medNames[x], med[[x]]))

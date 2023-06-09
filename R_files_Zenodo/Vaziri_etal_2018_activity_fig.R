require(ggplot2)
require(dplyr)

getwd()

obs<-read.csv("Vaziri_etal_2018_data_for_activity_figure.csv")
#find avg prop active for each treatment group by half hour since injection
head(obs)

unique(obs$bird_id_round)

act_av<-tapply(obs$propact, list(obs$comb_trt, obs$half_hr_since_inj),  mean, na.rm=T)
head(act_av)
dim(act_av)

# #find ses by treatment group by half hours  since injection
 se<-function(x) { sqrt((var(na.omit(x))/length(na.omit(x))))}   #calculates std error
 act_se<-tapply(obs$propact, list(obs$comb_trt, obs$half_hr_since_inj),  se)
 act_se<-as.vector(act_se)
 length(act_se)

#take row names from act_av and repeat the same number of times as there are columns in act-av
trt<-rep(dimnames(act_av)[[1]], 48)
length(trt)

#take column names from act_av (time interval in half hours since injection) and repeat the same number of times are there are rows (2)
hours_since_inj<-rep(dimnames(act_av)[[2]], each = nrow(act_av))
length(hours_since_inj)
hours_since_inj<-as.character(hours_since_inj)
hours_since_inj<-as.numeric(hours_since_inj)

#make act_av into a vector
act_av<-as.vector(act_av)
length(act_av)

#make a new dataframe
act_av_by_trt<-data.frame(hours_since_inj, act_av, act_se, trt)
str(act_av_by_trt)

head(act_av_by_trt)
#pull out treatment levels
act_av_by_trt$soc<-substr(act_av_by_trt$trt, 1,1)
act_av_by_trt$lps<-substr(act_av_by_trt$trt, 3,3)

# #create values for edges of upper and lower se ribbons
# act_av_by_trt$upper_se_act<-act_av_by_trt$act_av + act_av_by_trt$act_se
# act_av_by_trt$lower_se_act<-act_av_by_trt$act_av - act_av_by_trt$act_se

#create values for edges of upper and lower 90% confidence intervals
act_av_by_trt$upper_90<-act_av_by_trt$act_av + 1.645*act_av_by_trt$act_se
act_av_by_trt$lower_90<-act_av_by_trt$act_av - 1.645*act_av_by_trt$act_se

names(act_av_by_trt)
range(act_av_by_trt$hours_since_inj)
pdata<-act_av_by_trt[act_av_by_trt$hours_since_inj < 24,]
pdata$min_since_inj<-pdata$hours_since_inj*60
names(pdata)
levels(pdata$trt)

pdata$trt<-factor(pdata$trt, labels = c("Mixed Control (no LPS)","All LPS","Mixed LPS-treated"))


pdata$trt<-ordered(pdata$trt, levels =c("Mixed Control (no LPS)", "Mixed LPS-treated","All LPS"))

head(pdata)
# actplot<-pdata %>%
#   select(min_since_inj, act_av, upper_se_act, lower_se_act, soc, lps,trt) %>%
#   filter(min_since_inj< 520)%>%
#   na.omit() %>%
#   ggplot() +
#   #facet_wrap(~lps, scales = "free")+
#   #geom_rect(aes(xmin=8.33, xmax=16.75, ymin=-Inf, ymax=Inf), fill="#A9A9A9", alpha=0.01, inherit.aes = F)+
#   geom_line(aes(x = min_since_inj, y = act_av, group= trt, linetype = trt)) +
#   geom_ribbon(aes(x = min_since_inj, ymin = lower_se_act, ymax= upper_se_act, fill = trt), alpha=0.35)+
#   ylim(0,1)+
#   scale_linetype_manual(values = c(1,2,3))+
#   scale_fill_manual(values = c("#6f8d42","#9559b7","#c05347")) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 14),
#         legend.title = element_blank(),
#         legend.direction = "vertical",
#         legend.position = c(0.18,0.9),
#         legend.background = element_blank(),
#         axis.line = element_line(),
#         axis.title = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())+
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=0.8)+
#   labs(x="Time since injection (minutes)", y = "Proportion of time spent active")
# 
# actplot

CI.actplot<-pdata %>%
  select(hours_since_inj, act_av, upper_90, lower_90, soc, lps,trt) %>%
  filter(hours_since_inj< 520)%>%
  na.omit() %>%
  ggplot() +
  geom_rect(data=band_dat,aes(xmin=7.4,xmax=8,ymin=-Inf,ymax=Inf),
            fill="slategray3", alpha = 0.01)+
  geom_ribbon(aes(x = hours_since_inj, ymin = lower_90, ymax= upper_90, fill = trt), alpha=0.35)+
  geom_line(aes(x = hours_since_inj, y = act_av, group= trt, linetype = trt),size =1.1) +
  ylim(0,1)+
  scale_linetype_manual(values = c(1,2,3))+
  scale_fill_manual(values = c("#6f8d42","#9559b7","#c05347")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = c(0.2,0.88),
        legend.key.size = unit(0.8,"cm"),
        legend.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.background = element_blank())+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=0.8)+
  labs(x="Time since injection (hours)", y = "Proportion of time spent active")
  

CI.actplot


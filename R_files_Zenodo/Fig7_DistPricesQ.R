########################################################################
#Figure 7: Distributions of prices and quantities.
########################################################################
tab  <- recent[which(recent$ev=="half"&recent$cost=="future"&recent$load==2045&recent$scen==3),]
tab <- transform(tab, category=ifelse(tab$category=="fossil","Fossil",paste(tab$category)))
tab <- transform(tab, category=ifelse(tab$category=="rps_100","100% Clean",paste(tab$category)))
tab <- transform(tab, category=ifelse(tab$category=="free","Unconstrained",paste(tab$category)))
tab$category<- factor(tab$category, levels = c("Fossil", "Unconstrained","100% Clean"))
tab$pricing<- factor(tab$pricing, levels = c("flat", "dynamic"))

mu <- ddply(tab, c("category","pricing"), summarise, 
            grp.mean=mean(net.final.q))
p1 <- ggplot(tab, aes(x=net.final.q,weight=ts_scale_to_period.x)) +
  geom_histogram(aes(y=..density..,color=factor(pricing),alpha=factor(pricing),fill=factor(pricing)),position="identity",bins=20)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=pricing),
             linetype="dashed",size=0.8,show.legend=F)+
  facet_grid(.~category)+
  xlab("Hourly quantity (MWh)") +
  ylab("Relative Frequency") + theme_minimal() +
  scale_fill_manual("Pricing", labels=c("Flat pricing","RTP - Pessimistic"),
                      values = c("gray80",brewer.pal(6, "Paired")[4:4]))+
  scale_color_manual("", labels=c("Flat pricing","RTP - Pessimistic"),
                    values = c("gray80",brewer.pal(6, "Paired")[4:4]), guide='none')+
  scale_alpha_discrete(range = c(0.75, 0.3), guide='none')+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  theme(legend.position = c(0.62, 0.80),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        legend.text = element_text(colour = "black", size = 15), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size=14),axis.title = element_text(size=14))
p1 
ggsave(
  paste("plot/Fig7a_HISTOQFuture2045EVhalfPessimistic.pdf", sep = ""),height=3.1,width=10
)


################################################################################
pricebins <- within(tab,quartile <- as.integer(cut(net.final.q,
            quantile(net.final.q,probs=0:5/5),include.lowers=T)))
pricebins <- na.omit(pricebins)
ggplot(pricebins, aes(y=net.final.price, x=factor(quartile),color=pricing)) +
  geom_boxplot()+
  facet_grid(.~category)+
  xlab("Hourly Quantity Quartile") +
  ylab("Hourly Price ($/MWh)") + theme_minimal()+
  scale_color_manual("", labels=c("Flat pricing","RTP - Pessimistic"),
                     values = c("gray80",brewer.pal(6, "Paired")[4:4]))+
  theme(legend.direction="vertical",
        legend.title = element_text(colour = "black", size = 15),
        legend.text = element_text(colour = "black", size = 15),
        legend.position = c(0.8, 0.90), axis.ticks=element_blank(), 
        plot.title = element_text(hjust = -10, size = 15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size=14),axis.title = element_text(size=14))
ggsave(
  paste("plot/Fig7b_BoxplotPriceFuture2045EVhalfPessimistic.pdf", sep = ""),height=3.1,width=10)
################################END#############################################




library(MASS)
library(ggplot2)
library(ggmosaic)
library(ggthemes)
library(dplyr)
library(reshape2)
countsToCases <- function(x, countcol = "Freq") {
    # Get the row indices to pull from x
    idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
    # Drop count column#
    x[[countcol]] <- NULL
    # Get the rows from x
    x[idx, ]
}
dcounts <- data.frame(
    farc=c("absent", "present", "absent", "present","absent", "present","absent", "present", "absent", "present","absent", "present"), 
    pa=c("inside","inside","inside","inside", "5 km buffer", "5 km buffer","5 km buffer", "5 km buffer", "outside", "outside", "outside", "outside"),
    year=c("2017", "2017", "2018", "2018","2017", "2017", "2018", "2018","2017", "2017", "2018", "2018"), 
    Freq=c(263, 136, 241, 870, 218, 157,378,389,8096,877,15070,1602) )
dcases<-countsToCases(dcounts)
dcases$factorC <-with(dcases, interaction(farc,  pa))
dc<-table(dcases$year, dcases$factorC)
pdf("in_and_out_bar.pdf", w=10, h=5)
barplot(as.matrix(dc), beside = TRUE, space = c(0.1, 0.3),
	las = 1, xlab = "Guerrilla and Protected Area", ylab = "Frequency",
	col = c( "goldenrod1","firebrick"), legend.text = rownames(dc), 
	args.legend = list(x = 0.3, y = max(dc), xjust = 0))
dev.off()
dc2<-dc[1:2,1:4]
pdf("in_bar.pdf", w=7, h=5)
barplot(as.matrix(dc2), beside = TRUE, space = c(0.1, 0.3),
	las = 1, xlab = "Guerrilla and Protected Area", ylab = "Frequency",
	col = c( "goldenrod1","firebrick"), legend.text = rownames(dc2), 
	args.legend = list(x = 0.3, y = max(dc2), xjust = 0))
dev.off()
dcounts$factorC <-with(dcounts, interaction(farc,  pa))
ggplot(data = dcounts) + geom_mosaic(aes(weight = Freq, x = product(factorC), fill=factor(year)), na.rm=TRUE) +labs(x="Guerrilla and Protected Area", y="proportion") + guides(fill=guide_legend(title = "Year"))+theme_pander()+scale_fill_manual(values=c( "goldenrod1","firebrick"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("in_and_out.pdf", h=5, w=7)
dc1<-subset(dcounts,pa!="outside")
dc1$factorC<-droplevels(dc1)$factorC
ggplot(data = dc1) + geom_mosaic(aes(weight = Freq, x = product(factorC), fill=factor(year)), na.rm=TRUE) +labs(x="Guerrilla and Protected Area", y="proportion") + guides(fill=guide_legend(title = "Year"))+theme_pander()+scale_fill_manual(values=c( "goldenrod1","firebrick"))+ theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("in.pdf", h=5, w=6)
dc3<-as.data.frame(dcounts %>% group_by(factorC) %>% mutate(fraction=Freq/sum(Freq)))
dc3$year<- factor(dc3$year, levels = c("2018", "2017"))
ggplot() + geom_bar(aes(y = fraction, x = factorC, fill = year), data = dc3, stat="identity")+labs(x="Guerrilla and Protected Area", y="proportion") + guides(fill=guide_legend(title = "Year"))+theme_pander()+scale_fill_manual(values=c( "firebrick","goldenrod1"))+ theme(axis.text.x = element_text(angle = 35, hjust = 1))                           
ggsave("stacked1.pdf", h=5, w=6)
dc3$pa<- factor(dc3$pa, levels = c("outside", "5 km buffer", "inside"))
ggplot() + geom_bar(aes(y = fraction, x = farc, fill = year), data = dc3, stat="identity")+labs(x="Guerrillas", y="proportion", subtitle = "Location") + guides(fill=guide_legend(title = "Year"))+theme_pander()+scale_fill_manual(values=c( "firebrick","goldenrod1"))+ theme(axis.text.x = element_text(angle = 35, hjust = 1))+ facet_grid(~ pa, margins=F)
ggsave("stacked2.pdf", h=5, w=6)
time<-read.csv("time.csv")
t<-ts(data=time[,2:4], start=c(2001), end=c(2018), frequency=1)
t.melt<-melt(time,id.vars='Year', variable.name = 'location', value.name = 'fire_density')
t.melt$location<- factor(t.melt$location, levels = c("inside", "buffer", "outside"))
ggplot(data = t.melt, aes(x = Year, y = fire_density))+geom_line(aes(colour = location))+labs(x="Year", y="fire density (fires per km square)")+scale_color_few()+theme_pander()
ggsave("time_series.pdf", h=5, w=6)
mod0 <- glm(Freq ~ pa + year + farc, data=dcounts, family=poisson())
mod1 <- glm(Freq ~ pa * year + farc, data= dcounts, family=poisson())
mod2 <- glm(Freq ~ (pa + year) * farc, data=dcounts, family=poisson())
sink("contingency_in_and_out.txt")
#Are pa associated with year?
print("Are pa associated with year?")
print(anova(mod0, mod1, test="Chi")) #Test of joint independence
#Are pa independent of year, given the effect of farc?
print("Are pa independent of year, given the effect of farc?")
print(anova(mod0, mod2, test="Chi")) #Test of conditional independence
sink()
mod3 <- glm(Freq ~ pa*year*farc - pa:year:farc, data=dcounts, family=poisson())
mod4 <- glm(Freq ~ pa*year*farc, data=dcounts, family=poisson())
#Does the relationship between year and pa depend on the farc?
sink("contingency_in_and_out.txt", append=T)
print("Does the relationship between year and pa depend on the farc?")
print(anova(mod3, mod4, test="Chi"))
print(pchisq(deviance(mod3), df.residual(mod3), lower=F))
sink()
mod0i <- glm(Freq ~ pa + year + farc, data=dc1, family=poisson())
mod1i <- glm(Freq ~ pa * year + farc, data= dc1, family=poisson())
mod2i <- glm(Freq ~ (pa + year) * farc, data=dc1, family=poisson())
sink("contingency_in.txt")
#Are pa associated with year?
print("Are pa associated with year?")
print(anova(mod0i, mod1i, test="Chi")) #Test of joint independence
#Are farc independent of pa, given the effect of year?
print("Are pa independent of year, given the effect of farc?")
print(anova(mod0i, mod2i, test="Chi")) #Test of conditional independence
sink()
mod3i <- glm(Freq ~ pa*year*farc - pa:year:farc, data=dc1, family=poisson())
mod4i <- glm(Freq ~ pa*year*farc, data=dc1, family=poisson())
#Does the relationship between year and pa depend on the farc?
sink("contingency_in.txt", append=T)
print("Does the relationship between year and pa depend on the farc?")
print(anova(mod3i, mod4i, test="Chi"))
print(pchisq(deviance(mod3), df.residual(mod3i), lower=F))
sink()
sink("model_parameters.txt")
#Are pa associated with year?
print("Are pa associated with year?")
print(summary(mod1)) #Test of joint independence
#Are farc independent of pa, given the effect of year?
print("Are pa independent of year, given the effect of farc?")
print(summary(mod2)) #Test of conditional independence
#Does the relationship between year and pa depend on the farc?
print("Does the relationship between year and pa depend on the farc?")
print(summary(mod3)) 
sink()
sink("model_parameters_in.txt")
#Are pa associated with year?
print("Are pa associated with year?")
print(summary(mod1i)) #Test of joint independence
#Are farc independent of pa, given the effect of year?
print("Are pa independent of year, given the effect of farc?")
print(summary(mod2i)) #Test of conditional independence
#Does the relationship between year and pa depend on the farc?
print("Does the relationship between year and pa depend on the farc?")
print(summary(mod3i)) 
sink()
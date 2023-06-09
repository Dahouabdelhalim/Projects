library(ggplot2)
library(doBy) #for summaryBy() function to get mean/sd/se/n


d <- Lesions_2020
d <- as.data.frame(d)

fun <- function(x){c(m=mean(x), sd=sd(x),se=sd(x)/sqrt(length(x)), n=length(x))}

sum <- summaryBy(Lesion_Cover ~ Day_Year, data=d, FUN=fun)
sum_readable <- round(sum, digits = 1)

sum_readable$year <- c("2020", "2019", "2020", "2019", "2020", "2019", "2020", "2019", "2020", "2020")

ggplot(sum_readable, aes(x=Day_Year, y=Lesion_Cover.m, color=year))+
  geom_line(size=1.2)+
  geom_errorbar(aes(ymax=Lesion_Cover.m+Lesion_Cover.se,
                    ymin=Lesion_Cover.m-Lesion_Cover.se),
                width=.05)+
  scale_color_manual(values = c("darkolivegreen4", "darkolivegreen3"))+
  theme_classic()+
  ylab("Average % Lesion Cover")+ xlab("")+
  geom_point(size=3)+
  scale_x_continuous(breaks=c(175, 200, 225, 250), labels=c("June","July","August","Sept"))+
  theme(text=element_text(size=12))
ggsave("Field_Lesions_nomean_presentation.png", dpi=300, units="in", width=6, height=3)


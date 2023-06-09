library(dplyr)

### Read raw data from .csv file ###
eggdata <- read.csv("https://zenodo.org/record/4014998/files/1_Data_Eggshell_Thickness_1972-2019.csv")
fragdata <- subset(eggdata, Fragment=="+")
wholedata <- subset(eggdata, Fragment=="-")

### Correction for missing membrane (e.g. eggdata$Membrane == "-") ###
corr <- 0.071
eggdata$shell_c <- ifelse(eggdata$Membrane == "-", (eggdata$Shell + corr), eggdata$Shell)
fragdata$shell_c <- ifelse(fragdata$Membrane == "-", (fragdata$Shell + corr), fragdata$Shell)
wholedata$shell_c <- ifelse(wholedata$Membrane == "-", (wholedata$Shell + corr), wholedata$Shell)

### Clutch means ###
grouped_yl <- eggdata %>% group_by(Year, Nest_Site, Fragment)
grouped_yl_plot <- grouped_yl %>% summarise(tyk = mean(shell_c))
grouped_yl_year <- eggdata %>% group_by(Year, Nest_Site)
grouped_yl_year_plot <- grouped_yl_year %>% summarise(tyk = mean(shell_c))

### Check that count of fragments >= 20 for each clutch ###
grouped_yl_length_frag <- subset(grouped_yl, Fragment=="+") %>% summarise(tyk = length(shell_c))
grouped_yl_length_whole <- subset(grouped_yl, Fragment=="-") %>% summarise(tyk = length(shell_c))

### Year means, including linear regression line ###
grouped_year <- grouped_yl_year_plot %>% group_by(Year) # fragments and whole eggs combined for each clutch
grouped_year_plot <- grouped_year %>% summarise(tyk = mean(tyk))
linreg <- lm(grouped_year_plot$tyk ~ grouped_year_plot$Year)

### 95% confidence intervals for year means ###
conf_interval<-predict(linreg, newdata=data.frame(x=grouped_year_plot$Year), interval="confidence", level = 0.95)

### 95% prediction intervals for year means ###
pred_interval<-predict(linreg, newdata=data.frame(x=grouped_year_plot$Year), interval="prediction", level = 0.95)

### Calculate average increase per year (estimate of grouped_year_plot$Year) ### 
linreg_log <- lm(log(grouped_year_plot$tyk) ~ grouped_year_plot$Year)
summary(linreg_log)

### Plot data ###
par(fig=c(0, 0.80, 0, 1), cex.axis=1)
x_lim <- c(1971.5, 2030)
y_lim <- c(0.21, 0.43)
plot(fragdata$Year, fragdata$shell_c, xlab="Year", ylab="Thickness (mm) including membrane", main="Shell  thickness  of  Peregrine  Falcon  eggs  from  Greenland", type="p", col="green", xlim=x_lim, ylim=y_lim, pch=16, cex=0.6)
par(new=TRUE)
plot(wholedata$Year, wholedata$shell_c, xlab=NA, ylab=NA, type="p", col="orange", xlim=x_lim, ylim=y_lim, pch=16, cex=0.6, axes=F)
matlines(grouped_year_plot$Year, conf_interval[,2:3], col = "black", lty=2, lwd=1)
matlines(grouped_year_plot$Year, pred_interval[,2:3], col="black", lty=3, lwd=1)
polygon(c(rev(grouped_year_plot$Year), grouped_year_plot$Year), c(rev(pred_interval[ ,3]), pred_interval[ ,2]), col = rgb(200, 200, 200, max = 255, alpha = 50) , border = NA)
polygon(c(rev(grouped_year_plot$Year), grouped_year_plot$Year), c(rev(conf_interval[ ,3]), conf_interval[ ,2]), col = rgb(200, 200, 200, max = 255, alpha = 100) , border = NA)
par(new=TRUE)
plot(grouped_yl_plot$Year, grouped_yl_plot$tyk, xlab=NA, ylab=NA, type="p", col="blue", pch=16, xlim=x_lim, ylim=y_lim, cex=1, axes = F)
par(new=TRUE)
plot(grouped_year_plot$Year, grouped_year_plot$tyk, xlab=NA, ylab=NA, type="p", col="red", pch=16, xlim=x_lim, ylim=y_lim, cex=1.5, axes = F)
#plot(grouped_year_plot$Year, grouped_year_plot$tyk, xlab=NA, ylab=NA, type="p", col="blue", pch=16, xlim=x_lim, ylim=y_lim, cex=.5, axes = F)
abline(linreg, lty=1, col="red", lwd=2)

### Draw critical line (red) and pre-DDT line (blue) ###
abline(h=0.279, lty=2, col="red")
abline(h=0.336, lty=2, col="blue")

### Show legends ###
legend("topleft", bg=rgb(1,1,1,0.60), bty="n", cex=0.75, c("Individual fragment measurements", "Individual whole egg measuremets", "Clutch means", "Year means (based on clutch means)", paste("Linear regression (year means), slope =", format(coef(linreg)["grouped_year_plot$Year"], digits=3)), "95% confidence limits", "95% prediction limits", "Pre-DDT thickness (0.336mm)", "Critical thickness reduction (17%)"), pt.cex=c(0.6, 0.6, 1, 1.5), col=c("green", "orange", "blue", "red", "red", "black", "black", "blue", "red"), pch=c(16, 16, 16, 16, NA, NA, NA, NA, NA), lty=c(NA, NA, NA, NA, 1, 2, 3, 2, 2), lwd=c(NA, NA, NA, NA, 2, 1, 1, 1, 1))

### Boxplots of clutch means for each decade ###
par(fig=c(0.73, 1, 0, 1), new=TRUE, cex.axis=0.75)
int70 <- (subset(grouped_yl_plot, Year>=1970 & Year<1980))
int80 <- (subset(grouped_yl_plot, Year>=1980 & Year<1990))
int90 <- (subset(grouped_yl_plot, Year>=1990 & Year<2000))
int00 <- (subset(grouped_yl_plot, Year>=2000 & Year<2010))
int10 <- (subset(grouped_yl_plot, Year>=2010 & Year<2020))
boxplot(int70$tyk, int80$tyk, int90$tyk, int00$tyk, int10$tyk, ylim=y_lim, ylab=NA, names=NA, yaxt='n', xlab="Decades", cex=0.75)
axis(1, at=1:5, labels=c("70-79", "80-89", "90-99", "00-09", "10-19"), las=2)
legend("topleft", bg=rgb(1,1,1,0.60), bty="n", cex=0.75, c("Boxplots of", "clutch means", "for each decade"))


setwd("/Users/dbolnick/Documents/Am Nat EiC/AmNat history/DataAnalysis")
dat <- read.csv("AmNatTheoryData.csv")
dat[is.na(dat)] <- 0
head(dat)
d <- dat[,-c(1,2)]
d.rel <- d/rowSums(d)
npapers_per_issue <- rowSums(d)/3



{
##### Figure 1A.
par(mar = c(8,5,2,1), mfrow = c(1,1))
summary <- summary[order(summary, decreasing = T)]
barplot(summary, horiz = F, beside = F, legend.text = F, col = "grey70", las = 2, ylab = "Proportion of papers, 1980-2017", ylim = c(0,0.3))
box()
}

{##### Figure 1B.
#2017; excuding synthses:
par(mar = c(0.2,2,0.2,6), oma = c(0,0,0,0), par(mfrow = c(1,1)))
library(VennDiagram)
draw.pairwise.venn(area1 = 62, area2 = 103, cross.area = 22, category = c( "Theoretical", "Empirical"), cat.cex = 1.4, fill = rgb(0:1,0.3,1:0,0.5), cat.pos = c(0,0), cex = 1.4)
}




# Figure 3
par(mfrow = c(1,2), mar = c(5,5,3, 1))
{
  plot(d.rel$Paleontology ~ dat$Year, pch = 16, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Paleontology & Geology", cex.lab = 1.5)
  spl <-  smooth.spline(dat$Year, d.rel$Paleontology, df = 5)
  lines(spl, lwd = 3)
  text(1875, 0.29, "A", cex = 2)
  
  points( dat$Year,d.rel$Geology , pch = 1, col = "black", cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Geology", cex.lab = 1.5)
  spl <-  smooth.spline(dat$Year, d.rel$Geology, df = 5)
  lines(spl, lwd = 2, col = "black", lty = 2)

}


{
  plot(d.rel$Anthropology_archaeology ~ dat$Year, pch = 16, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Anthropology", cex.lab = 1.5, ylim = c(0,0.25))
  spl <-  smooth.spline(dat$Year, d.rel$Anthropology_archaeology, df = 5)
  lines(spl, lwd = 3, lty = 1)
  text(1875, 0.29, "B", cex = 2)
  
}



{
# Figure 4
par( mar = c(5,5,3,01), mfrow = c(1,2))

{
  plot(d.rel$Zoology ~ dat$Year, pch = 16, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Descriptive Zoology & Botany", cex.lab = 1.5)
  spl <-  smooth.spline(dat$Year, d.rel$Zoology, df = 5)
  lines(spl, lwd = 3)
  
  points(dat$Year, d.rel$Botany ,  pch = 1, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Botany", col = "black")
  spl <-  smooth.spline(dat$Year, d.rel$Botany, df = 5)
  lines(spl, col = "black", lwd = 2, lty = 2)
}


{
  plot(d.rel$Genetics ~ dat$Year, pch = 16, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Genetics & Ecology", , cex.lab = 1.5)
spl <-  smooth.spline(dat$Year, d.rel$Genetics, df = 5)
lines(spl, lwd = 3)
#lines(dat$Year, d.rel$Genetics)

points( dat$Year, d.rel$Ecology , pch = 1, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Ecology", col = "black", cex.lab = 1.5)
spl <-  smooth.spline(dat$Year, d.rel$Ecology, df = 5)
lines(spl, col = "black", lwd = 3)
#lines(dat$Year, d.rel$Ecology)
}
  }




{
#Figure 5
par( mar = c(5,5,3,01), mfrow = c(1,1))

plot(d.rel$Evolution ~ dat$Year, pch = 16, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Evolution", cex.lab = 1.5)
spl <-  smooth.spline(dat$Year, d.rel$Evolution, df = 6)
lines(spl, lwd = 3)

}






sex <-read.csv("SexRatioFirstIssueByYear.csv")
sex$ratio <- sex$N_women / (sex$N_women + sex$N_men)
plot(ratio ~ Year, sex, ylab = "Proprtion of women first authors", pch= 16, cex.lab = 1.5)
spl <-  smooth.spline(sex$Year, sex$ratio)
lines(spl, lwd = 3)
abline(h = 0.5)




plot(d.rel$Ecology ~ d.rel$Genetics)



{
  par(mar = c(5,5,2,1))
  include <- dat$Year >1950 & dat$Year <1990
  plot(d.rel$Genetics[include] ~ dat$Year[include], pch = 16, cex = 1, xlab = "Year", ylab = "Proportion of papers" , main = "Genetics (black) & Ecology (red)", cex.lab = 1.5)
  spl <-  smooth.spline(dat$Year[include], d.rel$Genetics[include], df = 5)
  lines(spl, lwd = 3)
  
  points(dat$Year[include], d.rel$Ecology[include] ,  pch = 16, cex = 1, col = "red")
  spl <-  smooth.spline(dat$Year[include], d.rel$Ecology[include], df = 5)
  lines(spl, col = "red", lwd = 3)
}



#### 
# Theory figure
data_dir <- file.path("AmNatTheoryData.csv")

# Import data
dat <- read.csv(file = data_dir, stringsAsFactors = F)

# Modify data
# One issue was a special one an wasn't numeric, but change here for analysis; gives warning
dat$Number[is.na(as.numeric(dat$Number))] <- 2
# Since entered in order, make all numbers 1, 2, 3 use as a decimal instead of counting years 3 times
dat$NumberMod <- rep(x = 1:3, times = length(unique(dat$Year)))
dat$YearNo <- dat$Year + (dat$NumberMod/10)

# Create new data
# dataframe with mean per year
year_mean <- tapply(X = dat$PropTheory, INDEX = dat$Year, FUN = mean, na.rm = T)
years <- as.numeric(names(year_mean))
df <- data.frame(Year = years, TotalTheory_mean = year_mean, row.names = NULL)

par(mar = c(5,5,1,1))
plot(x = dat$Year, y = dat$PropTheory, type = "p", pch = 16, las = 1, xlab = "Year", ylab = "Proportion of papers with quantiative theory", cex.lab = 1.4 , xlim = c(1920, 2020))
spl <-  smooth.spline(dat$Year, y = dat$PropTheory, df = 6)
lines(spl, lwd = 3)


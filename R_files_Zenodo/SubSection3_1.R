#This source code reproduces the results shown in subsection
#3.1 "Context for comparison"

#author: Andrea Falcon-Cortes

#Necessary packages
library(dplyr)
library(tidyverse)

#Data relative to Federal Budget from 2013 to 2020 in USD PPP

pf13 = 3956361600000/7.884
pf14 = 4467225800000/8.045
pf15 = 4694677400000/8.328
pf16 = 4763874000000/8.446
pf17 = 4888892500000/8.914
pf18 = 5279667000000/9.174
pf19 = 5838059700000/9.285
pf20 = 6107732400000/9.380

pf = c(pf13, pf14, pf15, pf16, pf17, pf18, pf19, pf20)

#Data relative to contracts made from 2013 to 2020

#Upload the source list
data <- read.csv("~/RPS_Cont_2013-2020.csv", 
                               stringsAsFactors=FALSE)

#Separate by years using filter function
dataEPN_13=data %>% filter(Year==2013)
dataEPN_14=data %>% filter(Year==2014)
dataEPN_15=data %>% filter(Year==2015)
dataEPN_16=data %>% filter(Year==2016)
dataEPN_17=data %>% filter(Year==2017)
dataEPN_18=data %>% filter(Year==2018)
dataAMLO_19=data %>% filter(Year==2019)
dataAMLO_20=data %>% filter(Year==2020)

#Separate each year by classes
#NC
dataEPN_13.nc = dataEPN_13 %>% filter(Status == "NC")
dataEPN_14.nc = dataEPN_14 %>% filter(Status == "NC")
dataEPN_15.nc = dataEPN_15 %>% filter(Status == "NC")
dataEPN_16.nc = dataEPN_16 %>% filter(Status == "NC")
dataEPN_17.nc = dataEPN_17 %>% filter(Status == "NC")
dataEPN_18.nc = dataEPN_18 %>% filter(Status == "NC")
dataAMLO_19.nc = dataAMLO_19 %>% filter(Status == "NC")
dataAMLO_20.nc = dataAMLO_20 %>% filter(Status == "NC")

#PCS
dataEPN_13.pcs = dataEPN_13 %>% filter(Status == "PCS")
dataEPN_14.pcs = dataEPN_14 %>% filter(Status == "PCS")
dataEPN_15.pcs = dataEPN_15 %>% filter(Status == "PCS")
dataEPN_16.pcs = dataEPN_16 %>% filter(Status == "PCS")
dataEPN_17.pcs = dataEPN_17 %>% filter(Status == "PCS")
dataEPN_18.pcs = dataEPN_18 %>% filter(Status == "PCS")
dataAMLO_19.pcs = dataAMLO_19 %>% filter(Status == "PCS")
dataAMLO_20.pcs = dataAMLO_20 %>% filter(Status == "PCS")

#EFOS
dataEPN_13.efos = dataEPN_13 %>% filter(Status == "EFOS")
dataEPN_14.efos = dataEPN_14 %>% filter(Status == "EFOS")
dataEPN_15.efos = dataEPN_15 %>% filter(Status == "EFOS")
dataEPN_16.efos = dataEPN_16 %>% filter(Status == "EFOS")
dataEPN_17.efos = dataEPN_17 %>% filter(Status == "EFOS")
dataEPN_18.efos = dataEPN_18 %>% filter(Status == "EFOS")
dataAMLO_19.efos = dataAMLO_19 %>% filter(Status == "EFOS")
dataAMLO_20.efos = dataAMLO_20 %>% filter(Status == "EFOS")

#Data relative to spending in public procurement per year 
#in USD PPP
gasto13 = sum(dataEPN_13$Spending)
gasto14 = sum(dataEPN_14$Spending)
gasto15 = sum(dataEPN_15$Spending)
gasto16 = sum(dataEPN_16$Spending)
gasto17 = sum(dataEPN_17$Spending)
gasto18 = sum(dataEPN_18$Spending)
gasto19 = sum(dataAMLO_19$Spending)
gasto20 = sum(dataAMLO_20$Spending)

gasto = c(gasto13, gasto14, gasto15, gasto16, gasto17, gasto18, gasto19, gasto20)

#Spending in public procurement relative to the federal
#budget in USD PPP 
pf.gs = gasto/pf

#Contracts per year
ncont13 = nrow(dataEPN_13)
ncont14 = nrow(dataEPN_14)
ncont15 = nrow(dataEPN_15)
ncont16 = nrow(dataEPN_16)
ncont17 = nrow(dataEPN_17)
ncont18 = nrow(dataEPN_18)
ncont19 = nrow(dataAMLO_19)
ncont20 = nrow(dataAMLO_20)

ncont = c(ncont13, ncont14, ncont15, ncont16, 
          ncont17, ncont18, ncont19, ncont20)

#Some parameters for graphs
anio = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)
col.n = c("darkgreen","darkgreen","darkgreen","darkgreen",
          "darkgreen","darkgreen",
          adjustcolor("purple", alpha.f = 0.5),
          adjustcolor("purple", alpha.f = 0.5))
min_pfgs = 0.04
max_pfgs = 0.17
norm_pfgs = (pf.gs - min_pfgs)/(max_pfgs-min_pfgs)
min_ncont = 137e3
max_ncont = 227e3
norm_ncont = (ncont - min_ncont)/(max_ncont - min_ncont)

#Figure 1 - Up
png("Fig1-Up.png", units = "in", width = 5, height = 5, res = 300)
par(oma=c(0,0,0,3))
plot(anio, norm_pfgs, las=1, axes=F,ylim=c(-0.01,1.01), 
     xlab="Year", ylab="", type="l",lty=2, lwd = 1, col = "cornflowerblue")
points(anio, norm_pfgs, col = col.n, pch=16, cex = 1.25)
lines(anio,norm_ncont, col="orange", lty = 2, lwd = 1)
points(anio,norm_ncont,col=col.n,pch=16, cex = 1.25)
axis(1)
axis(2, at = seq(0,1,0.2),las=1, 
     labels = seq( min_pfgs , max_pfgs,  length.out=6),
     col = "cornflowerblue", col.axis = "cornflowerblue")
mtext("TS/FB",side=2,adj=0.5,outer = T,padj = 1.5,col="cornflowerblue")
axis(4, at = seq(0,1,0.2),las=1, labels= seq(min_ncont/1e03, max_ncont/1e03,length.out=6),
     col="orange",col.axis="orange")
mtext("Total Contracts",side=4,adj=0.55,outer = T,padj = 1,col="orange")
legend("topright",col = c("darkgreen",adjustcolor("purple", alpha.f = 0.5)), legend = c("1st Period", "2nd Period"), 
       pch = c(16,16), bty = "n", cex = 1.0, pt.cex = c(1.25,1.25))
dev.off()

#Total spending and total of contracts per class per year
#Spending is reporting in USD PPP
#EFOS Contracts
efos13 = nrow(dataEPN_13.efos)
efos14 = nrow(dataEPN_14.efos)
efos15 = nrow(dataEPN_15.efos)
efos16 = nrow(dataEPN_16.efos)
efos17 = nrow(dataEPN_17.efos)
efos18 = nrow(dataEPN_18.efos)
efos19 = nrow(dataAMLO_19.efos)
efos20 = nrow(dataAMLO_20.efos)
#EFOS Spending
efos13.mont = sum(dataEPN_13.efos$Spending)
efos14.mont = sum(dataEPN_14.efos$Spending)
efos15.mont = sum(dataEPN_15.efos$Spending)
efos16.mont = sum(dataEPN_16.efos$Spending)
efos17.mont = sum(dataEPN_17.efos$Spending)
efos18.mont = sum(dataEPN_18.efos$Spending)
efos19.mont = sum(dataAMLO_19.efos$Spending)
efos20.mont = sum(dataAMLO_20.efos$Spending)

#PCS Contracts
pcs13 = nrow(dataEPN_13.pcs)
pcs14 = nrow(dataEPN_14.pcs)
pcs15 = nrow(dataEPN_15.pcs)
pcs16 = nrow(dataEPN_16.pcs)
pcs17 = nrow(dataEPN_17.pcs)
pcs18 = nrow(dataEPN_18.pcs)
pcs19 = nrow(dataAMLO_19.pcs)
pcs20 = nrow(dataAMLO_20.pcs)
#PCS Spending
pcs13.mont = sum(dataEPN_13.pcs$Spending)
pcs14.mont = sum(dataEPN_14.pcs$Spending)
pcs15.mont = sum(dataEPN_15.pcs$Spending)
pcs16.mont = sum(dataEPN_16.pcs$Spending)
pcs17.mont = sum(dataEPN_17.pcs$Spending)
pcs18.mont = sum(dataEPN_18.pcs$Spending)
pcs19.mont = sum(dataAMLO_19.pcs$Spending)
pcs20.mont = sum(dataAMLO_20.pcs$Spending)

#NC Contracts
nc13 = nrow(dataEPN_13.nc)
nc14 = nrow(dataEPN_14.nc)
nc15 = nrow(dataEPN_15.nc)
nc16 = nrow(dataEPN_16.nc)
nc17 = nrow(dataEPN_17.nc)
nc18 = nrow(dataEPN_18.nc)
nc19 = nrow(dataAMLO_19.nc)
nc20 = nrow(dataAMLO_20.nc)
#NC Spending
nc13.mont = sum(dataEPN_13.nc$Spending)
nc14.mont = sum(dataEPN_14.nc$Spending)
nc15.mont = sum(dataEPN_15.nc$Spending)
nc16.mont = sum(dataEPN_16.nc$Spending)
nc17.mont = sum(dataEPN_17.nc$Spending)
nc18.mont = sum(dataEPN_18.nc$Spending)
nc19.mont = sum(dataAMLO_19.nc$Spending)
nc20.mont = sum(dataAMLO_20.nc$Spending)

#Some parameters for graphs
year = c("2013", "2014", "2015", "2016", 
         "2017", "2018", "2019", "2020")

#Figure 1 - Center_Left
class.cont = data.frame(row.names = c("EFOS", "PCS", "NC"),
                        "2013" = c(efos13, pcs13, nc13),
                        "2014" = c(efos14, pcs14, nc14),
                        "2015" = c(efos15, pcs15, nc15),
                        "2016" = c(efos16, pcs16, nc16),
                        "2017" = c(efos17, pcs17, nc17),
                        "2018" = c(efos18, pcs18, nc18),
                        "2019" = c(efos19, pcs19, nc19),
                        "2020" = c(efos20, pcs20, nc20))
names(class.cont) = c ("2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")
col.b = c(adjustcolor("darkgreen", alpha.f = 0.05), 
          adjustcolor("darkgreen", alpha.f = 0.25),
          adjustcolor("darkgreen", alpha.f = 0.45),
          adjustcolor("darkgreen", alpha.f = 0.65),
          adjustcolor("darkgreen", alpha.f = 0.85),
          adjustcolor("darkgreen", alpha.f = 1.00),
          adjustcolor("purple", alpha.f = 0.5),
          adjustcolor("purple", alpha.f = 1.00))
png("Fig1-Center_Left.png", units = "in", width = 5, height = 5, res = 300)
barplot(t(as.matrix(class.cont)), beside=T, log = "y", col = col.b, las = 1, cex.axis = 0.85,
        cex.names = 1.0, ylab = "Contracts", xlab = "Class", cex.lab = 1.0,
        ylim = c(1e+00,1e07))
legend("topleft",legend = year, fill = col.b, bty = "n", cex = 1.0, ncol = 2)
dev.off()

#Figure 1 - Center_Right
class.mont = data.frame(row.names = c("EFOS", "PCS", "NC"),
                        "2013" = c(efos13.mont, pcs13.mont, nc13.mont),
                        "2014" = c(efos14.mont, pcs14.mont, nc14.mont),
                        "2015" = c(efos15.mont, pcs15.mont, nc15.mont),
                        "2016" = c(efos16.mont, pcs16.mont, nc16.mont),
                        "2017" = c(efos17.mont, pcs17.mont, nc17.mont),
                        "2018" = c(efos18.mont, pcs18.mont, nc18.mont),
                        "2019" = c(efos19.mont, pcs19.mont, nc19.mont),
                        "2020" = c(efos20.mont, pcs20.mont, nc20.mont))
names(class.mont) = c ("2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")
col.b = c(adjustcolor("darkgreen", alpha.f = 0.05), 
          adjustcolor("darkgreen", alpha.f = 0.25),
          adjustcolor("darkgreen", alpha.f = 0.45),
          adjustcolor("darkgreen", alpha.f = 0.65),
          adjustcolor("darkgreen", alpha.f = 0.85),
          adjustcolor("darkgreen", alpha.f = 1.00),
          adjustcolor("purple", alpha.f = 0.5),
          adjustcolor("purple", alpha.f = 1.00))
png("Fig1-Center_Right.png", units = "in", width = 5, height = 5, res = 300)
barplot(t(as.matrix(class.mont)), beside=T, log = "y", col = col.b, las = 1, cex.axis = 0.85,
        cex.names = 1.0, ylab = "Spending", xlab = "Class", cex.lab = 1.0,
        ylim = c(1e+06,1e12))
legend("topleft",legend = year, fill = col.b, bty = "n", cex = 1.0, ncol = 2)
dev.off()

#Figure 1 - Bottom_Left
class.cont.p = data.frame(row.names = c("EFOS", "PCS", "NC"),
                          "2013" = c(efos13/nrow(dataEPN_13),
                                     pcs13/nrow(dataEPN_13), 
                                     nc13/nrow(dataEPN_13)),
                          "2014" = c(efos14/nrow(dataEPN_14),
                                     pcs14/nrow(dataEPN_14), 
                                     nc14/nrow(dataEPN_14)),
                          "2015" = c(efos15/nrow(dataEPN_15), 
                                     pcs15/nrow(dataEPN_15), 
                                     nc15/nrow(dataEPN_15)),
                          "2016" = c(efos16/nrow(dataEPN_16),
                                     pcs16/nrow(dataEPN_16), 
                                     nc16/nrow(dataEPN_16)),
                          "2017" = c(efos17/nrow(dataEPN_17), 
                                     pcs17/nrow(dataEPN_17), 
                                     nc17/nrow(dataEPN_17)),
                          "2018" = c(efos18/nrow(dataEPN_18), 
                                     pcs18/nrow(dataEPN_18), 
                                     nc18/nrow(dataEPN_18)),
                          "2019" = c(efos19/nrow(dataAMLO_19), 
                                     pcs19/nrow(dataAMLO_19), 
                                     nc19/nrow(dataAMLO_19)),
                          "2020" = c(efos20/nrow(dataAMLO_20), 
                                     pcs20/nrow(dataAMLO_20), 
                                     nc20/nrow(dataAMLO_20)))
names(class.cont.p) = c ("2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")
col.b = c(adjustcolor("darkgreen", alpha.f = 0.05), 
          adjustcolor("darkgreen", alpha.f = 0.25),
          adjustcolor("darkgreen", alpha.f = 0.45),
          adjustcolor("darkgreen", alpha.f = 0.65),
          adjustcolor("darkgreen", alpha.f = 0.85),
          adjustcolor("darkgreen", alpha.f = 1.00),
          adjustcolor("purple", alpha.f = 0.5),
          adjustcolor("purple", alpha.f = 1.00))
png("Fig1-Bottom_Left.png", units = "in", width = 5, height = 5, res = 300)
barplot(t(as.matrix(class.cont.p*100)), beside=T, log = "y", col = col.b, las = 1, cex.axis = 0.85,
        cex.names = 1.0, ylab = "% of TC", xlab = "Class", cex.lab = 1.0,
        ylim = c(0.001,1e3))
legend("topleft",legend = year, fill = col.b, bty = "n", cex = 1.0, ncol = 2)
dev.off()

#Figure 1 -  Bottom_Right
class.mont.p = data.frame(row.names = c("EFOS", "PCS", "NC"),
                          "2013" = c(efos13.mont/gasto13,
                                     pcs13.mont/gasto13, 
                                     nc13.mont/gasto13),
                          "2014" = c(efos14.mont/gasto14,
                                     pcs14.mont/gasto14, 
                                     nc14.mont/gasto14),
                          "2015" = c(efos15.mont/gasto15, 
                                     pcs15.mont/gasto15, 
                                     nc15.mont/gasto15),
                          "2016" = c(efos16.mont/gasto16,
                                     pcs16.mont/gasto16, 
                                     nc16.mont/gasto16),
                          "2017" = c(efos17.mont/gasto17, 
                                     pcs17.mont/gasto17, 
                                     nc17.mont/gasto17),
                          "2018" = c(efos18.mont/gasto18, 
                                     pcs18.mont/gasto18, 
                                     nc18.mont/gasto18),
                          "2019" = c(efos19.mont/gasto19, 
                                     pcs19.mont/gasto19, 
                                     nc19.mont/gasto19),
                          "2020" = c(efos20.mont/gasto20, 
                                     pcs20.mont/gasto20, 
                                     nc20.mont/gasto20))
names(class.mont.p) = c ("2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")
col.b = c(adjustcolor("darkgreen", alpha.f = 0.05), 
          adjustcolor("darkgreen", alpha.f = 0.25),
          adjustcolor("darkgreen", alpha.f = 0.45),
          adjustcolor("darkgreen", alpha.f = 0.65),
          adjustcolor("darkgreen", alpha.f = 0.85),
          adjustcolor("darkgreen", alpha.f = 1.00),
          adjustcolor("purple", alpha.f = 0.5),
          adjustcolor("purple", alpha.f = 1.00))
png("Fig1-Bottom_Right.png", units = "in", width = 5, height = 5, res = 300)
barplot(t(as.matrix(class.mont.p*100)), beside=T, log = "y", col = col.b, las = 1, cex.axis = 0.85,
        cex.names = 1.0, ylab = "% of TS", xlab = "Class", cex.lab = 1.0,
        ylim = c(0.001,1e3))
legend("topleft",legend = year, fill = col.b, bty = "n", cex = 1.0, ncol = 2)
dev.off()
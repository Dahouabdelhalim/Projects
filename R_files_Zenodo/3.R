#### Calculating confidence intervals of proportions #####

#### Mammals #### 

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_mam<-c(3,4,0,2,2,2,3,2,6,3,4,8,1,3,6,6,7,3,17,14,6, 10,20,10,7,9,3,10,17,14,10,8.6,5,3, 4, 5,2, 5, 8,3,12,12, 6,6, 11, 4) # use Tracy et al. 2001 data averaged with Ward
length(successes_mam) # 46
plot(successes_mam) # does increase

total_mam<-c(78.54,133.137931, 98.8, 96.8,155,173.3333333,96.4, 99.66666667,164.5540541,133.890411,145,198.9333333, 175, 194.5666667,139.7142857, 168.3414634,230,  215,234.0208333, 238.8965517,270,198.1818182, 267.9428571,319.4018692,292.8348624,289.4956522,253.3414634,248.4615385,379.4864865,389.48,443.8289474, 330.9114,210.68,243.1851852,318.1714286,160.58,221.72,217.36,166.84, 250.9090909,213.28,192.6, 194.34, 178.88,225.8181818, 177.16)
length(total_mam) # 46
plot(total_mam) # 

plot(successes_mam/total_mam)

mam_data<-data.frame(year_all,successes_mam,total_mam)

head(mam_data)

confint<-binom.confint(mam_data$successes_mam,mam_data$total_mam,method="exact")

proport<-successes_mam/total_mam
lower<-confint$lower
upper<-confint$upper

mam_data2<-data.frame(year_all,proport,lower,upper)
head(mam_data2)

# PLOT 
require(ggplot2)
ggplot(mam_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower)) + ylim(c(0,0.25))


#### Fish #### 

# Data 1970-2001 is from Wood et al. 2010

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_fish<-c(10, 10,5,10,13,14,18,19,24,20, 15, 14, 11, 10,15,20,30,15,23,19,33, 31, 32,33,29,46,37,34,33, 41,50, 42,35.202,20.12,10.03,20,20,19.91,15.07,20.08,19.96,10.03,20.07,0,15.05,25) # 
length(successes_fish) # 46
plot(successes_fish) # 

total_fish<-c(160, 214,181,185,211,182, 231, 222,296,264,265,303,302,346,360,447, 438, 404,503,475, 536, 534,547, 555,596,645,  695,713,747,815,871,760.74,832.02, 1130.43, 752.06,781.18, 887.295,996.165,855.915, 1007.72, 794.4,1091.16,1025.73, 1275.47, 1183.88, 929.67)
length(total_fish) # 46
plot(total_fish) # 
plot(successes_fish/total_fish)

fish_data<-data.frame(year_all,successes_fish,total_fish)

head(fish_data)

confint<-binom.confint(fish_data$successes_fish,fish_data$total_fish,method="exact")

proport<-successes_fish/total_fish
lower<-confint$lower
upper<-confint$upper

fish_data2<-data.frame(year_all,proport,lower,upper)
head(fish_data2)

# PLOT
require(ggplot2)
ggplot(fish_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower)) + ylim(c(0,0.25))


##### Corals ######

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_cor<-c(0,0, 0,0,1,2,0,0, 1,1, 0,  0,1,  3,3, 1, 1,0,0,1, 0, 2,0, 0, 0,2, 3, 0,4,4,2,3.72, 15.3,12.06, 3.96,17.82, 10.68,11.76,19.6,15.8, 24, 16.92,15.96, 21.12, 19.68, 9.48)
length(successes_cor) # 46 

total_cor<-c(12,18.8,27.3,21.92,  42.56,36.6,17.2,27.72,   58.24,59.64,23.04, 78.89655172,109.5849057,  50, 65.82539683,123.6666667,  60.17910448,79.30434783, 69.78461538,55.45205479, 111.2394366,64.81690141,59.65217391,   104.1818182,74.61538462,129.754717,145,134.46, 124.5049505, 164.2142857,174.1176471, 249.3654514,454.4052288, 482.9064748, 440,438.8027211, 615, 490.6950355, 521.4054054, 626.4792899,  600,636.2512315,  621.3121693,700.195,634.3453608, 713.05)

cor_data<-data.frame(year_all,successes_cor,total_cor)
head(cor_data)

confint_cor<-binom.confint(cor_data$successes_cor,cor_data$total_cor,method="exact")

proport_cor<-successes_cor/total_cor
lower_cor<-confint_cor$lower
upper_cor<-confint_cor$upper

cor_data2<-data.frame(year_all,proport_cor,lower_cor,upper_cor)
head(cor_data2)

# PLOT
require(ggplot2)
ggplot(cor_data2, aes(x = year_all, y = proport_cor)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper_cor, ymin = lower_cor)) + ylim(c(0,0.25))


##### Urchins - 0s in denominator make ugly plots, but still improved

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_urch<-c(1, 1,0,0, 0,0,1,1, 0,0,0,0, 0,1,  4,1,0,1,2,2,0,0,1,0,3,1,0,2,1,1,1,0.5,1,0,3,0,1,1,2,1,2,2,2,1,2,1)
length(successes_urch) # 46
plot(successes_urch) # 

total_urch<-c(104.64,92.12, 0, 0,0,0,158, 209.72,0,0, 0,0, 0,344.0142857, 254.0377358, 361.1066667,355, 298.0645161,325,175.78, 0,0,188.16, 0,168.56,198,160.72,232,200.64,182.4,160.32, 48.41,55.04,49.92,108.48,63.36,91.84,77.76,72.9,61.62,78.3,61.92,104.92,105.64,81.84,63.84)
length(total_urch) # 46
plot(total_urch) # 

urch_data<-data.frame(year_all,successes_urch,total_urch)

head(urch_data)

confint<-binom.confint(urch_data$successes_urch,urch_data$total_urch,method="exact")

proport<-successes_urch/total_urch
lower<-confint$lower
upper<-confint$upper

urch_data2<-data.frame(year_all,proport,lower,upper)
head(urch_data2)

# PLOT
require(ggplot2)
ggplot(urch_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))


##### Seagrass #####

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_sg<-c(0,0,0, 0,0, 0, 0,  0, 1, 0, 0,0,0, 0, 2,1, 0,0,1,1,0,3,  0,0, 0,1,0,0,0, 1,0,0,1,0,1,0,0,0,0,1,0,1,1,1,0,0) 
length(successes_sg) # 46

total_sg<-c(3, 8, 9,4,12,6,9, 16, 29, 44.28,37,39.5, 62.16,0,30.36,57.66,   0,51.68,45.26,39.68, 45.44,49.84, 40.8, 0,0,67.8,0,0,0,89.28,0,46.4,83.2,72.72,50.84,59.64,49.8,73.26,65.96,50.56,40.2,58.08,50.4,83.3,60.72,69.3)
length(total_sg) # 46

sg_data<-data.frame(year_all,successes_sg,total_sg)

head(sg_data)

confint<-binom.confint(sg_data$successes_sg,sg_data$total_sg,method="exact")

proport<-successes_sg/total_sg
lower<-confint$lower
upper<-confint$upper

sg_data2<-data.frame(year_all,proport,lower,upper)
head(sg_data2)

# PLOT
require(ggplot2)
ggplot(sg_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))

#### Turtles #### 

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_turt<-c(0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,1,1,1,2,1, 3, 3, 1,3,2,2,3,2,5,2,1,2,4,4,8,2,6,3,1,2,3,4)
length(successes_turt) # 46

total_turt<-c(3,1,4,4, 7,9,11,7, 13,24, 25,24,19,21,16,26,22,23,15,22,34, 24,28,28,34,32,23,30,34,32,45,64.74, 120.54,142.08,90.16, 95.06, 161,148.96,150,146.02,156.8,167, 200,163,188.16,155)
length(total_turt) # 46

turt_data<-data.frame(year_all,successes_turt,total_turt)

head(turt_data)

confint<-binom.confint(turt_data$successes_turt,turt_data$total_turt,method="exact")

proport<-successes_turt/total_turt
lower<-confint$lower
upper<-confint$upper

turt_data2<-data.frame(year_all,proport,lower,upper)
head(turt_data2)

# PLOT
require(ggplot2)
ggplot(turt_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))

#### DECAPODS ###

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_dec<-c(4, 2,5,7,8,10,6, 5,7,7,4,11,5,8,8,7, 8,9,10,4,7,8, 10,7, 7,5,12,4,16,6,9,10,9.92,7.76,16.48,10.88,18.9,21.56,12.64,25.2,4.44,15.12,10.04,4.98,10,10.04)
length(successes_dec) # 46

total_dec<-c(42.30, 75.00,115.00,105.00,155.42,160.81,155.77,100.48, 245.98,295.00,215.46,318.33,330.56,235.00, 265.00,331.09,419.35,299.50, 298.22,288.39,305.00, 259.57,145.00, 209.34,145.00,139.55,155.00,274.62,413.47,324.16,249.06,351.5241,370.7956989, 585.8357143,285.686747,370.7789474,461.6428571,459.0980392,479.1504425,550,480.8135593,440.7457627,481.536,646.4742857,760, 573.3211679)             
length(total_dec) # 46

dec_data<-data.frame(year_all,successes_dec,total_dec)

head(dec_data)

confint<-binom.confint(dec_data$successes_dec,dec_data$total_dec,method="exact")

proport<-successes_dec/total_dec
lower<-confint$lower
upper<-confint$upper

dec_data2<-data.frame(year_all,proport,lower,upper)
head(dec_data2)

# PLOT
require(ggplot2)
ggplot(dec_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))  + ylim(c(0,0.25))

#### Sharks & Rays ###

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_elasmo<-c(3, 0,5,1, 3,3, 3,5,4,3,2,4,1,4, 4,5,3,6,1, 3,5,3,0, 5,1,2,2,2,5,5, 6,5.5,6,3,2,0,3,0,3,3,4,3,2,3,3,2)
length(successes_elasmo) # 46

total_elasmo<-c(55.04, 55.46,72.9,53.1,79.2,62.72,79.2,88.2,72.16,87.12, 82.8,85.36,93.24,105.3,113.52,112.64,  109.48,111.86,84.8,73.92,109.22, 79.18,120.32,118.68, 100.8,95.04,121.36,105.6, 101.08, 116.18,103.32, 89.18, 126,152.88, 108.08, 114.26,95.2,144.3,168,192.36, 135.72,161.7,188.8125,218.5806452, 159.4754098, 234.3188406)             
length(total_elasmo) # 46

plot(successes_elasmo/total_elasmo)

elasmo_data<-data.frame(year_all,successes_elasmo,total_elasmo)

head(elasmo_data)

confint<-binom.confint(elasmo_data$successes_elasmo,elasmo_data$total_elasmo,method="exact")

proport<-successes_elasmo/total_elasmo
lower<-confint$lower
upper<-confint$upper

elasmo_data2<-data.frame(year_all,proport,lower,upper)
head(elasmo_data2)

# PLOT
require(ggplot2)
ggplot(elasmo_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower)) + ylim(c(0,0.25))

###### Molluscs #####
year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978, 1979, 1980,  1981,1982,1983,1984, 1985,1986,  1987,1988,1989,1990,1991,1992,1993,  1994,1995, 1996,1997,1998,1999, 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_moll<-c(1,5, 4,3,1, 2, 6, 7, 6,6,4,9,8,5,8,8,8,12.24,13.92,7,6, 9,  7, 12.96,8.82,3, 17.6,8.68,23.1,19.8,14,19.24,19.8,7.04,10.02,24.78,14,13.92,9.9644,4.6,4.68,4.92,10.03,9.9,0,19.875)
length(successes_moll) # 46

total_moll<-c(129.6486486, 170,183.752809,220,214.5612245,259.0877193,265,289.1714286,409,403.1899441,394.0595238,469.365,472.29,416.955,522.345,457.015,526.95,532.68,621.72,406.8,455.52,407.895,587.52, 487.035,492.96,523.55,596.845,537.61,532.335,530.075, 615.81,637.4325,489.4060606,581.5263158, 373.828125,340,344.4926471,395.5683453,918.68,524.4502618,540.5806452,543.8989899,626.9, 545.0,499.4413408,610.075)             
length(total_moll) # 46

plot(successes_moll/total_moll)

moll_data<-data.frame(year_all,successes_moll,total_moll)

head(moll_data)

confint<-binom.confint(moll_data$successes_moll,moll_data$total_moll,method="exact")

proport<-successes_moll/total_moll
lower<-confint$lower
upper<-confint$upper

moll_data2<-data.frame(year_all,proport,lower,upper)
head(moll_data2)

# PLOT
require(ggplot2)
ggplot(moll_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower)) + ylim(c(0,0.25))


###### Echinoderms #####
year_all<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

successes_echino<-c(4,1,2,5,1,1,3,5,8,8,2,3,4,4,4)

length(successes_echino) # 15

total_echino<-c(314.11,147.84,187.68,173.28,118.44, 105.12, 144.32, 153.64, 120.96,145.14,114.66,127.4, 195.36, 126.54, 124)             
length(total_echino) # 15

plot(successes_echino/total_echino)

echino_data<-data.frame(year_all,successes_echino,total_echino)

head(echino_data)

confint<-binom.confint(echino_data$successes_echino,echino_data$total_echino,method="exact")

proport<-successes_echino/total_echino
lower<-confint$lower
upper<-confint$upper

echino_data2<-data.frame(year_all,proport,lower,upper)
head(echino_data2)

# PLOT
require(ggplot2)
ggplot(echino_data2, aes(x = year_all, y = proport)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower)) + ylim(c(0,0.25))

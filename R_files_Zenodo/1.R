##### Correlations between normalized disease reports and time ####

#################
#### 1. Mammals #####
#################

# 2001-2013 ###
mmdis_3yearFWD<-c(0.017375097,0.016213597,0.018681748,0.017576451,0.021053609,0.026657943,0.027636655,0.038723573,0.043508628,0.049814363,0.042240354,0.037709174,0.034944085)

year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

cor.test(mmdis_3yearFWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_mamm<-data.frame(year13,mmdis_3yearFWD)

par(mfrow=c(1,1))

require(ggplot2)  # Plot the value * 100 (percent)
ggplot(data_mamm, aes(x = year13, y = mmdis_3yearFWD*100)) +
  geom_point(size = 4) + ylab("Mammal %") + xlab("Year") +
  ylim(0,8)

# 1970-2013
year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

mm_Disease_45years_3yearFWD<-c(0.022747042,0.016901729,0.011188128,0.015034281,0.018520673,0.020908561,0.029216468,0.026311819,0.028818258,0.030069023,0.02450499, 0.020449214,0.021359317,0.031335169,0.03634047,0.026676705,0.039010459,0.048399788,0.051156033,0.043761236,0.049107906,0.052136672,0.043285186,0.028767109,0.022278178,0.027725984,0.032295594,0.04033014,0.034424648,0.027417894,0.023346998, 0.019948689,0.016213597,0.018681748,0.017576451, 0.021053609,0.026657943,0.027636655, 0.038723573, 0.043508628, 0.049814363,0.042240354, 0.037709174,0.034944085)

cor.test(mm_Disease_45years_3yearFWD,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)


# PLOTS 
data_mamm2<-data.frame(year_all,mm_Disease_45years_3yearFWD)

mean_mamm<-mean(mm_Disease_45years_3yearFWD) #
sd_mamm<-sd(mm_Disease_45years_3yearFWD)

mean_mamm-sd_mamm # 0.0194607
mean_mamm+sd_mamm # 0.04138178

require(ggplot2)  #
ggplot(data_mamm2, aes(x = year_all, y = mm_Disease_45years_3yearFWD*100)) +
  geom_point(size = 4) + ylab("Mammals %") + xlab('Year') +ylim(0,8)+
  geom_hline(yintercept = 0.0194607*100) +
  geom_hline(yintercept = 0.04138178*100) # lines show peaks


#################
#### 2. Fish #####
#################

## 2001-2013 ####
year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

fidis_3yearFWD<-c(0.033857588,0.024479321, 0.018910033,0.020491353,0.022711424, 0.020044706,0.019174108,0.020885714, 0.018083543, 0.017962516, 0.009587056,0.010760068, 0.013201974) 

cor.test(fidis_3yearFWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_fish<-data.frame(year13,fidis_3yearFWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_fish, aes(x = year13, y = fidis_3yearFWD*100)) +
  geom_point(size = 4) + xlab("Year") + ylab("Fish %") +
  ylim(0,8)

### 1970-2013 ##
# N.B. fish data 1970-2001 is from Wood et al. 2010

year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

fiDisease_45years_3yearFWD<-c(0.04561776,0.042802445,0.047763246,  0.064196168,0.072152176,0.08014358,0.081529582,0.080808081,0.071147477,0.05952199,0.046410745, 0.037176732,0.035664081,0.038437043,0.051634182,0.050121531,0.05044917,0.040951453,0.049097603, 0.053206533,0.059373504,0.058670936,0.055539364, 0.059811669,0.057737653,0.057413691,0.04836665,0.047389763,0.050629579, 0.05409323,0.051427701,0.038224759,0.024479321,0.018910033,0.020491353,0.022711424,0.020044706,0.019174108,0.020885714,0.018083543,0.017962516,0.009587056,0.010760068, 0.013201974)

cor.test(fiDisease_45years_3yearFWD,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_fish2<-data.frame(year_all,fiDisease_45years_3yearFWD)

mean_fish<-mean(fiDisease_45years_3yearFWD) #
sd_fish<-sd(fiDisease_45years_3yearFWD)

mean_fish-sd_fish # 0.02483619
mean_fish+sd_fish # 0.06397299

require(ggplot2)  #
ggplot(data_fish2, aes(x = year_all, y = fiDisease_45years_3yearFWD*100)) +
  geom_point(size = 4) + ylab("Fish %") + xlab('Year') +
  ylim(0,8)+
  geom_hline(yintercept = 0.02483619*100) +
  geom_hline(yintercept = 0.06397299*100) # These lines show peaks


#################
#### 3. Coral #####
#################

## 2001-2013 ####
year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

cordis_3yearFWD<-c(0.025139969,0.022548056, 0.024861428,0.022325453, 0.027314121,0.026307523,0.028925675, 0.034270339,0.030604525,0.03076028, 0.027481288,0.028958235, 0.024827379)

cor.test(cordis_3yearFWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)


# plots
data_cor<-data.frame(year13,cordis_3yearFWD)

ggplot(data_cor, aes(x = year13, y = cordis_3yearFWD*100)) + 
  ylab("Coral %") + xlab("Year") +
  geom_point(size = 4) + ylim(0,8)

### 1970-2013 ##
year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

Disease_45years_3yearFWD<-c(0, 0,0.00783208,0.026047017, 0.026047017,0.018214936,0.005723443, 0.011312533,0.011312533,0.00558909, 0.003041781, 0.023041781, 0.038233486,0.037887123,0.023426144,0.008234439, 0.005539021,0.006011199, 0.006011199,0.016296582,0.010285383, 0.010285383,0,0.005137899,0.012034451,0.012034451,0.017605631, 0.018828551, 0.02265738,
                            
                            0.016654991,0.019758981,0.024254745, 0.022548056,0.024861428,0.022325453,0.027314121,0.026307523, 0.028925675,0.034270339,  0.030604525, 0.03076028, 0.027481288, 0.028958235,0.024827379)

cor.test(Disease_45years_3yearFWD,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95) 

# Check with slight changes so not exact ties
Disease_45years_3yearFWD_2<-c(0, 0.00000001,0.00783208,0.026047017, 0.02604702,0.018214936,0.005723443, 0.0113125,0.011312533,0.00558909, 0.003041781, 0.023041781, 0.038233486,0.037887123,0.023426144,0.008234439, 0.005539021,0.006011199, 0.0060112,0.016296582,0.010285383, 0.0102854,0.000000000001,0.005137899,0.012034451,0.01203445,0.017605631, 0.018828551, 0.02265738,0.016654991,0.019758981,0.024254745, 0.022548056,0.024861428,0.022325453,0.027314121,0.026307523, 0.028925675,0.034270339,  0.030604525, 0.03076028, 0.027481288, 0.028958235,0.024827379)

cor.test(Disease_45years_3yearFWD_2,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95) 

# formerly: rho = 0.551, p = 0.000106*
# now: rho=0.551, p= 0.000136 - ok

# Check with exact=FALSE
cor.test(Disease_45years_3yearFWD,year_all, method = "spearman",
         continuity = FALSE,exact=FALSE,
         conf.level = 0.95) # p = 0.000106, rho = 0.551: same 

# PLOTS 
data_cor2<-data.frame(year_all,Disease_45years_3yearFWD)

mean_cor<-mean(Disease_45years_3yearFWD) #
sd_cor<-sd(Disease_45years_3yearFWD)

mean_cor-sd_cor # 0.007077434
mean_cor+sd_cor # 0.02831

require(ggplot2)  #
ggplot(data_cor2, aes(x = year_all, y = Disease_45years_3yearFWD*100)) +
  geom_point(size = 4) + ylab("Coral %") + xlab("Year") +
  ylim(0,8) +
  geom_hline(yintercept = 0.007077434*100) +
  geom_hline(yintercept = 0.02831*100) # These lines show peaks 


#################
# 4. Urchins
##############

year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

urdis_3yearFWD<-c(0.009499016,0.015274491,0.009218289,  0.01284779, 0.007916195, 0.017061142, 0.018841141,0.023068708,0.024690341, 0.025634889, 0.020275999,0.017655394, 0.016522733)

cor.test(urdis_3yearFWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_urch<-data.frame(year13,urdis_3yearFWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_urch, aes(x = year13, y = urdis_3yearFWD*100)) +
  geom_point(size = 4) + ylab("Urchin %") + xlab("Year")  + ylim(0,8) 

### 1970-2013 ##
year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

urDisease_45years_3yearFWD<-c(0.006803994,0.003618469,0,0,0.002109705,0.003699125,0.003699125,0.001589421,0,0,0,0.000968952, 0.006217516,0.007140604, 0.006171652, 0.002041414, 0.003169608,0.006962228, 0.005843902,0.00379262, 0.001771542, 0.001771542,0.007704148,0.007616107,0.007616107, 0.004557065,0.004534914,0.006362399,0.005568011,0.005628068, 0.009856784, 0.007777609, 0.015274491,  0.009218289,  0.01284779,  0.007916195,   0.017061142,  0.018841141,0.023068708,  0.024690341, 0.025634889, 0.020275999, 0.017655394, 0.016522733)
 
cor.test(urDisease_45years_3yearFWD,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_urch2<-data.frame(year_all,urDisease_45years_3yearFWD)

mean_urch<-mean(urDisease_45years_3yearFWD) #
sd_urch<-sd(urDisease_45years_3yearFWD)

mean_urch-sd_urch # 0.0007857823
mean_urch+sd_urch # 0.01483239

require(ggplot2)  #
ggplot(data_urch2, aes(x = year_all, y = urDisease_45years_3yearFWD*100)) +
  geom_point(size = 4) + ylab("Urchin %") + xlab("Year") +
  ylim(0,8)+
  geom_hline(yintercept = 0.0007857823*100) +
  geom_hline(yintercept = 0.01483239*100) # These lines show peaks


#################
# 5. Echinoderms 
################

year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

echino_3FWD<-c(0.010051597,0.015425179,0.015984854,0.015603688,0.01291439,0.020947895,0.039822771,0.05126679,0.046233212,0.03203665,0.020488592,0.025211153,0.028114548)

cor.test(echino_3FWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_echino<-data.frame(year13,echino_3FWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_echino, aes(x = year13, y = echino_3FWD*100)) +
  geom_point(size = 4) +  ylim(0,8) + ylab("Echinoderm %") + xlab("Year")
 

###############
# 6. Seagrass 
##############

year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

sgdis_3yearFWD<-c(0.00400641,0.010562927,0.006556517,0.006556517,0, 0, 0.006592827,0.006592827, 0.012332037,0.012352967, 0.016354568, 0.010615357, 0.004001601)
 
cor.test(sgdis_3yearFWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95) 

# PLOTS 
data_sg<-data.frame(year13,sgdis_3yearFWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_sg, aes(x = year13, y = sgdis_3yearFWD*100)) +
  ylab("Seagrass %") + xlab("Year") +
  geom_point(size = 4) + ylim(0,8)


### 1970-2013 ##
year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

sg_Disease_45years_3yearFWD<-c(0, 0,0,0,0,0, 0.011494253, 0.011494253,0.011494253, 0,0, 0,0.021958718, 0.027739733, 0.027739733, 0.005781015,0.007364855, 0.015765393,0.015765393,0.028464743, 0.020064205,0.020064205,0,0.004916421,0.004916421, 0.004916421,0,0.003733572,0.003733572,0.003733572,0.00400641,0.00400641,0.010562927, 0.006556517,  0.006556517, 0,0,0.006592827,0.006592827, 0.012332037,  0.012352967,0.016354568,   0.010615357, 0.004001601)

cor.test(sg_Disease_45years_3yearFWD,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# Plots
data_sg2<-data.frame(year_all,sg_Disease_45years_3yearFWD)

par(mfrow=c(1,1))

mean_sg<-mean(sg_Disease_45years_3yearFWD) #
sd_sg<-sd(sg_Disease_45years_3yearFWD)

mean_sg-sd_sg # -0.0002743467
mean_sg+sd_sg # 0.01625942

require(ggplot2)  #
ggplot(data_sg2, aes(x = year_all, y = sg_Disease_45years_3yearFWD*100)) +
  geom_point(size = 4) + ylim(0,8) + ylab("Seagrass %") + xlab("Year") + 
  geom_hline(yintercept = -0.0002743467*100) +
  geom_hline(yintercept = 0.01625942*100) #  lines show peaks

##############
# 7. Turtles
#############
year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

turt_3FWD<-c(0.026227773,0.022215992,0.015402438, 0.018991819, 0.024245637,0.0350103, 0.031294311,0.035098464,0.023308711,0.020409793,0.011744670,0.011071272,0.018006756)

cor.test(turt_3FWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_turt<-data.frame(year13,turt_3FWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_turt, aes(x = year13, y = turt_3FWD*100)) +
  ylab("Turtle %") + xlab("Year") +
  geom_point(size = 4) + ylim(0,8)


### 1970-2013 ##
year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

turt_45<-c(0,0,0,0,0.03030303, 0.03030303, 0.03030303,0,0,0.013888889, 0.031432749, 0.031432749,0.01754386,0.012820513,0.027972028,0.042464782,0.029644269, 0.029644269, 0.024955437,0.038844326, 0.047502334, 0.049603175,0.06512605,0.072566527, 0.075154518, 0.079076087, 0.06743393,0.07377451,0.062663399,0.054661949,0.047655285,0.030125255,0.022215992, 0.015402438, 0.018991819,0.024245637,0.0350103,0.031294311,0.035098464, 0.023308711, 0.020409793,0.01174467,0.011071272, 0.018006756)

cor.test(turt_45,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_turt2<-data.frame(year_all,turt_45)

mean_turt<-mean(turt_45) #
sd_turt<-sd(turt_45)

mean_turt-sd_turt # 0.009229683
mean_turt+sd_turt # 0.0536656

require(ggplot2)  # Graph for if it's NS
ggplot(data_turt2, aes(x = year_all, y = turt_45*100)) +
  geom_point(size = 4) +  ylab("Turtle %") + xlab("Year") +
  geom_hline(yintercept = 0.009229683*100) +
  geom_hline(yintercept = 0.0536656*100) # These lines show peaks

##############
# 8. Racoon rabies
#############

# Records from CDC not included because not publicly available. 

### Literature reports: 2001-2013

year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)
lit<-c(0.088311199,0.046644532,0.051359751,0.022222222,0.032026144,0.052320928,0.052320928,0.06335034,0.029513889, 0.045616948, 0.045359747,0.045359747,0.029256687) # 13

cor.test(lit,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_lit<-data.frame(year13,lit)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_lit, aes(x = year13, y = lit)) +ylab("Rabies %") + xlab("Year") +
  geom_point(size = 4) 


### 2001-2013 records vs. lit code and result (not runable because CDC records not included)

cor.test(lit,records, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)
# p = 0.042, rho = 0.57


### Literature reports 1970-2013
Year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

rabies_all_lit<-c(0.034100597,0.036714976,0.022222222,0.022222222,0,0,0,0.015873016,  0.015873016,0.028693529,0.012820513,0.036630037,0.034562212,0.047853015,0.034460158,0.036053148,0.045750851,0.035334185,0.053291536,0.062561095,0.092864125,0.062561095,0.052525253,0.031746032,0.031746032,0.033333333,0.06547619,0.09047619,0.092982456,0.092982456,0.085489459,0.088311199,0.046644532, 0.051359751, 0.022222222,0.032026144,0.052320928, 0.052320928, 0.06335034,  0.029513889,0.045616948, 0.045359747, 0.045359747,0.029256687)

cor.test(rabies_all_lit,Year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)


### Correlation between records & lit 1970-2013 (not runable because CDC records not included)

cor.test(rabies_all_lit,records_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)
# rho = 0.58, p = 3.881e-05

# PLOTS 
data_lit45<-data.frame(Year_all,rabies_all_lit)

mean_lit<-mean(rabies_all_lit) #
sd_lit<-sd(rabies_all_lit)

mean_lit-sd_lit # 0.01862608
mean_lit+sd_lit # 0.06823129

require(ggplot2)  #
ggplot(data_lit45, aes(x = year_all, y = rabies_all_lit*100)) +
  geom_point(size = 4) + ylab("Rabies %") + xlab("Year") +
  geom_hline(yintercept = 0.01862608*100) +
  geom_hline(yintercept = 0.06823129*100) # These lines show peaks


##############
# 9. Decapods 
#############
Year_13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

dec_3FWD<-c(0.018161943,0.032561624,0.033425075,0.042656643,0.039082005, 0.038094137,0.039719951,0.027144185,0.029786007,0.021463262,0.020953504,0.013904305,0.0127888) # 

cor.test(dec_3FWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_dec<-data.frame(year13,dec_3FWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_dec, aes(x = year13, y = dec_3FWD*100)) +
  ylab("Decapod %") + xlab("Year") +
  geom_point(size = 4) + ylim(0,8)


### 1970-2013 ##
Year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

dec_45<-c(0.054902525,0.045603865,0.053872904,0.06010853, 0.050726127, 0.050154929,0.038912416,0.033982201, 0.023583854,0.025616418,0.022748772,0.027907907,0.02645237,0.028457821,0.023469292,0.023423234,0.027553172,0.025817569,0.023451007,0.022546912,0.040912031, 0.0444077, 0.050226402,0.039180876,0.053841385,0.042604567,0.043560474,0.023923894,0.031113927,0.028706908,0.031454796,0.023824971,0.032561624,  0.033425075, 0.042656643,0.039082005, 0.038094137,0.039719951,0.027144185,0.029786007,0.021463262,0.020953504,0.013904305, 0.0127888) # 

cor.test(dec_45,Year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_dec2<-data.frame(Year_all,dec_45)

mean_dec<-mean(dec_45) #
sd_dec<-sd(dec_45)

mean_dec-sd_dec # 0.02228636
mean_dec+sd_dec # 0.04565133

require(ggplot2)  #
ggplot(data_dec2, aes(x = Year_all, y = dec_45*100)) +
  geom_point(size = 4) + ylab("Decapod %") + xlab("Year") +
  ylim(0,8)+
  geom_hline(yintercept = 0.02228636*100) +
  geom_hline(yintercept = 0.04565133*100) # These lines show peaks

##############
# 10. Elasmobranchs 
#############
year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

shark_3FWD<-c(0.039324267,0.028582364,0.012709348,0.016672472,0.010504202, 0.016456583,0.011150967,0.020975115, 0.021207026, 0.019539279,  0.014290102,0.014376369,0.013690656)  

cor.test(shark_3FWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_elasmo<-data.frame(year13,shark_3FWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_elasmo, aes(x = year13, y = shark_3FWD*100)) + ylab("Elasmobranch %") + xlab("Year") +
  geom_point(size = 4) + ylim(0,8)

### 1970-2013 ##
Year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

shark_45<-c(0.041030973,0.029139832,0.041766095, 0.034847604,0.041196403,0.047466588,0.050000168,0.048852326, 0.038007408,0.035150069, 0.027246652, 0.031857357,0.027982599,0.039203997,0.035675851, 0.041809982,0.030944398,0.035338448, 0.032718677, 0.041417311,0.027889172, 0.026672818,0.017350244,0.024364835,0.015814767, 0.01882102,0.028295019,0.037147277,0.050191482,0.054689307,0.056216767,0.043400509,0.028582364,0.012709348, 0.016672472,0.010504202,0.016456583, 0.011150967,0.020975115, 0.021207026, 0.019539279,0.014290102,0.014376369,0.013690656)# 

cor.test(shark_45,year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_elasmo2<-data.frame(Year_all,shark_45)

mean_elasmo<-mean(shark_45) #
sd_elasmo<-sd(shark_45)

mean_elasmo-sd_elasmo # 0.01810295
mean_elasmo+sd_elasmo # 0.04338162

require(ggplot2)  #
ggplot(data_elasmo2, aes(x = Year_all, y = shark_45*100)) +
  geom_point(size = 4) + ylab("Elasmobranch %") + xlab("Year") +
  ylim(0,8)+
  geom_hline(yintercept = 0.01810295*100) +
  geom_hline(yintercept = 0.04338162*100) # These lines show peaks


##############
# 11. Molluscs 
#############
year13<-c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)

mollusc_3FWD<-c(0.03279862,0.026455679,0.037264063, 0.046775196, 0.049570567,0.028892332,0.018269537,0.00942536, 0.008830769, 0.011242179,0.014411438, 0.011390152, 0.016914366)   

cor.test(mollusc_3FWD,year13, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_moll<-data.frame(year13,mollusc_3FWD)

par(mfrow=c(1,1))

require(ggplot2)  #
ggplot(data_moll, aes(x = year13, y = mollusc_3FWD*100)) +
  ylab("Mollusc %") + xlab("Year") +
  geom_point(size = 4) + ylim(0,8)


### 1970-2013 ##
Year_all<-c(1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010,2011,2012,2013)

moll_45<-c(0.019631098,0.021605501,0.013355138, 0.008672144, 0.011673859,0.018189332,0.020506177,0.017919448,0.013234,0.014735639, 0.015421446, 0.016035097, 0.014748665, 0.014937382,0.016000717,0.018554917, 0.020183118, 0.020858374,0.017589578,0.017481245,0.015716917,0.020196329, 0.018805468,0.016744009, 0.017703474, 0.017121346, 0.029675885, 0.032297488,0.034493739,0.030285465,0.031320131,0.027777393,0.026455679,0.037264063,0.046775196,0.049570567, 0.028892332,0.018269537,0.009425365,0.008830769,0.011242179,0.014411438,0.011390152,0.016914366) #

cor.test(moll_45,Year_all, method = "spearman",
         continuity = FALSE,
         conf.level = 0.95)

# PLOTS 
data_moll<-data.frame(Year_all,moll_45)

mean_moll<-mean(moll_45) #
sd_moll<-sd(moll_45)

mean_moll-sd_moll # 0.01128301
mean_moll+sd_moll # 0.02975846

require(ggplot2)  #
ggplot(data_moll, aes(x = Year_all, y = moll_45*100)) +
  geom_point(size = 4) + ylab("Mollusc %") + xlab("Year") +
  ylim(0,8)+
  geom_hline(yintercept = 0.01128301*100) +
  geom_hline(yintercept = 0.02975846*100) # These lines show peaks


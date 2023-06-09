##### Libraries and Preliminaries####

# Load the required libraries
library(tidyverse)


### Viscosity ####

#calculating sucrose solution viscosities
#create a tibble of the water data

water.dat <- read_csv("Data/water_dat.csv",col_types=list(col_double(),col_double(),
                                                          col_double(),col_factor()))
# Make a table with temps from -30C to 25C
suc.dat <- tribble(~temp,-30,-25,-20,-10,-5,0,5,10,15,20,25)

#same as previous line  with a wider range
suc.dat <- tibble(temp = -40:25)

#The temperature relation. Phi is used to calculate viscosity at diff. temperatures made above
suc.dat$phi <- (30-suc.dat$temp)/(91+suc.dat$temp)
## Changed to final concentration- 17.6 in summer and 35.2 in winter

#Calculate the mole fraction of sucrose in a 17% solution
mol <- 17.1/342.29648
mol.w <- 82.9/18.01528
mol.suc16<- mol/(mol+mol.w)
#calculate the mole fraction of sucrose in a 34% solution
mol.s30 <- 34.2/342.29648
mol.w30 <- 65.8/18.01528
mol.suc30<- mol.s30/(mol.s30+mol.w30)

mol.s25 <- 25/342.29648
mol.w25 <- 75/18.01528
mol.suc25<- mol.s25/(mol.s25+mol.w25)


suc.dat$vis.16 <- 10^(22.46*mol.suc16-0.114 + suc.dat$phi*(1.1 + 43.2*(mol.suc16^1.25)))

suc.dat$vis.30 <- 10^(22.46*mol.suc30-0.114 + suc.dat$phi*(1.1 + 43.2*(mol.suc30^1.25)))

#suc.dat$vis.25 <- 10^(22.46*mol.suc25-0.114 + suc.dat$phi*(1.1 + 43.2*(mol.suc25^1.25)))


#rename data columns
suc.dat <- rename(suc.dat, "vis.16"=vis.16, "vis.30"=vis.30)#, "25"=vis.25)
#convert to long from wide
#suc.long <- gather(suc.dat, "15", "30", "25", key="concentration",value= "viscosity")
suc.long <- gather(suc.dat, "vis.16", "vis.30", key="concentration",value= "viscosity")

#merge water and sucrose data set
suc.long <- full_join(suc.long,water.dat)

#save the data file

write_csv(suc.long, "Data/suc_plot_data_long_granular.csv")

#### RESISTANCE SECTION ####

#### Single Element Resistance 17% ####
#calculate the lumen resistance of one element for each location (top, mid, base)
# formula is 8*viscosity*Length/pi* radius^4

suc.dat$lumenres.base <- ((8 * suc.dat$vis.16 * 438.1581) / (pi*(27.97065^4)))

suc.dat$lumenres.mid <- ((8 * suc.dat$vis.16 * 345.5351) / (pi*(24.06255^4)))

suc.dat$lumenres.tip <- ((8 * suc.dat$vis.16 * 261.7154) / (pi*(9.9448^4)))


# Calculate the plate parameters for resistance

beta.base <- 1.24/0.21
beta.mid <- 1.37/0.25
beta.tip <- 0.55/0.08

alpha.base <- ((8*1)/(3*pi*1.24))
alpha.mid <- ((8*1)/(3*pi*1.37))
alpha.tip <- ((8*1)/(3*pi*0.55))

suc.dat$rp.base <- (((3*suc.dat$vis.16)/1.24^3)*1/208.7294068)*((1/(1+3*beta.base^2)))+(alpha.base/(1+6*beta.base^2+3*beta.base^4))
suc.dat$rp.mid <- (((3*suc.dat$vis.16)/1.37^3)*1/195.3674326)*((1/(1+3*beta.mid^2)))+(alpha.mid/(1+6*beta.mid^2+3*beta.mid^4))
suc.dat$rp.tip <-  (((3*suc.dat$vis.16)/0.55^3)*1/192.9410794)*((1/(1+3*beta.tip^2)))+(alpha.tip/(1+6*beta.tip^2+3*beta.tip^4))


#### Total Resistance for 17% sucrose ####

#total resistance in a single tube -- this is lumen resistance + plate resistance

suc.dat$r.base <- suc.dat$lumenres.base + suc.dat$rp.base
suc.dat$r.mid <- suc.dat$lumenres.mid + suc.dat$rp.mid
suc.dat$r.tip <- suc.dat$lumenres.tip + suc.dat$rp.tip

# total resistance in 1m this is the total element resistance * the # of cells in 1m (from jessica's data)

suc.dat$rbase.1m <- suc.dat$r.base*2282.281213
suc.dat$rmid.1m <- suc.dat$r.mid*2894.061993
suc.dat$rtip.1m <- suc.dat$r.tip*3820.94443

# total resistance for a 27m Quercus rubra assuming that radius is constant for each third
# Each 1m resistance calculated above is multiplied by 9

suc.dat$total.res <- (suc.dat$rbase.1m*9) + (suc.dat$rmid.1m*9) + (suc.dat$rtip.1m*9)



##### Resistance with 17% sucrose and 50% pore occlusion #####

# now do the same thing for resistance calculation at 16%, BUT pores have been partially clogged with extra callose in the system
# calculate the lumen resistance of one element for each location (top, mid, base)
# formula is 8*viscosity*Length/pi* radius^4


beta.baseOc <- (1.24*.5)/0.21
beta.midOc <- (1.37*.5)/0.25
beta.tipOc <- (0.55*.5)/0.08

alpha.baseOc <- ((8*1)/(3*pi*(1.24*.5)))
alpha.midOc <- ((8*1)/(3*pi*(1.37*.5)))
alpha.tipOc <- ((8*1)/(3*pi*(0.55*.5)))

suc.dat$rp.baseOc <- (((3*suc.dat$vis.16)/(1.24*.5)^3)*1/208.7294068)*((1/(1+3*beta.baseOc^2)))+(alpha.baseOc/(1+6*beta.baseOc^2+3*beta.baseOc^4))
suc.dat$rp.midOc <- (((3*suc.dat$vis.16)/(1.37*.5)^3)*1/195.3674326)*((1/(1+3*beta.midOc^2)))+(alpha.midOc/(1+6*beta.midOc^2+3*beta.midOc^4))
suc.dat$rp.tipOc <-  (((3*suc.dat$vis.16)/(0.55*.5)^3)*1/192.9410794)*((1/(1+3*beta.tipOc^2)))+(alpha.tipOc/(1+6*beta.tipOc^2+3*beta.tipOc^4))

#total resistance in a single element -- this is lumen resistance + plate resistance

suc.dat$r.baseOc <- suc.dat$lumenres.base + suc.dat$rp.baseOc
suc.dat$r.midOc <- suc.dat$lumenres.mid + suc.dat$rp.midOc
suc.dat$r.tipOc <- suc.dat$lumenres.tip + suc.dat$rp.tipOc

# total resistance in 1m this is the total element resistance * the # of cells in 1m (from jessica's data)

suc.dat$rbase.1mOc <- suc.dat$r.baseOc*2282.281213
suc.dat$rmid.1mOc <- suc.dat$r.midOc*2894.061993
suc.dat$rtip.1mOc <- suc.dat$r.tipOc*3820.94443

# total resistance for a 27m Quercus rubra assuming that radius is constant for each third
# Each 1m resistance calculated above is multiplied by 9

suc.dat$total.resOc <- (suc.dat$rbase.1mOc*9) + (suc.dat$rmid.1mOc*9) + (suc.dat$rtip.1mOc*9)



##### Resistance for 30% sucrose solution #####

# Calculate the lumen resistances using a 30 percent sucrose solution - all other parameters are the same

suc.dat$lumenres.base30 <- ((8 * suc.dat$vis.30 * 438.1581) / (pi*(27.97065^4)))

suc.dat$lumenres.mid30 <- ((8 * suc.dat$vis.30 * 345.5351) / (pi*(24.06255^4)))

suc.dat$lumenres.tip30 <- ((8 * suc.dat$vis.30 * 261.7154) / (pi*(9.9448^4)))


# Calculate the plate parameters for resistance

# alpha and beta are the same for this part as only viscosity is changing

suc.dat$rp.base30 <- (((3*suc.dat$vis.30)/1.24^3)*(1/208.7294068))*(((1/(1+3*beta.base^2)))+(alpha.base/(1+6*beta.base^2+3*beta.base^4)))
suc.dat$rp.mid30 <- (((3*suc.dat$vis.30)/1.37^3)*(1/195.3674326))*(((1/(1+3*beta.mid^2)))+(alpha.mid/(1+6*beta.mid^2+3*beta.mid^4)))
suc.dat$rp.tip30 <-  (((3*suc.dat$vis.30)/0.55^3)*(1/192.9410794))*(((1/(1+3*beta.tip^2)))+(alpha.tip/(1+6*beta.tip^2+3*beta.tip^4)))

# Total Resistance for 30% sucrose #
# total resistance in a single tube -- this is lumen resistance + plate resistance

suc.dat$r.base30 <- suc.dat$lumenres.base30 + suc.dat$rp.base30
suc.dat$r.mid30 <- suc.dat$lumenres.mid30 + suc.dat$rp.mid30
suc.dat$r.tip30 <- suc.dat$lumenres.tip30 + suc.dat$rp.tip30

# total resistance in 1m this is the total element resistance * the # of cells in 1m (from jessica's data)

suc.dat$rbase30.1m <- suc.dat$r.base30*2282.281213
suc.dat$rmid30.1m <- suc.dat$r.mid30*2894.061993
suc.dat$rtip30.1m <- suc.dat$r.tip30*3820.94443

# total resistance for a 27m Quercus rubra assuming that radius is constant for each third
# Each 1m resistance calculated above is multiplied by 9

suc.dat$total.res30 <- (suc.dat$rbase30.1m*9) + (suc.dat$rmid30.1m*9) + (suc.dat$rtip30.1m*9)


#### Resistance at 34% sucrose and 50% occluded pores #####


suc.dat$rp.baseOc30 <- (((3*suc.dat$vis.30)/(1.24*.5)^3)*1/208.7294068)*((1/(1+3*beta.baseOc^2)))+(alpha.baseOc/(1+6*beta.baseOc^2+3*beta.baseOc^4))
suc.dat$rp.midOc30 <- (((3*suc.dat$vis.30)/(1.37*.5)^3)*1/195.3674326)*((1/(1+3*beta.midOc^2)))+(alpha.midOc/(1+6*beta.midOc^2+3*beta.midOc^4))
suc.dat$rp.tipOc30 <-  (((3*suc.dat$vis.30)/(0.55*.5)^3)*1/192.9410794)*((1/(1+3*beta.tipOc^2)))+(alpha.tipOc/(1+6*beta.tipOc^2+3*beta.tipOc^4))



#total resistance in a single tube -- this is lumen resistance + plate resistance

suc.dat$r.baseOc30 <- suc.dat$lumenres.base30 + suc.dat$rp.baseOc30
suc.dat$r.midOc30 <- suc.dat$lumenres.mid30 + suc.dat$rp.midOc30
suc.dat$r.tipOc30 <- suc.dat$lumenres.tip30 + suc.dat$rp.tipOc30

# total resistance in 1m this is the total element resistance * the # of cells in 1m (from jessica's data)
base=1000000/438.1581
mid=1000000/345.5351
tip=9000000/261.754

suc.dat$rbase.1mOc30 <- suc.dat$r.baseOc30*2282.281213
suc.dat$rmid.1mOc30 <- suc.dat$r.midOc30*2894.061993
suc.dat$rtip.1mOc30 <- suc.dat$r.tipOc30*3820.94443

# total resistance for a 27m Quercus rubra assuming that radius is constant for each third
# Each 1m resistance calculated above is multiplied by 9

suc.dat$total.resOc30 <- (suc.dat$rbase.1mOc30*9) + (suc.dat$rmid.1mOc30*9) + (suc.dat$rtip.1mOc30*9)


#### Flow Rate ####

#pressure from Savage et al 2017 is 0.7MPa, but resistance is in mPa, so convert to mPa 0.7MPa = 0.7 * 1000000000

# Calculate flow rate using Q = deltaP/R, where R is the sum of the resistances calculated for the lumen and the plates

#Volumetric flow rate 
suc.dat$Q.16 = ((0.7*1000000000)/suc.dat$total.res) #17% sucrose
suc.dat$Q.16winter = ((0.7*1000000000)/4)/suc.dat$total.res #17% sucrose, occluded pores, 1/4 pressure gradient
suc.dat$Q.16occ = ((0.7*1000000000)/suc.dat$total.resOc)#17% sucrose with occluded pores
#flow rate for 30% viscosities - no radius change
suc.dat$Q.30 = (((0.7*1000000000)/suc.dat$total.res30))#34% sucrose, no occlusion
suc.dat$Q.30winter = (((0.7*1000000000)/4)/suc.dat$total.resOc30)#34% sucrose, occluded pores, 1/4 pressure gradient
suc.dat$Q.30occ = (0.7*1000000000)/suc.dat$total.resOc30 #34% sucrose, occluded pores


#Data converted into linear sap velocity by multiplying by the area of a midpoint sieve element

suc.dat$v.16 = ((0.7 * 1000000000)/((suc.dat$total.res)*(pi*24.06255^2))) #17% sucrose
suc.dat$v.16winter = ((0.7 * 1000000000)/4)/((suc.dat$total.res)*(pi*24.06255^2)) #17% sucrose, occluded pores, 1/4 pressure gradient
suc.dat$v.16occ = ((0.7 * 1000000000)/(suc.dat$total.resOc)*(pi*24.06255^2))#17% sucrose with occluded pores
#flow rate for 30% viscosities - no radius change
suc.dat$v.30 = ((0.7 * 1000000000)/((suc.dat$total.res30)*(pi*24.06255^2)))#34% sucrose, no occlusion
suc.dat$v.30winter = (((0.7 * 1000000000)/4)/((suc.dat$total.resOc30)*(pi*24.06255^2)))#34% sucrose, occluded pores, 1/4 pressure gradient
suc.dat$v.30occ = ((0.7 * 1000000000)/((suc.dat$total.resOc30)*(pi*24.06255^2)))#34% sucrose, occluded pores




#write the data to a data file

write_csv(suc.dat, "~/desktop/sucdat.csv")




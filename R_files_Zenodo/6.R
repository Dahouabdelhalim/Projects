####### Repeat authors 2001-2013 #####

### Corals ### 

year<-c(2001,2002,2003,2004, 2005,2006, 2007, 2008,2009,2010,2011, 2012,2013)

corals_full_3yr<-c(0.025139969,0.022548056,0.024861428, 0.022325453,0.027314121, 0.026307523, 0.028925675,0.034270339, 0.030604525, 0.03076028, 0.027481288, 0.028958235,0.024827379)
  
corals_no_common_3yr<-c(0.023892918, 0.021301004, 0.023357335, 0.018891821, 0.022739251,0.021446712, 0.025994402,0.032480306, 0.029127121, 0.029282876, 0.026003885, 0.028958235,0.024827379)

## Tests
cor.test(year,corals_full_3yr,method="spearman",conf.level=0.95) # NS after holm's bonferonni (below)

cor.test(year,corals_no_common_3yr,method="spearman",conf.level=0.95) #NS after holm's bonferonni

# Bonferonni with coral common author removed
shapiro_p<- c(0.0006517, 0.0003435, 0.01832, 0.01958, 0.016, 0.1692, 0.2094, 0.06109, 0.3633, 0.01115)
# Names: Mammals ,Fish,Coral,Urchins,Echinoderms,Seagrass,Turtles,Decapods,Elasmobranchs,Molluscs)

p.adjust(shapiro_p, method = "holm", n = 10) 
# corals: 0.11

##### Molluscs ####
year<-c(2001,2002,2003,2004, 2005,2006, 2007, 2008,2009,2010,2011, 2012,2013)

moll_full_3yr<-c(0.03279862,0.026455679,0.037264063,0.046775196,0.049570567,  0.028892332, 0.018269537, 0.009425365, 0.008830769, 0.011242179, 0.014411438, 0.011390152,  0.016914366)

moll_no_common_3yr<-c(0.03279862,0.026455679, 0.037264063, 0.046775196, 0.046638077, 0.025959842,0.015337047, 0.009425365,   0.008830769,0.011242179,0.008356393,0.005335106,0.00814449)

moll_3common_3yr<-c(0.016697942,0.015447511,0.035246384,0.046775196,0.046638077,0.024151901,0.013529106,0.007617424, 0.008830769, 0.011242179, 0.008356393, 0.005335106,0.00814449)

## Tests
cor.test(year,moll_full_3yr,method="spearman",conf.level=0.95) # rho = -0.692, p=0.01115

cor.test(year,moll_no_common_3yr,method="spearman",conf.level=0.95) # rho = -0.879, p= 1.901 e=05

cor.test(year,moll_3common_3yr,method="spearman",conf.level=0.95) # rho= -0.780, P=0.002621

# Bonferonni with common author removed
shapiro_p<- c(0.0006517, 0.0003435, 0.08122, 0.01958, 0.016, 0.1692, 0.2094, 0.06109, 0.3633, 1.901e-05)
# Names: Mammals ,Fish,Coral,Urchins,Echinoderms,Seagrass,Turtles,Decapods,Elasmobranchs,Molluscs)

p.adjust(shapiro_p, method = "holm", n = 10) # molluscs = 1.9e-04 

# Bonferonni with 3 most common removed
shapiro_p<- c(0.0006517, 0.0003435, 0.08122, 0.01958, 0.016, 0.1692, 0.2094, 0.06109, 0.3633,  0.002621)
# Names: Mammals ,Fish,Coral,Urchins,Echinoderms,Seagrass,Turtles,Decapods,Elasmobranchs,Molluscs

p.adjust(shapiro_p, method = "holm", n = 10) # molluscs = 0.021


### Elasmobranchs

year<-c(2001,2002,2003,2004, 2005,2006, 2007, 2008,2009,2010,2011, 2012,2013)

elasm_full_3yr<-c(0.039324267,0.028582364,0.012709348,0.016672472,0.010504202, 0.016456583,0.011150967,0.020975115, 0.021207026, 0.019539279,  0.014290102,0.014376369,0.013690656) 

elasm_no_common_3yr<-c(0.029178163,0.028582364,0.012709348,0.013171072,0.007002801,0.012955182,0.011150967, 0.020975115, 0.017084164, 0.011885578,0.006636401, 0.01084553, 0.013690656) # 2 most common authors removed (tied)

## Tests
cor.test(year,elasm_full_3yr,method="spearman",conf.level=0.95) # p = 0.3633, rho = -0.275

cor.test(year,elasm_no_common_3yr,method="spearman",conf.level=0.95) # p=0.1459, rho = -0.4285714 

# Bonferonni with 2 most common removed
shapiro_p<-  c(0.0006517, 0.0003435, 0.08122, 0.01958, 0.016, 0.1692, 0.2094, 0.06109, 0.1459, 0.01115)
# Names: Mammals ,Fish,Coral,Urchins,Echinoderms,Seagrass,Turtles,Decapods,Elasmobranchs,Molluscs

p.adjust(shapiro_p, method = "holm", n = 10) # molluscs = 0.023589
# now p = 0.44

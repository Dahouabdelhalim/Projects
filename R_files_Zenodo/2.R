#### Holm's bonferonni adjustment for results of correlation tests ###

#### 2001-2015 ####
# n= 10 
shapiro_p<- c(0.0006517, 0.0003435, 0.08122, 0.01958, 0.016, 0.1692, 0.2094, 0.06109, 0.3633, 0.01115)
# Order: Mammals ,Fish,Coral,Urchins,Echinoderms,Seagrass,Turtles,Decapods,Elasmobranchs,Molluscs

p.adjust(shapiro_p, method = "holm", n = 10) 
# Result: still significant are fish, mammals


#### 1970-2015 #####
# n=9, no echinoderms

shapiro.p2<-c(0.03616, 0.000001973, 0.000106, 4.33E-12, 0.1825, 0.1229, 0.01771, 0.00003511, 0.267)
# Order: Mammals ,Fish,Coral,Urchins,Seagrass,Turtles,Decapods,Elasmobranchs,Mollusc

p.adjust(shapiro.p2, method = "holm", n = 9) 


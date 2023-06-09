#African wild dog feeding priority

#do African wild dogs mainly approach carcasses when they have a chance at joining?
#data derived from "Approach by relativePFQ.xls"

#using all data (from 179 approach clusters)
prop.test (c(86,90,57), c(600, 207, 1024), conf.level = 0.95)
#higher up queue (+1)=better access. 86/600
#same age most access 90/207
#this could be because they are rejoining having just chased another away
#3-sample test for equality of proportions without continuity correction
#data:  c(86, 90, 57) out of c(600, 207, 1024)
#X-squared = 224.92, df = 2, p-value < 2.2e-16
#alternative hypothesis: two.sided
#sample estimates:
#  prop 1     prop 2     prop 3 
#0.14333333 0.43478261 0.05566406 
#(+1, higher up q) (0 same) (-1 lower in q)

#comparing those lower down queue (-1) and higher up (+1)
#using all data
prop.test (c(86,57), c(600, 1024), conf.level = 0.95, correct = FALSE)
#2-sample test for equality of proportions without continuity correction
#data:  c(86, 57) out of c(600, 1024)
#X-squared = 36.211, df = 1, p-value = 1.771e-09
#alternative hypothesis: two.sided
#95 percent confidence interval:
# 0.05631097 0.11902757
#sample estimates:
#  prop 1     prop 2 
#0.14333333 0.05566406  
#(+1, higher up q) (-1 lower in q)
#Those lower in queue are sig less likely to approach


#Checking now whether this because we have data from periods after these individuals have fed
#going to look at spread now in terms of carcass condition
#LATE CARCASS <25% REMAINING
#MID CARCASS 74-25% REMAINING
#EARLY CARCASS 75-100% REMAiNING

# LATE	40	423	37	105	41	522
# MID	5	37	29	62	9	268
# EARLY	41	140	20	35	7	190

#LATE carcass (all categories)
prop.test (c(40, 37, 41), c(423, 105, 522), conf.level = 0.95, correct = FALSE)
#3-sample test for equality of proportions without continuity correction
#data:  c(40, 37, 41) out of c(423, 105, 522)
#X-squared = 67.968, df = 2, p-value = 1.741e-15
#alternative hypothesis: two.sided
#sample estimates:
#  prop 1     prop 2     prop 3 
#0.09456265 0.35238095 0.07854406 
#Sig difference, prob driven by same PFQ joining

##LATE carcass (+1 higher up vs -1 lower down)
prop.test (c(40, 41), c(423, 522), conf.level = 0.95, correct = FALSE)
#2-sample test for equality of proportions with continuity correction
#data:  c(40, 41) out of c(423, 522)
#X-squared = 0.5743, df = 1, p-value = 0.4486
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  -0.02231764  0.05435482
#sample estimates:
# prop 1     prop 2 
#0.09456265 0.07854406 
#(+1, higher up q) (-1 lower in q)
#no sig diff, both approach at very low rates


#mid carcass (all categories)
prop.test (c(5,29,9), c(37,62,268), conf.level = 0.95, correct = TRUE)
#3-sample test for equality of proportions without continuity correction
#data:  c(5, 29, 9) out of c(37, 62, 268)
#X-squared = 91.884, df = 2, p-value < 2.2e-16
#alternative hypothesis: two.sided
#sample estimates:
#  prop 1     prop 2     prop 3 
#0.13513514 0.46774194 0.03358209 

#mid carcass (+1 higher up vs -1 lower down)
prop.test (c(5,9), c(37,268), conf.level = 0.95, correct = TRUE)
#2-sample test for equality of proportions with continuity correction
#data:  c(5, 9) out of c(37, 268)
#X-squared = 5.5127, df = 1, p-value = 0.01888
#alternative hypothesis: two.sided
#95 percent confidence interval:
# -0.02607314  0.22917923
#sample estimates:
#  prop 1     prop 2 
#0.13513514 0.03358209 
#(+1, higher up q) (-1 lower in q)
#Higher up queue sig more likely to approach (4x rate)


#EARLY carcass (all categories): EARLY	41	140	20	35	7	190
prop.test (c(41,20,7), c(140,35,190), conf.level = 0.95, correct = TRUE)
#3-sample test for equality of proportions without continuity correction
#data:  c(41, 20, 7) out of c(140, 35, 190)
#X-squared = 72.728, df = 2, p-value < 2.2e-16
#alternative hypothesis: two.sided
#sample estimates:
#  prop 1     prop 2     prop 3 
#0.29285714 0.57142857 0.03684211 

#EARLY carcass (+1 higher up vs -1 lower down)
prop.test (c(41,7), c(140,190), conf.level = 0.95, correct = TRUE)
#2-sample test for equality of proportions with continuity correction
#data:  c(41, 7) out of c(140, 190)
#X-squared = 40.47, df = 1, p-value = 1.997e-10
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  0.1698131 0.3422170
#sample estimates:
# prop 1     prop 2 
#0.29285714 0.03684211 
#(+1, higher up q) (-1 lower in q)
#Higher up queue sig more likely to approach (~8x rate)


#Lower down in queue no more likely to approach in early or late carcass?
#EARLY=7/190, LATE=41/522
prop.test (c(41,7), c(522,190), conf.level = 0.95, correct = TRUE)
  
#late stage comparison of higher up queue versus lower down queue
prop.test (c(41,40), c(522,423), conf.level = 0.95, correct = TRUE)


# script for pairwise comparisons for the significant three way interaction for 
# 500ms lick cluster size data. 

# Tail at postshift 1 

t.test(LCS ~ Contrast, data = TailPairwise500Postshift1) 

#Welch Two Sample t-test

#data:  LCS by Contrast
#t = -4.3336, df = 61.147, p-value = 5.56e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -9.130097 -3.364903
#sample estimates:
#  mean in group SNC mean in group SNCCon 
#13.76812             20.01562 

# Tail at postshift 2 

t.test(LCS ~ Contrast, data = TailPairwise500Postshift2) 

#Welch Two Sample t-test

#data:  LCS by Contrast
#t = -4.4609, df = 61.87, p-value = 3.513e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.062913 -3.072712
#sample estimates:
#  mean in group SNC mean in group SNCCon 
#14.81000             20.37781 

# tunnel at postshift 1

t.test(LCS ~ Contrast, data = TunnelPairwise500Postshift1)

# Welch Two Sample t-test
# 
# data:  LCS by Contrast
# t = -5.3285, df = 50.868, p-value = 2.275e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -10.425311  -4.719064
# sample estimates:
#   mean in group SNC mean in group SNCCon 
# 16.06156             23.63375 



# tunnel at postshift 2 

t.test(LCS ~ Contrast, data = TunnelPairwise500Postshift2)

#Welch Two Sample t-test

#data:  LCS by Contrast
#t = -1.1357, df = 59.965, p-value = 0.2606
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -4.860759  1.340134
#sample estimates:
#  mean in group SNC mean in group SNCCon 
#19.95656             21.71688 
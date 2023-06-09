# script to look at the contrast groups at the preshift phase (last four trials) with their 
# controls. 

# This is for SNC tail handled with SNC Con tail handled animals. Data file is 
# PreshiftTtestsTailSNC.csv

# want to do a unpaired t test between each of the groups to determine whether they had the 
# predicted differences at the preshift phase (ie because they are actually on different concentrations)

t.test(LCS500 ~ Contrast, data = PreshiftTtestsTailSNC)

# Welch Two Sample t-test
# 
# data:  LCS500 by Contrast
# t = 4.1316, df = 40.685, p-value = 0.0001747
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   4.342111 12.649764
# sample estimates:
#   mean in group SNC mean in group SNCCon 
# 25.74719             17.25125 

# This is for SNC tunnel handled with SNC Con tunnel handled animals. Data file is 
# PreshiftTtestsTunnelSNC.csv

t.test(LCS500 ~ Contrast, data = PreshiftTtestsTunnelSNC)

# Welch Two Sample t-test
# 
# data:  LCS500 by Contrast
# t = 0.94872, df = 58.462, p-value = 0.3467
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.710462  7.596087
# sample estimates:
#   mean in group SNC mean in group SNCCon 
# 25.00875             22.56594 

# This is for SPC tail handled with SPC Con tail handled animals. Data file is 
# PreshiftTtestsTailSPC.csv

t.test(LCS500 ~ Contrast, data = PreshiftTtestsTailSPC)

# Welch Two Sample t-test
# 
# data:  LCS500 by Contrast
# t = -1.4491, df = 53.908, p-value = 0.1531
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -7.263276  1.168901
# sample estimates:
#   mean in group SPC mean in group SPCCon 
# 16.93281             19.98000 

#  This is for SPC tunnel handled with SPC Con tunnel handled animals. Data file is 
# PreshiftTtestsTunnelSPC.csv

t.test(LCS500 ~ Contrast, data = PreshiftTtestsTunnelSPC)

# Welch Two Sample t-test
# 
# data:  LCS500 by Contrast
# t = 0.055812, df = 45.959, p-value = 0.9557
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -7.436754  7.860907
# sample estimates:
#   mean in group SPC mean in group SPCCon 
# 26.09645             25.88437 

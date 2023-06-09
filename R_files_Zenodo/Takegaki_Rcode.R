R_code

1) Tukey HSD test

data1<-read.csv("Cannibalism_data.csv")
Comparison of male body size among male status (courtship phase, parental care phase, cannibal males, non-cannibal males)
Code: TukeyHSD(aov(data1$size?data1$status))



2) Steel Dwass test

data1<-read.csv("Cannibalism_data.csv")
Source: source("http://aoki2.si.gunma-u.ac.jp/R/src/Steel-Dwass.R", encoding="euc-jp")
Code: Steel.Dwass(data1$kt, data1$group)



3) t test

Welch Two Sample t-test

data1<-read.csv("Cannibalism_data.csv")
Comparison of the number of eggs among male status (cannibal and non-cannibal)
Code: t.test(eggs ? status, data = data1)



Two Sample t-test

data2<-read.csv("EggDNA_data.csv")
Comparison of the DNA content per egg among male status (cannibal and non-cannibal)
Code: t.test(DNA ? status, data = data2, var=T)




4) ƒÔ2 test

data3<-read.csv("Remating_data.csv")
Comparison of remating success rate between males whose eggs were removed on the spawning day and males whose eggs were removed 4 to 5 days after spawning. 
Code: chisq.test(data3$timing, data3$remating)


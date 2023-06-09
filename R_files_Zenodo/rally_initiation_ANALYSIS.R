#African wild dog rally initiations in Jordan et al.
#Chi-square test for given probabilities testing whether PFQ (position in feeding queue) 
#is related to propensity to initiate pre-hunt rallies in African wild dogs
#46 rally initiations observed by individuals in PFQs 2-7

#See file "Rally initiator data.xls" for raw data
rally <- c(12, 18, 6, 10)
res <- chisq.test(rally, p = c(0.210300429, 0.384120172, 0.184549356, 0.221030043))
res

#Chi-squared test for given probabilities
#data:  rally
#X-squared = 1.2982, df = 3, p-value = 0.7296

#PFQ is unrelated to rally initiation propensity
#Microcosm Analyses
#First to Raid

#Size
Microcosm_Size <-
  matrix(c(6,4,4,6),
         nrow = 2,
         dimnames = list(Size = c("50", "10"),
                         Truth = c("Raided", "Not Raided")))


fet_size<-fisher.test(Microcosm_Size)

#Distance
Microcosm_Dist <-
  matrix(c(8,2,2,8),
         nrow = 2,
         dimnames = list(Distance = c("Near", "Far"),
                         Truth = c("Raided", "Not Raided")))


fet_dist<-fisher.test(Microcosm_Dist)



#All
Microcosm_All <-
  matrix(c(5,5,1,9,3,7,1,9),
         nrow = 4,
         dimnames = list(NestType = c("50 N","50 F", "10 N","10 F"),
                         Truth = c("Raided", "Not Raided")))

fet_all<-fisher.test(Microcosm_All)
#p-value = 0.8183
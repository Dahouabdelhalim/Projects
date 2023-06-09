### CONTAGION vs. AMPLIFICATION vs. NO EFFECT :  |pure_A - pure_B| vs. |mixed_A - mixed_B|

typeDAT<-read.csv('[pathname]/AGETypeSpecific.csv',h=T) # change AGE to GEN1, GEN2, or IW for other experiments
names(typeDAT)

typeDAT_pureX <- subset(typeDAT,Treat=='O') # where X is the type with the higher mean rmsd (W for IW, B for GEN1-2, O for AGE)
typeDAT_pureY <- subset(typeDAT,Treat=='Y') # where Y is the type with the lower mean rmsd (I for IW, A for GEN1-2, Y for AGE)
typeDAT_mixed <- subset(typeDAT,Treat=='YO') # "IW" for IW, "AB, for GEN1-2, "YO" for AGE

pureX <- typeDAT_pureX$RMSD
pureY <- typeDAT_pureY$RMSD
mixedX <- typeDAT_mixed$RMSD[typeDAT_mixed$Type=='O'] # where X is the type with the higher mean rmsd (W for IW, B for GEN1-2, O for AGE)
mixedY <- typeDAT_mixed$RMSD[typeDAT_mixed$Type=='Y'] # where Y is the type with the lower mean rmsd (I for IW, A for GEN1-2, Y for AGE)

diff_pure <- (pureX - pureY) 
diff_mixed <- (mixedX - mixedY)

shapiro.test(diff_pure)
shapiro.test(diff_mixed)

t.test(diff_pure,diff_mixed) 
# wilcox.test(diff_pure,diff_mixed) 

### CONTAGION ASYMMETRY: |Ap-Am| vs. |Bp-Bm|

Xpm<-abs(mixedX-pureX) 
Ypm<-abs(mixedY-pureY)

shapiro.test(Xpm)
shapiro.test(Ypm)

t.test(Xpm,Ypm) 
# wilcox.test(Xpm,Ypm) 
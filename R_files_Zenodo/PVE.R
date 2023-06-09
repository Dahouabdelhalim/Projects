# For PC
setwd("D:/PhD/Crossing/Analyses/Final_Heritability")
getwd()

data = read.csv("PVE.csv")

F2 <- subset(data, Generation=="F2")
North <- subset(data, Generation=="North")
South <- subset(data, Generation=="South")

# # REFLECTANCE # # ----------------------------
# -----------------------------------------
Ref.Plast_N <- (North$Reflectance.Plasticity)
Ref.Plast_S <- (South$Reflectance.Plasticity)
Ref.Plast_F2 <- (F2$Reflectance.Plasticity)

Ref.Plast_vS <- var(Ref.Plast_S, use="complete.obs")
Ref.Plast_vN <- var(Ref.Plast_N, use="complete.obs")
Ref.Plast_vF2 <- var(Ref.Plast_F2, use="complete.obs")

Ref.Cold_N <- (North$Reflectance.Cold)
Ref.Cold_S <- (South$Reflectance.Cold)
Ref.Cold_F2 <- (F2$Reflectance.Cold)

Ref.Cold_vS <- var(Ref.Cold_S, use="complete.obs")
Ref.Cold_vN <- var(Ref.Cold_N, use="complete.obs")
Ref.Cold_vF2 <- var(Ref.Cold_F2, use="complete.obs")

Ref.Warm_N <- (North$Reflectance.Warm)
Ref.Warm_S <- (South$Reflectance.Warm)
Ref.Warm_F2 <- (F2$Reflectance.Warm)

Ref.Warm_vS <- var(Ref.Warm_S, use="complete.obs")
Ref.Warm_vN <- var(Ref.Warm_N, use="complete.obs")
Ref.Warm_vF2 <- var(Ref.Warm_F2, use="complete.obs")

### Heritability
Ref.Plast_h2 <- (Ref.Plast_vF2-(sqrt(Ref.Plast_vN * Ref.Plast_vS)))/Ref.Plast_vF2
print(Ref.Plast_h2)

Ref.Cold_h2 <- (Ref.Cold_vF2-(sqrt(Ref.Cold_vN * Ref.Cold_vS)))/Ref.Cold_vF2
print(Ref.Cold_h2)

Ref.Warm_h2 <- (Ref.Warm_vF2-(sqrt(Ref.Warm_vN * Ref.Warm_vS)))/Ref.Warm_vF2
print(Ref.Warm_h2)

# # FLOWERING TIME # # ----------------------------
# -----------------------------------------
FT.Plast_N <- (North$FT.Plasticity)
FT.Plast_S <- (South$FT.Plasticity)
FT.Plast_F2 <- (F2$FT.Plasticity)

FT.Plast_vS <- var(FT.Plast_S, use="complete.obs")
FT.Plast_vN <- var(FT.Plast_N, use="complete.obs")
FT.Plast_vF2 <- var(FT.Plast_F2, use="complete.obs")

FT.Cold_N <- (North$FT.Cold)
FT.Cold_S <- (South$FT.Cold)
FT.Cold_F2 <- (F2$FT.Cold)

FT.Cold_vS <- var(FT.Cold_S, use="complete.obs")
FT.Cold_vN <- var(FT.Cold_N, use="complete.obs")
FT.Cold_vF2 <- var(FT.Cold_F2, use="complete.obs")

FT.Warm_N <- (North$FT.Warm)
FT.Warm_S <- (South$FT.Warm)
FT.Warm_F2 <- (F2$FT.Warm)

FT.Warm_vS <- var(FT.Warm_S, use="complete.obs")
FT.Warm_vN <- var(FT.Warm_N, use="complete.obs")
FT.Warm_vF2 <- var(FT.Warm_F2, use="complete.obs")

### Heritability
FT.Plast_h2 <- (FT.Plast_vF2-(sqrt(FT.Plast_vN * FT.Plast_vS)))/FT.Plast_vF2
print(FT.Plast_h2)

FT.Cold_h2 <- (FT.Cold_vF2-(sqrt(FT.Cold_vN * FT.Cold_vS)))/FT.Cold_vF2
print(FT.Cold_h2)

FT.Warm_h2 <- (FT.Warm_vF2-(sqrt(FT.Warm_vN * FT.Warm_vS)))/FT.Warm_vF2
print(FT.Warm_h2)


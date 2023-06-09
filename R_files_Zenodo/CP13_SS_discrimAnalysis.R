# CP13_SS_discrimAnalysis.R

#### Version 5: 2 November, 2016
#
# Reflects combining Leptospermum and Baeckea pollen counts.
#
# Reflects replacing the Astelia:Lomatia ratio with each taxon individually
# - PEH
#
####

# Clear and set work space ------------------------------------------------
rm(list=ls())
cat("\\14")
graphics.off()
library(MASS)
library(lattice)
setwd('L:\\\\1_projectsData\\\\Tasmania_Australia\\\\Analysis\\\\surf_samp\\\\charcoal')

# Load data ---------------------------------------------------------------
data = read.csv('TAS_SS_16_11_11.csv', header = T);   # pollen counts

burnCode = data$burn_index;   # Burn code: 0 = unburned, 1 = mixed, 2 = burned.
polSum = data$TOTAL_terrestrial  # Pollen sum
siteID = data$site_ID  # Site name
charCon = data$total_char_con_pieces_cm.3 # [#/cm3] charcoal concentration 


# Calculate pollen percentages --------------------------------------------
polIn = c(7:56)        # Index for pollen counts, to use in sum.
polData = data[,polIn]    # [counts] Raw pollen counts
taxa = colnames(polData)  # Taxa names for polData

#### COMBINE LEPTOSPERMUM AND BAECKEA
idx1 = which(taxa == 'Leptospermum')
idx2 = which(taxa == 'Baeckea')
polData[,idx1] = polData[,idx1] + polData[,idx2] # Lep. + Baeckea
polData[,idx2] = 0
taxa[idx1] = 'Leptospermum-Baeckea'
####

dimPolData = dim(polData) # Sample sizes.
nTaxa = dimPolData[2]     # Number of taxa
nSamples = dimPolData[1]  # Number of samples

polPer = matrix(NA,nSamples,nTaxa);    # Space for pollen percentages
for (i in 1:nSamples){
  polPer[i,] = t(polData[i,] / polSum[i])    # [%] Percent of pollen sum
}

# Modify taxa names
idx = which(taxa == 'Astelia_alpina')
taxa[idx] = 'Astelia alpina'


# Calculate pollen ratios -------------------------------------------------
    # Ratio suggestions: [a - b] / [a + b]
    #
    # x. a : b  
    # 1. Athrotaxis : Poaceae
    # 2. Athrotaxis : Poaceae+Juncus
    # 3. Athrotaxis+Astelia : Poaceae+Juncus  (I think this might be a good)
    # 4. Astelia : Lomatia
    # 5. Athrotaxis : Baeckea


ratios = c('Athrotaxis:Poaceae','Ath:Poa+Jun','Ath+Astelia:Poa+Jun','Astelia:Lomatia',
           'Ath:Baeckea')

polRatios = matrix(NA,nSamples,length(ratios));    # Space for pollen percentages

# Calculate Ratios:
Ath.in = which(taxa == "Athrotaxis.Diselma.Type")
Poa.in = which(taxa == "Poaceae")
Jun.in = which(taxa == "Juncus")
Astelia.in = which(taxa == "Astelia alpina")
Lomatia.in = which(taxa == "Lomatia")
Baeckea.in = which(taxa == "Baeckea")


polRatios[,1] = (polData[,Ath.in] - polData[,Poa.in]) / 
                (polData[,Ath.in] + polData[,Poa.in])
polRatios[,2] = (polData[,Ath.in] - (polData[,Poa.in] + polData[,Jun.in])) /
                (polData[,Ath.in] + (polData[,Poa.in] + polData[,Jun.in]))
polRatios[,3] = ((polData[,Ath.in] + polData[,Astelia.in]) - (polData[,Poa.in] + polData[,Jun.in])) /
                ((polData[,Ath.in] + polData[,Astelia.in]) + (polData[,Poa.in] + polData[,Jun.in]))
polRatios[,4] = (polData[,Astelia.in] - polData[,Lomatia.in]) / 
                (polData[,Astelia.in] + polData[,Lomatia.in])                                 
polRatios[,5] = (polData[,Ath.in] - polData[,Baeckea.in]) / 
                (polData[,Ath.in] + polData[,Baeckea.in])                                 

  #### NOTE: Ratio #1 and #4 are really the only two that are not highly correlated
# 
#     win.graph()
#     plot(polRatios[,1],polRatios[,4],pch = 21, bg = 'white',
#          xlab = ratios[1],
#          ylab = ratios[4])
#     #xlim = c(-5,5),
#     #ylim = c(-3.2,3.2))
#     
#     idx = which(burnCode == 1)
#     points(polRatios[idx,1],polRatios[idx,4],pch = 21, bg = 'lightgrey')
#     
#     idx = which(burnCode == 2)
#     points(polRatios[idx,1],polRatios[idx,4],pch = 21, bg = 'black')
#     
#     legend(0.5,0.3,c('Unburned','Mixed','Burned'),
#            xjust = 0,
#            pch = 21,
#            x.intersp = 1,
#            cex = 0.8,
#            pt.bg = c('white','lightgrey','black'))


# Add pollen ratios to taxa --------------------------------------
taxa2 = matrix(NA,1,length(taxa) + 1)    # Space for pollen percentages
taxa2[1,1:length(taxa)] = taxa
taxa2[1,seq(length(taxa)+1,length(taxa2))] = ratios[c(1)]

taxa = taxa2

# Plot data ---------------------------------------------------------------

### (1) All a priori relevant taxa + independent ratios:
#taxaPlot = c(1, 3, 4, 6, 11, 13, 15, 20, 21, 22, 26, 31, 32, 34, 39, 40)

### (2) All a priori relevant taxa, minus those in idenpendent ratios + independent ratios
#taxaPlot = c(1, 4, 11, 13, 15, 20, 21, 22, 32, 34, 39, 40)

### (3) Same as #2, but without taxa that did not differ among burn classes.
#taxaPlot = c(13, 15, 22, 34, 51, 52)

### (4) Same as #3, but taking Astelia and Lomatia out of the ratio and displaying individually
taxaPlot = c(3, 13, 15, 22, 31, 34, 51)

polDA = cbind(sqrt(0.5+100*polPer[,taxaPlot[1:(length(taxaPlot)-1)]]), polRatios[,1])  # Square-root transformed pollen data
#polDA =  polRatios # Square-root transformed pollen data
#taxaPlot = c(1:5)

#### Add charcoal concentration to pollen, for plotting and DA
polDA = cbind(polDA,log(charCon+1))  # Add charcoal concentrations

nPlots = length(taxaPlot)+1
pVal_taxa = matrix(NA,length(taxaPlot)+1,1)

## Kruskal Wallis rank sum test, for each taxon:
for (i in 1:nPlots){
  kw.results = kruskal.test(polDA[,i] ~ burnCode)
  pVal_taxa[i] = kw.results$p.value
}
sortResults = sort(pVal_taxa, index.return = T)
plotIn = sortResults$ix

win.graph()
par(mfrow=c(2,4),omi=c(0,0,0,0),mar = c(3,4,4,2))  
for (i in 1:nPlots){    
  boxplot(polDA[,plotIn[i]] ~ burnCode, col = c("white","lightgray","darkgray"),
          names=c("Unburned","Mixed","Burned"),
          notch = T,
          xlim = c(0.5,3.5),
          main = ' ',
          axes = F)

  axis(2)
  if(i < 2){
    title(ylab = 'Pollen ratio', font.lab = 2)
  }
  if((i > 1) & (i < nPlots)){
    title(ylab = 'Sqrt.-trans. pollen %', font.lab = 2)
  }
  if(i == nPlots){
    title(ylab = expression(paste("Log-trans. charcoal concentration (#/cm" ^"3", " )")), font.lab = 2)
  }
    axis(side = 1, at = c(1,2,3), label = c('Unburned','Mixed','Burned'))  
  if(i < nPlots){
    title(main = paste(taxa[taxaPlot[plotIn[i]]]),font.main = 2,cex.main = 0.9)
  }  else {
    title(main = paste('Macroscopic charcoal'),font.main = 2,cex.main = 0.9)
  }
}


# Linear discriminant analysis --------------------------------------------

#### Redefine "polDA" to exclude charcoal concentrations, for DA
polDA = cbind(sqrt(0.5+100*polPer[,taxaPlot[1:(length(taxaPlot)-1)]]), polRatios[,1])  # Square-root transformed pollen data
#polDA =  polRatios # Square-root transformed pollen data

## This index does nothing if charcoal samples are not used in DA, but code below utilized idx1.
idx1 = which(polDA[,1] > -999) 
burnCode.da = burnCode[idx1]

####

fit = lda(polDA[idx1,], burnCode.da, method = 'moment')
fit.cv = lda(polDA[idx1,], burnCode.da, CV = T)
da.values = predict(fit)

win.graph()
par(mfrow=c(1,2))
plot(da.values$x[,1],da.values$x[,2],pch = 21, bg = 'white',
     xlab = 'Discriminant function 1 (78%)',
     ylab = 'Discriminant function 2 (21%)')
     #xlim = c(-5,5),
     #ylim = c(-3.2,3.2))

idx = which(burnCode.da == 1)
points(da.values$x[idx,1],da.values$x[idx,2],pch = 21, bg = 'lightgrey')

idx = which(burnCode.da == 2)
points(da.values$x[idx,1],da.values$x[idx,2],pch = 21, bg = 'black')

legend(2,-2,c('Unburned','Mixed','Burned'),
       xjust = 0,
       pch = 21,
       x.intersp = 1,
       cex = 0.8,
       pt.bg = c('white','lightgrey','black'))

plot(fit$scaling[,1],fit$scaling[,2],pch = 3,
     xlab = 'Discriminant function 1 (77%)',
     ylab = 'Discriminant function 2 (22%)')
#xlim = c(-3,1),
#ylim = c(-2.5,0.75))

for (i in seq(1:length(taxaPlot)+1)){
  text(fit$scaling[i,1],fit$scaling[i,2],taxa[taxaPlot[i]],
       cex = 0.5, pos = 4)
}

#### Calculate classification rates for full DA

#TPR = matrix(NA,1,3) 
#FPR = matrix(NA,1,3)
acc = matrix(NA,1,3)
acc.xval = matrix(NA,1,3)

for (i in c(0,1,2)){
  # Accuracy, raw
    idx = which(da.values$class == i)
    P = length(idx)
    TP =  length(which(burnCode.da[idx] == i))  # 
    
    idx = which(da.values$class != i)
    N = length(idx)
    TN = length(which(burnCode.da[idx] != i))
    
    acc[i+1] = (TP + TN) / (P+N) 

    # Accuracy, cross-validation
    idx = which(fit.cv$class == i)
    P = length(idx)
    TP =  length(which(burnCode.da[idx] == i))  # 
    
    idx = which(fit.cv$class != i)
    N = length(idx)
    TN = length(which(burnCode.da[idx] != i))
    
    acc.xval[i+1] = (TP + TN) / (P+N) 
}

####
# plot(fit.cv$class,burnCode+1,xlab = 'Predicted',ylab = 'Observed',
#      axes = F,
#      main = 'Cross validation results')
# axis(side = 1, at = c(1,2,3), label = c('Unburned','Mixed','Burned'))
# axis(side = 2, at = c(1,2,3), label = c('Unburned','Mixed','Burned'))

# Apply linear discriminant function to unknown samples ------------------------

setwd('L:\\\\1_projectsData\\\\Tasmania_Australia\\\\Analysis\\\\1_Sites\\\\CP13_10\\\\pollen')
testData = read.csv('CP13_10A_pollen_16_11_10.csv', header = T);  # Test pollen data
setwd('L:\\\\1_projectsData\\\\Tasmania_Australia\\\\Analysis\\\\surf_samp')

testPolData = testData[,polIn]    # [counts] Raw pollen counts
testTaxa = colnames(testPolData)  # Taxa names for polData

dimTestPolData = dim(testPolData) # Sample sizes.
nTestTaxa = dimTestPolData[2]     # Number of taxa
nTestSamples = dimTestPolData[1]  # Number of samples

testPolSum = testData$TOTAL_terrestrial  # Pollen sum

testPolPer = matrix(NA,nTestSamples,nTestTaxa);    # Space for pollen percentages
for (i in 1:nTestSamples){
  testPolPer[i,] = t(testPolData[i,] / testPolSum[i])    # [%] Percent of pollen sum
}

# Add on pollen ratios

testPolRatio = matrix(NA,nTestSamples,1)

testPolRatio[,1] = (testPolData[,Ath.in] - testPolData[,Poa.in]) / 
  (testPolData[,Ath.in] + testPolData[,Poa.in])
# testPolRatio[,2] = (testPolData[,Astelia.in] - testPolData[,Lomatia.in]) / 
#  (testPolData[,Astelia.in] + testPolData[,Lomatia.in])         

testPolDA = cbind(sqrt(0.5+100*testPolPer[,taxaPlot[1:(length(taxaPlot)-1)]]),testPolRatio)  # Square-root transformed pollen data

# Plot test data
win.graph()
par(mfrow=c(2,4),omi=c(0,0,0,0))
for (i in 1:nPlots){
  boxplot(testPolDA[,plotIn[i]],
          notch = T,
          xlim = c(0.5,1.5),
          main = ' ',
          axes = F)
  axis(2)
  if(i == 9){
    title(ylab = 'Square-root transformed pollen %', font.lab = 2)
  }
  title(main = paste(taxa[taxaPlot[plotIn[i]]]), font.main = 1, cex.main = 0.9)
}


## Apply DA function: plot for full fossil record

predResults = predict(fit,testPolDA)

## Plot DA training and test datasets:
# win.graph()
# plot(da.values$x[,1],da.values$x[,2],pch = 21, bg = 'white',
#      xlab = 'Discriminant function 1 (77%)',
#      ylab = 'Discriminant function 2 (23%)',
#      xlim = c(-5,6.5),
#      ylim = c(-4.5,3))
# idx = which(burnCode == 1)
# points(da.values$x[idx,1],da.values$x[idx,2],pch = 21, bg = 'grey')
# 
# idx = which(burnCode == 2)
# points(da.values$x[idx,1],da.values$x[idx,2],pch = 21, bg = 'black')
# 
# points(predResults$x[,1],predResults$x[,2],pch = 4, col = 28)
# lines(predResults$x[,1],predResults$x[,2],pch = 4, col = 28)
# 
# for (i in 1:nTestSamples){
#   text(predResults$x[,1],predResults$x[,2],testData[,3]/1000)
# }
# 
# legend(4,-3,c('Unburned','Mixed','Burned','Unknown'),
#        xjust = 0,
#        pch = c(21,21,21,4),
#        col = c(24,24,24,28),
#        x.intersp = 1,
#        cex = 0.8,
#        pt.bg = c('white','lightgrey','black'))

#################################################
## Apply DA function: plot for 5 fossil samples

# sampIn = c(3, 5, 7, 10, 11)  # Index for samples used in SS paper. 
sampIn = c(2,5,8,10,13)  # Index for samples used in SS paper. 
predResults = predict(fit,testPolDA[sampIn,])

# Plot DA training and test datasets:
win.graph(9,7)

#### Plot all sites
plot(da.values$x[,1],da.values$x[,2],pch = 21, bg = 'darkgreen',
     xlab = 'Discriminant function 1 (79%)',
     ylab = 'Discriminant function 2 (21%)',
     xlim = c(-6.5,4.5), # c(-5,5),
     ylim = c(-3.8,3.8), # c(-4,2.5),
     cex = 2,
     axis = F,
     main = ' ')

#### Plot mixed sites
idx = which(burnCode.da == 1)
points(da.values$x[idx,1],da.values$x[idx,2],pch = 21, bg = 'grey',
       cex = 2)

#### Plot burned sites
idx = which(burnCode.da == 2)
points(da.values$x[idx,1],da.values$x[idx,2],pch = 21, bg = 'red',
       cex = 2)

#### Plot misclassified samples
text(-1.35,  1.6,'SS21',cex = 0.75)
text(1.4, -0.2444870,'SS2',cex = 0.75)

#### Plot species locations
offSet.x = c(-0.5,-0.8,-1,0,0.2,-1.5,-0.4)
offSet.y = c(0.3,0.2,0.20,-0.15,0.2,0.3,-0.3)

for (i in seq(1:7)){
  text(fit$scaling[i,1]+offSet.x[i],
       fit$scaling[i,2]+offSet.y[i],
       taxa[taxaPlot[i]],
       cex = 0.8, pos = 4)
  arrows(0,0,fit$scaling[i,1],fit$scaling[i,2],
         length = 0.1)
}
# text(fit$scaling[i,1]+offSet.x[i],
#      fit$scaling[i,2]+offSet.y[i],
#      cex = 0.8, pos = 4)

#### Plot fossil samples
points(predResults$x[,1],predResults$x[,2],pch = 4, cex = 2)

legend(-6.5,3.8,c('Unburned','Mixed','Burned','Unknown'),
       xjust = 0,
       pch = c(21,21,21,4),
       col = c(24,24,24,1),
       x.intersp = 1,
       cex = 1,
       pt.bg = c('darkgreen','grey','red'))


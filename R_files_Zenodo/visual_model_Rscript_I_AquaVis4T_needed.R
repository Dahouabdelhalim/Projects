######################################
#
#     AquaVis4T.R set of functions - contrast using QC only
#       Tetrachromat
#
#####################################
# Karen Carleton July 2021          
#
# This set of functions are useful for modeling in an aquatic enviroment. It includes
# funcitons for calculating visual pigments, absorptance, quantum catch and JND.
# It will then calculate JNDs between a pair of targets using a visual system with 
# UV, Short, Medium, and Long wavelength sensitive cones using 
# the Receptor Noise Limited model. 
# Refs: Vorobyev and Osorio 1998 for RNL  
# Champ et al 2016 and Olsson et al 2015 for the absolute calculations
# The functions contained in this version include: 
#   This includes 
#     Opsin_gov to calculate visual pigment templates 
#     Absorptance_calc to determine absorptance
#     Step_data calculates a single peak or single step color
#     Qcatch_calc to calculate quantum catch of U, S, M and L cones for a set of colors
#     Qcatchlight_calc calculates quantum catch of an illuminant (without reflecting off surface) such as sidewelling irradiance
#     QC2illum_calc calculates quantum catch but uses 2 illuminants: 1 for objects illuminated 
#            e.g. downwelling illum of rocks or algae and 1 for von Kries normalization
#     JND4one_calc to calculate JND from just one pair of colors for tetrachromat
#     JND4_calc to calculate JNDs for a set of colors for tetrachromat
# This assumes a tetrachromat with 2 single and 2 double cones.  These can coexpress visual pigments.
#
#================================================================================
#   FUNCTION: Opsin_gov 
#      Calculate visual pigments absorption using Govardovskii et al 2000 templates
#           lbeg, lend = wavelength range
#           lmax = peak wavelength
#           A1_chrom = % of A1 chromophore
#--------------------------------------------------------------------------------
Opsin_gov<-function(lbeg, lend, lmax, A1_chrom)  {
  # fit A1 peak
  a1<-0.8795+0.0459*exp(-(lmax-300)^2/11940)
  lmb1<-189+0.315*lmax
  p1<--40.5+0.195*lmax  #Govardovskii calls this b but this confuses it with the alpha peak equation
  i<-lbeg:lend
  x<-lmax/i
  Salpha_A1<-1/((exp(69.7*(a1-x))+exp(28*(0.922-x))+exp(-14.9*(1.104-x))+0.674))
  Sbeta_A1<-0.26*exp(-((i-lmb1)/p1)^2)
  Speak_A1<-Salpha_A1+Sbeta_A1
  # fit A2 peak
  a2<-0.875+0.0268*(exp((lmax-665)/40.7))  
  A2<- 62.7 + 1.834*exp((lmax-625)/54.2)
  lmb2<-216.7+0.287*lmax
  p2<-317-1.149*lmax+0.00124*lmax*lmax  #Govardovskii calls this b but this confuses it with the alpha peak equation
  Salpha_A2<-1/((exp(A2*(a2-x))+exp(20.85*(0.9101-x))+exp(-10.37*(1.1123-x))+0.5343))
  Sbeta_A2<-0.26*exp(-((i-lmb1)/p1)^2)
  Speak_A2<-Salpha_A2+Sbeta_A2
  # weight by chromophore, sum and normalize
  Speaktot<-A1_chrom*Speak_A1+(100-A1_chrom)*Speak_A2
  Snorm<-Speaktot/max(Speaktot)
  Gov_results<-data.frame(i,Snorm)
  return(Gov_results)}
#================================================================================
#   FUNCTION: Absorptance_calc 
#      Calculate visual pigments absorptance which is more exact and doesn't require
#         small absorption assumption
#           Vispig = matrix of visual pigment absorbances (U, S, M and L)
#			      Abs_peak = peak vis pig absorbances
#           Lengths = photoreceptor lengths (l)
#--------------------------------------------------------------------------------
Absorptance_calc<-function(Vispig, Abs_peak, Length) {
  colnames(Vispig)<-c("Wavelength","VPu","VPs","VPm","VPl")
  AbsorptanceU<-1-exp(-Vispig$VPu*Abs_peak[1]*Length[1])
  AbsorptanceS<-1-exp(-Vispig$VPs*Abs_peak[1]*Length[1])
  AbsorptanceM<-1-exp(-Vispig$VPm*Abs_peak[2]*Length[2])
  AbsorptanceL<-1-exp(-Vispig$VPl*Abs_peak[3]*Length[3])
  Absorptance_results<-data.frame(Vispig$Wavelength, AbsorptanceU, AbsorptanceS, AbsorptanceM, AbsorptanceL)
  colnames(Absorptance_results)<-c("Wavelength", "VPu", "VPs", "VPm", "VPl")
  return(Absorptance_results) }
#
#==============================================================================
#  FUNCTION: Step_data
#     Calculate ONE Gaussian stepped color (variable peak and rise rate)
#          lbeg, lend = wavelength range
#          peak = first stepped color (e.g. 400 nm)
#          peaksd = standard deviation of peak
#          peakoffset = offset from baseline
# ------------------------------------------------------------------------------
Step_data<-function(lbeg, lend, peak, peaksd, peakoffset) {
  wavetotal<-lend-lbeg+1
  waves<-matrix(lbeg:lend)
  peak_data<-matrix(nrow=wavetotal)
  for (i in 1:wavetotal) {
    coeff<-sqrt(2*pi)*peaksd*(1-peakoffset)+peakoffset
    peak_data[i]<-coeff*dnorm(waves[i], peak, peaksd)+peakoffset
    if (waves[i]>peak)  {
      peak_data[i]<-1
    }
  }   # i loop
  rownames(peak_data)<-waves
  return(peak_data)   }
#==============================================================================
#   FUNCTION:  Qcatch_calc  
#     Calculate the stimulation of the U, S, M and L cones in terms of quantum catch
#         light = the illuminant
#         fishlens = the lens transmission
#         Vispig = 3 visual pigments which can be coexpresing and involve A1 or A2
#         targets = the set of color targets
#     This version works on target file where there is no wavelength column
# ----------------------------------------------------------------------
Qcatch_calc<-function(light, fishlens, Vispig, targets) {
  target_labels<-colnames(targets)
  number_colors<-length(target_labels)   # number of colors in set
  #
  # Calculate the von Kries correction for accommodation to the illuminant
  colnames(Vispig)<-c("Wavelength","VPu","VPs","VPm","VPl")
  vonK_illum_QCu=0
  vonK_illum_QCs=0
  vonK_illum_QCm=0
  vonK_illum_QCl=0
  vonK_illum_QCu<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPu)
  vonK_illum_QCs<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPs)
  vonK_illum_QCm<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPm)
  vonK_illum_QCl<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPl)
  qcatchu=NULL
  qcatchs=NULL
  qcatchm=NULL
  qcatchl=NULL
  for (j in 1:number_colors) {
    qcatchu[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPu)/vonK_illum_QCu
    qcatchs[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPs)/vonK_illum_QCs
    qcatchm[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPm)/vonK_illum_QCm
    qcatchl[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPl)/vonK_illum_QCl
  }
  Qcatchdata_results<-data.frame(qcatchu,qcatchs,qcatchm,qcatchl)
  return(Qcatchdata_results)  }
#==============================================================================
#   FUNCTION:  Qcatchlight_calc  
#     Calculate the stimulation of the U, S, M and L cones in terms of quantum catch
#         light = the illuminant, side welling IRRadiance
#         fishlens = the lens transmission
#         Vispig = 3 visual pigments which can be coexpresing and involve A1 or A2
#         targets = the set of side welling irradiances that are themselves the "targets"
#     This version works on target file where there is no wavelength column
# ----------------------------------------------------------------------
Qcatchlight_calc<-function(light, fishlens, Vispig, targets) {
  target_labels<-colnames(targets)
  number_colors<-length(target_labels)   # number of colors in set
  #
  # Calculate the von Kries correction for accommodation to the illuminant
  colnames(Vispig)<-c("Wavelength","VPu","VPs","VPm","VPl")
  vonK_illum_QCu=0
  vonK_illum_QCs=0
  vonK_illum_QCm=0
  vonK_illum_QCl=0
  vonK_illum_QCu<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPu)
  vonK_illum_QCs<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPs)
  vonK_illum_QCm<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPm)
  vonK_illum_QCl<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPl)
  qcatchu=NULL
  qcatchs=NULL
  qcatchm=NULL
  qcatchl=NULL
  for (j in 1:number_colors) {   # the side welling irradiance is the target so dont multiply by illum
    qcatchu[j]<-sum(targets[,j]*fishlens$Transmission*Vispig$VPu)/vonK_illum_QCu
    qcatchs[j]<-sum(targets[,j]*fishlens$Transmission*Vispig$VPs)/vonK_illum_QCs
    qcatchm[j]<-sum(targets[,j]*fishlens$Transmission*Vispig$VPm)/vonK_illum_QCm
    qcatchl[j]<-sum(targets[,j]*fishlens$Transmission*Vispig$VPl)/vonK_illum_QCl
  }
  Qcatchlightdata_results<-data.frame(qcatchu,qcatchs,qcatchm,qcatchl)
  return(Qcatchlightdata_results)  }
#
#==============================================================================
#  FUNCTION : QC2illum_calc  
#     Calculate quantum catch 
#     This has 2 illuminants where one can illum object and other is for von Kries and scatter
#     This needs:
#         IRRlight1 = the illuminant, e.g. down welling IRRadiance for substrate
#         IRRlight2 = the illuminant that is used for von Kries
#         RADlight = the sidewelling RADiance
#         fishlens = the lens transmission
#         Vispig = 4 visual pigments which can be coexpresing and involve A1 or A2
#         targets = the set of side welling irradiance targets
#         atten = attenuation prop of water
#         dist = distance of colored target from viewer
# ----------------------------------------------------------------------
QC2illum_calc<-function(IRRlight1, IRRlight2, fishlens, Vispig, targets) {
  target_labels<-colnames(targets)
  number_colors<-length(target_labels)   # number of colors in set
  #
  # Calculate the von Kries correction for accommodation to illuminant 2
  colnames(Vispig)<-c("Wavelength","VPu","VPs","VPm","VPl")
  vonK_illum_QCu=0
  vonK_illum_QCs=0
  vonK_illum_QCm=0
  vonK_illum_QCl=0
  vonK_illum_QCu<-sum(IRRlight2$IRRAD*fishlens$Transmission*Vispig$VPu)
  vonK_illum_QCs<-sum(IRRlight2$IRRAD*fishlens$Transmission*Vispig$VPs)
  vonK_illum_QCm<-sum(IRRlight2$IRRAD*fishlens$Transmission*Vispig$VPm)
  vonK_illum_QCl<-sum(IRRlight2$IRRAD*fishlens$Transmission*Vispig$VPl)
  qcdistu=NULL
  qcdists=NULL
  qcdistm=NULL
  qcdistl=NULL
#  QC using illum1 but von Kries from illum2
  for (j in 1:number_colors) {
    qcdistu[j]<-sum((IRRlight1$IRRAD*targets[,j])*fishlens$Transmission*Vispig$VPu)/vonK_illum_QCu
    qcdists[j]<-sum((IRRlight1$IRRAD*targets[,j])*fishlens$Transmission*Vispig$VPs)/vonK_illum_QCs
    qcdistm[j]<-sum((IRRlight1$IRRAD*targets[,j])*fishlens$Transmission*Vispig$VPm)/vonK_illum_QCm
    qcdistl[j]<-sum((IRRlight1$IRRAD*targets[,j])*fishlens$Transmission*Vispig$VPl)/vonK_illum_QCl
  }
  QCdist_results<-data.frame(qcdistu,qcdists,qcdistm,qcdistl)
  colnames(QCdist_results)<-c("qcatchu", "qcatchs", "qcatchm", "qcatchl")
  return(QCdist_results)  
}
#==============================================================================
#  FUNCTION : JND4one_calc   
#     Calculate ONE JND value between two targets
#
# ----------------------------------------------------------------------
JND4one_calc<-function(QC1, QC2, Webers)  {
  wu<-Webers[1]
  ws<-Webers[2]
  wm<-Webers[3]
  wl<-Webers[4]
  denom<-(wu*ws*wm)^2+(wu*ws*wl)^2+(wu*wm*wl)^2+(ws*wm*wl)^2
  Dfu<-log(QC1[1]/QC2[1])        # note log calculates natural log
  Dfs<-log(QC1[2]/QC2[2])
  Dfm<-log(QC1[3]/QC2[3])
  Dfl<-log(QC1[4]/QC2[4])
  JND4one_result<-sqrt(((wu*ws)^2*(Dfm-Dfl)^2+(wu*wm)^2*(Dfs-Dfl)^2+(wu*wl)^2*(Dfs-Dfm)^2+(ws*wm)^2*(Dfu-Dfl)^2+(ws*wl)^2*(Dfu-Dfm)^2+(wm*wl)^2*(Dfu-Dfs)^2)/denom)
  return(JND4one_result)  }
#==============================================================================
#  FUNCTION : JND4_calc   
#     Calculate JND values between two sets of targets
#
# ----------------------------------------------------------------------
JND4_calc<-function(QCset1, colorlabel1, QCset2, colorlabel2, Webers)  {
  number_colors1<-length(colorlabel1)
  number_colors2<-length(colorlabel2)
  JND4_results=data.frame(matrix(nrow=number_colors1, ncol=number_colors2))
  wu<-Webers[1]
  ws<-Webers[2]
  wm<-Webers[3]
  wl<-Webers[4]
  denom<-(wu*ws*wm)^2+(wu*ws*wl)^2+(wu*wm*wl)^2+(ws*wm*wl)^2
  for (i in 1:number_colors1) {
    for (j in 1:number_colors2) {
      Dfu<-log(QCset1[i,1]/QCset2[j,1])
      Dfs<-log(QCset1[i,2]/QCset2[j,2])
      Dfm<-log(QCset1[i,3]/QCset2[j,3])
      Dfl<-log(QCset1[i,4]/QCset2[j,4])
      JND4_results[i,j]<-sqrt(((wu*ws)^2*(Dfm-Dfl)^2+(wu*wm)^2*(Dfs-Dfl)^2+(wu*wl)^2*(Dfs-Dfm)^2+(ws*wm)^2*(Dfu-Dfl)^2+(ws*wl)^2*(Dfu-Dfm)^2+(wm*wl)^2*(Dfu-Dfs)^2)/denom)
      } } 
  rownames(JND4_results)<-colorlabel1
  colnames(JND4_results)<-colorlabel2
  return(JND4_results)  }
#
#==============================================================================
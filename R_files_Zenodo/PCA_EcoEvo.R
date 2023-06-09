#Ecology and Evolution Final Pretty Code
#By the awesome Chloe Fouilloux

#First let's load some packages
library(vegan)
library(dplyr)
library(tidyr)
library(MASS)
library(sjPlot)
library(ggplot2)
library(AICcmodavg)

#### DATA----
#Let's load in that sweet, sweet 2019 data.
frog1<- read.csv("frogpool19_ECOEVO.csv", 
                 header = T,fileEncoding="UTF-8-BOM") #fromg Andrius

measure<- read.csv("FP_measures.csv", header = T,fileEncoding="UTF-8-BOM", sep = ",")
measure<- measure %>%
  mutate(pool_id = Pool_ID) %>%
  mutate(s_area = (Length/2*Width/2)*pi, #divide length and width by 2 for radii
         sa_depth = s_area/Depth) %>%
  dplyr::select(sa_depth, pool_id)

frog2<- dplyr::full_join(frog1, measure, by = "pool_id") #now has sa_depth ratio

################
#Tidy Data----
###############
#PCA must be numeric so let's just subset the numeric variables.

numfrog<- frog2 %>%
  mutate(total_other = oophagus_tp + tp_femo_hahneli,
         other_y_n = ifelse (oophagus_tp > 0, 1, 
                             ifelse(tp_femo_hahneli > 0, 1, 0)),
         tinc_Y_N = ifelse( tadpoles_tinc > 0, 1, 0),
         species = ifelse (tp_femo_hahneli> 0, "Femoralis", 
                           ifelse(oophagus_tp >0, "Oophagus", "Tinctorius")),
         pred_size_ave = predator_size_sum/predator_count,
         species = as.factor(species)) 

numfrog$pred_size_ave[is.nan(numfrog$pred_size_ave)]<-0

numfrog <- numfrog %>%
  tidyr::drop_na()

numfrog1<- numfrog %>%
  dplyr::select(- pool_type, - date, - pool_id, - tree,
                -pH_strip, -tadpoles_other, -species, -other_y_n, #remove cat. factors 
                -presence_by_size, -oophagus_tp, -tp_femo_hahneli, -predator_size_sum)  #remove double counts w/Janne, keep sa_depth

#remove tinctorius data
random_frog<- numfrog1 %>%
  dplyr::select(-tinc_Y_N, -tadpoles_tinc) %>%#remove response variable here
  rename("Height" = height,
         "Water_Cap" = water_capacity,
         "Leaf_litter"= leaflit_vol_ml,
         "Hardness" = hardness,
         "NO3" = no3,
         "Salinity" = salinity,
         "Invert_Div"= invert_diversity,
         "Invert_Dens" = dens_invert_total,
         "Pred_count" = predator_count,
         "Amphib_Div"= amphib_diversity,
         "Total_Other" = total_other) 


######## ######## ######## ########
# PCA MODELS ----
######## ######## ######## ########
#first we have to get all this data together
#remember our observation based PCA obs_PCA
dataset_obj = random_frog
obs_PCA = prcomp(dataset_obj,scale=T, center = T)
summary(obs_PCA)
axes <- predict(obs_PCA, newdata = random_frog) #
head(axes, 4)

PC_load<- as.data.frame(obs_PCA$rotation[,1:3]) #this is just loading for each variable, predict is better

dat <- cbind(numfrog, axes)
dat

#modelling----
#for binomial y/n
m1<-glm.nb(tinc_Y_N ~ PC1 + PC2 + PC3, data = dat) #a negative binomial is a better fit here

#trying interactions, this is to see if any of the components interact
#let's try out every possible combination here
m2<-glm.nb(tinc_Y_N ~ PC1 * PC2 + PC3, data = dat)
m3<-glm.nb(tinc_Y_N ~ PC1 * PC3 + PC2, data = dat)
m4<-glm.nb(tinc_Y_N ~ PC2 * PC3 + PC1, data = dat)
m5<-glm.nb(tinc_Y_N ~ PC1 * PC2* PC3, data = dat)

#AIC MODEL SELECTION----
aic.cand<- list( 
  m1, 
  m2, 
  m3, 
  m4, 
  m5)
Cand.names <- c( "m1",
                 "m2", 
                 "m3", 
                 "m4",
                 "m5")
aictab(aic.cand, Cand.names, sort = T, second.ord = T) #this is based on AICc

#looking at the lowest AIC model, we can see that the interaction is
#not significant, so let's drop it.
summary(m2)

#that leaves us with the first model, where the components are added
summary(m1)

tab_model(m1, 
          transform = NULL,
          show.obs = F, 
          show.r2 = F, 
          show.est = T, 
          string.est = "Estimates", 
          dv.labels = "Tadpole presence (Y/N)", 
          pred.labels = c("PC1", "PC2", "PC3"))

####
#PCA ANALYSIS#
####

#################### ########### ##########
### PCA COMPONENT ANALYSIS---- ####
#################### ########### ##########

######################## ########### ###########
#### Index of Loadings (original data) ---- ####
################ ########### ########### ###########
#ILij = uij^2 * lambai^2
#where uij is jth **loading** of the ith **eigenvector**
#saved in "loadings" for first 3 components
### Function for permuting a data frame by swapping values in a column ###
data_permute = function(x){
  for (i in 1:ncol(x)){
    x[,i]=base::sample(x[,i],size=length(x[,i]),replace=F)}
  return(x)
}

### IL function for prcomp ###
# From Björklund (2019)
get_IL = function(x,num_PCs){
  eigs = ((x$sdev^2)^2)[1:num_PCs] 
  loads = ((data.frame(x$rotation))^2)[1:num_PCs]
  ILs = apply(loads,MARGIN = 1, FUN = function(x){x*eigs})
  return(ILs)
}

### Funtion list to array ###
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}
##############################################################################

### SET UP ###################################################################

#Alright here we're taking random_frog, because we've subsetted all of the tinc data out of it.
dataset_obj = random_frog
PCs_to_retain = 3
p_value_correction_method = "BH" 
#holm, bonferroni, hochberg, hommel, BY all attempted. 
number_randomizations = 1000
#number of traits
p = 14 #14 traits (see random_frog data)

##############################################################################

### ANALYSIS #################################################################

# True IL

obs_PCA = prcomp(dataset_obj,scale=T, center = T)
summary(obs_PCA) #okay, identical to pca from original
IL_obs = get_IL(obs_PCA,PCs_to_retain) #this is the idex loading from our original pca (only first three components)

# Nulls ILs #Let's randomize this bitch!

nulls = list()
for (i in 1:number_randomizations){
  nulls[[i]] = prcomp(data_permute(dataset_obj),scale=T, center = T)
}


## Now, let's calculate psi and phi to make sure that the sure that the 
#data is acturally structured!

##############################################################################

###PSI ----
psi1 = function (x) {
  psi = sum((x$sdev^2-1)^2)
  return(psi)
}

psi_obs<- psi1(obs_PCA)
psi_obs #10.2

psi1(nulls)

psi_nulls = lapply(nulls,psi1)
psi_nulls= as.numeric(psi_nulls)
psi_nulls

mean(psi_nulls > psi_obs) #nulls never overlap with observed
mean(psi_nulls)
##############################################################################
#phi----
#######
phi1 = function (x) {
  phi =  sqrt(((sum((x$sdev^2)^2)) - (p))/ ((p)*(p-1)))
  return(phi)
}

phi_obs<- phi1(obs_PCA) 
phi_obs #0.2368

phi_nulls = lapply(nulls,phi1)
phi_nulls= as.numeric(phi_nulls)
phi_nulls

mean(phi_nulls) #0.122

mean(phi_nulls > phi_obs) #The null is never bigger than our observed°!

##############

IL_nulls = lapply(nulls,get_IL,num_PCs=PCs_to_retain )
IL_nulls_array = list2ary(IL_nulls) # Save this to calculate p-values

pvalues_out = expand.grid(PC = rownames(IL_nulls_array[,,1]), trait = colnames(IL_nulls_array[,,1]))
pvals = c()

IL_mat<- as.matrix(IL_nulls_array)
sum(IL_mat)

for (i in 1:nrow(pvalues_out)){
  trait = rownames(IL_nulls_array[,,1])==pvalues_out[i,1]
  PC = colnames(IL_nulls_array[,,1])==pvalues_out[i,2]
  null_vals = IL_nulls_array[trait,PC,]
  
  obs_val = IL_obs[trait,PC]
  
  pvals[[i]] = mean(null_vals >= obs_val)
}


pvalues_out$p_values = pvals
pvalues_out$p_values_adj = p.adjust(pvals,method=p_value_correction_method) 
# adjust p values this is with BH may be a little intense
#pvalues_out$p_values_adj_1 = p.adjust(pvals,method=p_value_correction_method1) 


pvalues_out = pvalues_out[order(pvalues_out$PC),] ### THIS IS THE DATAFRAME WITH PVALUES

pvalues_sig <- pvalues_out %>%
  filter(p_values <= 0.05, 
         PC == "PC1") %>% #PC1 is all thats good for yn
  dplyr::select(-p_values_adj)

pvalues_sig

##############################################################################

### PLOTTING PCA ################################################################

# custom quantiles function for making ggplot boxplot
f = function(x) {
  r = quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 0.95))
  names(r) = c("ymin", "lower", "middle", "upper", "ymax")
  r
  
}

# create data frame for ggplotting

#from randomized data
IL_nulls_plot = data.frame(do.call(rbind,IL_nulls))
IL_nulls_plot$trait = rownames(IL_nulls_plot)
IL_nulls_plot = tidyr::gather(IL_nulls_plot,"PC","IL_score",1:(ncol(IL_nulls_plot)-1)) # 1 to number of PCs (last column with trait names should be excluded)

IL_nulls_PC1 <- IL_nulls_plot[grep("PC1", IL_nulls_plot$trait), ]
IL_nulls_PC1 <- IL_nulls_PC1 %>%
  rename("PC1" = trait ,
         "Trait" = PC)

#real observed data
IL_obs = data.frame(PC = rownames(IL_obs),IL_obs)
IL_obs_plot = tidyr::gather(IL_obs,"trait","IL_score",2:ncol(IL_obs))
IL_obs_PC1<- IL_obs_plot %>%
  filter(PC == "PC1")

#FOR PC1

ggplot() +
  stat_summary(data= IL_nulls_PC1, fun.data = f, 
               geom="boxplot",aes(y= IL_score,x= Trait))+ #random 95box
  annotate("rect", xmin =2.5, xmax = 3.5, ymin = 0, ymax = Inf, fill = "#DeF7E9", alpha = 0.4)+
  annotate("rect", xmin =5.5, xmax = 6.5, ymin = 0, ymax = Inf, fill = "#ffb8d2", alpha = 0.3)+
  annotate("rect", xmin =1.5, xmax = 2.5, ymin = 0, ymax = Inf, fill = "#ffb8d2", alpha = 0.3)+
  annotate("rect", xmin =11.5, xmax = 12.5, ymin = 0, ymax = Inf, fill = "#ffb8d2", alpha = 0.3)+
  annotate("rect", xmin =4.5, xmax = 5.5, ymin = 0, ymax = Inf, fill = "#7ae2ff", alpha = 0.3)+
  geom_point(data=IL_obs_PC1, aes(x = trait, y= IL_score), color = "darkgreen", shape = 4, size = 4)+
  ylab("Index Loading for PC1")+
  xlab("")+
  theme_bw()+
  theme(text = element_text(size = 14))

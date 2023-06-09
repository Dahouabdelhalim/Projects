######################################################################
# SUPPLEMENTARY DATA -  CODE: R FUNCTION SCRIPT FILE (SCRIPT TO RUN SNA OF THE MAIN TEXT)
#
#   This is the R function required to run the SNA code present in 'Familiarity, dominance, sex and season
# shape common waxbill social networks' (Gomes et al.)
#
# This function is adapted from the code in Franks et al. (2021,  DOI: 10.1111/2041-210X.13429)
#
#   Notes:
#    The codes in this script file must be ran in R before running 'SNA_script' file 
#    Please read the SNA sections of the main text and of the supplementary data, with the same names as the 
# section in this file, to guide you through this code.
#    R version 4.0.0 and 4.1.0
#    'lme4' packages is required
#    Created by Ana Cristina R. Gomes (ana.gomes@cibio.up.pt) with input from remaining authors
#
######################################################################


#### SNA: PREDICTORS OF NETWORK CENTRALITY - MAIN ANALYSIS, DSP-BC procedure ####

# this code is adapted to include random effects to GLMM

lmer.dsp_main <- function(formula, random_effect, data = NULL, nperm = 1000){
  #
  #
  # This function is adapted from the code in Franks et al. (2021,  DOI: 10.1111/2041-210X.13429)
  ## it allows to perform the double-semi partialling permutation procedure in bias-corrected models
  ## (DSP-BC), considering general linear mixed models, which includes random effects
  #
  # ARGUMENTS.
  # "formula" is the formula for the GLMM following the terminology of the function lmer of the lme4 R package,
  ## without stating the random effect variable
  # "random_effect" is the name of the column where values for random effect are present
  # "data" is the name of the dataframe where the variables for the GLMM are present
  #
  # VALUES. it returns:
  ## 'summ_table', returns a matrix with the GLMM results and p-values obtained from the DSP-BC procedure.
  ### 
  #
  #
  
  mf <- model.frame(formula, data=data) #get the model frame
  
  if(ncol(mf) < 3){
    stop("Double-semi-partialling only useful for multiple regression. Use X or Y permutation instead.")
  }
  
  random_effect_formula<-paste("+(1|", random_effect, ")", sep="")
  
  formula_lmer<-paste(colnames(mf)[1],"~",
                      paste(colnames(mf)[2:length(mf)], collapse ="+"), random_effect_formula, sep="")
  
  library(lme4) # version 1.1-27.1
  # run the original GLMM, which is a bias-corrected model
  orig_fit_model <- lme4::lmer(formula_lmer, data=data, REML=T) 
  t.perm <- matrix(ncol = (ncol(mf)-(1)), nrow = nperm) #set up matrix to hold permuted t-values
  
  
  for(i in 1:(ncol(mf)-1)){#for each predictor
    # DSP procedure
    
    # 1) run a GLMM similar to the original but removing one predictor, at a time, and the dependent variable, 
    ## and placing that predictor as dependent variable instead 
    y <- mf[,1] #get the response
    pred <- mf[,-1] #get all predictors
    npred<-length(pred)
    xi <- as.numeric(pred[,i]) #the target predictor
    zi <- pred[,-i] #other predictors
    ri <-data[,c(random_effect)]
    mfi <- as.data.frame(cbind(xi,zi,ri), stringsAsFactors = F) #predictor model frame
    mfi$ri<-as.factor(mfi$ri)
    
    formula_residual<-paste(colnames(mfi)[1],"~",
                            paste(colnames(mfi)[2:(length(colnames(mfi))-1)], collapse ="+"),
                            "+(1|", colnames(mfi)[length(colnames(mfi))],")", sep="")
    
    rm(eps)
    # 2) calculate model residuals, for this specific model
    eps <- try(residuals(lme4::lmer(formula_residual, data=mfi, REML=F)), silent=T) #get residuals
    eps<-try(scale(eps), silent=T)
    
    # 3) permute those residuals
    if(is.numeric(eps)==T){
      for(j in 1:nperm){
        eps.p <- sample(eps) #permute the residuals
        
        mfij <- as.data.frame(cbind(y,eps.p,zi,ri)) #permuted model frame
        mfij$ri<-as.factor(mfij$ri)
        
        
        formula_sample<-paste(colnames(mfij)[1],"~",
                              paste(colnames(mfij)[2:(length(colnames(mfij))-1)], collapse ="+"),
                              "+(1|", colnames(mfij)[length(colnames(mfij))],")", sep="")
        
        # 4) re-run the original GLMM, replacing the target predictor by one set of permuted residuals
        fit_p_model <- lme4::lmer(formula_sample, data=mfij, REML=T) #permuted linear model
        # 5) save the value of t corresponding to the effect of the targeted predictor
        t.perm[j,i] <- summary(fit_p_model)$coef[2,3] #extract T values for predictor
        
        print(paste(i, "effect_",j,"permutation",sep = "_"))
      }
    }}
  
  # register the coefficients of the original model
  summ_table <- cbind(summary(orig_fit_model)$coef[2:(ncol(mf)),],NA)
  colnames(summ_table)[ncol(summ_table)] <- c("Pr(Permuted)")
  
  # function to compute p-values from permuted values
  get_pval <- function(orig,rand){
    all_val <- c(orig,rand)
    ls <- mean(all_val <= orig)
    gr <- mean(all_val >= orig)
    min(c(ls,gr))*2
  }
  
  # compute significance of each predictor in the model
  summ_table[,ncol(summ_table)] <- sapply(1:ncol(t.perm),function(z){
    get_pval(summary(orig_fit_model)$coef[(z+1),3],t.perm[,z])
  })
  
  return(summ_table)
}
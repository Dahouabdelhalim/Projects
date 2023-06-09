#.........................................................
## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla

#.........................................................
#       ** BEHAVIOURAL DATA **====
#.........................................................  

    # SUMMARY OF THE SCRIPT:====

# For all Behavioural test:
        #...Normality test and transformation of data (log or sqrt)
        #...(optional for some test)Calculation of "Activity level" variable
        #    with a Principal Component Analyses 
        #...(optional for some test)Calculation of "Relative body turns" variable == "RBT"

    # Acronyms
        
        # Pi......Immediate response to the predator presence
        # Pd......Delayed response to the predator presence
#.........................................................

library(dplyr)
library(vegan)
library(ggplot2)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(PerformanceAnalytics ) 


#..IMPORT DATA "SYAT_MS_BH_ANIBEHdata.txt"........====

MS_SYAT_BH <- read.delim("SYAT_MS_BH_ANIBEHdata.txt", 
                            na = ".") 

str(MS_SYAT_BH)

#..STEP 1: DATA DESCRIPTION AND VARIABLE PROCESSING ........====
  
  # RISK ASSESSMENT BEHAVIOURAL DATA

  # ... Exploratory behaviour DATA  ====
  #......................................
   unique(MS_SYAT_BH$TEST)

    MS_SYAT_BH_EXPL <- filter(MS_SYAT_BH,
                              TEST == "EXPLORATION")

    #...Normality test and transformation of data
    
    EXPL_var <- c("T_INITIAL", "T_SECOND", "T_THIRD",
                  "T_FOURTH", "T_FIFTH", "T_SIXTH",  
                  "N_TURNS", "N_MOV","T_CROSS",
                  "RICH_POS")
    
    EXPL_Norm <- MS_SYAT_BH_EXPL %>%
        shapiro_test(EXPL_var) 
    
    EXPL_Norm
    
    chart.Correlation(MS_SYAT_BH_EXPL[,EXPL_var])
    
    # Log or sqrt transformation and test

    EXPL_var_toT <- c("T_INITIAL", "T_SECOND", "T_THIRD",
                      "T_FOURTH", "T_FIFTH", "T_SIXTH",  
                      "N_MOV","T_CROSS")

    for(i in unique (EXPL_var_toT)){
        
        MS_SYAT_BH_EXPL[,paste(i,"_log",sep="")] <- log10(MS_SYAT_BH_EXPL[,i]+1)
        MS_SYAT_BH_EXPL[,paste(i,"_sqrt",sep="")] <- sqrt(MS_SYAT_BH_EXPL[,i]+1)
    }    

    names_log <- names(MS_SYAT_BH_EXPL)[grepl("_log",names(MS_SYAT_BH_EXPL))]
    names_sqrt <- names(MS_SYAT_BH_EXPL)[grepl("_sqrt",names(MS_SYAT_BH_EXPL))]
    
    EXPL_T_Norm <- MS_SYAT_BH_EXPL %>%
        shapiro_test(names_log,names_sqrt)     
    
    EXPL_T_Norm # Test
    
    # Variables selected

    EXPL_var_T <- c("N_MOV_sqrt","N_TURNS",
                    "T_INITIAL_log", "T_SECOND_log",
                    "T_THIRD_log","T_FOURTH_log",
                    "T_FIFTH_log","T_SIXTH_log",
                    "T_CROSS_log",
                    "RICH_POS")

    chart.Correlation(MS_SYAT_BH_EXPL[,EXPL_var_T])
    
    EXPL_T_Norm <- MS_SYAT_BH_EXPL %>%
        shapiro_test(EXPL_var_T) 
    EXPL_T_Norm # Test

    #...Calculation of "Activity level" variable with a Principal
    # Component Analyses

    EXPL_PCAvar <- c("T_INITIAL_log",
                     "N_MOV_sqrt",
                     "N_TURNS", 
                     "RICH_POS")
    
    EXPL_PCA <- MS_SYAT_BH_EXPL[,EXPL_PCAvar]
    
    EXPL_PCA <- na.omit(EXPL_PCA)
    
    chart.Correlation(MS_SYAT_BH_EXPL[,EXPL_PCAvar]) 
    
    # Standardized the variable prior PCA, set center and scale equal to TRUE 

    EXPL_PCA.STZ <- scale(EXPL_PCA,
                          center=TRUE, scale=TRUE)
    
    # PCA 
    #........How many factors..................................................
    ev <- eigen(cor(EXPL_PCA.STZ)) # get eigenvalues
    ap <- parallel(subject=nrow(EXPL_PCA.STZ),
                   var=ncol(EXPL_PCA.STZ),
                   rep=100,cent=.05)
    nS <- nScree(x=ev$values,
                 aparallel=ap$eigen$qevpea)
    plotnScree(nS) 
    
    #...........................................................................
    
    PCA.EXPL <- principal(EXPL_PCA.STZ,
                          scores=TRUE, 
                          rotate="varimax",
                          nfactors=1)
    PCA.EXPL
    summary(PCA.EXPL)
    names(PCA.EXPL)
    
    PCA.EXPL.loadings<-as.data.frame(PCA.EXPL$loadings[])
    PCA.EXPL.loadings

    PCA.EXPL.Scores<-as.data.frame(PCA.EXPL$scores)
    
    # plot results
    
    biplot.psych(PCA.EXPL)
    
    # Inclusion of factor scores in the original MS_SYAT_BH_EXPL database
    names(MS_SYAT_BH_EXPL)

    MS_SYAT_BH_EXPL.1<-merge(MS_SYAT_BH_EXPL,
                           PCA.EXPL.Scores, 
                           by="row.names", 
                           all.x = TRUE,
                           no.dups = TRUE)
    
    MS_SYAT_BH_EXPL.1 <- dplyr::rename(MS_SYAT_BH_EXPL.1,
                              Activity= "PC1")
    
    MS_SYAT_BH_EXPL <- MS_SYAT_BH_EXPL.1 # original name of the database
    
    #...Calculation of "Relative body turns" variable == "RBT"
    
    RBT_EXPL_lm <- lm(N_TURNS ~ N_MOV_sqrt,
                 data=MS_SYAT_BH_EXPL)
    
    par(mfrow=c(2,2))
    plot(RBT_EXPL_lm)
    
    summary(RBT_EXPL_lm)
    
    
    RBT <- data.frame(RBT =residuals(RBT_EXPL_lm))
    
    RBT.match <- match(row.names(RBT),row.names(MS_SYAT_BH_EXPL))
    
    MS_SYAT_BH_EXPL[RBT.match,"RBT"] <- RBT$RBT

    #................................................
    # ... Foraging DATA ====
    #............................
  
    MS_SYAT_BH_FORAG <- filter(MS_SYAT_BH,
                               TEST == "FORAGING")
    
    #...Normality test and transformation of data
    
    # Inclusion of a dichotomus variable for eat/no eat
    
    MS_SYAT_BH_FORAG <- mutate(MS_SYAT_BH_FORAG,
                               EAT_BIN = ifelse(T_1_EAT < 3090,
                                                1,
                                                0))
    # Select eating individuals to asess:
    
    FORAG_var <- c("T_1_EAT", "N_PECKS")
    
    FORAG_Norm <- MS_SYAT_BH_FORAG %>%
        filter(EAT_BIN==1)%>%
        shapiro_test(FORAG_var) 
    
    FORAG_Norm
    
    # Log of the number of pecks
    names(MS_SYAT_BH_FORAG)
    MS_SYAT_BH_FORAG <- mutate(MS_SYAT_BH_FORAG,
                               N_PECKS_log = log10(N_PECKS +1))
    MS_SYAT_BH_FORAG %>%
        filter(EAT_BIN==1) %>%
        shapiro_test(N_PECKS_log)
    
    
    chart.Correlation(MS_SYAT_BH_FORAG[MS_SYAT_BH_FORAG$EAT_BIN==1,
                                       c("T_1_EAT",
                                          "N_PECKS_log")])

    # ... Playback conspecific alarm call test DATA ====    
    #..................................................
    
    MS_SYAT_BH_PLAYB <- filter(MS_SYAT_BH,
                               TEST == "PLAYBACK")
    
    #### Add an unique code by individual of the stimulus order (1=SYAT-PHYBON)
    #### or (2=PHYBON-SYAT)."STIM_ID_ORDER_f" indicate the order of 
    #### conspecific (SYAT) playback
 
    MS_SYAT_BH_PLAYB$STIM_ID_ORDER_f<- NA 
    
    for (i in unique(MS_SYAT_BH_PLAYB$ID)){
        
        Datos.i <- filter(MS_SYAT_BH_PLAYB,
                          ID==i)
        
        STIM_ID_ORDER_f <- unique(Datos.i[Datos.i$STIMULUS =="SYLATR",
                                 "STIM_ORDER"])
        
        MS_SYAT_BH_PLAYB[MS_SYAT_BH_PLAYB$ID==i, 
                        "STIM_ID_ORDER_f"] <- STIM_ID_ORDER_f
    }

    MS_SYAT_BH_PLAYB$STIM_ID_ORDER_f <- as.factor(MS_SYAT_BH_PLAYB$STIM_ID_ORDER_f)
    
    # Two subsets:
    
    MS_SYAT_BH_PLAYB_i <- filter(MS_SYAT_BH_PLAYB,
                                 TIME_FRAME == "DIRECT") %>%
        mutate(KEY_BEHAV = factor(KEY_BEHAV))# reset levels
    
    MS_SYAT_BH_PLAYB_d <- filter(MS_SYAT_BH_PLAYB,
                                      TIME_FRAME != "DIRECT")

    #...Preliminar description of the data
    
    ## Frequency of individuals that did not displace by stimulus
    ## and time frame
    
    MS_SYAT_BH_PLAYB_d <- mutate(MS_SYAT_BH_PLAYB_d,
                                   Displ_BIN_d = ifelse(N_MOV >0,
                                                          1,0)
    )
    
    FREQ_NO_DISPL <- xtabs(~Displ_BIN_d + STIMULUS + 
                               TIME_FRAME,
                           data=MS_SYAT_BH_PLAYB_d)
    
    FREQ_NO_DISPL

    #... Normality test and transformation of data
    
    PLAYB_DPLY_vars <- c("T_INITIAL","N_MOV",
                         "N_TURNS","RICH_POS")
    
    PLAYB_DPLY_Norm <- MS_SYAT_BH_PLAYB_d %>%
        group_by(TIME_FRAME)%>% 
        shapiro_test(PLAYB_DPLY_vars) 
    
    PLAYB_DPLY_Norm
    
    chart.Correlation(MS_SYAT_BH_PLAYB_d[
        MS_SYAT_BH_PLAYB_d$TIME_FRAME== "DELAYED_B",
        PLAYB_DPLY_vars]) # For behaviours recorded 2 minutes before stimulus
    
    chart.Correlation(MS_SYAT_BH_PLAYB_d[
        MS_SYAT_BH_PLAYB_d$TIME_FRAME== "DELAYED_A",
        PLAYB_DPLY_vars]) # For behaviours recorded 2 minutes after stimulus
    
    # Log or sqrt transformation and test
    
    PLAYB_var_toT <- PLAYB_DPLY_vars # All
    
    for(i in unique (PLAYB_var_toT)){
        
        MS_SYAT_BH_PLAYB_d[,paste(i,"_log",sep="")] <- log10(MS_SYAT_BH_PLAYB_d[,i]+1)
        MS_SYAT_BH_PLAYB_d[,paste(i,"_sqrt",sep="")] <- sqrt(MS_SYAT_BH_PLAYB_d[,i]+1)
    }    
    
    names_log <- names(MS_SYAT_BH_PLAYB_d)[grepl("_log",names(MS_SYAT_BH_PLAYB_d))]
    names_sqrt <- names(MS_SYAT_BH_PLAYB_d)[grepl("_sqrt",names(MS_SYAT_BH_PLAYB_d))]
    
    PLAYB_T_Norm <- MS_SYAT_BH_PLAYB_d %>%
        group_by(TIME_FRAME)%>%
        shapiro_test(names_log,names_sqrt,
                     PLAYB_DPLY_vars)%>%
        print(n = nrow(.)) # Test by time frame
    
    PLAYB_T_Norm_all <- MS_SYAT_BH_PLAYB_d %>%
        shapiro_test(names_log,names_sqrt,
                     PLAYB_DPLY_vars)%>%
        print(n = nrow(.)) # Test by stimulus

    #...Calculation of "Activity level" variable with a Principal Component Analyses
    
    # **Caution note**: Variables selected were not normal, distribution with highest
    # P value was chosen and check normality of the residuals of the
    # models
    
    PLAYB_PCAvar <- c("N_MOV_log","N_TURNS_sqrt",
                     "RICH_POS_sqrt")
    
    PLAYB_PCA <- MS_SYAT_BH_PLAYB_d[,PLAYB_PCAvar]
    
    PLAYB_PCA <- na.omit(PLAYB_PCA)
    
    chart.Correlation(PLAYB_PCA[,PLAYB_PCAvar])
  
    # Standardize the variable prior PCA, set center and scale equal to TRUE 

    PLAYB_PCA.STZ <- scale(PLAYB_PCA, 
                           center=TRUE, scale=TRUE)
    # PCA 
    #........How many factors..................................................
    ev <- eigen(cor(PLAYB_PCA.STZ)) # get eigenvalues
    ap <- parallel(subject=nrow(PLAYB_PCA.STZ),
                   var=ncol(PLAYB_PCA.STZ),
                   rep=100,cent=.05)
    nS <- nScree(x=ev$values,
                 aparallel=ap$eigen$qevpea)
    plotnScree(nS) 
    
    #...........................................................................
    
    PCA.PLAYB <- principal(PLAYB_PCA.STZ,
                          scores=TRUE, 
                          rotate="varimax",
                          nfactors=1)
    PCA.PLAYB
    summary(PCA.PLAYB)
    
    PCA.PLAYB.loadings<-as.data.frame(PCA.PLAYB$loadings[])
    PCA.PLAYB.loadings
    
    PCA.PLAYB.Scores<-as.data.frame(PCA.PLAYB$scores)
    
    # plot results
    
    biplot.psych(PCA.PLAYB)
    
    # Inclusion of factor scores in the original MS_SYAT_BH_EXPL database
    
    MS_SYAT_BH_PLAYB_d.1<-merge(MS_SYAT_BH_PLAYB_d,
                             PCA.PLAYB.Scores, 
                             by="row.names", 
                             all.x = TRUE,
                             no.dups = TRUE)
    
    MS_SYAT_BH_PLAYB_d.1 <- dplyr::rename(MS_SYAT_BH_PLAYB_d.1,
                                Activity= "PC1")
    
    MS_SYAT_BH_PLAYB_d <- MS_SYAT_BH_PLAYB_d.1 # original name of the database    
    
    #...Calculation of "Relative body turns" variable == "RBT"
    
    RBT_PLAYB_lm <- lm(N_TURNS_sqrt ~ N_MOV_log,
                 data=MS_SYAT_BH_PLAYB_d)
    
    par(mfrow=c(2,2))
    plot(RBT_PLAYB_lm)
    
    summary(RBT_PLAYB_lm)

    RBT <- data.frame(RBT =residuals(RBT_PLAYB_lm))
    
    RBT.match <- match(row.names(RBT),row.names(MS_SYAT_BH_PLAYB_d))
    
    MS_SYAT_BH_PLAYB_d[RBT.match,"RBT"] <- RBT$RBT

    # ... Sparrowhawk presence test DATA ====    
    #..........................................
 
    MS_SYAT_BH_PRED.0 <- filter(MS_SYAT_BH,
                               TEST == "PREDATOR_PRESENCE")
    
    # ID with no predation test
    
    No_Ptest <-unique(MS_SYAT_BH_PRED.0[which(is.na(MS_SYAT_BH_PRED.0$STIM_ORDER)),
                    "ID"])
    
    # Database with ID with predation test (Ptest)
    
    MS_SYAT_BH_PRED <- filter(MS_SYAT_BH_PRED.0,
                              ID != No_Ptest)
    
    #### Add an unique code by individual of the stimulus order (1=SPARROW HAWK-BOTTLE)
    #### or (2=BOTTLE-SPARROW HAWK)."STIM_ID_ORDER_f" indicate the order of 
    #### SPARROW HAWK presence
    #### Stimulus codes are PRED for SPARROW HAWK and CTRL for BOTTLE.

    MS_SYAT_BH_PRED$STIM_ID_ORDER_f<- NA 
    
    for (i in unique(MS_SYAT_BH_PRED$ID)){
        
        Datos.i <- filter(MS_SYAT_BH_PRED,
                          ID==i)
        
        STIM_ID_ORDER_f <- Datos.i[Datos.i$STIMULUS =="PRED",
                                   "STIM_ORDER"]
        
        MS_SYAT_BH_PRED[MS_SYAT_BH_PRED$ID==i, 
                        "STIM_ID_ORDER_f"] <- STIM_ID_ORDER_f
    }
    
    MS_SYAT_BH_PRED$STIM_ID_ORDER_f <- as.factor(MS_SYAT_BH_PRED$STIM_ID_ORDER_f)
    
    # Two subsets:
    
    MS_SYAT_BH_Pi <- filter(MS_SYAT_BH_PRED,
                                      TIME_FRAME == "DIRECT")
    
    MS_SYAT_BH_Pd<- filter(MS_SYAT_BH_PRED,
                                   TIME_FRAME != "DIRECT")
    
    #...Preliminar description of the data
    
    ## Frequency of individuals that did not displace by stimulus
    ## and time frame
    
    MS_SYAT_BH_PRED <- mutate(MS_SYAT_BH_PRED,
                              Displ_BIN = ifelse(N_MOV >0,
                                                 1,0))
    
    FREQ_NO_DISPL_PRED <- xtabs(~Displ_BIN + STIMULUS + 
                               TIME_FRAME,
                           data=MS_SYAT_BH_PRED)
    
    FREQ_NO_DISPL_PRED
    
    # ...... (Pi) Immediate responses DATA====
    #..............................
    
    #... Normality test and transformation of data
    
    Pi_vars <- c("T_INITIAL","N_MOV",
                         "N_TURNS","RICH_POS")
    
    Pi_Norm <- MS_SYAT_BH_Pi %>% 
        shapiro_test(Pi_vars) 
    
    Pi_Norm
    
    chart.Correlation(MS_SYAT_BH_Pi[,
            Pi_vars]) 

    # Log or sqrt transformation and test
    
    Pi_var_toT <- Pi_vars # All
    
    for(i in unique (Pi_var_toT)){
        
        MS_SYAT_BH_Pi[,paste(i,"_log",sep="")] <- log10(MS_SYAT_BH_Pi[,i]+1)
        MS_SYAT_BH_Pi[,paste(i,"_sqrt",sep="")] <- sqrt(MS_SYAT_BH_Pi[,i]+1)
    }    
    
    names_log <- names(MS_SYAT_BH_Pi)[grepl("_log",names(MS_SYAT_BH_Pi))]
    names_sqrt <- names(MS_SYAT_BH_Pi)[grepl("_sqrt",names(MS_SYAT_BH_Pi))]
    
    Pi_T_Norm <- MS_SYAT_BH_Pi %>%
        shapiro_test(names_log,names_sqrt,
                     Pi_vars)%>%
        print(n = nrow(.)) # Test all

    #...Calculation of "Activity level" variable with a Principal
    # Component Analyses
 
    Pi_PCAvar <- c("T_INITIAL_log","N_MOV_sqrt",
                        "N_TURNS_sqrt", "RICH_POS_sqrt")
    
    Pi_PCA <- MS_SYAT_BH_Pi[,Pi_PCAvar]
    
    Pi_PCA <- na.omit(Pi_PCA)
    
    chart.Correlation(Pi_PCA[,Pi_PCAvar])
    
    # Standardize the variable prior PCA, set center and scale equal to TRUE 

    Pi_PCA.STZ <- scale(Pi_PCA, 
                           center=TRUE, scale=TRUE)
    # PCA 
    #........How many factors..................................................
    ev <- eigen(cor(Pi_PCA.STZ)) # get eigenvalues
    ap <- parallel(subject=nrow(Pi_PCA.STZ),
                   var=ncol(Pi_PCA.STZ),
                   rep=100,cent=.05)
    nS <- nScree(x=ev$values,
                 aparallel=ap$eigen$qevpea)
    plotnScree(nS) 
    
    #...........................................................................
    
    PCA.Pi <- principal(Pi_PCA.STZ,
                           scores=TRUE, 
                           rotate="varimax",
                           nfactors=1)
    PCA.Pi
    summary(PCA.Pi)
    
    PCA.Pi.loadings<-as.data.frame(PCA.Pi$loadings[])
    PCA.Pi.loadings
 
    PCA.Pi.Scores<-as.data.frame(PCA.Pi$scores)
    
    # plot results
    
    biplot.psych(PCA.Pi)
    
    # Inclusion of factor scores in the original MS_SYAT_BH_EXPL database
    
    MS_SYAT_BH_Pi.1<-merge(MS_SYAT_BH_Pi,
                                  PCA.Pi.Scores, 
                                  by="row.names", 
                                  all.x = TRUE,
                                  no.dups = TRUE)
    
    MS_SYAT_BH_Pi.1 <- dplyr::rename(MS_SYAT_BH_Pi.1,
                                     Activity= "PC1")
  
    MS_SYAT_BH_Pi <- MS_SYAT_BH_Pi.1 # original name of the database    

    #...Calculation of "Relative body turns" variable == "RBT"
    
    RBT_Pi_lm <- lm(N_TURNS_sqrt ~ N_MOV_sqrt,
                       data=MS_SYAT_BH_Pi)
    
    par(mfrow=c(2,2))
    plot(RBT_Pi_lm)
    
    summary(RBT_Pi_lm)
    
    RBT <- data.frame(RBT =residuals(RBT_Pi_lm))
    
    RBT.match <- match(row.names(RBT),row.names(MS_SYAT_BH_Pi))
    
    MS_SYAT_BH_Pi[RBT.match,"RBT"] <- RBT$RBT
    
    # ...... (Pd) Delayed responses DATA====
    #..............................
    
    #... Normality test and transformation of data
    
    Pd_vars <- c("N_MOV","N_TURNS",
                       "RICH_POS")
    
    Pd_Norm <- MS_SYAT_BH_Pd %>%
        group_by(TIME_FRAME)%>% 
        shapiro_test(Pd_vars) 
    
    Pd_Norm
    
    chart.Correlation(MS_SYAT_BH_Pd[
        MS_SYAT_BH_Pd$TIME_FRAME== "DELAYED_B",
        Pd_vars]) # For behaviours recorded 2 minutes before stimulus
    
    chart.Correlation(MS_SYAT_BH_Pd[
        MS_SYAT_BH_Pd$TIME_FRAME== "DELAYED_A",
        Pd_vars]) # For behaviours recorded 2 minutes after stimulus
    
    # Log or sqrt transformation and test
    
    Pd_var_toT <- Pd_vars # All
    
    for(i in unique (Pd_var_toT)){
        
        MS_SYAT_BH_Pd[,paste(i,"_log",sep="")] <- log10(MS_SYAT_BH_Pd[,i]+1)
        MS_SYAT_BH_Pd[,paste(i,"_sqrt",sep="")] <- sqrt(MS_SYAT_BH_Pd[,i]+1)
    }    
    
    names_log <- names(MS_SYAT_BH_Pd)[grepl("_log",names(MS_SYAT_BH_Pd))]
    names_sqrt <- names(MS_SYAT_BH_Pd)[grepl("_sqrt",names(MS_SYAT_BH_Pd))]
    
    Pd_T_Norm <- MS_SYAT_BH_Pd %>%
        group_by(TIME_FRAME)%>%
        shapiro_test(names_log,names_sqrt,
                     Pd_vars)%>%
        print(n = nrow(.)) # Test by time frame
    
    Pd_T_Norm_all <- MS_SYAT_BH_Pd %>%
        shapiro_test(names_log,names_sqrt,
                     Pd_vars)%>%
        print(n = nrow(.)) # Test without grouping 
    
    #...Calculation of "Activity level" variable with a Principal
    # Component Analyses
    
    # **Caution note**: Variables selected were not normal, distribution with highest
    # P value was chosen and check normality of the residuals of the
    # models
    
    Pd_PCAvar <- c("N_MOV_sqrt","N_TURNS_sqrt",
                      "RICH_POS_sqrt")
    
    Pd_PCA <- MS_SYAT_BH_Pd[,Pd_PCAvar]
    
    Pd_PCA <- na.omit(Pd_PCA)
    
    chart.Correlation(Pd_PCA[,Pd_PCAvar])
    
    # Standardize the variable prior PCA, set center and scale equal to TRUE 

    Pd_PCA.STZ <- scale(Pd_PCA, 
                           center=TRUE, scale=TRUE)
    # PCA 
    #........How many factors..................................................
    ev <- eigen(cor(Pd_PCA.STZ)) # get eigenvalues
    ap <- parallel(subject=nrow(Pd_PCA.STZ),
                   var=ncol(Pd_PCA.STZ),
                   rep=100,cent=.05)
    nS <- nScree(x=ev$values,
                 aparallel=ap$eigen$qevpea)
    plotnScree(nS) 
    
    #...........................................................................
    
    PCA.Pd <- principal(Pd_PCA.STZ,
                           scores=TRUE, 
                           rotate="varimax",
                           nfactors=1)
    PCA.Pd
    summary(PCA.Pd)
    
    PCA.Pd.loadings<-as.data.frame(PCA.Pd$loadings[])
    PCA.Pd.loadings
    
    PCA.Pd.Scores<-as.data.frame(PCA.Pd$scores)
    
    # plot results
    
    biplot.psych(PCA.Pd)
    
    # Inclusion of factor scores in the original MS_SYAT_BH_Pd database
    
    MS_SYAT_BH_Pd.1<-merge(MS_SYAT_BH_Pd,
                                  PCA.Pd.Scores, 
                                  by="row.names", 
                                  all.x = TRUE,
                                  no.dups = TRUE)
    
    MS_SYAT_BH_Pd.1 <- dplyr::rename(MS_SYAT_BH_Pd.1,
                                     Activity= "PC1")
    
    MS_SYAT_BH_Pd <- MS_SYAT_BH_Pd.1 # original name of the database    
    
    #...Calculation of "Relative body turns" variable == "RBT"
    
    RBT_Pd_lm <- lm(N_TURNS_sqrt ~ N_MOV_sqrt,
                       data=MS_SYAT_BH_Pd)
    
    par(mfrow=c(2,2))
    plot(RBT_Pd_lm)
    
    summary(RBT_Pd_lm)
    
    RBT <- data.frame(RBT =residuals(RBT_Pd_lm))
    
    RBT.match <- match(row.names(RBT),row.names(MS_SYAT_BH_Pd))
    
    MS_SYAT_BH_Pd[RBT.match,"RBT"] <- RBT$RBT
    
    # ... ALARM CALLS ====
    #..........................
    
    MS_SYAT_BH_Acall <- filter(MS_SYAT_BH,
                               TEST %in% c("FORAGING",
                                           "EXPLORATION"))%>%
      dplyr::select(ID, TEST, N_A_CALL) %>%
      group_by (ID) %>%
      dplyr::summarise(RA_N_A_CALL = sum(N_A_CALL, na.rm = T),
                       NODATA_N_A_CALL =sum(is.na(N_A_CALL))
                       )
    
    # Bird07 was not monitored until the end of the experiment and did not
    # performed any alarm call. It was not included in the analyses of likelihood
    # of display alarm calls during the whole trial
    
    MS_SYAT_BH_Acall <- filter(MS_SYAT_BH_Acall,
                               !ID=="Bird07")
    
 #...Normality test and transformation of data
    
    # Inclusion of a dichotomous variable for displaying or not alarm calls
    
    MS_SYAT_BH_Acall <- mutate(MS_SYAT_BH_Acall,
                               Acall_BIN = ifelse(RA_N_A_CALL >0,
                                                  1,
                                                  0))
    # Select individuals displaying alarm calls to asses:"N_A_CALL"
    
    Acall_Norm <- MS_SYAT_BH_Acall %>%
      filter(Acall_BIN==1)%>%
      shapiro_test(RA_N_A_CALL) 
    
    Acall_Norm
    
    # Log of the number of pecks
    
    MS_SYAT_BH_Acall <- mutate(MS_SYAT_BH_Acall,
                               RA_N_A_CALL_log = log10(RA_N_A_CALL +1))
    MS_SYAT_BH_Acall %>%
      filter(Acall_BIN==1) %>%
      shapiro_test(RA_N_A_CALL_log)
    
    MS_SYAT_BH_Acall <- data.frame(MS_SYAT_BH_Acall)
    
    par(mfrow=c(1,1))
    
    hist(MS_SYAT_BH_Acall[MS_SYAT_BH_Acall$Acall_BIN==1,
                          "RA_N_A_CALL_log"],
         breaks=15)
#.........................................................
    
# ** INDIVIDUAL PARASITEMIA AND BODY CONDITION DATA **====
#.........................................................
    
    # SUMMARY:====
    
    #..IMPORT DATA...SYAT_MS_BH_BIOL_ANIBEHdata.txt
    #..STEP 1: DATA DESCRIPTION AND VARIABLE PROCESSING
    #...Calculation of a variable of percentage of body mass loss during experiment
    #...Normality test and transformation of data
    #...Calculation of "Size" variable with a Principal Component Analyses
    #...Calculation of body condition of individual "BC" 
    #..STEP 2: PARASITEMIA RESULTS
    #.........................................................
    
    #..IMPORT DATA...MS_SYAT_BH_BIOL_ANIBEHdata.txt........====
    
    MS_SYAT_BH_BIOL <- read.delim("SYAT_MS_BH_BIOL_ANIBEHdata.txt", 
                                  na = ".") 
    
    str(MS_SYAT_BH_BIOL)
    
    #..STEP 1: DATA DESCRIPTION AND VARIABLE PROCESSING ........====
    
    #...Calculation of a variable of percentage of body mass loss during experiment
    
    MS_SYAT_BH_BIOL <- mutate(MS_SYAT_BH_BIOL,
                              BM_Loss = (((BM_test_0-BM_test_end)/
                                            BM_test_0)*100)
    )
    
    #...Normality test and transformation of data
    
    BIOL_var <- c("BM_test_0", "BM_test_end",
                  "BM_Loss",
                  "Tarsus_L","Wing_L")
    
    BIOL_Norm <- MS_SYAT_BH_BIOL %>%
      shapiro_test(BIOL_var) 
    
    BIOL_Norm
    
    chart.Correlation(MS_SYAT_BH_BIOL[,BIOL_var])
    
    # Log or sqrt transformation and test
    
    BIOL_var_toT <- c("BM_test_0", "BM_test_end")
    
    for(i in unique (BIOL_var_toT)){
      
      MS_SYAT_BH_BIOL[,paste(i,"_log",sep="")] <- log10(MS_SYAT_BH_BIOL[,i]+1)
      MS_SYAT_BH_BIOL[,paste(i,"_sqrt",sep="")] <- sqrt(MS_SYAT_BH_BIOL[,i]+1)
    }    
    
    names_log <- names(MS_SYAT_BH_BIOL)[grepl("_log",names(MS_SYAT_BH_BIOL))]
    names_sqrt <- names(MS_SYAT_BH_BIOL)[grepl("_sqrt",names(MS_SYAT_BH_BIOL))]
    
    BIOL_T_Norm <- MS_SYAT_BH_BIOL %>%
      shapiro_test(BIOL_var,
                   names_log,
                   names_sqrt)     
    
    BIOL_T_Norm # Test
    
    # Variables selected
    
    BIOL_var_T <- c("BM_test_0_log", "BM_test_end_log",
                    "BM_Loss",
                    "Tarsus_L","Wing_L")
    
    chart.Correlation(MS_SYAT_BH_BIOL[,BIOL_var_T])
    
    BIOL_T_Norm <- MS_SYAT_BH_BIOL %>%
      shapiro_test(BIOL_var_T) 
    BIOL_T_Norm # Test
    
    #...Calculation of "Size" variable with a Principal Component Analyses
    
    BIOL_PCAvar <- c("Tarsus_L","Wing_L")
    
    BIOL_PCA <- MS_SYAT_BH_BIOL[,BIOL_PCAvar]
    
    BIOL_PCA <- na.omit(BIOL_PCA)
    
    chart.Correlation(MS_SYAT_BH_BIOL[,BIOL_PCAvar]) 
    
    # Standardized the variable prior PCA, set center and scale equal to TRUE 
    
    BIOL_PCA.STZ <- scale(BIOL_PCA,
                          center=TRUE, scale=TRUE)
    
    # PCA 
    
    PCA.BIOL <- principal(BIOL_PCA.STZ,
                          scores=TRUE, 
                          rotate="none",
                          nfactors=1)
    PCA.BIOL
    summary(PCA.BIOL)
    names(PCA.BIOL)
    
    PCA.BIOL.loadings<-as.data.frame(PCA.BIOL$loadings[])
    PCA.BIOL.loadings
    
    PCA.BIOL.Scores<-as.data.frame(PCA.BIOL$scores)
    
    # plot results
    
    biplot.psych(PCA.BIOL)
    
    # Inclusion of factor scores in the original MS_SYAT_BH_BIOL database
    names(MS_SYAT_BH_BIOL)
    
    MS_SYAT_BH_BIOL.1<-merge(MS_SYAT_BH_BIOL,
                             PCA.BIOL.Scores, 
                             by="row.names", 
                             all.x = TRUE,
                             no.dups = TRUE)
    
    MS_SYAT_BH_BIOL.1 <- dplyr::rename(MS_SYAT_BH_BIOL.1,
                                Size= "PC1")
    
    MS_SYAT_BH_BIOL <- MS_SYAT_BH_BIOL.1 # original name of the database
    
    #...Calculation of body condition of individual "BC" using
    # residuals of the regression of size against body mass at the starting
    # of the experiment
    
    BC_lm <- lm(BM_test_0 ~ Size, data = MS_SYAT_BH_BIOL)
    summary(BC_lm)
    par(mfrow=c(2,2))
    plot(BC_lm)
    
    # **Caution note**  There was a potential influential data (individual) but 
    # results did not change removing it so it was kept. Alternative tests using
    # Pexp0_log did not improve residuals.
    
    # Plot of the relationship
    ggplot(aes(y=BM_test_0, x=Size),
           data=MS_SYAT_BH_BIOL) +
      geom_point() +
      geom_smooth(method=lm,se=TRUE)+
      ggtitle("BC")
    
    # Inclusion of the residuals in the database
    BC <- data.frame(BC =residuals(BC_lm))
    
    BC.match <- match(row.names(BC),row.names(MS_SYAT_BH_BIOL))
    
    MS_SYAT_BH_BIOL[BC.match,"BC"] <- BC$BC        
    
    #.........................................................
    #..STEP 2: PARASITEMIA RESULTS ........====
    
    # Sample size of birds in each group with respect INFECTION AND MEDICATION
    # TREATMENT
    
    xtabs(~INF_BIN +TREAT,data=MS_SYAT_BH_BIOL)
    xtabs(~ MULT_INF +TREAT,data=MS_SYAT_BH_BIOL)

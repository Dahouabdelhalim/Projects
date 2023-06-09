#.........................................................
## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla
#...................................................................

    ## SUMMARY OF THE SCRIPT:----

# **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
# TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

# ...For PREDATOR EXPOSURE analyses use: 

# (B) MS_SYAT_BM_Pd =For delayed behaviour to predator presence

# In each database and for each variable two type of test:
    # ... EFFICACY OF THE TEST (APPENDIX A3):
    #for "order of the STIMULUS_f playback" and "type of STIMULUS_f"
    # ... PARASITISM EFFECT when the behavioural variable passed ok
    # the "efficacy of the test".

#...................................................................

library(dplyr)
library(vegan)
library(ggplot2)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(lsmeans)
library(afex)
library(blmeco)
library(LMERConvenienceFunctions)

library(PerformanceAnalytics )
library(ggpubr)
library(effects)
library(gridExtra)
library(Rmisc)

library(broom)
    
    set.seed(123)
#...................................................................
    ## PRELIMINARY PREPARATION OF DATABASES FOR ANALYSES----
#...................................................................

    #### PREDATOR EXPOSURE DELAYED BEHAVIOUR (Pd)
    
    ### Merge DELAYED Behaviour and Biological database

    PRED_d_BIOL_DB <- merge (MS_SYAT_BH_Pd, 
                           MS_SYAT_BH_BIOL,
                           by.x ="ID", 
                           by.y="ID",
                           all.x = TRUE,
                           no.dups=FALSE)
    str(PRED_d_BIOL_DB)

    ### Variable structure for analyses

    PRED_d_BIOL_DB <- mutate(PRED_d_BIOL_DB,
                           INF_BIN_f = factor(INF_BIN),
                           MULT_INF_f = factor(MULT_INF),
                           STIM_ORDER_f = factor(STIM_ORDER),
                           STIMULUS_f =factor(STIMULUS),
                           KEY_BEHAV = factor(KEY_BEHAV)# reset levels of factor
        )

#............................
    ##...PREDATOR EXPOSURE analyses ----
#...................................................................
    # **Caution note**: When we had warnings in glmer models because
    # model failed to converged we include the following method (bobyqa)
    # to improve convergence controlling errors

    mi.control <- glmerControl(optimizer="bobyqa", 
                           optCtrl=list(maxfun=100000))
#...................................................................
    #### DELAYED BEHAVIOUR TO PREDATOR EXPOSURE ####
#...................................................................

    ###....(1) Activity level ####

    # ...Previous tests

    # Identify outliers

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f, TIME_FRAME) %>%
        identify_outliers(Activity)

    # Normality tests

    PRED_d_BIOL_DB %>%
        group_by(INF_BIN_f, TREAT, 
                STIMULUS_f, TIME_FRAME) %>%
        shapiro_test(Activity)
    
    ggqqplot(PRED_d_BIOL_DB, "Activity",
             facet.by = c("STIMULUS_f"))
    
    ggplot (PRED_d_BIOL_DB,
            aes(x=Activity))+
        geom_histogram()+
        facet_wrap(vars(STIMULUS_f))

    # Check normality also in residuals of the models

    # Levene's tests

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f, TIME_FRAME) %>%
        levene_test(Activity ~ STIM_ID_ORDER_f,
                    center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f, TIME_FRAME) %>%
        levene_test(Activity ~ TREAT*INF_BIN_f,
                    center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f, TIME_FRAME) %>%
        levene_test(Activity ~ TREAT,
                    center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f, TIME_FRAME) %>%
        levene_test(Activity ~ INF_BIN_f,
                 center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f, TIME_FRAME) %>%
        levene_test(Activity ~ MULT_INF_f,
                    center=mean) 

# Homogeneity of covariances assumption
# the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
# are unequal, we ignored it.

    box_m(PRED_d_BIOL_DB[,
                          "Activity",
                          drop = FALSE],PRED_d_BIOL_DB$TREAT) 

    PRED_d_BIOL_DB_MULT_G <- PRED_d_BIOL_DB %>%
        filter(!is.na(MULT_INF_f))# select not NA values

    box_m(PRED_d_BIOL_DB_MULT_G[,
                                 "Activity",
                                 drop = FALSE], 
          PRED_d_BIOL_DB_MULT_G$MULT_INF_f)

    box_m(PRED_d_BIOL_DB[,
                          "Activity",
                          drop = FALSE],
          PRED_d_BIOL_DB$STIM_ID_ORDER_f)

    #........ EFFICACY OF THE TEST #### 
    
    #........ Previous visualization====

    #........Differences between groups
    
    ggplot(data= PRED_d_BIOL_DB,
           aes(y=Activity, 
               x= TIME_FRAME))  +
        geom_boxplot(aes(fill=as.factor(STIMULUS_f))) + 
        geom_point() +
        theme_classic()+
        facet_grid(STIM_ID_ORDER_f ~ STIMULUS_f)+
        ggtitle("PREDATOR ORDER")+
        theme(legend.position = "none")

    #........Individual response

    ggplot(data=PRED_d_BIOL_DB, 
           aes(x=TIME_FRAME,
               y=Activity,
               group=ID)) +
        geom_point(aes(fill=STIMULUS_f, 
                       size=1.5,alpha=0.5))+
        geom_line(aes(color=ID))+
        facet_wrap(vars(STIM_ID_ORDER_f,STIMULUS_f))+
        theme(legend.position = "none")

    #........ Repeated measures test ====

    Act_OR_PRED_DLY_0 = aov_car(Activity ~ STIMULUS_f*STIM_ID_ORDER_f*
                                    TIME_FRAME +
                           Error(ID/STIMULUS_f*TIME_FRAME),
                       data=PRED_d_BIOL_DB,
                       type=3)

    summary(Act_OR_PRED_DLY_0) # No significant
    Act_OR_PRED_DLY_0

    # Bird activity change was not affected by the order of the predator
    
    Act_OR_PRED_DLY = aov_car(Activity ~ STIMULUS_f*TIME_FRAME +
                                    Error(ID/STIMULUS_f*TIME_FRAME),
                                data=PRED_d_BIOL_DB,
                                type=3)
    
    summary(Act_OR_PRED_DLY) # No significant
    Act_OR_PRED_DLY
    
    # Birds did not change activity after any stimulus.
    
    # **Cautious note** We also tested potential influence of parasitism
    # in the efficacy of the test but was not significant and therefore 
    # we do not included these analyses.

    ### ........ PARASITISM EFFECT: ####
    
    # * NOT TESTED because test did not worked* See above efficacy of the test*

    
#...................................................................

    ###.... (2) Relative body turns (RBT) ####

    # ...Previous tests

    # Identify outliers

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f,TIME_FRAME) %>%
            identify_outliers(RBT)

    # Normality tests

    PRED_d_BIOL_DB %>%
        group_by(INF_BIN_f, TREAT, 
             STIMULUS_f, TIME_FRAME) %>%
        shapiro_test(RBT)

    ggqqplot(PRED_d_BIOL_DB, "RBT",
         facet.by = c("STIMULUS_f","TIME_FRAME"))

    ggplot (PRED_d_BIOL_DB,
            aes(x=RBT))+
        geom_histogram()+
        facet_wrap(vars(STIMULUS_f,TIME_FRAME))

    # Check normality also in residuals of the models

    # Levene's tests

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f,TIME_FRAME) %>%
        levene_test(RBT ~ STIM_ID_ORDER_f,
                center=mean)

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f,TIME_FRAME) %>%
        levene_test(RBT ~ TREAT*INF_BIN_f,
                center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f,TIME_FRAME) %>%
        levene_test(RBT ~ TREAT,
                center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f,TIME_FRAME) %>%
        levene_test(RBT ~ INF_BIN_f,
                center=mean) 

    PRED_d_BIOL_DB %>%
        group_by(STIMULUS_f,TIME_FRAME) %>%
        levene_test(RBT ~ MULT_INF_f,
                center=mean) 

    # Homogeneity of covariances assumption
    # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
    # are unequal, we ignored it.

    box_m(PRED_d_BIOL_DB[,
                          "RBT",
                          drop = FALSE],PRED_d_BIOL_DB$TREAT) 

    PRED_d_BIOL_DB_MULT_G <- PRED_d_BIOL_DB %>%
        filter(!is.na(MULT_INF_f))# select not NA values

    box_m(PRED_d_BIOL_DB_MULT_G[,
                                 "RBT",
                                 drop = FALSE],
          PRED_d_BIOL_DB_MULT_G$MULT_INF_f)

    box_m(PRED_d_BIOL_DB[,
                          "RBT",
                          drop = FALSE],
          PRED_d_BIOL_DB$STIM_ID_ORDER_f)

    ###........  EFFICACY OF THE TEST ####
    
    #........ Previous visualization====
    
    #........Differences between groups
    
    ggplot(data= PRED_d_BIOL_DB,
           aes(y=RBT, 
               x= TIME_FRAME))  +
        geom_boxplot(aes(fill=as.factor(STIMULUS_f))) + 
        geom_point() +
        theme_classic()+
        facet_grid(STIM_ID_ORDER_f ~ STIMULUS_f)+
        ggtitle("PREDATOR ORDER")+
        theme(legend.position = "none")
    
    #........Individual response
    
    ggplot(data=PRED_d_BIOL_DB, 
           aes(x=TIME_FRAME,
               y=RBT,
               group=ID)) +
        geom_point(aes(fill=STIMULUS_f, 
                       size=1.5,alpha=0.5))+
        geom_line(aes(color=ID))+
        facet_wrap(vars(STIM_ID_ORDER_f,STIMULUS_f))+
        theme(legend.position = "none")
    
    #........ Repeated measures test ====
    
    RBT_OR_PRED_DLY_0 = aov_car(RBT ~ STIMULUS_f*STIM_ID_ORDER_f*
                                    TIME_FRAME +
                                    Error(ID/STIMULUS_f*TIME_FRAME),
                                data=PRED_d_BIOL_DB,
                                type=3)
    
    summary(RBT_OR_PRED_DLY_0) # No significant
    RBT_OR_PRED_DLY_0
    
    # Bird RBT change was not affected by the order of the predator
    
    RBT_OR_PRED_DLY = aov_car(RBT ~ STIMULUS_f*TIME_FRAME +
                                  Error(ID/STIMULUS_f*TIME_FRAME),
                              data=PRED_d_BIOL_DB,
                              type=3)
    
    summary(RBT_OR_PRED_DLY) # No significant
    RBT_OR_PRED_DLY
    
    # Birds did not change RBT after any stimulus.
    
    # **Cautious note** We also tested potential influence of parasitism
    # in the efficacy of the test but was not significant and therefore 
    # we do not included these analyses.
    
    ###........ PARASITISM EFFECT ####
    
    # * NOT TESTED because test did not worked* See above efficacy of the test*
    
 
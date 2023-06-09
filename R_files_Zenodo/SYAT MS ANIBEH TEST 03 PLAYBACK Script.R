#................................................................

## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla

#................................................................
        # SUMMARY OF THE SCRIPT:====

# **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
# TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

# ...For PLAYBACK CONSPECIFIC CALL analyses use:====
        ## MS_SYAT_BH_PLAYB_i database ====
        ## MS_SYAT_BH_PLAYB_d database ====

# In each database and for each variable two types of test:
        # ... EFFICACY OF THE TEST (APPENDIX A3):
    #for "order of the STIMULUS_f playback" and "type of STIMULUS_f"
        # ... PARASITISM EFFECT when the behavioural variable passed ok
    # the "efficacy of the test".

#......................................................

library(dplyr)
library(vegan)
library(ggplot2)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(elrm)
library(lsmeans)
library(emmeans)
library(afex)
library(blmeco)
library(LMERConvenienceFunctions)
library(DHARMa)
        
library(PerformanceAnalytics )
library(ggpubr)
library(effects)
library(gridExtra)
library(Rmisc)

library(jtools)
library(multcomp)
        
        set.seed(123)
#......................................................
 
    ## PRELIMINARY PREPARATION OF DATABASES FOR ANALYSES ----
#......................................................

        ### There are no empty rows:
        
       sum(ifelse(is.na(MS_SYAT_BH_PLAYB_d$T_INITIAL),1,0))


        ### Merge Behaviour (IMMEDIATE and DELAYED) and Biological database
       
            ####(A) PLAYBACK CONSPECIFIC CALL IMMEDIATE BEHAVIOUR
       
        PLAYB_i_BIOL_DB <- merge (MS_SYAT_BH_PLAYB_i, 
                                   MS_SYAT_BH_BIOL,
                                   by.x ="ID", 
                                   by.y="ID",
                                   all.x = TRUE,
                                   no.dups=FALSE)
        str(PLAYB_i_BIOL_DB)
       
            ####(B) PLAYBACK CONSPECIFIC CALL DELAYED BEHAVIOUR
       
        PLAYB_d_BIOL_DB <- merge (MS_SYAT_BH_PLAYB_d, 
                                MS_SYAT_BH_BIOL,
                                by.x ="ID", 
                                by.y="ID",
                                all.x = TRUE,
                                no.dups=FALSE)
        str(PLAYB_d_BIOL_DB)
        
        ### Variable structure for analyses
       
            ####(A) PLAYBACK CONSPECIFIC CALL IMMEDIATE BEHAVIOUR
        
        PLAYB_i_BIOL_DB <- mutate(PLAYB_i_BIOL_DB,
                    INF_BIN_f = factor(INF_BIN),
                    MULT_INF_f = factor(MULT_INF),
                    STIM_ORDER_f = factor(STIM_ORDER),
                    STIMULUS_f =factor(STIMULUS)
                    )
        
            ####(B) PLAYBACK CONSPECIFIC CALL DELAYED BEHAVIOUR
        
        PLAYB_d_BIOL_DB <- mutate(PLAYB_d_BIOL_DB,
                    INF_BIN_f = factor(INF_BIN),
                    MULT_INF_f = factor(MULT_INF),
                    STIM_ORDER_f = factor(STIM_ORDER),
                    STIMULUS_f =factor(STIMULUS),
                    # Dichotomous variable for some behaviors
                    Lat_BIN = ifelse(T_INITIAL_log < max(T_INITIAL_log),
                                     0,# bird displaced
                                     1), # bird did not displace
                    Activ_BIN = ifelse(Activity > min(Activity),
                                       1, 
                                       0)# no activity
        )

#......................................................
        
    ##...PLAYBACK CONSPECIFIC CALL analyses ----
#......................................................
    # **Caution note**: When we had warnings in glmer models because
        # model failed to converged we included the following method (bobyqa)
        # to improve convergence controlling errors

        mi.control <- glmerControl(optimizer="bobyqa", 
                                   optCtrl=list(maxfun=100000))
        
    ####(A) PLAYBACK CONSPECIFIC CALL BEHAVIOUR ####
  #......................................................
  
        ###.... Freezing behaviour ACCORDING TO RISK ####
        
        #........ EFFICACY OF THE TEST ####
        
        # Most individuals did not freeze to any stimulus. Only 10 birds displayed 
        # freezing behaviour immediately after an acoustic stimulus.
        
        xtabs(~KEY_BEHAV+STIMULUS_f
              ,data=PLAYB_i_BIOL_DB)
        
        ### We classified freezing behaviour of individuals in two categories:
        
        ###   (1) when individuals froze only with the conspecific playback
        ###   (0) Other alternatives
        
        # Transformation from long to wide format:
        
        PLAYB_i_BIOL_DB_w <-reshape(PLAYB_i_BIOL_DB, 
                                   direction='wide',
                                   timevar="STIMULUS",
                                   v.names= "KEY_BEHAV",
                                   idvar="ID") 
        
        # We created new dichotomous variable:
        
        PLAYB_i_BIOL_DB_w <- mutate(PLAYB_i_BIOL_DB_w,
                                   FROZE_ID = ifelse(KEY_BEHAV.SYLATR==1 &
                                                       KEY_BEHAV.PHYBON==0,
                                                     1,0))
        
        # Test only with the order of the stimulus
        
        xtabs(~FROZE_ID+STIM_ID_ORDER_f
              ,data=PLAYB_i_BIOL_DB_w)
        
        FROZE_ORS_PLAYB_glm = glm(FROZE_ID ~ STIM_ID_ORDER_f, 
                              family=binomial,
                              data=PLAYB_i_BIOL_DB_w)
        
        drop1(FROZE_ORS_PLAYB_glm,
              test="Chisq")
        
        # The order of the STIMULUS_f did not influence the likelihood of a bird
        # using freezing behaviour according to risk
        
        # ... Effect of the type of stimulus alone
        
        xtabs(~KEY_BEHAV+STIMULUS_f
              ,data=PLAYB_i_BIOL_DB)
        
        FROZE_STIM_PLAYB_glmer = glmer(KEY_BEHAV ~ STIMULUS_f +
                                     (1 | ID), 
                                   family=binomial,
                                   data=PLAYB_i_BIOL_DB)
        
        drop1(FROZE_STIM_PLAYB_glmer,
              test="Chisq")
        
        plot(allEffects (FROZE_STIM_PLAYB_glmer))
        
        # overdispersion test
        
        testDispersion(FROZE_STIM_PLAYB_glmer)
        
        # Birds froze more likely in response to a conspecific playback than to
        # the control song.
        # (LRT=4.33, P =0.04, Only 10 individuals displayed freezing behaviour
        # when challenged with an acoustic stimulus, 8 of which responded to 
        # low-pitched conspecific alarm calls and 2 to control stimulus ).
       
        # Correlation with body condition (BC)
        
        FROZE_BC_PLAYB_glm = glm(FROZE_ID ~ BC, 
                             family=binomial,
                             data=PLAYB_i_BIOL_DB_w)
        
        drop1(FROZE_BC_PLAYB_glm,
              test="Chisq") # not significant
        
        # overdispersion
        phi <- sum((residuals(FROZE_BC_PLAYB_glm, 
                              type="pearson"))^2)/FROZE_BC_PLAYB_glm$df.residual
        print(c("Pearson overdispersion =", 
                round(phi, 3)), quote=FALSE)
        
        testDispersion(FROZE_BC_PLAYB_glm)
        
        #........ PARASITISM EFFECT ####
        
        #........ STATUS OF INFECTION ====
        
        FX1 <- xtabs(~FROZE_ID+interaction(INF_BIN_f,TREAT)
                    ,data=PLAYB_i_BIOL_DB_w) 
        
        FX1dat <- data.frame(INF_BIN_f = rep(0:1, 2), 
                             TREAT = rep(c("PQ","Water"),
                                         each=2), 
                             FROZE_ID = c(FX1[2, ]),
                            n = colSums(FX1))
        FX1dat  # view collapsed data set
        
        # Exact logistic regression because some cells were empty
        
        FX1exact = elrm(FROZE_ID/n ~ INF_BIN_f +
                         TREAT +
                         INF_BIN_f:TREAT, 
                       interest = ~INF_BIN_f +
                         TREAT +
                         INF_BIN_f:TREAT,
                       iter= 1e+05, 
                       burnIn=2000, data=FX1dat,
                       r=2)
        
        F1exact2 <- update(FX1exact, iter = 1e+05, 
                            burnIn = 5000)
        
        summary(F1exact2) # interaction not significant
        
        # plot(F1exact2)
        
        # Testing simple effects without interaction
        
        FX1exact3 = elrm(FROZE_ID/n ~ INF_BIN_f +
                          TREAT , 
                        interest = ~INF_BIN_f +
                          TREAT,
                        iter= 1e+05, 
                        burnIn=2000, data=FX1dat,
                        r=2)
        summary(FX1exact3) # infection status not significant
                           # treatment almost significant
        
        # plot(FX1exact3)
        
        #........ MULTIPLE STATUS OF INFECTION ====
        
        # Not tested because too few birds exhibiting the freezing behaviour
        
        FMLGX1 <- xtabs(~FROZE_ID+interaction(MULT_INF,TREAT)
                     ,data=PLAYB_i_BIOL_DB_w) 
        
    ####(B) PLAYBACK CONSPECIFIC CALL DELAYED BEHAVIOUR ####
    
        #...............................................................              
        ###.... (1) Activity  #### 
        
        # ...Previous tests
        
        # Identify outliers
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          identify_outliers(Activity)
        
        # Normality tests
        
        PLAYB_d_BIOL_DB %>%
          group_by(INF_BIN_f, TREAT, 
                   STIMULUS_f, TIME_FRAME) %>%
          shapiro_test(Activity)
        
        ggqqplot(PLAYB_d_BIOL_DB, "Activity",
                 facet.by = c("STIMULUS_f","TIME_FRAME"))
        
        ggplot (PLAYB_d_BIOL_DB,
                aes(x=Activity))+
          geom_histogram()+
          facet_wrap(vars(STIMULUS_f, TIME_FRAME))
        # There are many birds with no activity after any playback
        
        # Check normality also in residuals of the models
        
        # Levene's tests
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ STIM_ID_ORDER_f,
                      center=mean)
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ TREAT*INF_BIN_f,
                      center=mean)
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ TREAT,
                      center=mean)
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ INF_BIN_f,
                      center=mean)
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ MULT_INF_f,
                      center=mean)
        
        # Homogeneity of covariances assumption
        # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
        # are unequal, we ignored it.
        
        box_m(PLAYB_d_BIOL_DB[,
                              "Activity",
                              drop = FALSE],PLAYB_d_BIOL_DB$TREAT) 
        
        PLAYB_d_BIOL_DB_MULT_G <- PLAYB_d_BIOL_DB %>%
          filter(!is.na(MULT_INF_f))# select not NA values
        
        box_m(PLAYB_d_BIOL_DB_MULT_G[,
                                     "Activity",
                                     drop = FALSE], PLAYB_d_BIOL_DB_MULT_G$MULT_INF_f) 
        
        box_m(PLAYB_d_BIOL_DB[,
                              "Activity",
                              drop = FALSE], PLAYB_d_BIOL_DB$STIM_ID_ORDER_f)
        
        # ** Cautious note** Heterocedasticity issue. High frequency of motionless 
        # individuals in some tests made variable distribution too skewed to meet
        # linear general model assumptions. Double approach for analyses: 
        
        # Visualization plot of individual changes of activity
        
        ggplot(data=PLAYB_d_BIOL_DB, aes(x=TIME_FRAME, 
                                         y=Activity,
                                         group=ID)) +
          geom_point(aes(fill=STIMULUS_f, size=1.5,alpha=0.5))+
          geom_line(aes(color=ID))+
          facet_wrap(vars(STIM_ID_ORDER_f,STIMULUS_f))+
          theme(legend.position = "none")
        
        #...... (a) Binomial analyses: dichotomous variable depending on birds 
        #were active or not
        
        #...... (b) Separated activity analyses with subset of only active birds
        
        #...... (a) Activity(YES/NO)...binomial analyses ====
        
        #........ EFFICACY OF THE TEST ####
        
        #...Repeated measures test
        
        Act_BIN_PLAYB_order.0 = glmer(Activ_BIN ~ STIMULUS_f +
                                        TIME_FRAME +
                                        STIM_ID_ORDER_f +
                                        STIMULUS_f:TIME_FRAME +
                                        STIMULUS_f:STIM_ID_ORDER_f+
                                        STIM_ID_ORDER_f:TIME_FRAME+
                                        STIM_ID_ORDER_f:STIMULUS_f:TIME_FRAME +
                                        (1 | ID), 
                                      family=binomial,
                                      data=PLAYB_d_BIOL_DB)
        
        drop1(Act_BIN_PLAYB_order.0, 
              test="Chisq")# Not significant
        # stimulus order
        
        # overdispersion
        dispersion_glmer(Act_BIN_PLAYB_order.0)
        
        # Model without stimulus order
        
        Act_BIN_PLAYB_eff = glmer(Activ_BIN ~ STIMULUS_f +
                                    TIME_FRAME +
                                    STIMULUS_f:TIME_FRAME +
                                    (1 | ID), 
                                  family=binomial,
                                  data=PLAYB_d_BIOL_DB)
        
        drop1(Act_BIN_PLAYB_eff,
              test="Chisq") #
        
        # There are NO differences in the delayed response between types
        # of stimulus. We cannot assume that the behaviour is induced by test.
        
        Act_BIN_PLAYB_eff_time = glmer(Activ_BIN ~ STIMULUS_f +
                                         TIME_FRAME +
                                         (1 | ID), 
                                       family=binomial,
                                       data=PLAYB_d_BIOL_DB)
        
        drop1(Act_BIN_PLAYB_eff_time,
              test="Chisq") 
        plot(allEffects(Act_BIN_PLAYB_eff_time))
        
        # Birds decreased activity equally after any playback test.
        
        # **Cautious note** We also tested potential influence of parasitism
        # in the efficacy of the test but was not significant and therefore 
        # we do not included these analyses.
        
        #........ PARASITISM EFFECT TEST ####
        
        # * NOT TESTED because test did not worked* See above efficacy of the test*
        
        #...... (b) Activity level in active birds====
        
        ### SUBSET DATA:
        
        # We included individuals that had always any activity
        
        ID.NO.Activ <- filter (PLAYB_d_BIOL_DB,
                               Activ_BIN==0) %>%
          dplyr::select(ID,STIM_ID_ORDER_f, STIMULUS_f,
                        TIME_FRAME )%>%
          mutate(ID = factor(ID)) # reset level
        
        length(unique(ID.NO.Activ$ID)) # 20 birds with no activity
        # on any occasion
        
        ID.No.Activ.ID <- levels(ID.NO.Activ$ID)
        
        PLAYB_d_BIOL_DB_ACTIV <- filter(PLAYB_d_BIOL_DB,
                                        !ID %in% ID.No.Activ.ID)%>%
          mutate(ID=factor(ID))
        
        length(levels(PLAYB_d_BIOL_DB_ACTIV$ID))# sample size:23
        
        # ...Previous tests
        
        # Identify outliers
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          identify_outliers(Activity)
        
        # 2 potential OUTLIERS
        
        # Normality test
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME, INF_BIN_f,TREAT) %>%
          shapiro_test(Activity) 
        
        ggqqplot(PLAYB_d_BIOL_DB_ACTIV, "Activity",
                 facet.by = c("STIMULUS_f","TIME_FRAME"))
        
        ggplot (PLAYB_d_BIOL_DB_ACTIV,
                aes(x=Activity))+
          geom_histogram()+
          facet_grid(TIME_FRAME~ STIMULUS_f*INF_BIN_f,
                     labeller = "label_both")
        
        # Levene's test
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ STIM_ID_ORDER_f,
                      center=mean)
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ TREAT*INF_BIN_f,
                      center=mean)
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ TREAT,
                      center=mean) 
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ INF_BIN_f,
                      center=mean) 
        
        PLAYB_d_BIOL_DB_ACTIV %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(Activity ~ MULT_INF_f,
                      center=mean)
        
        # Homogeneity of covariances assumption
        
        # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
        # are unequal, ignore it. However, if significant and you have unequal sample sizes, 
        # the test is not robust (https://en.wikiversity.org/wiki/Box%27s_M, Tabachnick & Fidell, 2001).
        
        box_m(PLAYB_d_BIOL_DB_ACTIV[,
                                    "Activity",
                                    drop = FALSE], PLAYB_d_BIOL_DB_ACTIV$TREAT)
        
        box_m(PLAYB_d_BIOL_DB_ACTIV[,
                                    "Activity",
                                    drop = FALSE], PLAYB_d_BIOL_DB_ACTIV$INF_BIN_f)   
        
        PLAYB_d_BIOL_DB_ACTIV_MULT_G <- PLAYB_d_BIOL_DB_ACTIV %>%
          filter(!is.na(MULT_INF_f))# select not NA values
        
        box_m(PLAYB_d_BIOL_DB_ACTIV_MULT_G[,
                                           "Activity",
                                           drop = FALSE], PLAYB_d_BIOL_DB_ACTIV_MULT_G$MULT_INF_f)
        
        box_m(PLAYB_d_BIOL_DB_ACTIV[,
                                    "Activity",
                                    drop = FALSE], PLAYB_d_BIOL_DB_ACTIV$STIM_ID_ORDER_f)
        
        #........ EFFICACY OF THE TEST ####
        
        #........ Previous visualization ====
        
        #........Differences between groups : EFFICACY OF THE TEST
        
        ggplot(data= PLAYB_d_BIOL_DB_ACTIV,
               aes(y=Activity, 
                   x= TIME_FRAME))  +
          geom_boxplot(aes(fill=as.factor(STIMULUS_f))) + 
          geom_point() +
          theme_classic()+
          facet_grid(STIM_ID_ORDER_f  ~ STIMULUS_f)+
          ggtitle("STIM_ID_ORDER_f")+
          theme(legend.position = "none")
        
        #.........Individual response
        
        ggplot2::ggplot(data=PLAYB_d_BIOL_DB_ACTIV, aes(x=ID, 
                                                        y=Activity,
                                                        group=ID)) +
          
          theme_classic()+
          geom_point(aes(color=TIME_FRAME, size=1.5,alpha=0.5))+
          geom_line()+
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 90))+
          facet_grid(STIMULUS_f~INF_BIN_f)
        
        #........ Repeated measures analyses ====
        
        Activ_PLAYB_order.0 =aov_car(Activity ~ STIMULUS_f*TIME_FRAME*STIM_ID_ORDER_f +
                                       Error(ID/STIMULUS_f*TIME_FRAME),
                                     data=PLAYB_d_BIOL_DB_ACTIV,
                                     type=3)
        
        summary(Activ_PLAYB_order.0) # significant influence of the order 
        # of stimulus
        
        # Residual visualization
        plots <- list()
        for (i in c(1:4)){
          p1 =ggqqplot(as.numeric(residuals(Activ_PLAYB_order.0$lm)[,i]), 
                       title=paste(colnames(residuals(Activ_PLAYB_order.0$lm))[i]))
          plots[[i]] <- p1 
        }
        multiplot(plotlist = plots, 
                  cols = 2)
        
        # Pos-hoc analyses to explore further effect of order of the conspecific
        # stimulus
        
        Activ.two.way <- PLAYB_d_BIOL_DB_ACTIV %>%
          dplyr::select(STIM_ID_ORDER_f,Activity,
                        STIMULUS_f, TIME_FRAME, ID) %>%
          group_by(STIM_ID_ORDER_f) %>%
          anova_test(dv = Activity, 
                     wid = ID, 
                     within = c(STIMULUS_f, TIME_FRAME))
        Activ.two.way
        
        # **When the conspecific stimulus is presented in the second place **
        # **activity differed between type of stimulus**.
        
        # POST-HOC to check time-frame differences in each type of stimulus
        
        # Subset sample size for individuals where conspecific alarm call
        # was played at second:
        
        PLAYB_d_BIOL_DB_ACTIV2 <- filter(PLAYB_d_BIOL_DB_ACTIV,
                                         STIM_ID_ORDER_f==2)%>%
          mutate(Activity_stdz = scale(Activity, 
                                       center=TRUE,scale=TRUE))
        
        PLAYB_d_BIOL_DB_ACTIV2 %>%
          group_by(STIMULUS_f)%>%
          pairwise_t_test(Activity ~ TIME_FRAME,
                          paired = TRUE, 
                          p.adjust.method = "bonferroni" )
        
        # ** Birds reduced activity only after the conspecific playback**
        # Sample size 
        length(unique(PLAYB_d_BIOL_DB_ACTIV2$ID)) # n = 11 
        
        #........ PARASITISM EFFECT TEST ####
        
        # * NOT TESTED because the subset of individuals responding was too small 
        # to allow statistical tests 
        
        #...................................................................              
        ###.... (2) Relative body turns (RBT) #### 
        
        # ...Previous tests
        
        # Identify outliers
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          identify_outliers(RBT)
        
        # 21 potential outliers but no one is extreme so we kept all.
        
        # Normality test
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          shapiro_test(RBT) # PHYBON before significant
        
        # Check normality also in residuals of the models
        
        ggqqplot(PLAYB_d_BIOL_DB, "RBT",
                 facet.by = c("STIMULUS_f","TIME_FRAME")) 
        
        ggplot(PLAYB_d_BIOL_DB, aes(x=RBT))+
          geom_histogram()+
          facet_wrap(vars(STIMULUS_f, TIME_FRAME))
        
        # Levene's test
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(RBT ~ STIM_ID_ORDER_f,
                      center=mean)
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(RBT ~ TREAT,
                      center=mean) 
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(RBT ~ INF_BIN_f,
                      center=mean) 
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(RBT ~ TREAT*INF_BIN_f,
                      center=mean) 
        
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(RBT ~ MULT_INF_f,
                      center=mean) 
        PLAYB_d_BIOL_DB %>%
          group_by(STIMULUS_f, TIME_FRAME) %>%
          levene_test(RBT ~ TREAT*MULT_INF_f,
                      center=mean) 
        
        # Homogeneity of covariances assumption
        # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
        # are unequal, we ignored it.
        
        box_m(PLAYB_d_BIOL_DB[,
                              "RBT",
                              drop = FALSE], PLAYB_d_BIOL_DB$TREAT)
        
        box_m(PLAYB_d_BIOL_DB[,
                              "RBT",
                              drop = FALSE], PLAYB_d_BIOL_DB$INF_BIN_f)   
        
        box_m(PLAYB_d_BIOL_DB[,
                              "RBT",
                              drop = FALSE], PLAYB_d_BIOL_DB$MULT_INF_f)
        
        box_m(PLAYB_d_BIOL_DB[,
                              "RBT",
                              drop = FALSE], PLAYB_d_BIOL_DB$STIM_ID_ORDER_f)
        
        #........ EFFICACY OF THE TEST ####
        
        #........ Previous visualization ====
        
        #........Differences between groups : EFFICACY OF THE TEST
        
        ggplot(data= PLAYB_d_BIOL_DB,
               aes(y=RBT, 
                   x= TIME_FRAME))  +
          geom_boxplot(aes(fill=as.factor(STIMULUS_f))) + 
          geom_point() +
          theme_classic()+
          facet_grid(STIM_ID_ORDER_f  ~ STIMULUS_f)+
          ggtitle("STIM_ID_ORDER_f")+
          theme(legend.position = "none")
        
        #.........Individual response
        
        ggplot2::ggplot(data=PLAYB_d_BIOL_DB, aes(x=ID,
                                                  y=RBT,
                                                  group=ID)) +
          theme_classic()+
          geom_point(aes(color=TIME_FRAME, size=1.5,alpha=0.5))+
          geom_line()+
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 90))+
          facet_grid(STIMULUS_f~INF_BIN_f)
        
        #........ Repeated measures analyses====
        
        RBT_ORS_PLAYB.0 =aov_car(RBT ~ STIMULUS_f*TIME_FRAME*STIM_ID_ORDER_f +
                                   Error(ID/STIMULUS_f*TIME_FRAME),
                                 data=PLAYB_d_BIOL_DB,
                                 type=3)
        
        summary(RBT_ORS_PLAYB.0) # not significant
        
        # Residual visualization
        plots <- list()
        for (i in c(1:4)){
          p1 =ggqqplot(as.numeric(residuals(RBT_ORS_PLAYB.0$lm)[,i]), 
                       title=paste(colnames(residuals(RBT_ORS_PLAYB.0$lm))[i]))
          plots[[i]] <- p1 
        }
        multiplot(plotlist = plots, 
                  cols = 2)
        
        # **Caution note**: Bad residuals for RBT responses after PHYBON playback
        # Deletion of some potential outliers did not change results. All data remained.
        
        # Model without stimulus order
        
        RBT_ORS_PLAYB_eff =aov_car(RBT ~ STIMULUS_f*TIME_FRAME +
                                     Error(ID/STIMULUS_f*TIME_FRAME),
                                   data=PLAYB_d_BIOL_DB,
                                   type=3)
        
        summary(RBT_ORS_PLAYB_eff) # not significant
        
        # Residual visualization
        plots <- list()
        for (i in c(1:4)){
          p1 =ggqqplot(as.numeric(residuals(RBT_ORS_PLAYB_eff$lm)[,i]), 
                       title=paste(colnames(residuals(RBT_ORS_PLAYB_eff$lm))[i]))
          plots[[i]] <- p1 
        }
        multiplot(plotlist = plots, 
                  cols = 2)
        # **Caution note**: Bad residuals for RBT responses after PHYBON playback.
        
        # There are NO differences in the delayed response between types
        # of stimulus. We cannot assume that the behaviour is induced by test.
        
        # **Cautious note** We also tested potential influence of parasitism
        # in the efficacy of the test but was not significant and therefore 
        # we do not included these analyses.
        
        #........ PARASITISM EFFECT TEST ####
        
        # * NOT TESTED because test did not worked* See above efficacy of the test*
            
#......................................................
    ###.... (3) Initial latency ####
    
    ### Preliminary description of the data
        
    # For conspecific playback latency analyses we discarded 
    # individuals which did not move in the two minutes previous 
    # to the STIMULUS_f because it is not possible to measure latency 
    # to stop displacement when they are already motionless 
      
        ID.NO.MOVE.Bf <- filter ( PLAYB_d_BIOL_DB,
                                  TIME_FRAME=="DELAYED_B",
                                  Lat_BIN==1) %>%
            dplyr::select(ID,STIM_ID_ORDER_f, STIMULUS_f, TIME_FRAME )%>%
            mutate(ID = factor(ID))
        
        length(unique(ID.NO.MOVE.Bf$ID)) # 16 birds were discarded
        
        ID.No.Bf <- levels(ID.NO.MOVE.Bf$ID)
        
        xtabs( ~ STIM_ID_ORDER_f + STIMULUS_f,
               data=ID.NO.MOVE.Bf)
        
        # Discard these individuals
        
        PLAYB_Lat_Bf_sel <- dplyr::filter(PLAYB_d_BIOL_DB,
                                          !ID %in% ID.No.Bf)%>%
            mutate(ID =factor(ID), # reset levels
                   T_INITIAL_log_stdz = scale(T_INITIAL_log,
                                              center=TRUE, scale=TRUE)
            )
        
        # Sample size
        length(levels(PLAYB_Lat_Bf_sel$ID)) # 27 individuals for analysis
        
        # ...Previous tests
        
        # Identify outliers
        
        PLAYB_Lat_Bf_sel %>%
            group_by(STIMULUS_f, TIME_FRAME) %>%
            identify_outliers(T_INITIAL_log)
        
        # Normality test
        
        PLAYB_Lat_Bf_sel %>%
            group_by(INF_BIN_f, TREAT, 
                     STIMULUS_f, TIME_FRAME) %>%
            shapiro_test(T_INITIAL_log)
        
        ggqqplot(PLAYB_Lat_Bf_sel, "T_INITIAL_log",
                 facet.by = c("STIMULUS_f","TIME_FRAME"))
        ggplot (PLAYB_Lat_Bf_sel,
                aes(x=T_INITIAL_log))+
            geom_histogram()+
            facet_wrap(vars(STIMULUS_f, TIME_FRAME))
        
        # Check normality also in residuals of the models

        # Levene's test

        PLAYB_Lat_Bf_sel %>%
            group_by(STIMULUS_f, TIME_FRAME) %>%
            levene_test(T_INITIAL_log ~ STIM_ID_ORDER_f,
                        center=mean)
        
        PLAYB_Lat_Bf_sel %>%
            group_by(STIMULUS_f, TIME_FRAME) %>%
            levene_test(T_INITIAL_log ~ TREAT*INF_BIN_f,
                        center=mean)
        
        PLAYB_Lat_Bf_sel %>%
            group_by(STIMULUS_f, TIME_FRAME) %>%
            levene_test(T_INITIAL_log ~ TREAT,
                        center=mean) 
        PLAYB_Lat_Bf_sel %>%
            group_by(STIMULUS_f, TIME_FRAME) %>%
            levene_test(T_INITIAL_log ~ INF_BIN_f,
                        center=mean)  

        PLAYB_Lat_Bf_sel %>%
            group_by(STIMULUS_f, TIME_FRAME) %>%
            levene_test(T_INITIAL_log ~ MULT_INF_f,
                        center=mean) 
        
        # Homogeneity of covariances assumption
        # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
        # are unequal, we ignored it.
        
        box_m(PLAYB_Lat_Bf_sel[,
                            "T_INITIAL_log",
                            drop = FALSE], PLAYB_Lat_Bf_sel$TREAT)

        PLAYB_Lat_Bf_sel_MULT_G <- PLAYB_Lat_Bf_sel %>%
            filter(!is.na(MULT_INF_f))# select not NA values
        
        box_m(PLAYB_Lat_Bf_sel_MULT_G[,
                            "T_INITIAL_log",
                            drop = FALSE], PLAYB_Lat_Bf_sel_MULT_G$MULT_INF_f)

        box_m(PLAYB_Lat_Bf_sel[,
                            "T_INITIAL_log",
                            drop = FALSE], PLAYB_Lat_Bf_sel$STIM_ID_ORDER_f)

        #........ EFFICACY OF THE TEST ####
        
        #........ Previous visualization ====
        
        #........Differences between groups : EFFICACY OF THE TEST
        
        ggplot(data= filter(PLAYB_Lat_Bf_sel),
               aes(y=T_INITIAL_log, 
                   x= TIME_FRAME))  +
            geom_boxplot(aes(fill=as.factor(STIM_ID_ORDER_f))) + 
            geom_point() +
            theme_classic()+
            facet_wrap(vars(STIM_ID_ORDER_f, STIMULUS_f))+
            ggtitle("SYAT PLAYBACK ORDER") +
            theme(axis.text.x = element_text(angle = 90)) 

        #........Individual response
        
        ggplot(data=PLAYB_Lat_Bf_sel, aes(x=ID, 
                                                y=T_INITIAL_log,
                                                group=ID)) +
            theme_classic()+
            geom_point(aes(color=TIME_FRAME, size=1.5,alpha=0.5))+
            geom_line()+
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 90))+
            facet_wrap("STIMULUS_f")
        
        #........ Repeated measures test ====
        
        Lat_PLAYB_order.0 =aov_car(T_INITIAL_log ~ STIMULUS_f*TIME_FRAME*STIM_ID_ORDER_f +
                                      Error(ID/STIMULUS_f*TIME_FRAME),
                                  data=PLAYB_Lat_Bf_sel,
                                  type=3)
        
        summary(Lat_PLAYB_order.0) # almost significant p=0.051
        
        # Post-hoc analyses to explore further effect of order of the conspecific
        # stimulus
        
        Lat.two.way <- PLAYB_Lat_Bf_sel %>%
            dplyr::select(STIM_ID_ORDER_f,T_INITIAL_log,
                          STIMULUS_f, TIME_FRAME, ID) %>%
            group_by(STIM_ID_ORDER_f) %>%
            anova_test(dv = T_INITIAL_log , 
                       wid = ID, 
                       within = c(STIMULUS_f, TIME_FRAME))
        Lat.two.way
        
        # ORDER 1 = SYAT playback displayed first # not significant
        # ORDER 2 = SYAT playback displayed second # significant
        
# **When the conspecific stimulus is presented in the second place **
# **latency differed between type of stimulus**.
        
    # POST-HOC to check time-frame differences in each type of stimulus
        
        # Subset sample size for individuals where conspecific alarm call
        # was played at second:
        
        PLAYB_Lat_Bf_sel_02=filter(PLAYB_Lat_Bf_sel,
                                   STIM_ID_ORDER_f==2) %>%
            mutate(ID=factor(ID))
        
        PLAYB_Lat_Bf_sel_02 %>%
            group_by(STIMULUS_f)%>%
            pairwise_t_test(T_INITIAL_log ~ TIME_FRAME,
                            paired = TRUE, 
                            p.adjust.method = "bonferroni" )
        
# ** Birds delayed response only after the conspecific playback**
    
        #........ PARASITISM EFFECT TEST ####
        
        # * NOT TESTED because the subset of individuals responding was too small 
        # to allow statistical tests 
        
        # sample size
        length(unique(PLAYB_Lat_Bf_sel_02$ID)) 
        
    # 14 individuals
        

#................................................................

## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla

#................................................................

    ## SUMMARY OF THE SCRIPT:----

# **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
# TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

    # ...For ALARM CALLS analyses use = MS_SYAT_BH_Acall database

#................................................................

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
library(PerformanceAnalytics )
library(ggpubr)

library(lsmeans)
library(afex)
library(jtools)
library(multcomp)
    
    set.seed(123)
#...................................................................
    ## PRELIMINARY PREPARATION OF DATABASES FOR ANALYSES ----
#...................................................................

    ##...ALARM CALL BEHAVIOUR

    Acall.DATA <- !is.na(MS_SYAT_BH_Acall$Acall_BIN)

    MS_SYAT_BH_Acall_f <- MS_SYAT_BH_Acall[Acall.DATA,]
    
    ### Merge Behavioural and Biological database
    
    Acall_BIOL_DB <- merge (MS_SYAT_BH_Acall_f, 
                            MS_SYAT_BH_BIOL,
                            by.x ="ID", 
                            by.y="ID",
                            all.x = TRUE,
                            no.dups=FALSE)
    str(Acall_BIOL_DB)
    
    ### Variable structure for analyses
    
    Acall_BIOL_DB <- mutate(Acall_BIOL_DB,
                            INF_BIN_f = as.factor(INF_BIN),
                            MULT_INF = as.factor(MULT_INF),
                            Acall_BIN_f = as.factor(Acall_BIN))

    ### Subset of dataframe for individuals displaying alarm calls
    
    Acall_BIOL_DB_SEL <- filter (Acall_BIOL_DB,
                                 Acall_BIN == 1)%>%
        mutate(RA_N_A_CALL_log_stdz = scale(RA_N_A_CALL_log,
                                    center=TRUE, scale=TRUE))
    
    # Differences between birds that elicit or not alarm calls with BC
    
    BC_ALARM.glm = glm(Acall_BIN  ~ BC, 
                     family=binomial,
                     data=Acall_BIOL_DB)
    
    # overdispersion test
    phi <- sum((residuals(BC_ALARM.glm, type="pearson"))^2)/BC_ALARM.glm$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)),quote=FALSE)
    
    drop1(BC_ALARM.glm,~.,
          test="Chisq") # no differences
#...................................................................
    ##...ALARM CALL analyses ----
#...................................................................
    
    ###...(1) ALARM CALL RISK ####
    
    # Number of individuals displaying alarm calls during 
    # risk-taking assessment time-frame
    
    length(unique(Acall_BIOL_DB_SEL$ID)) # 16

    #........ STATUS OF INFECTION ====
    
    A1 <- xtabs(~Acall_BIN_f+interaction(INF_BIN_f,TREAT)
                ,data=Acall_BIOL_DB)
    A1 # frequency table
    
    Acall_G_glm.0 = glm(Acall_BIN_f ~ INF_BIN_f*TREAT, 
                        family=binomial,
                        data=Acall_BIOL_DB)
    
    drop1(Acall_G_glm.0,~.,
          test="Chisq") # Significant interaction
    
    par(mfrow=c(2,3))
    plot(Acall_G_glm.0,
         ask = F, which = 1:6)
    
    # overdispersion test
    
    phi <- sum((residuals(Acall_G_glm.0, type="pearson"))^2)/Acall_G_glm.0$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)

    # POST-HOC
    
    # **For INFECTED**
    
    Acall_Risk_G_INF <- glm(Acall_BIN_f ~ TREAT,
                            family=binomial,
                            data=filter(Acall_BIOL_DB,
                                        INF_BIN_f==1)
    )
    
    drop1( Acall_Risk_G_INF, test="Chisq") # Significant differences
    
    # Frequencies
    
    xtabs(~ TREAT + Acall_BIN_f , 
          data=filter(Acall_BIOL_DB,
                      INF_BIN_f==1))
    
    # General infected individuals untreated were less likely to display 
    # alarm calls than individuals treated with PQ
   
    # **For NOT INFECTED**
    
    Acall_Risk_G_notINF <- glm(Acall_BIN_f ~ TREAT,
                               family=binomial,
                               data=filter(Acall_BIOL_DB,
                                           INF_BIN_f==0)
    )
    
    drop1( Acall_Risk_G_notINF, test="Chisq") # Not significant differences
    summary(Acall_Risk_G_notINF)
    
    # PLOT WITH RESULTS
    
    Acall_BIOL_DB <- mutate(Acall_BIOL_DB,
                            INF_G_TREAT = paste(INF_BIN_f, 
                                                TREAT, sep=""))
    
    ggplot(Acall_BIOL_DB, 
           aes(x = TREAT,
               fill = Acall_BIN_f)) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent)+
      theme_classic()+
      facet_wrap(vars(INF_BIN_f))
    
    p.adjust(c(0.034,0.388), method="bonferroni")
    
    #........ MULTIPLE STATUS OF INFECTION ====
 
    AMLG1 <- xtabs(~Acall_BIN_f+interaction(MULT_INF,TREAT)
                ,data=Acall_BIOL_DB)
    
    AMLG1dat <- data.frame(MULT_INF = rep(0:2, 2), 
                        TREAT = rep(c("PQ","Water"),
                                    each = 3), 
                        Acall_BIN_f = AMLG1[2, ], 
                        n = colSums(AMLG1))
    AMLG1dat  # view collapsed data set
    
    # Exact logistic regression because some cells were empty
    
    AMLG1exact1.0 = elrm(Acall_BIN_f/n ~ MULT_INF +
                     TREAT +
                      MULT_INF:TREAT, 
                   interest = ~MULT_INF +
                     TREAT +
                     MULT_INF:TREAT,
                   iter= 1e+05, 
                   burnIn=2000, data=AMLG1dat,
                   r=2)
    
    AMLG1exact2.0 <- update(AMLG1exact1.0, iter = 5e+06, 
                       burnIn = 4000)
    
    summary(AMLG1exact2.0) # not significant interaction
    # plot(AMLG1exact2.0)
 
    # Without interaction, there are no empty cells
    
    AMLG1_MULT_INF <- xtabs(~Acall_BIN_f+MULT_INF
                          ,data=Acall_BIOL_DB)
    AMLG1_TREAT <- xtabs(~Acall_BIN_f+TREAT
                       ,data=Acall_BIOL_DB)
    
    # normal logistic regression
    
    Acall_Risk_MLG <- glm(Acall_BIN_f ~ MULT_INF + TREAT,
                               family=binomial,
                               data=Acall_BIOL_DB)
    
    drop1( Acall_Risk_MLG, test="Chisq") 
    summary(Acall_Risk_MLG)
    
    plot(allEffects(Acall_Risk_MLG))
    
    library(multcomp)# for glht, mcp
    summary(glht(Acall_Risk_MLG, mcp(MULT_INF="Tukey")), 
            test=adjusted(type="bonferroni"))
    
    ###...(1) RA_N_A_CALL_log #### 
    
    # only data where number of calls could be measured during the trial
    Acall_BIOL_DB_SELdata <- filter(Acall_BIOL_DB_SEL,
                                    NODATA_N_A_CALL==0)
    # Sample size 
    
    length(unique(Acall_BIOL_DB_SELdata$ID)) #n= 15
    # **Caution note**: low sample size to test influence of infection and treatment
   
    ###...FDR CORRECTION BY MULTIPLE TESTING DURING THE RISK-TAKING ASSESSMENT TESTS  ----
    
    ALARM_p_unadjust <-c(0.045, # Acall_BIN_f ~ INF_BIN_f*TREAT 
                         0.03,0.32) # Acall_BIN_f ~ MULT_INF + TREAT
    
    test <- paste("test", 1:length(ALARM_p_unadjust), 
                  sep = "_")
    
    FDR_ALARM_results <- data.frame(test, ALARM_p_unadjust)
    
    
    FDR_ALARM_results <- arrange(FDR_ALARM_results, ALARM_p_unadjust)
    
    FDR_ALARM_results$FDR_ALARM_p_adjust <-p.adjust(FDR_ALARM_results$ALARM_p_unadjust,
                                                    method="fdr")
    
    FDR_ALARM_results
    
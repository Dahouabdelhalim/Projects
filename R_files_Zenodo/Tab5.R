#Table 5: Main Results: Comparison of prices, quantities, and surplus with at and RTP pricing.
###############################################################################
#Making latex tables
joinall <- transform(forlatex, category=ifelse(forlatex$category=="free","unconst.",paste(forlatex$category)))
joinall <- transform(joinall, category=ifelse(joinall$category=="rps_100","rps",paste(joinall$category)))
joinall <- joinall[ order(joinall[,1], joinall[,2],joinall[,3],joinall[,4],joinall[,5],joinall[,6],rev(joinall[,7])), ]
df <- joinall[which(joinall$scen!=2),]

################################################################################
##Table S1 main results in Appendix
###############################################################################
df$pricing<- factor(df$pricing, levels = c("flat", "dynamic"))
df <- df[order(df[,1],df[,2],df[,3],df[,4],df[,5],df[,6],rank(df[,7])), ]

df$Renewable_Share <- format(round(as.numeric(df$Renewable_Share),0), nsmall=0, big.mark=",") 
for (i in c("dcs","evcost","dps","tots","high","mid","inf","dyntot")){
  df[[i]] <- format(round(as.numeric(df[[i]]),1), nsmall=1, big.mark=",") 
}
for (i in c("mean","meanq","sd")){
  df[[i]] <- format(round(as.numeric(df[[i]]),0), nsmall=0, big.mark=",") 
}

#scen only one 
df <- transform(df, scen=ifelse(df$scen==1,"Optimistic", ifelse(df$scen==2,"Moderate","Pessimistic")))
df$scen<- factor(df$scen, levels = c("Optimistic", "Moderate","Pessimistic"))
scen1lastrow <- c("\\\\multirow{-2}{*}{Optimistic}")
scen2lastrow <- c("\\\\multirow{-2}{*}{Moderate}")
scen3lastrow <- c("\\\\multirow{-2}{*}{Pessimistic}")
df$scen <- factor(df$scen, levels = c(levels(df$scen), ""))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen1lastrow))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen2lastrow))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen3lastrow))

len <- nrow(df)
for (i in seq(1,len,1)){
  #find the last row to replace with multicolumn
  firstrow <- as.numeric(which(df[i,6] == df[i+1,6]))
  lastrow <- as.numeric(which(df[i,6] != df[i+1,6] | df[i,6] == df[len,6]))
  ifelse(length(firstrow)== 1,{
    df[i,6] <- ""}, 
    ifelse(length(lastrow)==1, {
      df[i,6] <- ifelse(df[i,6]=="Optimistic",scen1lastrow,ifelse(df[i,6]=="Moderate",scen2lastrow,scen3lastrow))},
      next))
}

##baseline 0
replacetext <- c("\\\\multicolumn{7}{c}{------------------------------ B  a  s  e  l  i  n  e  ------------------------------}\\\\tabularnewline")
ftx <- 0
index <- which(as.numeric(df$dcs)==ftx & as.numeric(df$evcost)==ftx 
               & as.numeric(df$dps)==ftx & as.numeric(df$tots)==ftx 
               & as.numeric(df$high)==ftx & as.numeric(df$mid)==ftx 
               & as.numeric(df$inf)==ftx & as.numeric(df$dyntot)==ftx, arr.ind=TRUE)
df[index,c("dcs")] <- replacetext
df[index,c("dps")] <- ""
df[index,c("evcost")] <- ""
df[index,c("tots")] <- ""
df[index,c("high")] <- ""
df[index,c("mid")] <- ""
df[index,c("inf")] <- ""

##dynprice 
ftx <- 0
index0 <- which(as.numeric(df$dyntot)==ftx, arr.ind=TRUE)
index1 <- index0+1
dynline <- paste("\\\\multirow{-2}{*}{",df[index1,19],"}\\\\tabularnewline", sep = "")
df[index1,c("dyntot")] <- paste(dynline,"\\\\cline{5-9}", sep = "")
df[index0,c("dyntot")] <- ""
df$pricing <-  sapply(df$pricing, str_to_title)

#100%renewable rotate
rpslastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering  100\\\\% Clean}}}")
fossillastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering Fossil}}}")
unconst.lastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering Unconstrained}}}")
df$category<- factor(df$category, levels = c("fossil", "rps","unconst."))
df$category <- factor(df$category, levels = c(levels(df$category), ""))
df$category <- factor(df$category, levels = c(levels(df$category), fossillastrow))
df$category <- factor(df$category, levels = c(levels(df$category), rpslastrow))
len <- nrow(df)
for (v in c("fossil", "rps")){
  for (i in seq(1,len,1)){
    #find the last row to replace with multicolumn
    firstrow <- as.numeric(which(df[i+1,4] == v & df[i,4] == v))
    lastrow <- as.numeric(which((df[i+1,4] != v | df[len,4] == v) & df[i,4] == v))
    ifelse(length(firstrow)== 1,{
      df[i,4] <- ""}, 
      ifelse(length(lastrow)==1, {
        prevval <- substr(df[i,c("inf")],1,38)
        ifelse(v=="unconst.",{
          df[i,4] <- eval(parse(text = paste(v,'lastrow', sep= "" )))
          df[i,c("inf")] <- prevval},{
            df[i,4] <- eval(parse(text = paste(v,'lastrow', sep= "" )))
            df[i,c("inf")] <- paste(prevval,"\\\\midrule",sep="")})},
        next))
  }
}

#cost rotate
df$cost<- factor(df$cost, levels = c("current", "future"))
currentlastrow <- paste("{\\\\multirow{-4}{*}{{2016}}}", sep = "")
futurelastrow <- paste("{\\\\multirow{-4}{*}{{2045}}}", sep = "")
df$cost <- factor(df$cost, levels = c(levels(df$cost), ""))
df$cost <- factor(df$cost, levels = c(levels(df$cost), currentlastrow))
df$cost <- factor(df$cost, levels = c(levels(df$cost), futurelastrow))
len <- nrow(df)
for (v in c("current","future")){
  for (i in seq(1,len,1)){
    #find the last row to replace with multicolumn
    firstrow <- as.numeric(which(df[i,5]==v & df[i,5] == df[i+1,5]))
    lastrow <- as.numeric(which(df[i,5]==v & (df[i,5] != df[i+1,5]| df[i,5] == df[len,5])))
    ifelse(length(firstrow)== 1,{
      df[i,5] <- ""}, 
      ifelse(length(lastrow)==1, {
        prevval <- substr(df[i,c("inf")],1,38)
        costlastrow <- paste(v,'lastrow', sep= "" )
        ifelse(costlastrow=="currentlastrow",{
          df[i,5] <- eval(parse(text = costlastrow))
          prevval <- df[i,c("inf")]
          df[i,c("inf")] <- paste(prevval,"\\\\cline{3-9}",sep="")
        },
        {df[i,5] <- eval(parse(text = costlastrow))}
        )
      },
      next))
  }
}
for (ev in c("2016","half","full")) {
  for (load in c("2045","2007")) {
    if ((  thetaval==0.1 & ev=="half" & load==2007  |
           thetaval==0.1 & ev=="2016" & load==2045  |
           thetaval==0.1 & ev=="full" & load==2045  |
           thetaval==0.1 & ev=="half" & load==2045  |
           thetaval==0.5 & ev=="half" & load==2045  | 
           thetaval==2 & ev=="half" & load==2045 )){
      
      num1 <- ifelse(thetaval==0.1 & ev=="half" & load==2007,"S2",
                     ifelse(thetaval==0.1 & ev=="2016" & load==2045,"S3",
                            ifelse(thetaval==0.1 & ev=="full" & load==2045,"S4",
                                   ifelse(thetaval==0.5 & ev=="half" & load==2045,"S5",
                                          ifelse(thetaval==2 & ev=="half" & load==2045,"S6","S1")))))
      
      df0 <- df[which(df$theta==thetaval & df$ev==ev & df$load==load),4:19]
      colnames(df0)<- c("Policy Objective","Demand Flexibility","Cost Assumptions","Pricing","Clean Shares",
                        "Price(mean) in $/MWh","Quantity(mean) in MWh","Price(SD) in $/MWh","Delta CS (%)",
                        "Delta EVChargeCost (%)", "Delta PS (%)","Delta TS (%)",
                        "Delta CS Highflex (%)","Delta CS Midflex (%)$",
                        "Delta CS Inflex (%)","Delta TS from Dynamic pricing (%)")
      latextab <- xtable(df0,caption=paste("Summary",ev,load,thetaval,"Relative to Baseline Expenditure"),sep="")
      print(latextab,
            only.contents=TRUE,
            type="latex", hline.after = NULL,
            include.rownames=FALSE, include.colnames=FALSE,
            sanitize.text.function = function(x){x},
            file=paste("tables/Tab",num1,"_summary_ev=",ev,"_load=",load,"_theta=",thetaval,"withbaselineexpenditure.tex"))
    } else {
      print("skip")
    }
  }
}

################################################################################
##main table#
# Smaller table with below differences compared to the main table in the appendix:
# - eliminate unconstrained
# 
# 9 columns:  
# - labels + Mean P, Mean Q, SD P, Delta TS, Delta TS RTP.
#
# Another table that focuses on distribution, columns 9, 11, 13-15.
################################################################################
#Main table 1
df <- joinall[which(joinall$scen!=2),]
df <- df[which(df$category!="unconst."),]
###############################################################################

df$pricing<- factor(df$pricing, levels = c("flat", "dynamic"))
df <- df[order(df[,1],df[,2],df[,3],df[,4],df[,5],df[,6],rank(df[,7])), ]

df$Renewable_Share <- format(round(as.numeric(df$Renewable_Share),0), nsmall=0, big.mark=",") 
for (i in c("tots","dyntot")){
  df[[i]] <- format(round(as.numeric(df[[i]]),1), nsmall=1, big.mark=",") 
}
for (i in c("mean","meanq","sd")){
  df[[i]] <- format(round(as.numeric(df[[i]]),0), nsmall=0, big.mark=",") 
}

#scen only one 
df <- transform(df, scen=ifelse(df$scen==1,"Optimistic", ifelse(df$scen==2,"Moderate","Pessimistic")))
df$scen<- factor(df$scen, levels = c("Optimistic", "Moderate","Pessimistic"))
scen1lastrow <- c("\\\\multirow{-2}{*}{Optimistic}")
scen2lastrow <- c("\\\\multirow{-2}{*}{Moderate}")
scen3lastrow <- c("\\\\multirow{-2}{*}{Pessimistic}")
df$scen <- factor(df$scen, levels = c(levels(df$scen), ""))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen1lastrow))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen2lastrow))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen3lastrow))

len <- nrow(df)
for (i in seq(1,len,1)){
  #find the last row to replace with multicolumn
  firstrow <- as.numeric(which(df[i,6] == df[i+1,6]))
  lastrow <- as.numeric(which(df[i,6] != df[i+1,6] | df[i,6] == df[len,6]))
  ifelse(length(firstrow)== 1,{
    df[i,6] <- ""}, 
    ifelse(length(lastrow)==1, {
      df[i,6] <- ifelse(df[i,6]=="Optimistic",scen1lastrow,ifelse(df[i,6]=="Moderate",scen2lastrow,scen3lastrow))},
      next))
}

##baseline 0
replacetext <- c("\\\\multicolumn{1}{c}{Baseline}\\\\tabularnewline")
ftx <- 0
index <- which(as.numeric(df$tots)==ftx 
               & as.numeric(df$dyntot)==ftx, arr.ind=TRUE)
df[index,c("tots")] <- replacetext

##dynprice 
ftx <- 0
index0 <- which(as.numeric(df$dyntot)==ftx, arr.ind=TRUE)
index1 <- index0+1
dynline <- paste("\\\\multirow{-2}{*}{",df[index1,19],"}\\\\tabularnewline", sep = "")
df[index1,c("dyntot")] <- paste(dynline,"\\\\cline{5-9}", sep = "")
df[index0,c("dyntot")] <- ""
df$pricing <-  sapply(df$pricing, str_to_title)

#100%renewable rotate
rpslastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering  100\\\\% Clean}}}")
fossillastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering Fossil}}}")
unconst.lastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering Unconstrained}}}")
df$category<- factor(df$category, levels = c("fossil", "rps","unconst."))
df$category <- factor(df$category, levels = c(levels(df$category), ""))
df$category <- factor(df$category, levels = c(levels(df$category), fossillastrow))
df$category <- factor(df$category, levels = c(levels(df$category), rpslastrow))
len <- nrow(df)
for (v in c("fossil", "rps")){
  for (i in seq(1,len,1)){
    #find the last row to replace with multicolumn
    firstrow <- as.numeric(which(df[i+1,4] == v & df[i,4] == v))
    lastrow <- as.numeric(which((df[i+1,4] != v | df[len,4] == v) & df[i,4] == v))
    ifelse(length(firstrow)== 1,{
      df[i,4] <- ""}, 
      ifelse(length(lastrow)==1, {
        prevval <- substr(df[i,c("inf")],1,38)
        ifelse(v=="unconst.",{
          df[i,4] <- eval(parse(text = paste(v,'lastrow', sep= "" )))
          df[i,c("inf")] <- prevval},{
            df[i,4] <- eval(parse(text = paste(v,'lastrow', sep= "" )))
            df[i,c("inf")] <- paste(prevval,"\\\\midrule",sep="")})},
        next))
  }
}

#cost rotate
df$cost<- factor(df$cost, levels = c("current", "future"))
currentlastrow <- paste("{\\\\multirow{-4}{*}{{2016}}}", sep = "")
futurelastrow <- paste("{\\\\multirow{-4}{*}{{2045}}}", sep = "")
df$cost <- factor(df$cost, levels = c(levels(df$cost), ""))
df$cost <- factor(df$cost, levels = c(levels(df$cost), currentlastrow))
df$cost <- factor(df$cost, levels = c(levels(df$cost), futurelastrow))
len <- nrow(df)
for (v in c("current","future")){
  for (i in seq(1,len,1)){
    #find the last row to replace with multicolumn
    firstrow <- as.numeric(which(df[i,5]==v & df[i,5] == df[i+1,5]))
    lastrow <- as.numeric(which(df[i,5]==v & (df[i,5] != df[i+1,5]| df[i,5] == df[len,5])))
    ifelse(length(firstrow)== 1,{
      df[i,5] <- ""}, 
      ifelse(length(lastrow)==1, {
        prevval <- substr(df[i,c("inf")],1,38)
        costlastrow <- paste(v,'lastrow', sep= "" )
        ifelse(costlastrow=="currentlastrow",{
          df[i,5] <- eval(parse(text = costlastrow))
          prevval <- df[i,c("inf")]
          df[i,c("inf")] <- paste(prevval,"\\\\cline{3-9}",sep="")
        },
        {df[i,5] <- eval(parse(text = costlastrow))}
        )
      },
      next))
  }
}
for (ev in c("half")) {
  for (load in c("2045")) {
    if ((load=="2007" & ev!="half" | thetaval!=0.1)){
      print("skip")
    } else {
      df0 <- df[which(df$theta==thetaval & df$ev==ev & df$load==load),c(4,5,6,7,9,10,11,15,19)]
      colnames(df0)<- c("Policy Objective","Demand Flexibility","Cost Assumptions","Pricing",
                        "Price(mean) in $/MWh","Quantity(mean) in MWh","Price(SD) in $/MWh",
                        "Delta TS (%)",
                        "Delta TS from Dynamic pricing (%)")
      latextab <- xtable(df0,caption=paste("Summary",ev,load,thetaval,"Relative to Baseline Expenditure"),sep="")
      #align(latextab) <- "p{0cm}p{3cm}p{2.5cm}p{2.5cm}p{2cm}p{0.6cm}p{0.6cm}p{0.6cm}p{0.6cm}p{0.6cm}p{0.6cm}" 
      #latextab <-  str_replace_all(latextab, pattern = "+\\\\+", " ")
      print(latextab,
            only.contents=TRUE,
            type="latex", hline.after = NULL,
            include.rownames=FALSE, include.colnames=FALSE,
            sanitize.text.function = function(x){x},
            file=paste("tables/Tab_5a_summary_ev=",ev,"_load=",load,"_theta=",thetaval,"withbaselineexpenditure.tex"))
    }
  }
}
################################################################################
#Main table 2
df <- joinall[which(joinall$scen!=2),]
df <- df[which(df$category!="unconst."),]
###############################################################################

df$pricing<- factor(df$pricing, levels = c("flat", "dynamic"))
df <- df[order(df[,1],df[,2],df[,3],df[,4],df[,5],df[,6],rank(df[,7])), ]

df$Renewable_Share <- format(round(as.numeric(df$Renewable_Share),0), nsmall=0, big.mark=",") 
for (i in c("tots","dyntot")){
  df[[i]] <- format(round(as.numeric(df[[i]]),1), nsmall=1, big.mark=",") 
}
for (i in c("mean","meanq","sd")){
  df[[i]] <- format(round(as.numeric(df[[i]]),0), nsmall=0, big.mark=",") 
}

#scen only one 
df <- transform(df, scen=ifelse(df$scen==1,"Optimistic", ifelse(df$scen==2,"Moderate","Pessimistic")))
df$scen<- factor(df$scen, levels = c("Optimistic", "Moderate","Pessimistic"))
scen1lastrow <- c("\\\\multirow{-2}{*}{Optimistic}")
scen2lastrow <- c("\\\\multirow{-2}{*}{Moderate}")
scen3lastrow <- c("\\\\multirow{-2}{*}{Pessimistic}")
df$scen <- factor(df$scen, levels = c(levels(df$scen), ""))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen1lastrow))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen2lastrow))
df$scen <- factor(df$scen, levels = c(levels(df$scen), scen3lastrow))

len <- nrow(df)
for (i in seq(1,len,1)){
  #find the last row to replace with multicolumn
  firstrow <- as.numeric(which(df[i,6] == df[i+1,6]))
  lastrow <- as.numeric(which(df[i,6] != df[i+1,6] | df[i,6] == df[len,6]))
  ifelse(length(firstrow)== 1,{
    df[i,6] <- ""}, 
    ifelse(length(lastrow)==1, {
      df[i,6] <- ifelse(df[i,6]=="Optimistic",scen1lastrow,ifelse(df[i,6]=="Moderate",scen2lastrow,scen3lastrow))},
      next))
}

##baseline 0
replacetext <- c("\\\\multicolumn{4}{c}{--------------------- B  a  s  e  l  i  n  e  ---------------------}\\\\tabularnewline")
ftx <- 0
index <- which(as.numeric(df$dcs)==ftx  & as.numeric(df$dps)==ftx  
               & as.numeric(df$high)==ftx & as.numeric(df$mid)==ftx 
               & as.numeric(df$inf)==ftx, arr.ind=TRUE)
df[index,c("dcs")] <- replacetext
df[index,c("dps")] <- ""
df[index,c("high")] <- ""
df[index,c("mid")] <- ""
df[index,c("inf")] <- ""

##dynprice 
ftx <- 0
index0 <- which(as.numeric(df$dyntot)==ftx, arr.ind=TRUE)
index1 <- index0+1
dynline <- paste(df[index1,18],"\\\\tabularnewline", sep = "")
df[index1,c("inf")] <- paste(dynline,"\\\\cline{5-9}", sep = "")
df$pricing <-  sapply(df$pricing, str_to_title)

#100%renewable rotate
rpslastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering  100\\\\% Clean}}}")
fossillastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering Fossil}}}")
unconst.lastrow <- c("{\\\\multirow{-8}{*}{{\\\\centering Unconstrained}}}")
df$category<- factor(df$category, levels = c("fossil", "rps","unconst."))
df$category <- factor(df$category, levels = c(levels(df$category), ""))
df$category <- factor(df$category, levels = c(levels(df$category), fossillastrow))
df$category <- factor(df$category, levels = c(levels(df$category), rpslastrow))
len <- nrow(df)
for (v in c("fossil", "rps")){
  for (i in seq(1,len,1)){
    #find the last row to replace with multicolumn
    firstrow <- as.numeric(which(df[i+1,4] == v & df[i,4] == v))
    lastrow <- as.numeric(which((df[i+1,4] != v | df[len,4] == v) & df[i,4] == v))
    ifelse(length(firstrow)== 1,{
      df[i,4] <- ""}, 
      ifelse(length(lastrow)==1, {
        prevval <- substr(df[i,c("inf")],1,38)
        ifelse(v=="unconst.",{
          df[i,4] <- eval(parse(text = paste(v,'lastrow', sep= "" )))
          df[i,c("inf")] <- prevval},{
            df[i,4] <- eval(parse(text = paste(v,'lastrow', sep= "" )))
            df[i,c("inf")] <- paste(prevval,"\\\\midrule",sep="")})},
        next))
  }
}

#cost rotate
df$cost<- factor(df$cost, levels = c("current", "future"))
currentlastrow <- paste("{\\\\multirow{-4}{*}{{2016}}}", sep = "")
futurelastrow <- paste("{\\\\multirow{-4}{*}{{2045}}}", sep = "")
df$cost <- factor(df$cost, levels = c(levels(df$cost), ""))
df$cost <- factor(df$cost, levels = c(levels(df$cost), currentlastrow))
df$cost <- factor(df$cost, levels = c(levels(df$cost), futurelastrow))
len <- nrow(df)
for (v in c("current","future")){
  for (i in seq(1,len,1)){
    #find the last row to replace with multicolumn
    firstrow <- as.numeric(which(df[i,5]==v & df[i,5] == df[i+1,5]))
    lastrow <- as.numeric(which(df[i,5]==v & (df[i,5] != df[i+1,5]| df[i,5] == df[len,5])))
    ifelse(length(firstrow)== 1,{
      df[i,5] <- ""}, 
      ifelse(length(lastrow)==1, {
        prevval <- substr(df[i,c("inf")],1,38)
        costlastrow <- paste(v,'lastrow', sep= "" )
        ifelse(costlastrow=="currentlastrow",{
          df[i,5] <- eval(parse(text = costlastrow))
          prevval <- df[i,c("inf")]
          df[i,c("inf")] <- paste(prevval,"\\\\cline{3-9}",sep="")
        },
        {df[i,5] <- eval(parse(text = costlastrow))}
        )
      },
      next))
  }
}
for (ev in c("half")) {
  for (load in c("2045")) {
    if ((load=="2007" & ev!="half" | thetaval!=0.1)){
      print("skip")
    } else {
      df0 <- df[which(df$theta==thetaval & df$ev==ev & df$load==load),c(4,5,6,7,12,14,16,17,18)]
      colnames(df0)<- c("Policy Objective","Demand Flexibility","Cost Assumptions","Pricing",
                        "Delta CS (%)", "Delta PS (%)",
                        "Delta CS Highflex (%)","Delta CS Midflex (%)$",
                        "Delta CS Inflex (%)")
      latextab <- xtable(df0,caption=paste("Summary",ev,load,thetaval,"Distribution"),sep="")
      print(latextab,
            only.contents=TRUE,
            type="latex", hline.after = NULL,
            include.rownames=FALSE, include.colnames=FALSE,
            sanitize.text.function = function(x){x},
            file=paste("tables/Tab_5b_summary_ev=",ev,"_load=",load,"_theta=",thetaval,"withbaselineexpenditure.tex"))
    }
  }
}

# R-code: Analyses r-scripts
# manuscrpt titled "Plant-insect chemical communication in ecological communities: an information theory perspective"
# Authors: Pengjuan Zu???, Reinaldo Garc??ia-Garc??ia, Meredith Schuman, Serguei Saavedra, Carlos J. Meli??an
# correspondence to:  Pengjuan.zu@gmail.com
# cite this paper: 

#########
library(poweRlaw)
library(tidyverse)

#####
### load VOC frequency and rank data for the four community papers and the review paper (in different categories)
df_VOC_all <- readRDS("df_combine_VOC_Freq.RDATA")
df_VOC_all[1:3, ]
Target_all <- levels(df_VOC_all$Target)
Category_all <- levels(df_VOC_all$Category)

### different distribution methods in "poweRlaw" package
# we select two to test
Fun_test <- function(Method_sel){
  E1 <- displ #discrete power law
  E2 <- dislnorm # discrete log normal
  E <- list(E1,E2)
  return(E[[Method_sel]])
}

method_all <- c("pl", "lnorm")

###############################################
### (I): plot
par(mfrow = c(4,4), mar = c(3, 3, 3, 3))
for (i in 1:length(Category_all)) {
  # select dataset
  Category_sel <- Category_all[i]
  
  df_review <- df_VOC_all %>%
    filter(Category == Category_sel)
  df_review[1:3,]
  Target_sel <- unique(df_review$Target)
  print(paste(Target_sel, Category_sel, sep = "; "))
  
  # extract Freq
  df_test <- df_review$Freq
  
  ######## test with power law and plot the best-of-fit lines
  ## functions in "poweRlaw" package
  m_pl <- displ$new(df_test)
  est_pl <- estimate_xmin(m_pl)
  m_pl$setXmin(est_pl)
  
  plot(m_pl, pch = 16, cex = 1.5, col =  "cadetblue3",  main = paste(Target_sel, Category_sel, sep = "; "),  
       panel.first = grid(col = "grey80") ,
       xlab = "", ylab = "")
  #   lines(m_pl, col = "black", lwd = 1.5) ## add power law fitting line
  ## ! note that scale for x- and y-axis increase in log-scale
}


###############################################
##### (II) FULL range power law: test whether whole dataset follow power law distribution: x_min set to be 1
output_all_min1 <- NULL
for (i in 1:length(Category_all)) {
  # select dataset
  Category_sel <- Category_all[i]
  
  df_review <- df_VOC_all %>%
    filter(Category == Category_sel)
  Target_sel <- unique(df_review$Target)
  print(paste(Target_sel, Category_sel, sep = "; "))
  
  # extract Freq
  df_test <- df_review$Freq
  
  # (II)-1. power law
  m_pl <- displ$new(df_test)
  m_pl$setXmin(1) # set Xmin to be 1
  est_pl <- estimate_pars(m_pl)
  m_pl$setPars(est_pl)
  
  x_min <- 1
  Alpha <- est_pl$pars
  
  ### (II)-2. Calculate the goodness-of-fit between the data and the power law.
  ### if p-value is greater than 0.1, the power law is a plausible hypothesis
  bs_p = bootstrap_p(m_pl, no_of_sims = 100, threads = 1) # test with small simulations
  # bs_p = bootstrap_p(m_pl, no_of_sims = 10000, threads = 4) # we did 10'000 simulations
  
  p_value <- bs_p$p
  gof <- bs_p$gof
  # plot(bs_p)
  
  output <- cbind.data.frame(Target_sel, Category_sel, "Power law", x_min, Alpha, p_value, gof)
  names(output) <- c("Target", "Category", "Test_name", "X_min", "alpha", "p_value", "gof")
  output_all_min1 <- rbind.data.frame(output_all_min1, output)
}


###############################################
##### (III) cut-off power law: tests with best estimated x_min
output_all <- NULL
# for each of the 16 cases
for (i in 1:length(Category_all)) {
  # select dataset
  Category_sel <- Category_all[i]
  
  df_review <- df_VOC_all %>%
    filter(Category == Category_sel)
  Target_sel <- unique(df_review$Target)
  print(paste(Target_sel, Category_sel, sep = "; "))
  
  # extract Freq
  df_test <- df_review$Freq
  
  ### test
  ##  estimate parameters x_min (cut-off) and alpha (exponent) of the power law model and other models (log normal)
  for(m in 1:length(method_all)){
    M_name <- paste0("m_", method_all[m])
    print(M_name)
    m_test <- Fun_test(m)$new(df_test)
    
    ### (III)-1. estimate
    # estimate x_min
    est_pl <- estimate_xmin(m_test)
    # test
    m_test$setXmin(est_pl)
    x_min <- m_test$xmin
    Alpha <- m_test$pars
    
    ### (III)-2. Calculate the goodness-of-fit between the data and the power law.
    ### if p-value is greater than 0.1, the power law is a plausible hypothesis
    bs_p = bootstrap_p(m_test, no_of_sims = 100, threads = 1) # test with small simulations
    #  bs_p = bootstrap_p(m_test, no_of_sims = 10000, threads = 4) # we did 10'000 simulations
    p_value <- bs_p$p
    gof <- bs_p$gof
    # plot(bs_p)
    
    ### write output
    output <- cbind.data.frame(Target_sel, Category_sel, M_name, x_min, Alpha, p_value, gof)
    names(output) <- c("Target", "Category", "Test_name", "X_min", "alpha", "p_value", "gof")
    output_all <- rbind.data.frame(output_all, output)
  }
}


###############################################
## (IV) detail look at the power law 
## select these ones that have small x_min and compare with other methods (by setting same x_min)
output_all %>% 
  filter(Test_name %in% "m_pl") %>%
  filter(p_value> 0.1)
m_pl_compare <- c("PV_k", "NS", "overall", "insect", "Orchidaceae")

### Compare the power law with alternative hypotheses (log normal) via a likelihood ratio test
# To compare two distributions, each distribution must have the same lower threshold. 
# So we first set the log normal distribution object to have the same xmin as the power law object
compare_output_all <- NULL
for(i in 1:length(m_pl_compare)){
  cat_sel <- m_pl_compare[i]
  print(paste0("Category selected = ", cat_sel))
  df_set <- df_VOC_all %>%
    filter(Category %in% cat_sel)
  df_test <- df_set$Freq
  
  # power law
  m_pl <- displ$new(df_test)
  est_pl <- estimate_xmin(m_pl)
  m_pl$setXmin(est_pl)
  
  est_pl <- estimate_pars(m_pl)
  m_pl$setPars(est_pl)
  
  # log normal
  m_ln <- dislnorm$new(df_test)
  # set Xmin to be the same as power law
  m_ln$setXmin(m_pl$getXmin())
  est_ln = estimate_pars(m_ln)
  m_ln$setPars(est_ln)
  
  # compare
  comp = compare_distributions(m_pl, m_ln)
  compare_distributions(m_pl, m_ln)$p_one_sided
  compare_distributions(m_ln, m_pl)$p_one_sided
  comp$p_two_sided
  comp$test_statistic
  compare_output <- cbind.data.frame("Category" = cat_sel ,"test_statistic" =comp$test_statistic,
                                     "p_two_sided" = comp$p_two_sided, "p_one_sided" = comp$p_one_sided,
                                     "p_other_sided" = compare_distributions(m_ln, m_pl)$p_one_sided)
  compare_output_all <- rbind.data.frame(compare_output_all, compare_output)
}


###### SIR GROUPS + DISEASE MODEL ###### 

# state variables: 
# g = mean number of uninfected individuals in uninfected groups
# s = mean number of uninfected individuals in infected groups
# i = mean number of infected individuals in infected groups
# r = mean number of recovered individuals in infected groups
# sR = mean number of uninfected individuals in recovered groups
# rR = mean number of recovered individuals in recovered groups
# G = number of uninfected groups
# I = number of infected groups
# R = number of recovered groups

# group parameters:
# Bf = number of breeding females
# b = birth 
# bf (just b in text) = number of offspring/breeding female
# d = death
# f = fission
# c = fusion
# a = aggression/attack (intraspecific mortality)
# e = group extinction rate

# disease parameters:
# beta = transmission rate between G/R and I groups
# betaW = within group transmission rate
# sigma = recovery rate
# alpha = disease-induced mortality
# phi = rate that immunity wanes

GrpDynSIR <- function(t, x, params){

  ## create vectors to hold states through time
  G <- x[1]
  I <- x[2]
  g <- x[3]
  s <- x[4]
  i <- x[5]
  r <- x[6]
  R <- x[7]
  sR <- x[8]
  rR <- x[9]

  ## model
  with(as.list(params),{

    ## equation 13 in text
    dg <- Bf*bf - (d+a)*g - (g*g*f)/((gf+g)*(d/(c+d))) - e*g*(G+I+R)
    ## equation 14 in text
    ds <- Bf*bf*((s+r)/(0.001+s+r)) - (d+a)*s + phi*r - s*f*(s+i+r)/((gf+s+i+r)*(d/(c+d))) - e*s*(G+I+R) - (betaW*s*i)/(1+s+i+r)
    ## equation 15 in text
    di <- (betaW*s*i)/(1+s+i+r) - (d+sigma+alpha+a)*i - e*i*(G+I+R)
    ## equation 16 in text
    dr <- sigma*i - (d+phi+a)*r - r*f*(s+i+r)/((gf+s+i+r)*(d/(c+d))) - e*r*(G+I+R)
    ## equation 17 in text
    dsR <- Bf*bf*((sR+rR)/(0.001+sR+rR)) - (d+a)*sR + phi*rR - ((sR+rR)*f*sR)/((gf+sR+rR)*(d/(c+d))) - e*sR*(G+I+R)
    ## equation 18 in text
    drR <- sigma*i/(1+i) - (d+phi+a)*rR - ((sR+rR)*rR*f)/((gf+sR+rR)*(d/(c+d))) - e*rR*(G+I+R)
    
    ## equation 19 in text
    dG <- (G*c*f*g*g)/((gf+g)*(d+c)) - a*G - e*G*(G+I+R) + R*phi/(1+rR) - beta*I*G/(G+I+R)
    ## equation 20 in text
    dI <- beta*I*(G+(R*(sR/(sR+rR))))/(G+I+R) - a*I - e*I*(G+I+R) - I*((alpha+sigma)/(1+i))
    ## equation 21 in text
    dR <- I*sigma/(1+i) + (c*f*R*((sR+rR)^2))/((gf+sR+rR)*(c+d)) - a*R - e*R*(G+I+R) - R*phi/(1+rR) - beta*I*(R*(sR/(sR+rR)))/(G+I+R)
    
    ## output
    res <- c(dG, dI, dg, ds, di, dr, dR, dsR, drR)
    list(res)
  })
}


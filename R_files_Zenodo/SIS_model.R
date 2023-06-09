###### SIS GROUPS + DISEASE MODEL ###### 

# state variables: 
# g = mean number of uninfected individuals in uninfected groups
# s = mean number of uninfected individuals in infected groups
# i = mean number of infected individuals in infected groups
# G = number of uninfected groups
# I = number of infected groups

# group parameters:
# Bf = number of breeding females
# b = birth 
# bf (just b in text) = number of offspring/breeding female
# d = death
# f = fission
# c = fusion
# a = aggression/attack (intraspecific mortality)
# e = group extinction rate
# A = allee threshold

# disease parameters:
# beta = transmission rate between G and I groups
# betaW = within group transmission rate
# sigma = recovery rate
# alpha = disease-induced mortality


GrpDynSI <- function(t, x, params){
  
  ## create vectors to hold states through time
  G <- x[1]  
  I <- x[2]
  g <- x[3]
  s <- x[4]
  i <- x[5]
  
  ## model
  with(as.list(params),{
    
    ## equation 8 in text
    dg <- Bf*bf*(g/(A+g)) - d*g - (g*g*f)/((gf+g)*(d/(c+d))) - a*g - e*g*(G+I)
    ## equation 9 in text
    ds <- ((s*Bf)/(s+0.001))*b*((s+i)/(A+s+i)) - d*s + sigma*i - ((s+i)*s*f/(gf+(s+i)))*(d/(c+d)) - a*s - e*s*(G+I) - (betaW*s*i)/(1+s+i)
    ## equation 10 in text
    di <- ((betaW*s*i)/(1+s+i)) - (d+alpha+sigma+a)*i - ((s+i)*i*f/((gf+s+i)*(d/(c+d)))) - e*i*(G+I)
    
    ## equation 11 in text
    dG <- (c*f*g*g*G)/((gf+g)*(d+c)) - a*G - e*G*(G+I) + I*sigma/(1+i) - beta*I*G/(G+I)
    ## equation 12 in text
    dI <- beta*I*G/(G+I) + (f*c*I*(s+i)^2)/((gf+s+i)*(d+c)) - a*I - e*I*(G+I) - I*((sigma+alpha)/(1+i))
    
    ## output  
    res <- c(dG, dI, dg, ds, di)
    list(res)
  })
}


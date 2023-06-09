library(foreign)
afdropidx <- c("Congo, Dem. Rep." , "Cape Verde" , "Kenya" , "Sao Tome and Principe" , "Seychelles") # ne 'estimate' data on these countries
bigdat <- read.dta(file="~/Dropbox/quantwork/WarDemocAIDS/statsworking/WarDemHIV2007.dta")
bigdat <- subset(bigdat, year >= 1990)
africa <- subset(bigdat, africayes == 1) 
afYest <- africa
for (i in afdropidx) {
    afYest <- subset(afYest, country != i)
}

devgnoaf <- subset(bigdat, africayes == 0 & devgyes == 1)
devdropidx <- c("Afghanistan", "Albania", "Bosnia and Herzegovina","Bulgaria","Bangladesh","Dominica","East Timor","Egypt, Arab Rep.","Grenada","Iraq","Jordan","Korea, Dem. Rep.","Libya","Macedonia, FYR","Maldives","Marshall Islands","Micronesia, Fed. Sts.","Palau","Philippines","Samoa","Solomon Islands","Sri Lanka","St. Kitts and Nevis","St. Lucia","St. Vincent and the Grenadines","Syrian Arab Republic","Tonga","Turkey","Turkmenistan","Vanuatu","Venezuela, RB","West Bank and Gaza","Yemen, Rep.")
devgYest <- devgnoaf
for (i in devdropidx){
    devgYest <- subset(devgYest, country != i)
}

library(ggplot2)
#using the africa data
ggplot() +
    geom_line(aes(y = estimate,x = year),data=afYest,alpha = 0.5) +
#    geom_density(aes(x = year),data=afYest,size = 0.25) +
    facet_wrap(facets = ~country) + #
    scale_y_continuous(limits = c(0,30)) +
    labs(title="HIV prevalence estimates, Sub-Saharan Africa, by country") +
    xlab("Year") + ylab("HIV Prevalence")

#selected country plots
# pullcountry <- c("Angola", "Botswana")
p <- ggplot() +
    geom_line(aes(y = estimate,x = year),data=subset(afYest, country=="Angola"),alpha = 0.5) +
    #    geom_density(aes(x = year),data=afYest,size = 0.25) +
#     facet_wrap(facets = ~pullcountry) + #
    scale_y_continuous(limits = c(0,4)) +
    labs(title="HIV prevalence estimates, Angola") +
    xlab("Year") + ylab("HIV Prevalence")
p + geom_abline(intercept=-238.2209838, slope=0.1198142, color="red")

coef1 <- coef(lm(estimate~year, data=subset(afYest, country=="Angola" & year <2002)))
coef2 <- coef(lm(estimate~year, data=subset(afYest, country=="Angola" & year >= 2002)))
p + geom_abline(intercept=coef1[[1]], slope=coef1[[2]], color="green") + geom_abline(intercept=coef2[[1]], slope=coef2[[2]], color="blue")

# Non-African Dev'g Countries
ggplot() +
    geom_line(aes(y = estimate,x = year),data=devgYest) + #,alpha = 0.5) +
#     geom_density(aes(x = estimate),data=devgYest,size = 0.25) +
    facet_wrap(facets = ~country) +
    scale_y_continuous(limits = c(0,4)) +
    labs(title="HIV prevalence estimates, Non-SSA Developing Countries, by country") +
    xlab("Year") + ylab("HIV Prevalence")


## as per the final editing suggestions suggestions, try to split this across two pp.
devtools::install_github("guiastrennec/ggplus")
library(ggplus)

pp2 <- ggplot(devgYest, aes(y = estimate,x = year) ) +
    geom_line() +       #,alpha = 0.5) +
    facet_wrap(facets = ~country) +
    scale_y_continuous(limits = c(0,4)) +
    labs(title="HIV prevalence estimates, Non-SSA Developing Countries, by country") +
    xlab("Year") + ylab("HIV Prevalence")

# Plot on multiple pages (output plot to R/Rstudio)
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('thomasp85/ggforce')
library(ggforce)

for (i in 1:2) {
gg <- ggplot(devgYest, aes(y = estimate,x = year) ) +
    geom_line() +       #,alpha = 0.5) +
    facet_wrap_paginate(facets = ~country, nrow=5, ncol=6, page=i) +
    scale_y_continuous(limits = c(0,4)) +
    labs(title="HIV prevalence estimates, Non-SSA Developing Countries, by country") +
    xlab("Year") + ylab("HIV Prevalence")
    print(gg)
}

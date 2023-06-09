library(foreign)
afdropidx <- c("Congo, Dem. Rep." , "Cape Verde" , "Kenya" , "Sao Tome and Principe" , "Seychelles") # ne 'estimate' data on these countries
bigdat <- read.dta(file="~/Dropbox/quantwork/WarDemocAIDS/statsworking/WarDemHIV2007.dta")
africa <- subset(bigdat, africayes == 1 & year >= 1990) 
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
	geom_histogram(aes(y = ..density..,x = estimate),data=afYest,alpha = 0.5) +
	geom_density(aes(x = estimate),data=afYest,size = 0.25) +
	facet_wrap(facets = ~country) +
	scale_y_continuous(limits = c(0,1.5)) +
	opts(title="Density-Histograms for Sub-Saharan Africa, by country") +
	xlab("HIV Prevalence")
 
# All in one
ggplot() +
    geom_histogram(aes(y = ..density..,x = estimate),data=afYest,alpha = 0.5) +
    geom_density(aes(x = estimate),data=afYest,size = 0.25) +
	opts(title="Density-Histograms for Sub-Saharan Africa") +
	xlab("HIV Prevalence")

# Non-African Dev'g Countries
ggplot() +
    geom_histogram(aes(y = ..density..,x = estimate),data=devgYest,alpha = 0.5) +
    geom_density(aes(x = estimate),data=devgYest,size = 0.25) +
    facet_wrap(facets = ~country) +
    scale_y_continuous(limits = c(0,1.5)) +
	opts(title="Density-Histograms for Non-SSA Developing Countries, by country") +
	xlab("HIV Prevalence")

ggplot() +
    geom_histogram(aes(y = ..density..,x = estimate),data=devgYest,alpha = 0.5) +
    geom_density(aes(x = estimate),data=devgYest,size = 0.25) +
    scale_y_continuous(limits = c(0,1.75)) +
	opts(title="Density-Histograms for Non-SSA Developing Countries, by country") +
	xlab("HIV Prevalence")

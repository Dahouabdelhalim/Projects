#Code for Taylor et al 2017 "Interacting effects of genetic variation for seed dormancy and flowering time on phenology, life history, and fitness of experimental Arabidopsis thaliana populations over multiple generations in the field"

#set working directory for where all data files are
setwd('/Users/mrktaylor531/Desktop/Projects/Demography_cage_expts/Data/dryad deposited data')

library(ggplot2)  # for plotting
library(plyr)     # for summarySE function
library(reshape2) # for converting from long to wide formata
library(scales)   # for setting axes labels
library(knitr)    
library(MCMCglmm) # for GLMM model fitting
library(coda)     # for visualizing Bayesian model fitting
library(foreign)
library(scatterplot3d)
library(psych)

#read in data
rawdata = read.table(file='census_data.txt', header=T)

#reformat data: each row = a date
mdata <- melt(rawdata, id=c("stock_no","block", "cage", "Generation")) 
names(mdata)[5] = 'date' #change date col name
mdata$date = as.Date(mdata$date, "X%m.%d.%y") #put dates into R format
mdata$value = as.numeric(mdata$value) #make census counts numeric class
mdata$logvalue = log10(mdata$value+1) #log+1 transform census counts

#insepct formatted phenotype data
str(mdata)
head(mdata)

#this function gives standard errors for some level of data (block, generation, genotype) to use as error bars
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length(na.omit(xx[[col]])),
                     mean = mean   (na.omit(xx[[col]])),
                     sd   = sd     (na.omit(xx[[col]]))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
str(mdata)

dfc = summarySE(mdata, measurevar='logvalue', groupvars=c('stock_no', 'Generation', 'date'))
names(dfc)[5] = 'value'

############################################################
############################################################
############################################################
dfc = summarySE(mdata, measurevar='logvalue', groupvars=c('stock_no', 'Generation', 'date'))
names(dfc)[5] = 'value'

head(dfc)
unique(dfc$stock_no)
#add strong vs weak DOG & Flowering to these bad boys
lesnames=c('Flowering','Dormancy')
dfc.1 = do.call(rbind, lapply(1:nrow(dfc), function(i){
  #i=4
  if(dfc$stock_no[i]=='Col'){df=data.frame('Early-flowering','Weak-dormancy') 
  names(df)=lesnames
  data.frame(dfc[i,],df)}
  else if(dfc$stock_no[i]=='Kas'){df=data.frame('Late-flowering','Strong-dormancy') 
  names(df)=lesnames
  data.frame(dfc[i,],df)}
  else if(dfc$stock_no[i]=='EF-wD'){df=data.frame('Early-flowering','Weak-dormancy') 
  names(df)=lesnames
  data.frame(dfc[i,],df)}
  else if(dfc$stock_no[i]=='EF-sD'){df=data.frame('Early-flowering','Strong-dormancy') 
  names(df)=lesnames
  data.frame(dfc[i,],df)}
  else if(dfc$stock_no[i]=='LF-sD'){df=data.frame('Late-flowering','Strong-dormancy') 
  names(df)=lesnames
  data.frame(dfc[i,],df)}
  else if(dfc$stock_no[i]=='LF-wD'){df=data.frame('Late-flowering','Weak-dormancy') 
  names(df)=lesnames
  data.frame(dfc[i,],df)}
}))

head(dfc.1)
number_ticks <- function(n) {function(limits) pretty(limits, n)} #function to set number of tick marks on axes of plot
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #set plotting color palette

dfc.2 = dfc.1[dfc.1$Generation=='Cotyledons' |dfc.1$Generation=='Reproductive', ] #subset to germinants and flowering only
ck = dfc.2[dfc.2$stock_no=='Col' | dfc.2$stock_no=='Kas',] #subset to ecotypes Col and Kas only

#add a flowering time or parental variable to facet on
head(dfc.2)
nrow(dfc.2)
dfc.3 = do.call(rbind, lapply(1:nrow(dfc.2), function(i){
  #i=624
  if(dfc.2$stock_no[i]=='Col') {var='parent'} else if(dfc.2$stock_no[i]=='Kas') {var='parent'} else {var=dfc.2$Flowering[i]}
  data.frame(dfc.2[i,],var)
}))
head(dfc.3)

####################################
#GERMINATION ONLY
ck = dfc.3[dfc.3$Generation=='Cotyledons',]
unique(ck$stock_no)
head(ck)
tail(ck)

ggplot() + 
  geom_line(data=ck[!is.na(ck$value),], aes(x = date, y = value, linetype=Dormancy), alpha=1, size=1) +
  geom_point(data=ck, aes(x = date, y = value, shape=Dormancy), size=3, alpha=1) +
  facet_grid(var ~ ., scales = "free_x") +
  geom_errorbar(data=ck, aes(date, ymin=value-se, ymax=value+se), alpha=0.8) + #cots error bar
  #xlab('date') +
  xlab('') +
  #ylab(bquote(~log[10]~ N*'')) +
  ylab('') +
  ggtitle('') +
  theme(plot.title = element_text(size = 0, vjust=0)) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  #theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) +
  #theme(strip.text.y = element_text(size = 20, colour = "black", angle = 0)) +
  scale_x_date(breaks=number_ticks(20),labels = date_format("%b")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust=-0.0005)) +
  scale_y_continuous(breaks=number_ticks(7), limits = c(0, 3)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_shape_manual(values=c(1,19), name="RIL Dormancy\\nlevel") +
  scale_linetype_manual(values=c("dashed","solid"), name="Dormancy level") +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  #theme(legend.position="none") +
  theme(strip.text.x = element_text(size=0, angle=75),
        strip.text.y = element_text(size=0, face="bold"),
        strip.background = element_rect(colour="white", fill="white"))

####################################
#FLOWERING ONLY
ck = dfc.3[dfc.3$Generation=='Reproductive',]
unique(ck$stock_no)

ggplot() + 
  geom_line(data=ck[!is.na(ck$value),], aes(x = date, y = value, linetype=Dormancy), alpha=1, size=1) +
  geom_point(data=ck, aes(x = date, y = value, shape=Dormancy), size=3, alpha=1) +
  facet_grid(var ~ ., scales = "free_x") +
  geom_errorbar(data=ck, aes(date, ymin=value-se, ymax=value+se), alpha=0.8) + #cots error bar
  #xlab('date') +
  xlab('') +
  #ylab(bquote(~log[10]~ N*'')) +
  ylab('') +
  ggtitle('') +
  theme(plot.title = element_text(size = 0, vjust=0)) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  #theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) +
  #theme(strip.text.y = element_text(size = 20, colour = "black", angle = 0)) +
  scale_x_date(breaks=number_ticks(20),labels = date_format("%b")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust=-0.0005)) +
  scale_y_continuous(breaks=number_ticks(7), limits = c(0, 3)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_shape_manual(values=c(1,19), name="RIL Dormancy\\nlevel") +
  scale_linetype_manual(values=c("dashed","solid"), name="Dormancy level") +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme(legend.position="none") +
  theme(strip.text.x = element_text(size=0, angle=75),
        strip.text.y = element_text(size=0, face="bold"),
        strip.background = element_rect(colour="white", fill="white"))

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
#FITNESS ANALYSIS
totalfec = read.table(file='fecundity.txt', header=T) #read in fecundity data
generations = unique(totalfec$generation) #get generations

#fitness in certain seasons
tosum = totalfec

summed = aggregate(tosum[,c(5,6,8)], by=list(tosum$genotype, tosum$block, tosum$generation, tosum$season), FUN=sum, na.rm=T) #sum per cage mass, silique #, and plant count
names(summed)[1] = 'genotype'
names(summed)[2] = 'block'
names(summed)[3] = 'year'
names(summed)[4] = 'season'
summed$avg_mass = summed$mass / summed$plant_count
summed$avg_silique = summed$silique_count / summed$plant_count

#divide fitness graphs by season
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #set plotting color palette

#add flowering and dormancy status to summed
head(summed)
lesnames = c('flowering','dormancy')
agmoss = summed
agmoss

summed.1 = do.call(rbind, lapply(1:nrow(agmoss), function(i){
  #i=4
  if(agmoss$genotype[i]=='Col'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='Kas'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-weakDOG'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-strongDOG'){df=data.frame('EF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-strongDOG'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-weakDOG'){df=data.frame('LF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
}))
summed = summed.1

summed.2 = do.call(rbind, lapply(1:nrow(summed.1), function(i) {
  #i=4
  if(summed.1$genotype[i]=='Col') {geno_name='Col'} else if(summed.1$genotype[i]=='EarlyFlower-strongDOG') {geno_name='EF-sD'} else if(summed.1$genotype[i]=='EarlyFlower-weakDOG') {geno_name='EF-wD'} else if(summed.1$genotype[i]=='LateFlower-strongDOG') {geno_name='LF-sD'} else if(summed.1$genotype[i]=='LateFlower-weakDOG') {geno_name='LF-wD'} else if(summed.1$genotype[i]=='Kas') {geno_name='Kas'}  
  data.frame(summed.1[i,],geno_name)
}))
summed.1 = summed.2
head(summed.1)

#############
#############
#############
#############
#############
#Fig. 4
#plot only for after spring 2012
head(summed.1)
fall2012 = summed.1[summed.1$season=='fall' & summed.1$year==2012,]
f2012 = ggplot(data=fall2012, aes(x = geno_name, y=log10(mass+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(a) fall 2012') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  #theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=24,color='black'), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20), axis.text.y=element_text(colour='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,1.6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
spring2013 = summed.1[summed.1$season=='spring' & summed.1$year==2013,]
s2013 = ggplot(data=spring2013, aes(x = geno_name, y=log10(mass+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(b) spring 2013') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20)) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,1.6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
fall2013 = summed.1[summed.1$season=='fall' & summed.1$year==2013,]
f2013 = ggplot(data=fall2013, aes(x = geno_name, y=log10(mass+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(c) fall 2013') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20)) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,1.6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
spring2014 = summed.1[summed.1$season=='spring' & summed.1$year==2014,]
s2014 = ggplot(data=spring2014, aes(x = geno_name, y=log10(mass+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(d) spring 2014') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20)) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,1.6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

fall2012 = summed.1[summed.1$season=='fall' & summed.1$year==2012,]
f2012.1 = ggplot(data=fall2012, aes(x = geno_name, y=log10(plant_count+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(N+1)')) +
  xlab('') +
  ggtitle('(e) fall 2012') +
  theme(plot.title = element_text(size = 30)) +
  #theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=24,color='black'), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20), axis.text.y=element_text(colour='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,3)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
spring2013 = summed.1[summed.1$season=='spring' & summed.1$year==2013,]
s2013.1 = ggplot(data=spring2013, aes(x = geno_name, y=log10(plant_count+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(f) spring 2013') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20)) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,3)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
fall2013 = summed.1[summed.1$season=='fall' & summed.1$year==2013,]
f2013.1 = ggplot(data=fall2013, aes(x = geno_name, y=log10(plant_count+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(g) fall 2013') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20)) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,3)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
spring2014 = summed.1[summed.1$season=='spring' & summed.1$year==2014,]
s2014.1 = ggplot(data=spring2014, aes(x = geno_name, y=log10(plant_count+1))) + 
  geom_boxplot(size=2) +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  xlab('') +
  ggtitle('(h) spring 2014') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 30)) +
  theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=20)) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(10),limits=c(0,3)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")
multiplot(f2012, f2012.1, s2013, s2013.1, f2013, f2013.1, s2014, s2014.1, cols=4)

#############
#############
#############
#############
#############
#Fig. 3

founders = summed.1[summed.1$plant_count <= 5 & summed.1$year==2012 & summed.1$season=='spring',]

found = ggplot(data=founders, aes(x = geno_name, y=log10(mass+1))) + 
  geom_boxplot(size=2, color='black') +
  #facet_grid(season ~ year) +
  xlab('genotype') +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  ggtitle('(a) founder fitness') +
  #scale_shape_discrete(name="genotype") +
  #scale_linetype_discrete(name="genotype") +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 38)) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=23), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=36), axis.text.y = element_text(angle=0, hjust=1, color='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(4)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  theme(panel.background = element_rect(colour = "black")) +
  scale_linetype_manual(values=c(5,1), name="dormancy class") +
  theme(legend.position="none")

#get only spring2012 flowering plants 
spring2012 = summed.1[summed.1$plant_count > 4 & summed.1$year==2012 & summed.1$season=='spring',]
dumrows = read.table(file='/Users/mrktaylor531/Desktop/Projects/Demography_cage_expts/Data/Fecundity_data/dummy_rows.txt', header=T)
spring2012 = rbind(spring2012, dumrows)

spring2012n = ggplot(data=spring2012, aes(x = geno_name, y=log10(plant_count+1))) + 
  geom_boxplot(size=2, color='black') +
  ylab(bquote(~log[10]~'(N+1)')) +
  ggtitle('(b) post-founder spring, 2012 N') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 38)) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=23), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=36), axis.text.y = element_text(angle=0, hjust=1, color='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

spring2012fit = ggplot(data=spring2012, aes(x = geno_name, y=log10(mass+1))) + 
  geom_boxplot(size=2, color='black') +
  ylab(bquote(~log[10]~'(biomass[g]+1)')) +
  ggtitle('(c) post-founder spring, 2012 fitness') +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 38)) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=23), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=36), axis.text.y = element_text(angle=0, hjust=1, color='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

multiplot(found, spring2012n, spring2012fit, cols=1)
####################################################################################################################
####################################################################################################################
###############################################################################################################
###############################################################################################################
#fitness proxy : biomass silique regression
#actual individual fitness from the third generation
thirdgen = totalfec[totalfec$generation==2014,]

#do this regression droppping top 2 mass points
head(thirdgen)

nrow(thirdgen)
#SI Fig. 1
ggplot(data=thirdgen, aes(x=log10(mass), y=log10(silique_count))) +
  geom_point(data=thirdgen, aes(shape=genotype), size=2.2) +
  stat_smooth(method=lm, se=F, color='black') +
  xlab(bquote(~log[10]~ 'mass(g)'*'')) +
  ylab(bquote(~log[10]~ 'silique count'*'')) +
  ggtitle('biomass fitness proxy') +
  theme(axis.text=element_text(size=18,color='black')) +
  theme(plot.title = element_text(size = 24)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(legend.text=element_text(size=20)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black'))+
  theme(panel.background = element_rect(colour = "black"))

thirdgen[match(max(thirdgen$mass), thirdgen$mass),] #find highest plant mass
thirdgen[match(max(thirdgen$silique_count), thirdgen$silique_count),] #find highest silique, i guess believable

#make new thirdgen df without NA of Inf values to run lm on
thirdgen.1 = data.frame(log10(thirdgen$silique_count),log10(thirdgen$mass))
thirdgen.1 = na.omit(thirdgen.1)
names(thirdgen.1)[1] = 'silique_count'
names(thirdgen.1)[2] = 'mass'
is.na(thirdgen.1) <- do.call(cbind,lapply(thirdgen.1, is.infinite)) #exclude infinite values

summary(lm(log10(silique_count) ~ log10(mass), na.action=na.omit,data = thirdgen.1))

#########################################################################################################################################
#########################################################################################################################################
###Interspecific and intraspecific competition

agdata = aggregate(totalfec[,c(5,6,8)], by=list(totalfec[,1], totalfec[,3], totalfec[,7], totalfec[,9]), 
                   FUN=sum, na.rm=TRUE)

names(agdata)[1] = 'block'
names(agdata)[2] = 'genotype'
names(agdata)[3] = 'generation'
names(agdata)[4] = 'season'

#read in interspecific competition data (called mossdat since it mainly consisted of moss)
mossdat = read.table(file='competition_data.txt', header=T)

#make rowid's that match between agdata and mossdat
agdata$rowid = paste(agdata$block, agdata$genotype, agdata$generation, sep='_')
mossdat$rowid = paste(mossdat$block, mossdat$genotype, mossdat$generation, sep='_')

#find and rbind matching rows
agmoss = do.call(rbind, lapply(1:nrow(mossdat), function(i) {
  #i=1
  tryCatch({
    data.frame(agdata[agdata$rowid==mossdat$rowid[i],], mossdat$score[i])
  }, error=function(e){})
}))
names(agmoss)[ncol(agmoss)] = 'compscore'
agmoss$massperplant = agmoss$mass / agmoss$plant_count

######################################################################################
#individual graphs for intraspecific competition versus interspecific competition

agmoss.1 = agmoss1[agmoss1$generation != 2012,]
head(agmoss.1)
#intraspecific
intraspp.pop = ggplot(data = agmoss.1, aes(x=log10(plant_count+1), y = log10(mass+1))) +
  geom_point(data = agmoss.1, aes(shape=genotype), size=4) +
  stat_smooth(method=lm, color='black', se = F) +
  ggtitle('') +
  xlab(bquote(~log[10]~ '(density+1)'*'')) +
  ylab(bquote(~log[10]~ 'population fitness: biomass (g)'*'')) +
  theme(axis.text=element_text(size=22,color='black')) +
  ggtitle('POPULATION: intraspecific') +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

attach(agmoss.1)
summary(lm(log10(mass+1) ~ log10(plant_count+1)))

#interspecific competition
interspp.pop = ggplot(data = agmoss.1, aes(x=compscore, y = log10(mass+1))) +
  geom_point(data = agmoss.1, aes(shape=genotype), size=4) +
  stat_smooth(method=lm, color='black', se = FALSE) +
  ggtitle('POPULATION: interspecific') +
  xlab('interspecific competition score') +
  ylab(bquote(~log[10]~ 'population fitness: biomass (g)'*'')) +
  theme(axis.text=element_text(size=22,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24, color='black')) +
  theme(text = element_text(size=24, color='black')) +
  theme(axis.text = element_text(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

agmoss.1.1 = agmoss.1[!is.na(agmoss.1$compscore),]
summary(lm(log10(mass+1) ~ compscore, data = agmoss.1.1))

####################
#plot average individual fitness biomass
lesnames=c('Flowering','Dormancy')
agmoss1 = do.call(rbind, lapply(1:nrow(agmoss), function(i){
  #i=4
  if(agmoss$genotype[i]=='Col'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='Kas'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-weakDOG'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-strongDOG'){df=data.frame('EF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-strongDOG'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-weakDOG'){df=data.frame('LF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
}))
agmoss.1.1 = agmoss1[!is.na(agmoss1$massperplant),]
agmoss.12014 = agmoss.1.1[agmoss.1.1$generation==2014,]

agmosslf = agmoss.12014[agmoss.12014$Flowering=='LF',]
agmossef = agmoss.12014[agmoss.12014$Flowering=='EF',]

#plot EF & LF together
intraspp.ind = ggplot() +
  geom_point(data = agmoss.12014, aes(x=log10(plant_count+1), y = log10(massperplant+1), shape=genotype), size=4) +
  stat_smooth(data = agmossef, aes(x=log10(plant_count+1), y = log10(massperplant+1)), method=lm, color='black', se = F) +
  stat_smooth(data = agmosslf, aes(x=log10(plant_count+1), y = log10(massperplant+1)), method=lm, linetype='dashed', color='black', se = F) +
  xlab(bquote(~log[10]~ '(N+1)'*'')) +
  ylab(bquote(~log[10]~ '(average per plant biomass[g])'*'')) +
  ggtitle('(a) intraspecific competition') +
  theme(plot.title = element_text(size = 32)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24, color='black')) +
  theme(text = element_text(size=28, color='black')) +
  theme(axis.text = element_text(size=28,colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks=number_ticks(10)) +
  scale_x_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.key = element_rect(fill = "white", colour = "white")) +
  theme(legend.position="none")

summary(lm(log10(massperplant+1) ~ log10(plant_count+1) + Flowering + log10(plant_count+1)*Flowering, data = agmoss.12014))

intraspp.ind.ef = ggplot(data = agmossef, aes(x=log10(plant_count+1), y = log10(massperplant+1))) +
  geom_point(data = agmossef, aes(shape=genotype), size=4) +
  stat_smooth(method=lm, color='black', se = F) +
  xlab(bquote(~log[10]~ 'N+1'*'')) +
  ylab(bquote(~log[10]~ 'average per plant biomass(g)'*'')) +
  theme(axis.text=element_text(size=22,color='black')) +
  ggtitle('intraspecific competition') +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24, color='black')) +
  theme(text = element_text(size=24, color='black')) +
  theme(axis.text = element_text(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks=number_ticks(10)) +
  scale_x_continuous(breaks=number_ticks(10)) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  #theme(legend.position="none")
  
  head(agmoss.1.1)
summary(lm(log10(massperplant+1) ~ log10(plant_count+1),data = agmossef))

intraspp.ind.lf = ggplot(data = agmosslf, aes(x=log10(plant_count+1), y = log10(massperplant+1))) +
  geom_point(data = agmosslf, aes(shape=genotype), size=4) +
  stat_smooth(method=lm, color='black', se = F) +
  xlab(bquote(~log[10]~ 'N+1'*'')) +
  ylab(bquote(~log[10]~ 'average per plant biomass(g)'*'')) +
  theme(axis.text=element_text(size=22,color='black')) +
  ggtitle('intraspecific competition') +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24, color='black')) +
  theme(text = element_text(size=24, color='black')) +
  theme(axis.text = element_text(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks=number_ticks(10)) +
  scale_x_continuous(breaks=number_ticks(10)) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  #theme(legend.position="none")
  
  summary(lm(log10(massperplant+1) ~ log10(plant_count+1),data = agmosslf))

multiplot(intraspp.ind.ef, intraspp.ind.lf, cols=2)

#interspecific competition
interspp.ind = ggplot(data = agmoss.1.1, aes(x=compscore, y = log10(massperplant+1))) +
  geom_point(data = agmoss.1, aes(shape=genotype), size=4) +
  stat_smooth(method=lm, color='black', se = FALSE) +
  ggtitle('INDIVIDUAL: interspecific') +
  xlab('interspecific competition score') +
  ylab(bquote(~log[10]~ '(average per plant biomass[g])'*'')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24, color='black')) +
  theme(text = element_text(size=24, color='black')) +
  theme(axis.text = element_text(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  #theme(legend.key = element_rect(fill = "white")) + #uncomment this to get the legend
  theme(legend.position="none")

summary(lm(log10(massperplant+1) ~ compscore, data = agmoss.1.1))
multiplot(intraspp.pop, interspp.pop, intraspp.ind, interspp.ind, cols=2)

##############################################################################################################
##############################################################################################################
#Assymetric competition
totalfec = read.table(file='/Users/mrktaylor531/Desktop/Projects/Demography_cage_expts/Data/Fecundity_data/combined_biomass_fec_2.txt', header=T) #read in fecundity data
third_gen = totalfec[totalfec$generation=='2014',]
third_gen$id = paste(third_gen$block, third_gen$genotype, sep='_')
cages = unique(third_gen$id)

cvs = do.call(rbind, lapply(1:length(cages), function(i) { #this builds a dataframe of all sorts of fun thangs
  #i=1
  block.df = third_gen[third_gen$id==cages[i],]
  cv.biom = sd(block.df$mass) / mean(block.df$mass)
  cv.sil = sd(block.df$silique_count) / mean(block.df$silique_count)
  plant_no = sum(block.df$plant_count)
  sil_sum = sum(block.df$silique_count)
  sil_gm = geometric.mean(block.df$silique_count) #geometric average of siliques per plant per cage
  sil_mean = mean(block.df$silique_count) #arithmatic avg of siliques per plant per cage
  df = data.frame(block.df$generation[1],block.df$genotype[1], plant_no, sil_sum, cv.biom, cv.sil, sil_gm, sil_mean)
  names(df)[1] = 'generation'
  names(df)[2] = 'genotype'
  df
}))

#plot asymmetric competition
asplot = ggplot(data = cvs, aes(x=log10(plant_no), y = cv.biom)) +
  geom_point(data = cvs, aes(x=log10(plant_no), y = cv.biom, shape=genotype), size=4) +
  stat_smooth(method=lm, color='black', se = F)+ #, aes(fill=generation)) +
  ggtitle('(b) asymmetric competition') +
  xlab(bquote(~log[10]~ '(N+1)'*'')) +
  ylab('CV for individual plant biomass(g)') +
  theme(plot.title = element_text(size = 32)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24, color='black')) +
  theme(text = element_text(size=28, color='black')) +
  theme(axis.text = element_text(size=28,colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks=number_ticks(10)) +
  scale_x_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.key = element_rect(fill = "white", colour = "white")) +
  theme(legend.position="none")

summary(lm(cv.biom ~ log10(plant_no), data = cvs))

#Fig. 5
multiplot(intraspp.ind, asplot, cols=2)

############################################################
############################################################
############################################################
#Visualize overdispersion in variables
#biomass
biom = ggplot(agdata, aes(mass)) + 
  geom_histogram(color='black') + 
  ggtitle('total biomass') +
  xlab('mass (g)') +
  theme(axis.text=element_text(size=22,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(7)) +
  theme(panel.background = element_rect(colour = "black")) 

#biomass per plant in each cage
biompp = ggplot(agmoss, aes(massperplant)) + 
  geom_histogram(color='black') + 
  ggtitle('biomass per plant') +
  xlab('average per plant mass (g)') +
  theme(axis.text=element_text(size=22,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(7)) +
  theme(panel.background = element_rect(colour = "black"))

#slique count
sil = ggplot(agdata, aes(silique_count/1000)) + 
  geom_histogram(color='black') + 
  ggtitle('fitness') +
  xlab('average per plant silique count (x1000)') +
  theme(axis.text=element_text(size=22,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(7)) +
  theme(panel.background = element_rect(colour = "black"))

#plant density
density = ggplot(agdata, aes(log10(plant_count+1))) + 
  geom_histogram(color='black') + 
  ggtitle('intraspecific competition') +
  xlab(bquote(~log[10]~ 'N'*'')) +
  theme(axis.text=element_text(size=22,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(7)) +
  theme(panel.background = element_rect(colour = "black"))

#competition overdispersion
comp = ggplot(agmoss, aes(compscore)) + 
  geom_histogram(color='black') + 
  ggtitle('interspecific competition') +
  xlab('competition score') +
  theme(axis.text=element_text(size=22,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(7)) +
  theme(panel.background = element_rect(colour = "black"))

#SI Fig. 2
multiplot(biom, density, comp, cols=1)
##############################################################################################################
################################################################################################################################################
#loop to add the flowering and seed dormancy status to each genotype
lesnames=c('Flowering','Dormancy')
agmoss1 = do.call(rbind, lapply(1:nrow(agmoss), function(i){
  #i=4
  if(agmoss$genotype[i]=='Col'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='Kas'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-weakDOG'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-strongDOG'){df=data.frame('EF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-strongDOG'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-weakDOG'){df=data.frame('LF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
}))

#convert response variables to non-negative integers to allow for poisson-distributed assumptions
#multiply by 1000 and round
agmoss1$mass.int = round(agmoss1$mass * 1000)
agmoss1$massperplant.int = round(agmoss1$massperplant*1000)

#define prior structure
prior2 = list(R = list(V = diag(2)/3, n = 2),
              G = list(G1 = list(V = diag(2)/3, n = 2),
                       G2 = list(V = diag(2)/3, n = 2),
                       G3 = list(V = diag(2)/3, n = 2)))
##################################################################################################################
#population-fitness MCMC model only
#response variable = population fitness only
#exclude the founders
agmoss1 = na.omit(agmoss1) #get rid of cages lacking interspp competition information
agmoss1.1 = agmoss1[agmoss1$plant_count != 4,] #must exclude founding generation from this because we controlled its population size (but not its fruiting)
agmoss1.1$genseas = paste(agmoss1.1$generation, agmoss1.1$season, sep='_')
agmoss1.1.2 = agmoss1.1[agmoss1.1$genseas != '2012_spring',]
m1<-MCMCglmm(mass.int ~ Dormancy + Flowering + Dormancy:Flowering + compscore + block + season + Flowering:season + Dormancy:season + season:Dormancy:Flowering,
             random = ~ block + compscore, rcov = ~ units,
             family = "poisson", nitt = 450000, burnin = 10000,
             thin=25, data = agmoss1.1.2)
#Table 2a
summary(m1)
#SI3.3
plot(m1)

#run GLMM on the plant abundance within a cage
m2<-MCMCglmm(plant_count ~ Dormancy + Flowering + Dormancy:Flowering + compscore + block + season + Flowering:season + Dormancy:season + season:Dormancy:Flowering,
             random = ~ block + compscore, rcov = ~ units,
             family = "poisson", nitt = 1000000, burnin = 10000,
             thin=25, data = agmoss1.1.2)

##Table 2b
summary(m2)
#SI3.3.2
plot(m2)

##################################################################################################################
##############################################################################################################
##########################
#CALCULATE LAMBDA GROWTH RATE
generations = unique(totalfec$generation) #get generations

#fitness in certain seasons
tosum = totalfec

summed = aggregate(tosum[,c(5,6,8)], by=list(tosum$genotype, tosum$block, tosum$generation, tosum$season), FUN=sum, na.rm=T) #sum per cage mass, silique #, and plant count
names(summed)[1] = 'genotype'
names(summed)[2] = 'block'
names(summed)[3] = 'year'
names(summed)[4] = 'season'
summed$massperplant = summed$mass / summed$plant_count
summed$avg_silique = summed$silique_count / summed$plant_count

#add flowering and dormancy status to summed
lesnames = c('flowering','dormancy')
agmoss = summed

#Calculate lambda growth rate
agmoss$silperplant = agmoss$silique_count / agmoss$plant_count
agmoss$cageid = paste(agmoss$block, agmoss$genotype, sep='_')
cageids = unique(agmoss$cageid)

#add up over seasons within each row
aggmoss = aggregate(agmoss[,c("plant_count","silique_count","massperplant")], by=list(agmoss$block,agmoss$genotype,agmoss$year,agmoss$cageid), "sum")
names(aggmoss)[1:4] = c('block','genotype','generation','cageid')
agmoss = aggmoss #rename so I can do the graphing more easily 

#check this whole damn thing because I do not think this is working : no block 5 and also repeating cage_ids each thing
r0 = do.call(rbind, lapply(1:length(cageids), function(i) {
  #i=25
  tryCatch({ #suppress error in names because last generation only has 3 blocks
    third = agmoss[agmoss$cageid == cageids[i] & agmoss$generation==2014,]  
    second = agmoss[agmoss$cageid == cageids[i] & agmoss$generation==2013,]
    first = agmoss[agmoss$cageid == cageids[i] & agmoss$generation==2012,]
    #third
    #second
    r013.totsil = log(third$silique_count+1) / log(first$silique_count+1) #calculate instantaneous rate of increase=r=ln(lambda) , lambda = Nt+1/Nt
    r012.avgmass =  second$massperplant / first$massperplant
    r023.avgmass = third$massperplant / second$massperplant
    r013.avgmass = third$massperplant / first$massperplant
    r013.totplant = log10(third$plant_count+1) / log10(4+1)
    r023.totplant = log10(third$plant_count+1) / log10(second$plant_count+1)
    r012.totplant = log10(second$plant_count+1) / log10(first$plant_count+1)
    data.frame(first[,c(1:4)], r013.totsil,r012.avgmass,r023.avgmass,r013.avgmass,r013.totplant,r023.totplant,r012.totplant)
  }, error=function(e){})
}))
r0
################################################
################################################
#now set up dataframe to do GLMM population growth

#GLMM on r0
#R0 ~ dormancy + flowring + dormancy:flowering + interspecific competition
#DOG and MAF interaction
#loop to add the flowering and seed dormancy status to each genotype
lesnames=c('flowering','dormancy')
head(r0)
agmoss = r0
r0.1 = do.call(rbind, lapply(1:nrow(agmoss), function(i){
  #i=4
  if(agmoss$genotype[i]=='Col'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='Kas'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-weakDOG'){df=data.frame('EF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='EarlyFlower-strongDOG'){df=data.frame('EF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-strongDOG'){df=data.frame('LF','sD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
  else if(agmoss$genotype[i]=='LateFlower-weakDOG'){df=data.frame('LF','wD') 
  names(df)=lesnames
  data.frame(agmoss[i,],df)}
}))

head(r0.1)
#reproductive individuals third year / first year
#visualize lambda

#add the geno_name column
r0.11 = do.call(rbind, lapply(1:nrow(r0.1), function(i) {
  #i=4
  if(r0.1$genotype[i]=='Col') {geno_name='Col'} else if(r0.1$genotype[i]=='EarlyFlower-strongDOG') {geno_name='EF-sD'} else if(r0.1$genotype[i]=='EarlyFlower-weakDOG') {geno_name='EF-wD'} else if(r0.1$genotype[i]=='LateFlower-strongDOG') {geno_name='LF-sD'} else if(r0.1$genotype[i]=='LateFlower-weakDOG') {geno_name='LF-wD'} else if(r0.1$genotype[i]=='Kas') {geno_name='Kas'}  
  data.frame(r0.1[i,],geno_name)
}))


head(r0.11)
#Fig. 6a
adult_lambda = ggplot(r0.11, aes(x = geno_name, y = r013.totplant)) +
  geom_boxplot(size=2) +
  ggtitle('(a)   2012-2014 r (adult-to-adult)') +
  xlab('') +
  ylab(expression(paste("ln(",lambda,")"))) +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 32)) +
  #theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=28,color='black'), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=28), axis.text.y=element_text(colour='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(6),limits=c(0,4)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

#Fig. 6b
#SEED-TO-SEED LAMBDA
head(r0.11)
seed_lambda = ggplot(r0.11, aes(x = geno_name, y = r013.totsil)) +
  geom_boxplot(size=2) +
  ggtitle('(b)   2012-2014 r (seed-to-seed)') +
  xlab('') +
  ylab(expression(paste("ln(",lambda,")"))) +
  theme(text = element_text(size=18)) +
  theme(plot.title = element_text(size = 32)) +
  #theme(axis.text.y = element_blank(),axis.title.y=element_blank()) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(text = element_text(size=28,color='black'), axis.text.x = element_text(angle=0, hjust=0.5, color='black', size=28), axis.text.y=element_text(colour='black')) +
  theme(axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=0)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(limits=c('Col','EF-sD','EF-wD','Kas','LF-sD','LF-wD')) + #chnge order of x-axis lbels
  scale_y_continuous(breaks=number_ticks(6),limits=c(0,4)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(legend.position="none")

#Fig. 6
multiplot(adult_lambda, seed_lambda, cols=1)

head(r0.1)
#add competition
#find and rbind matching rows
head(mossdat)
r0.1$generation = 2014
r0.1$rowid = paste(r0.1$block,r0.1$genotype,r0.1$generation,sep='_')

r0.11 = do.call(rbind, lapply(1:nrow(r0.1), function(i) {
  #i=1
  tryCatch({
    data.frame(r0.1[i,],mossdat[r0.1$rowid[i]==mossdat$rowid,5])
  }, error=function(e){})
}))
r0.11
names(r0.11)[ncol(r0.11)] = 'compscore'

####################################################################################
#GLMM for differences in population growth rate based on reproductive adults
head(r0.11)
m4<-MCMCglmm(r013.totplant ~ genotype + block + compscore,
             random = ~ block, rcov = ~ units,
             family = "gaussian", nitt = 600000, burnin = 10000,
             thin=25, data = r0.11)

#SI Table 2 
summary(m4)
#SI Fig. 3
plot(m4)
####################################################################################
#GLMM for differences in population growth rate based on SEED-to-SEED
head(r0.11)
m5<-MCMCglmm(r013.totsil ~ genotype + block + compscore,
             random = ~ block, rcov = ~ units,
             family = "gaussian", nitt = 600000, burnin = 10000,
             thin=25, data = r0.11)

#SI Table 2 
summary(m5)
#SI Fig. 3
plot(m5)











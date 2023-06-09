# title: "Do bumblebees use (social) cues more readily when the task is difficult?"
# author: "Marco Smolla"
# date: "10 November 2015"

library(ggplot2)
library(gridExtra)
library(knitr)
library(reshape2)
library(cowplot)
library(lme4)

p <- ggplot() +
	theme_bw() +
	xlab('Individual') +
	theme(axis.title.x = element_text(size=15,family='Helvetica'),
				axis.title.y = element_text(size=15, angle=90,family='Helvetica'),
				axis.text.x = element_text(size=13,family='Helvetica'),
				axis.text.y = element_text(size=13,family='Helvetica'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				legend.title = element_text(size=13, face='bold'),
				legend.text = element_text(size=13, family='Helvetica'),
				legend.key = element_blank(),
				plot.title = element_text(face='bold', lineheight=.8, size=10,family='Helvetica'),
				strip.text.x = element_text(face='bold', lineheight=.8, size=13,family='Helvetica'),
				strip.text.y = element_text(face='bold', lineheight=.8, size=13,family='Helvetica')
	)


q <-  ggplot() +
					theme_bw() +
					scale_y_continuous(lim=c(0,1.02), expand=c(0,0) ,breaks=seq(from = 0.2,to = 1.0,by=0.2)) +
					scale_fill_manual(values=c("black", "white")) +
					xlab('') +
					theme(
						panel.border = element_blank(),
						axis.line.x=element_line(),
						axis.line.y=element_line(),
							axis.title.x = element_text(size=8,family='Helvetica'),
							axis.title.y = element_text(size=8, angle=90,family='Helvetica'),
							axis.text.x = element_text(size=6,family='Helvetica'),
							axis.text.y = element_text(size=6,family='Helvetica'),
							panel.grid.major = element_blank(),
							panel.grid.minor = element_blank(),
							legend.title = element_text(size=6, face='bold'),
							legend.text = element_text(size=6, family='Helvetica'),
							legend.key = element_blank(),
							plot.title = element_text(face='bold', lineheight=.8, size=4,family='Helvetica'),
							strip.text.x = element_text(face='bold', lineheight=.8, size=6,family='Helvetica'),
							strip.text.y = element_text(face='bold', lineheight=.8, size=6,family='Helvetica')
						)

# Standard error
se <- function(x,na.rm=F){
	return( sd(x,na.rm = na.rm) / sqrt(length(x)) )
}

# Function to check whether chip landed on is a chip with a cue
llandingsF <- function(l,c){
	tmp <- list()
	for(i in 1:length(l)){
		tmp[[i]] <- l[[i]]%in%c[[i]]
	}
	return(tmp)
}

# Function to extract basic values from the landings: total landings, total landings on cue, first choice, first 4, first 10, total cue prop
data_add_infoF <- function(dat, ll){
	dat$total <- unlist(lapply(ll, length))
	dat$landedOnCue <- unlist(lapply(ll, sum))
	dat$totalCueProp <- unlist(lapply(ll, function(x) return(sum(x)/length(x))))
	dat$first10CueProp <-unlist(lapply(ll, function(x) {ifelse(length(x)<10, yes = sum(x)/length(x), no = sum(x[1:10])/10)}))
	dat$first4CueProp <- unlist(lapply(ll, function(x) {ifelse(length(x)<4 , yes = sum(x)/length(x), no = sum(x[1:4]) /4 )}))
	dat$firstChoice <- unlist(lapply(ll, '[',1))
	return(dat)
}

data_add_infoF2 <- function(dat, l, c, r){
	dat$uniqueFlowers <- unlist(lapply(lapply(l, unique), length))
	dat$uniqueCue <- unlist( lapply(c(1:length(l)), function(i) length(unique(l[[i]][l[[i]]%in%c[[i]]]))) )
	dat$collectedReward <- unlist( lapply(c(1:length(l)), function(i) length(unique(l[[i]][l[[i]]%in%r[[i]]]))) )

	visitsTillReward <- lapply(c(1:length(l)), function(i){
		first <- which(l[[i]]%in%r[[i]])[1]
		rr <- r[[i]][ r[[i]] != l[[i]][first] ]
		second <- which(l[[i]]==rr)[1]
		return(c(first, second))
	})
	dat$visitsTillFirstReward <- unlist(lapply(visitsTillReward, '[', 1))
	dat$visitsTillSecondReward <- unlist(lapply(visitsTillReward, '[', 2))

	return(dat)
}

### Functions for cue training
landingsF <- function(dat){
		lapply(dat$landed, function(x) as.numeric(unlist(strsplit(as.character(x), split=' '))))}

cueF <- function(dat){
		lapply(dat$cue, function(x) as.numeric(unlist(strsplit(as.character(x), split=' '))))}

rewardF <- function(dat){
		lapply(dat$reward, function(x) as.numeric(unlist(strsplit(as.character(x), split=' '))))}





##### Read in data #####
data <- DATA <- read.csv("~/Desktop/beeLandings.csv",header=T,sep=',')
### Adding results from the simultions
load('~/Desktop/simulation1')
DF_sim_1 <- DF
load('~/Desktop/simulation2')
DF_sim_2 <- DF
### Transform empirical data set
landings <- lapply(data$landed, function(x) as.numeric(unlist(strsplit(as.character(x), split=' '))))
cue <- lapply(data$cue, function(x) as.numeric(unlist(strsplit(as.character(x), split=' '))))
### Check length of input
if(length(landings)!=length(cue)) stop('objects have different length')
### Create landings list
llandings <- llandingsF(l=landings, c=cue)
### Add info to data.frame
data <- data_add_infoF(dat=data, ll=llandings)
### Setting the right world (diff=high-varaince / easy=no-varaince)
data$world <- factor(ifelse(test = data$group=='test_easy', yes = 'easy', no = 'diff'))
### Extracting worlds for last cue(=dummy) training
dummy_worlds <- do.call(rbind,lapply(unique(data$id), function(id){
	data.frame(id=id, world=data[data$id==id & data$group!="dummy_last", "world"])
}))
### Writing the worlds to the according id's
sort_id <- match(x=dummy_worlds$id, table=data$id[data$group=="dummy_last"])
data$world[data$group=="dummy_last"] <- dummy_worlds$world[sort_id]
### Order data and fill with -1 to bring all vectors the same length
longest <- max(unlist(lapply(llandings, length)))
r <- data.frame(do.call(rbind, lapply(llandings, function(x)c(x,rep(-1,longest-length(x)) ))))



##### Extract flower choices #####
### Bee dummy
x <- data$cueType=="clay" & data$group%in%c("test_diff","test_easy")
### First Flower Choice
diffFirstChoice_clay <- data$firstChoice[data$world=='diff' & x ]
easyFirstChoice_clay <- data$firstChoice[data$world=='easy' & x]
### First 4 choices
diffFirst4Avg_clay <- data$first4CueProp[data$world=='diff' & x]
easyFirst4Avg_clay <- data$first4CueProp[data$world=='easy' & x]
### First 10 choices
diffFirst10Avg_clay <- data$first10CueProp[data$world=='diff' & x]
easyFirst10Avg_clay <- data$first10CueProp[data$world=='easy' & x]

### Craft Foam dummy
y <- data$cueType=="green" & data$group%in%c("test_diff","test_easy")
### First Flower Choice
diffFirstChoice_green <- data$firstChoice[data$world=='diff' & y]
easyFirstChoice_green <- data$firstChoice[data$world=='easy' & y]
### First 4 choices
diffFirst4Avg_green <- data$first4CueProp[data$world=='diff' & y]
easyFirst4Avg_green <- data$first4CueProp[data$world=='easy' & y]
### First 10 choices
diffFirst10Avg_green <- data$first10CueProp[data$world=='diff' & y]
easyFirst10Avg_green <- data$first10CueProp[data$world=='easy' & y]





##### Main Text Stistical tests #####
#### Results from Test ####
### First flower choice different from chance (AS REPORTED IN MAIN TEXT)?
### Bee dummy - high-variance
binom.test(x=c(sum(diffFirstChoice_clay==1),sum(diffFirstChoice_clay==0)), p=1/3)#, alternative="greater")
### Bee dummy - no-variance
binom.test(x=c(sum(easyFirstChoice_clay==1),sum(easyFirstChoice_clay==0)), p=1/3)#, alternative="less")
### Craft foam - high-variance
binom.test(x=c(sum(diffFirstChoice_green==1),sum(diffFirstChoice_green==0)), p=1/3)#, alternative="greater")
### Craft foam - no-variance
binom.test(x=c(sum(easyFirstChoice_green==1),sum(easyFirstChoice_green==0)), p=1/3)#, alternative="less")



#### Logistic Regression ####
### Selecting test data
df <- data[data$group!="dummy_last",]
df$cueType <- factor(df$cueType, levels=c("clay","green"))
df$world <- factor(df$world, levels=c("diff","easy"))

#### Setting 1's and 0's for cue type and reward distribution
# 1== bee dummy, 0== craft foam
# 1== high-variance, 0== no-variance
df$cueType <- factor(ifelse(df$cueType=='clay',yes=1,no=0))
df$world <- factor(ifelse(df$world=='diff',yes=1,no=0))

# All variables with colony as random effect (AS REPORTED IN MAIN TEXT)
mod_col <- glmer(data = df, family = "binomial", formula = firstChoice ~ cueType + world + cueType*world + (1|colony))
summary(mod_col)


### Splitting the data for the two different cues (AS REPORTED IN MAIN TEXT)
mod_clay <- glmer(data = df[df$cueType=="1",], family = "binomial", formula = firstChoice ~ world + (1|colony))
summary(mod_clay)
mod_green <- glmer(data = df[df$cueType=="0",], family = "binomial", formula = firstChoice ~ world + (1|colony))
summary(mod_green)

### Splitting the data for the two different reward distributions (not reported)
mod_uneven <- glmer(data = df[df$world==0,], family = "binomial", formula = firstChoice ~ cueType + (1|colony))
summary(mod_uneven)
mod_even <- glmer(data = df[df$world==1,], family = "binomial", formula = firstChoice ~ cueType + (1|colony))
summary(mod_even)












#### Comparisson of empirical data with results from simulations ####
### Data flower choice step by step ###
datas <- do.call(rbind, lapply(c("clay","green"), function(cueT){
	do.call(rbind, lapply(c('diff','easy'), function(w){
		do.call(rbind, lapply(c('firstChoice','first4CueProp','first10CueProp'), function(ch){
			d <- data[data$world==w & data$cueType==cueT & data$group%in%c("test_diff","test_easy"), colnames(data)==ch]
			data.frame(world=w, cueType=cueT,choice=ch,avg=mean(d),sd=sd(d),se=se(d),n=length(d))
			}))
		}))
	}))
### Transform data from simulations
sim_diff <- unlist(DF_sim_1$meanILprop)
mean_diff <- 1-mean(sim_diff)
sim_easy <- unlist(DF_sim_2$meanILprop)
mean_easy <- 1-mean(sim_easy)
### Adding a column to combine the test conditions (cueType and world trained)
datas_first <- datas[datas$choice=="firstChoice",]
datas_first$test <- paste(datas_first$world, as.character(datas_first$cueType), sep='_')
datas_first$test <- factor(datas_first$test, levels=c("diff_clay","easy_clay","diff_green","easy_green"))
### Adding simulation data
datas_first <- rbind(data.frame(world=c('diff','easy'), cueType='simulation', choice='first', avg=c(mean_diff,mean_easy), se=c(se(sim_diff),se(sim_easy)),sd=c(sd(sim_diff),sd(sim_easy)),n=25, test=c('diff_sim','easy_sim')), datas_first)

### Comparison to Bee Model (AS REPORTED IN MAIN TEXT)
binom.test(x=c(sum(diffFirstChoice_clay==1),sum(diffFirstChoice_clay==0)), p=mean_diff)
binom.test(x=c(sum(easyFirstChoice_clay==1),sum(easyFirstChoice_clay==0)), p=mean_easy)
### Comparison to Craft Foam (AS REPORTED IN MAIN TEXT)
binom.test(x=c(sum(diffFirstChoice_green==1),sum(diffFirstChoice_green==0)), p=mean_diff)
binom.test(x=c(sum(easyFirstChoice_green==1),sum(easyFirstChoice_green==0)), p=mean_easy)









#### Results from last training round ####
### Did bees differ in training success? — Last Training Data ###
### Subsetting raw data for training with bee models and craft foam
data_lastTraining <- data[data$group=='dummy_last', ]
### Transforming strings of landings and places of cues to numeric vectors for cue
landings_lastTraining <- landingsF(data_lastTraining)
cue_lastTraining <- cueF(data_lastTraining)
reward_lastTraining <- rewardF(data_lastTraining)
### Did the chips landed on were with cue? Returns vectors with T and F
llandings_lastTraining <- llandingsF(l=landings_lastTraining, c=cue_lastTraining)
### Extracting additional information from the landings and adding to the data.frames
data_lastTraining <- data_add_infoF(dat = data_lastTraining, ll=llandings_lastTraining)
### ... and more data that could be important: unique flowers visited, unique cues landed on, collected rewards, visits befpre first and second reward
data_lastTraining <- data_add_infoF2(dat=data_lastTraining, l=landings_lastTraining, c=cue_lastTraining, r=reward_lastTraining)
### Adding informaiton whether bees were tested later on in easy or diff
data_lastTraining$world <- data$world[ match(c(as.character(data_lastTraining$id)), as.character(data$id)) ]

### First Flower Choice last training ###
### Bee Model
x <- data_lastTraining$cueType=="clay"# & data_lastTraining$group=="dummy_last"
diffFirstChoice_clay_last <- data_lastTraining$firstChoice[data_lastTraining$world=='diff' & x ]
easyFirstChoice_clay_last <- data_lastTraining$firstChoice[data_lastTraining$world=='easy' & x ]
### Craft Foam Dummy
y <- data_lastTraining$cueType=="green"# & data_lastTraining$group=="dummy_last"
diffFirstChoice_green_last <- data_lastTraining$firstChoice[data_lastTraining$world=='diff' & y ]
easyFirstChoice_green_last <- data_lastTraining$firstChoice[data_lastTraining$world=='easy' & y ]


#### For First Choice (test of equal proportions) in the last cue learning round, do groups differ from each other? (NOT REPORTED)
### Bee Model
firstChoice_last <- matrix(c(sum(diffFirstChoice_clay_last),sum(!diffFirstChoice_clay_last), sum(easyFirstChoice_clay_last),sum(!easyFirstChoice_clay_last)), nrow=2)
colnames(firstChoice_last) <- c('diff','easy')
rownames(firstChoice_last) <- c('cue','nocue')
firstChoice_last
fisher.test(firstChoice_last)
### Craft Foam Dummy
firstChoice_last <- matrix(c(sum(diffFirstChoice_green_last),sum(!diffFirstChoice_green_last), sum(easyFirstChoice_green_last),sum(!easyFirstChoice_green_last)), nrow=2)
colnames(firstChoice_last) <- c('diff','easy')
rownames(firstChoice_last) <- c('cue','nocue')
firstChoice_last
fisher.test(firstChoice_last)


### Is the cue use from bees trained in high/no variance significantly different from chance (1/3) in last cue training round? (AS REPORTED IN MAIN TEXT)
# Bee Model - High-variance
binom.test(x=c(sum(diffFirstChoice_clay_last==1),sum(diffFirstChoice_clay_last==0)), p=1/3)
# Bee Model - No-variance
binom.test(x=c(sum(easyFirstChoice_clay_last==1),sum(easyFirstChoice_clay_last==0)), p=1/3)
# Craft Foam - high-variance
binom.test(x=c(sum(diffFirstChoice_green_last==1),sum(diffFirstChoice_green_last==0)), p=1/3)
# Craft Foam - no-variance
binom.test(x=c(sum(easyFirstChoice_green_last==1),sum(easyFirstChoice_green_last==0)), p=1/3)


















##### Main Text Figures #####
### Figure 1a
# Error bars
limits <- aes(x=world, y=avg, ymax = avg + se, ymin=avg - se)
q+geom_bar(data=datas_first[datas_first$cueType=='simulation',], aes(x=world, y=avg, fill=world, group=world), colour=c(NA,"black"), width=1, stat='identity') + geom_errorbar(data=datas_first[datas_first$cueType=='simulation',], limits, width=.1, size=.25, colour=c('grey25','black'))+ylab('')+guides(fill=FALSE) + theme(axis.ticks.x=element_blank())

### Figure 1b
# Confidence intervalls, using http://www.graphpad.com/quickcalcs/ConfInterval1.cfm
datas_first$CI_lower <- c(NA, NA, 0.6562, 0.1068, 0.3555, 0.2195)
datas_first$CI_upper <- c(NA, NA, 0.9773, 0.4365, 0.7805, 0.6445)
q + geom_bar(data=datas_first[datas_first$cueType!='simulation',], aes(x=cueType, y=avg, fill=world), colour="black", position='dodge', width=.8, stat='identity') + ylab('')+guides(fill=FALSE)+ geom_hline(yintercept=1/3, linetype='dashed', colour='black', size=.5)+ geom_errorbar(data=datas_first[datas_first$cueType!='simulation',], aes(x=cueType, y=avg, ymin=CI_lower, ymax=CI_upper, group=world), width=.1, position = position_dodge(0.8), colour='grey40', size=.25)











##### SI #####
#### Gini index for simulation ####
gini <- function(x, na.rm = FALSE){
	if (!is.numeric(x)){
		warning("'x' is not numeric; returning NA")
		return(NA)
	}
	if (!na.rm && any(na.ind <- is.na(x)))
		stop("'x' contain NAs")
	if (na.rm)
		x <- x[!na.ind]
	n <- length(x)
	mu <- mean(x)
	N <- n * (n + 1)
	ox <- x[order(x)]
	dsum <- drop(2*crossprod(1:n,  ox))
	(dsum / (mu * N)) - 1
}

MEAN <- 100/12 #(100µl in 12 flowers)
K <- (MEAN^2)/var(c(rep(0,10),rep(50,2)))
THETA <- MEAN/K
rand <- rgamma(n = 100000, shape = K, scale=THETA)
### Calculating Gini index of resource distribution used in our experiment to compare with Smolla, Gilman, Galla, and Shultz (2015)
gini(rand)



#### SI figures ####
### Print version for Figure S1a and S2a:

### Adding factor levels
datas$choice <- factor(c('first','first 4',' first 10'), levels=c(c('first','first 4',' first 10')))
### Setting standard error bars
limits <- aes(x=world, y=avg, ymax = avg + sd, ymin=avg - sd)
### Poprotion of landings in first, first four, and first ten landings on flowers with a cue
x <- 'clay' # Bee dummy
# x <- 'green' # Craft foam dummy
dat <- datas[datas$cueType==x,]
qchoices <- q + geom_bar(data=dat, aes(x=world, y=avg, fill=world), stat='identity', colour=c('black'), width=1)  + geom_errorbar(data=dat[dat$choice!='first',], limits, width=.1, size=.25, colour=c('grey25','black','grey25','black')) + xlab(label='') + geom_hline(yintercept=1/3, linetype='dashed', colour='grey25', size=.5) + facet_grid(.~choice) + ylab('') + guides(fill=FALSE) + theme(axis.ticks.x=element_blank(), panel.margin = unit(.1, "lines"), strip.text=element_blank()) + theme(axis.title.x = element_text(size=15,family='Helvetica'),
				axis.title.y = element_text(size=15, angle=90,family='Helvetica'),
				axis.text.x = element_text(size=13,family='Helvetica'),
				axis.text.y = element_text(size=13,family='Helvetica'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				legend.title = element_text(size=13, face='bold'),
				legend.text = element_text(size=13, family='Helvetica'))
### Landings per reward distribution, step-by-step
countData <-
	do.call(rbind, lapply(c("clay","green"), function(cueT){
		do.call(rbind, lapply(c(1:10), function(i){
		x <-  data$cueType==cueT & data$group%in%c("test_diff","test_easy")
		data.frame(cueType=cueT, step=i,
							 diff=sum(r[data$world=='diff' & x,i]==1),#/length(r[data$world=='diff' & x,i]),
							 easy=sum(r[data$world=='easy' & x,i]==1))#/length(r[data$world=='easy' & x,i]))
	}))
}))
countData$id <- 1:nrow(countData)
step_cue_landing_melted <- melt(countData,value.name=c('landingOnCue'), varnames='landing', id=c('id','step','cueType'))

### Print Version for Figure S1b and S2b:#SP15HJ
x <- 'clay' # Bee dummy
# x <- 'green' # Craft foam dummy
qsteps <- q+geom_bar(data=step_cue_landing_melted[step_cue_landing_melted$cueType==x,], aes(x=factor(step), y=landingOnCue, fill=variable), colour='black', size=.2, stat='identity', position='dodge') + xlab('Landing number') + ylab('') + guides(fill=FALSE) + scale_fill_manual(values=c("black", "white")) + scale_y_continuous(lim=c(0,15.2), expand=c(0,0.00)) + theme(axis.title.x = element_text(size=15,family='Helvetica'),
				axis.title.y = element_text(size=15, angle=90,family='Helvetica'),
				axis.text.x = element_text(size=13,family='Helvetica'),
				axis.text.y = element_text(size=13,family='Helvetica'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				legend.title = element_text(size=13, face='bold'),
				legend.text = element_text(size=13, family='Helvetica'))
### Plot figures
ggdraw() +
  draw_plot(qchoices, 0, 0, .33, 1) +
  draw_plot(qsteps, .33, 0, .66, 1) +
  draw_plot_label(c("a", "b"), fontface='italic', c(0, .33), c(1, 1), size = 15)



#### SI Statistics ####
### Extracting landings from the last first 4 and first 10 landings
first4 <- unlist(lapply(llandings, function(x) {ifelse(length(x)<4 , yes = sum(x), no = sum(x[1:4]) )}))
first10 <- unlist(lapply(llandings, function(x) {ifelse(length(x)<10 , yes = sum(x), no = sum(x[1:10]) )}))

x <- data$cueType=="clay" & data$group%in%c("test_diff","test_easy")
### First 4 choices
diffFirst4_clay <- first4[data$world=='diff' & x]
easyFirst4_clay <- first4[data$world=='easy' & x]
### First 10 choices
diffFirst10_clay <- first10[data$world=='diff' & x]
easyFirst10_clay <- first10[data$world=='easy' & x]

y <- data$cueType=="green" & data$group%in%c("test_diff","test_easy")
### First 4 choices
diffFirst4_green <- first4[data$world=='diff' & y]
easyFirst4_green <- first4[data$world=='easy' & y]
### First 10 choices
diffFirst10_green <- first10[data$world=='diff' & y]
easyFirst10_green <- first10[data$world=='easy' & y]


#### For first 4 and first 10 choices testing whether landings differ from random choice
### First four for bee dummy
binom.test(x=c(sum(diffFirst4_clay),(length(diffFirst4_clay)*4)-sum(diffFirst4_clay)), p=1/3)
binom.test(x=c(sum(easyFirst4_clay),(length(easyFirst4_clay)*4)-sum(easyFirst4_clay)), p=1/3)
### First ten for bee foam
binom.test(x=c(sum(diffFirst10_clay),(length(diffFirst10_clay)*10)-sum(diffFirst10_clay)), p=1/3)
binom.test(x=c(sum(easyFirst10_clay),(length(easyFirst10_clay)*10)-sum(easyFirst10_clay)), p=1/3)

### First four for craft foam
binom.test(x=c(sum(diffFirst4_green),(length(diffFirst4_green)*4)-sum(diffFirst4_green)), p=1/3)
binom.test(x=c(sum(easyFirst4_green),(length(easyFirst4_green)*4)-sum(easyFirst4_green)), p=1/3)
### First ten for craft foam
binom.test(x=c(sum(diffFirst10_green),(length(diffFirst10_green)*10)-sum(diffFirst10_green)), p=1/3)
binom.test(x=c(sum(easyFirst10_green),(length(easyFirst10_green)*10)-sum(easyFirst10_green)), p=1/3)


### Do groups trained in different reward distributions differ between each other for the first 4 and first 10 landings?
### Bee dummy first four landings
first4ChoiceM <- matrix(c(sum(diffFirst4_clay),(length(diffFirst4_clay)*4)-sum(diffFirst4_clay), sum(easyFirst4_clay),(length(easyFirst4_clay)*4)-sum(easyFirst4_clay)), nrow=2)
colnames(first4ChoiceM) <- c('diff','easy')
rownames(first4ChoiceM) <- c('cue','nocue')
first4ChoiceM
chisq.test(first4ChoiceM)

### Bee dummy first ten landings
first10ChoiceM <- matrix(c(sum(diffFirst10_clay),(length(diffFirst10_clay)*10)-sum(diffFirst10_clay), sum(easyFirst10_clay),(length(easyFirst10_clay)*10)-sum(easyFirst10_clay)), nrow=2)
colnames(first10ChoiceM) <- c('diff','easy')
rownames(first10ChoiceM) <- c('cue','nocue')
chisq.test(first10ChoiceM)

### Craft foam first four landings
first4ChoiceM <- matrix(c(sum(diffFirst4_green),(length(diffFirst4_green)*4)-sum(diffFirst4_green), sum(easyFirst4_green),(length(easyFirst4_green)*4)-sum(easyFirst4_green)), nrow=2)
colnames(first4ChoiceM) <- c('diff','easy')
rownames(first4ChoiceM) <- c('cue','nocue')
chisq.test(first4ChoiceM)

### Craft foam first ten landings
first10ChoiceM <- matrix(c(sum(diffFirst10_green),(length(diffFirst10_green)*10)-sum(diffFirst10_green), sum(easyFirst10_green),(length(easyFirst10_green)*10)-sum(easyFirst10_green)), nrow=2)
colnames(first10ChoiceM) <- c('diff','easy')
rownames(first10ChoiceM) <- c('cue','nocue')
chisq.test(first10ChoiceM)

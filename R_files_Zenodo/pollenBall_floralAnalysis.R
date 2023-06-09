# James, Geber, and Toews 
## Code for Molecular Assays of Pollen Use Consistently Reflect Pollinator Visitation Patterns in a System of Flowering Plants

# clear 
#rm(list=ls())

# packages 
require(reshape2)
require(ggplot2)
require(tidyr)
require(rstatix)
require(dplyr)
require(bipartite)
require(data.table)
require(vegan)

# setwd as desired
setwd('~/Desktop/revisions MER 2021/')

# data 
d<-read.csv("qAMPdat.csv", header=T)

# format date
d$date<-as.Date(a$date, format=c("%m/%d/%y"), tryFormats=c("%m/%d/%Y"))

names(d)<-c("date", "nc", "ns", "nu", "nx", "bgenus",  'bsubgenus',"bspecies","id", "flid", "am", "site", "type", "cct", "sct", "uct", "xct", "crel","srel", "urel", "xrel") #rab stands for relative abundance

# get relative amount of flowers in community 
d$flct<-rowSums(d[,c(2:5)])
d$pc<-d$nc/d$flct
d$ps<-d$ns/d$flct
d$pu<-d$nu/d$flct
d$px<-d$nx/d$flct

head(d)

# add differences between community and pollen balls  
fl<-d[,(23:26)] # floral abundances 
rra<-d[,(18:21)] # relative read abundance
rra[rra < 0.05]<-0 # this is optional - adjust rra by giving it a cutoff
qamp<-d[,(14:17)] # qamp estimate 

fl.qamp<-qamp-fl # the difference in proportion between flowers and qampseq
fl.rra<-rra-fl # the difference in proportion between flowers and rra

# qamp diffs 
d$cdiff<-fl.qamp[,1]
d$sdiff<-fl.qamp[,2]
d$udiff<-fl.qamp[,3]
d$xdiff<-fl.qamp[,4]

# rra diffs
d$cdiff2<-fl.rra[,1]
d$sdiff2<-fl.rra[,2]
d$udiff2<-fl.rra[,3]
d$xdiff2<-fl.rra[,4]

# melt 
e<-d[,c(1,6:12,27:34)]

# remove single species sites 
e<-subset(e, site=c("BG"|"KS"|"LCG"|"LT"|"MC"|"S8"|"SC"|"UCG"))

# melt 
dat<-melt(e, id.vars=c("date","bgenus","bspecies","id","flid","am","site"), measure.vars=c("cdiff", "sdiff", "udiff", "xdiff", "cdiff2", "sdiff2", "udiff2", "xdiff2"))
diffs<-dat$variable[c(1:608)]
dat$diffs<-rep(diffs,2)
dat$Method<-c(rep("qAMPseq", 608), rep("RRA", 608))

## analysis for preference
rracut<-dat$value[which(dat$Method=='RRA')]

## t test between datasets
nudat<-cbind(dat[which(dat$Method=='qAMPseq'),'value'], rracut)
t.test(nudat[,1], nudat[,2], paired=T)

# use qampseq only, because not different
half<-dat[which(dat$Method=='qAMPseq'),]
prefnov<-lm(value~diffs-1, data=half)

# complete cases data 
half<-half[which(complete.cases(half==T)),]

##
cols=c("#6F6CB2","#F9C2B1")

ggplot(half, aes(diffs, value))+geom_hline(yintercept=0, col='gray88', lty='dashed', lwd=1)+stat_summary(position=position_dodge(width=0.5), size=1.2, fun.data='mean_cl_boot')+theme_classic()+scale_color_manual(values=cols)+labs(x="Species", y="Preference")+theme(axis.text.x=element_text(size=13, , face="italic"), axis.text.y=element_text(size=13), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))+scale_x_discrete(labels=c("C. cylindrica", "C. speciosa", "C. unguiculata", "C. xantiana"))

### 
bp<-d[,c(2:5,18:21, 14:17)]
bp<-bp[which(complete.cases(bp)==T),]

## shaping the data for plotting 
## without threshholds for RRA bp
# add number of species in sample types and communities. 

## with 0% threshhold added to RRA 
ap<-bp
ap[ap != 0]<- 1
ap$sample<-seq(1:nrow(ap))
ap$community<-rowSums(ap[1:4])
ap$rra<-rowSums(ap[5:8])
ap$qamp<-rowSums(ap[9:12])

## with 5 and 10 % threshholds added to the RRA 
fivep<-bp[,c(5:8)]
fivep$crel[which(fivep$crel<0.05)]<-0
fivep$srel[which(fivep$srel<0.05)]<-0
fivep$urel[which(fivep$urel<0.05)]<-0
fivep$xrel[which(fivep$xrel<0.05)]<-0
fivep[fivep != 0]<- 1
ap$rra5<-rowSums(fivep)

tenp<-bp[,c(5:8)]
tenp$crel[which(tenp$crel<0.10)]<-0
tenp$srel[which(tenp$srel<0.10)]<-0
tenp$urel[which(tenp$urel<0.10)]<-0
tenp$xrel[which(tenp$xrel<0.10)]<-0
tenp[tenp != 0]<- 1
ap$rra10<-rowSums(tenp)

# remove samples where estimated qampseq = 0 (not enough to estimate proportions from)
ap<-ap[-c(138,144),]

# melt 
bq<-melt(ap, measure.vars=c("community", 'qamp', 'rra', 'rra5', 'rra10'))
bq$value<-as.factor(bq$value)
bq<-plyr::rename(bq, replace=c("value"="N.Species"))

# plot
bq$variable<-factor(bq$variable, levels=c("community", "qamp",  "rra10", "rra5", 'rra'))
cols1=c("white", "gray90", "gray50", "black")

ggplot(bq, aes(x=variable, group=N.Species, fill=as.factor(N.Species)))+geom_bar(stat='count', position='stack', color='black')+scale_fill_manual(values=cols1)+theme_bw()+scale_x_discrete(labels=c("Community", "qAMPseq", "RRA, 10%", "RRA, 5%", "RRA"))+labs(x="Method", y="Sample Count")+theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))

countcompare<-as.data.frame(bq %>% group_by(N.Species, variable) %>% tally())
#write.csv(countcompare, 'fig4_table.csv')

## probably want to do a chisq test 
countest<-countcompare[which(countcompare$variable=='community'|countcompare$variable=='qamp'),]
setorder(countest, cols='variable')

mat<-matrix(c(44,63,9,34,114,31,5,0), ncol=4, byrow=T)
rownames(mat)=c("comm","qamp")
colnames(mat)=c("1","2","3","4")

chisq.test(mat)

############ look at bees with diff flowers 
## optional figure  
bd<-d[,c(6:8,10,14:17)]
bd<-bd[which(complete.cases(bd)==T),]

cols2<-c("#6498CE","#FFF25F","#FF8C36","#FF392F")
bd2<-melt(bd, id.vars=c("bgenus","bsubgenus","bspecies","flid"))
pdat<-as.data.frame(bd2 %>% group_by(bgenus, bsubgenus, bspecies, variable) %>% summarize(sumb=sum(value, na.rm=T)))
pdat2<-pdat%>%unite("b", 1:3, sep=" ")

ggplot(pdat2, aes(x=b, y=sumb, fill=variable))+geom_col(color="#223447")+scale_fill_manual(values=cols2)+theme_classic()+theme(axis.text.x=element_text(size=15, face="italic", angle=30, hjust=1), axis.text.y=element_text(size=8), axis.title.x=element_text(size=25), axis.title.y=element_text(size=25))+labs(x="Bee species", y="Relative Pollen Contents")+labs(fill="Species")

########## look at diff from zero
anovadat<-dat[which(dat$Method=="qAMPseq"),]
summary.lm(aov(formula = value ~ diffs-1, data = anovadat))
TukeyHSD(aov(formula = value ~ diffs-1, data = anovadat))


#Make networks 
## from where bees were caught 
comm_net<-d[c(2:8,10)]
comm_net$bid<-paste(comm_net$bgenus, comm_net$bspecies)
comm_web<-table(comm_net$bid, comm_net$flid)
comm_web<-as.matrix(comm_web[c(1:7),])
H2fun(comm_web)

## from qamp
qamp_net<-d[c(6:8,14:17)]
c<-qamp_net[,c(4:7)]*100
c<-round(c)
qamp_net$bid<-paste(qamp_net$bgenus, qamp_net$bspecies)
qamp_net<-cbind(qamp_net, c)
qamp_net<-qamp_net[which(complete.cases(qamp_net)==T),]

qamp_net<-as.data.frame(qamp_net[,c(8:12)])

qnet<-qamp_net%>%group_by(bid)%>%summarise_each(sum)
qnet2<-as.data.frame(qnet[,c(2:5)])
rownames(qnet2)<-qnet$bid

## from RRA
rra_net<-d[c(6:8,18:21)]

# cutoff for rra 5%
rra_net$crel[which(rra_net $crel<0.05)]<-0
rra_net$srel[which(rra_net $srel<0.05)]<-0
rra_net$urel[which(rra_net $urel<0.05)]<-0
rra_net$xrel[which(rra_net $xrel<0.05)]<-0

# carry on as before!
c2<-rra_net[,c(4:7)]*100
c2<-round(c2)
rra_net$bid<-paste(rra_net$bgenus, rra_net$bspecies)
rra_net<-cbind(rra_net, c2)
rra_net<-rra_net[which(complete.cases(rra_net)==T),]
rra_net<-as.data.frame(rra_net[,c(8:12)])

rnet<-rra_net%>%group_by(bid)%>%summarise_each(sum)
rnet2<-as.data.frame(rnet[,c(2:5)])
rownames(rnet2)<-rnet$bid

## plot all networks
par(mfrow=c(1,3))
plotweb(comm_web, method='normal', col.interaction='#6498B0', col.high='#A498BC', col.low='#A498BC', bor.col.interaction='#000000')
plotweb(qnet2, method='normal', col.interaction='#6498B0', col.high='#A498BC', col.low='#A498BC', bor.col.interaction='#000000')
plotweb(rnet2, method='normal', col.interaction='#6498B0', col.high='#A498BC', col.low='#A498BC', bor.col.interaction='#000000')

## comparing abundances 
abuncompare<-melt(abun, measure.vars=c('nc','ns','nu', 'nx'))
acomp<-abuncompare[which(complete.cases(abuncompare==T)),]
amod<-aov(value~variable, data=acomp)
summary(amod)
TukeyHSD(amod)

# change entries? 
acomp$variable<-gsub("nc", "C. cylindrica", acomp$variable)
acomp$variable<-gsub("ns", "C. speciosa", acomp$variable)
acomp$variable<-gsub("nu", "C. unguiculata", acomp$variable)
acomp$variable<-gsub("nx", "C. xantiana", acomp$variable)

# plot
ggplot(acomp, aes(x=variable, y=value))+geom_jitter(alpha=0.05, size=2, width=0.01)+stat_summary(col='red', size=1)+theme_bw()+theme(axis.text.x=element_text(size=15, face="italic"), axis.text.y=element_text(size=18), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+labs(x="Species", y="log(Abundance)")+labs(fill="Species")+annotate(geom='text', x=1, y=7, label='ab', size=6)+annotate(geom='text', x=2, y=4.5, label='c', size=6)+annotate(geom='text', x=3, y=6, label='bc', size=6)+annotate(geom='text', x=4, y=5.5, label='c', size=6)

#### procrustes analysis to compare network topology (vegan package)
?procrustes #
?protest #tests the significance of the nonrandomness of two configurations (matrices). if significant, they are nonrandom and statistically indistinguishable :) uses permutations 
# matrices to use 
qnet2 # qampseq matrix
comm_web # visitation matrix
protest(qnet2, comm_web)
summary(protest(comm_web, qnet2)) # the matrix makes it pretty obvious that the scaling is almost exactly the identity matrix lol so it's not different. still can look at the output 

# matrices to use 
# rnet2 # qampseq matrix
# comm_web # visitation matrix
protest(rnet2, comm_web)
summary(protest(comm_web, rnet2)) 

# matrices to use 
rent2 # qampseq matrix
qnet2 # visitation matrix

protest(rnet2, qnet2)
summary(protest(qnet2, rnet2)) 



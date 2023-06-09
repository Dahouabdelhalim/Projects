##### trying to label the mixed egg masses
library(ggplot2)

#rawdata
map_pop<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/IndivPlast_MappingPops/Date_forR.csv", header=TRUE, sep=",")
SvN<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/SoCALvsNorCAL/SoCalvNorCal_Data.csv", header=TRUE, sep=",")

#combine into HvL

#get date into the right format
map_pop$date<-as.Date(map_pop$date, "%m/%d/%Y")
SvN$date<-as.Date(SvN$date, "%m/%d/%Y")

sel_all<-rbind(SvN[,c("Capsule", "original", "date", "Family", "ID", "eggmass")],map_pop[, c("Capsule", "original", "date", "Family", "ID", "eggmass")])

#if we aggregate by date, eggmass and ID, we should get a id of each eggmass 
#should be able to use aggregate. Want ID, by date, and eggmass

sel_all_mlp<-aggregate(sel_all$Capsule, list(sel_all$date, sel_all$eggmass, sel_all$ID, sel_all$original), function(x) ifelse(x > 0.15, "l", "p"))

#this returns a vector (x) that contains for each eggmass a id of each egg capsule that was measured. 

#rename columns

colnames(sel_all_mlp)<-c("date", "eggmass", "ID", "original", "x")

#we now want to ID the groups that have l and p, these are our mixed eggmasses
#can do this in two parts, first seeing which are l
          
sel_all_mlp$is_l<-grepl("l", sel_all_mlp$x)

#now which are p

sel_all_mlp$is_p<-grepl("p", sel_all_mlp$x)


#and now which are true for both l and p

sel_all_mlp$is_m<- sel_all_mlp$is_l == T & sel_all_mlp$is_p == T

table(sel_all_mlp$is_m)


#now make a vector that actually identifies each eggmass as l, p, or m


sel_all_mlp$eggmass_type<-ifelse(sel_all_mlp$is_m == T, "m", ifelse(sel_all_mlp$is_l == T, "l", "p"))

table(sel_all_mlp$eggmass_type)

#success

#look to see how many egg masses were laid in high versus low salinity
table(sel_all_mlp$eggmass_type,sel_all_mlp$original)



#lets just look at map pop to see if we have more mixed eggmasses after selection

map_pop_mlp<-aggregate(map_pop$Capsule, list(map_pop$date, map_pop$eggmass, map_pop$ID), function(x) ifelse(x > 0.15, "l", "p"))

colnames(map_pop_mlp)<-c("date", "eggmass", "ID", "x")
map_pop_mlp$is_l<-grepl("l", map_pop_mlp$x)
map_pop_mlp$is_p<-grepl("p", map_pop_mlp$x)
map_pop_mlp$is_m<- map_pop_mlp$is_l == T & map_pop_mlp$is_p == T
map_pop_mlp$eggmass_type<-ifelse(map_pop_mlp$is_m == T, "m", ifelse(map_pop_mlp$is_l == T, "l", "p"))
table(map_pop_mlp$eggmass_type)

#then merge with map_pop data to get the generation back in

map_pop_mlp_tog<-merge(map_pop_mlp, map_pop[,c("Family", "original","date", "eggmass", "ID", "gen")], by=c("date", "eggmass", "ID"))

table(map_pop_mlp_tog$eggmass_type, map_pop_mlp_tog$gen)


# p   F1   F2   F3
# l  299  366  156    0
# m  167   88   48    0
# p 2690  873   77    8


table(map_pop_mlp_tog$eggmass_type, map_pop_mlp_tog$gen:map_pop_mlp_tog$original)

# p:20C16ppt p:20C32ppt F1:20C16ppt F1:20C32ppt F2:20C16ppt F2:20C32ppt F3:20C16ppt F3:20C32ppt
# l         55        244          97         269          33         123           0           0
# m         62        105          37          51           4          44           0           0
# p       1171       1519         226         647           2          75           4           4

map_pop_mlp_tog$region<-ifelse(grepl("SC*", map_pop_mlp_tog$Family), "Long Beach", ifelse(grepl("W128|13[0-4]", map_pop_mlp_tog$Family), "Tomales", "Mill Valley"))


#proportion of m increases as does l by the F2 generation
# > sum(299, 167, 2690)
# [1] 3156
# > sum(366, 88, 873)
# [1] 1327
# > sum(156, 48, 77)
# [1] 281

# l
# > 299/3156
# [1] 0.09474018
# > 366/1327
# [1] 0.2758101
# > 156/281
# [1] 0.5551601

# m
# > 167/3156
# [1] 0.05291508
# > 88/1327
# [1] 0.066315
# > 48/281
# [1] 0.1708185



#Do all this for sel all data frame which contains all the data for all the things

sel_sal<-read.table("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/SoCALvsNorCAL/F1s/Sel_sal1.csv", header=TRUE, sep=",")

map_pop<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/IndivPlast_MappingPops/Date_forR.csv", header=TRUE, sep=",")

#merge them together

head(sel_sal)
head(map_pop)


colnames(F1_sal)<-c("date", "order", "Capsule", "mom_id", "mom_env", "ID", "eggmass", "original", "Family")


map_pop$date<-as.Date(map_pop$date, "%m/%d/%Y")
sel_sal$date<-as.Date(sel_sal$date, "%m/%d/%Y")

cols <- intersect(colnames(sel_sal[,-2]), colnames(map_pop[,-2]))
sel_all<- rbind(sel_sal[,cols], map_pop[,cols])

sel_all_mlp<-aggregate(sel_all$Capsule, list(sel_all$date, sel_all$eggmass, sel_all$ID, sel_all$original), function(x) ifelse(x > 0.15, "l", "p"))

#this returns a vector (x) that contains for each eggmass a id of each egg capsule that was measured. 
head(sel_all_mlp)
#rename columns

colnames(sel_all_mlp)<-c("date", "eggmass", "ID", "original", "x")

#we now want to ID the groups that have l and p, these are our mixed eggmasses
#can do this in two parts, first seeing which are l

sel_all_mlp$is_l<-grepl("l", sel_all_mlp$x)

#now which are p

sel_all_mlp$is_p<-grepl("p", sel_all_mlp$x)


#and now which are true for both l and p

sel_all_mlp$is_m<- sel_all_mlp$is_l == T & sel_all_mlp$is_p == T

table(sel_all_mlp$is_m)


#now make a vector that actually identifies each eggmass as l, p, or m


sel_all_mlp$eggmass_type<-ifelse(sel_all_mlp$is_m == T, "m", ifelse(sel_all_mlp$is_l == T, "l", "p"))

table(sel_all_mlp$eggmass_type)

#success

#look to see how many egg masses were laid in high versus low salinity

table(sel_all_mlp$eggmass_type,sel_all_mlp$original)
# 20C16ppt 20C32ppt
# l      183      465
# m       70      114
# p      784     1230

#to get high and low salinity by generation we will ahve to merge back with sel_all

sel_all_mlp_tog<-merge(sel_all_mlp, sel_all[,-2], by=c("date", "ID", "eggmass", "original"))


table(sel_all_mlp_tog$eggmass_type,sel_all_mlp_tog$original:sel_all_mlp_tog$gen)

# 20C16ppt:F1 20C16ppt:p 20C16ppt:F2 20C16ppt:F3 20C32ppt:F1 20C32ppt:p 20C32ppt:F2 20C32ppt:F3
# l         193        207          33           0         440        630         123           0
# m          80        128           4           0          82        233          44           0
# p         408       1824           2           4         821       2527          75           4


#divide by region

sel_all_mlp_tog$region<-ifelse(grepl("SC*", sel_all_mlp_tog$ID), "Long Beach", ifelse(grepl("W128|13[0-4]", sel_all_mlp_tog$ID), "Tomales", "Mill Valley"))


#do for the second experiment 

F1_sal<-read.table("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/SoCALvsNorCAL/F1s/Capsule_data.csv", header=TRUE, sep=",")
SvN<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/SoCALvsNorCAL/SoCalvNorCal_Data.csv", header=TRUE, sep=",")


head(F1_sal)
head(SvN)

colnames(F1_sal)<-c("date", "order", "Capsule", "mom_id", "mom_env", "ID", "eggmass", "original", "Family")
colnames(SvN)<-c("date", "order", "Capsule", "Family", "ID", "eggmass", "original", "mom_env", "region", "starved")

#add in generation
F1_sal$generation<-ifelse(F1_sal$mom_env == "20C16ppt", "low selected", "high selected")
SvN$generation<-as.vector(strrep("Parental", 1))

cols2<-intersect(colnames(F1_sal[,-2]), colnames(SvN[,-2]))
svn_sel_all<-rbind(F1_sal[,cols2], SvN[,cols2])

head(svn_sel_all)

svn_sel_all_mpl<-aggregate(svn_sel_all$Capsule, list(svn_sel_all$date, svn_sel_all$eggmass, svn_sel_all$ID, svn_sel_all$original, svn_sel_all$generation), function(x) ifelse(x > 0.15, "l", "p"))


colnames(svn_sel_all_mpl)<-c("date", "eggmass", "ID", "original", "generation", "x")

#we now want to ID the groups that have l and p, these are our mixed eggmasses
#can do this in two parts, first seeing which are l

svn_sel_all_mpl$is_l<-grepl("l", svn_sel_all_mpl$x)

#now which are p

svn_sel_all_mpl$is_p<-grepl("p", svn_sel_all_mpl$x)


#and now which are true for both l and p

svn_sel_all_mpl$is_m<- svn_sel_all_mpl$is_l == T & svn_sel_all_mpl$is_p == T

table(svn_sel_all_mpl$is_m)


#now make a vector that actually identifies each eggmass as l, p, or m


svn_sel_all_mpl$eggmass_type<-ifelse(svn_sel_all_mpl$is_m == T, "m", ifelse(svn_sel_all_mpl$is_l == T, "l", "p"))

table(svn_sel_all_mpl$eggmass_type)

#success

#look to see how many egg masses were laid in high versus low salinity

table(svn_sel_all_mpl$eggmass_type,svn_sel_all_mpl$original)
#   20C16ppt 20C32ppt
#l      122      268
#m       45       64
#p      411      581

#to get high and low salinity by generation we will ahve to merge back with sel_all
parental<-svn_sel_all_mpl[svn_sel_all_mpl$generation == "Parental",]

table(parental$eggmass_type, parental$original)
# 20C16ppt 20C32ppt
# l       75      186
# m       26       51
# p      328      502

kids_low<-svn_sel_all_mpl[svn_sel_all_mpl$generation == "low selected",]
kids_high<-svn_sel_all_mpl[svn_sel_all_mpl$generation == "high selected",]

table(kids_low$eggmass_type, kids_low$original)

# 20C16ppt 20C32ppt
# l       14       24
# m        8        8
# p       44       53

table(kids_high$eggmass_type, kids_high$original)

# 20C16ppt 20C32ppt
# l       33       58
# m       11        5
# p       39       26



####### ID trait lability ##############



lab<-table(sel_all_mlp$ID, sel_all_mlp$eggmass_type)

length(lab)
#1884

#first how many only laid l?

#which are only l?
l<-lab[c(which(lab[,2] == 0 & lab[,3] == 0)),]
length(lab[c(which(lab[,2] == 0 & lab[,3] == 0)),])
#192

#and which are only p?
p<-lab[c(which(lab[,1] == 0 & lab[,2] == 0)),]
length(lab[c(which(lab[,1] == 0 & lab[,2] == 0)),])
#1080


#are there any that laid only mixed?
m<-lab[c(which(lab[,1] == 0 & lab[,3] == 0)),]
length(lab[c(which(lab[,1] == 0 & lab[,3] == 0)),])
#21

#add those together to get the number that laid the same kind of egg mass

(192+1080+21)/length(lab)
#[1] 0.6863057


#which made lecitho and plankto
lp<-lab[c(which(lab[,1] > 0 & lab[,3] > 0 & lab[,2] == 0)),]

length(lp)
#234

#which made lecitho and mixed?
lm<-lab[c(which(lab[,1] > 0 & lab[,2] > 0 & lab[,3] == 0)),]
length(lm)
#63


#plankto and mixed
pm<-lab[c(which(lab[,2] > 0 & lab[,3] > 0 & lab[,1] == 0)),]
length(pm)
#147

#add em all up
192+1080+21+234+63+147

#test significnace of switching with salinity and maternal family
library(data.table)
#score the slugs that did not switch with a 0 and those that did with a 1
switcher<-setDT(as.data.frame.matrix(rbind(lm, lp, pm)), keep.rownames="ID")
switcher$s<-rep(1)

#get the non-switchers
const<-setDT(as.data.frame.matrix(rbind(l,m, p)), keep.rownames = "ID")
const$s<-rep(0)

tog<-rbind(switcher, const)

colnames(tog)

#now merge with dataframe that has family and salinity data

test<-unique(merge(tog[,c("ID", "s")], sel_all_mlp_tog[,c("ID", "original", "Family", "gen", "region")], by="ID"))

model_switch<-glm(s~original+gen+region, data=test, family="binomial")
summary(model_switch)

ggplot(data=test, aes(s))+geom_histogram()+facet_grid(gen~original)

# Code released under MIT license
require(ggplot2); require(ggbeeswarm)

# Rare and Declining lists
#####
rar_species = readRDS("/Users/Mark/Desktop/DecBees_FINAL_documents/rare_species.rds")
dec_species = readRDS("/Users/Mark/Desktop/DecBees_FINAL_documents/declining_species.rds")

threat_list = c(rar_species, dec_species)
#####

# Load in data for in-text results
#####
bee_abun = readRDS("/Users/Mark/Desktop/DecBees_FINAL_documents/BeeAbun_ByPlantSp.rds")

# wildflowers with sampling time correction
pr.f = bee_abun[[1]]
pt.f = bee_abun[[2]]
mf.f = bee_abun[[3]]

# crops (sampling time correction NA for crops)
bl.f = bee_abun[[4]]
cr.f = bee_abun[[5]]
wa.f = bee_abun[[6]]

## Abundance

# wildflowers without sampling time correction (useful for summary data)
pr.a2 = bee_abun[[7]]
pt.a2 = bee_abun[[8]]
mf.a2 = bee_abun[[9]]

# wildflowers with sampling time correction
pr.a = bee_abun[[10]]
pt.a = bee_abun[[11]]
mf.a = bee_abun[[12]]

# crops (sampling time correction NA for crops)
bl.a = bee_abun[[13]]
cr.a = bee_abun[[14]]
wa.a = bee_abun[[15]]
#####

# in-text results

# Data summary
#####
# Bee abundance - wildflowers
sum(pr.a2) + sum(pt.a2) + sum(mf.a2)

# Bee richness - wildflowers
asw.a = c(rownames(pr.a), rownames(pt.a), rownames(mf.a))
length(unique(asw.a))

# Bee abundance - crops
sum(bl.a) + sum(cr.a) + sum(wa.a)

# Bee richness - crops
asc.a = union(rownames(bl.a), rownames(cr.a))
asc.a = union(asc.a, rownames(wa.a))
length(unique(asc.a))

# Bee richness - total
length(unique(c(asw.a, asc.a)))

# How many species are rare?
as.a = union(asw.a, asc.a)
as.a.rare = intersect(as.a, rar_species)
length(as.a.rare)

# How many species were declining?
as.a.dec = intersect(as.a, dec_species)
length(as.a.dec)

# How many species were both rare and declining?
length(intersect(as.a.rare, as.a.dec))
#####

# Q1.1 What is the mean abundance of rare and declining bee species across site-years?
#####

# Setup 
# P reptans
pr.a.rd = pr.a[which(rownames(pr.a) %in% threat_list),] # rare + dec
pr.a.r = pr.a[which(rownames(pr.a) %in% rar_species),] # just rare
pr.a.d = pr.a[which(rownames(pr.a) %in% dec_species),] # just declining

pr.Q1a = colSums(pr.a.rd) / colSums(pr.a) # relativize
pr.Q1a.r = colSums(pr.a.r) / colSums(pr.a) # relativize
pr.Q1a.d = colSums(pr.a.d) / colSums(pr.a) # relativize

# P tanacetifolia 
# (comments as above for all plant species)
pt.a.rd = pt.a[which(rownames(pt.a) %in% threat_list),]
pt.a.r = pt.a[which(rownames(pt.a) %in% rar_species),]
pt.a.d = pt.a[which(rownames(pt.a) %in% dec_species),]

pt.Q1a = colSums(pt.a.rd) / colSums(pt.a)
pt.Q1a.r = colSums(pt.a.r) / colSums(pt.a)
pt.Q1a.d = colSums(pt.a.d) / colSums(pt.a)

# M fistulosa
mf.a.rd = mf.a[which(rownames(mf.a) %in% threat_list),]
mf.a.r = mf.a[which(rownames(mf.a) %in% rar_species),]
mf.a.d = mf.a[which(rownames(mf.a) %in% dec_species),]

mf.Q1a = colSums(mf.a.rd) / colSums(mf.a)
mf.Q1a.r = colSums(mf.a.r) / colSums(mf.a)
mf.Q1a.d = colSums(mf.a.d) / colSums(mf.a)

# blueberry
bl.a.rd = bl.a[which(rownames(bl.a) %in% threat_list),]
bl.a.r = bl.a[which(rownames(bl.a) %in% rar_species),]
bl.a.d = bl.a[which(rownames(bl.a) %in% dec_species),]

bl.Q1a = colSums(bl.a.rd) / colSums(bl.a)
bl.Q1a.r = colSums(bl.a.r) / colSums(bl.a)
bl.Q1a.d = colSums(bl.a.d) / colSums(bl.a)

# cranberry
cr.a.rd = cr.a[which(rownames(cr.a) %in% threat_list),]
cr.a.r = cr.a[which(rownames(cr.a) %in% rar_species),]
cr.a.d = cr.a[which(rownames(cr.a) %in% dec_species),]

cr.Q1a = colSums(cr.a.rd) / colSums(cr.a)
cr.Q1a.r = colSums(cr.a.r) / colSums(cr.a)
cr.Q1a.d = colSums(cr.a.d) / colSums(cr.a)

# watermelon
wa.a.rd = wa.a[which(rownames(wa.a) %in% threat_list),]
wa.a.r = wa.a[which(rownames(wa.a) %in% rar_species),]
wa.a.d = wa.a[which(rownames(wa.a) %in% dec_species),]

wa.Q1a = colSums(wa.a.rd) / colSums(wa.a)
wa.Q1a.r = colSums(wa.a.r) / colSums(wa.a)
wa.Q1a.d = colSums(wa.a.d) / colSums(wa.a)

# combine
as.Q1a = c(pr.Q1a, pt.Q1a, mf.Q1a, bl.Q1a, cr.Q1a, wa.Q1a)
as.Q1a.r = c(pr.Q1a.r, pt.Q1a.r, mf.Q1a.r, bl.Q1a.r, cr.Q1a.r, wa.Q1a.r)
as.Q1a.d = c(pr.Q1a.d, pt.Q1a.d, mf.Q1a.d, bl.Q1a.d, cr.Q1a.d, wa.Q1a.d)

# What is the mean abundance of rare and declining bee species across site-years?
mean(as.Q1a)
Q1a.ci = t.test(as.Q1a, conf.level = 0.95)
c(Q1a.ci$estimate, Q1a.ci$conf.int)

# as above but just wildflowers
as.Q1a.w = as.Q1a[1:144]; mean(as.Q1a.w) 
Q1a.w.ci = t.test(as.Q1a.w, conf.level = 0.95)
c(Q1a.w.ci$estimate, Q1a.w.ci$conf.int)

# as above but just crops
as.Q1a.c = as.Q1a[145:240]; mean(as.Q1a.c) 
Q1a.c.ci = t.test(as.Q1a.c, conf.level = 0.95)
c(Q1a.c.ci$estimate, Q1a.c.ci$conf.int)

# as above but just rare species
mean(as.Q1a.r) 
Q1a.ci.r = t.test(as.Q1a.r, conf.level = 0.95)
c(Q1a.ci.r$estimate, Q1a.ci.r$conf.int)

# as above but just declining species
mean(as.Q1a.d) 
Q1a.ci.d = t.test(as.Q1a.d, conf.level = 0.95)
c(Q1a.ci.d$estimate, Q1a.ci.d$conf.int)
#####

# Q1.2 What is the mean function of rare and declining bee species across site-years?
#####

# Setup
# (comments as for abundance)
# P reptans
pr.f.rd = pr.f[which(rownames(pr.f) %in% threat_list),]
pr.f.r = pr.f[which(rownames(pr.f) %in% rar_species),]
pr.f.d = pr.f[which(rownames(pr.f) %in% dec_species),]

pr.Q1f = colSums(pr.f.rd) / colSums(pr.f)
pr.Q1f.r = colSums(pr.f.r) / colSums(pr.f)
pr.Q1f.d = colSums(pr.f.d) / colSums(pr.f)

# P tanacetifolia
pt.f.rd = pt.f[which(rownames(pt.f) %in% threat_list),]
pt.f.r = pt.f[which(rownames(pt.f) %in% rar_species),]
pt.f.d = pt.f[which(rownames(pt.f) %in% dec_species),]

pt.Q1f = colSums(pt.f.rd) / colSums(pt.f)
pt.Q1f.r = colSums(pt.f.r) / colSums(pt.f)
pt.Q1f.d = colSums(pt.f.d) / colSums(pt.f)

# M fistulosa
mf.f.rd = mf.f[which(rownames(mf.f) %in% threat_list),]
mf.f.r = mf.f[which(rownames(mf.f) %in% rar_species),]
mf.f.d = mf.f[which(rownames(mf.f) %in% dec_species),]

mf.Q1f = colSums(mf.f.rd) / colSums(mf.f)
mf.Q1f.r = colSums(mf.f.r) / colSums(mf.f)
mf.Q1f.d = colSums(mf.f.d) / colSums(mf.f)

# Blueberry
bl.f.rd = bl.f[which(rownames(bl.f) %in% threat_list),]
bl.f.r = bl.f[which(rownames(bl.f) %in% rar_species),]
bl.f.d = bl.f[which(rownames(bl.f) %in% dec_species),]

bl.Q1f = colSums(bl.f.rd) / colSums(bl.f)
bl.Q1f.r = colSums(bl.f.r) / colSums(bl.f)
bl.Q1f.d = colSums(bl.f.d) / colSums(bl.f)

# Cranberry
cr.f.rd = cr.f[which(rownames(cr.f) %in% threat_list),]
cr.f.r = cr.f[which(rownames(cr.f) %in% rar_species),]
cr.f.d = cr.f[which(rownames(cr.f) %in% dec_species),]

cr.Q1f = colSums(cr.f.rd) / colSums(cr.f)
cr.Q1f.r = colSums(cr.f.r) / colSums(cr.f)
cr.Q1f.d = colSums(cr.f.d) / colSums(cr.f)

# Watermelon
wa.f.rd = wa.f[which(rownames(wa.f) %in% threat_list),]
wa.f.r = wa.f[which(rownames(wa.f) %in% rar_species),]
wa.f.d = wa.f[which(rownames(wa.f) %in% dec_species),]

wa.Q1f = colSums(wa.f.rd) / colSums(wa.f)
wa.Q1f.r = colSums(wa.f.r) / colSums(wa.f)
wa.Q1f.d = colSums(wa.f.d) / colSums(wa.f)

# combine
as.Q1f = c(pr.Q1f, pt.Q1f, mf.Q1f, bl.Q1f, cr.Q1f, wa.Q1f)
as.Q1f.r = c(pr.Q1f.r, pt.Q1f.r, mf.Q1f.r, bl.Q1f.r, cr.Q1f.r, wa.Q1f.r)
as.Q1f.d = c(pr.Q1f.d, pt.Q1f.d, mf.Q1f.d, bl.Q1f.d, cr.Q1f.d, wa.Q1f.d)

# What is the mean function of rare and declining bee species across site-years?
mean(as.Q1f)
Q1f.ci = t.test(as.Q1f, conf.level = 0.95)
c(Q1f.ci$estimate, Q1f.ci$conf.int)

# as above but just wildflowers
as.Q1f.w = as.Q1f[1:144]; mean(as.Q1f.w)
Q1f.w.ci = t.test(as.Q1f.w, conf.level = 0.95)
c(Q1f.w.ci$estimate, Q1f.w.ci$conf.int)

# as above but just crops
as.Q1f.c = as.Q1f[145:240]; mean(as.Q1f.c)
Q1f.c.ci = t.test(as.Q1f.c, conf.level = 0.95)
c(Q1f.c.ci$estimate, Q1f.c.ci$conf.int)

# as above but just rare species
mean(as.Q1f.r)
Q1f.ci.r = t.test(as.Q1f.r, conf.level = 0.95)
c(Q1f.ci.r$estimate, Q1f.ci.r$conf.int)

# as above but just declining species
mean(as.Q1f.d)
Q1f.ci.d = t.test(as.Q1f.d, conf.level = 0.95)
c(Q1f.ci.d$estimate, Q1f.ci.d$conf.int)
#####

# Q2 How many rare/dec species are important pollinators when the data are summed?

# P reptans
q2.pr = rev(sort(rowSums(pr.f)))
q2.pr = cumsum(q2.pr)/sum(q2.pr)
q2.pr = names(q2.pr[1:min(which(q2.pr>0.5))])

# PT
q2.pt = rev(sort(rowSums(pt.f)))
q2.pt = cumsum(q2.pt)/sum(q2.pt)
q2.pt = names(q2.pt[1:min(which(q2.pt>0.5))])

# MF
q2.mf = rev(sort(rowSums(mf.f)))
q2.mf = cumsum(q2.mf)/sum(q2.mf)
q2.mf = names(q2.mf[1:min(which(q2.mf>0.5))])

# Blueberry
q2.bl = rev(sort(rowSums(bl.f)))
q2.bl = cumsum(q2.bl)/sum(q2.bl)
q2.bl = names(q2.bl[1:min(which(q2.bl>0.5))])

# Cranberry
q2.cr = rev(sort(rowSums(cr.f)))
q2.cr = cumsum(q2.cr)/sum(q2.cr)
q2.cr = names(q2.cr[1:min(which(q2.cr>0.5))])

# Watermelon
q2.wa = rev(sort(rowSums(wa.f)))
q2.wa = cumsum(q2.wa)/sum(q2.wa)
q2.wa = names(q2.wa[1:min(which(q2.wa>0.5))])

# aggregated min set size, all species (not used)
q2.as = c(length(q2.pr), length(q2.pt), length(q2.mf), 
          length(q2.bl), length(q2.cr), length(q2.mf))
names(q2.as) = c("PR", "PT", "MF", "BL", "CR", "WA")
q2.as

# aggregated min set size, RD species
q2.rd = c(length(q2.pr[which(q2.pr %in% threat_list)]), 
          length(q2.pt[which(q2.pt %in% threat_list)]), 
          length(q2.mf[which(q2.mf %in% threat_list)]), 
          length(q2.bl[which(q2.bl %in% threat_list)]), 
          length(q2.cr[which(q2.cr %in% threat_list)]), 
          length(q2.wa[which(q2.wa %in% threat_list)]))
names(q2.rd) = c("PR", "PT", "MF", "BL", "CR", "WA")

# How many rare/dec species are important pollinators for >=1 plant sp. when the data are summed? 
mean(q2.rd)

# How many of these are unique?
unique.sp = unique(c(q2.pr, q2.pt, q2.mf, q2.bl, q2.cr, q2.wa)) # all important species
unique.rd.sp = unique.sp[which(unique.sp %in% threat_list)] # just rare/dec important sp
length(unique.rd.sp) 

#####

# 3. How many rare/dec species are important when the data are analyzed 
# for each site and year separately?
#####
# output of the genetic optimizer
q3.all = read.csv("/Users/Mark/Desktop/DecBees_FINAL_documents/Fig2_Data.csv")

# the output includes more rows than we need
q3 = q3.all[which(q3.all$threshold==0.5),] # just the 50% threshold

# choose max # of sites for each plant species
q3.pr = q3[which(q3$plant=="Polemonium" & q3$nsites==24),]$decl_plus_rare
q3.pt = q3[which(q3$plant=="Phacelia" & q3$nsites==24),]$decl_plus_rare
q3.mf = q3[which(q3$plant=="Monarda" & q3$nsites==24),]$decl_plus_rare
q3.bl = q3[which(q3$plant=="Blueberry" & q3$nsites==16),]$decl_plus_rare
q3.cr = q3[which(q3$plant=="Cranberry" & q3$nsites==16),]$decl_plus_rare
q3.wa = q3[which(q3$plant=="Watermelon" & q3$nsites==16),]$decl_plus_rare

# combine so we can take a mean
q3.vec = c(q3.pr, q3.pt, q3.mf, q3.bl, q3.cr, q3.wa)

# Average rare/dec important bee species (average taken across plant species)
mean(q3.vec)

minset.unique = read.csv("/Users/Mark/Desktop/DecBees_FINAL_documents/MinSet_Unique.csv")

# How many unique, important, rare/dec species across ALL plant species
length(unique(minset.unique$as[which(minset.unique$as %in% threat_list)]))
length(unique(minset.unique$as[which(minset.unique$as %in% rar_species)]))# rare
length(unique(minset.unique$as[which(minset.unique$as %in% dec_species)])) # declining
length(unique(minset.unique$as[which(minset.unique$as %in% dec_species &
                                       minset.unique$as %in% rar_species)])) # both

# one species is both rare and declining, which is why 25 + 8 != 32

# for each rare/dec bee species in the minimum set, find out how many plant species
# needed that rare/dec species

q3.part2 = summary(as.factor(unlist(
  apply(minset.unique[,2:7], 2, function(x){x[which(x %in% threat_list)]}))))

length(q3.part2) # 32, as above
length(which(q3.part2==1)) # how many rare/dec bee species important to just one plant sp

#####

# setup Fig 1
#####
fig1.df = readRDS("/Users/Mark/Desktop/DecBees_FINAL_documents/Fig1_Data.rds")
box.colors = c("goldenrod", "skyblue")
mytheme = readRDS("/Users/Mark/Desktop/DecBees_FINAL_documents/mytheme.rds")
#####

# Figure 1a - abundance
#####
ggplot(fig1.df, aes(plant_sp, pct_abun*100)) +
  
  geom_quasirandom(fill=fig1.df$color, size=5, alpha=0.4, shape=21, stroke=0.2, color="black") +
  geom_boxplot(color=c("gray30", rep("goldenrod2", 3), rep("skyblue2", 3)), 
               fill=NA, lwd=1.3, alpha=0.8) +
  stat_summary(fun=mean, geom="point", shape=21, size=5, stroke=1,
               fill=c("gray30", rep(box.colors[1], 3), rep(box.colors[2], 3))) +
  scale_y_continuous(name="Percent Abundance", limits=c(0,100)) +
  scale_x_discrete(name="", labels=c("All", "P. rep", "P. tan", "M. fis", "Blue", "Cran", "Wat")) +
  
  geom_hline(yintercept = 25, lwd=0.5, alpha=0.5, lty=2, color="black") +
  geom_hline(yintercept = 50, lwd=0.5, alpha=0.5, lty=2, color="black") +
  geom_hline(yintercept = 75, lwd=0.5, alpha=0.5, lty=2, color="black") +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  mytheme +
  theme(legend.position = "none")
#####

# Figure 1b - function
#####
ggplot(fig1.df, aes(plant_sp, pct_fn*100)) +
  
  geom_quasirandom(fill=fig1.df$color, size=5, alpha=0.4, shape=21, stroke=0.2, color="black") +
  geom_boxplot(color=c("gray30", rep("goldenrod3", 3), rep("skyblue3", 3)), 
               fill=NA, lwd=1.3, alpha=0.8) +
  stat_summary(fun=mean, geom="point", shape=21, size=5, stroke=1,
               fill=c("gray30", rep(box.colors[1], 3), rep(box.colors[2], 3))) +
  scale_y_continuous(name="Percent Pollination", limits=c(0,100)) +
  scale_x_discrete(name="", labels=c("All", "P. rep", "P. tan", "M. fis", "Blue", "Cran", "Wat")) +
  
  geom_hline(yintercept = 25, lwd=0.5, alpha=0.5, lty=2, color="black") +
  geom_hline(yintercept = 50, lwd=0.5, alpha=0.5, lty=2, color="black") +
  geom_hline(yintercept = 75, lwd=0.5, alpha=0.5, lty=2, color="black") +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  mytheme +
  theme(legend.position = "none")
#####

# Fig 2a - P reptans
#####
# renaming
fig2 = q3.all

fig2$pct.rd = fig2$decl_plus_rare/fig2$all

library(colorspace)
fig2.colors = diverge_hcl(12)

fig2a = fig2[which(fig2$plant=="Polemonium" & fig2$threshold==0.5),-1]
fig2a = fig2a[-which(fig2a$nsites==2),]

#fig2a = rbind(fig2a, fig2a[7,])
#fig2a[8,c(1,4:11)] = c(30,5,0,1,0,0,0,1,0)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  #geom_line(aes(x=c(1,24.7), y=rep(q2.rd[1] + 0.20, 2)), lty=2, lwd=1.3, col="purple", alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.r[1]        , 2)), lty=2, lwd=1.3, col="blue"  , alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.d[1]        , 2)), lty=2, lwd=1.3, col="red"   , alpha=0.7) +
  
  geom_errorbar(aes(x=fig2a$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, fig2a$decl_plus_rare-(1.96*fig2a$decl_plus_rare_sd)), 
                    ymax=fig2a$decl_plus_rare+(1.96*fig2a$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2a$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, fig2a$decl-(1.96*fig2a$decl_sd)), 
                    ymax=fig2a$decl+(1.96*fig2a$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2a$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, fig2a$rare-(1.96*fig2a$rare_sd)), 
                    ymax=fig2a$rare+(1.96*fig2a$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=fig2a$nsites + c(rep(-0.3, 5), 0), y=fig2a$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2a$nsites + c(rep( 0.3, 5), 0), y=fig2a$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2a$nsites + c(rep( 0.0, 5), 0), y=fig2a$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=fig2a$nsites + c(rep(-0.3, 5), 0), y=fig2a$decl_plus_rare), 
             color="black", fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2a$nsites + c(rep( 0.3, 5), 0), y=fig2a$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2a$nsites + c(rep( 0.0, 5), 0), y=fig2a$rare), 
             color="black", fill="skyblue", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,15))

#####

# Fig 2b - P tanacetifolia
#####
fig2b = fig2[which(fig2$plant=="Phacelia" & fig2$threshold==0.5),-1]
fig2b = fig2b[-which(fig2b$nsites==2),]

#fig2b = rbind(fig2b, fig2b[7,])
#fig2b[8,c(1,4:11)] = c(30,5,0,1,0,0,0,1,0)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  #geom_line(aes(x=c(1,24.7), y=rep(q2.rd[1] + 0.20, 2)), lty=2, lwd=1.3, col="purple", alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.r[1]        , 2)), lty=2, lwd=1.3, col="blue"  , alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.d[1]        , 2)), lty=2, lwd=1.3, col="red"   , alpha=0.7) +
  
  geom_errorbar(aes(x=fig2b$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, fig2b$decl_plus_rare-(1.96*fig2b$decl_plus_rare_sd)), 
                    ymax=fig2b$decl_plus_rare+(1.96*fig2b$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2b$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, fig2b$decl-(1.96*fig2b$decl_sd)), 
                    ymax=fig2b$decl+(1.96*fig2b$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2b$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, fig2b$rare-(1.96*fig2b$rare_sd)), 
                    ymax=fig2b$rare+(1.96*fig2b$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=fig2b$nsites + c(rep(-0.3, 5), 0), y=fig2b$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2b$nsites + c(rep( 0.3, 5), 0), y=fig2b$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2b$nsites + c(rep( 0.0, 5), 0), y=fig2b$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=fig2b$nsites + c(rep(-0.3, 5), 0), y=fig2b$decl_plus_rare), 
             color="black", fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2b$nsites + c(rep( 0.3, 5), 0), y=fig2b$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2b$nsites + c(rep( 0.0, 5), 0), y=fig2b$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,15))

#####

# Fig 2c - M fistulosa
#####
fig2c = fig2[which(fig2$plant=="Monarda" & fig2$threshold==0.5),-1]
fig2c = fig2c[-which(fig2c$nsites==2),]

#fig2c = rbind(fig2c, fig2c[7,])
#fig2c[8,c(1,4:11)] = c(30,5,0,1,0,0,0,1,0)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  #geom_line(aes(x=c(1,24.7), y=rep(q2.rd[1] + 0.20, 2)), lty=2, lwd=1.3, col="purple", alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.r[1]        , 2)), lty=2, lwd=1.3, col="blue"  , alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.d[1]        , 2)), lty=2, lwd=1.3, col="red"   , alpha=0.7) +
  
  geom_errorbar(aes(x=fig2c$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, fig2c$decl_plus_rare-(1.96*fig2c$decl_plus_rare_sd)), 
                    ymax=fig2c$decl_plus_rare+(1.96*fig2c$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2c$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, fig2c$decl-(1.96*fig2c$decl_sd)), 
                    ymax=fig2c$decl+(1.96*fig2c$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2c$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, fig2c$rare-(1.96*fig2c$rare_sd)), 
                    ymax=fig2c$rare+(1.96*fig2c$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=fig2c$nsites + c(rep(-0.3, 5), 0), y=fig2c$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2c$nsites + c(rep( 0.3, 5), 0), y=fig2c$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2c$nsites + c(rep( 0.0, 5), 0), y=fig2c$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=fig2c$nsites + c(rep(-0.3, 5), 0), y=fig2c$decl_plus_rare), 
             color="black", fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2c$nsites + c(rep( 0.3, 5), 0), y=fig2c$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2c$nsites + c(rep( 0.0, 5), 0), y=fig2c$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,15))

#####

# Fig 2d - Blueberry
#####
fig2d = fig2[which(fig2$plant=="Blueberry" & fig2$threshold==0.5),-1]
fig2d = fig2d[-which(fig2d$nsites==2),]

#fig2d = rbind(fig2d, fig2d[7,])
#fig2d[8,c(1,4:11)] = c(30,5,0,1,0,0,0,1,0)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  #geom_line(aes(x=c(1,24.7), y=rep(q2.rd[1] + 0.20, 2)), lty=2, lwd=1.3, col="purple", alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.r[1]        , 2)), lty=2, lwd=1.3, col="blue"  , alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.d[1]        , 2)), lty=2, lwd=1.3, col="red"   , alpha=0.7) +
  
  geom_errorbar(aes(x=fig2d$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, fig2d$decl_plus_rare-(1.96*fig2d$decl_plus_rare_sd)), 
                    ymax=fig2d$decl_plus_rare+(1.96*fig2d$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2d$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, fig2d$decl-(1.96*fig2d$decl_sd)), 
                    ymax=fig2d$decl+(1.96*fig2d$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2d$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, fig2d$rare-(1.96*fig2d$rare_sd)), 
                    ymax=fig2d$rare+(1.96*fig2d$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=fig2d$nsites + c(rep(-0.3, 4), 0), y=fig2d$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2d$nsites + c(rep( 0.3, 4), 0), y=fig2d$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2d$nsites + c(rep( 0.0, 4), 0), y=fig2d$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=fig2d$nsites + c(rep(-0.3, 4), 0), y=fig2d$decl_plus_rare), 
             color="black", fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2d$nsites + c(rep( 0.3, 4), 0), y=fig2d$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2d$nsites + c(rep( 0.0, 4), 0), y=fig2d$rare), 
             color="black", fill="skyblue", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,15))

#####

# Fig 2e - Cranberry
#####
fig2e = fig2[which(fig2$plant=="Cranberry" & fig2$threshold==0.5),-1]
fig2e = fig2e[-which(fig2e$nsites==2),]

#fig2e = rbind(fig2e, fig2e[7,])
#fig2e[8,c(1,4:11)] = c(30,5,0,1,0,0,0,1,0)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  #geom_line(aes(x=c(1,24.7), y=rep(q2.rd[1] + 0.20, 2)), lty=2, lwd=1.3, col="purple", alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.r[1]        , 2)), lty=2, lwd=1.3, col="blue"  , alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.d[1]        , 2)), lty=2, lwd=1.3, col="red"   , alpha=0.7) +
  
  geom_errorbar(aes(x=fig2e$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, fig2e$decl_plus_rare-(1.96*fig2e$decl_plus_rare_sd)), 
                    ymax=fig2e$decl_plus_rare+(1.96*fig2e$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2e$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, fig2e$decl-(1.96*fig2e$decl_sd)), 
                    ymax=fig2e$decl+(1.96*fig2e$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=fig2e$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, fig2e$rare-(1.96*fig2e$rare_sd)), 
                    ymax=fig2e$rare+(1.96*fig2e$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=fig2e$nsites + c(rep(-0.3, 4), 0), y=fig2e$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2e$nsites + c(rep( 0.3, 4), 0), y=fig2e$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2e$nsites + c(rep( 0.0, 4), 0), y=fig2e$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=fig2e$nsites + c(rep(-0.3, 4), 0), y=fig2e$decl_plus_rare), 
             color="black", fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2e$nsites + c(rep( 0.3, 4), 0), y=fig2e$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2e$nsites + c(rep( 0.0, 4), 0), y=fig2e$rare), 
             color="black", fill="skyblue", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,15))

#####

# Fig 2f - Watermelon
#####
fig2f = fig2[which(fig2$plant=="Watermelon" & fig2$threshold==0.5),-1]
fig2f = fig2f[-which(fig2f$nsites==2),]

#fig2f = rbind(fig2f, fig2f[7,])
#fig2f[8,c(1,4:11)] = c(30,5,0,1,0,0,0,1,0)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  #geom_line(aes(x=c(1,24.7), y=rep(q2.rd[1] + 0.20, 2)), lty=2, lwd=1.3, col="purple", alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.r[1]        , 2)), lty=2, lwd=1.3, col="blue"  , alpha=0.7) +
  #geom_line(aes(x=c(1,24.7), y=rep(q2.d[1]        , 2)), lty=2, lwd=1.3, col="red"   , alpha=0.7) +
  
  geom_errorbar(aes(x=fig2f$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, fig2f$decl_plus_rare-(1.96*fig2f$decl_plus_rare_sd)), 
                    ymax=fig2f$decl_plus_rare+(1.96*fig2f$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7) +
  geom_errorbar(aes(x=fig2f$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, fig2f$decl-(1.96*fig2f$decl_sd)), 
                    ymax=fig2f$decl+(1.96*fig2f$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7) +
  geom_errorbar(aes(x=fig2f$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, fig2f$rare-(1.96*fig2f$rare_sd)), 
                    ymax=fig2f$rare+(1.96*fig2f$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7) +
  
  geom_line(aes(x=fig2f$nsites + c(rep(-0.3, 4), 0), y=fig2f$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2f$nsites + c(rep( 0.3, 4), 0), y=fig2f$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=fig2f$nsites + c(rep( 0.0, 4), 0), y=fig2f$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=fig2f$nsites + c(rep(-0.3, 4), 0), y=fig2f$decl_plus_rare), 
             color="black", fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2f$nsites + c(rep( 0.3, 4), 0), y=fig2f$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=fig2f$nsites + c(rep( 0.0, 4), 0), y=fig2f$rare), 
             color="black", fill="skyblue", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,15))

#####

# Fig3 - All six plant species
#####
fig3 = read.csv("/Users/Mark/Desktop/DecBees_FINAL_documents/Fig3_Data.csv")

f3.list = vector("list", 6)

for(i in 1:6){
  f3.list[[i]] = fig3[,(i+1)]
  f3.list[[i]] = f3.list[[i]][-which(f3.list[[i]]=="")]
}

fig3.out = matrix(NA, 6, 10000)
fig3.rd.out = matrix(NA, 6, 10000)
fig3.r.out = matrix(NA, 6, 10000)
fig3.d.out = matrix(NA, 6, 10000)

for(j in 1:6){
  for(i in 1:10000){
    
    fig3.samp = sample(1:6, j) 
    hold = unique(unlist(f3.list[fig3.samp]))
    fig3.out[j,i] = length(hold)
    
    fig3.rd.out[j,i] = length(hold[which(hold %in% threat_list)])
    fig3.r.out[j,i] = length(hold[which(hold %in% rar_species)])
    fig3.d.out[j,i] = length(hold[which(hold %in% dec_species)])
    
  }
}

fig3.rd.final = rowMeans(fig3.rd.out)
fig3.rd.sd = apply(fig3.rd.out, 1, sd)

fig3.r.final = rowMeans(fig3.r.out)
fig3.r.sd = apply(fig3.r.out, 1, sd)

fig3.d.final = rowMeans(fig3.d.out)
fig3.d.sd = apply(fig3.d.out, 1, sd)

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept=10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept=30), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=1:6, ymin=pmax(0, fig3.rd.final-(1.96*fig3.rd.sd)), 
                    ymax=fig3.rd.final+(1.96*fig3.rd.sd)), 
                color="gray20", width=0.25, lwd=0.8, alpha=0.7) +
  geom_errorbar(aes(x=1:6+c(rep(0.12,5),0), ymin=pmax(0, fig3.r.final -(1.96*fig3.r.sd)),  
                    ymax=fig3.r.final +(1.96*fig3.r.sd )), 
                color="skyblue3", width=0.25, lwd=0.8, alpha=0.7) +
  geom_errorbar(aes(x=1:6+c(rep(0.24,5),0), ymin=pmax(0, fig3.d.final -(1.96*fig3.d.sd)),  
                    ymax=fig3.d.final +(1.96*fig3.d.sd )), 
                color="goldenrod3", width=0.25, lwd=0.8, alpha=0.7) +
  
  geom_line(aes(x=1:6, y=fig3.rd.final), color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=1:6+c(rep(0.12,5),0), y=fig3.r.final), color="skyblue3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=1:6+c(rep(0.24,5),0), y=fig3.d.final), color="goldenrod3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=1:6, y=fig3.rd.final), color="black", fill="gray30", 
             size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=1:6+c(rep(0.12,5),0), y=fig3.r.final), color="black",fill="skyblue", 
             size=4, pch=23, stroke=1, alpha=0.5) +
  geom_point(aes(x=1:6+c(rep(0.24,5),0), y=fig3.d.final), color="black",fill="goldenrod", 
             size=4, pch=22, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Plant Species", breaks = c(1:6)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,35))
#####

# Supplementary  

# Fig S5 25% 
#####
# panel A
figS2.25.a = fig2[which(fig2$plant=="Polemonium" & fig2$threshold==0.25),-1]
figS2.25.a = figS2.25.a[-which(figS2.25.a$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 2), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 4), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 6), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.25.a$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, figS2.25.a$decl_plus_rare-(1.96*figS2.25.a$decl_plus_rare_sd)), 
                    ymax=figS2.25.a$decl_plus_rare+(1.96*figS2.25.a$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.a$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, figS2.25.a$decl-(1.96*figS2.25.a$decl_sd)), 
                    ymax=figS2.25.a$decl+(1.96*figS2.25.a$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.a$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, figS2.25.a$rare-(1.96*figS2.25.a$rare_sd)), 
                    ymax=figS2.25.a$rare+(1.96*figS2.25.a$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.25.a$nsites + c(rep(-0.3, 5), 0), y=figS2.25.a$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.a$nsites + c(rep( 0.3, 5), 0), y=figS2.25.a$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.a$nsites + c(rep( 0.0, 5), 0), y=figS2.25.a$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.25.a$nsites + c(rep(-0.3, 5), 0), y=figS2.25.a$decl_plus_rare), 
             color="black",fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.a$nsites + c(rep( 0.3, 5), 0), y=figS2.25.a$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.a$nsites + c(rep( 0.0, 5), 0), y=figS2.25.a$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,7.5))


# panel B
figS2.25.b = fig2[which(fig2$plant=="Phacelia" & fig2$threshold==0.25),-1]
figS2.25.b = figS2.25.b[-which(figS2.25.b$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 2), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 4), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 6), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.25.b$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, figS2.25.b$decl_plus_rare-(1.96*figS2.25.b$decl_plus_rare_sd)), 
                    ymax=figS2.25.b$decl_plus_rare+(1.96*figS2.25.b$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.b$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, figS2.25.b$decl-(1.96*figS2.25.b$decl_sd)), 
                    ymax=figS2.25.b$decl+(1.96*figS2.25.b$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.b$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, figS2.25.b$rare-(1.96*figS2.25.b$rare_sd)), 
                    ymax=figS2.25.b$rare+(1.96*figS2.25.b$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.25.b$nsites + c(rep(-0.3, 5), 0), y=figS2.25.b$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.b$nsites + c(rep( 0.3, 5), 0), y=figS2.25.b$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.b$nsites + c(rep( 0.0, 5), 0), y=figS2.25.b$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.25.b$nsites + c(rep(-0.3, 5), 0), y=figS2.25.b$decl_plus_rare), 
             color="black",fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.b$nsites + c(rep( 0.3, 5), 0), y=figS2.25.b$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.b$nsites + c(rep( 0.0, 5), 0), y=figS2.25.b$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,7.5))

# panel C
figS2.25.c = fig2[which(fig2$plant=="Monarda" & fig2$threshold==0.25),-1]
figS2.25.c = figS2.25.c[-which(figS2.25.c$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 2), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 4), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 6), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.25.c$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, figS2.25.c$decl_plus_rare-(1.96*figS2.25.c$decl_plus_rare_sd)), 
                    ymax=figS2.25.c$decl_plus_rare+(1.96*figS2.25.c$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.c$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, figS2.25.c$decl-(1.96*figS2.25.c$decl_sd)), 
                    ymax=figS2.25.c$decl+(1.96*figS2.25.c$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.c$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, figS2.25.c$rare-(1.96*figS2.25.c$rare_sd)), 
                    ymax=figS2.25.c$rare+(1.96*figS2.25.c$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.25.c$nsites + c(rep(-0.3, 5), 0), y=figS2.25.c$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.c$nsites + c(rep( 0.3, 5), 0), y=figS2.25.c$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.c$nsites + c(rep( 0.0, 5), 0), y=figS2.25.c$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.25.c$nsites + c(rep(-0.3, 5), 0), y=figS2.25.c$decl_plus_rare), 
             color="black",fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.c$nsites + c(rep( 0.3, 5), 0), y=figS2.25.c$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.c$nsites + c(rep( 0.0, 5), 0), y=figS2.25.c$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,7.5))


# panel D
figS2.25.d = fig2[which(fig2$plant=="Blueberry" & fig2$threshold==0.25),-1]
figS2.25.d = figS2.25.d[-which(figS2.25.d$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 2), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 4), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 6), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.25.d$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, figS2.25.d$decl_plus_rare-(1.96*figS2.25.d$decl_plus_rare_sd)), 
                    ymax=figS2.25.d$decl_plus_rare+(1.96*figS2.25.d$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.d$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, figS2.25.d$decl-(1.96*figS2.25.d$decl_sd)), 
                    ymax=figS2.25.d$decl+(1.96*figS2.25.d$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.d$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, figS2.25.d$rare-(1.96*figS2.25.d$rare_sd)), 
                    ymax=figS2.25.d$rare+(1.96*figS2.25.d$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.25.d$nsites + c(rep(-0.3, 4), 0), y=figS2.25.d$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.d$nsites + c(rep( 0.3, 4), 0), y=figS2.25.d$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.d$nsites + c(rep( 0.0, 4), 0), y=figS2.25.d$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.25.d$nsites + c(rep(-0.3, 4), 0), y=figS2.25.d$decl_plus_rare), 
             color="black", fill="gray20", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.d$nsites + c(rep( 0.3, 4), 0), y=figS2.25.d$decl), 
             color="black", fill="goldenrod3", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.d$nsites + c(rep( 0.0, 4), 0), y=figS2.25.d$rare), 
             color="black", fill="skyblue3", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,7.5))


# panel E
figS2.25.e = fig2[which(fig2$plant=="Cranberry" & fig2$threshold==0.25),-1]
figS2.25.e = figS2.25.e[-which(figS2.25.e$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 2), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 4), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 6), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.25.e$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, figS2.25.e$decl_plus_rare-(1.96*figS2.25.e$decl_plus_rare_sd)), 
                    ymax=figS2.25.e$decl_plus_rare+(1.96*figS2.25.e$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.e$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, figS2.25.e$decl-(1.96*figS2.25.e$decl_sd)), 
                    ymax=figS2.25.e$decl+(1.96*figS2.25.e$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.e$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, figS2.25.e$rare-(1.96*figS2.25.e$rare_sd)), 
                    ymax=figS2.25.e$rare+(1.96*figS2.25.e$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.25.e$nsites + c(rep(-0.3, 4), 0), y=figS2.25.e$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.e$nsites + c(rep( 0.3, 4), 0), y=figS2.25.e$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.e$nsites + c(rep( 0.0, 4), 0), y=figS2.25.e$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.25.e$nsites + c(rep(-0.3, 4), 0), y=figS2.25.e$decl_plus_rare), 
             color="black", fill="gray20", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.e$nsites + c(rep( 0.3, 4), 0), y=figS2.25.e$decl), 
             color="black", fill="goldenrod3", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.e$nsites + c(rep( 0.0, 4), 0), y=figS2.25.e$rare), 
             color="black", fill="skyblue3", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,7.5))

# panel F
figS2.25.f = fig2[which(fig2$plant=="Watermelon" & fig2$threshold==0.25),-1]
figS2.25.f = figS2.25.f[-which(figS2.25.f$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 2), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 4), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 6), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.25.f$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, figS2.25.f$decl_plus_rare-(1.96*figS2.25.f$decl_plus_rare_sd)), 
                    ymax=figS2.25.f$decl_plus_rare+(1.96*figS2.25.f$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.f$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, figS2.25.f$decl-(1.96*figS2.25.f$decl_sd)), 
                    ymax=figS2.25.f$decl+(1.96*figS2.25.f$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.25.f$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, figS2.25.f$rare-(1.96*figS2.25.f$rare_sd)), 
                    ymax=figS2.25.f$rare+(1.96*figS2.25.f$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.25.f$nsites + c(rep(-0.3, 4), 0), y=figS2.25.f$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.f$nsites + c(rep( 0.3, 4), 0), y=figS2.25.f$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.25.f$nsites + c(rep( 0.0, 4), 0), y=figS2.25.f$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.25.f$nsites + c(rep(-0.3, 4), 0), y=figS2.25.f$decl_plus_rare), 
             color="black", fill="gray20", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.f$nsites + c(rep( 0.3, 4), 0), y=figS2.25.f$decl), 
             color="black", fill="goldenrod3", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.25.f$nsites + c(rep( 0.0, 4), 0), y=figS2.25.f$rare), 
             color="black", fill="skyblue3", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,7.5))
#####

# Fig S5 75% 
#####
# panel A
figS2.75.a = fig2[which(fig2$plant=="Polemonium" & fig2$threshold==0.75),-1]
figS2.75.a = figS2.75.a[-which(figS2.75.a$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0 ), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5 ), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 25), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.75.a$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, figS2.75.a$decl_plus_rare-(1.96*figS2.75.a$decl_plus_rare_sd)), 
                    ymax=figS2.75.a$decl_plus_rare+(1.96*figS2.75.a$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.a$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, figS2.75.a$decl-(1.96*figS2.75.a$decl_sd)), 
                    ymax=figS2.75.a$decl+(1.96*figS2.75.a$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.a$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, figS2.75.a$rare-(1.96*figS2.75.a$rare_sd)), 
                    ymax=figS2.75.a$rare+(1.96*figS2.75.a$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.75.a$nsites + c(rep(-0.3, 5), 0), y=figS2.75.a$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.a$nsites + c(rep( 0.3, 5), 0), y=figS2.75.a$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.a$nsites + c(rep( 0.0, 5), 0), y=figS2.75.a$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.75.a$nsites + c(rep(-0.3, 5), 0), y=figS2.75.a$decl_plus_rare), 
             color="black",fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.a$nsites + c(rep( 0.3, 5), 0), y=figS2.75.a$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.a$nsites + c(rep( 0.0, 5), 0), y=figS2.75.a$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,25))


# panel B
figS2.75.b = fig2[which(fig2$plant=="Phacelia" & fig2$threshold==0.75),-1]
figS2.75.b = figS2.75.b[-which(figS2.75.b$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0 ), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5 ), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 25), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.75.b$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, figS2.75.b$decl_plus_rare-(1.96*figS2.75.b$decl_plus_rare_sd)), 
                    ymax=figS2.75.b$decl_plus_rare+(1.96*figS2.75.b$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.b$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, figS2.75.b$decl-(1.96*figS2.75.b$decl_sd)), 
                    ymax=figS2.75.b$decl+(1.96*figS2.75.b$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.b$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, figS2.75.b$rare-(1.96*figS2.75.b$rare_sd)), 
                    ymax=figS2.75.b$rare+(1.96*figS2.75.b$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.75.b$nsites + c(rep(-0.3, 5), 0), y=figS2.75.b$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.b$nsites + c(rep( 0.3, 5), 0), y=figS2.75.b$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.b$nsites + c(rep( 0.0, 5), 0), y=figS2.75.b$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.75.b$nsites + c(rep(-0.3, 5), 0), y=figS2.75.b$decl_plus_rare), 
             color="black",fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.b$nsites + c(rep( 0.3, 5), 0), y=figS2.75.b$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.b$nsites + c(rep( 0.0, 5), 0), y=figS2.75.b$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,25))

# panel C
figS2.75.c = fig2[which(fig2$plant=="Monarda" & fig2$threshold==0.75),-1]
figS2.75.c = figS2.75.c[-which(figS2.75.c$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0 ), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5 ), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 25), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.75.c$nsites + c(rep(-0.3, 5), 0),
                    ymin=pmax(0, figS2.75.c$decl_plus_rare-(1.96*figS2.75.c$decl_plus_rare_sd)), 
                    ymax=figS2.75.c$decl_plus_rare+(1.96*figS2.75.c$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.c$nsites + c(rep( 0.3, 5), 0), 
                    ymin=pmax(0, figS2.75.c$decl-(1.96*figS2.75.c$decl_sd)), 
                    ymax=figS2.75.c$decl+(1.96*figS2.75.c$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.c$nsites + c(rep( 0.0, 5), 0), 
                    ymin=pmax(0, figS2.75.c$rare-(1.96*figS2.75.c$rare_sd)), 
                    ymax=figS2.75.c$rare+(1.96*figS2.75.c$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.75.c$nsites + c(rep(-0.3, 5), 0), y=figS2.75.c$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.c$nsites + c(rep( 0.3, 5), 0), y=figS2.75.c$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.c$nsites + c(rep( 0.0, 5), 0), y=figS2.75.c$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.75.c$nsites + c(rep(-0.3, 5), 0), y=figS2.75.c$decl_plus_rare), 
             color="black",fill="gray30", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.c$nsites + c(rep( 0.3, 5), 0), y=figS2.75.c$decl), 
             color="black", fill="goldenrod", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.c$nsites + c(rep( 0.0, 5), 0), y=figS2.75.c$rare), 
             color="black", fill="skyblue",  size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16,24)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,25))


# panel D
figS2.75.d = fig2[which(fig2$plant=="Blueberry" & fig2$threshold==0.75),-1]
figS2.75.d = figS2.75.d[-which(figS2.75.d$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0 ), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5 ), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 25), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.75.d$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, figS2.75.d$decl_plus_rare-(1.96*figS2.75.d$decl_plus_rare_sd)), 
                    ymax=figS2.75.d$decl_plus_rare+(1.96*figS2.75.d$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.d$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, figS2.75.d$decl-(1.96*figS2.75.d$decl_sd)), 
                    ymax=figS2.75.d$decl+(1.96*figS2.75.d$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.d$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, figS2.75.d$rare-(1.96*figS2.75.d$rare_sd)), 
                    ymax=figS2.75.d$rare+(1.96*figS2.75.d$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.75.d$nsites + c(rep(-0.3, 4), 0), y=figS2.75.d$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.d$nsites + c(rep( 0.3, 4), 0), y=figS2.75.d$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.d$nsites + c(rep( 0.0, 4), 0), y=figS2.75.d$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.75.d$nsites + c(rep(-0.3, 4), 0), y=figS2.75.d$decl_plus_rare), 
             color="black", fill="gray20", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.d$nsites + c(rep( 0.3, 4), 0), y=figS2.75.d$decl), 
             color="black", fill="goldenrod3", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.d$nsites + c(rep( 0.0, 4), 0), y=figS2.75.d$rare), 
             color="black", fill="skyblue3", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,25))


# panel E
figS2.75.e = fig2[which(fig2$plant=="Cranberry" & fig2$threshold==0.75),-1]
figS2.75.e = figS2.75.e[-which(figS2.75.e$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0 ), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5 ), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 25), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.75.e$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, figS2.75.e$decl_plus_rare-(1.96*figS2.75.e$decl_plus_rare_sd)), 
                    ymax=figS2.75.e$decl_plus_rare+(1.96*figS2.75.e$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.e$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, figS2.75.e$decl-(1.96*figS2.75.e$decl_sd)), 
                    ymax=figS2.75.e$decl+(1.96*figS2.75.e$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.e$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, figS2.75.e$rare-(1.96*figS2.75.e$rare_sd)), 
                    ymax=figS2.75.e$rare+(1.96*figS2.75.e$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.75.e$nsites + c(rep(-0.3, 4), 0), y=figS2.75.e$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.e$nsites + c(rep( 0.3, 4), 0), y=figS2.75.e$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.e$nsites + c(rep( 0.0, 4), 0), y=figS2.75.e$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.75.e$nsites + c(rep(-0.3, 4), 0), y=figS2.75.e$decl_plus_rare), 
             color="black", fill="gray20", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.e$nsites + c(rep( 0.3, 4), 0), y=figS2.75.e$decl), 
             color="black", fill="goldenrod3", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.e$nsites + c(rep( 0.0, 4), 0), y=figS2.75.e$rare), 
             color="black", fill="skyblue3", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,25))

# panel F
figS2.75.f = fig2[which(fig2$plant=="Watermelon" & fig2$threshold==0.75),-1]
figS2.75.f = figS2.75.f[-which(figS2.75.f$nsites==2),]

ggplot() +
  
  geom_hline(aes(yintercept= 0 ), lty=1, lwd=0.6, col="black", alpha=0.7) +
  geom_hline(aes(yintercept= 5 ), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 10), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 15), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 20), lty=2, lwd=0.3, col="black", alpha=0.3) +
  geom_hline(aes(yintercept= 25), lty=2, lwd=0.3, col="black", alpha=0.3) +
  
  geom_errorbar(aes(x=figS2.75.f$nsites + c(rep(-0.3, 4), 0),
                    ymin=pmax(0, figS2.75.f$decl_plus_rare-(1.96*figS2.75.f$decl_plus_rare_sd)), 
                    ymax=figS2.75.f$decl_plus_rare+(1.96*figS2.75.f$decl_plus_rare_sd)), 
                color="gray20", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.f$nsites + c(rep( 0.3, 4), 0), 
                    ymin=pmax(0, figS2.75.f$decl-(1.96*figS2.75.f$decl_sd)), 
                    ymax=figS2.75.f$decl+(1.96*figS2.75.f$decl_sd)), 
                color="goldenrod3", width=0.5, lwd=0.8, alpha=0.7)+
  geom_errorbar(aes(x=figS2.75.f$nsites + c(rep( 0.0, 4), 0), 
                    ymin=pmax(0, figS2.75.f$rare-(1.96*figS2.75.f$rare_sd)), 
                    ymax=figS2.75.f$rare+(1.96*figS2.75.f$rare_sd)), 
                color="skyblue3", width=0.5, lwd=0.8, alpha=0.7)+
  
  geom_line(aes(x=figS2.75.f$nsites + c(rep(-0.3, 4), 0), y=figS2.75.f$decl_plus_rare), 
            color="gray20", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.f$nsites + c(rep( 0.3, 4), 0), y=figS2.75.f$decl), 
            color="goldenrod3", alpha=0.7, lwd=1.5) +
  geom_line(aes(x=figS2.75.f$nsites + c(rep( 0.0, 4), 0), y=figS2.75.f$rare), 
            color="skyblue3", alpha=0.7, lwd=1.5) +
  
  geom_point(aes(x=figS2.75.f$nsites + c(rep(-0.3, 4), 0), y=figS2.75.f$decl_plus_rare), 
             color="black", fill="gray20", size=4, pch=21, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.f$nsites + c(rep( 0.3, 4), 0), y=figS2.75.f$decl), 
             color="black", fill="goldenrod3", size=4, pch=22, stroke=1, alpha=0.5) +
  geom_point(aes(x=figS2.75.f$nsites + c(rep( 0.0, 4), 0), y=figS2.75.f$rare), 
             color="black", fill="skyblue3", size=4, pch=23, stroke=1, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Number of Sites", breaks = c(1,4,8,12,16)) +
  scale_y_continuous("Minimum Number\\n of Bee Species Needed", limits=c(0,25))
#####


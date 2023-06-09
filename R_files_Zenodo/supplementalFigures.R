# # If you run the code you may encounter errors, please reach out to
# # the corresponding author, Lea Richardson lkr626@gmail.com for additional information.
# # Make Figure S1 ####
# # Make a script for visualizing the nmds plot associated with ord.Rdata
# # ord.Rdata is a saved ordination made in ff.R that used the whole season of 
# # ff data visualized in 2 dimensions, the axes from this analysis are 
# # the nmds axes in the ff.csv plot, the info for ord.Rdata is shown in comments 
# library(tidyverse)
# library(vegan)
# 
# ff <- read.csv('ff.csv', stringsAsFactors = F)
# set.seed(1)
# # add a columns for proportion of richness that is invasive vs native #
# # make a df to get native and invasive info for each species
# # spList <- read.csv('cofloweringPlants/spForMichele.csv', stringsAsFactors = F)
# # names(spList)[2] <- 'CoFlowSp'
# # ffLong <- read.csv("cofloweringPlants/fnc2011jan26.csv", stringsAsFactors = F)
# # spList <- merge(spList, unique(ffLong[, c('CoFlowSp', 'coflspid')]))
# # spList$sp <- paste('infl', spList$coflspid, sep = '.')
# # write.csv(spList[, c(4,3,1)], 'cofloweringPlants/natInv.csv', row.names = F)
# # get top10 most prevalent
# # g <- vector("numeric", 41L)
# # x <- names(ff[,c(13:53)])
# # for (i in 1:41){
# #   g[i] <- sum(ifelse(ff[, x[i]] == 0, 0, 1))
# # }
# # 
# # g <- as.data.frame(g)
# # g$sp <- x
# # 
# # g[order(g$g), ]
# top10 <- c('infl.ancy', "infl.copa", 'infl.dapu', 'infl.syal',
#            'infl.amca', 'infl.melu', 'infl.meal', 'infl.meof',
#            "infl.mesa", 'infl.ecan')
# 
# remRow <- rowSums(ff[, top10])
# 
# 
# 
# ord <- metaMDS(ff[!remRow %in% 0, top10], trymax = 10000, k=2 )
# 
# ff <- ff %>%
#   mutate(fncDay = as.factor(ff$fncDay)) %>%
#   mutate(infl.ecan = infl.ecan) %>%
#   mutate(inflCt = rowSums(ff[,c(13:53)])) %>%
#   mutate(logAbun = log(inflCt + 1))
# 
# spList <- read.csv('cofloweringPlants/natInv.csv', stringsAsFactors = F)
# ff$invRich <- specnumber(ff[, spList[spList$native %in% 0, 'sp']])
# ff$natRich <- specnumber(ff[, spList[spList$native %in% 1, 'sp']])
# # read in df with native (y/n) column
# # pdf(file = '/Users/learichardson/Desktop/ffOrd.pdf', width = 8, height = 8)
# # par(mar=c(5.1,4.1,4.1,2.1))
# # Begin ord visualisation 
# ffSmall <- ff[!(ff$natRich + ff$invRich) %in% 0 %in% NA, ]
# fig <- ordiplot(ord, type = 'none', bty = 'n', xlim = c(-2,3), ylim = c(-1,2.5))
# text(ord, display = "spec", cex=0.7, col="black")
# fit <- envfit(ord, ff[!remRow %in% 0 , c(2:8, 55:57 )], perm = 1000, na.rm = T)
# fit
# plot(fit, p.max = 0.05, cex=1, labels= list(factors = paste(' ')), 
#      arrow.mul = 4.15, col = 'black', lwd = 4,)
# ffSmall$propNat <- ffSmall$natRich/(ffSmall$natRich+ffSmall$invRich)
# ffSmall[is.nan(ffSmall$propNat), 'propNat'] <- 0
# ffSmall$grayCol <- gray(1-ffSmall$propNat)
# points(fig, "sites", pch=21, bg=ffSmall$grayCol)
# 
# # dev.off()
# 
# # Make Figure S2 ####
# source('pollinators/bb.R')
# 
# 
# #pdf('c:/users/thayes/dropbox/fnc2017/plots/bbBarPlotsSupplemental.pdf', width = 4, height = 8)
# par(mfrow = c(3,1))
# 
# # BB by Site 
# 
# # Add abundance and richness for site combination
# bbSite <- aggregate(bb.2$sp, by = list(bb.2$site, bb.2$sp), FUN = length)
# names(bbSite) <- c('site', 'sp', 'abun')
# 
# # Make a bbSiteTable that has abundance and richness by site (not just abun by species by site)
# bbSiteRich <- aggregate(bbSite$sp, by = list(bbSite$site), FUN = length)
# names(bbSiteRich) <- c('site', 'rich')
# bbSiteAbun <- aggregate(bbSite$abun, by = list(bbSite$site), FUN = sum)
# names(bbSiteAbun) <- c('site', 'abun')
# bbSiteTable <- merge(bbSiteRich, bbSiteAbun, by = 'site')
# cols = c('deepskyblue1', 'darkorange1')
# barplot(height = t(bbSiteTable[, c('rich', 'abun')]), beside = T,
#         names.arg = bbSiteTable$site, col = cols, border = cols,
#         ylim = c(0,32), legend.text=c('richness','abundance'),
#         args.legend=list(col=cols,border=cols,bty='n'),
#         ylab = '# species, # individuals', xlab = 'Site')
# 
# # BB by Day 
# 
# # Add abundance and richness for obsDay combination
# bbDay <- aggregate(bb.2$sp, by = list(bb.2$obsDay, bb.2$sp), FUN = length)
# names(bbDay) <- c('obsDay', 'sp', 'abun')
# 
# # Make a bbDayTable that has abundance and richness by obsDay
# bbDayRich <- aggregate(bbDay$sp, by = list(bbDay$obsDay), FUN = length)
# names(bbDayRich) <- c('obsDay', 'rich')
# bbDayAbun <- aggregate(bbDay$abun, by = list(bbDay$obsDay), FUN = sum)
# names(bbDayAbun) <- c('obsDay', 'abun')
# bbDayTable <- merge(bbDayRich, bbDayAbun, by = 'obsDay')
# cols = c('deepskyblue1', 'darkorange1')
# par(mar = c(5,5,2,2))
# barplot(height = t(bbDayTable[, c('rich', 'abun')]), beside = T,
#         names.arg = bbDayTable$obsDay, col = cols, border = cols,
#         ylim = c(0,70), legend.text=c('richness','abundance'),
#         args.legend=list(col=cols,border=cols,bty='n'),
#         ylab = '# species, # individuals', xlab = 'Observation day')
# 
# # BB by Site/Day 
# # Add abundance and richness for site combination
# bbBoth <- aggregate(bb.2$sp, by = list(bb.2$site, bb.2$sp, bb.2$obsDay), FUN = length)
# names(bbBoth) <- c('site', 'sp', 'obsDay', 'abun')
# bbBoth$siteDay <- paste(bbBoth$site, bbBoth$obsDay, sep = '')
# 
# # Make a bbSiteTable that has abundance and richness by both site and obsDay
# bbBothRich <- aggregate(bbBoth$sp, by = list(bbBoth$siteDay), FUN = length)
# names(bbBothRich) <- c('siteDay', 'rich')
# bbBothAbun <- aggregate(bbBoth$abun, by = list(bbBoth$siteDay), FUN = sum)
# names(bbBothAbun) <- c('siteDay', 'abun')
# bbBothTable <- merge(bbBothRich, bbBothAbun, by = 'siteDay')
# cols = c('deepskyblue1', 'darkorange1')
# barplot(height = t(bbBothTable[, c('rich', 'abun')]), beside = T,
#         names.arg = bbBothTable$siteDay, col = cols, border = cols,
#         ylim = c(0,13), legend.text=c('richness','abundance'),
#         args.legend=list(col=cols,border=cols,bty='n'),
#         ylab = '# species, # individuals', xlab = 'Site-day combination',
#         las = 2)
# 
# #dev.off()
# par(mfrow = c(1,1))
# 
# 
# # I need ff to get characteristics by site
# ff <- read.csv('cofloweringPlants/ff.csv', stringsAsFactors = F)
# siteNn2 <- aggregate(ff$nn2, by = list(ff$site), FUN = mean)
# siteRich <- aggregate(ff$rich, by = list(ff$site), FUN = mean)
# siteEnv <- merge(siteNn2, siteRich, by = 'Group.1')
# names(siteEnv) <- c('site', 'meanNn2', 'meanRich')
# 
# siteDapu <- aggregate(ff$infl.dapu, by = list(ff$site), FUN = mean)
# names(siteDapu) <- c('site', 'dapu')
# siteAncy <- aggregate(ff$infl.ancy, by = list(ff$site), FUN = mean)
# names(siteAncy) <- c('site', 'ancy')
# siteMelu <- aggregate(ff$infl.melu, by = list(ff$site), FUN = mean)
# names(siteMelu) <- c('site', 'melu')
# siteSyal <- aggregate(ff$infl.syal, by = list(ff$site), FUN = mean)
# names(siteSyal) <- c('site', 'syal')
# siteMesa <- aggregate(ff$infl.mesa, by = list(ff$site), FUN = mean)
# names(siteMesa) <- c('site', 'mesa')
# siteAmca <- aggregate(ff$infl.amca, by = list(ff$site), FUN = mean)
# names(siteAmca) <- c('site', 'amca')
# siteCopa <- aggregate(ff$infl.copa, by = list(ff$site), FUN = mean)
# names(siteCopa) <- c('site', 'copa')
# 
# m1 <- merge(siteDapu, siteAncy, by = 'site')
# m2 <- merge(m1, siteMelu, by = 'site')
# m3 <- merge(m2, siteSyal, by = 'site')
# m4 <- merge(m3, siteMesa, by = 'site')
# m5 <- merge(m4, siteAmca, by = 'site')
# m6 <- merge(m5, siteCopa, by = 'site')
# 
# siteEnv <- merge(siteEnv, m6, by = 'site')
# 
# 
# # Start ordination by site
# site1<- sort(rep(unique(bbSite$site), length(unique(bbSite$sp))))
# sp1 <- rep(sort(unique(bbSite$sp)), length(unique(bbSite$site)))
# 
# preBbSiteOrd <- as.data.frame(cbind(site1, sp1))
# names(preBbSiteOrd) <- c("site", "sp")
# 
# preBbSiteOrd2 <- merge(preBbSiteOrd, bbSite[ , c("site", "sp", 'abun')], by = c("site", "sp"), all.x= T)
# 
# #make empty dataframe formatted for NMDS
# 
# bbSiteOrd <- as.data.frame(matrix(ncol = length(unique(preBbSiteOrd2$sp)), nrow = 0))
# names(bbSiteOrd) <- sort(unique(preBbSiteOrd2$sp))
# spList <- sort(unique(preBbSiteOrd2$sp))
# 
# for (activesite in sort(unique(preBbSiteOrd2$site))) {
#   bbSiteOrd[nrow(bbSiteOrd)+1, ] <- NA
#   
#   for (activeSp in 1:length(spList)) {
#     bbSiteOrd[nrow(bbSiteOrd), activeSp] <- preBbSiteOrd2[preBbSiteOrd2$site %in% activesite & preBbSiteOrd2$sp %in% spList[activeSp], ]$abun
#   }
# }
# 
# bbSiteOrd[is.na(bbSiteOrd)] <- 0
# everything <- cbind(siteEnv, bbSiteOrd)
# 
# #try NMDS, k = 2
# 
# ord <- metaMDS(everything[, c(11:19)], zerodist="add", k = 2)
# 
# site.sc <- scores(ord, display = 'sites')
# 
# # pdf('c:/users/thayes/dropbox/fnc2017/plots/bbSiteNMDS.pdf', width = 7, height = 5)
# 
# plot(ord, type ="n", xlim = c(-1.5, 1.5))
# points(ord, display = "sites", cex = 1, pch=21, col="white",
#        bg="white")
# text(ord, display = "species", cex=0.8, col="blue")
# ordihull(ord, unique(everything$site), label = T)
# 
# 
# fit <- envfit(ord, everything[, c(1:10)], perm = 1000, na.rm = T)
# fit
# plot(fit, p.max = 0.05, cex=1, col = 'black')
# 
# # dev.off()
# 
# 
# stress <- as.data.frame(matrix(nrow = 0, ncol = 3))
# 
# names(stress) <- c('type', 'k2', 'k3')
# stress[1,'type'] <- 'bySite'
# stress[1,'k2'] <- ord$stress
# 
# # Make Figure S3 ####
# library(vegan) 
# library(MASS)
# 
# # Import and Clean Data
# pol <-read.csv('ppHet.csv')
# abunRich <-read.csv('ppCon.csv')
# 
# levels <- levels(abunRich$sp)
# levels[length(levels) + 1] <- "None"
# abunRich$sp <- factor(abunRich$sp, levels = levels)
# abunRich$sp <- as.character(abunRich$sp)
# abunRich[is.na(abunRich$sp), "sp"] <- "unidentified"
# 
# 
# # Remove non-pollinators and one's with very low abundance
# abunRich <-abunRich[abunRich['totAbun']!=0,]
# abunRich <-abunRich[abunRich['totAbun']!=1,] 
# abunRich <-abunRich[abunRich['sp']!='Coleoptera sp. 2',]
# abunRich <-abunRich[abunRich['sp']!='Coelioxys rufitarsis',]
# abunRich <-abunRich[abunRich['sp']!='unidentified',]
# abunRich <-abunRich[abunRich['sp']!='Megachile sp. 1',]
# abunRich <-abunRich[abunRich['sp']!='Megachile sp. 2',]
# 
# 
# 
# # add sp to pol data 
# pol <- merge(pol, abunRich,
#              by.y=c("specimenID"),
#              all.x=TRUE)
# 
# # remove NAs 
# pol <- na.omit(pol)
# str(pol)
# 
# # long to wide 
# pol <- reshape(pol, idvar = c('specimenID' , 'id' , 'sp' , 'notes'  
#                               , 'totAbun' , 'codeRichnessWithSpro',  'numSpro', 'hetAbun') , direction = 'wide',
#                timevar = 'polCode')
# 
# 
# # rename columns 
# names(pol) <- gsub("polCount.", "", names(pol), fixed = TRUE)
# 
# 
# # replace NAs w/ 0 
# str(pol)
# i <- sapply(pol, is.factor)
# pol[i] <- lapply(pol[i], as.character)
# str(pol)
# pol[is.na(pol)] <- 0
# pol[8:16] <- lapply(pol[8:16], as.integer)
# pol <- na.omit(pol)
# str(pol)
# 
# 
# # analysis 
# 
# xx <- pol[,1:7] # all the other variables
# xx [xx %in% NA]
# xx$sp <- factor(xx$sp)
# 
# pp <- pol[,9:16]  # just pollen types 
# # pp$sp <-pol[,3]
# pp[pp %in% NA]
# str(pp)
# str(xx$sp)
# pp.mds <- metaMDS(pp, zerodist="add")
# pp.mds
# summary(pp.mds)
# plot(pp.mds, type= "t")
# # 
# pp.mds$points
# vegdist(pp)
# 
# 
# # propSPRO column in xx 
# xx$prop <- (xx$numSpro / xx$totAbun)
# xx$col[xx$prop==0]="grey75"
# xx$col[xx$prop>0 & xx$prop<=0.2]="grey65"
# xx$col[xx$prop>0.2 & xx$prop<=0.4]="grey52"
# xx$col[xx$prop>0.4 & xx$prop<=0.6]="grey39"
# xx$col[xx$prop>0.6 & xx$prop<=0.8]="grey26"
# xx$col[xx$prop>0.8 & xx$prop<1]="grey20"
# xx$col[xx$prop==1]="grey0"
# 
# # Plot S3 NMDS
# fig <- ordiplot(pp.mds, type = 'none', bty = 'n')
# ordihull(fig, xx$sp, lwd = 3,
#          col = brewer.pal(name ='Dark2',n =8))
# points(fig, "sites", pch=as.integer(xx$sp), col=adjustcolor(xx$col, alpha.f = 0.8), 
#        cex=1.5 , lwd=2.5) 
# 
# names(xx)[3] <- 'richness of \\npollen categories'
# fit <- envfit(pp.mds, xx[,c(1:3)], perm = 1000)
# fit
# plot(fit, p.max = 0.05, col = 'black', cex=1)
# legend("topleft", title="Pollinator Taxa", pch = c(1:8),
#        c("Agopostemon virescens", 'Augochlorella striata',
#          "Ceratina calcarata/dupla","Dialictus spp.", 'Halictus ligatus',
#          'Halictus parallelus', 'Melissodes spp.',
#          'Pseudopanurgus spp.'), 
#        col = brewer.pal(name ='Dark2',n =8), cex=0.8, bty = 'n')
# 
# legend("bottomleft", title="Proportion of \\nconspecific pollen", pch = c(15),
#        c("0", '>0-0.2', '>0.2-0.4', '>0.4-0.6', '>0.6-0.8', '>0.8-0.99', '1'), 
#        col = c("grey75","grey65","grey52","grey39","grey26", "grey15","grey0"), 
#        cex=0.8, bty = 'n')
# title('Characterization of pollen community on pollinators 
#       \\nshaded by proportion of conspecific pollen')
# 

#### load packages and data ####
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

ff <- read.csv('ff.csv', stringsAsFactors = F)
bb <- read.csv('bb.csv', stringsAsFactors = F)
ppCon <-read.csv('ppCon.csv')
ppHet <-read.csv('ppHet.csv', stringsAsFactors = F)
ssHet <- read.csv('ssHet.csv', stringsAsFactors = F)
ssCon <- read.csv('ssCon.csv', stringsAsFactors = F)
rr <- read.csv('rr.csv', stringsAsFactors = F)

bff <- merge(bb, ff, by = c('plantObsID', 'site'))
bff$fncDay<- as.factor(bff$fncDay)
sproNOEA <- c('infl.acmi',
'infl.ciar',
'infl.cifl',
'infl.copa',
'infl.erst',
'infl.phpi',
'infl.somi')
polColSpro <- bff[bff$polCollected %in% 1, c('specimenID', sproNOEA)]
sum(rowSums(polColSpro[,sproNOEA]) >0)  # 24 insects were collected in areas with noea spro 
look <- polColSpro[rowSums(polColSpro[,sproNOEA]) >0, ]

expandlook <- bff[bff$specimenID %in% look$specimenID, ] 

abunRich <- ppCon
levels <- levels(abunRich$sp)
levels[length(levels) + 1] <- "None"
abunRich$sp <- factor(abunRich$sp, levels = levels)
abunRich$sp <- as.character(abunRich$sp)
abunRich[is.na(abunRich$sp), "sp"] <- "unidentified"
abunRich <- na.omit(abunRich)
pp1 <-abunRich[abunRich['totAbun']!=0,]
pp2 <-pp1[pp1['totAbun']!=1,]
pp3 <-pp2[pp2['sp']!='Coleoptera sp. 2',]
pp4 <-pp3[pp3['sp']!="Coelioxys rufitarsis",]
pp5 <-pp4[pp4['sp']!="Megachile sp. 1",]
pp6 <-pp5[pp5['sp']!="Megachile sp. 2",]
pp7 <-pp6[pp6['sp']!="unidentified",]
pp7$sp<-factor(pp7$sp)

expandPP <- merge(expandlook, pp7, by = 'specimenID') # 19 pols could have noEA spro



j <- merge(expandPP, ppHet, by = 'specimenID')

#### determine pollinators and styles at risk of mis-identified pollen ####
#### this includes all specimens from neighborhoods that have NOEA flowering plants ####
sum(j$polCode %in% 'NOEA') # 5 have id'd NOEA pollen, which leaves 14 remaining specimens with spro in the 
# neighborhood and no mention of NOEA, which suggests that 14/100 specimens in the pp dataset

specToRemove <- j[j$polCode %in% 'NOEA', 'specimenID'] 

specsToRemoveFinal <- expandPP[!expandPP$specimenID %in% specToRemove, ]

totalProbSpecs <- j[!j$specimenID %in% specToRemove, 'specimenID']
  
ppNoPolProbs <- pp7[!pp7$specimenID %in% totalProbSpecs, ]



# Creating ssffCon
str(ssCon) # 337 obs. of  6 variables
str(ff) # 224 obs. of  53 variables

ssffCon <- merge(ssCon, ff, by = 'plantObsID', all = T)
str(ssffCon) # 348 obs. of  58 variables
ssffCon$plantObsID <- as.factor(ssffCon$plantObsID)
ssffCon$styleID <- as.factor(ssffCon$styleID)
ssffCon$fncDay <- as.factor(ssffCon$fncDay)
ssffCon$site <- as.factor(ssffCon$site)
ssffCon$tag <- as.factor(ssffCon$tag)

str(ssffCon)
head(ssffCon)
names(ssffCon) # preStyle = 1 is pre-obs style
ssffCon0 <- ssffCon[ssffCon$preStyle == 0,]

specIDsToRemove <- j[j$polCode %in% 'NOEA', 'specimenID']

specsWithNoeaSpro<- expandPP[!expandPP$specimenID %in% specIDsToRemove,  ]


ssffHet <- merge(ssHet, ff, by = 'plantObsID', all = T)
str(ssffHet) # 402 obs. of  56 variables
ssffHet$plantObsID <- as.factor(ssffHet$plantObsID)
ssffHet$styleID <- as.factor(ssffHet$styleID)
ssffHet$fncDay <- as.factor(ssffHet$fncDay)
ssffHet$site <- as.factor(ssffHet$site)
ssffHet$tag <- as.factor(ssffHet$tag)



sum(na.omit(rowSums(ssffCon[,sproNOEA]) >0))  # 73 styles were collected in areas with non ech spro 
sum(ssffHet$polCode %in% 'NOEA') # 14 styles collected in areas with NOEA code 
stylesWNoea<- ssffHet[ssffHet$polCode %in% 'NOEA', 'styleID']
stylesToRem1 <- ssffCon[(rowSums(ssffCon[,sproNOEA]) >0), 'styleID']
stylesToRemFinal <- stylesToRem1[!stylesToRem1 %in% stylesWNoea]

sum(na.omit(rowSums(ssffCon[,sproNOEA]) >0)) - 
  sum(ssffHet$polCode %in% 'NOEA') # 59 styles that likely could have non ech spro 59/337= 18%

#### Fig 4a, Visualize proportion of pollen type by pollinator taxon ####

pol <- ppHet
abunRich <- ppCon
pol <- pol[pol$specimenID %in% pp7$specimenID, ]
abunRich <- abunRich[abunRich$specimenID %in% pp7$specimenID, ]
# rename unidentified rows in sp
levels <- levels(abunRich$sp)
levels[length(levels) + 1] <- "None"
abunRich$sp <- factor(abunRich$sp, levels = levels)
abunRich$sp <- as.character(abunRich$sp)
abunRich[is.na(abunRich$sp), "sp"] <- "unidentified"
# abun

# total abundance by sp 
abun <- aggregate(abunRich$totAbun, by=list(sp=abunRich$sp), FUN=sum)

colnames(abun)[2] <- "totAbun"

# create new df 
merge <- merge(pol, abunRich, by = c("specimenID"))

# combine pollen types by sp 
merge <- aggregate(polCount~sp+polCode, data=merge, FUN=sum) 

#create column with pollen count totals for each sp 
merge <- merge(merge, abun, by = c("sp"))

# create proportions of pollen types by sp 
merge<- transform(merge, prop =  polCount / totAbun )

# remove rows w/ zero values  
merge <- merge[!(merge$polCount==0 & merge$totAbun==0),]
merge[merge==""]<-NA
merge<-merge[complete.cases(merge),]

# remove non-pollinators & Megachile sp. 1 (bc it has only 1 pollen grain)
merge <-merge[merge['sp']!='Coleoptera sp. 2',]
merge <-merge[merge['sp']!='Coelioxys rufitarsis',]
merge <-merge[merge['sp']!='Megachile sp. 1',]

# stacked barchart, code from MKG 


table(merge$sp)

sp_order <- c('Halictus ligatus', 'Pseudopanurgus sp.', 'Megachile sp. 2',
              'Dialictus sp.', 'Agapostemon virescens', 
              'Ceratina calcarata/dupla',  'Melissodes sp.', 'Halictus parallelus',
              'Augochlorella striata', 'unidentified')

pol_order <- c('LISU', 'MOFI', 'NOEA', 
               'SMOS', 'SMPS', 'SMRO',
               'UNID', 'SPRO')

#value <- c('.0996','0.995', '0.995','0.913','.901','.893',  '.853','.842','.770', .140')
# for writing on the bars? 

merge <- transform(merge, sp = factor(sp, levels = sp_order))
merge <- transform(merge, polCode = factor(polCode, levels = pol_order))

merge <-merge[merge['sp']!='unidentified',]
merge <-merge[merge['sp']!='Megachile sp. 2',]

c <- ggplot(merge) +
  geom_bar(data=merge,
           aes(y = prop, x = sp, fill = polCode),
           stat="identity",
           position='stack') +
  scale_fill_brewer(palette = 'Paired', direction = -1)+
  scale_x_discrete(limits=c( "Melissodes sp.",
                             "Halictus parallelus",
                             "Agapostemon virescens",
                             "Halictus ligatus",
                             "Ceratina calcarata/dupla",
                             "Pseudopanurgus sp.",
                             'Augochlorella striata',
                             "Dialictus sp."),
                   labels = c( "8. Melissodes \\nspp.",
                               "7. Halictus \\nparallelus",
                               "6. Agapostemon \\nvirescens",
                               "5. Halictus \\nligatus",
                               "4. Ceratina \\ncalcarata/dupla",
                               "3. Pseudopanurgus \\nspp.",
                               '2. Augochlorella \\nstriata',
                               "1. Dialictus \\nspp."))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_flip() +
  ggtitle("Proportion of pollen by category on pollinators") +
  labs(x = "Pollinator taxon", y = "Proportion of pollen", fill = "Pollen \\ncategory")

#### Does pollen abundance differ by pollinator taxon? ####
# model abun
alt.modelA <- glm(totAbun ~ sp, data = pp7, family = quasipoisson)
null.modelA <- glm(totAbun ~ 1, data = pp7, family = quasipoisson)
anova(null.modelA, alt.modelA, test = "F")
summary(alt.modelA) 

# order bee taxa by body size

order <- ordered(pp7$sp, levels = c("Melissodes sp.", "Halictus parallelus", "Agopostemon virescens", 
                                    "Halictus ligatus", "Ceratina calcarata/dupla", "Pseudopanurgus sp.",
                                    'Augochlorella striata', "Dialictus sp."))


# predicted values, abun

pp7$fitA <- predict.glm(alt.modelA, type = 'response', se.fit = T)$fit
pp7$predictA <- predict.glm(alt.modelA, type = 'response', se.fit = T)$se.fit
pp7

fitted <- unique(pp7$fitA)
stderr <- unique(pp7$predictA)
fitplus <- fitted + stderr
fitminus <- fitted - stderr




#### Make Figure 4b ####

a <- ggplot(pp7, aes(sp, fitA)) +
  geom_point() +
  geom_errorbar(aes(ymin = fitA - predictA, 
                    ymax = fitA + predictA,
                    width=0.05)) +
  scale_x_discrete(limits=c("Dialictus sp.",
                            'Augochlorella striata',
                            "Pseudopanurgus sp.",
                            "Ceratina calcarata/dupla",
                            "Halictus ligatus",
                            "Agapostemon virescens",
                            "Halictus parallelus",
                            "Melissodes sp.")) +
  theme(axis.text.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = '', y = "Total pollen grains") +
  ggtitle("Pollen abundance by pollinator taxon (+/- SE)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 200))
# theme(axis.text.x = element_blank(),)
#### Does count of pollen types differ by pollinator taxon ####
alt.modelR <- glm(codeRichnessNoSpro ~ sp, data = pp7, family = quasipoisson)
null.modelR <- glm(codeRichnessNoSpro ~   1, data = pp7, family = quasipoisson)
anova(null.modelR, alt.modelR, test = "F")
summary(alt.modelR) # might be good enough to just report this table? IDK

# predicted values, rich
pp7$fitR <- predict.glm(alt.modelR, type = 'response', se.fit = T)$fit
pp7$predictR <- predict.glm(alt.modelR, type = 'response', se.fit = T)$se.fit
pp7

fittedR <- unique(pp7$fitR)
stderrR <- unique(pp7$predictR)
fitplusR <- fittedR + stderrR
fitminusR <- fittedR - stderrR

#### Make Figure 4c ####
b <- ggplot(pp7, aes(sp, fitR)) +
  geom_point() +
  geom_errorbar(aes(ymin = fitR - predictR, 
                    ymax = fitR + predictR,
                    width=0.05)) +
  scale_x_discrete(limits=c("Dialictus sp.",
                            'Augochlorella striata',
                            "Pseudopanurgus sp.",
                            "Ceratina calcarata/dupla",
                            "Halictus ligatus",
                            "Agapostemon virescens",
                            "Halictus parallelus",
                            "Melissodes sp."),
                   labels = c('1','2','3','4', '5','6', '7', '8')) +
  theme(axis.text.x = element_text(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(x = "Pollinator taxon", y = "Total pollen richness") +
  ggtitle("Pollen richness by pollinator taxon (+/- SE)") +
  theme(plot.title = element_text(hjust = 0.5)) 


#### Full Figure 4 ####
ppPlots <- ggarrange(c, ggarrange(a, b, labels = c("b", "c"),
                                  ncol = 1, nrow = 2,  align = 'v'), labels = c( 'a', ''),
                     ncol = 2, nrow = 1, widths = c(2,1))


ppPlots


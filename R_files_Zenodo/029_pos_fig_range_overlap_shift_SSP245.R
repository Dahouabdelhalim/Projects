#' ----
#' Range overlap SSP245
#' ---- 

# prepare r ---------------------------------------------------------------

options(digits = 3, scipen = 999)

setwd(here())

source('00_packages.R')

source('01_settings.R')

# heatmap present ---------------------------------------------------------

setwd('range_tables')

withself <- read.xlsx('range_overlap_present_R1.xlsx', sheetIndex = 1)

omatrix <- reshape(withself, direction="wide", idvar="sp2", timevar="sp1")

omatrix

withself$sp1 <- sub('_', ' ', stringr::str_to_sentence(withself$sp1))

withself$sp2 <- sub('_', ' ', stringr::str_to_sentence(withself$sp2))

col.range=c(10000,max(withself$area_overlapped))

#b=c(1200000, 440000, 162000 ,59000, 20000)

#limits=col.range,

oheat <- ggplot(data=withself, aes(y = sp1, x=sp2, fill =area_overlapped ) ) +
  geom_tile(alpha=0.75) + scale_fill_viridis( option = "D", direction=-1,  #trans = "log",
                                              limits=col.range) + #breaks=c ) +
  theme_bw()+
  theme(plot.title.position = "plot", axis.title=element_text(size=11),
        axis.text.y = element_text(size=13, face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 0.99, size=13,  face = "italic"))+
  ylab('') +xlab('') +labs(fill = "Area overlap (sqkm)        ") 

oheat

#fign <- paste0('Fig_overlap_present_', version_suffix, '.jpg')
#ggsave(filename = fign, oheat, width =33, height = 29, units = 'cm', dpi=600)

# heatmap future ---------------------------------------------------------

withselff <- read.xlsx('overlap_future_SSP245_R1.xlsx', sheetIndex = 1)
head(withselff)

omatrixf <- reshape(withself, direction="wide", idvar="sp2", timevar="sp1")

withselff$sp1 <- sub('_', ' ', stringr::str_to_sentence(withselff$sp1))

withselff$sp2 <- sub('_', ' ', stringr::str_to_sentence(withselff$sp2))

oheatfut <- ggplot(data=withselff, aes(y = sp1, x=sp2, fill =area_overlapped_future ) ) +
  geom_tile(alpha=0.75) + scale_fill_viridis( option = "D", direction=-1,  #trans = "log",
                                              limits=col.range) + 
  theme_bw()+
  theme(plot.title.position = "plot", axis.title=element_text(size=11),
        axis.text.y = element_text(size=13, face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 0.99, size=13,  face = "italic"))+
  ylab('') +xlab('') +labs(fill = "Future area overlap (sqkm)") 

oheatfut
getwd()

compo <- gridExtra::grid.arrange(oheat, oheatfut, nrow = 2)
compo

fign <- paste0('Fig_overlap_future_SSP245_R1_', version_suffix, '.jpg')
#ggsave(filename = fign, compo, width =26, height = 37, units = 'cm', dpi=600)

# Left join

withdif <- left_join(withself, withselff)

withdif$dif <- withdif$area_overlapped_future - withdif$area_overlapped
hist(withdif$dif)

summary(withdif$dif)
col.range1=c(min(withdif$dif),max(withdif$dif))

# Negative values is area overlap loss in the future
withdif$dif

oheatdif<- ggplot(data=withdif, aes(y = sp1, x=sp2, fill =scale(dif) ) ) +
  geom_tile(alpha=0.75) + scale_fill_viridis( option = "D", direction=-1               ) + 
  theme_bw()+
  theme(plot.title.position = "plot", axis.title=element_text(size=11),
        axis.text.y = element_text(size=13, face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 0.99, size=13,  face = "italic"))+
  ylab('') +xlab('') +labs(fill = "Shift in future area overlap (scaled sqkm)") 

oheatdif

oheatdifr <- ggplot(data=withdif, aes(y = sp1, x=sp2, fill =dif ) ) +
  geom_tile(alpha=0.75) + scale_fill_viridis( option = "D", direction=-1    ,
                                              breaks=c(400000, 0, -500000 ,-1000000, -1500000,-2000000) ) + 
  theme_bw()+
  theme(plot.title.position = "plot", axis.title=element_text(size=11),
        axis.text.y = element_text(size=13, face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 0.99, size=13,  face = "italic"))+
  ylab('') +xlab('') +labs(fill = "Shift in future area overlap (sqkm)") 

oheatdifr

#figr <- paste0('Fig_overlap_supplements_SSP245_R1', version_suffix, '.jpg')
#ggsave(filename = figr, oheatdifr, width =33, height = 29, units = 'cm', dpi=600)

# overlap present-future --------------------------------------------------

head(withdif)

withdif$difpct <- withdif$area_overlapped_future / withdif$area_overlapped 
withdif$difpct

require(tidyverse)
require(janitor)
nrow(withdif)

withdif$difpct

withdif$category <- ifelse(withdif$difpct == 'NaN', "Never overlaps",
                           ifelse(withdif$difpct == 'Inf', "Starts overlapping in the future", 
                                  ifelse(withdif$difpct == 0, "Stops overlaping in the future",
                                         ifelse(withdif$difpct <  1, "Overlaps less in the future",
                                                ifelse(withdif$difpct > 1, "Overlaps more in the future", NA)))))

table(withdif$difpct == 0)

0/0
1/0
0/1

withdif$difpct == 'Inf'
withdif

mypali <- viridisLite::viridis(5)

oheatdifpct <- ggplot(data=withdif, aes(y = sp1, x=sp2, fill =category ) ) +
  geom_tile(alpha=0.75) + scale_fill_manual(values=mypali )+
  theme_bw()+
  theme(plot.title.position = "plot", axis.title=element_text(size=11),
        axis.text.y = element_text(size=13, face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 0.99, size=13,  face = "italic"))+
  ylab('') +xlab('') +labs(fill = "Shift in area overlap") 

oheatdifpct

head(withdif)

fign <- paste0('Fig_overlap_shift_SSP245_R1_', version_suffix, '.jpg')
ggsave(filename = fign, oheatdifpct, width =33, height = 29, units = 'cm', dpi=600)

write.xlsx(withdif, 'overlap_shift_SSP245_R1.xlsx', row.names = FALSE)

# Summary table merge

head(withdif)

trange <- read.xlsx('Table_Range_results_40o_ss6_maxent_15rep.xlsx', sheetIndex = 1)

head(trange)

withdif$species <- withdif$sp1

suptable245 <- left_join(withdif, trange)

head(suptable245)

suptable245$pct_range_overlap_present <- suptable245$area_overlapped / suptable245$Potential.range

hist( suptable245$pct_range_overlap_present )

suptable245$pct_range_overlap_futureSSP245 <- suptable245$area_overlapped_future/ suptable245$Future.range.SSP.245.BCC

hist( suptable245$pct_range_overlap_futureSSP245  )

supexpo245 <- suptable245 %>% arrange(species)

options(digits = 3)

head(supexpo245)

write.xlsx(supexpo245, 'range_overlap_shift_SSP245.xlsx', row.names = FALSE)

#-------------------------------------------------------------------------------------------------------------
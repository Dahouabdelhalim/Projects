# Range overlap Future
# Building loop for individual maps for every scenario, GCM and period

options(digits = 3, scipen = 999)

setwd(projdir)

source('00_packages.R')

source('01_settings.R')

setwd(enmresultsdir)
setwd('Projection')
list.files()

setwd(paste0('BCC-CSM2-MR_ssp245_2081-2100_27kms/','/Ensemble/W_MEAN/MAX_TSS/'))

all <- stack(list.files(pattern='.tif$'))

setwd(aa)

all_aa <- stack(list.files(pattern='.tif$'))

allm <- mask(all, all_aa)

allm[values(allm)==0] <- NA

recebeself <- data.frame()

for(i in names(allm))
  for(j in rev(names(allm)) )
  {
    if (i == j){
      next
    }
    summed <- allm[[i]] + allm[[j]]
    
    values(summed)[values(summed) != 2] = NA
    
    ar <- area(summed, na.rm=TRUE)
    
    meui <- i
    meuj <- j
    areasoma <- sum(values(ar),  na.rm = TRUE)
    
    aux <- data.frame(sp1= meui, sp2=meuj ,area_overlapped_future= areasoma)
    recebeself <- rbind(recebeself, aux)
    
    print(i)
    print(j) 
  }

head(recebeself)

rexpo <- recebeself %>% arrange(desc(area_overlapped_future))

setwd(here())

#export---

setwd('range_tables')

write.xlsx(rexpo, 'overlap_future_SSP245_R1.xlsx', row.names = FALSE)

withself <- read.xlsx('overlap_future_SSP245_R1.xlsx', sheetIndex = 1)
head(withself)

omatrix <- reshape(withself, direction="wide", idvar="sp2", timevar="sp1")

omatrix

head(rexpo)

withself$sp1 <- sub('_', ' ', stringr::str_to_sentence(withself$sp1))

withself$sp2 <- sub('_', ' ', stringr::str_to_sentence(withself$sp2))

#_SSP245
oheatfut <- ggplot(data=withself, aes(y = sp1, x=sp2, fill =area_overlapped_future ) ) +
  geom_tile(alpha=0.75) + scale_fill_viridis( option = "D", direction=-1,  trans = "log") + 
  theme_bw()+
  theme(plot.title.position = "plot", axis.title=element_text(size=11),
        axis.text.y = element_text(size=13, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 0.99, size=13,  face = "italic"))+
  ylab('') +xlab('') +labs(fill = "Future area overlap (sqkm)") 

oheatfut

fign <- paste0('Fig_overlap_future_SSP245_R1_', version_suffix, '.jpg')
ggsave(filename = fign, oheatfut, width =33, height = 29, units = 'cm', dpi=600)
##############################################################################################################

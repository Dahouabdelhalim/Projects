library(tidyverse)
library(nlme)
library(multcomp)

#  irep data analysis

rp <- read.table("resubmission_files/irepvals_allsamples_v_all_filt_5kb_derep_98_genomes.txt", header=T, sep='\\t')
bin.tax <- read.table("resubmission_files/gtdbtk.bac120.summary.tsv", header=T, sep='\\t') %>%
  dplyr::rename(genome = user_genome)

rp.l <- rp %>%
  dplyr::select(-ExtractionID) %>%
  pivot_longer(!c(BeeID,Age_days,Colony), names_to="genome",values_to="iRep") %>%
  drop_na() %>%
  full_join(bin.tax, by="genome") %>%
  mutate(genus = classification) %>%
  separate(col=genus,sep=";",into=c("Domain",
                                    "Phylum",
                                    "Class",
                                    "Order",
                                    "Family",
                                    "Genus",
                                    "Species")) %>%
  mutate(Genus_Sp_Bin = paste(Genus,Species,genome))


# plot schmid and gilliamella

rp.l.filt_sg <- rp.l %>%
  filter(Genus %in% c("g__Gilliamella","g__Schmidhempelia")) %>%
  mutate(Genus = as.factor(gsub("g__","",Genus))) %>%
  mutate(Colony = as.factor(Colony))

rp.l.filt_sg$Genus <- factor(rp.l.filt_sg$Genus, c("Schmidhempelia","Gilliamella"))

## Fig. 5

irep.plot.sg <- ggplot(rp.l.filt_sg, aes(x=Age_days,y=iRep, fill=Colony)) +
  geom_point(size=4, alpha=0.7, shape=21) +
  scale_fill_manual(values=c("Y"="#ffeda0",
                              "B"="#9ecae1",
                              "W"="#f0f0f0")) +
  ylim(c(0.95,2.3)) +
  facet_wrap(~Genus, nrow=2,ncol=1, scales="fixed") +
  geom_hline(yintercept = 1,linetype="dashed",color="black") +
  theme_bw() +
  ylab("Replication index") + xlab("Age (days)") +
  theme(#legend.position = "none",
    strip.text = element_text(face = "italic"),
    panel.grid.minor = element_blank())
irep.plot.sg


##stats on iRep for Schmid and Gill only
rp.l.filt_s <- filter(rp.l.filt_sg, Genus == "Schmidhempelia")
rp.l.filt_g <- filter(rp.l.filt_sg, Genus == "Gilliamella")

##Schmid:

m.s.age_wInt <- lm(iRep ~ Age_days * Colony, data=rp.l.filt_s)
m.s.age_nInt <- lm(iRep ~ Age_days + Colony, data=rp.l.filt_s)
anova(m.s.age_wInt,m.s.age_nInt) # p = 0.86 -> use simpler model w/o intxn
summary(m.s.age_nInt)
qqnorm(resid(m.s.age_nInt))

#comparing levels within Colony
posthoc.colony.test <- glht(m.s.age_nInt, linfct = mcp(Colony = 'Tukey'))
summary(posthoc.colony.test)


##Gill:

m.g.age_wInt <- lm(iRep ~ Age_days * Colony, data=rp.l.filt_g)
m.g.age_nInt <- lm(iRep ~ Age_days + Colony, data=rp.l.filt_g)
anova(m.g.age_wInt,m.g.age_nInt) # p = 0.91 -> use simpler model w/o intxn
summary(m.g.age_nInt)
qqnorm(resid(m.g.age_nInt))

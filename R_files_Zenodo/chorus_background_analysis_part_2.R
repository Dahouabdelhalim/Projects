#Part 2: compare female choices for call rate in tests with vs. without a chorus background

#load data
data<-read.csv("cb_other_combined.csv", header = T)

data.choices<-droplevels(data[data$mate_choice!="NC",]) #just analyze choices for this analysis
data.choices<-droplevels(data.choices[is.na(data.choices$mate_choice)==F,]) #drop NA values
unique(data.choices$ID) #34 individuals have at least 1 measure of choice in absence of chorus background 
#reclass choices as factor (they were imported as character...)
data.choices$mate_choice <- as.factor(data.choices$mate_choice)

#Do 26 vs 37 call/min assays differ in choices from 31 vs 37 calls/min assays? 
#test difference of proportions using fisher's exact test
summary(data.choices[data.choices$treatment=="26v37",]$mate_choice) #17 F, 13S
summary(data.choices[data.choices$treatment=="31v37",]$mate_choice) #9 F, 10S
assay.matrix <- matrix(c(17,13,9,10), nrow = 2, ncol = 2, byrow = T, dimnames = list(c("26v37", "31v37"), c("F", "S")))
fisher.test(assay.matrix) #p = 0.5688, prefs do not differ between these assays.  


#collapsing all tests without a chorus background into 'noChorus' treatment:
data.choices2 <- data.choices
levels(data.choices2$treatment)[1:2] <- "noChorus"

#is there a preference in each treatment for the subset of females that were tested in 
#no chorus AND chorus background tests?
summary(data.choices[data.choices$treatment=="26v37",]$mate_choice) #17 F, 13S
binom.test(x = c(17,13)) #no sig. pref (p = 0.58)
summary(data.choices[data.choices$treatment=="31v37",]$mate_choice) #9 F, 10S
binom.test(x = c(9,10)) #no sig. pref (p = 1)
summary(data.choices2[data.choices2$treatment=="noChorus",]$mate_choice) #26 F, 23S
binom.test(x = c(26,23)) #no sig. pref (p = 0.78)
summary(data.choices2[data.choices2$treatment=="pureSM",]$mate_choice) #10F, 14S
binom.test(x = c(10,14)) #no sig. pref (p = 0.54)
summary(data.choices2[data.choices2$treatment=="mixed",]$mate_choice) #8F, 20S
binom.test(x = c(8,20)) #sig. pref (p = 0.036)


#model effect of chorus treatment on mate choice
library(lme4)
trt.mod <- glmer(formula = mate_choice ~ treatment + (1|ID), data = data.choices2, family = "binomial")

#post-hoc tests: does no-chorus differ from mixed chorus? Does no-chorus differ from pureSM?
library(multcomp)
summary(glht(trt.mod, linfct=mcp(treatment=c("noChorus - mixed = 0",
                                               "noChorus - pureSM = 0"))),
        test = adjusted("none")
        )

#post hoc tests using emmeans
library(emmeans)
emmeans(trt.mod, specs = "treatment", infer =c(T, T), type = "response")
#repeating post-hoc tests with p-value adjustment for 2 comparisons:
summary(glht(trt.mod, linfct=mcp(treatment=c("noChorus - mixed = 0",
                                             "noChorus - pureSM = 0")))
)

# new fig. panel SxA
bar.counts <- c(10,14,8,20,17,13,9,10)
chorus <- c("pure species", "pure species", "mixed species", "mixed species", 
            "no chorus\\26 v. 37", "no chorus\\26 v. 37",
            "no chorus\\31 v. 37", "no chorus\\31 v. 37")
chorus <- factor(chorus, levels = c("pure species", "mixed species", "no chorus\\26 v. 37", "no chorus\\31 v. 37"))
Stimulus <- c("fast", "slow","fast", "slow", "fast", "slow", "fast", "slow")
bar.df <- cbind.data.frame(bar.counts, chorus, Stimulus)

panelA <- ggplot(data = bar.df, aes(x=chorus, y=bar.counts, fill = Stimulus)) +
  geom_bar(stat="identity", color = "black", position = position_dodge())+
  theme_classic()+
  ylab("Number of females choosing stimulus")+
  scale_fill_manual(values = c("darkgrey", "white"))+
  xlab("Chorus background")+
  scale_x_discrete(labels = c("pure species\\n26 v. 37", "mixed species\\n26 v. 37", 
                              "no chorus\\n26 v. 37", "no chorus\\n31 v. 37"))
#new fig. panelSXB
#subset relevant columns so that this is less annoying to do:
data.fig <- data.choices2[,c(3,7,9)]
unique(data.fig$ID) #34 females
yax.displace <- seq(from = -0.16, to = 0.17, by = 0.01) #34 displacements, one for each female

#add displacement column to data.fig according to female ID
yax.jitter <- rep(NA, times = nrow(data.fig))
data.fig <- cbind.data.frame(data.fig, yax.jitter)
for (i in 1:34) {
  fem.ID <- unique(data.fig$ID)[i]
  data.fig[data.fig$ID == fem.ID,4] <- yax.displace[i]
}
as.numeric(data.fig$mate_choice)
#refactor mate choice so that S = 1, F = 2
data.fig$mate_choice <- factor(data.fig$mate_choice, levels = c("S", "F"))
data.fig$yax.val <- as.numeric(data.fig$mate_choice) + data.fig$yax.jitter
#refactor treatment so that no chorus = 1, pure SM = 2, mixed = 3
data.fig$treatment <- factor(data.fig$treatment, levels = c("noChorus", "pureSM", "mixed"))


panelSX <- ggplot(data = data.fig, mapping = aes(x = treatment, y = yax.val, group = ID))+
  #add points showing females' choices
  geom_point(shape = 1, size = 3)+
  geom_line()+
  theme_classic()+
  scale_x_discrete(labels = c("No Chorus", "Pure\\nSpecies", "Mixed\\nSpecies"), 
                     name = "Chorus Treatment")+
  scale_y_continuous(breaks = c(1,2), 
                     labels = c("Slow", "Fast"), name = "Female choice")+
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 14))

#Export fig
svg(filename="figSX.svg", 
    width=4, 
    height=3, 
    pointsize=12, 
    antialias = c("subpixel")
)
panelSX
dev.off()
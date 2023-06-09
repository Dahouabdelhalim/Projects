###############################
#      Process Raw Data       #
#   Catherine Triandafillou   #
###############################



# Setup -------------------------------------------------------------------

path_to_here <- "~/Dropbox (Drummond Lab)/cat/pH-induction-paper/final-version/"
path_to_raw_data <- "~/Dropbox (Drummond Lab)/cat/data/flow/"

library(data.table)
library(flowCore)
library(flowViz)
library(reshape2)
library(tidyr)
library(dplyr)
library(flownalysis)
library(ggplot2)
#setwd(path_to_raw_data)
setwd(path_to_here)

sigmoid <- function(params, x) {
  (params[1] / (1 + exp(-params[2] * (x-params[3])))) + params[4]
}


# Ionophore data ----------------------------------------------------------


# Data from January 09, 2017


raw.jan <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-01-07/", sep = ""),
                                      alter.names = TRUE), rename.cols = FALSE)

jan.buffer.values <- c(5.05, 5.59, 6.12, 6.51, 7.01, 7.52, 8.13)
jan.inputs <- c("5p0", "5p5", "6p0", "6p5", "7p0", "7p5", "8p0")

jan.bkgs.dist <- raw.jan %>%
  filter(grepl("controls", exp)) %>%
  separate(exp, c("condition", "strain", "shock", "media", "rep", "fcs")) %>%
  filter(FSC.A < 70000, FITC.A < 2000) %>%
  mutate(rel.red = PE.Texas.Red.A / FSC.A)

jan.bkgs <- jan.bkgs.dist %>%
  group_by(strain, media) %>%
  summarise(red.bkg = median(rel.red), FITC.bkg = median(FITC.A), BV.bkg = median(BV.510.A))

## January 09 calibration curve

jan.cc <- raw.jan %>%
  filter(grepl("cc", exp)) %>%
  separate(exp, c("condition", "pH", "fcs")) %>%
  mutate(pH.ratio = (BV.510.A - as.numeric(jan.bkgs[1,5])) / (FITC.A - as.numeric(jan.bkgs[1,5])), pH = as.numeric(plyr::mapvalues(pH, from=jan.inputs, to=jan.buffer.values))) %>%
  group_by(pH) %>%
  summarise(med.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))


x <- jan.cc$pH
y <- jan.cc$med.ratio

fitmodel <- nls(y~a/(1 + exp(-b * (x-c))) + d, start=list(a=.5, b=2, c=7, d=0.25))
params = coef(fitmodel)

x2 <- seq(5.0, 8.3, 0.01)
y2 <- sigmoid(params, x2)


q <- ggplot(data=NULL) + geom_point(aes(x=x, y=y)) + geom_line(aes(x=x2, y=y2), size=0.5) 
q <- q + geom_errorbar(data=jan.cc, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05)
q <- q + theme_minimal() + labs(
  x = "pH",
  y = "Ratio 405/488",
  title = "pHluorin calibration curve") + geom_text(aes(x=5.8, y=0.55, label="ratio==frac(a,1+e^(-b(x-c))) + d"), parse = T)
q

ratio.to.pH <- function(ratio, p=params, lims=y2){
  # ratio = [a / (1 + exp(-b(pH-c)))] + d
  # pH = [log(a / (ratio - d) - 1) / -b] + c
  min.ratio = min(y2) - 0.05*min(y2)
  max.ratio = max(y2) + 0.05*max(y2)
  ifelse(ratio < min.ratio, NA, ifelse(ratio < max.ratio, (log((p[1]/(ratio-p[4])) - 1)/-p[2]) + p[3], NA))
}

##

jan.ts <- raw.jan %>%
  filter(!grepl("cc", exp) & !grepl("controls", exp)) %>%
  separate(exp, c("timepoint", "condition", "pH.of.shock")) %>%
  mutate(rel.red = (PE.Texas.Red.A / FSC.A) / as.numeric(jan.bkgs[2,3]), pH.ratio = (BV.510.A - as.numeric(jan.bkgs[1,5])) / (FITC.A - as.numeric(jan.bkgs[1,4])),
         pH.m = ratio.to.pH(pH.ratio),
         timepoint = as.numeric(timepoint),
         date = "jan09")

percent.not.na <- nrow(filter(jan.ts, !is.nan(pH.m))) / nrow(jan.ts) * 100


# Data from January 16, 2017
raw.jan.16 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-01-16/", sep = ""),
                                         alter.names = T), rename.cols = FALSE)

jan.16.bkgs.dist <- raw.jan.16 %>%
  separate(exp, c("sample", "strain", "condition")) %>%
  filter(sample == "control") %>%
  mutate(population = ifelse(FSC.A > 100*FITC.A - 5000, "wt", "pH"), rel.red = PE.Texas.Red.A / FSC.A) %>%
  filter(FSC.A < 70000, FITC.A < 2000)


jan.16.bkgs <- jan.16.bkgs.dist %>%
  group_by(strain, population) %>%
  summarise(med.FITC = median(FITC.A), med.BV = median(BV.510.A), med.Texas = median(PE.Texas.Red.A), med.FSC = median(FSC.A), med.SSC = median(SSC.A), red.bkg = median(rel.red), sample.size = n())


to_join <- jan.16.bkgs %>%
  filter(strain == "mix") %>%
  select(population, red.bkg)

jan.16.ts.dist <- raw.jan.16 %>%
  separate(exp, c("timepoint", "condition", "shock.pH")) %>%
  filter(timepoint != "control") %>%
  mutate(population = ifelse(FSC.A > 100*FITC.A - 5000, "wt", "pH"), pH.ratio = BV.510.A / FITC.A, timepoint = as.numeric(timepoint)) %>%
  left_join(., to_join, by=c("population")) %>%
  select(-strain) %>%
  mutate(rel.red = (PE.Texas.Red.A / FSC.A) / red.bkg, date="jan16")

jan.16.ts <- jan.16.ts.dist %>%
  group_by(timepoint, condition, shock.pH, population) %>%
  summarise(med.red = median(rel.red), red.high = quantile(rel.red, 0.75),
            red.low = quantile(rel.red, 0.25),
            med.pH = median(pH.ratio),
            pH.high = quantile(pH.ratio, 0.75),
            pH.low = quantile(pH.ratio, 0.25),
            sample.size=n())


# Data from January 30


raw.jan.30 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-01-30/", sep = ""),
                                         alter.names = T), rename.cols = F)

jan.30.bkg.dist <- raw.jan.30 %>%
  separate(exp, c("exp", "sample", "condition", "media")) %>%
  filter(exp == "controls") %>%
  mutate(population = ifelse(FSC.A > 110*FITC.A - 5000, "wt", "pH"))


jan.30.bkg <- raw.jan.30 %>%
  separate(exp, c("exp", "sample", "condition", "media")) %>%
  filter(exp == "controls" & media == "sc") %>%
  mutate(population = ifelse(FSC.A > 100*FITC.A - 5000, "wt", "pH"), rel.red = PE.Texas.Red.A / FSC.A) %>%
  group_by(sample, population) %>%
  summarise(med.FITC = median(FITC.A),
            med.BV = median(BV.510.A),
            med.Texas = median(PE.Texas.Red.A),
            med.FSC = median(FSC.A),
            med.SSC = median(SSC.A),
            red.bkg = median(rel.red),
            sample.size = n())

wt.red.bkg <- as.numeric(jan.30.bkg[2,8])
pH.red.bkg <- as.numeric(jan.30.bkg[5,8])


## January 30 calibration curve

jan.30.buffer.values <- c(5.03, 5.57, 6.08, 6.61, 6.98, 7.46, 8.08)

cc.bkgs <- raw.jan.30 %>%
  filter(grepl("controls", exp)) %>%
  separate(exp, c("exp", "sample", "condition", "media")) %>%
  filter(media == "buffer") %>%
  summarise(med.FITC = median(FITC.A),
            med.BV = median(BV.510.A),
            med.Texas = median(PE.Texas.Red.A),
            med.FSC = median(FSC.A),
            med.SSC = median(SSC.A),
            sample.size = n())

jan.30.cc <- raw.jan.30 %>%
  filter(grepl("cc", exp)) %>%
  separate(exp, c("condition", "pH", "fcs")) %>%
  mutate(pH.ratio = (BV.510.A - as.numeric(cc.bkgs[1,2])) / (FITC.A - as.numeric(cc.bkgs[1,1])), pH = as.numeric(plyr::mapvalues(pH, from=jan.inputs, to=jan.30.buffer.values))) %>%
  group_by(pH) %>%
  summarise(med.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))


x <- jan.30.cc$pH
y <- jan.30.cc$med.ratio

fitmodel <- nls(y~a/(1 + exp(-b * (x-c))) + d, start=list(a=.5, b=2, c=7, d=0.25))
params = coef(fitmodel)

x2 <- seq(5.0, 8.3, 0.01)
y2 <- sigmoid(params, x2)


q <- ggplot(data=NULL) + geom_point(aes(x=x, y=y)) + geom_line(aes(x=x2, y=y2), size=0.5) 
q <- q+ geom_errorbar(data=jan.30.cc, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05)
q <- q + theme_minimal() + labs(
  x = "pH",
  y = "Ratio 405/488",
  title = "pHluorin calibration curve\\nJanuary 30 2017") + geom_text(aes(x=5.8, y=0.55, label="ratio==frac(a,1+e^(-b(x-c))) + d"), parse = T)
q

ratio.to.pH.30 <- function(ratio, p=params, lims=y2){
  # ratio = [a / (1 + exp(-b(pH-c)))] + d
  # pH = [log(a / (ratio - d) - 1) / -b] + c
  min.ratio = min(y2) - 0.1*min(y2)
  max.ratio = max(y2) + 0.1*max(y2)
  ifelse(ratio < min.ratio, NA, ifelse(ratio < max.ratio, (log((p[1]/(ratio-p[4])) - 1)/-p[2]) + p[3], NA)) # This is terrible
}


jan.30.ts.dist <- raw.jan.30 %>%
  separate(exp, c("timepoint", "sample", "shock.pH", "treatment")) %>%
  filter(timepoint != "controls" & timepoint != "cc") %>%
  mutate(timepoint = as.numeric(timepoint),
         pH.ratio = (BV.510.A - as.numeric(jan.30.bkg[2,4])) / (FITC.A - as.numeric(jan.30.bkg[2,3])),
         pH.m = ratio.to.pH.30(pH.ratio),
         rel.red = (PE.Texas.Red.A / FSC.A) / pH.red.bkg,
         population = ifelse(BV.421.A > 200, "bv+", ifelse(FSC.A > 100*FITC.A - 5000, "wt", "pH")),
         date="jan30")

percent.not.na.30 <- nrow(filter(jan.30.ts.dist, !is.nan(pH.m))) / nrow(jan.30.ts.dist) * 100

print("Percent of NA values for pH calibration curve January 30")
print(percent.not.na.30)


# Data from February 20

feb.20.raw <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-02-20/", sep = ""),
                                         alter.names = T), rename.cols = FALSE)

feb.bkg.dist <- feb.20.raw %>%
  separate(exp, c("controls", "treatment", "strain", "media")) %>%
  filter(controls == "controls") %>%
  mutate(rel.red = PE.Texas.Red.A / FSC.A,
         pH.ratio = BV.510.A / FITC.A,
         population = ifelse(BV.421.A > 200, "bv+", ifelse(FSC.A > 95*FITC.A - 7000, "wt", "pH")))

feb.bkgs <- feb.bkg.dist %>%
  group_by(strain, media, population) %>%
  summarise(med.FITC = median(FITC.A), med.BV = median(BV.510.A), med.rel.red = median(rel.red), med.pH.ratio = median(pH.ratio), sample.size = n())

# It's only okay that this is here because the ratio.to.pH function from Jan 30 can be used to analyze this data; this assumption is tested in the pH section below
feb.ts.dist <- feb.20.raw %>%
  separate(exp, c("timepoint", "sample", "shock.pH", "treatment")) %>%
  filter(timepoint != "cc" & timepoint != "controls") %>%
  mutate(timepoint = as.numeric(timepoint),
         pH.ratio = (BV.510.A - as.numeric(feb.bkgs[2,5])) / (FITC.A - as.numeric(feb.bkgs[2,4])),
         pH.m = ratio.to.pH.30(pH.ratio),
         rel.red = (PE.Texas.Red.A / FSC.A) / as.numeric(feb.bkgs[9,6]),
         population = ifelse(BV.421.A > 200, "bv+", ifelse(FSC.A > 95*FITC.A - 7000, "wt", "pH")),
         date="feb20")

percent.not.na.feb20 <- nrow(filter(feb.ts.dist, population == "pH" & !is.na(pH.m)))/nrow(filter(feb.ts.dist, population == "pH")) * 100

print("Percent of NA values for pH calibration curve February 20")
print(percent.not.na.feb20)


feb.cc <- feb.20.raw %>%
  separate(exp, c("condition", "pH")) %>%
  filter(condition == "cc") %>%
  mutate(pH.ratio = (BV.510.A - as.numeric(feb.bkgs[2,5])) / (FITC.A - as.numeric(feb.bkgs[2,4])),
         pH = as.numeric(plyr::mapvalues(pH, from=jan.inputs, to=jan.30.buffer.values))) %>%
  group_by(pH) %>%
  summarise(med.pH.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))

q + geom_point(data = feb.cc, aes(x = pH, y = med.pH.ratio), color = "red", alpha = 0.5) +
  geom_errorbar(data=feb.cc, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05, color = "red", alpha=0.5)


feb.ts <- feb.ts.dist %>%
  group_by(timepoint, sample, shock.pH, treatment, population) %>%
  summarise(med.pH = median(pH.m, na.rm = T),
            pH.high = quantile(pH.m, 0.75, na.rm = T),
            pH.low = quantile(pH.m, 0.25, na.rm = T),
            med.red = median(rel.red, na.rm = T),
            red.high = quantile(rel.red, 0.75, na.rm = T),
            red.low = quantile(rel.red, 0.25, na.rm = T),
            n = n())

# Data from August 10, 2017

raw.aug.10 <- read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-08-10/", sep = ""), alter.names = T)
aug.10.df <- merge_flowSet(raw.aug.10) %>%
  rename(PEDazzle.A = PEDazzle594.A)

aug.10.bkgs <- filter(aug.10.df) %>%
  filter(grepl("controls", exp))

aug.10.cc <- filter(aug.10.df) %>%
  filter(grepl("cc", exp))

cc.analysis(aug.10.cc, aug.10.bkgs, "aug10", xmin = 5.0, xmax = 8.35, start.list = list(a=-1, b=-1, c=7.3, d=3))


aug.10.ind.ts <- process.pH.induction(filter(aug.10.df, grepl("hs", exp) | grepl("controls", exp))) %>%
  filter(fitness == "ycgt028") %>%
  mutate(pH.m = convert.to.pH.aug10(FITC.A, BV510.A), population = ifelse(BV421.A > 0.5*FITC.A - 200, "bv", "pH"), date = "aug10")

colnames(aug.10.ind.ts) <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W", "FITC.A", "FITC.H", "FITC.W", "PE.Texas.Red.A", "PE.Texas.Red.H",
                             "PE.Texas.Red.W", "BV.421.A", "BV.421.H", "BV.421.W", "BV.510.A", "BV.510.H", "BV.510.W", "Time", "timepoint", "treatment",
                             "shock.pH", "sample", "rel.red", "pH.m", "population", "date")



# Data from September 28 2017

raw.sept.28 <- read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-09-28/", sep = ""), alter.names = T)
raw.sept.28.df <- merge_flowSet(raw.sept.28)
sept.28.bkgs <- filter(raw.sept.28.df, grepl("controls", exp))
sept.28.cc <- filter(raw.sept.28.df, grepl("cc", exp))

cc.analysis(sept.28.cc, sept.28.bkgs, "sept28", xmax = 8.4)

colnames(raw.sept.28.df) <- gsub("594", "", colnames(raw.sept.28.df))

red.bkg <- filter(sept.28.bkgs, exp == "controls_ycgt028-media.fcs") %>%
  summarise(med.red = median(PEDazzle594.A)) %>%
  as.numeric()

sept.28.ts <- filter(raw.sept.28.df, grepl("hs", exp)) %>%
  separate(exp, c("timepoint", "treatment", "shock.pH", "strain"), convert = T, extra = "drop") %>%
  mutate(pH.m = convert.to.pH.sept28(FITC.A, BV510.A),
         population = ifelse(BV421.A < 0.3*FITC.A - 50, "pH", "bv"),
         rel.red = PEDazzle.A / red.bkg,
         date = "sept28")

nrow(filter(sept.28.ts, !is.na(pH.m))) / nrow(sept.28.ts) # quality of calibration curve


# Data from October 23

raw.oct23 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-10-23/", sep = "")),
                           method = "new")
oct23.bkgs <- filter(raw.oct23, grepl("controls", exp))
oct23.cc <- filter(raw.oct23, grepl("cc", exp))

cc.analysis(oct23.cc, oct23.bkgs, id = "oct23", xmin = 4.9, xmax = 8.5, FITC.thresh = 400, BV.thresh = 400,
            start.list = c(a = .2, b = 2, c = 7, d = 1))


oct23.red.bkg.all <- filter(oct23.bkgs, exp == "controls_ycgt028_media_rep2_004.fcs") %>%
  mutate(pH.m = convert.to.pH.oct23(FITC.A, BV510.A))
#write resting pH dist to file
fwrite(oct23.red.bkg.all, "resting-pH.csv")

oct23.red.bkg <- oct23.red.bkg.all %>%
  summarise(med.red = median(PEDazzle.A / FSC.A)) %>%
  as.numeric()

oct23.ts <- filter(raw.oct23, grepl("hs", exp)) %>%
  separate(exp, c("timepoint", "treatment", "fitness", "shock.pH", "replicate"), extra = "drop", convert = T) %>%
  mutate(rel.red = (PEDazzle.A / FSC.A) / oct23.red.bkg,
         pH.m = convert.to.pH.oct23(FITC.A, BV510.A),
         population = ifelse(BV421.A > 500, "bv", ifelse(FITC.A > 300, "pH", "wt")),
         date = "oct23")


# Data from Nov 14

raw.nov14 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-11-14/", sep = "")), method = "new")

nov14.bkg <- filter(raw.nov14, grepl("controls", exp))

nov14.red.bkg <- filter(nov14.bkg, grepl("ycgt028", exp)) %>%
  mutate(rel.red = PEDazzle.A / FSC.A) %>%
  summarise(med.red = median(rel.red, na.rm = T)) %>%
  as.numeric()

nov14.red.w.bkg <- filter(nov14.bkg, grepl("ycgt028", exp)) %>%
  mutate(rel.red = PEDazzle.A / FSC.W) %>%
  summarise(med.red = median(rel.red, na.rm = T)) %>%
  as.numeric()

cc1 <- filter(raw.nov14, grepl("cc_", exp))
cc.analysis(cc1, nov14.bkg, "nov14.1", FITC.thresh = 450, BV.thresh = 450)

nov14.ind <- filter(raw.nov14, !grepl("cc", exp) & !grepl("controls", exp)) %>%
  separate(exp, c("timepoint", "treatment", "sample", "shock.pH", "replicate"), extra = "drop", convert = TRUE) %>%
  filter(shock.pH == "5p0") %>% # Filter out other non-ionophore experiments done on the same day
  mutate(rel.red = (PEDazzle.A/FSC.A) / nov14.red.bkg,
         rel.red.w = (PEDazzle.A/FSC.W) / nov14.red.w.bkg,
         pH.m = convert.to.pH.nov14.1(FITC.A, BV510.A),
         population = case_when(BV421.A > 500 ~ "bv",
                                FITC.A > 300 ~ "pH",
                                TRUE ~ "wt"),
         date = "nov14")


# Aggregate data

colnames(jan.16.ts.dist)[21] = "treatment"
jan.16.ts.dist$sample = "mix"

jan.ts$date = "jan07"
jan.ts$population = "pH"
colnames(jan.ts)[21] = "treatment"
colnames(jan.ts)[22] = "shock.pH"

colnames(sept.28.ts) <- gsub("BV", "BV.", colnames(sept.28.ts))
colnames(sept.28.ts) <- gsub("PEDazzle", "PE.Texas.Red", colnames(sept.28.ts))

colnames(oct23.ts) <- gsub("BV", "BV.", colnames(oct23.ts))
colnames(oct23.ts) <- gsub("PEDazzle", "PE.Texas.Red", colnames(oct23.ts))
colnames(oct23.ts) <- gsub("fitness", "sample", colnames(oct23.ts))


agg.data <- bind_rows(jan.ts, jan.16.ts.dist, jan.30.ts.dist, feb.ts.dist, aug.10.ind.ts, sept.28.ts, oct23.ts, nov14.ind) %>%
  filter(!grepl("native", shock.pH))
agg.data$replicate[is.na(agg.data$replicate)] <- "rep1"

#fwrite(agg.data, file = paste(path_to_here, "data/ionophore-induction-data.tsv", sep = ""), sep = "\\t", na = "NA")
#fwrite(agg.data, file = paste(path_to_here, "../ionophore-induction-data.tsv", sep = ""), sep = "\\t", na = "NA")
write.csv(agg.data, file = gzfile("all_ionophore_data.csv.gz"))


# Data for recovery in buffered SC pH 7.4

## May 24 2018

raw.may24 <- merge_flowSet(read.flowSet(path = paste0(path_to_raw_data,"Triandafillou_18-05-24/"), alter.names = T))

may24.bkgs <- filter(raw.may24, grepl("controls", exp))

may24.red.bkg <- #filter(raw.may24, exp == "30-mock_native-rep1.fcs" | exp == "30-mock_native-rep2.fcs") %>% # ratio to first mock timepoint
  filter(raw.may24, exp == "controls_ycgt028-media.fcs") %>%
  mutate(rel.red = PEDazzle594.A / FSC.A) %>%
  summarise(med.red = median(rel.red, na.rm = T)) %>%
  as.numeric()

may24.cc <- filter(raw.may24, grepl("cc", exp))

cc.analysis(may24.cc, may24.bkgs, "may24", start.list = c(a = 1.6, b = 2, c = 7, d = 0.25), xmin = 4.6, xmax = 8.6)

may24.ts <- filter(raw.may24, grepl("hs", exp) | grepl("mock", exp)) %>%
  separate(exp, c("timepoint", "treatment", "shock.pH", "replicate"), extra = "drop", convert = T) %>%
  mutate(population = ifelse(FITC.A > 600, "pH", ifelse(FITC.A < 400, "wt", "bv")),
         shock.pH.f = gsub("p", ".", shock.pH),
         shock.pH.f = gsub(".0", "", shock.pH.f),
         pH.m = convert.to.pH.may24(FITC.A, BV510.A),
         rel.red = (PEDazzle594.A / FSC.A) / may24.red.bkg,
         date = "may24")


## June 5 2018

raw.june5 <- merge_flowSet(read.flowSet(path=paste0(path_to_raw_data,"Triandafillou_18-06-05/"), alter.names = T))

june5.bkgs <- filter(raw.june5, grepl("controls", exp))

june5.red.bkg <- filter(raw.june5, exp == "controls_ycgt028-media.fcs") %>%
  mutate(rel.red = PEDazzle594.A / FSC.A) %>%
  summarise(med.red = median(rel.red, na.rm = T)) %>%
  as.numeric()

june5.cc <- filter(raw.june5, grepl("cc", exp))

cc.analysis(june5.cc, june5.bkgs, "june5", start.list = c(a = 1.6, b = 2, c = 7, d = 0.25), xmin = 4.6, xmax = 8.6)



june5.ts <- filter(raw.june5, grepl("hs", exp) | grepl("mock", exp)) %>%
  separate(exp, c("timepoint", "treatment", "replicate", "shock.pH"), extra = "drop", convert = T) %>%
  mutate(population = ifelse(log10(BV421.A) > log10(FITC.A) - 0.5, "bv", "pH"),
         shock.pH.f = gsub("p", ".", shock.pH),
         shock.pH.f = gsub(".0", "", shock.pH.f),
         pH.m = convert.to.pH.june5(FITC.A, BV510.A),
         rel.red = (PEDazzle594.A / FSC.A) / june5.red.bkg,
         date = "june5") %>%
  filter(!is.na(population))


## June 7 2018

raw.june7 <- merge_flowSet(read.flowSet(path = paste0(path_to_raw_data,"Triandafillou_18-06-07/"), alter.names = F), method = "new")

june7.bkgs <- filter(raw.june7, grepl("controls", exp))

june7.red.bkg <- filter(raw.june7, exp == "controls_ycgt028_media_003.fcs") %>%
  mutate(rel.red = PEDazzle.A / FSC.A) %>%
  summarise(med.red = median(rel.red, na.rm = T)) %>%
  as.numeric()

june7.cc <- filter(raw.june7, grepl("cc", exp))

cc.analysis(june7.cc, june7.bkgs, "june7", start.list = c(a = 1.6, b = 2, c = 7, d = 0.25), xmin = 4.6, xmax = 8.6, FITC.thresh = 400, BV.thresh = 400)

june7.ts <- filter(raw.june7, grepl("hs", exp) | grepl("mock", exp)) %>%
  separate(exp, c("timepoint", "treatment", "shock.pH", "replicate"), extra = "drop", convert = T) %>%
  mutate(population = ifelse(FITC.A < 350, "wt", ifelse(log10(BV421.A) < log10(FITC.A) - 0.6, "pH", "bv")),
         shock.pH.f = gsub("p", ".", shock.pH),
         shock.pH.f = gsub(".0", "", shock.pH.f),
         pH.m = convert.to.pH.june7(FITC.A, BV510.A),
         rel.red = (PEDazzle.A / FSC.A) / june7.red.bkg,
         date = "june7") %>% 
  filter(!is.na(population)) %>% # 0.4% not assigned due to low/negative signal in BV421 or FITC
  rename("PEDazzle594.A" = "PEDazzle.A", "PEDazzle594.H"= "PEDazzle.H", "PEDazzle594.W" = "PEDazzle.W") # to be consistent with other datasets


neutral.sc.recovery <- bind_rows(may24.ts, june5.ts, june7.ts)

write.csv(neutral.sc.recovery, file = gzfile("buffered-sc-recovery.csv.gz"))


# Neutral Ionophore data --------------------------------------------------



# Data from April 13

raw.apr.13 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-04-13/", sep = ""),
                                         alter.names = T), rename.cols = T)

apr.13.bkgs.dist <- raw.apr.13 %>%
  filter(grepl("controls", exp)) %>%
  separate(exp, c("experiment", "strain", "background")) %>%
  group_by(strain, background, experiment) %>%
  select(-PE.A, -PE.H, -PE.W) %>%
  mutate(rel.red = PEDazzle594.A / FSC.A)


apr.13.bkgs <- summarise_all(apr.13.bkgs.dist, median)

red.bkg <- as.numeric(filter(apr.13.bkgs, strain == "ycgt0028")$rel.red)
FITC.bkg <- as.numeric(filter(apr.13.bkgs, strain == "by4743" & background == "buffer")$FITC.A)
BV.bkg <- as.numeric(filter(apr.13.bkgs, strain == "by4743" & background == "buffer")$BV510.A)

FITC.bkg.sample <- as.numeric(filter(apr.13.bkgs, strain == "by4743" & background == "media")$FITC.A)
BV.bkg.sample <- as.numeric(filter(apr.13.bkgs, strain == "by4743" & background == "media")$BV510.A)

# April 13 Calibration Curve

buffer.values <- c(4.97, 5.5, 6.08, 6.45,7.05, 7.46, 8.04)
inputs <- c("5p0", "5p5", "6p0", "6p5", "7p0", "7p5", "8p0")

cc <- filter(raw.apr.13, grepl("cc", exp)) %>%
  filter(FITC.A > 800 & BV510.A > 800) %>% # High-quality points only, preserves 99% of data
  separate(exp, c("experiment", "pH")) %>%
  mutate(pH = as.numeric(plyr::mapvalues(pH, from=inputs, to=buffer.values)), pH.ratio = (BV510.A - BV.bkg) / (FITC.A - FITC.bkg)) %>%
  group_by(pH) %>%
  summarise(med.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))


x <- cc$pH
y <- cc$med.ratio

fitmodel <- nls(y~a/(1 + exp(-b * (x-c))) + d, start=list(a=.5, b=2, c=7, d=0.25))
params = coef(fitmodel)

x2 <- seq(4.9, 8.2, 0.01)
y2 <- sigmoid(params, x2)


q <- ggplot(data=NULL) + geom_point(aes(x=x, y=y)) + geom_line(aes(x=x2, y=y2), size=0.5) 
q <- q+ geom_errorbar(data=cc, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05)
q <- q + theme_minimal() + labs(
  x = "pH",
  y = "Ratio 405/488",
  title = "pHluorin calibration curve") + geom_text(aes(x=5.8, y=0.55, label="ratio==frac(a,1+e^(-b(x-c))) + d"), parse = T)
q

ratio.to.pH.apr13 <- function(ratio, p=params, lims=y2){
  # ratio = [a / (1 + exp(-b(pH-c)))] + d
  # pH = [log(a / (ratio - d) - 1) / -b] + c
  min.ratio = min(y2) - 0.1*min(y2)
  max.ratio = max(y2) + 0.1*max(y2)
  ifelse(ratio < min.ratio, NA, ifelse(ratio < max.ratio, (log((p[1]/(ratio-p[4])) - 1)/-p[2]) + p[3], NA)) # This is terrible
}

apr.13.ts.dist <- raw.apr.13 %>%
  filter(grepl("spike", exp)) %>%
  separate(exp, c("timepoint", "treatment", "fit", "shock.pH")) %>%
  select(-PE.A, -PE.H, -PE.W) %>%
  mutate(rel.red = (PEDazzle594.A / FSC.A) / red.bkg,
         population = ifelse(BV421.A > 450, "bv", ifelse(SSC.A > 30*FITC.A - 7000, "wt", "pH")),
         pH.ratio = ifelse(population == "pH", (BV510.A - BV.bkg.sample) / (FITC.A - FITC.bkg.sample), NA),
         pH.m = ratio.to.pH.apr13(pH.ratio),
         timepoint = as.numeric(timepoint),
         date="apr13")


# Data from April 17

raw.17 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-04-17/", sep = ""),
                                     alter.names = T))
raw.17.2 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-04-17/", sep = ""),
                                       alter.names = T))

raw.17.df <- full_join(raw.17, raw.17.2) #%>% select(-PE.A, -PE.H, -PE.W)

bkgs.17 <- filter(raw.17.df, grepl("controls", exp)) %>%
  separate(exp, c("experiment", "strain", "background", "condition", "fcs")) %>%
  mutate(rel.red = PEDazzle594.A / FSC.A) %>%
  group_by(strain, background, condition) %>%
  summarise(red.bkg = median(rel.red), FITC.bkg = median(FITC.A), BV.bkg = median(BV510.A))

BV.bkg.17 = as.numeric(bkgs.17[1,6])
FITC.bkg.17 = as.numeric(bkgs.17[1,5])

BV.bkg.media.17 = as.numeric(bkgs.17[2,6])
FITC.bkg.media.17 = as.numeric(bkgs.17[2,5])

red.bkg.17 = as.numeric(bkgs.17[3,4])

## April 17 calibration curve

buffer.values.17 <- c(5.07, 5.53, 6.04, 6.58, 7.09, 7.50, 8.07, 8.34)
inputs.17 <- c("5p0", "5p5", "6p0", "6p5", "7p1", "7p5", "8p0", "8p3")

cc.17 <- filter(raw.17.df, grepl("cc", exp)) %>%
  filter(FITC.A > 800 & BV510.A > 800) %>% # High-quality points only, preserves 99% of data
  separate(exp, c("experiment", "pH")) %>%
  filter(pH != "8p3") %>% # 8.3 sample was very noisy and is excluded from the calibration curve
  mutate(pH = as.numeric(plyr::mapvalues(pH, from=inputs.17, to=buffer.values.17)), pH.ratio = (BV510.A - BV.bkg.17) / (FITC.A - FITC.bkg.17)) %>%
  group_by(pH) %>%
  summarise(med.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))

x <- cc.17$pH
y <- cc.17$med.ratio

fitmodel <- nls(y~a/(1 + exp(-b * (x-c))) + d, start=list(a=.5, b=2, c=7, d=0.25))
params = coef(fitmodel)

x2 <- seq(4.9, 8.4, 0.01)
y2 <- sigmoid(params, x2)


q <- ggplot(data=NULL) + geom_point(aes(x=x, y=y)) + geom_line(aes(x=x2, y=y2), size=0.5) 
q <- q+ geom_errorbar(data=cc.17, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05)
q <- q + theme_minimal() + labs(
  x = "pH",
  y = "Ratio 405/488",
  title = "pHluorin calibration curve") + geom_text(aes(x=5.8, y=0.55, label="ratio==frac(a,1+e^(-b(x-c))) + d"), parse = T)
q

ratio.to.pH.17 <- function(ratio, p=params, lims=y2){
  # ratio = [a / (1 + exp(-b(pH-c)))] + d
  # pH = [log(a / (ratio - d) - 1) / -b] + c
  min.ratio = min(y2)
  max.ratio = max(y2)
  ifelse(ratio < min.ratio, NA, ifelse(ratio < max.ratio, (log((p[1]/(ratio-p[4])) - 1)/-p[2]) + p[3], NA)) # This is terrible
}


apr.17.ts.dist <- raw.17.df %>%
  filter(grepl("spike", exp)) %>%
  separate(exp, c("timepoint", "treatment", "fit", "shock.pH")) %>%
  mutate(timepoint = as.numeric(timepoint),
         pH.ratio = (BV510.A - BV.bkg.media.17) / (FITC.A - FITC.bkg.media.17),
         pH.m = ratio.to.pH.17(pH.ratio),
         population=ifelse(BV421.A > 2500, "bv", ifelse(BV510.A > 1200, "pH", "wt")),
         rel.red = (PEDazzle594.A / FSC.A) / red.bkg,
         date = "apr17") %>%
  filter(rel.red < 80) # Preserves 99.9% of data


# Data from August 3

raw.aug.03 <- read.flowSet(path=paste(path_to_raw_data,"Triandafillou_17-08-03/",sep=""), alter.names = F)
raw.aug.03.df <- merge_flowSet(raw.aug.03, method="new")
labeled.aug.03.df <- raw.aug.03.df %>% 
  mutate(population = ifelse(BV421.A > 600, "bv", ifelse(FITC.A > 270, "pH", "wt")))

aug.03.cc <- filter(raw.aug.03.df, grepl("cc", exp))
aug.03.bkgs <- filter(raw.aug.03.df, grepl("controls", exp))
cc.analysis(aug.03.cc, aug.03.bkgs, id = "aug03", xmin = 4.8, xmax = 8.7, FITC.thresh = 400, BV.thresh = 400)



aug.03.ind.ts.dist <- process.pH.induction(labeled.aug.03.df, summarise = F) %>%
  mutate(date = "aug03", pH.m = convert.to.pH.aug03(FITC.A, BV510.A), replicate = "rep1")
colnames(aug.03.ind.ts.dist) <- gsub("fitness", "fit", colnames(aug.03.ind.ts.dist))



agg.data.neutral <- bind_rows(apr.13.ts.dist, apr.17.ts.dist, aug.03.ind.ts.dist)


#fwrite(agg.data.neutral, file = paste(path_to_here, "data/near-seven-ionophore-data.tsv"), sep = "\\t", na = "NA")
write.csv(agg.data.neutral, file = gzfile("neutral_ionophore_data.csv.gz"))


# Native data -------------------------------------------------------------




# Native shock: Data from 3/09/2017

raw.native <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-03-07/", sep = "")),
                            method = "new")

pH.labels <- c("5p0", "6p5", "7p0", "7p5", "8p0")
mar.07.pH.values <- c(4.97, 6.45, 7.05, 7.44, 8.04)

bkgs.dist <- filter(raw.native, grepl("controls", exp)) %>%
  separate(exp, c("controls", "condition", "strain", "media", "sample.no", "fcs")) %>%
  filter(media != "voltages") %>%
  select(-controls, -condition, -sample.no, -fcs, -Time) %>%
  mutate(rel.red = PEDazzle.A / FSC.A)

bkgs <- bkgs.dist %>%
  group_by(strain, media) %>%
  summarise_all(median)

red.bkg <- as.numeric(bkgs[4,21])

induction.ts.dist <- raw.native %>%
  separate(exp, c("timepoint", "treatment", "shock.pH", "sample", "rep", "sample.no", "fcs")) %>%
  filter(timepoint != "controls" & timepoint != "cc") %>%
  mutate(timepoint = as.numeric(timepoint),
         rel.red = (PEDazzle.A / FSC.A) / red.bkg,
         population = ifelse(FSC.A > 120*FITC.A - 15000, "wt", "pH")) %>%
  filter(shock.pH == "native") %>%
  select(-fcs, -sample.no, -rep) %>%
  mutate(date = "mar07")



# Native shock: Data from 5/03/2017
raw.df <- merge_flowSet(read.flowSet(path=paste(path_to_raw_data, "Triandafillou_17-05-03/", sep = ""), alter.names = T))

bkgs.dist <- raw.df %>%
  filter(grepl("controls", exp)) %>%
  separate(exp, c("experiment", "strain", "background"), extra="drop", convert = TRUE) %>%
  group_by(strain, background, experiment) %>%
  mutate(rel.red = PEDazzle594.A / FSC.A)

bkgs <- summarise_all(bkgs.dist, median)

red.bkg <- as.numeric(filter(bkgs, strain == "ycgt028")$rel.red)
FITC.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "buffer")$FITC.A)
BV.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "buffer")$BV510.A)

FITC.bkg.sample <- as.numeric(filter(bkgs, strain == "by4743" & background == "media")$FITC.A)
BV.bkg.sample <- as.numeric(filter(bkgs, strain == "by4743" & background == "media")$BV510.A)

may.03.cc <- filter(raw.df, grepl("cc", exp))
may.03.bkgs <- filter(raw.df, grepl("controls", exp))
may.03.inputs <- c("5p0", "5p5", "6p0", "6p5", "7p1", "7p5", "8p0", "8p6")
may.03.values <- c(5.07, 5.55, 6.04, 6.58, 7.09, 7.50, 8.07, 8.58) # Double check these values 5.7.17!!

cc.analysis(may.03.cc, may.03.bkgs, "may.03", may.03.inputs, may.03.values)

induction <- filter(raw.df, grepl("spike", exp)) %>%
  separate(exp, c("fit", "treatment", "timepoint"), extra="drop", convert = T) %>%
  mutate(rel.red = (PEDazzle594.A / FSC.A) / red.bkg,
         pH.m = convert.to.pH.may.03(FITC.value=FITC.A, BV.value=BV510.A),
         population = ifelse(BV421.A > 2200, "bv", ifelse(BV510.A > 1000, "pH", "wt"))) %>%
  group_by(fit, treatment, population, timepoint) %>%
  summarise(med.pH = median(pH.m, na.rm = T),
            pH.high = quantile(pH.m, 0.75, na.rm = T),
            pH.low = quantile(pH.m, 0.25, na.rm = T),
            med.red = median(rel.red),
            red.high = quantile(rel.red, 0.75),
            red.low = quantile(rel.red, 0.25),
            n.cells = n())


induction.dist <- filter(raw.df, grepl("spike", exp)) %>%
  separate(exp, c("fit", "treatment", "timepoint"), extra="drop", convert = T) %>%
  mutate(rel.red = (PEDazzle594.A / FSC.A) / red.bkg,
         pH.m = convert.to.pH.may.03(FITC.value=FITC.A, BV.value=BV510.A),
         population = ifelse(BV421.A > 2200, "bv", ifelse(BV510.A > 1000, "pH", "wt")),
         date = "may03")



# Same day, pH recovery timecourse:

#pH.ts <- raw.df %>%
#  filter(grepl("classification", exp)) %>%
#  separate(exp, c("experiment", "strain", "treatment", "timepoint"), extra="drop", convert=T) %>%
#  mutate(rel.red = (PEDazzle594.A / FSC.A) / red.bkg,
#         pH.m = convert.to.pH.may.03(FITC.value=FITC.A, BV.value=BV510.A),
#         population = ifelse(BV421.A > 2200, "bv", ifelse(BV510.A > 1000, "pH", "wt"))) %>%
#  filter(population == "pH") %>%
#  mutate(timepoint.f = as.factor(timepoint), treatment.f = factor(treatment, levels = c("unshocked", "shocked")))

#fwrite(pH.ts, "../../pH-induction-paper/17-05-03_recovery-timecourse.csv")



# Data from July 29:

raw.jul.29 <- read.flowSet(path=paste(path_to_raw_data, "Triandafillou_17-07-29/", sep = ""))
jul.29.df <- merge_flowSet(raw.jul.29, method = "new")

jul.29.cc <- filter(jul.29.df, grepl("cc", exp))
jul.29.bkgs <- filter(jul.29.df, grepl("controls", exp))

cc.analysis(jul.29.cc, jul.29.bkgs, "july29", xmin = 4.9, xmax=8.5, FITC.thresh = 300, BV.thresh = 700)

jul.29.red.bkg <- filter(jul.29.df, exp == "controls_ycgt028_media_003.fcs") %>%
  summarise(med.red = median(PEDazzle.A / FSC.A, na.rm = T)) %>%
  as.numeric()

jul.29.ind <- filter(jul.29.df, !grepl("controls", exp) & !grepl("42c10", exp) & !grepl("cc", exp)) %>%
  separate(exp, c("timepoint", "strain"), convert = T, extra = "drop") %>%
  mutate(treatment = "hs", fitness = NA,
         population = gsub("yct028", "pH", strain),
         population = gsub("by4743", "wt", population),
         replicate = "rep1",
         date = "jul29",
         rel.red = (PEDazzle.A / FSC.A) / jul.29.red.bkg,
         pH.m = convert.to.pH.july29(FITC.A, BV510.A))



# Native induction w/ replicates: data from 8/04


raw.aug.04 <- read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-08-04/", sep = ""))
raw.aug.04.df <- merge_flowSet(raw.aug.04, method = "new")
aug.04.bkgs <- filter(raw.aug.04.df, grepl("controls", exp))
aug.04.cc <- filter(raw.aug.04.df, grepl("cc", exp))
cc.analysis(aug.04.cc, aug.04.bkgs, "aug.04", xmin = 5.0, xmax = 7.85, FITC.thresh = 500, BV.thresh = 500)

shock.ts <- filter(raw.aug.04.df, grepl("native", exp) | grepl("controls", exp)) %>%
  mutate(population = ifelse(BV421.A > 1300, "bv", ifelse(FITC.A > 250, "pH", "wt")),
         pH.m = convert.to.pH.aug.04(FITC.A, BV510.A),
         date = "aug04")

aug.04.ind.ts.dist <- process.pH.induction(shock.ts)



# Native induction with replicates and fitness:

raw.oct19 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-10-19/", sep = ""),
                                        alter.names = T), method = "old")
oct19.bkgs <- filter(raw.oct19, grepl("controls", exp))
oct19.cc <- filter(raw.oct19, grepl("cc", exp))

cc.analysis(oct19.cc, oct19.bkgs, id = "oct19", xmin = 4.9, xmax = 8.5, FITC.thresh = 600, BV.thresh = 600,
            start.list = c(a = 1, b = 2, c = 7, d = .5))

oct19.red.bkg <- filter(oct19.bkgs, exp == "controls_ycgt028-media.fcs") %>%
  summarise(med.red = median(PEDazzle594.A / FSC.A)) %>%
  as.numeric()

oct19.ts <- filter(raw.oct19, grepl("hs", exp)) %>%
  separate(exp, c("timepoint", "treatment", "shock.pH", "fitness", "replicate"), extra = "drop", convert = T) %>%
  mutate(rel.red = (PEDazzle594.A / FSC.A) / oct19.red.bkg,
         pH.m = convert.to.pH.oct19(FITC.A, BV510.A),
         population = ifelse(BV421.A > 300, "bv", ifelse(FITC.A > 500, "pH", "wt")),
         date = "oct19")


# Aggregate data

colnames(induction.dist) <- gsub("fit", "fitness", colnames(induction.dist))
colnames(induction.dist) <- gsub("594", "", colnames(induction.dist))
colnames(oct19.ts) <- gsub("594", "", colnames(oct19.ts))
induction.dist$shock.pH = "native"
colnames(induction.ts.dist) <- gsub("sample", "fitness", colnames(induction.ts.dist))

all.native <- bind_rows(induction.dist, induction.ts.dist, aug.04.ind.ts.dist, oct19.ts, jul.29.ind)
all.native["replicate"][is.na(all.native["replicate"])] <- "rep1"

#fwrite(all.native, paste(path_to_raw_data, "data/no-ionophore-induction.tsv"), sep = "\\t")
write.csv(all.native, file = gzfile("no_ionophore_data.csv.gz"))


# Dynamics data -----------------------------------------------------------



# pH dynamics: Data from July 29


together <- filter(jul.29.df, grepl("shock", exp) | grepl("recovery", exp)) %>%
  separate(exp, c("shock", "strain", "timepoint", "volume"), convert = T, extra = "drop") %>%
  mutate(pH.m = convert.to.pH.july29(FITC.A, BV510.A), timepoint = factor(timepoint, levels = c("shock", "recovery"))) %>%
  group_by(timepoint, volume) %>%
  sample_n(5000)


# pH dynamics: Data from October 16

raw.oct16 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_17-10-16/", sep = ""),
                                        alter.names = T))

bkgs <- filter(raw.oct16, grepl("controls", exp))
cc <- filter(raw.oct16, grepl("cc", exp))

cc.analysis(cc, bkgs, "oct16")

hs.ts <- filter(raw.oct16, grepl("hs-ts", exp)) %>%
  separate(exp, c("hs", "ts", "experiment", "media", "rep"), extra = "drop", convert = T) %>%
  mutate(pH.m = convert.to.pH.oct16(FITC.A, BV510.A)) %>%
  select(-hs, -ts, -experiment)








# QC data -----------------------------------------------------------------

# Data from Jan 30, 2018:

raw.jan.30 <- merge_flowSet(read.flowSet(path = paste(path_to_raw_data, "Triandafillou_18-01-30", sep = "")), method = "new")

jan.30.bkgs <- filter(raw.jan.30, grepl("controls", exp))

jan.30.red.bkg <- filter(jan.30.bkgs, grepl("ycgt028", exp)) %>%
  mutate(rel.red = PEDazzle.A / FSC.A) %>%
  summarise(med.red = median(rel.red, na.rm = T)) %>%
  as.numeric()

cc1 <- filter(raw.jan.30, grepl("cc_", exp))
cc.analysis(cc1, jan.30.bkgs, "jan30.1", FITC.thresh = 450, BV.thresh = 450, xmin = 4.8, xmax = 8.7)

iono.plus.heat <- filter(raw.jan.30, grepl("equilib", exp) | grepl("midstress", exp) | grepl("poststress", exp)) %>%
  separate(exp, c("timepoint", "time", "buffer.pH"), extra = "drop") %>%
  filter(FITC.A > 450 & BV510.A > 450) %>% # Filtering for high-quality cells, retains 97% of data
  mutate(pH.m = convert.to.pH.jan30.1(FITC.A, BV510.A),
         buffer.pH.num = as.numeric(gsub("p", ".", buffer.pH)))

reference.dist <- filter(jan.30.bkgs, grepl("ycgt028", exp)) %>%
  filter(FITC.A > 450 & BV510.A > 450) %>% # Filtering for high-quality cells
  separate(exp, c("timepoint", "strain", "background")) %>%
  mutate(pH.m = convert.to.pH.jan30.1(FITC.A, BV510.A))

ionophore.qc.data <- full_join(iono.plus.heat, reference.dist)

#fwrite(ionophore.qc.data, file = paste(path_to_here, "ionophore-qc-data.tsv", sep = ""), sep = "\\t", na = "NA")
write.csv(ionophore.qc.data, file = gzfile("ionophore_qc_data.csv.gz"))

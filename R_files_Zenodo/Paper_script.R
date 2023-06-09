#### This R-script contains the complete data analysis of the paper
#### Peplau, T., Schroeder, J., Gregorich, E., Poeplau, C. (2021): 
#### Long-term geothermal warming affects carbon, but not nitrogen, in a subarctic forest soil ####

library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(reshape2)
library(gridExtra)
library(vegan)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colorblind friendly color palette
setwd() # set to your working directory with the data files

results_HS <- read_excel("HS_data_final.xlsx", 
                      sheet = "sample_data") # the dataset of the soil samples ...
plots_HS <- read_excel("HS_data_final.xlsx", 
                       sheet = "plot_data") # ... and plot-wise data ...

teabags_HS <- read_excel("teabags_HS.xlsx") # the dataset of the teabags

##### Temperature data (Tinytags)--------
Sys.setlocale(locale = "US") # sets locale to US so that American english is used for abbreviations of the months in all plots
clist <- c("W1_10cm", "W1_50cm", "W2_10cm", "W2_50cm", "W3_10cm", "W3_50cm", "W4_10cm", "W4_50cm") # List with the sheet names from excel
sites_raw <- data_frame() # New dataframe for the complete data

for (i in clist) {
  sites_raw <- rbind(sites_raw, read_excel("temperature.xlsx"
                                           ,sheet = i, col_types = c("date","numeric","text")))
}                   # The for-loop iterates through the excel sheets and puts all data togehter, using the rbind-function

climate_daily_takhini <- read_csv("climate-daily-takhini.txt") # Daily climate data. Downloaded from https://climatedata.ca/download/ --> Station data --> Whitehorse Auto with the same dates as the tinytags

# Setting up the Hampel-filter
HampelFilter <- function (x, k,t0){
  n <- length(x)
  y <- x
  L <- 1.4826                                      # Estimate standard deviation of a Gaussian normal distribution 
  for (i in (k + 1):(n - k)) {
    x0 <- median(x[(i - k):(i + k)])               # Median calculation
    S0 <- L * median(abs(x[(i - k):(i + k)] - x0)) # MAD scale estimate
    if (abs(x[i] - x0) > t0 * S0) {
      y[i] <- x0
    }
  }
  data.frame(y)
}                   # The Hampel filter applies a moving window with the size k to the data x and filters the data
# for the threshold of the t-fold MAD (Median Absolute Deviation)

enddate = "2020-06-30"   # "cut-off date" of the timeseries. The loggers were dug out after the cod
startdate = "2019-07-01" # startdate is chosen to minimize the extreme values at the middle of june

sites_filter <- data.frame(sites_raw$Time,
                           Temperature = (HampelFilter(sites_raw$Temperature, 8,t0 = 2)),
                           sites_raw$Class)

sites <- (filter(sites_filter ,sites_raw.Time<enddate & sites_raw.Time>startdate) %>%
            ddply(. ,.(sites_raw.Time,sites_raw.Class), summarize, Temperature=mean(y)))

# After using the Hampel filter and cutting the time-series at cod, the ddply-function calculates
# day mean values for each day

colnames(sites) <- c("Time",
                     "Class",
                     "Temperature")


# Descriptive statistics

stat <- data.frame()

for (i in clist) {
  stat[i,1] <- filter(sites, Class == i) %>%
    summarise(. ,mean(Temperature))
  stat[i,2] <- round(filter(sites ,Temperature<0&
                              Class == i) %>%
                       nrow()/12, digits = 0) 
  stat[i,3] <- filter(sites ,Temperature<0&
                        Class == i) %>%
    summarise(. , mean(Temperature))
  stat[i,4] <- round(filter(sites ,Temperature>0&
                              Class == i) %>%
                       nrow()/12, digits = 0)
  stat[i,5] <- filter(sites ,Temperature>0&
                        Class == i) %>%
    summarise(. , mean(Temperature))
  stat[i,6] <- filter(sites ,Class == i) %>%
    summarise(. , min(Temperature))
  stat[i,7] <- filter(sites ,Class == i) %>%
    summarise(. , max(Temperature))
}  
colnames(stat) <- c("Total mean", 
                    "Days <0°C", 
                    "Mean <0°C", 
                    "Days >0°C", 
                    "Mean >0°C", 
                    "Minimum Temperature", 
                    "Maximum Temperature")

# Plotting the data
# Subplot 1: Temperature in 10 cm
Topsoil <- ggplot() +
  geom_line(data=(filter(sites, Class == "W1_10cm" | Class == "W2_10cm" | Class == "W3_10cm" | Class == "W4_10cm"  )), aes(x=Time, y=Temperature, color=Class))+
  geom_line(data=climate_daily_takhini, aes(x=climate_daily_takhini$LOCAL_DATE, 
                                            y=climate_daily_takhini$MEAN_TEMPERATURE),
            linetype="dotted")+
  ylab("Temperature in °C")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(name="Plot:",
                     values=c("#000000", "#009E73", "#56B4E9","#E69F00"),
                     labels=c("W1", "W2", "W3", "W4"))+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        axis.ticks.x= element_blank())+
  theme_bw()

# Subplot 2: Temperature in 50 cm
Subsoil <- ggplot() +
  geom_line(data=(filter(sites, Class == "W1_50cm" | Class == "W2_50cm" | Class == "W3_50cm" | Class == "W4_50cm"  )), aes(x=Time, y=Temperature, color=Class))+
  geom_line(data=climate_daily_takhini, aes(x=climate_daily_takhini$LOCAL_DATE, 
                                            y=climate_daily_takhini$MEAN_TEMPERATURE),
            linetype="dotted")+
  xlab("Time ")+
  ylab("Temperature in °C")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(name="Plot:",
                     values=c("#000000", "#009E73", "#56B4E9","#E69F00"),
                     labels=c("W1", "W2", "W3", "W4"))+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+
  theme_bw()

# arrange both plots to one
ggarrange(Topsoil, Subsoil, nrow=2, common.legend = TRUE, legend="bottom")

##### General C and N----------
### TOC
stats_conc <- results_HS %>% 
  group_by(Land_Use,Depth) %>%
  summarise(mean_TOC = mean(`TOC [%]`, na.rm=TRUE)*10,  # calculates the mean of each group
            sd_TOC = sd(`TOC [%]`, na.rm =TRUE)*10, # calculates the standard deviation of each group
            n_TOC = n(),  # calculates the sample size per group
            SE_TOC = sd(`TOC [%]`)/sqrt(n())*10) # calculates the standard error of each group


TOC = ggplot(data = stats_conc, aes(color = Land_Use))+
  geom_point(aes(x=mean_TOC, y=Depth))+
  geom_path((aes(x=mean_TOC, y=Depth, group = Land_Use)))+
  geom_errorbarh(aes(y=Depth,xmin=mean_TOC-SE_TOC, xmax=mean_TOC+SE_TOC, height = 2))+
  scale_color_manual(name="Warming",
                     values=c("#000000", "#009E73", "#56B4E9","#E69F00" ),
                     labels=c("+7.7 °C","+1.8 °C", "+0.5 °C", "Reference"))+
  scale_y_reverse(name = "Depth [cm]",
                  breaks=c(10,30,50,70))+
  xlab(bquote('SOC content ['*'g '*kg^-1*']'))+
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 39.5, y=60, label = "A")+
  guides(color = guide_legend(reverse = TRUE))


### C stocks
stats_C_stock <- results_HS %>% 
  group_by(Land_Use,Depth) %>%
  summarise(mean_stock = mean(C_cumu_uncorrected, na.rm=TRUE),  # calculates the mean of each group
            sd_stock = sd(C_cumu_uncorrected, na.rm =TRUE), # calculates the standard deviation of each group
            n_stock = n(),  # calculates the sample size per group
            SE_stock = sd(C_cumu_uncorrected)/sqrt(n())) # calculates the standard error of each group


C_stock = ggplot(data = stats_C_stock, aes(color = Land_Use))+
  geom_point(aes(x=mean_stock, y=Depth))+
  geom_path((aes(x=mean_stock, y=Depth)))+
  geom_errorbarh(aes(y=Depth,xmin=mean_stock-SE_stock, xmax=mean_stock+SE_stock, height = 2))+
  scale_y_reverse(breaks=c(0,20,40,60,80))+
  scale_color_manual(name="Warming",
                     values=c("#000000", "#009E73", "#56B4E9","#E69F00" ),
                     labels=c("+7.7 °C","+1.8 °C", "+0.5 °C", "Reference"))+
  xlab(bquote('Cumulative SOC stock ['*'Mg '*ha^-1*']'))+
  ylab("Depth [cm]")+
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 60, y=60, label = "C")+
  guides(color = guide_legend(reverse = TRUE))


### C:N ratio

stats_cn <- results_HS %>%
  group_by(Land_Use,Depth) %>%
  summarise(mean_CN = mean(`C/N`, na.rm=TRUE), # calculates the mean of each group
            sd_CN = sd(`C/N`, na.rm=TRUE), # calculates the standard deviation of each group
            n_CN = n(), # calculates the sample size per group
            SE_CN = sd(`C/N`)/sqrt(n())) # calculates the standard error of each group

cn =ggplot(data= stats_cn, aes(color = Land_Use))+
  geom_point(aes(x = mean_CN, y=Depth))+
  geom_path(aes(x = mean_CN, y=Depth, group = Land_Use))+
  geom_errorbarh(aes(y=Depth,xmin=mean_CN-SE_CN, xmax=mean_CN+SE_CN, height = 2))+
  scale_color_manual(name="Warming",
                     values=c("#000000", "#009E73", "#56B4E9","#E69F00" ),
                     labels=c("+7.7 °C","+1.8 °C", "+0.5 °C", "Reference"))+
  scale_y_reverse(name = "Depth [cm]",
                  breaks=c(10,30,50,70))+  
  xlab("C:N")+
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 22.5, y=60, label = "B")+
  guides(color = guide_legend(reverse = TRUE))

### N stocks
stats_N_stock <- results_HS %>% 
  group_by(Land_Use,Depth) %>%
  summarise(mean_stock = mean(N_cumu_uncorrected, na.rm=TRUE),  # calculates the mean of each group
            sd_stock = sd(N_cumu_uncorrected, na.rm =TRUE), # calculates the standard deviation of each group
            n_stock = n(),  # calculates the sample size per group
            SE_stock = sd(N_cumu_uncorrected)/sqrt(n())) # calculates the standard error of each group


N_stock = ggplot(data = stats_N_stock, aes(color = Land_Use))+
  geom_point(aes(x=mean_stock, y=Depth))+
  geom_line((aes(x=mean_stock, y=Depth)))+
  geom_errorbarh(aes(y=Depth,xmin=mean_stock-SE_stock, xmax=mean_stock+SE_stock, height = 2))+
  scale_y_reverse(breaks=c(0,20,40,60,80))+
  scale_color_manual(name="Warming",
                     values=c("#000000", "#009E73", "#56B4E9","#E69F00" ),
                     labels=c("+7.7 °C","+1.8 °C", "+0.5 °C", "Reference"))+
  xlab(bquote('Cumulative N stock ['*'Mg '*ha^-1*']'))+
  ylab("Depth [cm]")+
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 4.3, y=60, label = "D")+
  guides(color = guide_legend(reverse = TRUE))


#### Putting the subplots together---
p1 <- ggarrange(TOC,cn, C_stock, N_stock, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 
##### Plotting of the fractionation data  -----------------
Depthlabs <- c(
  '5'="0 - 10 cm",
  '15'="10 - 20 cm",
  '30'="20 - 40 cm",
  '50'="40 - 60 cm",
  '70'="60 - 80 cm"
) # preparing labels for the facet plot

# C fractions
# Proportional c content of the fractions, normalised on 100%
dfA <- select(results_HS, c(warming, Depth, POM_C_relative, DOC_relative, SA_C_relative, SC_C_relative, rSOC_relative)) # selects the relevant columns from the dataset

DF1A <- melt(dfA, id.var=c("warming", "Depth")) # converts the subset from wide to long format
DF2A <- aggregate(x=DF1A$value,
                  by=list(DF1A$warming, DF1A$Depth, DF1A$variable),
                  FUN = mean, na.rm=TRUE) #aggregates the data for plotting

relative <- DF2A %>%
  ggplot(., aes(x = as.factor(Group.1), y = x, fill = Group.3)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~Group.2, labeller = as_labeller(Depthlabs), ncol = 1)+
  scale_fill_manual(name="Fractions:",
                    values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                    labels=c("POM", "DOC", "S+A","S+C", "rSOC"))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  ylab("Proportions of C-Fractions")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle(label = "Relative C content")

#Absolute C content of the fractions in g/kg
dfB <- select(results_HS, c(warming, Depth, POM_C_gkg, DOC_gkg, SA_C_gkg, SC_C_gkg, rSOC_gkg)) # selects the relevant columns from the dataset

DF1B <- melt(dfB, id.var=c("warming", "Depth")) # converts the subset from wide to long format
DF2B <- aggregate(x=DF1B$value,
                  by=list(DF1B$warming, DF1B$Depth, DF1B$variable),
                  FUN = mean, na.rm=TRUE) # aggregates the data for plotting

absolute <- DF2B %>%
  ggplot(., aes(x = as.factor(Group.1), y = x, fill = Group.3)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~Group.2, labeller = as_labeller(Depthlabs), scales = "free_y", ncol = 1)+
  scale_fill_manual(name="Fractions:",
                    values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                    labels=c("POM", "DOC", "S+A","S+C", "rSOC"))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  ylab("Amount of C [g/kg] in fractions")+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle(label = "Absolute C content")

# Creating a common legend for both plots 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(absolute)
mm <- theme(plot.margin=unit(rep(1.5,4), "line"),
            legend.position = "none")

# Arranging the two subplots to one plot
p2 <- grid.arrange(arrangeGrob(absolute + mm,
                               relative + mm,
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))

# Proportional N content of the fractions, normalised on 100%

dfAn <- select(results_HS, c(warming, Depth, POM_N_relative, SA_N_relative, SC_N_relative, rSN_relative)) # selects the relevant columns from the dataset

DF1An <- melt(dfAn, id.var=c("warming", "Depth")) # converts the subset from wide to long format
DF2An <- aggregate(x=DF1An$value,
                  by=list(DF1An$warming, DF1An$Depth, DF1An$variable),
                  FUN = mean, na.rm=TRUE) # aggregates the data for plotting

relativeN <- DF2An %>%
  ggplot(., aes(x = as.factor(Group.1), y = x, fill = Group.3)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~Group.2, labeller = as_labeller(Depthlabs), ncol = 1)+
  scale_fill_manual(name="Fractions:",
                    values=c("#000000", "#56B4E9", "#009E73", "#F0E442"),
                    labels=c("POM", "S+A","S+C", "rSN"))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  ylab("Proportions of N-Fractions")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle(label = "Relative N content")


# Absolute N content of the fractions und mg/kg
dfBn <- select(results_HS, c(warming, Depth, POM_N_gkg, SA_N_gkg, SC_N_gkg, rSN_gkg))

DF1Bn <- melt(dfBn, id.var=c("warming", "Depth"))
DF2Bn <- aggregate(x=DF1Bn$value,
                  by=list(DF1Bn$warming, DF1Bn$Depth, DF1Bn$variable),
                  FUN = mean, na.rm=TRUE)

absoluteN <- DF2Bn %>%
  ggplot(., aes(x = as.factor(Group.1), y = x, fill = Group.3)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~Group.2, labeller = as_labeller(Depthlabs), scales = "free_y", ncol = 1)+
  scale_fill_manual(name="Fractions:",
                    values=c("#000000", "#56B4E9", "#009E73", "#F0E442"),
                    labels=c("POM", "S+A","S+C", "rSN"))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  ylab("Amount of N [g/kg] in fractions")+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle(label = "Absolute N content")

mylegendN<-g_legend(absoluteN)
mm <- theme(plot.margin=unit(rep(1.5,4), "line"),
            legend.position = "none")

p3 <- grid.arrange(arrangeGrob(absoluteN + mm,
                               relativeN + mm,
                               nrow=1),
                   mylegendN, nrow=2,heights=c(10, 1))

##### Overview plot Fractions/Topsoil/Subsoil--------------
# C fractions
fracC <- select(plots_HS, c(plot, Whole_Profile_POM_C, Whole_Profile_DOC,Whole_Profile_SA_C, Whole_Profile_SC_rSOC)) # selects the relevant columns from the dataset
fracC1 <- melt(fracC, id.vars = c("plot")) # converts the subset from wide to long format
fracC2 <- aggregate(x=fracC1$value,
                    by=list(fracC1$plot, fracC1$variable),
                    FUN = mean, na.rm=TRUE) # aggregates the data for plotting

C_frac <- fracC2 %>%
  ggplot(., aes(x = rev(as.factor(Group.1)), y = x, fill = Group.2)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(name="Fractions:",
                    values=c("#000000", "#E69F00", "#56B4E9", "#009E73"),
                    labels=c("POM", "DOC", "S+A","S+C"))+
  ylab(bquote('C stocks in fractions [ '*'Mg '*ha^-1*']'))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  theme_bw()+
  theme(legend.position = "none", legend.title = element_blank())+
  annotate("text", x = 4, y = 8, label = "A")

# N fractions 
fracN <- select(plots_HS, c(plot, Whole_Profile_POM_N, Whole_Profile_SA_N, Whole_Profile_SC_rSN)) # selects the relevant columns from the dataset
fracN1 <- melt(fracN, id.vars = c("plot")) # converts the subset from wide to long format
fracN2 <- aggregate(x=fracN1$value,
                    by=list(fracN1$plot, fracN1$variable),
                    FUN = mean, na.rm=TRUE) # aggregates the data for plotting

N_frac <- fracN2 %>%
  ggplot(., aes(x = rev(as.factor(Group.1)), y = x, fill = Group.2)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(name="Fractions:",
                    values=c("#000000",  "#56B4E9", "#009E73"),
                    labels=c("POM", "S+A","S+C"))+
  ylab(bquote('C stocks in fractions [ '*'Mg '*ha^-1*']'))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  theme_bw()+
  theme(legend.position = "none", legend.title = element_blank())+
  annotate("text", x = 4, y = 0.5, label = "B")

fractions <- ggarrange(C_frac, N_frac, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom") # arranges the two subplots to one plot

##### Topsoil / Subsoil
# C in topsoil and subsoil
topsub <- select(plots_HS, c(plot, Topsoil_C, Subsoil_C))  # selects the relevant columns from the dataset
topsub1 <- melt(topsub, id.vars = "plot") # converts the subset from wide to long format
topsub2 <- aggregate(x=topsub1$value,
                     by=list(topsub1$plot, topsub1$variable),
                     FUN=mean) # aggregates the data for plotting

topsub_C <- topsub2 %>%
  ggplot(., aes(x = as.factor(Group.1), y=rev(x), fill = Group.2))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name="Share of topsoil and subsoil in total C stock",
                    values=c("darkgreen","grey60"),
                    labels=c("Topsoil", "Subsoil"))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  ylab(bquote('C stocks in soil layers [ '*'Mg '*ha^-1*']'))+
  theme_bw()+
  theme(legend.position = "none", legend.title = element_blank())+
  annotate("text", x = 4, y = 8, label = "C")

# N in topsoil and subsoil
topsubN <- select(plots_HS, c(plot, Topsoil_N, Subsoil_N)) # selects the relevant columns from the dataset
topsub1N <- melt(topsubN, id.vars = "plot") # converts the subset from wide to long format
topsub2N <- aggregate(x=topsub1N$value,
                     by=list(topsub1$plot, topsub1$variable),
                     FUN=mean)  # aggregates the data for plotting

topsub_N <- topsub2N %>%
  ggplot(., aes(x = as.factor(Group.1), y=rev(x), fill = Group.2))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name="Share of topsoil and subsoil in total N stock",
                    values=c("darkgreen","grey60"),
                    labels=c("Topsoil", "Subsoil"))+
  scale_x_discrete(name = "Warming intensity",
                   labels = c("Reference", "+0.5 °C", "+1.8 °C", "+7.7 °C"))+
  ylab(bquote('N stocks in soil layers [ '*'Mg '*ha^-1*']'))+
  theme_bw()+
  theme(legend.position = "none", legend.title = element_blank())+
  annotate("text", x = 4, y = 0.5, label = "D")

Topsoil_Subsoil <- ggarrange(topsub_C, topsub_N, ncol =2, nrow = 1, common.legend = TRUE, legend="bottom") # arranges the two subplots to one plot

p4 <-ggarrange(fractions, Topsoil_Subsoil, ncol=1, nrow=2) # creates one plot with four subplots 

##### Litter decomposition with teabags-------

teabags_HS <- teabags_HS %>%
  mutate(., loss = weight_start-weight_end) %>%
  mutate(., loss_perc = (1-(weight_end/weight_start))*100)%>%
  mutate(., rest = (1-loss_perc)) %>% # adds columns with the percentage of loss and origin 
  mutate(., temperature = case_when(
    plot == "W1" & depth == "10" ~ 10,
    plot == "W1" & depth == "50" ~ 12.4,
    plot == "W2" & depth == "10" ~ 4.3,
    plot == "W2" & depth == "50" ~ 6.3,
    plot == "W3" & depth == "10" ~ 3.6,
    plot == "W3" & depth == "50" ~ 4.4,
    plot == "W4" & depth == "10" ~ 3.1,
    plot == "W4" & depth == "50" ~ 3.8
  ))

teabags_mean <- melt(teabags_HS, id.var=c("temperature", "depth", "sample", "rep")) %>% # prepares the dataframe for plotting
  filter(., .$variable=="loss_perc") 

Depthlabs <- c('10' = "10 cm",
               '50' = "50 cm") # labels for facetting


lm_eqn = function(teabags_mean){
  m = lm(log(as.numeric(value)) ~ temperature, teabags_mean);
  eq <- substitute(italic(y) == a + e^(b+t) %.% ","~~italic(r)^2~"="~r2, 
                   list(a = round(as.numeric(coef(m)[1]),digits = 3), 
                        b = round(as.numeric(coef(m)[2]), digits = 4), 
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));                 
} # a function that calculates and extracts the parameters of the linear model. Will be used to display the equation in a plot

eq <- ddply(teabags_mean,.(depth),lm_eqn) # applies the created function to the teabag data

teabags_mean %>%
  filter(.$variable == "loss_perc") %>%
  na.omit(.)%>%
  ggplot(., aes(x = temperature, y = as.numeric(value))) + 
  geom_point(size = 6)+
  geom_smooth(method = "lm", formula = y ~ log(x))+
  scale_y_continuous(name = "Decomposition [%]")+
  scale_x_continuous(name = "Mean annual soil temperature [°C]")+
  facet_wrap(~depth,labeller = as_labeller(Depthlabs), scales = "free")+
  geom_text(data=eq,aes(x = 7, y = 50,label=V1), parse = TRUE, inherit.aes=FALSE, size = 4)+
  theme_bw()+
  theme(text = element_text(size = 20))


#### Calculation of Q10 values
depth10 <- subset(teabags_mean, depth == 10)
exponential.model10 <- lm(log(as.numeric(value))~ temperature, data = depth10)
summary(exponential.model10)
q10_10 <- exp(10*(summary(exponential.model10)$coefficients[2,1]))



depth50 <- subset(teabags_mean, depth == 50)
exponential.model50 <- lm(log(as.numeric(value))~ temperature, data = depth50)
summary(exponential.model50)
q10_50 <- exp(10*(summary(exponential.model50)$coefficients[2,1]))
##### Statistics--------------

plots_HS <- plots_HS %>%
  mutate(., temperature_topsoil = case_when(
    plot == "W1" ~ 10,
    plot == "W2" ~ 4.3,
    plot == "W3" ~ 3.6,
    plot == "W4" ~ 3.1
  )) %>%
  mutate(., temperature_subsoil = case_when(
    plot == "W1" ~ 12.4,
    plot == "W2" ~ 6.3,
    plot == "W3" ~ 4.4,
    plot == "W4" ~ 3.8
  )) %>%
  mutate(., temperature_mean = case_when(
    plot == "W1" ~ 11.2,
    plot == "W2" ~ 5.3,
    plot == "W3" ~ 4,
    plot == "W4" ~ 3.5))

cortable <- select(plots_HS, c(temperature_topsoil, temperature_subsoil, temperature_mean, Whole_Profile_C, Topsoil_C, Subsoil_C)) %>%
  na.omit() %>%
  cor()

# Mean C and N in Topsoil and Subsoil, including standard error 

mean_stocks <- plots_HS %>% 
  group_by(plot) %>%
  summarise(mean_C_topsoil = mean(Topsoil_C, na.rm=TRUE),  # calculates the mean of each plot
            sd_C_topsoil = sd(Topsoil_C, na.rm =TRUE), # calculates the standard deviation of each plot
            n_plots = n(),  # calculates the sample size of each plot plot
            SE_C_topsoil = sd(Topsoil_C)/sqrt(n()), # calculates the standard deviation of each plot
            mean_C_subsoil = mean(Subsoil_C, na.rm=TRUE), 
            sd_C_subsoil = sd(Subsoil_C, na.rm =TRUE),
            SE_C_subsoil = sd(Subsoil_C)/sqrt(n()),
            mean_N_topsoil = mean(Topsoil_N, na.rm=TRUE),  
            sd_N_topsoil = sd(Topsoil_N, na.rm =TRUE), 
            SE_N_topsoil = sd(Topsoil_N)/sqrt(n()),
            mean_N_subsoil = mean(Subsoil_N, na.rm=TRUE),  
            sd_N_subsoil = sd(Subsoil_N, na.rm =TRUE), 
            SE_N_subsoil = sd(Subsoil_N)/sqrt(n()))

# Anosim
# C fracions
d <- 70 # for depth-wise ANOSIM. Chose the desired depth to investigate
veg <- select(results_HS, c(warming, Depth, POM_C_relative, DOC_relative, SA_C_relative, SC_C_relative, rSOC_relative)) %>%
  filter(., Depth == d) %>%
  filter(., warming == 7.7 | warming == 0) %>%
  select(., c(3,4,5,6,7)) %>%
  as.matrix(.) # brings the data in a matrix that can be used with the anosim function
grp <-select(results_HS, c(warming, Depth, POM_C_relative, DOC_relative, SA_C_relative, SC_C_relative, rSOC_relative)) %>% 
  filter(., Depth == d) %>%
  filter(., warming == 7.7 | warming == 0) # Needed for the grouping parameter in anosim. The second filter can be set to any desired combination of groups or can be ommitted if you wish to compare all groups
 
veggi <- anosim(x = veg, grouping = grp$warming, permutations = 999)
plot(veggi)


# N fractions
d <- 70 # for depth-wise ANOSIM. Chose the desired depth to investigate
veg <- select(results_HS, c(warming, Depth, POM_N_relative, SA_N_relative, SC_N_relative, rSN_relative)) %>%
  filter(., Depth == d) %>%
  filter(., warming == 7.7 | warming == 0) %>%
  select(., c(3,4,5,6)) %>%
  as.matrix(.)
grp <-select(results_HS, c(warming, Depth, POM_N_relative, SA_N_relative, SC_N_relative, rSN_relative)) %>% 
  filter(., Depth == d) %>%
  filter(., warming == 7.7 | warming == 0)

veggi <- anosim(x = veg, grouping = grp$warming, permutations = 999)
plot(veggi)

# linear model of fraction stocks, logarithmic data

coefLog <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefLog[1,1] <- sub$Depth[1]
coefLog[1,2] <- "SOC stock"
coefLog[1,3] <- summary(lm(log(SOC_corr_Topsoil) ~ warming, data=sub))$coefficients[1,1]
coefLog[1,4] <- summary(lm(log(SOC_corr_Topsoil) ~ warming, data=sub))$coefficients[2,1]
coefLog[1,5] <- summary(lm(log(SOC_corr_Topsoil) ~ warming, data=sub))$r.squared
coefLog[1,6] <- summary(lm(log(SOC_corr_Topsoil) ~ warming, data=sub))$coefficients[2,4]

coefLog[2,1] <- sub$Depth[1]
coefLog[2,2] <- "POM"
coefLog[2,3] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[2,4] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[2,5] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$r.squared
coefLog[2,6] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[3,1] <- sub$Depth[1]
coefLog[3,2] <- "SA_C_stock"
coefLog[3,3] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[3,4] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[3,5] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$r.squared
coefLog[3,6] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[4,1] <- sub$Depth[1]
coefLog[4,2] <- "SC_C_stock"
coefLog[4,3] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[4,4] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[4,5] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[4,6] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[5,1] <- sub$Depth[1]
coefLog[5,2] <- "rSOC_stock"
coefLog[5,3] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[5,4] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[5,5] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[5,6] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[6,1] <- sub$Depth[1]
coefLog[6,2] <- "DOC_C_stock"
coefLog[6,3] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[6,4] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[6,5] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[6,6] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefLog[7,1] <- sub$Depth[1]
coefLog[7,2] <- "SOC_corr_Topsoil"
coefLog[7,3] <- summary(lm(log(SOC_corr_Topsoil) ~ warming  , data=sub))$coefficients[1,1]
coefLog[7,4] <- summary(lm(log(SOC_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,1]
coefLog[7,5] <- summary(lm(log(SOC_corr_Topsoil) ~ warming  , data=sub))$r.squared
coefLog[7,6] <- summary(lm(log(SOC_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,4]

coefLog[8,1] <- sub$Depth[1]
coefLog[8,2] <- "POM_C_stock"
coefLog[8,3] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[8,4] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[8,5] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$r.squared
coefLog[8,6] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[9,1] <- sub$Depth[1]
coefLog[9,2] <- "SA_C_stock"
coefLog[9,3] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[9,4] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[9,5] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$r.squared
coefLog[9,6] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[10,1] <- sub$Depth[1]
coefLog[10,2] <- "SC_C_stock"
coefLog[10,3] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[10,4] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[10,5] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[10,6] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[11,1] <- sub$Depth[1]
coefLog[11,2] <- "rSOC_stock"
coefLog[11,3] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[11,4] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[11,5] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[11,6] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[12,1] <- sub$Depth[1]
coefLog[12,2] <- "DOC_C_stock"
coefLog[12,3] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[12,4] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[12,5] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[12,6] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefLog[13,1] <- sub$Depth[1]
coefLog[13,2] <- "SOC_corr"
coefLog[13,3] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLog[13,4] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLog[13,5] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$r.squared
coefLog[13,6] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLog[14,1] <- sub$Depth[1]
coefLog[14,2] <- "POM_C_stock"
coefLog[14,3] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[14,4] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[14,5] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$r.squared
coefLog[14,6] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[15,1] <- sub$Depth[1]
coefLog[15,2] <- "SA_C_stock"
coefLog[15,3] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[15,4] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[15,5] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$r.squared
coefLog[15,6] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[16,1] <- sub$Depth[1]
coefLog[16,2] <- "SC_C_stock"
coefLog[16,3] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[16,4] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[16,5] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[16,6] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[17,1] <- sub$Depth[1]
coefLog[17,2] <- "rSOC_stock"
coefLog[17,3] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[17,4] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[17,5] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[17,6] <- summary(lm(log(rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[18,1] <- sub$Depth[1]
coefLog[18,2] <- "DOC_C_stock"
coefLog[18,3] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[18,4] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[18,5] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[18,6] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefLog[19,1] <- sub$Depth[1]
coefLog[19,2] <- "SOC_corr"
coefLog[19,3] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLog[19,4] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLog[19,5] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$r.squared
coefLog[19,6] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLog[20,1] <- sub$Depth[1]
coefLog[20,2] <- "POM_C_stock"
coefLog[20,3] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[20,4] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[20,5] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$r.squared
coefLog[20,6] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[21,1] <- sub$Depth[1]
coefLog[21,2] <- "SA_C_stock"
coefLog[21,3] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[21,4] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[21,5] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$r.squared
coefLog[21,6] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[22,1] <- sub$Depth[1]
coefLog[22,2] <- "SC_C_stock"
coefLog[22,3] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[22,4] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[22,5] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[22,6] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[23,1] <- sub$Depth[1]
coefLog[23,2] <- "DOC_C_stock"
coefLog[23,3] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[23,4] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[23,5] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[23,6] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefLog[24,1] <- sub$Depth[1]
coefLog[24,2] <- "SOC_corr"
coefLog[24,3] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLog[24,4] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLog[24,5] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$r.squared
coefLog[24,6] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLog[25,1] <- sub$Depth[1]
coefLog[25,2] <- "POM_C_stock"
coefLog[25,3] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[25,4] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[25,5] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$r.squared
coefLog[25,6] <- summary(lm(log(POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[26,1] <- sub$Depth[1]
coefLog[26,2] <- "SA_C_stock"
coefLog[26,3] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[26,4] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[26,5] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$r.squared
coefLog[26,6] <- summary(lm(log(SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLog[27,1] <- sub$Depth[1]
coefLog[27,2] <- "SC_C_stock"
coefLog[27,3] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[27,4] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[27,5] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[27,6] <- summary(lm(log(SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]


coefLog[28,1] <- sub$Depth[1]
coefLog[28,2] <- "DOC_C_stock"
coefLog[28,3] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLog[28,4] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLog[28,5] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLog[28,6] <- summary(lm(log(DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefLog) <- c("Depth",
                    "Fraction",
                    "Intercept",
                    "Slope",
                    "R2",
                    "p")


# linear model of fraction stocks, untransformed data
coefLin <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefLin[1,1] <- sub$Depth[1]
coefLin[1,2] <- "SOC stock"
coefLin[1,3] <- summary(lm((SOC_corr_Topsoil) ~ warming, data=sub))$coefficients[1,1]
coefLin[1,4] <- summary(lm((SOC_corr_Topsoil) ~ warming, data=sub))$coefficients[2,1]
coefLin[1,5] <- summary(lm((SOC_corr_Topsoil) ~ warming, data=sub))$r.squared
coefLin[1,6] <- summary(lm((SOC_corr_Topsoil) ~ warming, data=sub))$coefficients[2,4]

coefLin[2,1] <- sub$Depth[1]
coefLin[2,2] <- "POM"
coefLin[2,3] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[2,4] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[2,5] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$r.squared
coefLin[2,6] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[3,1] <- sub$Depth[1]
coefLin[3,2] <- "SA_C_stock"
coefLin[3,3] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[3,4] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[3,5] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$r.squared
coefLin[3,6] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[4,1] <- sub$Depth[1]
coefLin[4,2] <- "SC_C_stock"
coefLin[4,3] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[4,4] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[4,5] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[4,6] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[5,1] <- sub$Depth[1]
coefLin[5,2] <- "rSOC_stock"
coefLin[5,3] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[5,4] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[5,5] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[5,6] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[6,1] <- sub$Depth[1]
coefLin[6,2] <- "DOC_C_stock"
coefLin[6,3] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[6,4] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[6,5] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[6,6] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefLin[7,1] <- sub$Depth[1]
coefLin[7,2] <- "SOC_corr_Topsoil"
coefLin[7,3] <- summary(lm((SOC_corr_Topsoil) ~ warming  , data=sub))$coefficients[1,1]
coefLin[7,4] <- summary(lm((SOC_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,1]
coefLin[7,5] <- summary(lm((SOC_corr_Topsoil) ~ warming  , data=sub))$r.squared
coefLin[7,6] <- summary(lm((SOC_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,4]

coefLin[8,1] <- sub$Depth[1]
coefLin[8,2] <- "POM_C_stock"
coefLin[8,3] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[8,4] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[8,5] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$r.squared
coefLin[8,6] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[9,1] <- sub$Depth[1]
coefLin[9,2] <- "SA_C_stock"
coefLin[9,3] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[9,4] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[9,5] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$r.squared
coefLin[9,6] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[10,1] <- sub$Depth[1]
coefLin[10,2] <- "SC_C_stock"
coefLin[10,3] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[10,4] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[10,5] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[10,6] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[11,1] <- sub$Depth[1]
coefLin[11,2] <- "rSOC_stock"
coefLin[11,3] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[11,4] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[11,5] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[11,6] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[12,1] <- sub$Depth[1]
coefLin[12,2] <- "DOC_C_stock"
coefLin[12,3] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[12,4] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[12,5] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[12,6] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefLin[13,1] <- sub$Depth[1]
coefLin[13,2] <- "SOC_corr"
coefLin[13,3] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLin[13,4] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLin[13,5] <- summary(lm((SOC_corr) ~ warming  , data=sub))$r.squared
coefLin[13,6] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLin[14,1] <- sub$Depth[1]
coefLin[14,2] <- "POM_C_stock"
coefLin[14,3] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[14,4] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[14,5] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$r.squared
coefLin[14,6] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[15,1] <- sub$Depth[1]
coefLin[15,2] <- "SA_C_stock"
coefLin[15,3] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[15,4] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[15,5] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$r.squared
coefLin[15,6] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[16,1] <- sub$Depth[1]
coefLin[16,2] <- "SC_C_stock"
coefLin[16,3] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[16,4] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[16,5] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[16,6] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[17,1] <- sub$Depth[1]
coefLin[17,2] <- "rSOC_stock"
coefLin[17,3] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[17,4] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[17,5] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[17,6] <- summary(lm((rSOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[18,1] <- sub$Depth[1]
coefLin[18,2] <- "DOC_C_stock"
coefLin[18,3] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[18,4] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[18,5] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[18,6] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefLin[19,1] <- sub$Depth[1]
coefLin[19,2] <- "SOC_corr"
coefLin[19,3] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLin[19,4] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLin[19,5] <- summary(lm((SOC_corr) ~ warming  , data=sub))$r.squared
coefLin[19,6] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLin[20,1] <- sub$Depth[1]
coefLin[20,2] <- "POM_C_stock"
coefLin[20,3] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[20,4] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[20,5] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$r.squared
coefLin[20,6] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[21,1] <- sub$Depth[1]
coefLin[21,2] <- "SA_C_stock"
coefLin[21,3] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[21,4] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[21,5] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$r.squared
coefLin[21,6] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[22,1] <- sub$Depth[1]
coefLin[22,2] <- "SC_C_stock"
coefLin[22,3] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[22,4] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[22,5] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[22,6] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[23,1] <- sub$Depth[1]
coefLin[23,2] <- "DOC_C_stock"
coefLin[23,3] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[23,4] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[23,5] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[23,6] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefLin[24,1] <- sub$Depth[1]
coefLin[24,2] <- "SOC_corr"
coefLin[24,3] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLin[24,4] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLin[24,5] <- summary(lm((SOC_corr) ~ warming  , data=sub))$r.squared
coefLin[24,6] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLin[25,1] <- sub$Depth[1]
coefLin[25,2] <- "POM_C_stock"
coefLin[25,3] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[25,4] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[25,5] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$r.squared
coefLin[25,6] <- summary(lm((POM_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[26,1] <- sub$Depth[1]
coefLin[26,2] <- "SA_C_stock"
coefLin[26,3] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[26,4] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[26,5] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$r.squared
coefLin[26,6] <- summary(lm((SA_C_stock) ~ warming  , data=sub))$coefficients[2,4]

coefLin[27,1] <- sub$Depth[1]
coefLin[27,2] <- "SC_C_stock"
coefLin[27,3] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[27,4] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[27,5] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[27,6] <- summary(lm((SC_C_stock) ~ warming  , data=sub))$coefficients[2,4]


coefLin[28,1] <- sub$Depth[1]
coefLin[28,2] <- "DOC_C_stock"
coefLin[28,3] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[1,1]
coefLin[28,4] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,1]
coefLin[28,5] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$r.squared
coefLin[28,6] <- summary(lm((DOC_C_stock) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefLin) <- c("Depth",
                       "Fraction",
                       "Intercept",
                       "Slope",
                       "R2",
                       "p")
# linear model of concentrations of the fractions, logarithmic data

coefLogConc <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefLogConc[1,1] <- sub$Depth[1]
coefLogConc[1,2] <- "SOC gkg"
coefLogConc[1,3] <- summary(lm(log(TOC_gkg) ~ warming, data=sub))$coefficients[1,1]
coefLogConc[1,4] <- summary(lm(log(TOC_gkg) ~ warming, data=sub))$coefficients[2,1]
coefLogConc[1,5] <- summary(lm(log(TOC_gkg) ~ warming, data=sub))$r.squared
coefLogConc[1,6] <- summary(lm(log(TOC_gkg) ~ warming, data=sub))$coefficients[2,4]

coefLogConc[2,1] <- sub$Depth[1]
coefLogConc[2,2] <- "POM"
coefLogConc[2,3] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[2,4] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[2,5] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[2,6] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[3,1] <- sub$Depth[1]
coefLogConc[3,2] <- "SA_C_gkg"
coefLogConc[3,3] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[3,4] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[3,5] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[3,6] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[4,1] <- sub$Depth[1]
coefLogConc[4,2] <- "SC_C_gkg"
coefLogConc[4,3] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[4,4] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[4,5] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[4,6] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[5,1] <- sub$Depth[1]
coefLogConc[5,2] <- "rSOC_gkg"
coefLogConc[5,3] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[5,4] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[5,5] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[5,6] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[6,1] <- sub$Depth[1]
coefLogConc[6,2] <- "DOC_gkg"
coefLogConc[6,3] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[6,4] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[6,5] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[6,6] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefLogConc[7,1] <- sub$Depth[1]
coefLogConc[7,2] <- "TOC_gkg"
coefLogConc[7,3] <- summary(lm(log(TOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[7,4] <- summary(lm(log(TOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[7,5] <- summary(lm(log(TOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[7,6] <- summary(lm(log(TOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[8,1] <- sub$Depth[1]
coefLogConc[8,2] <- "POM_C_gkg"
coefLogConc[8,3] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[8,4] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[8,5] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[8,6] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[9,1] <- sub$Depth[1]
coefLogConc[9,2] <- "SA_C_gkg"
coefLogConc[9,3] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[9,4] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[9,5] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[9,6] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[10,1] <- sub$Depth[1]
coefLogConc[10,2] <- "SC_C_gkg"
coefLogConc[10,3] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[10,4] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[10,5] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[10,6] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[11,1] <- sub$Depth[1]
coefLogConc[11,2] <- "rSOC_gkg"
coefLogConc[11,3] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[11,4] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[11,5] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[11,6] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[12,1] <- sub$Depth[1]
coefLogConc[12,2] <- "DOC_gkg"
coefLogConc[12,3] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[12,4] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[12,5] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[12,6] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefLogConc[13,1] <- sub$Depth[1]
coefLogConc[13,2] <- "SOC_corr"
coefLogConc[13,3] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[13,4] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[13,5] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$r.squared
coefLogConc[13,6] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[14,1] <- sub$Depth[1]
coefLogConc[14,2] <- "POM_C_gkg"
coefLogConc[14,3] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[14,4] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[14,5] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[14,6] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[15,1] <- sub$Depth[1]
coefLogConc[15,2] <- "SA_C_gkg"
coefLogConc[15,3] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[15,4] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[15,5] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[15,6] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[16,1] <- sub$Depth[1]
coefLogConc[16,2] <- "SC_C_gkg"
coefLogConc[16,3] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[16,4] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[16,5] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[16,6] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[17,1] <- sub$Depth[1]
coefLogConc[17,2] <- "rSOC_gkg"
coefLogConc[17,3] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[17,4] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[17,5] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[17,6] <- summary(lm(log(rSOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[18,1] <- sub$Depth[1]
coefLogConc[18,2] <- "DOC_gkg"
coefLogConc[18,3] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[18,4] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[18,5] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[18,6] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefLogConc[19,1] <- sub$Depth[1]
coefLogConc[19,2] <- "SOC_corr"
coefLogConc[19,3] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[19,4] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[19,5] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$r.squared
coefLogConc[19,6] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[20,1] <- sub$Depth[1]
coefLogConc[20,2] <- "POM_C_gkg"
coefLogConc[20,3] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[20,4] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[20,5] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[20,6] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[21,1] <- sub$Depth[1]
coefLogConc[21,2] <- "SA_C_gkg"
coefLogConc[21,3] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[21,4] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[21,5] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[21,6] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[22,1] <- sub$Depth[1]
coefLogConc[22,2] <- "SC_C_gkg"
coefLogConc[22,3] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[22,4] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[22,5] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[22,6] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[23,1] <- sub$Depth[1]
coefLogConc[23,2] <- "DOC_gkg"
coefLogConc[23,3] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[23,4] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[23,5] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[23,6] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefLogConc[24,1] <- sub$Depth[1]
coefLogConc[24,2] <- "SOC_corr"
coefLogConc[24,3] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[24,4] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[24,5] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$r.squared
coefLogConc[24,6] <- summary(lm(log(SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[25,1] <- sub$Depth[1]
coefLogConc[25,2] <- "POM_C_gkg"
coefLogConc[25,3] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[25,4] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[25,5] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[25,6] <- summary(lm(log(POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[26,1] <- sub$Depth[1]
coefLogConc[26,2] <- "SA_C_gkg"
coefLogConc[26,3] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[26,4] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[26,5] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[26,6] <- summary(lm(log(SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLogConc[27,1] <- sub$Depth[1]
coefLogConc[27,2] <- "SC_C_gkg"
coefLogConc[27,3] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[27,4] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[27,5] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[27,6] <- summary(lm(log(SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]


coefLogConc[28,1] <- sub$Depth[1]
coefLogConc[28,2] <- "DOC_gkg"
coefLogConc[28,3] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLogConc[28,4] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLogConc[28,5] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$r.squared
coefLogConc[28,6] <- summary(lm(log(DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefLogConc) <- c("Depth",
                       "Fraction",
                       "Intercept",
                       "Slope",
                       "R2",
                       "p")


# linear model of concentrations of the fractions, untransformed data
coefLinConc <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefLinConc[1,1] <- sub$Depth[1]
coefLinConc[1,2] <- "SOC gkg"
coefLinConc[1,3] <- summary(lm((TOC_gkg) ~ warming, data=sub))$coefficients[1,1]
coefLinConc[1,4] <- summary(lm((TOC_gkg) ~ warming, data=sub))$coefficients[2,1]
coefLinConc[1,5] <- summary(lm((TOC_gkg) ~ warming, data=sub))$r.squared
coefLinConc[1,6] <- summary(lm((TOC_gkg) ~ warming, data=sub))$coefficients[2,4]

coefLinConc[2,1] <- sub$Depth[1]
coefLinConc[2,2] <- "POM"
coefLinConc[2,3] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[2,4] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[2,5] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[2,6] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[3,1] <- sub$Depth[1]
coefLinConc[3,2] <- "SA_C_gkg"
coefLinConc[3,3] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[3,4] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[3,5] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[3,6] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[4,1] <- sub$Depth[1]
coefLinConc[4,2] <- "SC_C_gkg"
coefLinConc[4,3] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[4,4] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[4,5] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[4,6] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[5,1] <- sub$Depth[1]
coefLinConc[5,2] <- "rSOC_gkg"
coefLinConc[5,3] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[5,4] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[5,5] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[5,6] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[6,1] <- sub$Depth[1]
coefLinConc[6,2] <- "DOC_gkg"
coefLinConc[6,3] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[6,4] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[6,5] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[6,6] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefLinConc[7,1] <- sub$Depth[1]
coefLinConc[7,2] <- "TOC_gkg"
coefLinConc[7,3] <- summary(lm((TOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[7,4] <- summary(lm((TOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[7,5] <- summary(lm((TOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[7,6] <- summary(lm((TOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[8,1] <- sub$Depth[1]
coefLinConc[8,2] <- "POM_C_gkg"
coefLinConc[8,3] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[8,4] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[8,5] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[8,6] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[9,1] <- sub$Depth[1]
coefLinConc[9,2] <- "SA_C_gkg"
coefLinConc[9,3] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[9,4] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[9,5] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[9,6] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[10,1] <- sub$Depth[1]
coefLinConc[10,2] <- "SC_C_gkg"
coefLinConc[10,3] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[10,4] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[10,5] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[10,6] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[11,1] <- sub$Depth[1]
coefLinConc[11,2] <- "rSOC_gkg"
coefLinConc[11,3] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[11,4] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[11,5] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[11,6] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[12,1] <- sub$Depth[1]
coefLinConc[12,2] <- "DOC_gkg"
coefLinConc[12,3] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[12,4] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[12,5] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[12,6] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefLinConc[13,1] <- sub$Depth[1]
coefLinConc[13,2] <- "SOC_corr"
coefLinConc[13,3] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[13,4] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[13,5] <- summary(lm((SOC_corr) ~ warming  , data=sub))$r.squared
coefLinConc[13,6] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[14,1] <- sub$Depth[1]
coefLinConc[14,2] <- "POM_C_gkg"
coefLinConc[14,3] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[14,4] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[14,5] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[14,6] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[15,1] <- sub$Depth[1]
coefLinConc[15,2] <- "SA_C_gkg"
coefLinConc[15,3] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[15,4] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[15,5] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[15,6] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[16,1] <- sub$Depth[1]
coefLinConc[16,2] <- "SC_C_gkg"
coefLinConc[16,3] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[16,4] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[16,5] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[16,6] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[17,1] <- sub$Depth[1]
coefLinConc[17,2] <- "rSOC_gkg"
coefLinConc[17,3] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[17,4] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[17,5] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[17,6] <- summary(lm((rSOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[18,1] <- sub$Depth[1]
coefLinConc[18,2] <- "DOC_gkg"
coefLinConc[18,3] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[18,4] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[18,5] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[18,6] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefLinConc[19,1] <- sub$Depth[1]
coefLinConc[19,2] <- "SOC_corr"
coefLinConc[19,3] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[19,4] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[19,5] <- summary(lm((SOC_corr) ~ warming  , data=sub))$r.squared
coefLinConc[19,6] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[20,1] <- sub$Depth[1]
coefLinConc[20,2] <- "POM_C_gkg"
coefLinConc[20,3] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[20,4] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[20,5] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[20,6] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[21,1] <- sub$Depth[1]
coefLinConc[21,2] <- "SA_C_gkg"
coefLinConc[21,3] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[21,4] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[21,5] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[21,6] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[22,1] <- sub$Depth[1]
coefLinConc[22,2] <- "SC_C_gkg"
coefLinConc[22,3] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[22,4] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[22,5] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[22,6] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[23,1] <- sub$Depth[1]
coefLinConc[23,2] <- "DOC_gkg"
coefLinConc[23,3] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[23,4] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[23,5] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[23,6] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefLinConc[24,1] <- sub$Depth[1]
coefLinConc[24,2] <- "SOC_corr"
coefLinConc[24,3] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[24,4] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[24,5] <- summary(lm((SOC_corr) ~ warming  , data=sub))$r.squared
coefLinConc[24,6] <- summary(lm((SOC_corr) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[25,1] <- sub$Depth[1]
coefLinConc[25,2] <- "POM_C_gkg"
coefLinConc[25,3] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[25,4] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[25,5] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[25,6] <- summary(lm((POM_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[26,1] <- sub$Depth[1]
coefLinConc[26,2] <- "SA_C_gkg"
coefLinConc[26,3] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[26,4] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[26,5] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[26,6] <- summary(lm((SA_C_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefLinConc[27,1] <- sub$Depth[1]
coefLinConc[27,2] <- "SC_C_gkg"
coefLinConc[27,3] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[27,4] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[27,5] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[27,6] <- summary(lm((SC_C_gkg) ~ warming  , data=sub))$coefficients[2,4]


coefLinConc[28,1] <- sub$Depth[1]
coefLinConc[28,2] <- "DOC_gkg"
coefLinConc[28,3] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefLinConc[28,4] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefLinConc[28,5] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$r.squared
coefLinConc[28,6] <- summary(lm((DOC_gkg) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefLinConc) <- c("Depth",
                       "Fraction",
                       "Intercept",
                       "Slope",
                       "R2",
                       "p")

# write the coefficients to an excel file
model_output <- list("Linear stocks" = coefLin, "log_Transformed stocks" = coefLog)
write.xlsx(model_output, "D:/Tino/Auswertung/Model.xlsx")



# linear model of fraction stocks, logarithmic data

coefNLog <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefNLog[1,1] <- sub$Depth[1]
coefNLog[1,2] <- "N stock"
coefNLog[1,3] <- summary(lm(log(N_corr_Topsoil) ~ warming, data=sub))$coefficients[1,1]
coefNLog[1,4] <- summary(lm(log(N_corr_Topsoil) ~ warming, data=sub))$coefficients[2,1]
coefNLog[1,5] <- summary(lm(log(N_corr_Topsoil) ~ warming, data=sub))$r.squared
coefNLog[1,6] <- summary(lm(log(N_corr_Topsoil) ~ warming, data=sub))$coefficients[2,4]

coefNLog[2,1] <- sub$Depth[1]
coefNLog[2,2] <- "POM"
coefNLog[2,3] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[2,4] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[2,5] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[2,6] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[3,1] <- sub$Depth[1]
coefNLog[3,2] <- "SA_N_stock"
coefNLog[3,3] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[3,4] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[3,5] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[3,6] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[4,1] <- sub$Depth[1]
coefNLog[4,2] <- "SC_N_stock"
coefNLog[4,3] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[4,4] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[4,5] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[4,6] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[5,1] <- sub$Depth[1]
coefNLog[5,2] <- "rN_stock"
coefNLog[5,3] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[5,4] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[5,5] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$r.squared
coefNLog[5,6] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefNLog[7,1] <- sub$Depth[1]
coefNLog[7,2] <- "N_corr_Topsoil"
coefNLog[7,3] <- summary(lm(log(N_corr_Topsoil) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[7,4] <- summary(lm(log(N_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[7,5] <- summary(lm(log(N_corr_Topsoil) ~ warming  , data=sub))$r.squared
coefNLog[7,6] <- summary(lm(log(N_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[8,1] <- sub$Depth[1]
coefNLog[8,2] <- "POM_N_stock"
coefNLog[8,3] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[8,4] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[8,5] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[8,6] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[9,1] <- sub$Depth[1]
coefNLog[9,2] <- "SA_N_stock"
coefNLog[9,3] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[9,4] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[9,5] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[9,6] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[10,1] <- sub$Depth[1]
coefNLog[10,2] <- "SC_N_stock"
coefNLog[10,3] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[10,4] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[10,5] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[10,6] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[11,1] <- sub$Depth[1]
coefNLog[11,2] <- "rN_stock"
coefNLog[11,3] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[11,4] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[11,5] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$r.squared
coefNLog[11,6] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefNLog[13,1] <- sub$Depth[1]
coefNLog[13,2] <- "N_corr"
coefNLog[13,3] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[13,4] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[13,5] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$r.squared
coefNLog[13,6] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[14,1] <- sub$Depth[1]
coefNLog[14,2] <- "POM_N_stock"
coefNLog[14,3] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[14,4] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[14,5] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[14,6] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[15,1] <- sub$Depth[1]
coefNLog[15,2] <- "SA_N_stock"
coefNLog[15,3] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[15,4] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[15,5] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[15,6] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[16,1] <- sub$Depth[1]
coefNLog[16,2] <- "SC_N_stock"
coefNLog[16,3] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[16,4] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[16,5] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[16,6] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[17,1] <- sub$Depth[1]
coefNLog[17,2] <- "rN_stock"
coefNLog[17,3] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[17,4] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[17,5] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$r.squared
coefNLog[17,6] <- summary(lm(log(rSON_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefNLog[19,1] <- sub$Depth[1]
coefNLog[19,2] <- "N_corr"
coefNLog[19,3] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[19,4] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[19,5] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$r.squared
coefNLog[19,6] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[20,1] <- sub$Depth[1]
coefNLog[20,2] <- "POM_N_stock"
coefNLog[20,3] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[20,4] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[20,5] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[20,6] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[21,1] <- sub$Depth[1]
coefNLog[21,2] <- "SA_N_stock"
coefNLog[21,3] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[21,4] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[21,5] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[21,6] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[22,1] <- sub$Depth[1]
coefNLog[22,2] <- "SC_N_stock"
coefNLog[22,3] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[22,4] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[22,5] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[22,6] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefNLog[24,1] <- sub$Depth[1]
coefNLog[24,2] <- "N_corr"
coefNLog[24,3] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[24,4] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[24,5] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$r.squared
coefNLog[24,6] <- summary(lm(log(N_Corr) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[25,1] <- sub$Depth[1]
coefNLog[25,2] <- "POM_N_stock"
coefNLog[25,3] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[25,4] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[25,5] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[25,6] <- summary(lm(log(POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[26,1] <- sub$Depth[1]
coefNLog[26,2] <- "SA_N_stock"
coefNLog[26,3] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[26,4] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[26,5] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[26,6] <- summary(lm(log(SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLog[27,1] <- sub$Depth[1]
coefNLog[27,2] <- "SC_N_stock"
coefNLog[27,3] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLog[27,4] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLog[27,5] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLog[27,6] <- summary(lm(log(SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefNLog) <- c("Depth",
                       "Fraction",
                       "Intercept",
                       "Slope",
                       "R2",
                       "p")


# linear model of fraction stocks, untransformed data
coefNLin <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefNLin[1,1] <- sub$Depth[1]
coefNLin[1,2] <- "N stock"
coefNLin[1,3] <- summary(lm((N_corr_Topsoil) ~ warming, data=sub))$coefficients[1,1]
coefNLin[1,4] <- summary(lm((N_corr_Topsoil) ~ warming, data=sub))$coefficients[2,1]
coefNLin[1,5] <- summary(lm((N_corr_Topsoil) ~ warming, data=sub))$r.squared
coefNLin[1,6] <- summary(lm((N_corr_Topsoil) ~ warming, data=sub))$coefficients[2,4]

coefNLin[2,1] <- sub$Depth[1]
coefNLin[2,2] <- "POM"
coefNLin[2,3] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[2,4] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[2,5] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[2,6] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[3,1] <- sub$Depth[1]
coefNLin[3,2] <- "SA_N_stock"
coefNLin[3,3] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[3,4] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[3,5] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[3,6] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[4,1] <- sub$Depth[1]
coefNLin[4,2] <- "SC_N_stock"
coefNLin[4,3] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[4,4] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[4,5] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[4,6] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[5,1] <- sub$Depth[1]
coefNLin[5,2] <- "rN_stock"
coefNLin[5,3] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[5,4] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[5,5] <- summary(lm((rSON_stock) ~ warming  , data=sub))$r.squared
coefNLin[5,6] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefNLin[7,1] <- sub$Depth[1]
coefNLin[7,2] <- "N_corr_Topsoil"
coefNLin[7,3] <- summary(lm((N_corr_Topsoil) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[7,4] <- summary(lm((N_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[7,5] <- summary(lm((N_corr_Topsoil) ~ warming  , data=sub))$r.squared
coefNLin[7,6] <- summary(lm((N_corr_Topsoil) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[8,1] <- sub$Depth[1]
coefNLin[8,2] <- "POM_N_stock"
coefNLin[8,3] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[8,4] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[8,5] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[8,6] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[9,1] <- sub$Depth[1]
coefNLin[9,2] <- "SA_N_stock"
coefNLin[9,3] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[9,4] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[9,5] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[9,6] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[10,1] <- sub$Depth[1]
coefNLin[10,2] <- "SC_N_stock"
coefNLin[10,3] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[10,4] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[10,5] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[10,6] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[11,1] <- sub$Depth[1]
coefNLin[11,2] <- "rN_stock"
coefNLin[11,3] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[11,4] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[11,5] <- summary(lm((rSON_stock) ~ warming  , data=sub))$r.squared
coefNLin[11,6] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefNLin[13,1] <- sub$Depth[1]
coefNLin[13,2] <- "N_Corr"
coefNLin[13,3] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[13,4] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[13,5] <- summary(lm((N_Corr) ~ warming  , data=sub))$r.squared
coefNLin[13,6] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[14,1] <- sub$Depth[1]
coefNLin[14,2] <- "POM_N_stock"
coefNLin[14,3] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[14,4] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[14,5] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[14,6] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[15,1] <- sub$Depth[1]
coefNLin[15,2] <- "SA_N_stock"
coefNLin[15,3] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[15,4] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[15,5] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[15,6] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[16,1] <- sub$Depth[1]
coefNLin[16,2] <- "SC_N_stock"
coefNLin[16,3] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[16,4] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[16,5] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[16,6] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[17,1] <- sub$Depth[1]
coefNLin[17,2] <- "rN_stock"
coefNLin[17,3] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[17,4] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[17,5] <- summary(lm((rSON_stock) ~ warming  , data=sub))$r.squared
coefNLin[17,6] <- summary(lm((rSON_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefNLin[19,1] <- sub$Depth[1]
coefNLin[19,2] <- "N_Corr"
coefNLin[19,3] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[19,4] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[19,5] <- summary(lm((N_Corr) ~ warming  , data=sub))$r.squared
coefNLin[19,6] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[20,1] <- sub$Depth[1]
coefNLin[20,2] <- "POM_N_stock"
coefNLin[20,3] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[20,4] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[20,5] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[20,6] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[21,1] <- sub$Depth[1]
coefNLin[21,2] <- "SA_N_stock"
coefNLin[21,3] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[21,4] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[21,5] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[21,6] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[22,1] <- sub$Depth[1]
coefNLin[22,2] <- "SC_N_stock"
coefNLin[22,3] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[22,4] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[22,5] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[22,6] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefNLin[24,1] <- sub$Depth[1]
coefNLin[24,2] <- "N_Corr"
coefNLin[24,3] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[24,4] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[24,5] <- summary(lm((N_Corr) ~ warming  , data=sub))$r.squared
coefNLin[24,6] <- summary(lm((N_Corr) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[25,1] <- sub$Depth[1]
coefNLin[25,2] <- "POM_N_stock"
coefNLin[25,3] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[25,4] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[25,5] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[25,6] <- summary(lm((POM_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[26,1] <- sub$Depth[1]
coefNLin[26,2] <- "SA_N_stock"
coefNLin[26,3] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[26,4] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[26,5] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[26,6] <- summary(lm((SA_N_stock) ~ warming  , data=sub))$coefficients[2,4]

coefNLin[27,1] <- sub$Depth[1]
coefNLin[27,2] <- "SC_N_stock"
coefNLin[27,3] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[1,1]
coefNLin[27,4] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,1]
coefNLin[27,5] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$r.squared
coefNLin[27,6] <- summary(lm((SC_N_stock) ~ warming  , data=sub))$coefficients[2,4]


colnames(coefLin) <- c("Depth",
                       "Fraction",
                       "Intercept",
                       "Slope",
                       "R2",
                       "p")

coefNLogConc <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefNLogConc[1,1] <- sub$Depth[1]
coefNLogConc[1,2] <- "SOC gkg"
coefNLogConc[1,3] <- summary(lm(log(TN_gkg) ~ warming, data=sub))$coefficients[1,1]
coefNLogConc[1,4] <- summary(lm(log(TN_gkg) ~ warming, data=sub))$coefficients[2,1]
coefNLogConc[1,5] <- summary(lm(log(TN_gkg) ~ warming, data=sub))$r.squared
coefNLogConc[1,6] <- summary(lm(log(TN_gkg) ~ warming, data=sub))$coefficients[2,4]

coefNLogConc[2,1] <- sub$Depth[1]
coefNLogConc[2,2] <- "POM"
coefNLogConc[2,3] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[2,4] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[2,5] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[2,6] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[3,1] <- sub$Depth[1]
coefNLogConc[3,2] <- "SA_N_gkg"
coefNLogConc[3,3] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[3,4] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[3,5] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[3,6] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[4,1] <- sub$Depth[1]
coefNLogConc[4,2] <- "SC_N_gkg"
coefNLogConc[4,3] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[4,4] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[4,5] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[4,6] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[5,1] <- sub$Depth[1]
coefNLogConc[5,2] <- "rSN_gkg"
coefNLogConc[5,3] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[5,4] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[5,5] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[5,6] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefNLogConc[7,1] <- sub$Depth[1]
coefNLogConc[7,2] <- "TN_gkg"
coefNLogConc[7,3] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[7,4] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[7,5] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[7,6] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[8,1] <- sub$Depth[1]
coefNLogConc[8,2] <- "POM_N_gkg"
coefNLogConc[8,3] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[8,4] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[8,5] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[8,6] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[9,1] <- sub$Depth[1]
coefNLogConc[9,2] <- "SA_N_gkg"
coefNLogConc[9,3] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[9,4] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[9,5] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[9,6] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[10,1] <- sub$Depth[1]
coefNLogConc[10,2] <- "SC_N_gkg"
coefNLogConc[10,3] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[10,4] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[10,5] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[10,6] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[11,1] <- sub$Depth[1]
coefNLogConc[11,2] <- "rSN_gkg"
coefNLogConc[11,3] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[11,4] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[11,5] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[11,6] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefNLogConc[13,1] <- sub$Depth[1]
coefNLogConc[13,2] <- "TN_gkg"
coefNLogConc[13,3] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[13,4] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[13,5] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[13,6] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[14,1] <- sub$Depth[1]
coefNLogConc[14,2] <- "POM_N_gkg"
coefNLogConc[14,3] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[14,4] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[14,5] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[14,6] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[15,1] <- sub$Depth[1]
coefNLogConc[15,2] <- "SA_N_gkg"
coefNLogConc[15,3] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[15,4] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[15,5] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[15,6] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[16,1] <- sub$Depth[1]
coefNLogConc[16,2] <- "SC_N_gkg"
coefNLogConc[16,3] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[16,4] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[16,5] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[16,6] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[17,1] <- sub$Depth[1]
coefNLogConc[17,2] <- "rSN_gkg"
coefNLogConc[17,3] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[17,4] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[17,5] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[17,6] <- summary(lm(log(rSN_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefNLogConc[19,1] <- sub$Depth[1]
coefNLogConc[19,2] <- "TN_gkg"
coefNLogConc[19,3] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[19,4] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[19,5] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[19,6] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[20,1] <- sub$Depth[1]
coefNLogConc[20,2] <- "POM_N_gkg"
coefNLogConc[20,3] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[20,4] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[20,5] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[20,6] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[21,1] <- sub$Depth[1]
coefNLogConc[21,2] <- "SA_N_gkg"
coefNLogConc[21,3] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[21,4] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[21,5] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[21,6] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[22,1] <- sub$Depth[1]
coefNLogConc[22,2] <- "SC_N_gkg"
coefNLogConc[22,3] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[22,4] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[22,5] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[22,6] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefNLogConc[24,1] <- sub$Depth[1]
coefNLogConc[24,2] <- "TN_gkg"
coefNLogConc[24,3] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[24,4] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[24,5] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[24,6] <- summary(lm(log(TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[25,1] <- sub$Depth[1]
coefNLogConc[25,2] <- "POM_N_gkg"
coefNLogConc[25,3] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[25,4] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[25,5] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[25,6] <- summary(lm(log(POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[26,1] <- sub$Depth[1]
coefNLogConc[26,2] <- "SA_N_gkg"
coefNLogConc[26,3] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[26,4] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[26,5] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[26,6] <- summary(lm(log(SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLogConc[27,1] <- sub$Depth[1]
coefNLogConc[27,2] <- "SC_N_gkg"
coefNLogConc[27,3] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLogConc[27,4] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLogConc[27,5] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLogConc[27,6] <- summary(lm(log(SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefNLogConc) <- c("Depth",
                           "Fraction",
                           "Intercept",
                           "Slope",
                           "R2",
                           "p")


# linear model of concentrations of the fractions, untransformed data
coefNLinConc <- data.frame()

sub <-subset(results_HS, Depth == 5)
coefNLinConc[1,1] <- sub$Depth[1]
coefNLinConc[1,2] <- "SOC gkg"
coefNLinConc[1,3] <- summary(lm((TN_gkg) ~ warming, data=sub))$coefficients[1,1]
coefNLinConc[1,4] <- summary(lm((TN_gkg) ~ warming, data=sub))$coefficients[2,1]
coefNLinConc[1,5] <- summary(lm((TN_gkg) ~ warming, data=sub))$r.squared
coefNLinConc[1,6] <- summary(lm((TN_gkg) ~ warming, data=sub))$coefficients[2,4]

coefNLinConc[2,1] <- sub$Depth[1]
coefNLinConc[2,2] <- "POM"
coefNLinConc[2,3] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[2,4] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[2,5] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[2,6] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[3,1] <- sub$Depth[1]
coefNLinConc[3,2] <- "SA_N_gkg"
coefNLinConc[3,3] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[3,4] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[3,5] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[3,6] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[4,1] <- sub$Depth[1]
coefNLinConc[4,2] <- "SC_N_gkg"
coefNLinConc[4,3] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[4,4] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[4,5] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[4,6] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[5,1] <- sub$Depth[1]
coefNLinConc[5,2] <- "rSN_gkg"
coefNLinConc[5,3] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[5,4] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[5,5] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[5,6] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 15)
coefNLinConc[7,1] <- sub$Depth[1]
coefNLinConc[7,2] <- "TN_gkg"
coefNLinConc[7,3] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[7,4] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[7,5] <- summary(lm((TN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[7,6] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[8,1] <- sub$Depth[1]
coefNLinConc[8,2] <- "POM_N_gkg"
coefNLinConc[8,3] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[8,4] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[8,5] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[8,6] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[9,1] <- sub$Depth[1]
coefNLinConc[9,2] <- "SA_N_gkg"
coefNLinConc[9,3] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[9,4] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[9,5] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[9,6] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[10,1] <- sub$Depth[1]
coefNLinConc[10,2] <- "SC_N_gkg"
coefNLinConc[10,3] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[10,4] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[10,5] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[10,6] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[11,1] <- sub$Depth[1]
coefNLinConc[11,2] <- "rSN_gkg"
coefNLinConc[11,3] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[11,4] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[11,5] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[11,6] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 30)
coefNLinConc[13,1] <- sub$Depth[1]
coefNLinConc[13,2] <- "TN_gkg"
coefNLinConc[13,3] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[13,4] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[13,5] <- summary(lm((TN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[13,6] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[14,1] <- sub$Depth[1]
coefNLinConc[14,2] <- "POM_N_gkg"
coefNLinConc[14,3] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[14,4] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[14,5] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[14,6] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[15,1] <- sub$Depth[1]
coefNLinConc[15,2] <- "SA_N_gkg"
coefNLinConc[15,3] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[15,4] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[15,5] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[15,6] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[16,1] <- sub$Depth[1]
coefNLinConc[16,2] <- "SC_N_gkg"
coefNLinConc[16,3] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[16,4] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[16,5] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[16,6] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[17,1] <- sub$Depth[1]
coefNLinConc[17,2] <- "rSN_gkg"
coefNLinConc[17,3] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[17,4] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[17,5] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[17,6] <- summary(lm((rSN_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 50)
coefNLinConc[19,1] <- sub$Depth[1]
coefNLinConc[19,2] <- "TN_gkg"
coefNLinConc[19,3] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[19,4] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[19,5] <- summary(lm((TN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[19,6] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[20,1] <- sub$Depth[1]
coefNLinConc[20,2] <- "POM_N_gkg"
coefNLinConc[20,3] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[20,4] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[20,5] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[20,6] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[21,1] <- sub$Depth[1]
coefNLinConc[21,2] <- "SA_N_gkg"
coefNLinConc[21,3] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[21,4] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[21,5] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[21,6] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[22,1] <- sub$Depth[1]
coefNLinConc[22,2] <- "SC_N_gkg"
coefNLinConc[22,3] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[22,4] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[22,5] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[22,6] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

sub <-subset(results_HS, Depth == 70)
coefNLinConc[24,1] <- sub$Depth[1]
coefNLinConc[24,2] <- "TN_gkg"
coefNLinConc[24,3] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[24,4] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[24,5] <- summary(lm((TN_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[24,6] <- summary(lm((TN_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[25,1] <- sub$Depth[1]
coefNLinConc[25,2] <- "POM_N_gkg"
coefNLinConc[25,3] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[25,4] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[25,5] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[25,6] <- summary(lm((POM_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[26,1] <- sub$Depth[1]
coefNLinConc[26,2] <- "SA_N_gkg"
coefNLinConc[26,3] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[26,4] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[26,5] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[26,6] <- summary(lm((SA_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

coefNLinConc[27,1] <- sub$Depth[1]
coefNLinConc[27,2] <- "SC_N_gkg"
coefNLinConc[27,3] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[1,1]
coefNLinConc[27,4] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,1]
coefNLinConc[27,5] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$r.squared
coefNLinConc[27,6] <- summary(lm((SC_N_gkg) ~ warming  , data=sub))$coefficients[2,4]

colnames(coefNLinConc) <- c("Depth",
                           "Fraction",
                           "Intercept",
                           "Slope",
                           "R2",
                           "p")

model_output <- list("Linear C stocks" = coefLin, "log_Transformed C stocks" = coefLog,
                     "Linear C concentrations" = coefLinConc, "log_Transformed C concentrations" = coefLogConc,
                     "Linear N stocks" = coefNLin, "log_Transformed N stocks" = coefNLog,
                     "Linear N concentrations" = coefNLinConc, "log_Transformed N concentrations" = coefNLogConc)
write.xlsx(model_output, "Model.xlsx")
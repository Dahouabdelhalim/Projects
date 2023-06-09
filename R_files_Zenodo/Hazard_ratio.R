source ("Library.R")


####Figure 1 ####

a <- read_excel("copulation_rate.xlsx") %>%
  dplyr::mutate(ctime=ctime/60) %>%
  dplyr::mutate(event=ifelse(ctime>30, 0, 1)) %>%
  dplyr::filter(ctime >0) %>%
  dplyr::mutate(ID = row_number())
a$IPI <- as.factor(as.character(a$IPI))
a$species <- as.factor(a$species)

IPI_order <- c("35", "15", "55", "75", "95")
species_order <- c("mel", "sim")
a$IPI <- factor(a$IPI, levels = IPI_order)
a$species <- factor(a$species, levels = species_order)


surv_object <- Surv(time = a$ctime, event = a$event)
fit.coxph <- coxph(Surv(ctime, event) ~ IPI * species,
                   data = a)

fit.zph <- cox.zph(fit.coxph)
fit.zph


#Divide time to keep hazard proportionality

a.split <- survSplit(Surv(ctime, event)~.,
                     data = a, cut = c(7, 31),
                     episode = "tgroup",
                     id ="ID1")

a.split0 <- subset(a.split, tstart == 0)#0-7分
a.split7 <- subset(a.split, tstart == 7)#7-30分


fit.split0 <- coxph(Surv(ctime, event) ~ 
                      IPI*species,
                    data = a.split0)
cox.zph(fit.split0)

fit.split7 <- coxph(Surv(ctime, event) ~ 
                       IPI * species,
                     data = a.split7)
cox.zph(fit.split7)

sum_fit0 <- summary(fit.split0)
sum_fit7 <- summary(fit.split7)

HR_data = NULL
group_order = c("IPI(35-15ms)", "IPI(35-55ms)","IPI(35-75ms)", "IPI(35-95ms)",
                "melanogaster-simulans", "(35-15)*(mel-sim)", "(35-55)*(mel-sim)",
                "(35-75)*(mel-sim)", "(35-95)*(mel-sim)")



.HR_data = data_frame(Time = "0-7 min", 
                      group = group_order,
                      HR = sum_fit0$conf.int[c(1:9)],
                      lower.95 = sum_fit0$conf.int[c(19:27)],
                      upper.95 = sum_fit0$conf.int[c(28:36)],
                      p.value = sum_fit0$coefficients[c(37:45)])
..HR_data = data_frame(Time = "7-30 min", 
                       group = group_order,
                      HR = sum_fit7$conf.int[c(1:9)],
                      lower.95 = sum_fit7$conf.int[c(19:27)],
                      upper.95 = sum_fit7$conf.int[c(28:36)],
                      p.value = sum_fit7$coefficients[c(37:45)])
HR_data <- rbind(.HR_data, ..HR_data)

write.csv(HR_data, "Hazard_ratio.csv")

####Figure S1B ####

a <- read_excel("copulation_rate_S1B.xlsx") %>%
  dplyr::mutate(ctime=ctime/60) %>%
  dplyr::mutate(event=ifelse(ctime>30, 0, 1)) %>%
  dplyr::filter(ctime >0) %>%
  dplyr::mutate(ID = row_number())
a$IPI <- as.factor(as.character(a$IPI))
a$species <- as.factor(a$species)

IPI_order <- c("35", "15", "55")
species_order <- c("mel", "sim")
a$IPI <- factor(a$IPI, levels = IPI_order)
a$species <- factor(a$species, levels = species_order)


surv_object <- Surv(time = a$ctime, event = a$event)
fit.coxph <- coxph(Surv(ctime, event) ~ IPI * species,
                   data = a)

fit.zph <- cox.zph(fit.coxph)
fit.zph


#Divide time to keep hazard proportionality

a.split <- survSplit(Surv(ctime, event)~.,
                     data = a, cut = c(7, 31),
                     episode = "tgroup",
                     id ="ID1")

a.split0 <- subset(a.split, tstart == 0)#0-7分
a.split7 <- subset(a.split, tstart == 7)#7-30分


fit.split0 <- coxph(Surv(ctime, event) ~ 
                      IPI*species,
                    data = a.split0)
cox.zph(fit.split0)

fit.split7 <- coxph(Surv(ctime, event) ~ 
                       IPI * species,
                     data = a.split7)
cox.zph(fit.split7)

sum_fit0 <- summary(fit.split0)
sum_fit7 <- summary(fit.split7)

HR_data = NULL
group_order = c("IPI(35-15ms)", "IPI(35-55ms)",
                "melanogaster-simulans", "(35-15)*(mel-sim)", "(35-55)*(mel-sim)")



.HR_data = data_frame(Time = "0-7 min", 
                      group = group_order,
                      HR = sum_fit0$conf.int[c(1:5)],
                      lower.95 = sum_fit0$conf.int[c(11:15)],
                      upper.95 = sum_fit0$conf.int[c(16:20)],
                      p.value = sum_fit0$coefficients[c(21:25)])
..HR_data = data_frame(Time = "7-30 min", 
                       group = group_order,
                       HR = sum_fit7$conf.int[c(1:5)],
                       lower.95 = sum_fit7$conf.int[c(11:15)],
                       upper.95 = sum_fit7$conf.int[c(16:20)],
                       p.value = sum_fit7$coefficients[c(21:25)])
HR_data <- rbind(.HR_data, ..HR_data)

write.csv(HR_data, "Hazard_ratio_S1B.csv")





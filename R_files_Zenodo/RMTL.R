source ("Library.R")

#####Figure 1#####

RMTL = NULL
data <- read_excel("copulation_rate.xlsx") %>%
  dplyr::mutate(ctime=ctime/60) %>% 
  dplyr:: mutate(event=ifelse(ctime>30, 0, 1)) %>%
  dplyr:: filter(ctime > 0)

mel1535 <- data[(data$species == "mel" & (data$IPI ==15 | data$IPI ==35)),] %>% 
  dplyr:: mutate(sound=ifelse(IPI==15, 0, 1))
rfit_mel1535 <- rmst2(mel1535$ctime, mel1535$event, mel1535$sound)

.RMTL<- data.frame(species = "mel", IPI = 15, 
                   estimate = rfit_mel1535$RMST.arm0$rmtl[1],
                   lower.95 = rfit_mel1535$RMST.arm0$rmtl[3],
                   upper.95 = rfit_mel1535$RMST.arm0$rmtl[4],
                   se = rfit_mel1535$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "mel", IPI = 35, 
                     estimate = rfit_mel1535$RMST.arm1$rmtl[1],
                     lower.95 = rfit_mel1535$RMST.arm1$rmtl[3],
                     upper.95 = rfit_mel1535$RMST.arm1$rmtl[4],
                     se = rfit_mel1535$RMST.arm1$rmtl[2])
RMTL <- rbind(.RMTL, ..RMTL)

mel5575 <- data[(data$species == "mel" &(data$IPI ==55 | data$IPI ==75)),] %>% 
  dplyr:: mutate(sound=ifelse(IPI==55, 0, 1))
rfit_mel5575 <- rmst2(mel5575$ctime, mel5575$event, mel5575$sound)

.RMTL<- data.frame(species = "mel", IPI = 55, 
                   estimate = rfit_mel5575$RMST.arm0$rmtl[1],
                   lower.95 = rfit_mel5575$RMST.arm0$rmtl[3],
                   upper.95 = rfit_mel5575$RMST.arm0$rmtl[4],
                   se = rfit_mel5575$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "mel", IPI = 75, 
                     estimate = rfit_mel5575$RMST.arm1$rmtl[1],
                     lower.95 = rfit_mel5575$RMST.arm1$rmtl[3],
                     upper.95 = rfit_mel5575$RMST.arm1$rmtl[4],
                     se = rfit_mel5575$RMST.arm1$rmtl[2])
RMTL <- rbind(RMTL, .RMTL, ..RMTL)

mel95sim15 <- rbind(data[(data$species == "mel" &data$IPI ==95)|(data$species == "sim" & data$IPI ==15),]) %>% 
  dplyr:: mutate(sound=ifelse(IPI==95, 0, 1))
rfit_mel95sim15 <- rmst2(mel95sim15$ctime, mel95sim15$event, mel95sim15$sound)

.RMTL<- data.frame(species = "mel", IPI = 95, 
                   estimate = rfit_mel95sim15$RMST.arm0$rmtl[1],
                   lower.95 = rfit_mel95sim15$RMST.arm0$rmtl[3],
                   upper.95 = rfit_mel95sim15$RMST.arm0$rmtl[4],
                   se = rfit_mel95sim15$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "sim", IPI = 15, 
                     estimate = rfit_mel95sim15$RMST.arm1$rmtl[1],
                     lower.95 = rfit_mel95sim15$RMST.arm1$rmtl[3],
                     upper.95 = rfit_mel95sim15$RMST.arm1$rmtl[4],
                     se = rfit_mel95sim15$RMST.arm1$rmtl[2])
RMTL <- rbind(RMTL, .RMTL, ..RMTL)

sim3555 <- data[data$species == "sim" & (data$IPI ==35 | data$IPI ==55),] %>% 
  dplyr:: mutate(sound=ifelse(IPI==35, 0, 1))
rfit_sim3555 <- rmst2(sim3555$ctime, sim3555$event, sim3555$sound)

.RMTL<- data.frame(species = "sim", IPI = 35, 
                   estimate = rfit_sim3555$RMST.arm0$rmtl[1],
                   lower.95 = rfit_sim3555$RMST.arm0$rmtl[3],
                   upper.95 = rfit_sim3555$RMST.arm0$rmtl[4],
                   se = rfit_sim3555$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "sim", IPI = 55, 
                     estimate = rfit_sim3555$RMST.arm1$rmtl[1],
                     lower.95 = rfit_sim3555$RMST.arm1$rmtl[3],
                     upper.95 = rfit_sim3555$RMST.arm1$rmtl[4],
                     se = rfit_sim3555$RMST.arm1$rmtl[2])
RMTL <- rbind(RMTL, .RMTL, ..RMTL)

sim7595 <- data[data$species == "sim" & (data$IPI ==75 | data$IPI ==95),] %>% 
  dplyr:: mutate(sound=ifelse(IPI==75, 0, 1))
rfit_sim7595 <- rmst2(sim7595$ctime, sim7595$event, sim7595$sound)

.RMTL<- data.frame(species = "sim", IPI = 75, 
                   estimate = rfit_sim7595$RMST.arm0$rmtl[1],
                   lower.95 = rfit_sim7595$RMST.arm0$rmtl[3],
                   upper.95 = rfit_sim7595$RMST.arm0$rmtl[4],
                   se = rfit_sim7595$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "sim", IPI = 95, 
                     estimate = rfit_sim7595$RMST.arm1$rmtl[1],
                     lower.95 = rfit_sim7595$RMST.arm1$rmtl[3],
                     upper.95 = rfit_sim7595$RMST.arm1$rmtl[4],
                     se = rfit_sim7595$RMST.arm1$rmtl[2])
RMTL <- rbind(RMTL, .RMTL, ..RMTL)

write.csv(RMTL, "RMTL.csv",  row.names=FALSE)

#####Figure S1B#####


RMTL = NULL
data <- read_excel("copulation_rate_S1B.xlsx") %>%
  dplyr::mutate(ctime=ctime/60) %>% 
  dplyr:: mutate(event=ifelse(ctime>30, 0, 1)) %>%
  dplyr:: filter(ctime > 0)

mel1535 <- data[(data$species == "mel" & (data$IPI ==15 | data$IPI ==35)),] %>% 
  dplyr:: mutate(sound=ifelse(IPI==15, 0, 1))
rfit_mel1535 <- rmst2(mel1535$ctime, mel1535$event, mel1535$sound)

.RMTL<- data.frame(species = "mel", IPI = 15, 
                   estimate = rfit_mel1535$RMST.arm0$rmtl[1],
                   lower.95 = rfit_mel1535$RMST.arm0$rmtl[3],
                   upper.95 = rfit_mel1535$RMST.arm0$rmtl[4],
                   se = rfit_mel1535$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "mel", IPI = 35, 
                     estimate = rfit_mel1535$RMST.arm1$rmtl[1],
                     lower.95 = rfit_mel1535$RMST.arm1$rmtl[3],
                     upper.95 = rfit_mel1535$RMST.arm1$rmtl[4],
                     se = rfit_mel1535$RMST.arm1$rmtl[2])
RMTL <- rbind(.RMTL, ..RMTL)


mel55sim15 <- rbind(data[(data$species == "mel" &data$IPI ==55)|(data$species == "sim" & data$IPI ==15),]) %>% 
  dplyr:: mutate(sound=ifelse(IPI==55, 0, 1))
rfit_mel55sim15 <- rmst2(mel55sim15$ctime, mel55sim15$event, mel55sim15$sound)

.RMTL<- data.frame(species = "mel", IPI = 55, 
                   estimate = rfit_mel55sim15$RMST.arm0$rmtl[1],
                   lower.95 = rfit_mel55sim15$RMST.arm0$rmtl[3],
                   upper.95 = rfit_mel55sim15$RMST.arm0$rmtl[4],
                   se = rfit_mel55sim15$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "sim", IPI = 15, 
                     estimate = rfit_mel95sim15$RMST.arm1$rmtl[1],
                     lower.95 = rfit_mel95sim15$RMST.arm1$rmtl[3],
                     upper.95 = rfit_mel95sim15$RMST.arm1$rmtl[4],
                     se = rfit_mel55sim15$RMST.arm1$rmtl[2])
RMTL <- rbind(RMTL, .RMTL, ..RMTL)

sim3555 <- data[data$species == "sim" & (data$IPI ==35 | data$IPI ==55),] %>% 
  dplyr:: mutate(sound=ifelse(IPI==35, 0, 1))
rfit_sim3555 <- rmst2(sim3555$ctime, sim3555$event, sim3555$sound)

.RMTL<- data.frame(species = "sim", IPI = 35, 
                   estimate = rfit_sim3555$RMST.arm0$rmtl[1],
                   lower.95 = rfit_sim3555$RMST.arm0$rmtl[3],
                   upper.95 = rfit_sim3555$RMST.arm0$rmtl[4],
                   se = rfit_sim3555$RMST.arm0$rmtl[2])
..RMTL <- data.frame(species = "sim", IPI = 55, 
                     estimate = rfit_sim3555$RMST.arm1$rmtl[1],
                     lower.95 = rfit_sim3555$RMST.arm1$rmtl[3],
                     upper.95 = rfit_sim3555$RMST.arm1$rmtl[4],
                     se = rfit_sim3555$RMST.arm1$rmtl[2])
RMTL <- rbind(RMTL, .RMTL, ..RMTL)


write.csv(RMTL, "RMTL_S1B.csv",  row.names=FALSE)




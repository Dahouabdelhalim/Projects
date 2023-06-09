source ("Library.R")

#####melanogaster####

data= NULL
for (i in 1:8) {
  sname <- paste("Fly",i, sep="") 
  .data = read.xlsx("raw_data_imaging_3E_melanogaster.xlsx",sheetName= sname) %>%
    dplyr::mutate(Fly = paste(i)) %>%
    tidyr::gather(key = Frequency, value = response, -Frame, -Fly)
  .data$Frequency = str_sub(.data$Frequency, start = 2)
  background = .data %>%
    dplyr::filter(Frame < 20) %>%
    dplyr::group_by(Frequency) %>%
    dplyr::summarise(background = mean(response))
  ..data = full_join(.data, background) %>%
    mutate(norm_response = (response-background)/background)
  data = rbind (data, ..data) # bind the files
}
head(data)
data $Fly = as.numeric(data$Fly)

.data = data %>%
  dplyr::filter(grepl("Hz", Frequency)) %>%
  tidyr::separate(col = Frequency, into = c("Frequency", "Hz", "rep"), sep = "_") %>%
  dplyr::select(-Hz) 
head(.data)

.data$Fly = as.factor(.data$Fly)
Frequency_order = as.character(c(40, seq(100, 300, by=100)))
.data$Hz = factor(.data$Frequency, levels = Frequency_order)

fnrollmean <- function (x) {
  if (length(x) < 3) {
    rep(NA,length(x)) 
  } else {
    rollmean(x,3,align="center",na.pad=TRUE)
  }
}
.data <- .data %>% group_by(Fly,Frequency) %>% 
  mutate(rollavg=fnrollmean(norm_response))


.data.mean = .data %>% group_by(Frame, Fly, Frequency) %>%
  dplyr::summarise(mean.response = mean(rollavg), sd.response = sd(rollavg))
head(.data.mean)


.data.for.model = .data %>%
  filter(Frequency != "test") %>% 
  mutate(norm_response.for.model = ifelse((Frame < 20 |Frame >130), rollavg, NA)) #select Frame 1~9, 41~50
head(.data.for.model)
.data.for.model <- .data.for.model[!(.data.for.model$Frame==1),]

model = .data.for.model %>% group_by(Fly, Frequency, rep) %>%
  do(fit = nls(norm_response.for.model ~ a * Frame + b, data = ., start = list(a = 1, b = 0.1))) %>% #fit to lm
  mutate(tidys = list(broom::tidy(fit))) %>%
  unnest(tidys)
head(model)

model.coef = model %>% 
  dplyr::select(Fly, Frequency, rep, term, estimate) %>%
  tidyr::spread(key = term, value = estimate)
head(model.coef) 

data.fitting = .data %>%
  dplyr::filter(Frequency != "test") %>%
  dplyr::select(-response, -background) %>% 
  left_join(., model.coef, by = c("Fly", "Frequency", "rep")) %>%
  dplyr::mutate(fit_response = a * Frame + b) %>%
  dplyr::mutate(norm_norm_response = rollavg - fit_response)
head(data.fitting)


data.fitting.mean = data.fitting %>% dplyr::filter(Frame > 1, Frame < 150) %>% group_by(Fly, Frame, Frequency) %>%
  dplyr::summarise(average.response = mean(norm_norm_response), sd.response = sd(norm_norm_response))
data.fitting.mean =  dplyr::mutate(data.fitting.mean, time = Frame/10)
head(data.fitting.mean)
data.fitting.mean$Frequency = as.numeric(data.fitting.mean$Frequency)
data.fitting.mean[order(data.fitting.mean$Frequency, decreasing=F),]

write.csv(data.fitting.mean, "timetrace_imaging_3E_melanogaster.csv")


max.data.fitting.mean = data.fitting.mean %>%
  group_by(Fly, Frequency) %>%
  summarize(max = max(average.response, na.rm = TRUE))
integral = max.data.fitting.mean %>% dplyr::group_by(Fly) %>%
  dplyr::summarise(integral = sum(max))
max.data.fitting.mean = full_join(max.data.fitting.mean, integral) %>%
  mutate(normmax = max/integral) %>%
  group_by(Frequency) %>%
  mutate(rank = row_number(Fly)) %>%
  dplyr::select(ID = rank, dplyr::everything())

write.csv(max.data.fitting.mean, "normmax_imaging_3F_melanogaster.csv")

#####simulans####

data= NULL
for (i in 1:8) {
  sname <- paste("Fly",i, sep="") 
  .data = read.xlsx("raw_data_imaging_3E_simulans.xlsx",sheetName= sname) %>%
    dplyr::mutate(Fly = paste(i)) %>%
    tidyr::gather(key = Frequency, value = response, -Frame, -Fly)
  .data$Frequency = str_sub(.data$Frequency, start = 2)
  background = .data %>%
    dplyr::filter(Frame < 20) %>%
    dplyr::group_by(Frequency) %>%
    dplyr::summarise(background = mean(response))
  ..data = full_join(.data, background) %>%
    mutate(norm_response = (response-background)/background)
  data = rbind (data, ..data) # bind the files
}
head(data)
data $Fly = as.numeric(data$Fly)

.data = data %>%
  dplyr::filter(grepl("Hz", Frequency)) %>%
  tidyr::separate(col = Frequency, into = c("Frequency", "Hz", "rep"), sep = "_") %>%
  dplyr::select(-Hz) 
head(.data)

.data$Fly = as.factor(.data$Fly)
Frequency_order = as.character(c(40, seq(100, 300, by=100)))
.data$Hz = factor(.data$Frequency, levels = Frequency_order)

fnrollmean <- function (x) {
  if (length(x) < 3) {
    rep(NA,length(x)) 
  } else {
    rollmean(x,3,align="center",na.pad=TRUE)
  }
}
.data <- .data %>% group_by(Fly,Frequency) %>% 
  mutate(rollavg=fnrollmean(norm_response))


.data.mean = .data %>% group_by(Frame, Fly, Frequency) %>%
  dplyr::summarise(mean.response = mean(rollavg), sd.response = sd(rollavg))
head(.data.mean)


.data.for.model = .data %>%
  filter(Frequency != "test") %>% 
  mutate(norm_response.for.model = ifelse((Frame < 20 |Frame >130), rollavg, NA)) #select Frame 1~9, 41~50
head(.data.for.model)
.data.for.model <- .data.for.model[!(.data.for.model$Frame==1),]

model = .data.for.model %>% group_by(Fly, Frequency, rep) %>%
  do(fit = nls(norm_response.for.model ~ a * Frame + b, data = ., start = list(a = 1, b = 0.1))) %>% #fit to lm
  mutate(tidys = list(broom::tidy(fit))) %>%
  unnest(tidys)
head(model)

model.coef = model %>% 
  dplyr::select(Fly, Frequency, rep, term, estimate) %>%
  tidyr::spread(key = term, value = estimate)
head(model.coef) 

data.fitting = .data %>%
  dplyr::filter(Frequency != "test") %>%
  dplyr::select(-response, -background) %>% 
  left_join(., model.coef, by = c("Fly", "Frequency", "rep")) %>%
  dplyr::mutate(fit_response = a * Frame + b) %>%
  dplyr::mutate(norm_norm_response = rollavg - fit_response)
head(data.fitting)


data.fitting.mean = data.fitting %>% dplyr::filter(Frame > 1, Frame < 150) %>% group_by(Fly, Frame, Frequency) %>%
  dplyr::summarise(average.response = mean(norm_norm_response), sd.response = sd(norm_norm_response))
data.fitting.mean =  dplyr::mutate(data.fitting.mean, time = Frame/10)
head(data.fitting.mean)
data.fitting.mean$Frequency = as.numeric(data.fitting.mean$Frequency)
data.fitting.mean[order(data.fitting.mean$Frequency, decreasing=F),]

write.csv(data.fitting.mean, "timetrace_imaging_3E_simulans.csv")


max.data.fitting.mean = data.fitting.mean %>%
  group_by(Fly, Frequency) %>%
  summarize(max = max(average.response, na.rm = TRUE))
integral = max.data.fitting.mean %>% dplyr::group_by(Fly) %>%
  dplyr::summarise(integral = sum(max))
max.data.fitting.mean = full_join(max.data.fitting.mean, integral) %>%
  mutate(normmax = max/integral) %>%
  group_by(Frequency) %>%
  mutate(rank = row_number(Fly)) %>%
  dplyr::select(ID = rank, dplyr::everything())

write.csv(max.data.fitting.mean, "normmax_imaging_3E_simulans.csv")






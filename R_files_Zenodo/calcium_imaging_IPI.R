source ("Library.R")

#####melanoagster#####

data= NULL
for (i in 1:12) {
  sname <- paste("Fly",i, sep="") 
  .data = read.xlsx("raw_data_imaging_4B_melanogaster.xlsx",sheetName= sname) %>%
    dplyr::mutate(Fly = paste(i)) %>%
    tidyr::gather(key = IPI, value = response, -Frame, -Fly)
  background = .data %>%
    dplyr::filter(Frame < 20) %>%
    dplyr::group_by(IPI) %>%
    dplyr::summarise(background = mean(response))
  ..data = full_join(.data, background) %>%
    mutate(norm_response = (response-background)/background)
  data = rbind (data, ..data) # bind the files
}
head(data)
data $Fly = as.numeric(data$Fly)


data $Fly = as.numeric(data$Fly)

.data = data %>%
  dplyr::filter(grepl("IPI", IPI)) %>%
  tidyr::separate(col = IPI, into = c("value", "IPI", "rep"), sep = "_") %>%
  dplyr::select(-value) 
head(.data)

.data$Fly = as.factor(.data$Fly)
IPI_order = as.character(c(seq(15, 105, by=10)))
.data$IPI = factor(.data$IPI, levels = IPI_order)

fnrollmean <- function (x) {
  if (length(x) < 3) {
    rep(NA,length(x)) 
  } else {
    rollmean(x,3,align="center",na.pad=TRUE)
  }
}

.data <- .data %>% group_by(Fly,IPI) %>% 
  mutate(rollavg=fnrollmean(norm_response))
head(.data)


.data.mean = .data %>% group_by(Frame, Fly, IPI) %>%
  dplyr::summarise(mean.response = mean(rollavg), sd.response = sd(rollavg))
head(.data.mean)

.data.for.model = .data %>%
  mutate(norm_response.for.model = ifelse((Frame < 10 |Frame >40), rollavg, NA)) #select Frame 1~9, 41~50
head(.data.for.model)
.data.for.model <- .data.for.model[!(.data.for.model$Frame==1),]

model = .data.for.model %>% group_by(Fly, IPI, rep) %>%
  do(fit = nls(norm_response.for.model ~ a * Frame + b, data = ., start = list(a = 1, b = 0.1))) %>% #fit to lm
  mutate(tidys = list(broom::tidy(fit))) %>%
  unnest(tidys)
head(model)

model.coef = model %>% 
  dplyr::select(Fly, IPI, rep, term, estimate) %>%
  tidyr::spread(key = term, value = estimate)
head(model.coef) 

data.fitting = .data %>%
  dplyr::filter(IPI != "test") %>%
  dplyr::select(-response, -background) %>% 
  left_join(., model.coef, by = c("Fly", "IPI", "rep")) %>%
  dplyr::mutate(fit_response = a * Frame + b) %>%
  dplyr::mutate(norm_norm_response = rollavg - fit_response)
head(data.fitting)


data.fitting.mean = data.fitting %>% dplyr::filter(Frame > 1, Frame < 50) %>% group_by(Frame, Fly, IPI) %>%
  dplyr::summarise(average.response = mean(norm_norm_response), sd.response = sd(norm_norm_response))
data.fitting.mean =  dplyr::mutate(data.fitting.mean, time = Frame/10)
head(data.fitting.mean)

write.csv(data.fitting.mean, "timetrace_imaging_4B_melanogaster.csv")

max.data.fitting.mean = data.fitting.mean %>%
  group_by(Fly, IPI) %>%
  summarize(max = max(average.response, na.rm = TRUE))
integral = max.data.fitting.mean %>% dplyr::group_by(Fly) %>%
  dplyr::summarise(integral = sum(max))
max.data.fitting.mean = full_join(max.data.fitting.mean, integral) %>%
  mutate(normmax = max/integral) %>%
  group_by(IPI) %>%
  mutate(rank = row_number(Fly)) %>%
  dplyr::select(ID = rank, dplyr::everything())

write.csv(max.data.fitting.mean, "normmax_imaging_4C_melanogaster.csv")

#####simulans####

data= NULL
for (i in 1:13) {
  sname <- paste("Fly",i, sep="") 
  .data = read.xlsx("raw_data_imaging_4B_simulans.xlsx",sheetName= sname) %>%
    dplyr::mutate(Fly = paste(i)) %>%
    tidyr::gather(key = IPI, value = response, -Frame, -Fly)
  background = .data %>%
    dplyr::filter(Frame < 20) %>%
    dplyr::group_by(IPI) %>%
    dplyr::summarise(background = mean(response))
  ..data = full_join(.data, background) %>%
    mutate(norm_response = (response-background)/background)
  data = rbind (data, ..data) # bind the files
}
head(data)
data $Fly = as.numeric(data$Fly)


data $Fly = as.numeric(data$Fly)

.data = data %>%
  dplyr::filter(grepl("IPI", IPI)) %>%
  tidyr::separate(col = IPI, into = c("value", "IPI", "rep"), sep = "_") %>%
  dplyr::select(-value) 
head(.data)

.data$Fly = as.factor(.data$Fly)
IPI_order = as.character(c(seq(15, 105, by=10)))
.data$IPI = factor(.data$IPI, levels = IPI_order)

fnrollmean <- function (x) {
  if (length(x) < 3) {
    rep(NA,length(x)) 
  } else {
    rollmean(x,3,align="center",na.pad=TRUE)
  }
}

.data <- .data %>% group_by(Fly,IPI) %>% 
  mutate(rollavg=fnrollmean(norm_response))
head(.data)


.data.mean = .data %>% group_by(Frame, Fly, IPI) %>%
  dplyr::summarise(mean.response = mean(rollavg), sd.response = sd(rollavg))
head(.data.mean)

.data.for.model = .data %>%
  mutate(norm_response.for.model = ifelse((Frame < 10 |Frame >40), rollavg, NA)) #select Frame 1~9, 41~50
head(.data.for.model)
.data.for.model <- .data.for.model[!(.data.for.model$Frame==1),]

model = .data.for.model %>% group_by(Fly, IPI, rep) %>%
  do(fit = nls(norm_response.for.model ~ a * Frame + b, data = ., start = list(a = 1, b = 0.1))) %>% #fit to lm
  mutate(tidys = list(broom::tidy(fit))) %>%
  unnest(tidys)
head(model)

model.coef = model %>% 
  dplyr::select(Fly, IPI, rep, term, estimate) %>%
  tidyr::spread(key = term, value = estimate)
head(model.coef) 

data.fitting = .data %>%
  dplyr::filter(IPI != "test") %>%
  dplyr::select(-response, -background) %>% 
  left_join(., model.coef, by = c("Fly", "IPI", "rep")) %>%
  dplyr::mutate(fit_response = a * Frame + b) %>%
  dplyr::mutate(norm_norm_response = rollavg - fit_response)
head(data.fitting)


data.fitting.mean = data.fitting %>% dplyr::filter(Frame > 1, Frame < 50) %>% group_by(Frame, Fly, IPI) %>%
  dplyr::summarise(average.response = mean(norm_norm_response), sd.response = sd(norm_norm_response))
data.fitting.mean =  dplyr::mutate(data.fitting.mean, time = Frame/10)
head(data.fitting.mean)

write.csv(data.fitting.mean, "timetrace_imaging_4B_simulans.csv")

max.data.fitting.mean = data.fitting.mean %>%
  group_by(Fly, IPI) %>%
  summarize(max = max(average.response, na.rm = TRUE))
integral = max.data.fitting.mean %>% dplyr::group_by(Fly) %>%
  dplyr::summarise(integral = sum(max))
max.data.fitting.mean = full_join(max.data.fitting.mean, integral) %>%
  mutate(normmax = max/integral) %>%
  group_by(IPI) %>%
  mutate(rank = row_number(Fly)) %>%
  dplyr::select(ID = rank, dplyr::everything())

write.csv(max.data.fitting.mean, "normmax_imaging_4C_simulans.csv")

#####delta response (Figure 4D)#####


st_data <- read.csv("normmax_imaging_4C_melanogaster.csv") %>% mutate(species = "mel")
st_data = st_data[, colnames(st_data) %in% c("species","ID", "IPI", "normmax")]
.st_data <- read.csv("normmax_imaging_4C_simulans.csv") %>% mutate(species = "sim")
.st_data = .st_data[, colnames(.st_data) %in% c("species","ID", "IPI", "normmax")]
st_data <- rbind(st_data, .st_data)
head(st_data)


order_f <- c(15, seq(35, 105, by = 10))

st = st_data %>% filter (IPI == 25)
names(st)[ which( names(st)=="normmax" ) ] <- "max25"

st_mel = st %>% filter(species == "mel")
st_sim = st %>% filter(species == "sim")

st_gp = st_data %>% filter (IPI == 25)
names(st_gp)[ which( names(st_gp)=="normmax" ) ] <- "max25"
st_gp25 = st_gp[, colnames(st_gp) == "max25"]

st_gp = NULL
for(i in order_f){
  st_i = st_data %>% filter(IPI == i)
  st_i = st_i[, colnames(st_i) %in% c("species","ID", "IPI", "normmax")]
  st_ij = cbind(st_i, st_gp25) %>%
    mutate(df = normmax - st_gp25) 
  st_ij = st_ij[, colnames(st_ij) %in% c("species","ID", "IPI", "df")]
  st_gp = rbind(st_gp, st_ij)
}

write.csv(st_gp, "delta_response_4D.csv")



source ("Library.R")

# read data
data= NULL
for (i in 1:8) {
  sname <- paste("Fly",i, sep="") 
  .data = read.xlsx("raw_data_imaging_S4B.xlsx",sheetName= sname) %>%
    dplyr::mutate(Fly = paste(i)) %>%
    tidyr::gather(key = sound, value = response, -Frame, -Fly) 
  background = .data %>%
    dplyr::filter(Frame < 10) %>%
    dplyr::group_by(sound) %>%
    dplyr::summarise(background = mean(response))
  ..data = full_join(.data, background) %>%
    mutate(norm_response = (response-background)/background)
  data = rbind (data, ..data) 
}
head(data)

data $Fly = as.numeric(data$Fly)

.data = data %>%
  dplyr::filter(grepl("Hz", sound)) %>%
  tidyr::separate(col = sound, into = c("Frequency", "Hz","IPI","ms", "rep"), sep = "_") %>%
  dplyr::select(-Hz, -ms) 
head(.data)

.data$Fly = as.factor(.data$Fly)
Frequency_order = as.character(c(167, 167, 167, 333, 333, 333))
IPI_order = as.character(c(15, 25, 35, 15, 25, 35))
group_order = as.character(c("167_15","167_25","167_35", "333_15", "333_25", "333_35" ))
.data <- .data %>% mutate(group = paste(!!!rlang::syms(c("Frequency", "IPI")), sep="_"))

fnrollmean <- function (x) {
  if (length(x) < 3) {
    rep(NA,length(x)) 
  } else {
    rollmean(x,3,align="center",na.pad=TRUE)
  }
}

.data <- .data %>% group_by(Fly,group) %>% 
  mutate(rollavg=fnrollmean(norm_response))
head(.data)

.data.mean = .data %>% group_by(Frame, Fly, group) %>%
  dplyr::summarise(mean.response = mean(rollavg), sd.response = sd(rollavg))
head(.data.mean)


.data.for.model = .data %>%
  mutate(norm_response.for.model = ifelse((Frame < 10 |Frame >40), rollavg, NA)) #select Frame 1~9, 41~50
head(.data.for.model)
.data.for.model <- .data.for.model[!(.data.for.model$Frame==1),]

model = .data.for.model %>% group_by(Fly, group, rep) %>%
  do(fit = nls(norm_response.for.model ~ a * Frame + b, data = ., start = list(a = 1, b = 0.1))) %>% #fit to lm
  mutate(tidys = list(broom::tidy(fit))) %>%
  unnest(tidys)
head(model)

model.coef = model %>% 
  dplyr::select(Fly, group, rep, term, estimate) %>%
  tidyr::spread(key = term, value = estimate)
head(model.coef) 

data.fitting = .data %>%
  dplyr::filter(group != "test") %>%
  dplyr::select(-response, -background) %>% 
  left_join(., model.coef, by = c("Fly", "group", "rep")) %>%
  dplyr::mutate(fit_response = a * Frame + b) %>%
  dplyr::mutate(norm_norm_response = rollavg - fit_response)
head(data.fitting)


data.fitting.mean = data.fitting %>% dplyr::filter(Frame > 1, Frame < 50) %>% group_by(Fly, Frame, IPI, Frequency, group) %>%
  dplyr::summarise(average.response = mean(norm_norm_response), sd.response = sd(norm_norm_response))
data.fitting.mean =  dplyr::mutate(data.fitting.mean, time = Frame/10)
head(data.fitting.mean)

write.csv(data.fitting.mean, "timetrace_imaging_S4B.csv")


max.data.fitting.mean = data.fitting.mean %>%
  group_by(Fly, IPI, Frequency, group) %>%
  summarize(max = max(average.response, na.rm = TRUE))
integral = max.data.fitting.mean %>% dplyr::group_by(Fly, Frequency) %>%
  dplyr::summarise(integral = sum(max))
max.data.fitting.mean = full_join(max.data.fitting.mean, integral)
max.data.fitting.mean = max.data.fitting.mean %>% mutate(normmax = max/integral)
max.data.fitting.mean  <- max.data.fitting.mean %>%
  group_by(IPI, Frequency) %>%
  mutate(rank = row_number(Fly))
max.data.fitting.mean <- dplyr::select(max.data.fitting.mean, ID = rank, dplyr::everything())

head(max.data.fitting.mean)

write.csv(max.data.fitting.mean, "norm_max_mean_S4B.csv")



st_data <- read.csv("norm_max_mean_S4B.csv")
st_data = st_data[, colnames(st_data) %in% c("ID","group","Frequency", "IPI", "normmax")]
head(st_data)


order_f <- c(15, 35)

st = st_data %>% filter (IPI == 25)
names(st)[ which( names(st)=="normmax" ) ] <- "max25"

st_167 = st %>% filter(Frequency == 167)
st_333 = st %>% filter(Frequency == 333)

st_gp = st_data %>% filter (IPI == 25)
names(st_gp)[ which( names(st_gp)=="normmax" ) ] <- "max25"
st_gp25 = st_gp[, colnames(st_gp) == "max25"]

st_gp = NULL
for(i in order_f){
  st_i = st_data %>% filter(IPI == i)
  st_i = st_i[, colnames(st_i) %in% c("Frequency","ID", "IPI","group", "normmax")]
  st_ij = cbind(st_i, st_gp25) %>%
    mutate(df = normmax - st_gp25) 
  st_ij = st_ij[, colnames(st_ij) %in% c("Frequency","ID", "IPI", "group", "df")]
  st_gp = rbind(st_gp, st_ij)
}

write.csv(st_gp, "delta_response_S4B.csv")



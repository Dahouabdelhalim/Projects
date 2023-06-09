## Dominance Hierarchies and Temperature in Black-capped Chickadees - Simplified

# Adding fonts

#library("showtext")
#font_add_google(name = "Noto Sans", family = "Noto Sans")
#{
#  font_add(family = "Colus", regular = "C:/Users/joshk/Desktop/NewFonts/Colus-Regular.ttf")
#  font_add(family = "LibreB", regular = "C:/Users/joshk/Desktop/NewFonts/LibreBaskerville-Regular.otf")
#  font_add(family = "Gidole", regular = "C:/Users/joshk/Desktop/NewFonts/Gidole-Regular.ttf")
#  font_add(family = "Reforma1918", regular = "C:/Users/joshk/Desktop/NewFonts/Reforma1918-Gris.ttf")
#  font_add(family = "Reforma1969", regular = "C:/Users/joshk/Desktop/NewFonts/Reforma1969-Gris.ttf")
#  font_add(family = "Reforma2018", regular = "C:/Users/joshk/Desktop/NewFonts/Reforma2018-Gris.ttf")
#}
#showtext_auto()

# Pulling in required packages

library('devtools')

#install_github("joshuakrobertson/ReVuePro")
#install_github("samclifford/mgcv.helper")

library("easypackages")
libraries(
  "akima", "aphid", "Biostrings", "doParallel", "dplyr",
  "effects", "EloRating", "emmeans", "ggeffects", "ggplot2", "ggsignif","glmmTMB", "grid",
  "gridExtra", "gtable", "gtools", "itsadug","lme4", "lmerTest", 
  "mgcv", "mgcv.helper", "plotly", "pracma", "ReVuePro", "sjPlot", 
  "splitstackshape","stringr",
  "stringi", "tidyverse"
)

# font_import()
# load_fonts()

# my.theme <- theme(
#   panel.grid.minor = element_blank(),
#   axis.title = element_text(size = 18, family = "Noto Sans"),
#   axis.text = element_text(size = 14, colour = "black", family = "Noto Sans"),
#   axis.title.y = element_text(vjust = 0.5), panel.grid.major =
#     element_line(colour = "grey75"), legend.title = element_text(
#     size = 14,
#     colour = "black", family = "Noto Sans"
#   ), legend.text = element_text(
#     size = 12,
#     colour = "black", family = "Noto Sans"
#   )
# )

my.theme_2 <- theme(
  panel.grid.minor = element_blank(),
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 14, colour = "black"),
  axis.title.y = element_text(vjust = 0.5), panel.grid.major =
    element_line(colour = "grey75"), legend.title = element_text(
    size = 14,
    colour = "black", 
  ), legend.text = element_text(
    size = 12,
    colour = "black"
  )
)

## Additional Functions

{
pull.cat <- function(x) {
  bins <- 5
  increments <- (range(x)[2] - range(x)[1])/(bins - 1)
  to_return <- seq(range(x)[1], range(x)[2], increments)
  return(to_return)
}

up.cat <- function(new_bins) {
  up_bins = new_bins
  body(pull.cat)[[2]] <<- substitute(bins <- up_bins)
}

T2S <- function(x) {
  if (is.na(x) == TRUE) {
    return(NA)
  } else if (length(unlist(strsplit(as.character(x), ":"))) >= 2) {
    cut <- c(unlist(strsplit(as.character(x), ":")))
    Min <- as.numeric(cut[1])
    Sec <- as.numeric(cut[2])
    ABS <- Min * 60 + Sec
    return(ABS)
  } else if (length(unlist(strsplit(as.character(x), ":"))) < 2) {
    if (nchar(as.character(x)) == 6) {
      Hour <- as.numeric(substr(as.character(x), 0, 2))
      Min <- as.numeric(substr(as.character(x), 3, 4))
      Sec <- as.numeric(substr(as.character(x), 5, 6))
      ABS <- Hour * 3600 + Min * 60 + Sec
      return(ABS)
    } else if (nchar(as.character(x)) == 5) {
      PadTime <- str_pad(as.character(x), 6, side = "left", "0")
      Hour <- as.numeric(substr(PadTime, 0, 2))
      Min <- as.numeric(substr(PadTime, 3, 4))
      Sec <- as.numeric(substr(PadTime, 5, 6))
      ABS <- Hour * 3600 + Min * 60 + Sec
      return(ABS)
    } else if (nchar(as.character(x)) > 6) {
      Time_Cut <- c(unlist(strsplit(as.character(x), "_")))
      Base_Time <- Time_Cut[1]
      Add_Time <- Time_Cut[2]
      Hour <- as.numeric(substr(Base_Time, 0, 2))
      Min <- as.numeric(substr(Base_Time, 3, 4))
      Sec <- as.numeric(substr(Base_Time, 5, 6))
      ABS <- Hour * 3600 + Min * 60 + Sec
      Add_Seconds <- as.numeric(Add_Time)
      ABS <- ABS + Add_Seconds
      return(ABS)
    } else {
      return(x)
    }
  }
}

BinF <- function(vect, width, lowest = "TRUE", right = "FALSE") {
  F_breaks <- seq(from = min(vect), to = max(vect), by = width)
  F_labels <- c(as.numeric(as.character(F_breaks))[-1])
  F_bins <- cut(vect, F_breaks, include.lowest = lowest, right = right, labels = F_labels)
  return(F_bins)
}

BinSum <- function(data, bin_col, width, group, response, lowest = "TRUE", right = "FALSE",
                   conf = 95, na.rm = "TRUE", cut = "TRUE") {
  require("dplyr")
  require("ggplot2")
  conf.dec <- (100 - ((100 - conf) / 2)) / 100
  z.est <- qnorm(abs(conf.dec))
  bin_var <- c(with(data, get(bin_col)))
  frame <- as.data.frame(data)

  if (cut == "TRUE") {
    F_breaks <- seq(from = min(bin_var), to = max(bin_var), by = width)
    F_labels <- c(as.numeric(as.character(F_breaks))[-1])
    F_bins <- cut(bin_var, F_breaks, include.lowest = lowest, right = right, labels = F_labels)
    frame$vect <- F_bins
    if (na.rm == "TRUE") {
      frame_sum <- as.data.frame(frame %>% group_by(get(group), vect) %>%
        dplyr::summarise(
          Mean = mean(get(response), na.rm = T),
          Upper.CI = (mean(get(response), na.rm = T) +
            z.est * (sd(na.omit(get(response))) /
              sqrt(mean(get(response), na.rm = T)))),
          Lower.CI = (mean(get(response), na.rm = T) -
            z.est * (sd(na.omit(get(response))) /
              sqrt(mean(get(response), na.rm = T))))
        ))
    } else if (na.rm == "FALSE") {
      frame_sum <- as.data.frame(frame %>% group_by(get(group), vect) %>%
        dplyr::summarise(
          Mean = mean(get(response), na.rm = F),
          Upper.CI = (mean(get(response), na.rm = F) +
            z.est * (sd(get(response)) /
              sqrt(mean(get(response), na.rm = F)))),
          Lower.CI = (mean(get(response), na.rm = F) -
            z.est * (sd(get(response)) /
              sqrt(mean(get(response), na.rm = F))))
        ))
    }
  } else if (cut == "FALSE") {
    if (na.rm == "TRUE") {
      frame_sum <- as.data.frame(frame %>% group_by(get(group), get(bin_col)) %>%
        dplyr::summarise(
          Mean = mean(get(response), na.rm = T),
          Upper.CI = (mean(get(response), na.rm = T) +
            z.est * (sd(na.omit(get(response))) /
              sqrt(mean(get(response), na.rm = T)))),
          Lower.CI = (mean(get(response), na.rm = T) -
            z.est * (sd(na.omit(get(response))) /
              sqrt(mean(get(response), na.rm = T))))
        ))
    } else if (na.rm == "FALSE") {
      frame_sum <- as.data.frame(frame %>% group_by(get(group), get(bin_col)) %>%
        dplyr::summarise(
          Mean = mean(get(response), na.rm = F),
          Upper.CI = (mean(get(response), na.rm = F) +
            z.est * (sd(get(response)) /
              sqrt(mean(get(response), na.rm = F)))),
          Lower.CI = (mean(get(response), na.rm = F) -
            z.est * (sd(get(response)) /
              sqrt(mean(get(response), na.rm = F))))
        ))
    }
  }
  colnames(frame_sum)[1:length(group)] <- paste(group)
  colnames(frame_sum)[(1 + length(group))] <- paste(bin_col, "bin", sep = "_")
  return(frame_sum)
}

Smooth_CI <- function(data, group, xvar, Upper.Conf = "Upper.CI", Lower.Conf = "Lower.CI") {
  require("dplyr")
  dat_list <- split(data, f = with(data, get(group)))
  for (i in 1:length(dat_list)) {
    Upper.Mod <- loess(with(dat_list[[i]], get(Upper.Conf)) ~ with(dat_list[[i]], get(xvar)))
    Lower.Mod <- loess(with(dat_list[[i]], get(Lower.Conf)) ~ with(dat_list[[i]], get(xvar)))
    dat_list[[i]]$Smooth_Upper <- predict(Upper.Mod, type = "response")
    dat_list[[i]]$Smooth_Lower <- predict(Lower.Mod, type = "response")
  }
  dat_return <- bind_rows(dat_list)
  return(dat_return)
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

rotate <- function(x) t(apply(x, 2, rev))

Dom_trans <- function(x){
  m2d = as.data.frame(rotate(x))
  disem = vector('list', ncol(m2d))
  for (i in 1:ncol(m2d)){
    disem[[i]] = data.frame("ID" = as.character(colnames(m2d)[i]),
                            "Elo" = as.integer(m2d[,i]),
                            "Time" = c(1:length(m2d[,i])))
  }
  long_d = bind_rows(disem)
  return(long_d)
}

Elo_Pull = function(win, lose, s.denom = 100, random = TRUE){
  Elo = elo_scores(winners = win, losers = lose, sigmoid.param = I(1/s.denom), randomise = random)
  Trans_Elo = Dom_trans(Elo)
  Elo_Rep = estimate_uncertainty_by_repeatability(winners = win, losers = lose,
        sigmoid.param = I(1/s.denom))
  Elo_Plot = Trans_Elo %>% group_by(ID) %>%
                 dplyr::summarise(Mean = mean(Elo),
                 SE = sd(Elo)/sqrt(n()),
                 Upper.CI = Mean + 1.96*SE,
                 Lower.CI = Mean - 1.96*SE) %>%
                  ggplot(aes(x = ID, y = Mean, fill = ID)) +
                    geom_point(pch = 21, size = 4, colour = "black") +
                    geom_errorbar(aes(x = ID, ymin = Lower.CI, ymax = Upper.CI),
                      size = 1, width = 0.3, colour = "black") +
                      theme_bw()
   Rank_Est = as.data.frame(Trans_Elo %>% group_by(ID) %>%
              dplyr::summarise(Mean_Rank = mean(Elo)) %>%
              arrange(desc(Mean_Rank)) %>%
              mutate("Rank" = c(1:nrow(.))) %>%
              select(ID, Rank))
   Return_Dat = vector('list', 5)
   Return_Dat[[1]] = Elo
   Return_Dat[[2]] = Trans_Elo
   Return_Dat[[3]] = Elo_Rep
   Return_Dat[[4]] = Rank_Est
   Return_Dat[[5]] = Elo_Plot
   names(Return_Dat) = c("Elo.Matrix", "Elo.DF", "Intraclass.r", "Ranks", "Elo.Plot")
   return(Return_Dat)
}

Rand_Elo <- function(data, presence, kvals = NULL, sample_size = 20, niter = 100){
  Elo_means = vector('list', niter)
  if (is.null(kvals) == TRUE){
    for (i in 1:niter){
      while (TRUE){
        sample_dates = sample(unique(data$Date), size = sample_size, replace = F)
        sample_dat = subset(data, Date %in% sample_dates)
        sample_pres = subset(presence, Date >= min(sample_dates) & Date <= max(sample_dates))
        Elo_samples = try(elo.seq(winner = sample_dat$Winner,
                               loser = sample_dat$Loser,
                               Date = sample_dat$Date,
                               presence = sample_pres,
                               draw = sample_dat$Draw), silent = TRUE)
        if (!is(Elo_samples, 'try-error')) break
      }
    Elo_means[[i]] = as.data.frame(as.data.frame(Elo_samples$mat) %>%
                     dplyr::summarise_all(list(~ mean(., na.rm = T)))) %>%
                     mutate("Sample" = i)
    }
  } else if (is.null(kvals) == FALSE & is.list(kvals) == FALSE) {
    return("kvals must be a list containing the names and values of each interaction type.")
  } else {
    for (i in 1:niter){
      while (TRUE){
        sample_dates = sample(unique(data$Date), size = sample_size, replace = F)
        sample_dat = subset(data, Date %in% sample_dates)
        sample_pres = subset(presence, Date >= min(sample_dates) & Date <= max(sample_dates))
        Elo_samples = try(elo.seq(winner = sample_dat$Winner,
                               loser = sample_dat$Loser,
                               Date = sample_dat$Date,
                               presence = sample_pres,
                               draw = sample_dat$Draw,
                               intensity = sample_dat$Type,
                               k = kvals), silent = TRUE)
        if (!is(Elo_samples, 'try-error')) break
      }
    Elo_means[[i]] = as.data.frame(as.data.frame(Elo_samples$mat) %>%
                     dplyr::summarise_all(list(~ mean(., na.rm = T)))) %>%
                     mutate("Sample" = i)
      }
    }
    Elo_full = bind_rows(Elo_means)
    return(Elo_full)
  }
}

# Pulling in Data for Elo calculations

All <- read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Digital Bird With RT.csv")

Pens <- split(All, f = All$Pen_Short)

{
  Bird.IDs <- vector("list", 4)
  Bird.IDs[[1]] <- c("BdABd", "BdAO", "ABlBl", "AYY", "YAO")
  Bird.IDs[[2]] <- c("OOA", "BdABl", "YAO2", "ABlO", "ABdBd")
  Bird.IDs[[3]] <- c("YAR", "RAR", "ABlR", "ARO", "BdABl2")
  Bird.IDs[[4]] <- c("BlAR", "BdAO2", "AOR", "YAY", "AYBd")
}

Pens.Adj <- vector("list", 4)

for (i in 1:length(Pens)) {
  B1 <- Pens[[i]][, 1:16] %>%
    rename("PA" = "Bird.1") %>%
    mutate("Bird.ID" = Bird.IDs[[i]][1])
  B2 <- Pens[[i]][, c(1:15, 17)] %>%
    rename("PA" = "Bird.2") %>%
    mutate("Bird.ID" = Bird.IDs[[i]][2])
  B3 <- Pens[[i]][, c(1:15, 18)] %>%
    rename("PA" = "Bird.3") %>%
    mutate("Bird.ID" = Bird.IDs[[i]][3])
  B4 <- Pens[[i]][, c(1:15, 19)] %>%
    rename("PA" = "Bird.4") %>%
    mutate("Bird.ID" = Bird.IDs[[i]][4])
  B5 <- Pens[[i]][, c(1:15, 20)] %>%
    rename("PA" = "Bird.5") %>%
    mutate("Bird.ID" = Bird.IDs[[i]][5])
  Pens.Adj[[i]] <- rbind(B1, B2, B3, B4, B5)
  Keep <- c(which(Pens.Adj[[i]]$PA == "1"))
  Pens.Adj[[i]] <- Pens.Adj[[i]] %>% mutate("Feeding" = 0)
  Pens.Adj[[i]]$Feeding[Keep] <- as.character(Pens.Adj[[i]]$Bird.ID[Keep])
  rm("Keep")
  rm("B1", "B2", "B3", "B4", "B5")
}

for (i in 1:length(Pens.Adj)) {
  Pens.Adj[[i]]$Feeding[c(which(Pens.Adj[[i]]$Feeding == "0"))] <- NA
}

# Summation

PenSums <- vector("list", 4)
doubles <- vector("list", 4)

NAaFunc <- function(x) {
  unique(x[!is.na(x)])
}

for (i in 1:length(Pens.Adj)) {
  PenSums[[i]] <- as.data.frame(
    Pens.Adj[[i]] %>%
      group_by(Pen_Short, Treatment, Date.of.Photo, DOE, Absolute.Seconds, Mean.Amb) %>%
      summarise("Unique" = paste(NAaFunc(Feeding), collapse = ","))
  )
  doubles[[i]] <- c(grep(",", PenSums[[i]]$Unique))
}

# Assessing doubles

print(doubles) # quite many!

for (i in 1:length(PenSums)) {
  PenSums[[i]]$Unique[c(which(PenSums[[i]]$Unique == ""))] <- "0"
}

# Recombining

AllLong <- bind_rows(PenSums)

# Calculating average time spent feeding

Pres = AllLong[c(which(AllLong$Unique != "0")),]

Switch = c()

for (i in 2:nrow(Pres)){
  if (Pres$Absolute.Seconds[i]-1 == Pres$Absolute.Seconds[i-1]){
    if (Pres$Unique[i] == Pres$Unique[i-1]){
      Switch[i] = "Hold"
    } else {
      Switch[i] = "Trans"
    }
  } else {
    Switch[i] = "Trans"
  }
}

Switch[1] = "Hold"
Pres$Switch = Switch

Trans_points = c(as.numeric(which(Pres$Switch == "Trans")))
head(Trans_points)

Bout_length = c()

for (i in 2:length(Trans_points)){
  Bout_length[i] = Trans_points[i] - Trans_points[i-1]
}

Bout_length = Bout_length[-c(which(is.na(Bout_length)))]

mean(Bout_length)
getmode(Bout_length)
sd(Bout_length)
Bout_length

ggplot(subset(as.data.frame(Bout_length), Bout_length <=60), aes(Bout_length)) +
  geom_histogram(bins = 20, colour = "black", fill = "royalblue2") +
  theme_bw() + xlab("Duration of Feeding Bout (s)") +
  geom_vline(xintercept = 5, colour = "red", linetype = "longdash")

# A criteria of 4 seconds appears more reasonable than 3 seconds.

# Preparing for dominance assessment

Divide = split(AllLong, f = AllLong$Pen_Short)
for (i in 1:length(Divide)){
  Divide[[i]] = split(Divide[[i]], f = Divide[[i]]$Date.of.Photo)
}

for (k in 1:length(Divide)){
  for (j in 1:length(Divide[[k]])){
    Trans_Number_S1 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Number_S2 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Number_W1 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Number_W2 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Type_S1 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Type_S2 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Type_W1 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Type_W2 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_ID_S1 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_ID_S2 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_ID_W1 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_ID_W2 = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Number_D = vector("list", nrow(Divide[[k]][[j]]))
    Trans_ID_D = vector("list", nrow(Divide[[k]][[j]]))
    Trans_Type_D = vector("list", nrow(Divide[[k]][[j]]))

    for (i in 1:(nrow(Divide[[k]][[j]])-5)){
      if (Divide[[k]][[j]][i,"Unique"] != "0" &
          grepl(",", Divide[[k]][[j]][i, "Unique"]) == "FALSE" &
          Divide[[k]][[j]][i,"Unique"] != Divide[[k]][[j]][(i+1),"Unique"] &
          Divide[[k]][[j]][(i+1),"Unique"] != "0" &
          grepl(",", Divide[[k]][[j]][(i+1), "Unique"]) == "FALSE" &
          Divide[[k]][[j]][(i+1),"Unique"] == Divide[[k]][[j]][(i+2),"Unique"] &
          Divide[[k]][[j]][(i+1),"Unique"] == Divide[[k]][[j]][(i+3),"Unique"] &
          Divide[[k]][[j]][(i+1),"Unique"] == Divide[[k]][[j]][(i+4),"Unique"] &
          Divide[[k]][[j]][(i+1),"Unique"] == Divide[[k]][[j]][(i+5),"Unique"]
          ){
            Trans_Number_S1[[i]] <- c(i, i + 1, i + 2, i + 3, i + 4, i + 5)
            Trans_ID_S1[[i]] <- c("B1", "B2", "B2", "B2", "B2", "B2")
            Trans_Type_S1[[i]] <- c("Supplant", "Supplant", "Supplant", "Supplant", "Supplant", "Supplant")
      } else if (Divide[[k]][[j]][i,"Unique"] != "0" &
                grepl(",", Divide[[k]][[j]][i, "Unique"]) == "FALSE" &
                Divide[[k]][[j]][i,"Unique"] != Divide[[k]][[j]][(i+1),"Unique"] &
                grepl(",", Divide[[k]][[j]][(i+1), "Unique"]) == "TRUE" &
                Divide[[k]][[j]][i,"Unique"] != Divide[[k]][[j]][(i+2),"Unique"] &
                Divide[[k]][[j]][(i+2),"Unique"] != "0" &
                grepl(",", Divide[[k]][[j]][(i+2), "Unique"]) == "FALSE" &
                Divide[[k]][[j]][(i+2),"Unique"] == Divide[[k]][[j]][(i+3),"Unique"] &
                Divide[[k]][[j]][(i+2),"Unique"] == Divide[[k]][[j]][(i+4),"Unique"] &
                Divide[[k]][[j]][(i+2),"Unique"] == Divide[[k]][[j]][(i+5),"Unique"] &
                Divide[[k]][[j]][(i+2),"Unique"] == Divide[[k]][[j]][(i+6),"Unique"]
                ){
                  Trans_Number_S2[[i]] <- c(i, i + 1, i + 2, i + 3, i + 4, i + 5, i + 6)
                  Trans_ID_S2[[i]] <- c("B1", "GAP", "B2", "B2", "B2", "B2", "B2")
                  Trans_Type_S2[[i]] <- c("Supplant", "Supplant", "Supplant", "Supplant", "Supplant", "Supplant", "Supplant")
      } else if (Divide[[k]][[j]][i,"Unique"] != "0" &
               grepl(",", Divide[[k]][[j]][i, "Unique"]) == "FALSE" &
               Divide[[k]][[j]][(i+1),"Unique"] == "0" &
               Divide[[k]][[j]][i,"Unique"] != Divide[[k]][[j]][(i+2),"Unique"] &
               Divide[[k]][[j]][(i+2),"Unique"] != "0" &
               grepl(",", Divide[[k]][[j]][(i+2), "Unique"]) == "FALSE" &
               Divide[[k]][[j]][(i+2),"Unique"] == Divide[[k]][[j]][(i+3),"Unique"] &
               Divide[[k]][[j]][(i+2),"Unique"] == Divide[[k]][[j]][(i+4),"Unique"]){
                 Trans_Number_W1[[i]] <- c(i, i + 1, i + 2, i + 3, i + 4)
                 Trans_ID_W1[[i]] <- c("B1", "GAP", "B2", "B2", "B2")
                 Trans_Type_W1[[i]] <- c("Wait", "Wait", "Wait", "Wait", "Wait")
      } else if (Divide[[k]][[j]][i,"Unique"] != "0" &
               grepl(",", Divide[[k]][[j]][i, "Unique"]) == "FALSE" &
               Divide[[k]][[j]][(i+1),"Unique"] == "0" &
               Divide[[k]][[j]][(i+2),"Unique"] == "0" &
               Divide[[k]][[j]][i,"Unique"] != Divide[[k]][[j]][(i+3),"Unique"] &
               Divide[[k]][[j]][(i+3),"Unique"] != "0" &
               grepl(",", Divide[[k]][[j]][(i+3), "Unique"]) == "FALSE" &
               Divide[[k]][[j]][(i+3),"Unique"] == Divide[[k]][[j]][(i+4),"Unique"] &
               Divide[[k]][[j]][(i+3),"Unique"] == Divide[[k]][[j]][(i+5),"Unique"]
               ){
                 Trans_Number_W2[[i]] <- c(i, i + 1, i + 2, i + 3, i + 4, i + 5)
                 Trans_ID_W2[[i]] <- c("B1", "GAP", "GAP", "B2", "B2", "B2")
                 Trans_Type_W2[[i]] <- c("Wait", "Wait", "Wait", "Wait", "Wait", "Wait")
      } else if (Divide[[k]][[j]][i,"Unique"] != "0" &
               grepl(",", Divide[[k]][[j]][i, "Unique"]) == "FALSE" &
               Divide[[k]][[j]][i,"Unique"] != Divide[[k]][[j]][(i+1),"Unique"] &
               grepl(",", Divide[[k]][[j]][(i+1), "Unique"]) == "TRUE" &
               Divide[[k]][[j]][(i+1),"Unique"] == Divide[[k]][[j]][(i+2),"Unique"]){
                 Trans_Number_D[[i]] <- c(i, i + 1, i + 2)
                 Trans_ID_D[[i]] <- c("B1", "B2", "B2")
                 Trans_Type_D[[i]] <- c("Draw", "Draw", "Draw")
    }
    Trans_S1 <- data.frame(
      "Row" = c(unlist(Trans_Number_S1)), "Type" = c(unlist(Trans_ID_S1)),
      "Trans_Factor" = c(unlist(Trans_Type_S1))
    )
    Trans_S2 <- data.frame(
      "Row" = c(unlist(Trans_Number_S2)), "Type" = c(unlist(Trans_ID_S2)),
      "Trans_Factor" = c(unlist(Trans_Type_S2))
    )
    Trans_W1 <- data.frame(
      "Row" = c(unlist(Trans_Number_W1)), "Type" = c(unlist(Trans_ID_W1)),
      "Trans_Factor" = c(unlist(Trans_Type_W1))
    )
    Trans_W2 <- data.frame(
      "Row" = c(unlist(Trans_Number_W2)), "Type" = c(unlist(Trans_ID_W2)),
      "Trans_Factor" = c(unlist(Trans_Type_W2))
    )
    Trans_D <- data.frame(
      "Row" = c(unlist(Trans_Number_D)), "Type" = c(unlist(Trans_ID_D)),
      "Trans_Factor" = c(unlist(Trans_Type_D))
    )

    Divide[[k]][[j]]$Transitions <- "0"
    Divide[[k]][[j]]$Trans_Factor <- "0"
    Divide[[k]][[j]]$Transitions <- as.character(Divide[[k]][[j]]$Transitions)
    Divide[[k]][[j]]$Trans_Factor <- as.character(Divide[[k]][[j]]$Trans_Factor)

    for (i in 1:nrow(Trans_S1)) {
      Divide[[k]][[j]]$Transitions[Trans_S1$Row[i]] <- as.character(Trans_S1$Type[i])
      Divide[[k]][[j]]$Trans_Factor[Trans_S1$Row[i]] <- as.character(Trans_S1$Trans_Factor[i])
    }
    for (i in 1:nrow(Trans_S2)) {
      Divide[[k]][[j]]$Transitions[Trans_S2$Row[i]] <- as.character(Trans_S2$Type[i])
      Divide[[k]][[j]]$Trans_Factor[Trans_S2$Row[i]] <- as.character(Trans_S2$Trans_Factor[i])
    }
    for (i in 1:nrow(Trans_W1)) {
      Divide[[k]][[j]]$Transitions[Trans_W1$Row[i]] <- as.character(Trans_W1$Type[i])
      Divide[[k]][[j]]$Trans_Factor[Trans_W1$Row[i]] <- as.character(Trans_W1$Trans_Factor[i])
    }
    for (i in 1:nrow(Trans_W2)) {
      Divide[[k]][[j]]$Transitions[Trans_W2$Row[i]] <- as.character(Trans_W2$Type[i])
      Divide[[k]][[j]]$Trans_Factor[Trans_W2$Row[i]] <- as.character(Trans_W2$Trans_Factor[i])
    }
    for (i in 1:nrow(Trans_D)) {
      Divide[[k]][[j]]$Transitions[Trans_D$Row[i]] <- as.character(Trans_D$Type[i])
      Divide[[k]][[j]]$Trans_Factor[Trans_D$Row[i]] <- as.character(Trans_D$Trans_Factor[i])
    }
   }
  }
}

recombined <- vector("list", length(Divide))
for (i in 1:length(Divide)) {
  recombined[[i]] <- bind_rows(Divide[[i]])
}
Full <- as.data.frame(bind_rows(recombined))
Transitions <- subset(Full, Trans_Factor == "Supplant" | Trans_Factor == "Wait" | Trans_Factor == "GAP" | Trans_Factor == "Draw")
Rm_Gaps <- subset(Transitions, Transitions != "GAP")

Displaced <- c()
Displacer <- c()
Waitee <- c()
Waiter <- c()
Disp_keep <- c()
Wait_keep <- c()
Draw_1 <- c()
Draw_2 <- c()
Draw_keep <- c()

for (i in 1:nrow(Rm_Gaps)) {
  if (Rm_Gaps$Trans_Factor[i] == "Supplant" & Rm_Gaps$Transitions[i] == "B1"){
      Displaced[i] <- as.character(Rm_Gaps$Unique[i])
      Displacer[i] <- as.character(Rm_Gaps$Unique[i + 1])
      Disp_keep[i] <- i
  } else if (Rm_Gaps$Trans_Factor[i] == "Wait" & Rm_Gaps$Transitions[i] == "B1"){
      Waitee[i] <- as.character(Rm_Gaps$Unique[i])
      Waiter[i] <- as.character(Rm_Gaps$Unique[i + 1])
      Wait_keep[i] <- i
  } else if (Rm_Gaps$Trans_Factor[i] == "Draw" & Rm_Gaps$Transitions[i] == "B1"){
      Draw_1[i] <- as.character(Rm_Gaps$Unique[i])
      Draw_2[i] <- as.character(gsub(",", "", gsub(Rm_Gaps$Unique[i], "", Rm_Gaps$Unique[i+1])))
      Draw_keep[i] <- i
  }
}

Displacement <- rbind(data.frame("Pen_Short" = Rm_Gaps[Disp_keep,1], "Treatment" = Rm_Gaps[Disp_keep,2],
                                 "Date.of.Photo" = Rm_Gaps[Disp_keep,3], "DOE" = Rm_Gaps[Disp_keep,4],
                                 "Absolute.Seconds" = Rm_Gaps[Disp_keep,5], "Mean.Amb" = Rm_Gaps[Disp_keep,6],
                                 "Dominant" = Displacer, "Subordinate" = Displaced, "Type" = "Supplant"),
                      data.frame("Pen_Short" = Rm_Gaps[Wait_keep,1], "Treatment" = Rm_Gaps[Wait_keep,2],
                                 "Date.of.Photo" = Rm_Gaps[Wait_keep,3], "DOE" = Rm_Gaps[Wait_keep,4],
                                 "Absolute.Seconds" = Rm_Gaps[Wait_keep,5], "Mean.Amb" = Rm_Gaps[Wait_keep,6],
                                 "Dominant" = Waitee, "Subordinate" = Waiter, "Type" = "Wait"),
                      data.frame("Pen_Short" = Rm_Gaps[Draw_keep,1], "Treatment" = Rm_Gaps[Draw_keep,2],
                                 "Date.of.Photo" = Rm_Gaps[Draw_keep,3], "DOE" = Rm_Gaps[Draw_keep,4],
                                 "Absolute.Seconds" = Rm_Gaps[Draw_keep,5], "Mean.Amb" = Rm_Gaps[Draw_keep,6],
                                 "Dominant" = Draw_1, "Subordinate" = Draw_2, "Type" = "Draw")
                      )

Dominance = na.omit(Displacement)
testing = Dominance %>% arrange(Pen_Short,Date.of.Photo)
print(testing)

#write.csv(testing, "/home/joshk/git_repositories/BCCH_Dominance/Data/Dominance_Dat5.csv")

testing = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Dominance_Dat5.csv")
summed <- as.data.frame(testing %>% group_by(Pen_Short, Dominant, Subordinate) %>% summarise("Count" = n()))
print(summed)

# Formatting data sets for calculation of Elo Ratings.

testing$Draw = "FALSE"
testing$Draw = as.character(testing$Draw)
testing$Draw[c(which(testing$Type == "Draw"))] = "TRUE"
testing$Draw = as.logical(testing$Draw)

EloDat_min <- testing %>%
  dplyr::select(
    "Date" = Date.of.Photo, "Time" = Absolute.Seconds,
    "Winner" = Dominant, "Loser" = Subordinate, Type,
    "Draw" = Draw, "Pen" = Pen_Short
  ) %>%
  arrange(Pen, Date, Time)

head(EloDat_min)

# Adjusting Date Format

New_date <- c()

for (i in 1:nrow(EloDat_min)) {
  New_date[i] <- paste(substr(EloDat_min$Date[i], 0, 4),
    substr(EloDat_min$Date[i], 5, 6),
    substr(EloDat_min$Date[i], 7, 8),
    sep = "-"
  )
}

EloDat_min$Date <- as.character(New_date)

# Splitting by Pen and Naming Lists

EloDat_min$Date <- as.Date(EloDat_min$Date)

# Double-checking that there are no oddities.

unique(EloDat_min$Winner)
unique(EloDat_min$Loser)

# Yes, one oddity (OOAABLO). Removing. 

EloDat_min = EloDat_min[-c(which(EloDat_min$Loser == "OOAABlO")),]

EloDat_list <- split(EloDat_min, f = EloDat_min$Pen)
names(EloDat_list) <- c(as.character(unique(EloDat_min$Pen)))

# Quantifying number of interactions experienced by an individual.

IperBird = function(x){
  IDs = unique(c(as.character(x$Winner), as.character(x$Loser)))
  Ints = c()
  for (i in 1:length(IDs)){
    Ints[i] = length(which(x$Winner == IDs[i])) +
      length(which(x$Loser == IDs[i]))
  }
  return(data.frame("Bird.ID" = IDs, "nInt" = Ints))
}

bind_rows(lapply(EloDat_list, IperBird)) %>%
  summarise(Mean = mean(nInt), SD = sd(nInt), Total = sum(nInt))

# Average number of interactions per bird = 50.300 +- 23.414. Reasonable sample size. 

# Pulling ratings according to raw, unrandomised Elo scores.

Elo_Check <- vector("list", 4)

for (i in 1:length(Elo_Check)) {
  Elo_Check[[i]] <- seqcheck(
    winner = as.factor(EloDat_list[[i]]$Winner),
    loser = as.factor(EloDat_list[[i]]$Loser),
    Date = EloDat_list[[i]]$Date,
    draw = EloDat_list[[i]]$Draw
  )
}

Elo_Check

## Adding presence data.

{
  Years <- c(rep("2018", (7 + 31 + 19)))
  Months <- c(rep("04", 7), rep("05", 31), rep("06", 19))
  Days <- c(c(24:30), c(1:31), c(1:19))
  Days <- str_pad(Days, 2, side = "left", pad = "0")
  Possible_Dates <- paste(Years, Months, Days, sep = "-")
  Possible_Dates_BP <- vector("list", 4)

  min_Date <- c()
  max_Date <- c()
  for (i in 1:length(Possible_Dates_BP)) {
    min_Date[i] <- as.character(range(EloDat_list[[i]]$Date)[1])
    max_Date[i] <- as.character(range(EloDat_list[[i]]$Date)[2])
    Possible_Dates_BP[[i]] <- Possible_Dates[c(which(Possible_Dates == min_Date[i]):
    which(Possible_Dates == max_Date[i]))]
  }
}

{
  Pen_NW <- data.frame(
    "Date" = as.character(Possible_Dates_BP[[1]]), "BdAO" = "1", "ABlBl" = "1",
    "BdABd" = "1", "AYY" = "1", "YAO" = "1"
  )
  Pen_NE <- data.frame(
    "Date" = as.character(Possible_Dates_BP[[2]]), "ABlO" = "1", "OOA" = "1", "YAO2" = "1",
    "ABdBd" = "1", "BdABl" = "1"
  )
  Pen_SW <- data.frame(
    "Date" = as.character(Possible_Dates_BP[[3]]), "RAR" = "1", "BdABl2" = "1",
    "ARO" = "1", "ABlR" = "1", "YAR" = "1"
  )
  Pen_SE <- data.frame(
    "Date" = as.character(Possible_Dates_BP[[4]]), "BlAR" = "1", "BdAO2" = "1",
    "AOR" = "1", "AYBd" = "1", "YAY" = "1"
  )
}

Presence <- vector("list", 4)

{
  Presence[[1]] <- Pen_NW
  Presence[[1]]$Date <- as.Date(as.character(Presence[[1]]$Date))
  Presence[[2]] <- Pen_NE
  Presence[[2]]$Date <- as.Date(as.character(Presence[[2]]$Date))
  Presence[[3]] <- Pen_SW
  Presence[[3]]$Date <- as.Date(as.character(Presence[[3]]$Date))
  Presence[[4]] <- Pen_SE
  Presence[[4]]$Date <- as.Date(as.character(Presence[[4]]$Date))
}

for (j in 1:4) {
  for (i in 2:6) {
    Presence[[j]][, i] <- as.numeric(Presence[[j]][, i])
  }
}

# Testing hierarchies per flight enclosure with raw Ela ratings.

Elo_Results <- vector("list", 4)
Elo_Stab <- vector("list", 4)
Elo_Vals <- vector("list", 4)
for (i in 1:length(Elo_Vals)) {
  Elo_Vals[[i]] <- vector("list", length(unique(EloDat_list[[i]]$Date)))
}

for (i in 1:length(Elo_Results)) {
  Elo_Results[[i]] <- elo.seq(
    winner = EloDat_list[[i]]$Winner,
    loser = EloDat_list[[i]]$Loser,
    Date = EloDat_list[[i]]$Date,
    draw = EloDat_list[[i]]$Draw,
    presence = Presence[[i]]
  )
  Elo_Stab[[i]] <- stab_elo(Elo_Results[[i]])
}

Pens <- c("NW", "NE", "SW", "SE")
for (j in 1:length(Elo_Results)) {
  for (i in 1:length(unique(EloDat_list[[j]]$Date))) {
    Elo_Vals[[j]][[i]] <- as.data.frame(extract_elo(Elo_Results[[j]],
      extractdate = unique(EloDat_list[[j]]$Date)[i],
      IDs = c(unique(as.character(EloDat_list[[j]]$Winner)))
    ))
    Elo_Vals[[j]][[i]]$Date <- unique(EloDat_list[[j]]$Date)[i]
    Elo_Vals[[j]][[i]]$Bird.ID <- rownames(Elo_Vals[[j]][[i]])
    colnames(Elo_Vals[[j]][[i]])[1] <- "Elo"
  }
  Elo_Vals[[j]] <- bind_rows(Elo_Vals[[j]])
  Elo_Vals[[j]]$Pen <- Pens[[j]]
}

Elo_Vals <- bind_rows(Elo_Vals)

# Assigning ranks by raw, unrandomised Elo scores

Elo_Rank <- vector('list', 4)
for (j in 1:length(Elo_Results)) {
    Elo_Rank[[j]] <- as.data.frame(extract_elo(Elo_Results[[j]]))
    Elo_Rank[[j]]$Bird.ID <- c(row.names(as.data.frame(extract_elo(Elo_Results[[j]]))))
    Elo_Rank[[j]]$Pen <- Pens[j]
    Elo_Rank[[j]]$E.Rank <- as.factor(as.character(c(1:5)))
}
Elo_dat <- bind_rows(Elo_Rank)
colnames(Elo_dat)[1] = "Elo_Score"
Elo_dat = Elo_dat %>% dplyr::select(Bird.ID, Pen, Elo_Score, E.Rank)

### Using custom randomisation function to calculate randomised Elo Ratings from randomly sampled
# days of observation.

R.Elos = vector('list', 4)
Stab = vector('list', 4)

for (i in 1:4){
  output_list = vector('list', 2)
  sampled_Elo = Rand_Elo(data = EloDat_list[[i]], presence = Presence[[i]], kvals = list(Supplant = 50, Wait = 25, Draw = 0), sample_size = 20, niter = 1000)
  for_Stab = as.data.frame(sampled_Elo %>% gather(IDs, Elo_Score, colnames(sampled_Elo)[1]:colnames(sampled_Elo)[5]) %>%
    arrange(Sample, desc(Elo_Score)) %>% mutate("Ranks" = rep(c(1:5),1000), "Status" = rep(c(1,1,2,2,2), 1000)))
  summed_Elo = sampled_Elo %>% gather(IDs, Elo_Score, colnames(sampled_Elo)[1]:colnames(sampled_Elo)[5]) %>%
    group_by(IDs) %>% dplyr::summarise(Mean_Elo = mean(Elo_Score)) %>%
    arrange(desc(Mean_Elo)) %>% mutate("Rank" = c(1:nrow(.)))
  sampled_plot = sampled_Elo %>% gather(IDs, Elo_Score, colnames(sampled_Elo)[1]:colnames(sampled_Elo)[5]) %>%
    group_by(IDs) %>% dplyr::summarise(Mean_Elo = mean(Elo_Score),
                                       SE = sd(Elo_Score)/sqrt(n()),
                                       Lower.CI = Mean_Elo - 1.96*SE,
                                       Upper.CI = Mean_Elo + 1.96*SE) %>%
    ggplot(aes(x = IDs, y = Mean_Elo, fill = IDs)) +
      geom_point(pch = 21, colour = "black", size = 4) +
      geom_errorbar(mapping = aes(x = IDs, ymin = Lower.CI, ymax = Upper.CI),
        colour = "black", size = 1, width = 0.3) + theme_bw()
    output_list[[1]] = summed_Elo
    output_list[[2]] = sampled_plot
    R.Elos[[i]] = output_list
    Stab[[i]] = for_Stab
}

# Quick plot 

grid.arrange(R.Elos[[1]][[2]], R.Elos[[2]][[2]],
  R.Elos[[3]][[2]], R.Elos[[4]][[2]], nrow = 2)

# Testing that customised function is functioning correctly by plotting
# randomised Elo ratings by raw, unrandomised Elo ratings.

Random_Elo = bind_rows(as.data.frame(R.Elos[[1]][[1]] %>%
    mutate("Pen" = "NW")), as.data.frame(R.Elos[[2]][[1]] %>%
    mutate("Pen" = "NE")), as.data.frame(R.Elos[[3]][[1]] %>%
    mutate("Pen" = "SW")), as.data.frame(R.Elos[[4]][[1]] %>%
    mutate("Pen" = "SE"))
    ) %>% rename("R_Elo_Rank" = Rank, "Bird.ID" = IDs,
    "R_Elo_Value" = Mean_Elo) %>% 
    mutate(R_Elo_Rank = as.integer(R_Elo_Rank))

Raw_Elo = Elo_dat %>% rename("NR_Elo_Rank" = E.Rank,
    "NR_Elo_Value" = Elo_Score) %>% 
    mutate(NR_Elo_Rank = as.integer(NR_Elo_Rank))

Validation = left_join(Random_Elo, Raw_Elo, by = c("Bird.ID", "Pen"))

ggplot(Validation, aes(x = NR_Elo_Value, y = R_Elo_Value)) + 
    geom_point(size = 3, pch = 21, fill = "seagreen", 
    colour = "black", alpha = 0.3) + 
    geom_smooth(method = "lm", size = 1, colour = "black") + 
    theme_bw() + xlab("Non-random Elo Score") + 
    ylab("Random Elo Score")

summary(lm(R_Elo_Value ~ NR_Elo_Value, data = Validation))

ggplot(Validation, aes(x = NR_Elo_Rank, y = R_Elo_Rank)) + 
    geom_point(size = 3, pch = 21, fill = "seagreen", 
    colour = "black", alpha = 0.3) + 
    geom_smooth(method = "lm", size = 1, colour = "black") + 
    theme_bw() + xlab("Non-random Elo Rank") + 
    ylab("Random Elo Rank")

summary(lm(R_Elo_Rank ~ NR_Elo_Rank, data = Validation))

# Significantly predictive in both cases, however, with some variation.
# Proceeding with randomised Elo ratings according to Sanchez et al. (2018).

# Producing plots for supplemental information.

Plot.Data = vector('list', 4)

for (i in 1:4){
  sampled_Elo = Rand_Elo(data = EloDat_list[[i]], presence = Presence[[i]], kvals = list(Supplant = 50, Wait = 25, Draw = 0), sample_size = 20, niter = 1000)
  plot_dat = as.data.frame(sampled_Elo %>% gather(IDs, Elo_Score, colnames(sampled_Elo)[1]:colnames(sampled_Elo)[5]) %>%
    group_by(IDs) %>% dplyr::summarise(Mean_Elo = mean(Elo_Score),
                                       SE = sd(Elo_Score)/sqrt(n()),
                                       Lower.CI = Mean_Elo - 2.58*SE,
                                       Upper.CI = Mean_Elo + 2.58*SE))
  Plot.Data[[i]] = plot_dat
}

for (i in 1:4){
  Plot.Data[[i]] = Plot.Data[[i]][c(order(Plot.Data[[i]][,2], decreasing = T)),]
  Plot.Data[[i]]$Status = c("Dominant", "Dominant", rep("Subordinate", 3))
  Plot.Data[[i]]$Sex = c(rep("Female",5))
  Plot.Data[[i]]$Sex.S = c(rep("\\u2640",5))
}

Plot.Data[[1]]$Sex[1:4] = "Male"
Plot.Data[[2]]$Sex[2:3] = "Male"
Plot.Data[[3]]$Sex[3:4] = "Male"
Plot.Data[[4]]$Sex[c(1,5)] = "Male"
Plot.Data[[1]]$Sex.S[1:4] = "\\u2642"
Plot.Data[[2]]$Sex.S[2:3] = "\\u2642"
Plot.Data[[3]]$Sex.S[3:4] = "\\u2642"
Plot.Data[[4]]$Sex.S[c(1,5)] = "\\u2642"

Plot.Data[[1]]$Bird.Lab = sample(1:5, 5)
Plot.Data[[2]]$Bird.Lab = sample(1:5, 5)
Plot.Data[[3]]$Bird.Lab = sample(1:5, 5)
Plot.Data[[4]]$Bird.Lab = sample(1:5, 5)

for (j in 1:4){
  Plot.Data[[j]]$Status = factor(Plot.Data[[j]]$Status)
  Plot.Data[[j]]$Sex = factor(Plot.Data[[j]]$Sex)
  Plot.Data[[j]]$Bird.Lab = as.character(Plot.Data[[j]]$Bird.Lab)
  Plot.Data[[j]]$Group = paste(Plot.Data[[j]]$Status, Plot.Data[[j]]$Sex, sep = " ")
}

for (j in 1:4){
  Plot.Data[[j]]$Group = factor(Plot.Data[[j]]$Group, levels = c("Dominant Female", "Dominant Male",
   "Subordinate Female", "Subordinate Male"))
}

Plot.Data[[2]]

Final.Plots = vector('list', 4)
xpos = c(-Inf, Inf, -Inf, Inf)
ypos = c(Inf, Inf, Inf, Inf)
hpos = c(-0.75, 1.75, -0.75, 1.75)
vpos = c(1.75,1.75,1.75,1.75)

min_break = c()
max_break = c()
interv = c()

for (i in 1:4){
  min_break[i] = min(Plot.Data[[i]]$Lower.CI) - 5
  max_break[i] = max(Plot.Data[[i]]$Upper.CI) + 5
  interv[i] = ((max(Plot.Data[[i]]$Upper.CI) + 5) - (min(Plot.Data[[i]]$Lower.CI) - 5))/4
}

min_break = c(960,960,950,940)
max_break = c(1030,1030, 1070, 1050)
intervs = vector('list', 4)
intervs[[1]] = seq(960, 1035, by = 15)
intervs[[2]] = seq(960, 1035, by = 15)
intervs[[3]] = seq(950, 1075, by = 25)
intervs[[4]] = seq(935, 1060, by = 25)

for (i in 1:4){
  Final.Plots[[i]] = Plot.Data[[i]] %>%
    ggplot(aes(x = Bird.Lab, y = Mean_Elo, fill = Status, shape = Sex)) +
      geom_point(size = 3, colour = "black") +
      geom_errorbar(mapping = aes(x = Bird.Lab, ymin = Lower.CI, ymax = Upper.CI),
        width = 0.3, size = 0.75, colour = "black") +
      theme_bw() +
      scale_fill_viridis_d() +
      scale_shape_manual(values = c(21,24)) +
      theme(axis.ticks = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), ) +
      scale_x_discrete(name ="") + guides(fill = "Social Ranks") +
      scale_y_continuous(breaks = intervs[[i]]) + theme(axis.text.x = element_text(size = 20, family = "Noto Sans"))
      annotate("text", x = xpos[i], y = ypos[i], hjust = hpos[i], vjust = vpos[i], size = 8, label = i)
}

Legend_Plot = Plot.Data[[2]] %>%
  ggplot(aes(x = Bird.Lab, y = Mean_Elo, fill = Group, shape = Group)) +
    geom_point(size = 3, colour = "black") +
    geom_errorbar(mapping = aes(x = Bird.Lab, ymin = Lower.CI, ymax = Upper.CI),
      width = 0.3, size = 0.75, colour = "black") +
    theme_bw() +
    scale_fill_manual(name = "Social Rank and Sex",
                   labels = c("Dominant Female", "Dominant Male", "Subordinate Female", "Subordinate Male"),
                   values = c(viridis::viridis(n = 2)[1], viridis::viridis(n = 2)[1],
                    viridis::viridis(n = 2)[2], viridis::viridis(n = 2)[2])) +
                   my.theme +
    scale_shape_manual(name = "Social Rank and Sex",
                   labels = c("Dominant Female", "Dominant Male", "Subordinate Female", "Subordinate Male"),
                   values = c(21,24,21,24)) +
    theme(axis.ticks = element_blank(), axis.title.y = element_blank(), ) +
    scale_x_discrete(name ="")

legend = gtable_filter(ggplotGrob(Legend_Plot), "guide-box")
yaxis_lab = grid.text("Randomized Elo Score", gp=gpar(fontfamily = "Noto Sans", cex = 1.75), rot = 90)
bottom_lab = grid.text("Bird Identity", gp=gpar(fontfamily = "Noto Sans", cex = 1.75))

Status.Plots = grid.arrange(arrangeGrob(Final.Plots[[1]] + theme(legend.position="none"),
  Final.Plots[[2]] + theme(legend.position="none"),
  Final.Plots[[3]] + theme(legend.position="none"),
  Final.Plots[[4]] + theme(legend.position="none"),
  nrow = 2,
  left = yaxis_lab,
  bottom = bottom_lab), legend,
  widths=unit.c(unit(1, "npc") - legend$width, legend$width),
  nrow = 1)

showtext_opts(dpi = 2000)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Status_Plots.jpeg", Status.Plots, dpi = 2000,
  width = 9.67, height = 6.67, units = "in")

## Assessing stability of ranks and status

Pens = c("NW", "NE", "SW", "SE")
Stab.per.pen = vector('list',4)

for (j in 1:4){
  Stab.per.bird = vector('list', 5)
  for (i in 1:5){
    bird = unique(Stab[[j]]$IDs)[i]
    mode_rank = getmode(subset(Stab[[j]], IDs == bird)$Ranks)
    mode_status = getmode(subset(Stab[[j]], IDs == bird)$Status)
    stab_rank = length(which(subset(Stab[[j]], IDs == bird)$Ranks == mode_rank))/1000
    stab_status = length(which(subset(Stab[[j]], IDs == bird)$Status == mode_status))/1000
    Stab.per.bird[[i]] = data.frame("Pen" = Pens[j], "Bird.ID" = bird, "Rank.Mode" = mode_rank, "Status.Mode" = mode_status,
    "Rank.Stability" = stab_rank, "Status.Stability" = stab_status)
  }
  Stab.per.pen[[j]] = bind_rows(Stab.per.bird)
}

Stability.all = bind_rows(Stab.per.pen)
#write.csv(Stability.all, "/home/joshk/git_repositories/BCCH_Dominance/Data/Stability_data.csv")

# Stability data can be read in here

Stability.all = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Stability_data.csv")

mean(Stability.all$Rank.Stability)
sd(Stability.all$Rank.Stability)

mean(1-Stability.all$Status.Stability)
sd(Stability.all$Status.Stability)

summary(aov(Stability.all$Status.Stability ~ Stability.all$Pen))

t.test(Stability.all$Rank.Stability, Stability.all$Status.Stability)

# Saving final ranks

Validation = Validation %>% rename("Elo_Rank" = NR_Elo_Rank, "RS.Rank" = R_Elo_Rank)
Validation$E.Status = "1"
Validation$E.Status[c(which(Validation$Elo_Rank >=3))] = "2"
Validation$RSE.Status = "1"
Validation$RSE.Status[c(which(Validation$RS.Rank >=3))] = "2"
Validation = Validation %>% mutate(E.Status = factor(E.Status),
    RSE.Status = factor(RSE.Status))

#write.csv(Validation, "/home/joshk/git_repositories/BCCH_Dominance/Data/Rank_Summaries.csv")

######################################################################################
# Final rank data saved. Fusing with remainder of data (eye temperature and feeding).#
######################################################################################

E_join = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Rank_Summaries.csv")

All_Ranks = as.data.frame(E_join %>% dplyr::select(Bird.ID, Pen, Elo_Rank,
    "RSE_Rank" = RS.Rank, E.Status, RSE.Status))

# Expanding and applying to feeding data frame to assess the effects of 
# social status on access to food.

{
  All_IDs <- read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/All Combined Thermal Data.csv")
  All_IDs$Maximum.Eye.Temp <- as.numeric(as.character(All_IDs$Maximum.Eye.Temp))
  All_IDs$Amb.Temp <- as.numeric(as.character(All_IDs$Amb.Temp))

  ResFrame <- data.frame(
    "Row" = 1:nrow(All_IDs),
    "Residuals" = residuals(lm(Maximum.Eye.Temp ~ Amb.Temp,
                               data = All_IDs, na.action = na.exclude
    ))
  )

  to_remove <- c(
    c(which(All_IDs$Maximum.Eye.Temp <=
              mean(na.omit(All_IDs$Maximum.Eye.Temp)) -
              4 * sd(na.omit(All_IDs$Maximum.Eye.Temp)))),
    c(which(ResFrame$Residuals <=
              mean(na.omit(ResFrame$Residuals)) -
              4 * sd(na.omit(ResFrame$Residuals)))),
    c(which(ResFrame$Residuals >=
              mean(na.omit(ResFrame$Residuals)) +
              4 * sd(na.omit(ResFrame$Residuals))))
  )

  length(to_remove)

  # Removing 10 extreme and unlikely readings

  All_IDs = All_IDs[-to_remove, ]
  
  All_IDs$Treatment <- ordered(All_IDs$Treatment, levels = c("Rest", "Stress"))
  All_IDs$Bird.ID <- factor(All_IDs$Bird.ID)
  All_IDs$Date.of.Photo <- factor(All_IDs$Date.of.Photo)
  All_IDs$Pen_Short <- factor(All_IDs$Pen_Short)

  ## Adding bird capture locales

  Clean <- All_IDs

  Capture.Locale <- c()
  Locale.Type <- c()
  Sex <- c()

  for (i in 1:nrow(Clean)) {
    if (is.na(Clean$Bird.ID[i])) {
      Capture.Locale[i] = NA
      Locale.Type[i] = NA
      Sex[i] = NA
    } else if (Clean$Bird.ID[i] == "BdABd") {
      Capture.Locale[i] = "Corwhin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "BdAO1") {
      Capture.Locale[i] = "Erin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "ABlBl") {
      Capture.Locale[i] = "Ruthven.Park"
      Locale.Type[i] = "Urban"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "AYY") {
      Capture.Locale[i] = "Brantford"
      Locale.Type[i] = "Urban"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "YAO1") {
      Capture.Locale[i] = "Guelph"
      Locale.Type[i] = "Urban"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "OOA") {
      Capture.Locale[i] = "Ruthven.Park"
      Locale.Type[i] = "Rural"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "BdABl1") {
      Capture.Locale[i] = "Ruthven.Park"
      Locale.Type[i] = "Rural"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "YAO2") {
      Capture.Locale[i] = "Corwhin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "ABlO") {
      Capture.Locale[i] = "Brantford"
      Locale.Type[i] = "Urban"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "ABdBd") {
      Capture.Locale[i] = "Guelph"
      Locale.Type[i] = "Urban"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "YAR") {
      Capture.Locale[i] = "Cambridge"
      Locale.Type[i] = "Urban"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "RAR") {
      Capture.Locale[i] = "Corwhin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "ABlR") {
      Capture.Locale[i] = "Brantford"
      Locale.Type[i] = "Urban"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "ARO") {
      Capture.Locale[i] = "Corwhin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "BdABl2") {
      Capture.Locale[i] = "Corwhin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "BlAR") {
      Capture.Locale[i] = "Corwhin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "BdAO2") {
      Capture.Locale[i] = "Guelph"
      Locale.Type[i] = "Urban"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "AOR") {
      Capture.Locale[i] = "Guelph"
      Locale.Type[i] = "Urban"
      Sex[i] = "Female"
    } else if (Clean$Bird.ID[i] == "YAY") {
      Capture.Locale[i] = "Guelph"
      Locale.Type[i] = "Urban"
      Sex[i] = "Male"
    } else if (Clean$Bird.ID[i] == "AYBd") {
      Capture.Locale[i] = "Erin"
      Locale.Type[i] = "Rural"
      Sex[i] = "Female"
    } else {
      Capture.Locale[i] = NA
      Locale.Type[i] = NA
    }
  }

  Clean$Capture.Locale <- factor(Capture.Locale)
  Clean$Locale.Type <- factor(Locale.Type)
  Clean$Sex <- factor(Sex)

  # Removing those with poor image quality

  Clean$Maximum.Eye.Temp[c(which(Clean$Eye.Quality != "Good"))] <- NA

  Feeding.Dat <- read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Digital Bird With RT.csv")

  Feeding <- Feeding.Dat %>% dplyr::select(Date.of.Photo,
                                           Pen = Pen_Full, Treatment,
                                           DOE, Absolute.Seconds, Amb.Temp = Amb_Temp,
                                           Maximum.Eye.Temp,
                                           Time.of.Photo,
                                           Bird.1, Bird.2, Bird.3, Bird.4, Bird.5
  )

  F_by_HH <- as.data.frame(Feeding %>% group_by(Date.of.Photo, Pen, Treatment) %>%
                             summarize(
                               Bird.1.count = stri_count_fixed(paste(Bird.1, collapse = ""), "01"),
                               Bird.1 = sum(Bird.1),
                               Bird.2.count = stri_count_fixed(paste(Bird.2, collapse = ""), "01"),
                               Bird.2 = sum(Bird.2),
                               Bird.3.count = stri_count_fixed(paste(Bird.3, collapse = ""), "01"),
                               Bird.3 = sum(Bird.3),
                               Bird.4.count = stri_count_fixed(paste(Bird.4, collapse = ""), "01"),
                               Bird.4 = sum(Bird.4),
                               Bird.5.count = stri_count_fixed(paste(Bird.5, collapse = ""), "01"),
                               Bird.5 = sum(Bird.5),
                               All.Birds = sum(Bird.1, Bird.2, Bird.3, Bird.4, Bird.5),
                               Mean.Temp = mean(Amb.Temp, na.rm = T),
                               Time.of.Day = min(Time.of.Photo, na.rm = T)
                             ))

  F_by_HH$Pen_Short <- gsub("[[:space:]].*", "", F_by_HH$Pen)
  F_by_HH$Pen_Short <- factor(F_by_HH$Pen_Short)

  # Converting to long form
  {
    Bird.1 <- F_by_HH[, c(1:5, 15:17)]
    Cols <- colnames(Bird.1)
    Cols <- c(Cols[1:3], "Count", "TSF", Cols[6:8])
    colnames(Bird.1) <- Cols
    Bird.2 <- F_by_HH[, c(1:3, 6, 7, 15:17)]
    colnames(Bird.2) <- Cols
    Bird.3 <- F_by_HH[, c(1:3, 8, 9, 15:17)]
    colnames(Bird.3) <- Cols
    Bird.4 <- F_by_HH[, c(1:3, 10, 11, 15:17)]
    colnames(Bird.4) <- Cols
    Bird.5 <- F_by_HH[, c(1:3, 12, 13, 15:17)]
    colnames(Bird.5) <- Cols
    All_Birds <- F_by_HH[, c(1:3, 14, 14, 15:17)]
    colnames(All_Birds) <- Cols
  }
  Long_FbyHH <- rbind(Bird.1, Bird.2, Bird.3, Bird.4, Bird.5, All_Birds)

  Bird.IDs <- c(
    rep("Bird.1", nrow(F_by_HH)),
    rep("Bird.2", nrow(F_by_HH)),
    rep("Bird.3", nrow(F_by_HH)),
    rep("Bird.4", nrow(F_by_HH)),
    rep("Bird.5", nrow(F_by_HH)),
    rep("All", nrow(F_by_HH))
  )

  Long_FbyHH$Bird.ID <- Bird.IDs
  colnames(Long_FbyHH)[4] <- "Feeding.Bouts"
  colnames(Long_FbyHH)[5] <- "Feeding.Time"

  Long_FbyHH$Bird.ID <- factor(Long_FbyHH$Bird.ID)
  Long_FbyHH$Treatment <- ordered(Long_FbyHH$Treatment)
  Long_FbyHH$Date.of.Photo <- factor(Long_FbyHH$Date.of.Photo)
  Long_FbyHH$Feeding.Time <- as.integer(Long_FbyHH$Feeding.Time)

  Capture.Locale <- c()
  Locale.Type <- c()
  Bands <- c()
  Sex <- c()

  for (i in 1:nrow(Long_FbyHH)) {
    if (is.na(Long_FbyHH$Pen_Short[i])) {
      Capture.Locale[i] <- NA
      Locale.Type[i] <- NA
      Bands[i] <- NA
      Sex[i] <- NA
    } else if (Long_FbyHH$Pen_Short[i] == "NW") {
      if (Long_FbyHH$Bird.ID[i] == "Bird.1") {
        Capture.Locale[i] <- "Corwhin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "BdABd"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.2") {
        Capture.Locale[i] <- "Erin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "BdAO1"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.3") {
        Capture.Locale[i] <- "Ruthven.Park"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "ABlBl"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.4") {
        Capture.Locale[i] <- "Brantford"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "AYY"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.5") {
        Capture.Locale[i] <- "Guelph"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "YAO1"
        Sex[i] <- "Female"
      } else {
        Capture.Locale[i] <- NA
        Locale.Type[i] <- NA
        Bands[i] <- NA
        Sex[i] <- NA
      }
    } else if (Long_FbyHH$Pen_Short[i] == "NE") {
      if (Long_FbyHH$Bird.ID[i] == "Bird.1") {
        Capture.Locale[i] <- "Ruthven.Park"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "OOA"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.2") {
        Capture.Locale[i] <- "Ruthven.Park"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "BdABl1"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.3") {
        Capture.Locale[i] <- "Corwhin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "YAO2"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.4") {
        Capture.Locale[i] <- "Brantford"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "ABlO"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.5") {
        Capture.Locale[i] <- "Guelph"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "ABdBd"
        Sex[i] <- "Female"
      } else {
        Capture.Locale[i] <- NA
        Locale.Type[i] <- NA
        Sex[i] <- NA
      }
    } else if (Long_FbyHH$Pen_Short[i] == "SW") {
      if (Long_FbyHH$Bird.ID[i] == "Bird.1") {
        Capture.Locale[i] <- "Cambridge"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "YAR"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.2") {
        Capture.Locale[i] <- "Corwhin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "RAR"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.3") {
        Capture.Locale[i] <- "Brantford"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "ABlR"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.4") {
        Capture.Locale[i] <- "Corwhin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "ARO"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.5") {
        Capture.Locale[i] <- "Corwhin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "BdABl2"
        Sex[i] <- "Male"
      } else {
        Capture.Locale[i] <- NA
        Locale.Type[i] <- NA
        Bands[i] <- NA
        Sex[i] <- NA
      }
    } else if (Long_FbyHH$Pen_Short[i] == "SE") {
      if (Long_FbyHH$Bird.ID[i] == "Bird.1") {
        Capture.Locale[i] <- "Corwhin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "BlAR"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.2") {
        Capture.Locale[i] <- "Guelph"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "BdAO2"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.3") {
        Capture.Locale[i] <- "Guelph"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "AOR"
        Sex[i] <- "Female"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.4") {
        Capture.Locale[i] <- "Guelph"
        Locale.Type[i] <- "Urban"
        Bands[i] <- "YAY"
        Sex[i] <- "Male"
      } else if (Long_FbyHH$Bird.ID[i] == "Bird.5") {
        Capture.Locale[i] <- "Erin"
        Locale.Type[i] <- "Rural"
        Bands[i] <- "AYBd"
        Sex[i] <- "Female"
      } else {
        Capture.Locale[i] <- NA
        Locale.Type[i] <- NA
        Bands[i] <- NA
      }
    } else {
      Capture.Locale[i] <- NA
      Locale.Type[i] <- NA
      Bands[i] <- NA
    }
  }

  Long_FbyHH$Capture.Locale <- Capture.Locale
  Long_FbyHH$Locale.Type <- Locale.Type
  Long_FbyHH$Bands <- Bands
  Long_FbyHH$Sex <- Sex
  Long_FbyHH$Capture.Locale <- factor(Long_FbyHH$Capture.Locale)
  Long_FbyHH$Locale.Type <- ordered(Long_FbyHH$Locale.Type)
  Long_FbyHH$Bands <- factor(Long_FbyHH$Bands)
  Long_FbyHH$Sex <- factor(Long_FbyHH$Sex)
  Long_FbyHH$Abs.Sec <- c(unlist(lapply(Long_FbyHH$Time.of.Day, T2S)))

  Long_FbyHH$Hour <- round(Long_FbyHH$Abs.Sec / 3600) * 3600
  Long_FbyHH$Hour <- (Long_FbyHH$Hour / 3600)

  # Preparing to remove observations where an individual was unlikely to be feeding,
  # and more likely to be resting.

  high <- mean(na.omit(subset(Long_FbyHH, Bird.ID != "All")$Feeding.Time)) +
    3 * sd(na.omit(subset(Long_FbyHH, Bird.ID != "All")$Feeding.Time))
  # high = 187.893 seconds

  PCFeel.All <- with(Long_FbyHH, data.frame(Feeding.Time, Feeding.Bouts))
  PCFeel.AllOut <- prcomp(na.omit(PCFeel.All), center = TRUE, scale = TRUE)
  cumsum((PCFeel.AllOut$sdev)^2) / sum(PCFeel.AllOut$sdev^2)

  # Removing NAs

  to_rm <- c(
    c(which(is.na(Long_FbyHH$Feeding.Bouts))),
    c(which(is.na(Long_FbyHH$Feeding.Time)))
  )
  
  AllFeed.NoNA <- Long_FbyHH[-c(to_rm), ]
  AllFeed.NoNA$FeedingPC <- PCFeel.AllOut$x[, 1]
  AllFeed.NoNA$FeedingPC.Scale <- c(((AllFeed.NoNA$FeedingPC - min(AllFeed.NoNA$FeedingPC)) /
    (max(AllFeed.NoNA$FeedingPC) - min(AllFeed.NoNA$FeedingPC))) * 100)

  Eye_Mean_Dat <- as.data.frame(Clean %>% group_by(Date.of.Photo, Pen, Treatment, Bird.ID) %>%
                                  summarize(Mean.Eye = mean(Maximum.Eye.Temp, na.rm = T)))

  Long_FbyHH <- subset(Long_FbyHH, Bird.ID != "All")
  Eye_Feed <- left_join(Long_FbyHH, Eye_Mean_Dat, 
    by = c("Date.of.Photo", "Pen", "Treatment", "Bands" = "Bird.ID"))

  Eye_Feed$Treatment <- factor(Eye_Feed$Treatment, ordered = T)
  Eye_Feed$Date.of.Photo <- factor(Eye_Feed$Date.of.Photo)
  Eye_Feed$Bands <- factor(Eye_Feed$Bands)
  Eye_Feed$Pen_Short <- factor(Eye_Feed$Pen_Short)

  Feed.PCA.Dat <- with(Eye_Feed, data.frame(Feeding.Time, Feeding.Bouts))
  Feed.PCA <- prcomp(na.omit(Feed.PCA.Dat), center = TRUE, scale = TRUE)
  cumsum((Feed.PCA$sdev)^2) / sum(Feed.PCA$sdev^2)

  to_rm <- c(
    c(which(is.na(Eye_Feed$Feeding.Bouts))),
    c(which(is.na(Eye_Feed$Feeding.Time)))
  )
  
  Eye.Feed <- Eye_Feed[-c(to_rm), ]
  Eye.Feed$FeedingPC <- Feed.PCA$x[, 1]

  Eye.Feed$FeedingPC <- c(((Eye.Feed$FeedingPC - min(Eye.Feed$FeedingPC)) /
                             (max(Eye.Feed$FeedingPC) - min(Eye.Feed$FeedingPC))) * 100)
}

# Correcting column names for simplicity in analyses.

Eye.Feed <- Eye.Feed %>% rename(
  "Bird.Number" = "Bird.ID", "Bird.ID" = "Bands",
  "Pen_Incorrect" = "Pen", "Pen" = "Pen_Short",
  "Date" = "Date.of.Photo"
)

# Calculating mean feeding metrics per individual.

Eye.Feed.Day <- as.data.frame(Eye.Feed %>% group_by(
  Date, Treatment, Pen_Incorrect,
  Capture.Locale, Locale.Type, Sex,
  Bird.ID
  ) %>%
  summarise(
    Mean.PC = mean(FeedingPC, na.rm = T),
    Mean.Feed.Time = mean(Feeding.Time, na.rm = T),
    Mean.Feed.Bouts = mean(Feeding.Bouts, na.rm = T),
    Mean.Ambient.Temp = mean(Mean.Temp, na.rm = T),
    Mean.Eye.Temp = mean(Mean.Eye, na.rm = T),
    Mean.Time = mean(Hour, na.rm = T)
))

Eye.Feed.Day$Date <- as.Date(as.character(paste(substr(Eye.Feed.Day$Date, 0, 4),
                                                substr(Eye.Feed.Day$Date, 5, 6),
                                                substr(Eye.Feed.Day$Date, 7, 8),
                                                sep = "-"
)))


Eye.Feed.Day$Pen = gsub("[[:space:]].*", "", Eye.Feed.Day$Pen_Incorrect)

# Binding social status metrics.

{
  All_Ranks$Bird.ID = as.character(All_Ranks$Bird.ID)  
  All_Ranks$Bird.ID[c(which(All_Ranks$Bird.ID == "BdAO"))] = "BdAO1"
  All_Ranks$Bird.ID[c(which(All_Ranks$Bird.ID == "YAO"))] = "YAO1"
  All_Ranks$Bird.ID[c(which(All_Ranks$Bird.ID == "BdABl"))] = "BdABl1"
  All_Ranks$Bird.ID = as.factor(All_Ranks$Bird)
  All_Ranks = All_Ranks %>% dplyr::rename("RSE.Rank" = RSE_Rank)
  All_Ranks = All_Ranks %>% dplyr::rename("E.Rank" = Elo_Rank)
}

Feed_Rank <- left_join(Eye.Feed.Day, All_Ranks, by = c("Pen", "Bird.ID"))
Feed_Rank$RSE.Rank = as.integer(as.character(Feed_Rank$RSE.Rank))

# Plotting Rank by Feeding

ggplot(subset(Feed_Rank, !is.na(Treatment)), aes(x = RSE.Rank, y = Mean.PC)) +
  geom_smooth(method = "gam", formula = y~s(x, k= 3, bs = "cr"),
              colour = "black") +
  theme_bw() + xlab("Social Rank") + ylab("Feeding Score (PC1)") +
  facet_wrap(~Treatment)

ggplot(subset(Feed_Rank, !is.na(Treatment)), aes(x = RSE.Rank, y = Mean.PC)) +
  geom_smooth(method = "lm", colour = "black") +
  theme_bw() + xlab("Social Rank") + ylab("Feeding Score (PC1)")

# Running a quick linear model - probing.

summary(lm(Mean.PC ~ RSE.Rank, data = Feed_Rank))

# Significant negative correlation

# Categorizing as dominant or subbordinate and replotting

Statuses = rep("Subordinate", nrow(Feed_Rank))
Statuses[c(which(Feed_Rank$RSE.Rank <= 2))] = "Dominant"

Feed_Rank %>% mutate("Status" = Statuses) %>%
  group_by(Status) %>% dplyr::summarise(Feed_Score = mean(Mean.PC, na.rm = T),
                              SE = sd(na.omit(Mean.PC))/sqrt(n()),
                              Upper.CI = Feed_Score + 1.96*SE,
                              Lower.CI = Feed_Score - 1.96*SE) %>%
  ggplot(aes(x = Status, y = Feed_Score, fill = Status)) + geom_point(size = 4, pch = 21,
    colour = "black") + geom_errorbar(mapping = aes(x = Status, ymin = Lower.CI,
      ymax = Upper.CI), width = 0.3, size = 1, colour = "black") + theme_bw() +
      scale_fill_manual(values = c("royalblue3", "navajowhite"))

t.test(Mean.PC ~ Statuses, data = Feed_Rank)

# Modeling feeding rate by social status.

# Adjusting predictors into necessary categies.

{
Feed_Rank$Mean.Time = as.integer(Feed_Rank$Mean.Time)
Feed_Rank$Mean.Feed.Bouts = as.integer(Feed_Rank$Mean.Feed.Bouts)
Feed_Rank$RSE.Status = factor(Feed_Rank$RSE.Status)
Feed_Rank$Date = factor(Feed_Rank$Date)
Feed_Rank$Pen = factor(Feed_Rank$Pen)
Feed_Rank$Bird.ID = factor(Feed_Rank$Bird.ID)
Feed_Rank$O.Treatment = ordered(Feed_Rank$Treatment)
Feed_Rank$Sex = factor(Feed_Rank$Sex)
}

# Using bam (R package mgcv) to permit implementation of non-gaussian family.

Feed_gam = bam(Mean.Feed.Bouts ~ O.Treatment + Sex + RSE.Status +
  s(Mean.Time, k = 4, bs = "cr") +
  s(Mean.Time, k = 4, by = O.Treatment, bs = "cr", m = 1) +
  s(Mean.Ambient.Temp, k = 4, bs = "cr") +
  s(Mean.Ambient.Temp, by = O.Treatment, bs = "cr", m = 1) +
  s(Date, bs = "re") + s(Pen, bs = "re") +
  s(Bird.ID, bs = "re") + 
  s(Capture.Locale, bs = "re"),
  method = "REML", family = nb(),
  data = Feed_Rank
)

# Assessing suitability of knots for time and ambient temperature.

gam.check(Feed_gam)

# Time appears largely linear. Replacing time as a parametric predictor.

Feed_gam2 = bam(Mean.Feed.Bouts ~ O.Treatment*Mean.Time + Sex + RSE.Status +
  s(Mean.Ambient.Temp, k = 4, bs = "cr") +
  s(Mean.Ambient.Temp, k = 4, by = O.Treatment, bs = "cr", m = 1) +
  s(Date, bs = "re") + s(Pen, bs = "re") +
  s(Bird.ID, bs = "re") + 
  s(Capture.Locale, bs = "re"),
  method = "REML", family = nb(),
  data = Feed_Rank
)

# Assessing distribution of residuals

{
p1 = ggplot(data = data.frame("Res" = resid_gam(Feed_gam2, incl_na = T)),
  aes(x = 1:length(resid_gam(Feed_gam2, incl_na = T)), y = Res)) + my.theme + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalized Residuals")

p2 = ggplot(data = data.frame("Res" = resid_gam(Feed_gam2, incl_na = T),
  "Fit" = predict(Feed_gam2, type = "response")),
  aes(x = Fit, y = Res)) + my.theme + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalized Residuals")

p3 = ggplot(data.frame("Res" = resid_gam(Feed_gam2, incl_na = T)), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(resid_gam(Feed_gam2)), sd = sd(resid_gam(Feed_gam2)))) +
  theme_bw() + my.theme + xlab("Normalized Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(resid_gam(Feed_gam2)) - 3*sd(resid_gam(Feed_gam2))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(resid_gam(Feed_gam2)) + 3*sd(resid_gam(Feed_gam2))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = resid_gam(Feed_gam2, incl_na = T), "Fit" = predict(Feed_gam2)),
  aes(sample=Res)) + stat_qq(colour = "gold") + stat_qq_line() + my.theme + theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Model Residuals")
}

# Quite clean, except for residuals by predicted values. Testing whether status should be 
# included as an interacting predictor.

Feed_Int_gam = bam(Mean.Feed.Bouts ~ O.Treatment*Mean.Time + 
    Sex + O.Treatment*RSE.Status +
      s(Mean.Ambient.Temp, k = 4, bs = "cr") +
      s(Mean.Ambient.Temp, k = 4, by = O.Treatment, bs = "cr", m = 1) +
      s(Date, bs = "re") + s(Pen, bs = "re") +
      s(Bird.ID, bs = "re") +
      s(Capture.Locale, bs = "re"),
      method = "REML", family = nb(),
      data = Feed_Rank
)

summary(Feed_Int_gam)
logLik(Feed_Int_gam)
logLik(Feed_gam2)

test_stat = -2 * (logLik(Feed_gam2) - logLik(Feed_Int_gam))
pchisq(as.numeric(test_stat), df = (74 - 72), lower.tail = FALSE)

# No significant improvement in model. Proceeding to assess residual distribution

{
p1 = ggplot(data = as.data.frame(Feed_Rank %>% 
    mutate("Res" = resid_gam(Feed_gam2, incl_na = T))),
  aes(x = Mean.Time, y = Res)) + geom_point(size = 3, 
    colour = "cornflowerblue", alpha = 0.5) + 
    xlab("Hour") + ylab("Normalized Residuals") + 
    theme_bw() + my.theme

p2 = ggplot(data = as.data.frame(Feed_Rank %>% 
    mutate("Res" = resid_gam(Feed_gam2, incl_na = T))),
    aes(x = Mean.Ambient.Temp, y = Res)) + geom_point(size = 3, 
    colour = "cornflowerblue", alpha = 0.5) + 
    xlab("Ambient Temperature (Degrees C)") + 
    ylab("Normalized Residuals") + 
    theme_bw() + my.theme

p3 = ggplot(data = as.data.frame(Feed_Rank %>% 
    mutate("Res" = resid_gam(Feed_gam2, incl_na = T))),
    aes(x = Sex, y = Res)) + geom_boxplot( 
    colour = "slateblue") + 
    xlab("Sex") + 
    ylab("Normalized Residuals") + 
    theme_bw() + my.theme

p4 = ggplot(data = as.data.frame(Feed_Rank %>% 
    mutate("Res" = resid_gam(Feed_gam2, incl_na = T))),
    aes(x = O.Treatment, y = Res)) + geom_boxplot( 
    colour = "darkorchid") + 
    xlab("Treatment") + 
    ylab("Normalized Residuals") + 
    theme_bw() + my.theme

p5 = ggplot(data = as.data.frame(Feed_Rank %>% 
    mutate("Res" = resid_gam(Feed_gam2, incl_na = T))),
    aes(x = RSE.Status, y = Res)) + geom_boxplot( 
    colour = "black") + 
    xlab("Status") + 
    ylab("Normalized Residuals") + 
    theme_bw() + my.theme

p6 = ggplot(data = as.data.frame(Feed_Rank %>% 
    mutate("Res" = resid_gam(Feed_gam2, incl_na = T))),
  aes(x = Mean.Time, y = Res, fill = O.Treatment)) + 
    geom_point(size = 3, alpha = 0.5, pch = 21) + 
    xlab("Hour") + ylab("Normalized Residuals") + 
    theme_bw() + my.theme + facet_grid(~O.Treatment) + 
    scale_fill_manual(values = c("navajowhite", "slateblue"))

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, top = "Model Residuals")
}

# Residuals appear satisfactorily distributed.

# Summarising and plotting model results.

summary(Feed_gam2)

# Mean feeding rate per social category.

emmeans(Feed_gam2, specs = pairwise ~ RSE.Status, 
  cov.reduce = mean, type = "response",
  nesting = NULL, nesting.order = FALSE,
  p.adjust.method = "bonferroni",
  data = Feed_Rank)

BMD = emmip(Feed_gam2, RSE.Status ~ Mean.Time,
  data = Feed_Rank, nesting = NULL, nesting.order = FALSE,
  type = "response", plot = FALSE, cov.reduce = FALSE,
  CIs = T, at = list(
    Mean.Time = c(8,10,12,14,16)))

BMD$Status = "Subordinate"
BMD$Status[c(which(BMD$RSE.Status == 1))] = "Dominant"
BMD$Status = factor(BMD$Status)

Feed_Plot = ggplot(BMD, aes(x = Mean.Time, y = yvar,
  linetype = Status, fill = Status)) +
  geom_smooth(method = "gam", formula = y~s(x, k = 5, bs = "cr"),
    size = 0.75, se = F, colour = "black") +
  geom_ribbon(mapping = aes(x = Mean.Time, ymin = LCL, ymax = UCL),
    alpha = 0.5, colour = NA) + theme_bw() +
  scale_fill_viridis_d() +
    my.theme_2 + ylab("Feeding Rate (visits/per hour)") + xlab("Time (hour)") +
    labs(fill = "Social Rank", linetype = "Social Rank")

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Feed_Plot_New.pdf", Feed_Plot, dpi = 800,
    width = 8.67, height = 6.67, units = "in")

# With Treatment

BMD = emmip(Feed_gam2, RSE.Status ~ Mean.Time | O.Treatment,
  data = Feed_Rank, nesting = NULL, nesting.order = FALSE, 
  type = "response", cov.reduce = FALSE, plot = FALSE,
  CIs = T, at = list(
    Mean.Time = c(8,10,12,14,16)))

BMD$Status = "Subordinate"
BMD$Status[c(which(BMD$RSE.Status == 1))] = "Dominant"
BMD$Status = factor(BMD$Status)

facet_names_new <- c(`Stress` = "Stress-Exposed",
                    `Rest` = "Control"
                    )

Feed_Plot = ggplot(BMD, aes(x = Mean.Time, y = yvar,
  linetype = Status, fill = Status)) +
  geom_smooth(method = "gam", formula = y~s(x, k = 5, bs = "cr"),
    size = 0.75, se = F, colour = "black") +
  geom_ribbon(mapping = aes(x = Mean.Time, ymin = LCL, ymax = UCL),
    alpha = 0.5, colour = NA) + theme_bw() +
  scale_fill_viridis_d() +
    my.theme_2 + ylab("Feeding Rate (visits/per hour)") + xlab("Time (hour)") +
    labs(fill = "Social Rank", linetype = "Social Rank") +
    facet_wrap(~O.Treatment, labeller = as_labeller(facet_names_new)) +
    theme(strip.text.x = element_text(size = 16),
    strip.background = element_rect(fill="gray90"))

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Feed_Plot_Facet_New.pdf",
  Feed_Plot, dpi = 800, width = 9.67, height = 6.67, units = "in")

## Pulling in Condition/Mass Data for analysis by rank and treatment.

{
  Cond_dat <- data.frame(
    "Bird.ID" = c(
      "BdABd", "BdAO1", "ABlBl", "AYY",
      "YAO1", "OOA", "BdABl1", "YAO2",
      "ABlO", "ABdBd", "YAR", "RAR",
      "ABlR", "ARO", "BdABl2", "BlAR",
      "BdAO2", "AOR", "YAY", "AYBd"
    ),
    "Pen" = c(rep("NW", 5), rep("NE", 5), rep("SW", 5), rep("SE", 5)),
    "Capture.Locale" = c(
      "Corwhin", "Erin", "Ruthven.Park",
      "Brantford", "Guelph", "Ruthven.Park",
      "Ruthven.Park", "Corwhin", "Brantford",
      "Guelph", "Cambridge", "Corwhin",
      "Brantford", "Corwhin", "Corwhin",
      "Corwhin", "Guelph", "Guelph",
      "Guelph", "Erin"
    ),
    "Locale.Type" = c(
      "Rural", "Rural", "Rural", "Urban",
      "Urban", "Rural", "Rural", "Rural",
      "Urban", "Urban", "Urban", "Rural",
      "Urban", "Rural", "Rural", "Rural",
      "Urban", "Urban", "Urban", "Rural"
    ),
    "Sex" = c("Male", "Male", "Male", "Male", "Female",
    "Female", "Female", "Male", "Male", "Female","Male",
    "Female", "Female", "Female", "Male","Male", "Female", 
    "Female", "Male", "Female"
    ),
    "Wing.Chord" = c(
      66.0, 63.0, NA, 64.5, 66.0, 64.5, 65.0,
      62.0, 69.0, 63.0, 69.0, 64.0, 66.0, 62.0,
      68.0, 68.0, 63.0, 64.0, 68.0, 63.0
    ),
    "Tarsus" = c(
      18.6, 18.2, 18.8, 17.9, 18.0, 19.0, 18.1,
      17.8, 18.9, 17.0, 20.0, 18.2, 18.0, 19.1,
      18.0, 19.0, 19.5, 19.0, 19.1, 17.0
    ),
    "Mass" = c(
      11.6, 11.3, 10.7, 10.7, 11.2, 9.9,
      11.3, NA, 11.4, 10.9,
      12.1, 9.9, 10.3, 10.2, 12.1, 11.7,
      10.8, 9.8, 11.5, 10.4
    ),
    "Measure.Time" = c(rep("Rest.Start", 20))
  )

  Rest_End_Cond_dat <- Cond_dat
  Rest_End_Cond_dat$Measure.Time <- c(rep("Rest.End", 20))
  Rest_End_Cond_dat$Mass <- c(
    11.9, 11.1, 10.5, 10.7, 11.8,
    9.9, 11.4, NA, NA, 10.6,
    12.1, 10.0, 10.5, 9.8, NA,
    12.1, 10.2, 8.7, 11.1, 10.0
  )

  Stress_Start_Cond_dat <- Cond_dat
  Stress_Start_Cond_dat$Measure.Time <- c(rep("Stress.Start", 20))
  Stress_Start_Cond_dat$Mass <- c(
    11.9, 11.1, 10.5, 10.7, 11.8,
    10.3, 11.6, 10.9, 11.6, 10.6,
    11.9, 9.7, 10.5, 10.4, 11.8,
    12.1, 10.2, 8.7, 11.1, 10.0
  )

  Stress_End_Cond_dat <- Cond_dat
  Stress_End_Cond_dat$Measure.Time <- c(rep("Stress.End", 20))
  Stress_End_Cond_dat$Mass <- c(
    NA, 11.0, 9.9, 10.2, 11.4,
    9.9, 11.3, NA, 11.4, 10.9,
    12.1, 9.9, 10.3, 10.2, 12.1,
    NA, 9.9, 9.1, 11.1, 10.1
  )

  Cond_dat <- rbind(
    Cond_dat,
    Rest_End_Cond_dat,
    Stress_Start_Cond_dat,
    Stress_End_Cond_dat
  )

  Cond_dat$Condition <- residuals(lm(Mass ~ Wing.Chord, data = Cond_dat, na.action = na.exclude))
  summary(lm(Mass ~ Wing.Chord, data = Cond_dat, na.action = na.exclude))
  # Normalizing condition metric

  Cond_dat$Condition <- (Cond_dat$Condition - min(Cond_dat$Condition, na.rm = T)) /
    (max(Cond_dat$Condition, na.rm = T) - min(Cond_dat$Condition, na.rm = T))

  Start_Cond = Cond_dat[c(grep("Start", Cond_dat$Measure.Time)),]
  Start_Cond$Treatment = "Rest"
  Start_Cond$Treatment[c(grep("Stress", Start_Cond$Measure.Time))] = "Stress"
  End_Cond = Cond_dat[c(grep("End", Cond_dat$Measure.Time)),]
  End_Cond$Treatment = "Rest"
  End_Cond$Treatment[c(grep("Stress", End_Cond$Measure.Time))] = "Stress"
  Start_Cond$Treatment = factor(Start_Cond$Treatment)
  End_Cond$Treatment = factor(End_Cond$Treatment)
  Start_Cond = Start_Cond %>%
    dplyr::select("Bird.ID", "Pen", "Treatment", "Capture.Locale", "Locale.Type", "Sex", "Wing.Chord", "Mass", "Condition")
  End_Cond = End_Cond %>%
    dplyr::select("Bird.ID", "Pen", "Treatment", "Wing.Chord", "Mass", "Condition")

  Cond_delta = left_join(End_Cond, Start_Cond, by = c("Bird.ID", "Pen", "Treatment"))

  Means = c()
  for (i in 1:nrow(Cond_delta)){
    if (is.na(Cond_delta$Condition.x[i]) & !is.na(Cond_delta$Condition.y[i])){
    Means[i] = Cond_delta$Condition.y[i]
    } else if (!is.na(Cond_delta$Condition.x[i]) & is.na(Cond_delta$Condition.y[i])){
    Means[i] = Cond_delta$Condition.x[i]
    } else if (is.na(Cond_delta$Condition.x[i]) & is.na(Cond_delta$Condition.y[i])){
    Means[i] = NA
    } else {
    Means[i] = mean(Cond_delta$Condition.x[i], Cond_delta$Condition.y[i])
    }
  }
  Cond_delta$Condition = Means
}

# Fusing data with social status metrics for plotting.

Cond_Bound = left_join(Cond_delta, All_Ranks, by = c("Pen", "Bird.ID"))

ggplot(subset(Cond_Bound, !is.na(Treatment)), aes(x = RSE.Rank, y = Condition)) +
  geom_smooth(method = "lm", colour = "black", fill = "royalblue3") +
  theme_bw() + xlab("Social Rank") + ylab("Condition")

ggplot(subset(Cond_Bound, !is.na(Treatment)), aes(x = E.Rank, y = Condition)) +
  geom_smooth(method = "lm", colour = "black", fill = "royalblue3") +
  theme_bw() + xlab("Social Rank") + ylab("Condition")

# Testing according to binomial status

Cond.Statuses = rep("Subordinate", nrow(Cond_Bound))
Cond.Statuses[c(which(Cond_Bound$RSE.Rank <= 2))] = "Dominant"

Cond_Bound %>% mutate("Status" = Cond.Statuses) %>%
    group_by(Status) %>% dplyr::summarise(Mean.Condition = mean(Condition, na.rm = T),
                                SE = sd(na.omit(Condition))/sqrt(n()),
                                Upper.CI = Mean.Condition + 1.96*SE,
                                Lower.CI = Mean.Condition - 1.96*SE) %>%
    ggplot(aes(x = Status, y = Mean.Condition, fill = Status)) + geom_point(size = 4, pch = 21,
      colour = "black") + geom_errorbar(mapping = aes(x = Status, ymin = Lower.CI,
        ymax = Upper.CI), width = 0.3, size = 1, colour = "black") + theme_bw() +
        scale_fill_manual(values = c("royalblue3", "navajowhite"))

# Mild difference

# Formally modelling the effects of social status on mass
# and the change in mass of individuals across treatments.

Sex = c()

for (i in 1:nrow(Cond_Bound)) {
  if (is.na(Cond_Bound$Bird.ID[i])) {
    Sex[i] = NA
  } else if (Cond_Bound$Bird.ID[i] == "BdABd") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "BdAO1") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "ABlBl") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "AYY") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "YAO1") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "OOA") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "BdABl1") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "YAO2") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "ABlO") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "ABdBd") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "YAR") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "RAR") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "ABlR") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "ARO") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "BdABl2") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "BlAR") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "BdAO2") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "AOR") {
    Sex[i] = "Female"
  } else if (Cond_Bound$Bird.ID[i] == "YAY") {
    Sex[i] = "Male"
  } else if (Cond_Bound$Bird.ID[i] == "AYBd") {
    Sex[i] = "Female"
  } else {
  }
}

Cond_Bound$Sex = factor(Sex)

# First, however, saving the data frame for future usage. 

#write.csv(Cond_Bound, "/home/joshk/git_repositories/BCCH_Dominance/Data/Condition_Final.csv")
Cond_Bound = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Condition_Final.csv")

Counts = Cond_Bound %>% dplyr::select(Bird.ID, Treatment, Mass.x, Mass.y) %>% 
group_by(Treatment)

# Aqcuiring sample sizes per grouping, and preparing summarised data.

length(which(is.na(subset(Counts, Treatment == "Stress")$Mass.x)))
length(which(is.na(subset(Counts, Treatment == "Rest")$Mass.x)))
length(which(is.na(subset(Counts, Treatment == "Stress")$Mass.y)))
length(which(is.na(subset(Counts, Treatment == "Rest")$Mass.y)))

CB_NaFree = as.data.frame(Cond_Bound %>% 
    dplyr::select(Mass.x, Mass.y, Treatment, RSE.Status, Bird.ID, Pen, Sex) %>% 
    na.omit(.)) %>% mutate(Delta = Mass.x - Mass.y)
CB_NaFree$RSE.Status = factor(CB_NaFree$RSE.Status)

nrow(subset(CB_NaFree, Treatment == "Stress"))
nrow(subset(CB_NaFree, Treatment == "Rest"))

nrow(subset(CB_NaFree, Treatment == "Stress" & RSE.Status == "1"))
nrow(subset(CB_NaFree, Treatment == "Stress" & RSE.Status == "2"))

nrow(subset(CB_NaFree, Treatment == "Rest" & RSE.Status == "1"))
nrow(subset(CB_NaFree, Treatment == "Rest" & RSE.Status == "2"))

# Modeling the change in mass by status and treatment.

mod = glmmTMB(Delta ~ Treatment*RSE.Status + Sex + (1|Pen) + (1|Bird.ID),
  weights = as.integer(Treatment), data = CB_NaFree)

summary(mod)

# Flight enclosure explains negligable variance in data. Removing.

mod = glmmTMB(Delta ~ Treatment*RSE.Status + Sex + (1|Bird.ID),
  weights = as.integer(Treatment), data = CB_NaFree)

Diag_plot = sjPlot::plot_model(mod, "diag")
Diag_plot[[1]]
Diag_plot[[2]]
Diag_plot[[3]]
Diag_plot[[4]]

# Residuals appear acceptably normal and heterogenous. One tentative outlier,
# however, no logical reason to remove this individual. Plotting
# residuals by predictors.

{
p1 = na.omit(CB_NaFree) %>% mutate("Res" = residuals(mod, type = "pearson"),
    "Fit" = predict(mod, newdata=na.omit(CB_NaFree))) %>%
    ggplot(aes(x = Treatment, y = Res)) +
    geom_boxplot(colour = "seagreen") + 
    theme_bw() + ylab("Pearson Residuals") + facet_grid(~RSE.Status)

p2 = na.omit(CB_NaFree) %>% mutate("Res" = residuals(mod, type = "pearson"),
    "Fit" = predict(mod, newdata=na.omit(CB_NaFree))) %>%
    ggplot(aes(x = RSE.Status, y = Res)) +
    geom_boxplot(colour = "darkorchid") + 
    theme_bw() + ylab("Pearson Residuals")

p3 = na.omit(CB_NaFree) %>% mutate("Res" = residuals(mod, type = "pearson"),
    "Fit" = predict(mod, newdata=na.omit(CB_NaFree))) %>%
    ggplot(aes(x = Sex, y = Res)) +
    geom_boxplot(colour = "grey60") + 
    theme_bw() + ylab("Pearson Residuals")

grid.arrange(p1, p2, p3, nrow = 2, top = "Model Residuals")
}

# Again, distributions appear acceptable.

# Assessing and plotting trends across groupings.

emmeans(mod, specs = pairwise ~ RSE.Status | Treatment, 
  cov.reduce = mean, type = "response",
  nesting = NULL, nesting.order = FALSE,
  p.adjust.method = "bonferroni",
  data = CB_NaFree)

cond_em = emmip(mod, RSE.Status~Treatment, nesting = NULL, 
  nesting.order = FALSE, plot = F,
  data = CB_NaFree, CIs = T, cov.reduce = F)
cond_em$Status = "Dominant"
cond_em$Status[c(which(cond_em$RSE.Status == "2"))] = "Subordinate"
cond_em$Treatment_New = "Control"
cond_em$Treatment_New[c(which(cond_em$Treatment == "Stress"))] = "Stress"

Mass_Plot = cond_em %>% 
    ggplot(aes(x = Treatment_New, y = yvar, fill = Status)) + 
    geom_errorbar(mapping = aes(x = Treatment_New, ymin = LCL, ymax = UCL),
        width = 0.4, fill = 1, colour = "black", 
        position = position_dodge(width = 0.5)) +
    geom_point(size = 4, colour = "black", pch = 21,
        position = position_dodge(width = 0.5)) + 
    theme_bw() + scale_fill_viridis_d() + 
    xlab("Treatment") + ylab("Change in Mass (g)") + my.theme_2 + geom_signif(
        xmin=1.75,
        xmax=2.25,
        y_position = 0.36, annotations = c("*"), size = 1,
        textsize = 8
      ) + ylim(c(-0.4,0.4))

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Delta_Mass_Rev.pdf",
  Mass_Plot, dpi = 800, width = 9.67, height = 6.67, units = "in")

summary(mod)

# Running planned comparison under stress.

emmeans(mod, specs = pairwise ~ RSE.Status|Treatment, nesting = NULL, 
  nesting.order = FALSE, data = CB_NaFree, CIs = T, adjust = "bonferroni")

# Quantifying change.

emmeans(mod, specs = ~RSE.Status, at = list(Treatment = "Stress"),
  data = CB_NaFree, CIs = T, adjust = "bonferroni")

# Quick test of final mass across treatments

Cond_Bound$RSE.Status = factor(Cond_Bound$RSE.Status)

Condition_mod = glmmTMB(Mass.x ~ Treatment*RSE.Status + Sex + (1|Bird.ID),
  weights = as.integer(Treatment), data = Cond_Bound)

# Counting sample size 

Cond_Bound %>% na.omit(.) %>% group_by(RSE.Status, Treatment) %>% dplyr::summarize(n())

condition_em = emmip(Condition_mod, RSE.Status~Treatment, nesting = NULL, 
  nesting.order = FALSE, plot = F,
  data = Cond_Bound, CIs = T, cov.reduce = F)

condition_em$Status = "Dominant"
condition_em$Status[c(which(condition_em$RSE.Status == "2"))] = "Subordinate"
condition_em$Treatment_New = "Control"
condition_em$Treatment_New[c(which(condition_em$Treatment == "Stress"))] = "Stress"

Cond_Plot = condition_em %>% mutate(Response = yvar, UCL = UCL, LCL = LCL) %>%
    ggplot(aes(x = Treatment_New, y = Response, fill = Status)) + 
    geom_errorbar(mapping = aes(x = Treatment_New, ymin = LCL, ymax = UCL),
        width = 0.4, size = 1, colour = "black", 
        position = position_dodge(width = 0.5)) +
    geom_point(size = 4, colour = "black", pch = 21,
        position = position_dodge(width = 0.5)) + 
    theme_bw() + scale_fill_viridis_d() + 
    xlab("Treatment") + ylab("Mass (g)") + my.theme_2

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Conditon_Plot_Final.pdf",
  Cond_Plot, dpi = 2000, width = 8.67, height = 6.67, units = "in")

emmeans(Condition_mod, specs = pairwise ~ RSE.Status | Treatment, 
  cov.reduce = mean, type = "response",
  nesting = NULL, nesting.order = FALSE,
  p.adjust.method = "bonferroni",
  data = Cond_Bound)

# No different between groups here.

summary(Condition_mod)

###---- Rank and Eye Temperature ------- ##

# Reloading eye data

{
  All_IDs = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/All Combined Thermal Data.csv")
  All_IDs$Maximum.Eye.Temp = as.numeric(as.character(All_IDs$Maximum.Eye.Temp))
  All_IDs$Amb.Temp = as.numeric(as.character(All_IDs$Amb.Temp))

  ResFrame = data.frame(
    "Row" = 1:nrow(All_IDs),
    "Residuals" = residuals(lm(Maximum.Eye.Temp ~ Amb.Temp,
      data = All_IDs, na.action = na.exclude
      ))
    )

    to_remove = c(
        c(which(All_IDs$Maximum.Eye.Temp <=
          mean(na.omit(All_IDs$Maximum.Eye.Temp)) -
          4 * sd(na.omit(All_IDs$Maximum.Eye.Temp)))),
        c(which(ResFrame$Residuals <=
          mean(na.omit(ResFrame$Residuals)) -
          4 * sd(na.omit(ResFrame$Residuals)))),
        c(which(ResFrame$Residuals >=
          mean(na.omit(ResFrame$Residuals)) +
          4 * sd(na.omit(ResFrame$Residuals))))
    )

    nrow(All_IDs[to_remove, ])

    # Minimal, and therefore keeping extremes.

    All_IDs$Treatment = ordered(All_IDs$Treatment, levels = c("Rest", "Stress"))
    All_IDs$Bird.ID = factor(All_IDs$Bird.ID)
    All_IDs$Date.of.Photo = factor(All_IDs$Date.of.Photo)
    All_IDs$Pen_Short = factor(All_IDs$Pen_Short)

    Clean = All_IDs

    Capture.Locale = c()
    Locale.Type = c()
    Sex = c()

    for (i in 1:nrow(Clean)) {
      if (is.na(Clean$Bird.ID[i])) {
        Capture.Locale[i] = NA
        Locale.Type[i] = NA
        Sex[i] = NA
      } else if (Clean$Bird.ID[i] == "BdABd") {
        Capture.Locale[i] = "Corwhin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "BdAO1") {
        Capture.Locale[i] = "Erin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "ABlBl") {
        Capture.Locale[i] = "Ruthven.Park"
        Locale.Type[i] = "Urban"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "AYY") {
        Capture.Locale[i] = "Brantford"
        Locale.Type[i] = "Urban"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "YAO1") {
        Capture.Locale[i] = "Guelph"
        Locale.Type[i] = "Urban"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "OOA") {
        Capture.Locale[i] = "Ruthven.Park"
        Locale.Type[i] = "Rural"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "BdABl1") {
        Capture.Locale[i] = "Ruthven.Park"
        Locale.Type[i] = "Rural"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "YAO2") {
        Capture.Locale[i] = "Corwhin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "ABlO") {
        Capture.Locale[i] = "Brantford"
        Locale.Type[i] = "Urban"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "ABdBd") {
        Capture.Locale[i] = "Guelph"
        Locale.Type[i] = "Urban"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "YAR") {
        Capture.Locale[i] = "Cambridge"
        Locale.Type[i] = "Urban"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "RAR") {
        Capture.Locale[i] = "Corwhin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "ABlR") {
        Capture.Locale[i] = "Brantford"
        Locale.Type[i] = "Urban"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "ARO") {
        Capture.Locale[i] = "Corwhin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "BdABl2") {
        Capture.Locale[i] = "Corwhin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "BlAR") {
        Capture.Locale[i] = "Corwhin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "BdAO2") {
        Capture.Locale[i] = "Guelph"
        Locale.Type[i] = "Urban"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "AOR") {
        Capture.Locale[i] = "Guelph"
        Locale.Type[i] = "Urban"
        Sex[i] = "Female"
      } else if (Clean$Bird.ID[i] == "YAY") {
        Capture.Locale[i] = "Guelph"
        Locale.Type[i] = "Urban"
        Sex[i] = "Male"
      } else if (Clean$Bird.ID[i] == "AYBd") {
        Capture.Locale[i] = "Erin"
        Locale.Type[i] = "Rural"
        Sex[i] = "Female"
      } else {
        Capture.Locale[i] = NA
        Locale.Type[i] = NA
      }
    }

    Clean$Capture.Locale = factor(Capture.Locale)
    Clean$Locale.Type = factor(Locale.Type)
    Clean$Sex = factor(Sex)
    Clean$Maximum.Eye.Temp[c(which(Clean$Eye.Quality != "Good"))] = NA
}

Clean_Ranks = left_join(Clean, All_Ranks, by = c("Bird.ID"))
Clean_Ranks$Hour = floor(Clean_Ranks$Absolute.Seconds/3600)

### Effect of rank on eye temperature responses to stress 

#write.csv(Clean_Ranks, "C:/Users/joshk/Documents/BCCH_Dominance/BCCH-Dominance-Analyses/Clean_Ranks.csv")
Clean_Ranks = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Clean_Ranks.csv")

Clean_Ranks$Treatment.O = factor(Clean_Ranks$Treatment, ordered = "T")
Clean_Ranks$RSE.Status.O = factor(Clean_Ranks$RSE.Status, ordered = "T")
Clean_Ranks$Date.of.Photo = factor(Clean_Ranks$Date.of.Photo)
Clean_Ranks$RSE.Status = factor(Clean_Ranks$RSE.Status)

# Collecting samples sizes

Clean_Ranks %>% filter(Eye.Quality != "Poor Eye" & !is.na(Maximum.Eye.Temp) & !is.na(Bird.ID)) %>%
  group_by(Treatment.O, RSE.Status.O) %>% summarise(Count = n())

Clean_Ranks %>% filter(Eye.Quality != "Poor Eye" & !is.na(Maximum.Eye.Temp) & !is.na(Bird.ID)) %>%
  group_by(Treatment.O, RSE.Status.O) %>% summarise(Count = n()) %>% 
  summarise(Sum = sum(Count))

Clean_Ranks %>% filter(Eye.Quality != "Poor Eye" & !is.na(Maximum.Eye.Temp) & !is.na(Bird.ID)) %>%
  group_by(Treatment.O, RSE.Status.O) %>% summarise(Count = n()) %>% 
  ungroup() %>% summarise(Sum = sum(Count))

# Averaging by hour

hour_means = as.data.frame((Clean_Ranks) %>% filter(Eye.Quality != "Poor Eye") %>%
  group_by(Date.of.Photo,Hour,Pen_Short,Treatment,Bird.ID, Capture.Locale,Sex,RSE.Status) %>%
  dplyr::summarize(Mean.Eye = mean(Maximum.Eye.Temp, na.rm = T),
    Mean.Amb = mean(Amb.Temp, na.rm = T)))

# Correcting factor types

hour_means$Treatment.O = factor(hour_means$Treatment, ordered = "T")
hour_means$RSE.Status.O = factor(hour_means$RSE.Status, ordered = "T")

# Eye temperature modeling

# Note that capture locale has been excluded because we have previously shown
# that it does not influence surface temperature (Robertson et al, 2020; JEB).

# Adding orientation of pens to account for change in solar radiation by time.

hour_means = hour_means %>% mutate("Direction" = "West")
hour_means$Direction[c(which(hour_means$Pen_Short == "NE" | hour_means$Pen_Short == "SE"))] = "East"
hour_means$Direction = factor(hour_means$Direction)

Eye_Gamm <- mgcv::gamm(Mean.Eye ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1),
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  method = "REML", na.action = na.exclude, data = hour_means)
)

hour_means %>% filter(!is.na(Bird.ID)) %>% group_by(Treatment.O) %>%
  summarise(Mean_Temp = mean(Mean.Eye, na.rm = T),
    SD_Temp = sd(Mean.Eye, na.rm = T))

hour_means %>% filter(!is.na(Bird.ID)) %>% group_by(RSE.Status.O) %>%
  summarise(Mean_Temp = mean(Mean.Eye, na.rm = T),
    SD_Temp = sd(Mean.Eye, na.rm = T))

hour_means %>% filter(!is.na(Bird.ID)) %>% group_by(Sex) %>%
  summarise(Mean_Temp = mean(Mean.Eye, na.rm = T),
    SD_Temp = sd(Mean.Eye, na.rm = T))

# Checking predictive capacity 

summary(lm(hour_means$Mean.Eye~predict(Eye_Gamm$gam, type = "response", newdata = hour_means)))

Eye_Gamm <- mgcv::gamm(Mean.Eye ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1),
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  method = "REML", na.action = na.omit, data = hour_means)

Null_Gamm <- mgcv::gamm(Mean.Eye ~ 1, 
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  data=hour_means)

(sum(as.numeric(residuals(Null_Gamm$gam, type = "deviance"))^2) - 
sum(as.numeric(residuals(Eye_Gamm$gam, type = "deviance"))^2))/
(sum(as.numeric(residuals(Null_Gamm$gam, type = "deviance"))^2))

# R-squared is reasonable, deviance explained high, and relationship between 
# fitte values and response is very close to 1.
# Assessing residuals, but first re-running model with NAs included.

Eye_Gamm <- mgcv::gamm(Mean.Eye ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1),
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  method = "REML", na.action = na.exclude, data = hour_means)


{
p1 = ggplot(data = data.frame("Res" = residuals(Eye_Gamm$lme, type = "normalized")),
  aes(x = 1:length(residuals(Eye_Gamm$lme, type = "normalized")), y = Res)) + my.theme + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalized Residuals")

p2 = ggplot(data = data.frame("Res" = residuals(Eye_Gamm$lme, type = "normalized"), 
  "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))),
  aes(x = Fit, y = Res)) + my.theme + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalized Residuals")

p3 = ggplot(data.frame("Res" = residuals(Eye_Gamm$lme, type = "normalized")), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(residuals(Eye_Gamm$lme, type = "normalized")), 
    sd = sd(residuals(Eye_Gamm$lme, type = "normalized")))) +
  theme_bw() + my.theme + xlab("Normalized Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(residuals(Eye_Gamm$lme, type = "normalized")) - 
    3*sd(residuals(Eye_Gamm$lme, type = "normalized"))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(residuals(Eye_Gamm$lme, type = "normalized")) + 
    3*sd(residuals(Eye_Gamm$lme, type = "normalized"))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = residuals(Eye_Gamm$lme, type = "normalized"), 
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Model Residuals")
}

# Fair, with very few outliers. Some left-skew to residuals, but histogram looks
# reasonable. Checking for patterns across predictors.

{
p1 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Mean.Amb, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "slateblue", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~Treatment)

p2 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Hour, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "aquamarine", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~Treatment)

p3 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Mean.Amb, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "slateblue", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(Treatment~RSE.Status.O)

p4 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Hour, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "aquamarine", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(Treatment~RSE.Status.O)

p5 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Treatment, y = Res)) +
    geom_boxplot(fill = "navajowhite", colour = "black") + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~RSE.Status.O)

p6 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Pen_Short, y = Res)) +
    geom_boxplot(fill = "navajowhite", colour = "black") + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~RSE.Status.O)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, top = "Model Residuals")
}

# No clear trends. Assessing autocorrelation. 

acf(residuals(Eye_Gamm$lme, type = "normalized"))

# Autocorrelation present, but minimal. This is 
# likely due to uncontrolled environmental parameters 
# (i.e. humidity, rainfall). Testing whether there is 
# temperoral autocorrelation within individuals by 
# re-ordering data frame by bird.

hour_means_ordered = as.data.frame(hour_means %>% 
  arrange(., Bird.ID, Date.of.Photo, Hour))

Eye_Gamm_ordered <- mgcv::gamm(Mean.Eye ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1), 
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1),  
  method = "REML", na.action = na.exclude, data = hour_means_ordered)
)

acf(residuals(Eye_Gamm_ordered$lme, type = "normalized"))
acf(residuals(Eye_Gamm_ordered$lme, type = "normalized"), plot = F)$acf[2]

# Autocorrelation no longer present (rho = -0.151). External, environmental parameters
# are likely at play here, and should be influencing each treatment
# and status equivalently. Testing whether this reordering is likely
# to influence results.

summary(Eye_Gamm$gam)
summary(Eye_Gamm_ordered$gam)

# Significance levels remain the same. Autocorrelation is therefore being ignored, 
# given that correlation is both across individuals and across treatments.
# Now checking for colinearity between sex and status using cars VIF function, with slight 
# modifications to function with gamm.

vif_gamm = function (mod, ...) {
    if (any(is.na(fixef(mod$lme)))) 
        stop("there are aliased coefficients in the model")
    v <- as.matrix(vcov(mod$lme))
    assign <- c(0:(nrow(summary(mod$lme)$tTable[-c(grep("Fx", rownames(summary(mod$lme)$tTable))),])-1))
    if (grep("Intercept", names(fixef(mod$lme)[1])) == 1) {
        v <- v[-1, -1]
        assign <- assign[-1]
    }
    else warning("No intercept: vifs may not be sensible.")
    to_scan <- gsub(".L", "", names(fixef(mod$lme))[-c(grep("Fx", names(fixef(mod$lme))))])
    terms <- c()

    for (i in 1:length(labels(terms(mod$gam)))){
      if (identical(grep(labels(terms(mod$gam))[i], to_scan), 
      integer(0)) == FALSE){
      terms <- c(terms, labels(terms(mod$gam))[i])
      }
    }

    n.terms <- length(terms)
    if (n.terms < 2) 
        stop("model contains fewer than 2 terms")
    R <- cov2cor(v)
    detR <- base::det(R)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    for (term in 1:n.terms) {
        subs <- which(assign == term)
        result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, 
            -subs]))/detR
        result[term, 2] <- length(subs)
    }
    if (all(result[, 2] == 1)){ 
        result <- result[, 1]
    } else { result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    }
    result
}

vif_gamm(Eye_Gamm)

# All VIF values are quite low, suggesting very limited likelihood of multicollinearity.

# Checking whether predicted values are biased by treatment and status.

p1 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Fit, y = Mean.Eye)) +
    geom_point(size = 3, fill = "navajowhite", colour = "black", alpha = 0.5, pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = F) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(RSE.Status.O~Treatment.O)

p2 = na.omit(hour_means) %>% mutate("Res" = residuals(Eye_Gamm$lme, type = "normalized"),
    "Fit" = predict(Eye_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Fit, y = Res)) +
    geom_point(size = 3, fill = "slateblue", colour = "black", alpha = 0.5, pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = F) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(RSE.Status.O~Treatment.O)

grid.arrange(p1, p2, top = "Model Residuals")

# No detectable bias.
# Proceeding to model summary and preliminary plots.

EM_Grid = ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = seq(2.5, 38.5, by = I(35/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

emmip(EM_Grid, Treatment.O ~ Mean.Amb | RSE.Status.O,
    cov.reduce = F, CIs = T, data = hour_means)

summary(Eye_Gamm$gam)

# Coefficient summary 

Coefs = data.frame("Predictors" = names(Eye_Gamm$gam$coefficients),
  "Coefficients" = as.numeric(Eye_Gamm$gam$coefficients),
  "E.D.F." = as.numeric(Eye_Gamm$gam$edf),
  "S.E.M" = as.numeric(sqrt(diag(vcov(Eye_Gamm$gam, unconditional = TRUE))))
  )

print(Coefs)

df.residual(Eye_Gamm$gam) # 1005.593

# Averages summaries across splines

Coefs$Predictor_Combined = gsub("\\\\.[[:digit:]].*", "", Coefs[,1])
Spline_Sum = Coefs %>% group_by(Predictor_Combined) %>%
  summarise(Mean.Beta = mean(Coefficients, na.rm = T),
    Mean.SE = mean(S.E.M, na.rm = T),
    edf = sum(E.D.F.))

print(Spline_Sum)

# Running planned comparison by temperature zone

{
    Sub_TNZ =  ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = seq(2.5, 13.99, by = I((13.99-2.5)/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

    TNZ =  ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = seq(14, 30, by = I((30-14)/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

    Supra_TNZ =  ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = seq(30.01, 38.5, by = I((38.5-30.01)/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

    Sub_TNZ_Cont = emmeans(Sub_TNZ, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

    TNZ_Cont = emmeans(TNZ, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

    Supra_TNZ_Cont = emmeans(Supra_TNZ, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)
}

print(Sub_TNZ_Cont)
print(TNZ_Cont)
print(Supra_TNZ_Cont)

# Between status at given temperature zones

{
Sub_TNZ_Cont_Rev = emmeans(Sub_TNZ, specs = pairwise ~ RSE.Status.O | Treatment.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

TNZ_Cont_Rev = emmeans(TNZ, specs = pairwise ~ RSE.Status.O | Treatment.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

Supra_TNZ_Cont_Rev = emmeans(Supra_TNZ, specs = pairwise ~ RSE.Status.O | Treatment.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)
}

print(Sub_TNZ_Cont_Rev)
print(TNZ_Cont_Rev)
print(Supra_TNZ_Cont_Rev)

# Interesting; no trends observed. 

# Assessing means at minimum and maximum observed temperatures.

Max_Temp = max(Clean_Ranks$Amb.Temp, na.rm = T)
Min_Temp = min(Clean_Ranks$Amb.Temp, na.rm = T)

EM_Grid_Lower_Limit = ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = Min_Temp), 
    cov.reduce = mean, type = "response",
    nesting = NULL, nesting.order = FALSE,
    data = hour_means)

EM_Grid_Upper_Limit = ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = Max_Temp), 
    cov.reduce = mean, type = "response",
    nesting = NULL, nesting.order = FALSE,
    data = hour_means)

emmeans(EM_Grid_Lower_Limit, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

emmeans(EM_Grid_Upper_Limit, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

# Pulling marginal means by treatment and plotting more comprehensively

EM_Grid = ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = seq(2.5, 38.5, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

MMeans = emmip(Eye_Gamm, Treatment.O ~ Mean.Amb | RSE.Status.O,
  at = list(Mean.Amb = seq(2.5, 38.5, by = I(36/4))), nesting = NULL,
  nesting.order = FALSE,
  data = hour_means, type = "response", CIs = T, 
  cov.reduce = FALSE, plot = F)

MMeans$RSE.Status.O = factor(MMeans$RSE.Status.O)
MMeans$Status = "Subordinate"
MMeans$Status[c(which(MMeans$RSE.Status.O == "1"))] = "Dominant"
MMeans$Treatment = "Stress-Exposed"
MMeans$Treatment[c(which(MMeans$Treatment.O == "Rest"))] = "Control"
MMeans$Status = factor(MMeans$Status)
MMeans$Treatment = factor(MMeans$Treatment)

Raw_Points = as.data.frame(hour_means %>% 
    filter(!is.na(RSE.Status.O)) %>% 
    mutate("Status" = "Dominant")) %>%
    mutate(Status = as.character(Status), 
    Treatment = as.character(Treatment))

Raw_Points$Status[c(which(Raw_Points$RSE.Status.O == "2"))] = "Subordinate"
Raw_Points$Treatment[c(which(Raw_Points$Treatment == "Rest"))] = "Control"
Raw_Points$Treatment[c(which(Raw_Points$Treatment == "Stress"))] = "Stress-Exposed"

Raw_Points = Raw_Points %>% 
  mutate(Status = as.factor(Status), Treatment = as.factor(Treatment))

# Without dots

MMean_plot = MMeans %>% ggplot(aes(x = Mean.Amb, y = yvar, fill = Treatment, linetype = Treatment)) +
  geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cs"), size = 1, colour = "black", se = F) +
  geom_ribbon(mapping = aes(x = Mean.Amb, ymin = LCL, ymax = UCL), alpha = 0.5, colour = NA) +
  theme_bw() + scale_fill_manual(values = c(viridis::viridis(n = 6)[2], "mediumaquamarine")) + facet_wrap(~Status) +
  xlab("Ambient Temperature (°C)") + ylab("Maximum Eye Temperature (°C)") +
  my.theme_2 + theme(strip.text.x = element_text(size = 16),
  strip.background = element_rect(fill="gray90")) + 
  annotate("rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf, fill = "grey20", alpha = 0.4) + 
  annotate("text", label = "TNZ", size = 8, x = 22, y = 39, colour = "black")

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Means_Plot_No_Dots.pdf", MMean_plot, dpi = 2000,
  width = 9, height = 6.67, units = "in")

# With dots 

MMean_plot = MMeans %>% ggplot(aes(x = Mean.Amb, y = yvar, fill = Treatment, linetype = Treatment)) +
  geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cs"), size = 1, colour = "black", se = F) +
  geom_ribbon(mapping = aes(x = Mean.Amb, ymin = LCL, ymax = UCL), alpha = 0.5, colour = NA) +
  geom_point(data = Raw_Points, aes(x = Mean.Amb, y = Mean.Eye),
  pch = 21, size = 2, alpha = 0.3, colour = "black") + 
  theme_bw() + scale_fill_manual(values = c(viridis::viridis(n = 6)[2], "mediumaquamarine")) + facet_wrap(~Status) +
  xlab("Ambient Temperature (°C)") + ylab("Maximum Eye Temperature (°C)") +
  my.theme_2 + theme(strip.text.x = element_text(size = 16),
  strip.background = element_rect(fill="gray90")) + 
  annotate("rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf, fill = "grey20", alpha = 0.4) + 
  annotate("text", label = "TNZ", size = 8, x = 22, y = 39, colour = "black")

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Means_Plot_New_Dots.pdf", MMean_plot, dpi = 800,
  width = 9, height = 6.67, units = "in")

# By temperature zone 

EM_Grid_Wide = ref_grid(Eye_Gamm,  
    at = list(Mean.Amb = seq(2.5, 38.5, by = 0.1)), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

Grid_Wide_Plot = emmip(EM_Grid_Wide, Treatment.O ~ Mean.Amb | RSE.Status.O,
  data = hour_means, type = "response", CIs = T, 
  cov.reduce = FALSE, plot = F)

Grid_Wide_Plot$Zone = "Thermoneutrality"
Grid_Wide_Plot$Zone[c(which(Grid_Wide_Plot$Mean.Amb < 14))] = "Below\\nThermoneutrality"
Grid_Wide_Plot$Zone[c(which(Grid_Wide_Plot$Mean.Amb > 30))] = "Above\\nThermoneutrality"
Grid_Wide_Plot$Zone = factor(Grid_Wide_Plot$Zone, levels = c("Below\\nThermoneutrality", 
"Thermoneutrality", "Above\\nThermoneutrality"))

hour_means$Zone = "Thermoneutrality"
hour_means$Zone[c(which(hour_means$Mean.Amb < 14))] = "Below\\nThermoneutrality"
hour_means$Zone[c(which(hour_means$Mean.Amb > 30))] = "Above\\nThermoneutrality"
hour_means$Zone = factor(hour_means$Zone, levels = c("Below\\nThermoneutrality", 
"Thermoneutrality", "Above\\nThermoneutrality"))

Grid_Wide_Plot %>% group_by(Zone, Treatment.O, RSE.Status.O) %>%
  summarise(Mean_Eye = mean(yvar), SE = sd(yvar)/sqrt(n()), 
  LCL = Mean_Eye - 1.96*SE, UCL = Mean_Eye + 1.96*SE) %>%
ggplot(aes(x = Zone, y = Mean_Eye, fill = Treatment.O)) + 
  geom_errorbar(aes(x = Zone, ymin = LCL, ymax = UCL), size = 1, 
  width = 0.5, colour = "black", position = position_dodge(width = 0.5)) + 
  geom_point(pch = 21, size = 5, position = position_dodge(width = 0.5),
  colour = "black") + 
  geom_point(data = subset(hour_means, !is.na(RSE.Status.O)), aes(x = Zone, y = Mean.Eye),
  size = 1, alpha = 0.2, position = position_dodge(width = 0.5),
  colour = "black") + theme_bw() + facet_wrap(~RSE.Status.O)

# As boxplot

facet_names_new <- c(`1` = "Dominant",
                    `2` = "Subordinate"
                    )
Stars = data.frame("Treatment.O" = c("Rest","Stress","Stress"), 
    "RSE.Status.O" = c("1","2","2"), 
    "Zone" = c("Above\\nThermoneutrality", "Below\\nThermoneutrality",
    "Above\\nThermoneutrality"), "yvar" = c(41,32,41))

Stars$Treatment.O = factor(Stars$Treatment.O)    
Stars$RSE.Status.O = factor(Stars$RSE.Status.O)    
Stars$Zone = factor(Stars$Zone)    

EM_Boxplot = Grid_Wide_Plot %>% 
ggplot(aes(x = Zone, y = yvar, colour = Treatment.O)) + 
  geom_boxplot(size = 0.75, position = position_dodge(width = 0.75)) + 
  geom_jitter(data = as.data.frame(hour_means %>% filter(!is.na(RSE.Status.O)) %>%
  mutate("yvar" = Mean.Eye)), aes(x = Zone, y = yvar, colour = Treatment.O),
  size = 1, alpha = 0.4, position = position_dodge(width = 0.75)) + 
  theme_bw() + facet_wrap(~RSE.Status.O, labeller = as_labeller(facet_names_new)) + 
  scale_colour_manual(values = c(viridis::viridis(n = 6)[2], "mediumaquamarine"),
  name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
  ylab("Maximum Eye Temperature (°C)") + my.theme_2 + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
  strip.text.x = element_text(size = 12, colour = "black"), 
  panel.grid.major = element_blank(), legend.position = "bottom", 
  axis.title.x = element_blank()) + 
  geom_text(data = Stars, label = "*", size = 10, colour = "black")

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/EM_Boxplot_Low_New.pdf", EM_Boxplot, dpi = 2000,
  height = 6.67, width = 10, units = "in")

# Conducting heat transfer calculations using hourly means. Constructing heat transfer function to begin.

{
dim = 0.011
Area = ((1.1/2)*(1.0/2)*pi)*0.0001
Area = Area*2 # For two eyes 

qrad_Clean = c()
for (i in 1:nrow(hour_means)) {
  qrad_Clean[i] = Area * (5.67 * 10^-8) * 0.95 * 0.95 *
    ((hour_means$Mean.Eye[i] + 273.15)^4 - (hour_means$Mean.Amb[i] + 273.15)^4)
}

Kt_clean = c()
v_clean = c()
alpha_clean = c()
Gr_clean = c()
Nu_clean = c()
Hc_clean = c()
qconv_clean = c()

for (i in 1:nrow(hour_means)) {
  Kt_clean[i] = airtconductivity(hour_means$Mean.Amb[i])
  v_clean[i] = kinematic.viscosity(hour_means$Mean.Amb[i], 101.325)
  alpha_clean[i] = 1 / (hour_means$Mean.Amb[i] + 273.15)
}

for (i in 1:nrow(hour_means)) {
  Gr_clean[i] = ((alpha_clean[i]) * (9.81) * (dim)^3 *
                   (hour_means$Mean.Eye[i] - hour_means$Mean.Amb[i])) / (v_clean[i])^2
}

for (i in 1:nrow(hour_means)) {
  Nu_clean[i] = sign(Gr_clean[i]) * 0.50 * (abs(Gr_clean[i]))^0.25
}

for (i in 1:nrow(hour_means)) {
  Hc_clean[i] = Nu_clean[i] * (Kt_clean[i] / dim)
}

for (i in 1:nrow(hour_means)) {
  qconv_clean[i] = Area * (Hc_clean[i]) * (hour_means$Mean.Eye[i] - hour_means$Mean.Amb[i])
}

qtot_clean = c()

for (i in 1:nrow(hour_means)) {
  qtot_clean[i] = qconv_clean[i] + qrad_Clean[i]
}

hour_means$qtot = qtot_clean
}

# Converting W to mW and extending to two eye.

hour_means$mW = hour_means$qtot*1000

# Assessing means per category

hour_means %>% filter(!is.na(Bird.ID)) %>% group_by(Treatment.O) %>%
  summarise(Mean_mW = mean(mW, na.rm = T),
    SD_mW = sd(mW, na.rm = T))

hour_means %>% filter(!is.na(Bird.ID)) %>% group_by(RSE.Status.O) %>%
  summarise(Mean_mW = mean(mW, na.rm = T),
    SD_mW = sd(mW, na.rm = T))
    
hour_means %>% filter(!is.na(Bird.ID)) %>% group_by(Sex) %>%
  summarise(Mean_mW = mean(mW, na.rm = T),
    SD_mW = sd(mW, na.rm = T))

# Testing evidence for polynomial relationship between ambient temperature and heat transfer

summary(lm(mW ~ poly(Mean.Amb, 2), data = hour_means))
ggplot(hour_means, aes(x = Mean.Amb, y = mW)) + geom_point(pch = 21, size = 3, alpha = 0.5, fill = "magenta") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, colour = "black", linetype = "solid") +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Within this ambient temperature range, a linear trend appears sufficient. Proceeding with analysis

hour_means[c(which(hour_means$mW <= 0)),]

ggplot(hour_means, aes(x = Mean.Amb, y = mW, linetype = Treatment)) + geom_point(pch = 21, size = 3, alpha = 0.5, fill = "magenta") +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, colour = "black") +
  theme_bw()

ggplot(hour_means, aes(mW)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
      aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
      args = list(mean = mean(hour_means$mW, na.rm = T), sd = sd(hour_means$mW, na.rm = T))) +
  theme_bw() + my.theme + xlab("Normalized Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(hour_means$mW, na.rm = T) - 3*sd(hour_means$mW, na.rm = T)),
      size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(hour_means$mW, na.rm = T) + 3*sd(hour_means$mW, na.rm = T)),
      size = 1, linetype = "dashed")

head(hour_means)

# Ensuring data is organised for mixed effects model with auroregressive 
# correlation structure, as predicted given necessity in surface temperature model.

str(hour_means)

HT_Gamm <- mgcv::gamm(mW ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1), 
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  method = "REML", na.action = na.exclude, data = hour_means)
)

# Checking model residuals.

{
p1 = ggplot(data = data.frame("Res" = residuals(HT_Gamm$lme, type = "normalized")),
  aes(x = 1:length(residuals(HT_Gamm$lme, type = "normalized")), y = Res)) + my.theme + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalized Residuals")

p2 = ggplot(data = data.frame("Res" = residuals(HT_Gamm$lme, type = "normalized"), 
  "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))),
  aes(x = Fit, y = Res)) + my.theme + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalized Residuals")

p3 = ggplot(data.frame("Res" = residuals(HT_Gamm$lme, type = "normalized")), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(residuals(HT_Gamm$lme, type = "normalized")), 
    sd = sd(residuals(HT_Gamm$lme, type = "normalized")))) +
  theme_bw() + my.theme + xlab("Normalized Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(residuals(HT_Gamm$lme, type = "normalized")) - 
    3*sd(residuals(HT_Gamm$lme, type = "normalized"))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(residuals(HT_Gamm$lme, type = "normalized")) + 
    3*sd(residuals(HT_Gamm$lme, type = "normalized"))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = residuals(HT_Gamm$lme, type = "normalized"), 
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Model Residuals")
}

# Quite similar to surface temperature model - nicely distributed but with slight
# left skew. Assessing patterns by predictors.

{
p1 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Mean.Amb, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "slateblue", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~Treatment)

p2 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Hour, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "aquamarine", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~Treatment)

p3 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Mean.Amb, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "slateblue", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(Treatment~RSE.Status.O)

p4 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Hour, y = Res)) +
    geom_point(size = 3, pch = 21, fill = "aquamarine", colour = "black", alpha = 0.5) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(Treatment~RSE.Status.O)

p5 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Treatment, y = Res)) +
    geom_boxplot(fill = "navajowhite", colour = "black") + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~RSE.Status.O)

p6 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Pen_Short, y = Res)) +
    geom_boxplot(fill = "navajowhite", colour = "black") + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(~RSE.Status.O)

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, top = "Model Residuals")
}

# Few observations made late in the day, but the do not appear misplaced. Assessing autocorrelation again.

acf(residuals(HT_Gamm$lme, type = "pearson"))
acf(residuals(HT_Gamm$lme, type = "normalized"))

# Autocorrelation structure corrects to some degree, but residual environmental effects are
# likely remaining at play. Continuing on to predictive capacity and bias or fitted values.

summary(lm(hour_means$mW~predict(HT_Gamm$gam, type = "response", newdata = hour_means)))

HT_Gamm <- mgcv::gamm(mW ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1), 
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  method = "REML", na.action = na.omit, data = hour_means)


Null_Gamm <- mgcv::gamm(mW ~ 1, 
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  data=hour_means)

(sum(as.numeric(residuals(Null_Gamm$gam, type = "deviance"))^2) - 
sum(as.numeric(residuals(HT_Gamm$gam, type = "deviance"))^2))/
(sum(as.numeric(residuals(Null_Gamm$gam, type = "deviance"))^2))
# Highly predictive (beta near zero) and r-squared is quite high (0.8325).

# Re-running model with NAs included for further inspection

HT_Gamm <- mgcv::gamm(mW ~ Treatment.O*RSE.Status.O + Sex +
  s(Mean.Amb, k = 4, bs = "cr") +
  ti(Mean.Amb, Treatment.O, k = 4, bs = c("cr","fs")) +
  ti(Mean.Amb, Treatment.O, k = 4, by = RSE.Status.O, bs = c("cr", "fs"), m = 1) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs")) +
  te(Hour, Direction, k = 4, bs = c("cs", "fs"), by = Treatment.O, m = 1), 
  random = list(Bird.ID = ~1, Date.of.Photo = ~1, Pen_Short = ~1), 
  method = "REML", na.action = na.exclude, data = hour_means)
)

{
p1 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Fit, y = mW)) +
    geom_point(size = 3, fill = "navajowhite", colour = "black", alpha = 0.5, pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = F) + 
    theme_bw() + ylab("mW") + facet_grid(RSE.Status.O~Treatment.O)

p2 = na.omit(hour_means) %>% mutate("Res" = residuals(HT_Gamm$lme, type = "normalized"),
    "Fit" = predict(HT_Gamm$gam, newdata=na.omit(hour_means))) %>%
    ggplot(aes(x = Fit, y = Res)) +
    geom_point(size = 3, fill = "slateblue", colour = "black", alpha = 0.5, pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = F) + 
    theme_bw() + ylab("Normalized Residuals") + facet_grid(RSE.Status.O~Treatment.O)

grid.arrange(p1, p2, top = "Model Residuals")
}

# No apparent bias, though a very subtle fanning can be observed in the stress-exposed 
# treatments, regardless of status. This fan is likely too small to be influential, however.

# Checking for collinearity due to inclusion of both sex and status as parametric predictors.

vif_gamm(HT_Gamm)

# Summarising model output

Coefs = data.frame("Predictors" = names(HT_Gamm$gam$coefficients),
  "Coefficients" = as.numeric(HT_Gamm$gam$coefficients),
  "E.D.F." = as.numeric(HT_Gamm$gam$edf),
  "S.E.M" = as.numeric(sqrt(diag(vcov(HT_Gamm$gam, unconditional = TRUE))))
  )

print(Coefs)

df.residual(HT_Gamm$gam) # 1004.824

# Averages summaries across splines

Coefs$Predictor_Combined = gsub("\\\\.[[:digit:]].*", "", Coefs[,1])
Spline_Sum = Coefs %>% group_by(Predictor_Combined) %>%
  summarise(Mean.Beta = mean(Coefficients, na.rm = T),
    Mean.SE = mean(S.E.M, na.rm = T),
    edf = sum(E.D.F.))

print(Spline_Sum)

summary(HT_Gamm$gam)

{
    Sub_TNZ =  ref_grid(HT_Gamm,  
    at = list(Mean.Amb = seq(2.5, 13.99, by = I((13.99-2.5)/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

    TNZ =  ref_grid(HT_Gamm,  
    at = list(Mean.Amb = seq(14, 30, by = I((30-14)/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

    Supra_TNZ =  ref_grid(HT_Gamm,  
    at = list(Mean.Amb = seq(30.01, 38.5, by = I((38.5-30.01)/9))), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

    Sub_TNZ_Cont = emmeans(Sub_TNZ, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

    TNZ_Cont = emmeans(TNZ, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

    Supra_TNZ_Cont = emmeans(Supra_TNZ, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)
}

print(Sub_TNZ_Cont)
print(TNZ_Cont)
print(Supra_TNZ_Cont)

# Results are highly similar to those derived from surface temperature modeling
# (as expected). 

# Assessing means at minimum and maximum observed temperatures.

Max_Temp = max(Clean_Ranks$Amb.Temp, na.rm = T)
Min_Temp = min(Clean_Ranks$Amb.Temp, na.rm = T)

EM_Grid_Lower_Limit = ref_grid(HT_Gamm,  
    at = list(Mean.Amb = Min_Temp), 
    cov.reduce = mean, type = "response",
    nesting = NULL, nesting.order = FALSE,
    data = hour_means)

EM_Grid_Upper_Limit = ref_grid(HT_Gamm,  
    at = list(Mean.Amb = Max_Temp), 
    cov.reduce = mean, type = "response",
    nesting = NULL, nesting.order = FALSE,
    data = hour_means)

emmeans(EM_Grid_Lower_Limit, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

emmeans(EM_Grid_Upper_Limit, specs = pairwise ~ Treatment.O | RSE.Status.O, 
    cov.reduce = mean, type = "response",
    p.adjust.method = "bonferroni",
    data = hour_means)

# Plotting full model.

HT_EMM_Grid = ref_grid(HT_Gamm,  
    at = list(Mean.Amb = seq(2.5, 38.5, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = hour_means)

emmip(HT_EMM_Grid, Treatment.O ~ Mean.Amb | RSE.Status.O,
  data = hour_means, type = "response", CIs = T, 
  cov.reduce = FALSE, plot = T)

HT_EMM = emmip(HT_EMM_Grid, Treatment.O ~ Mean.Amb | RSE.Status.O,
  data = hour_means, type = "response", CIs = T, 
  cov.reduce = FALSE, plot = F)

#write.csv(HT_EMM, "/home/joshk/git_repositories/BCCH_Dominance/Data/Heat_T_EMM.csv")
HT_EMM = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Heat_T_EMM.csv")

# First plotting boxplot

HT_EMM$Zone = "Thermoneutrality"
HT_EMM$Zone[c(which(HT_EMM$Mean.Amb < 14))] = "Below\\nThermoneutrality"
HT_EMM$Zone[c(which(HT_EMM$Mean > 30))] = "Above\\nThermoneutrality"
HT_EMM$Zone = factor(HT_EMM$Zone, levels = c("Below\\nThermoneutrality", 
"Thermoneutrality", "Above\\nThermoneutrality"))

facet_names_new <- c(`1` = "Dominant",
                    `2` = "Subordinate"
                    )

HT_EM_Boxplot = HT_EMM %>% rename("mW" = yvar) %>%
ggplot(aes(x = Zone, y = mW, colour = Treatment.O)) + 
  geom_boxplot(size = 0.75, position = position_dodge(width = 0.75)) + 
  geom_jitter(data = as.data.frame(hour_means %>% filter(!is.na(RSE.Status.O))), 
  aes(x = Zone, y = mW, colour = Treatment.O),
  size = 1, alpha = 0.4, position = position_dodge(width = 0.75)) + 
  theme_bw() + facet_wrap(~RSE.Status.O, labeller = as_labeller(facet_names_new)) + 
  scale_colour_manual(values = c(viridis::viridis(n = 6)[2], "mediumaquamarine"),
  name = "Treatment") + ylab("Heat Transfer Rate (mW)") + my.theme_2 + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
  strip.text.x = element_text(size = 12, colour = "black"), 
  panel.grid.major = element_blank(), legend.position = "bottom", 
  axis.title.x = element_blank())

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/HT_EM_Boxplot_Low_New.pdf", HT_EM_Boxplot, dpi = 800,
  height = 6.67, width = 10, units = "in")

# Looping through confidence intervals for graded bands; it appears graded ribbons are not possible to construct otherwise (at least simply).

HT_EMM = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Heat_T_EMM.csv")
vir_select = viridis::viridis(n = 8)[4]

RRy = max(subset(HT_EMM, RSE.Status.O == "1" & Treatment.O == "Rest")$yvar)
SRy = max(subset(HT_EMM, RSE.Status.O == "2" & Treatment.O == "Rest")$yvar)
RSy = max(subset(HT_EMM, RSE.Status.O == "1" & Treatment.O == "Stress")$yvar)
SSy = max(subset(HT_EMM, RSE.Status.O == "2" & Treatment.O == "Stress")$yvar)
x1 = min(HT_EMM$Mean.Amb) + 1

Sub_seg_dat1 = data.frame(x=x1,y=SRy,xend=11,yend=17,
                      RSE.Status.O="2")
Sub_seg_dat2 = data.frame(x=x1,y=SSy,xend=11,yend=17,
                      RSE.Status.O="2")
Dom_seg_dat1 = data.frame(x=x1,y=RRy,xend=11,yend=17,
                      RSE.Status.O="1")
Dom_seg_dat2 = data.frame(x=x1,y=RSy,xend=11,yend=17,
                      RSE.Status.O="1")

HT_EMM$Treatment = "Control"
HT_EMM$Treatment[c(which(HT_EMM$Treatment.O == "Stress"))] = "Stress-Exposed"
HT_EMM$Treatment = factor(HT_EMM$Treatment)

Fade_plot = ggplot(HT_EMM, aes(x = Mean.Amb, y = yvar, linetype = Treatment,
    fill = Treatment)) + 
    geom_ribbon(aes(x = Mean.Amb, ymin = LCL, ymax = UCL),
    colour = NA, alpha = 0.5) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"),
      size = 1, colour = "black", se = FALSE) +
  annotate("rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf,
        fill = "grey20", alpha = 0.5) +
  theme_bw() + my.theme_2 + xlab("Ambient Temperature (°C)") +
  ylab("Heat Transfer Rate (mW)") +
  annotate("text", x = 22, y = 52, label = "TNZ", size = 8) +
  theme(strip.text.x = element_text(size = 16),
    strip.background = element_rect(fill="gray90")) +
  geom_segment(data = Sub_seg_dat1, aes(x=x,y=y,yend=66,xend=xend),
    inherit.aes=FALSE, size=1, color="black", lineend = "round",
    linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
    alpha = 0.7) +
  geom_segment(data = Sub_seg_dat2, aes(x=x,y=y,yend=66,xend=xend),
    inherit.aes=FALSE, size=1, color="black", lineend = "round",
    linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
    alpha = 0.7) +
  geom_segment(data = Dom_seg_dat1, aes(x=x,y=y,yend=66,xend=xend),
      inherit.aes=FALSE, size=1, color="black", lineend = "round",
      linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
      alpha = 0.7) +
   geom_segment(data = Dom_seg_dat2, aes(x=x,y=y,yend=66,xend=xend),
      inherit.aes=FALSE, size=1, color="black", lineend = "round",
      linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
      alpha = 0.7) +
   annotate("text", x = 12.4, y = 67, label = expression(Delta), size = 8,
    parse = T,
    alpha = 0.8) + scale_fill_manual(values = c(viridis::viridis(n = 6)[2], "mediumaquamarine")) + 
  facet_wrap(~RSE.Status.O, labeller = as_labeller(facet_names_new))

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Heat_Fade_Low_New.pdf", Fade_plot, dpi = 800,
  width = 9.67, height = 6.67, units = "in")

# Plotting with raw data points

hour_means$Treatment = "Control"
hour_means$Treatment[c(which(hour_means$Treatment.O == "Stress"))] = "Stress-Exposed"
hour_means$Treatment = factor(hour_means$Treatment)

Fade_plot_points = ggplot(HT_EMM, aes(x = Mean.Amb, y = yvar, 
    linetype = Treatment, fill = Treatment)) + 
  geom_ribbon(aes(x = Mean.Amb, ymin = LCL, ymax = UCL),
    colour = NA, alpha = 0.5) + 
  geom_point(data = as.data.frame(hour_means %>% 
    filter(!is.na(RSE.Status.O))), aes(x = Mean.Amb, y = mW),
  pch = 21, size = 2, alpha = 0.3, colour = "black") + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4, bs = "cs"),
      size = 1, colour = "black", se = FALSE) +
  annotate("rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf,
        fill = "grey20", alpha = 0.5) +
  theme_bw() + my.theme_2 + xlab("Ambient Temperature (°C)") +
  ylab("Heat Transfer Rate (mW)") +
  annotate("text", x = 22, y = 52, label = "TNZ", size = 8) +
  theme(strip.text.x = element_text(size = 16),
    strip.background = element_rect(fill="gray90")) +
  geom_segment(data = Sub_seg_dat1, aes(x=x,y=y,yend=66,xend=xend),
    inherit.aes=FALSE, size=1, color="black", lineend = "round",
    linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
    alpha = 0.7) +
  geom_segment(data = Sub_seg_dat2, aes(x=x,y=y,yend=66,xend=xend),
    inherit.aes=FALSE, size=1, color="black", lineend = "round",
    linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
    alpha = 0.7) +
  geom_segment(data = Dom_seg_dat1, aes(x=x,y=y,yend=66,xend=xend),
      inherit.aes=FALSE, size=1, color="black", lineend = "round",
      linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
      alpha = 0.7) +
   geom_segment(data = Dom_seg_dat2, aes(x=x,y=y,yend=66,xend=xend),
      inherit.aes=FALSE, size=1, color="black", lineend = "round",
      linejoin = "round", arrow=arrow(ends = "first", length = unit(0.3, "cm")),
      alpha = 0.7) +
   annotate("text", x = 12.4, y = 67, label = "Delta", size = 8,
    parse = T,
    alpha = 0.8) + scale_fill_manual(values = c(viridis::viridis(n = 6)[2], "mediumaquamarine")) + 
  facet_wrap(~RSE.Status.O, labeller = as_labeller(facet_names_new))

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Heat_Fade_Low_New_points.pdf", 
    Fade_plot_points, dpi = 800, width = 9.67, height = 6.67, units = "in")

# Integrating across difference curves (that is, the difference in heat transfer between stress-exposed
# and control treatments, within a given social status).

library("pracma")
Areas_It = c()
Areas_All = vector('list',2)

for (k in 1:2){
  B_sample = subset(HT_EMM, RSE.Status.O == k)
  for (j in 1:1000){
    curve_est_Stress = c()
    curve_est_Rest = c()
    for (i in 1:length(unique(B_sample$Mean.Amb))){
      curve_est_Stress[i] = rnorm(n = 1, mean = subset(B_sample, Treatment.O == "Stress")$yvar[i], 
        sd = subset(B_sample, Treatment.O == "Stress")$SE[i])
      curve_est_Rest[i] = rnorm(n = 1, mean = subset(B_sample, Treatment.O == "Rest")$yvar[i], 
        sd = subset(B_sample, Treatment.O == "Rest")$SE[i])
    }
    f1 = approxfun(unique(B_sample$Mean.Amb), curve_est_Stress-curve_est_Rest)
    f2 = function(x) abs(f1(x))
    Area_Est_Low = integral(f2, 2.5, 14)
    Area_Est_High = integral(f2, 30, 38.5)
    Areas_It[j] = Area_Est_Low + Area_Est_High
  }
  Areas_All[[k]] = data.frame("Area" = Areas_It, "Status" = paste(k))
}

Areas_DF = bind_rows(Areas_All)
Areas_DF$RSE.Status.O = "Dominant"
Areas_DF$RSE.Status.O[c(which(Areas_DF$Status == 2))] = "Subordinate"

Areas_DF %>% group_by(RSE.Status.O) %>%
  dplyr::summarize(Mean = mean(Area), SD = sd(Area))

t.test(Areas_DF$Area ~ Areas_DF$RSE.Status.O)

# Ensuring directionality is correct

Areas_It_1 = c()
Areas_All_1 = vector('list',2)
Areas_It_2 = c()
Areas_All_2 = vector('list',2)
Curve_store = vector('list',1000)

for (k in 1:2){
  B_sample = subset(HT_EMM, RSE.Status.O == k)
  for (j in 1:1000){
    curve_est_Stress = c()
    curve_est_Rest = c()
    for (i in 1:length(unique(B_sample$Mean.Amb))){
      curve_est_Stress[i] = rnorm(n = 1, mean = subset(B_sample, Treatment.O == "Stress")$yvar[i], 
        sd = subset(B_sample, Treatment.O == "Stress")$SE[i])
      curve_est_Rest[i] = rnorm(n = 1, mean = subset(B_sample, Treatment.O == "Rest")$yvar[i], 
        sd = subset(B_sample, Treatment.O == "Rest")$SE[i])
    }
    #f1 = approxfun(B_sample$Mean.Amb, curve_est_Rest-curve_est_Stress)
    Temp_Range = unique(B_sample$Mean.Amb)
    mod = gam(I(curve_est_Rest-curve_est_Stress) ~ s(Temp_Range, k = 4), method = "REML")
    f1 = function(x) {
      data = data.frame("Temp_Range" = x)
      predict.gam(mod, newdata = data)
    }
    Curve_store[[j]] = f1
    f2 = function(x) abs(f1(x))
    Area_Est_Low_1 = integral(f1, 2.5, 14)
    Area_Est_High_1 = integral(f1, 30, 38.5)*-1
    Areas_It_1[j] = Area_Est_Low_1 + Area_Est_High_1
    Area_Est_Low_2 = integral(f2, 2.5, 14)
    Area_Est_High_2 = integral(f2, 30, 38.5)
    Areas_It_2[j] = Area_Est_Low + Area_Est_High
  }
  Areas_All_1[[k]] = data.frame("Area" = Areas_It_1, "Status" = paste(k))
  Areas_All_2[[k]] = data.frame("Area" = Areas_It_2, "Status" = paste(k))
}

Areas_DF_1 = bind_rows(Areas_All_1)
Areas_DF_1$RSE.Status.O = "Dominant"
Areas_DF_1$RSE.Status.O[c(which(Areas_DF_1$Status == 2))] = "Subordinate"
Areas_DF_2 = bind_rows(Areas_All_2)
Areas_DF_2$RSE.Status.O = "Dominant"
Areas_DF_2$RSE.Status.O[c(which(Areas_DF_2$Status == 2))] = "Subordinate"

# Plotting mean curve

Curve_dat = vector('list', 2)
for (i in 1:2){
  Curve_dat[[i]] = vector('list', 1000)
}

for (k in 1:2){
  B_sample = subset(HT_EMM, RSE.Status.O == k)

  for (j in 1:1000){
    curve_est_Stress = c()
    curve_est_Rest = c()

    for (i in 1:length(B_sample$Mean.Amb)){
        curve_est_Stress[i] = rnorm(n = 1, mean = subset(B_sample, Treatment.O == "Stress")$yvar[i], 
        sd = subset(B_sample, Treatment.O == "Stress")$SE[i])
        curve_est_Rest[i] = rnorm(n = 1, mean = subset(B_sample, Treatment.O == "Rest")$yvar[i], 
        sd = subset(B_sample, Treatment.O == "Rest")$SE[i])
      }
      mod_frame_low = data.frame("Mean.Amb" = unique(B_sample$Mean.Amb), "Diff" = c(curve_est_Rest-curve_est_Stress))
      mod_frame_high = data.frame("Mean.Amb" = unique(B_sample$Mean.Amb), "Diff" = c(curve_est_Stress-curve_est_Rest))
      mod_frame_low = subset(mod_frame_low, Mean.Amb <= 29.5)
      mod_frame_high = subset(mod_frame_high, Mean.Amb >= 14.5)
      mod_frame = rbind(mod_frame_low, mod_frame_high)

    mod = gam(Diff ~ s(Mean.Amb, k = 4), data = mod_frame, method = "REML")
    dat = data.frame("Mean.Amb" = unique(B_sample$Mean.Amb))
    Curve_dat[[k]][[j]] = data.frame("Mean.Amb" = unique(B_sample$Mean.Amb), 
        "Pred" = predict.gam(mod, newdata = dat))
    }
}

for (i in 1:2){
  Curve_dat[[i]] = bind_rows(Curve_dat[[i]])
  Curve_dat[[i]] = Curve_dat[[i]] %>% group_by(Mean.Amb) %>%
  dplyr::summarize(Mean = mean(Pred), SE = sd(Pred)/sqrt(20), UCL = Mean + 1.96*SE, LCL = Mean - 1.96*SE)
}

Curve_dat[[1]]$Status = "Dominant"
Curve_dat[[2]]$Status = "Subordinate"
Curve_dat = bind_rows(Curve_dat)
Curve_dat = as.data.frame(Curve_dat)
vir_select = viridis::viridis(n = 8)[4]

# Creating fills for each plot, again with absurd looping through confidence bands

p_func_dom = with(subset(Curve_dat, Status == "Dominant"), approxfun(Mean.Amb, Mean))
p_func_sub = with(subset(Curve_dat, Status == "Subordinate"), approxfun(Mean.Amb, Mean))
T_Range = seq(min(Curve_dat$Mean.Amb), max(Curve_dat$Mean.Amb), by = 0.01)

Poly_All = rbind(data.frame("X" = T_Range,
  "Y" = p_func_dom(T_Range), "Status" = "Dominant"),
    data.frame("X" = T_Range,
    "Y" = p_func_sub(T_Range), "Status" = "Subordinate")
)

Poly_All$Ymin = 0
Poly_All$Ymax = 0
Poly_All$Ymin[c(which(Poly_All$Y<=0))] = Poly_All$Y[c(which(Poly_All$Y<=0))]
Poly_All$Ymax[c(which(Poly_All$Y>0))] = Poly_All$Y[c(which(Poly_All$Y>0))]
Poly_All$Status = factor(Poly_All$Status)
Poly_All$Mean = 0

Text_pos = data.frame("X" = c(10,10), "Y" = c(0.5, 1.5), "Status" = c("Dominant", "Subordinate"))
Text_pos$Status = factor(Text_pos$Status)

Curve_plot = ggplot(Curve_dat, aes(x = Mean.Amb, y = Mean, group = Status)) +
  geom_ribbon(data = Poly_All, mapping = aes(x = X, ymin = Ymin, ymax = Ymax),
    fill = "grey60", alpha = 0.5) +
  geom_ribbon(mapping = aes(x = Curve_dat$Mean.Amb, ymin = LCL, ymax = UCL), 
    alpha = 0.5, colour = NA, fill = vir_select) + 
    geom_smooth(mapping = aes(x = Curve_dat$Mean.Amb, y = Mean),
      colour = "black", size = 1, se = FALSE) +
    annotate("rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf,
        fill = "grey20", alpha = 0.5) +
    geom_text(data = Text_pos, aes(x = X, y = Y), colour = "black", size = 6,
      family = "Noto Sans", label = "h(x)") +
      annotate("text", x = 22, y = -1.25, label = "TNZ", size = 5, 
        family = "Noto Sans", colour = "black") +
      geom_segment(mapping=aes(x=14.5, y=-0.75, xend=29.5, yend=-0.75),
                   arrow=arrow(ends = "both"), size=1.25, color="black",
                   lineend = "round", linejoin = "round") +
        theme_bw() + my.theme + xlab("Ambient Temperature (°C)") +
        ylab(expression(Delta*" Heat Transfer (mW)")) +
    theme(strip.text.x = element_text(size = 16, family = "Noto Sans"),
      strip.background = element_rect(fill="gray90")) +
  facet_wrap(~Status)

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Curve_Diff_Plot_New.jpeg", Curve_plot, dpi = 800,
  width = 9.67, height = 6.67, units = "in")

# Plotting sums across sub and supra-thermoneutrality; non-absolute values.

Heat_dotplot_NA = Areas_DF_1 %>% group_by(RSE.Status.O) %>%
  dplyr::summarize(Mean = mean(Area),
    SE = sd(Area)/sqrt(10),
    UCL = Mean + 1.96*SE,
    LCL = Mean - 1.96*SE) %>%
    ggplot(aes(x = RSE.Status.O, y = Mean, fill = RSE.Status.O)) +
      geom_errorbar(mapping = aes(x = RSE.Status.O, ymin = LCL, ymax = UCL),
      colour = "black", size = 1, width = 0.2) +
      geom_point(size = 5, pch = 21, colour = "black") +
      theme_bw() + my.theme + scale_fill_viridis_d() + ylab(expression(Delta*" Total Heat Transfer (mW x °C)")) +
      theme(axis.title.x = element_blank(), legend.position = "none", axis.ticks.y.right = element_blank()) +
      geom_signif(
        comparisons = list(c("Dominant", "Subordinate")),
        y_position = 85, annotations = c("*"), size = 1,
        textsize = 8
      ) +
      scale_y_continuous(limits = c(-5, 90), sec.axis = dup_axis(~ ./10, name = "Relative Energy Conservation",
                                             labels = c("-", "", "", "", "+")))


ggsave("/home/joshk/git_repositories/BCCH_Dominance/Figures/Heat_DotPlot_New.jpeg", Heat_dotplot_NA, dpi = 800,
  width = 7.67, height = 6.67, units = "in")

# Using absolute values

# Heat_dotplot = Areas_DF %>% group_by(RSE.Status.O) %>%
#   dplyr::summarize(Mean = mean(Area),
#     SE = sd(Area)/sqrt(10),
#     UCL = Mean + 1.96*SE,
#     LCL = Mean - 1.96*SE) %>%
#     ggplot(aes(x = RSE.Status.O, y = Mean, fill = RSE.Status.O)) +
#       geom_errorbar(mapping = aes(x = RSE.Status.O, ymin = LCL, ymax = UCL),
#       colour = "black", size = 1, width = 0.3) +
#       geom_point(size = 5, pch = 21, colour = "black") +
#       theme_bw() + my.theme + scale_fill_viridis_d() + ylab(expression(Delta*" Total Heat Transfer (mW x °C)")) +
#       theme(axis.title.x = element_blank(), legend.position = "none", axis.ticks.y.right = element_blank()) +
#       geom_signif(
#         comparisons = list(c("Dominant", "Subordinate")),
#         y_position = 98, annotations = c("*"), size = 1,
#         textsize = 8
#       ) +
#       scale_y_continuous(limits = c(40, 105), sec.axis = dup_axis(~ ./10, name = "Relative Energy Conservation",
#                                              labels = c("-", "","","+")))


# showtext_opts(dpi = 2000)
# showtext_auto()

# ggsave("C:/Users/joshk/Documents/BCCH_Dominance/BCCH-Dominance-Analyses/Heat_DotPlot.jpeg", Heat_dotplot, dpi = 2000,
#   width = 7.67, height = 6.67, units = "in")

Areas_DF$Status = factor(Areas_DF$Status)

summary(gls(Area ~ Status, weights = varIdent(form = ~1|Status), data = Areas_DF))

# Plotting attribute-ordered networks

testing = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Dominance_Dat5.csv")
summed <- as.data.frame(testing %>% group_by(Pen_Short, Dominant, Subordinate) %>% summarise("Count" = n()))
E_join = read.csv("/home/joshk/git_repositories/BCCH_Dominance/Data/Rank_Summaries.csv")

summed = summed[-c(which(summed$Subordinate == "OOAABlO")),]
E_join$ID = as.character(E_join$Bird.ID)
summed$Dominant = as.character(summed$Dominant)
summed$Subordinate = as.character(summed$Subordinate)

Rank_d = c()
Rank_s = c()

for (i in 1:nrow(summed)){
    Rank_d[i] = E_join$RS.Rank[c(which(E_join$ID == summed$Dominant[i]))]
    Rank_s[i] = E_join$RS.Rank[c(which(E_join$ID == summed$Subordinate[i]))]
}

summed$Dominant_Ranks = Rank_d
summed$Subordinate_Ranks = Rank_s
summed$Rank_Diff = with(summed, Dominant_Ranks - Subordinate_Ranks)

Attrib = summed %>% dplyr::select("Pen" = Pen_Short, "ID" = Dominant, "Rank" = Dominant_Ranks) %>%
  rbind(., as.data.frame(summed %>% dplyr::select("Pen" = Pen_Short, "ID" = Subordinate, "Rank" = Subordinate_Ranks))) %>%
  unique(.) %>% mutate(x = 1) %>% dplyr::select(Pen, ID, x, Rank)

summed = summed %>% dplyr::select("Pen" = Pen_Short, "ID1" = Dominant, "ID2" = Subordinate, "Wins" = Count, Rank_Diff)

At_list = split(Attrib, f = Attrib$Pen)
for (i in 1:length(At_list)){
  At_list[[i]] = At_list[[i]] %>% dplyr::select(-Pen)
}

Edge_list = split(summed, f = summed$Pen)
for (i in 1:length(Edge_list)){
  Edge_list[[i]] = Edge_list[[i]] %>% dplyr::select(-Pen)
}

libraries('igraph','ggraph')

Hier_Plot = vector('list', 4)
HP_Out = vector('list', 4)
att_ord = vector('list', 4)
tiewidth = vector('list', 4)

for (i in 1:4){
  Hier_Plot[[i]] = graph.data.frame(Edge_list[[i]], directed = TRUE, vertices = At_list[[i]])
  E(Hier_Plot[[i]])[Rank_Diff<0]$color = viridis::viridis(n = 2)[1]
  E(Hier_Plot[[i]])[Rank_Diff>0]$color = viridis::viridis(n = 2)[2]

  # Setting axes

  x = V(Hier_Plot[[i]])$x
  y = V(Hier_Plot[[i]])$Rank*-1
  att_ord[[i]] = cbind(x, y)

  # Setting node sized and plotting

  nodesize = 10
  tiewidth[[i]] = E(Hier_Plot[[i]])$Wins/2
}

tiff(filename="C:/Users/joshk/Documents/BCCH_Dominance/BCCH-Dominance-Analyses/Figures/NE_Hier.tiff", width = 5, height = 7, units = "in", res = 2000)

plot(Hier_Plot[[1]], layout = att_ord[[1]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[1]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[1]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

dev.off()

tiff(filename="C:/Users/joshk/Documents/BCCH_Dominance/BCCH-Dominance-Analyses/Figures/NW_Hier.tiff", width = 5, height = 7, units = "in", res = 2000)

plot(Hier_Plot[[2]], layout = att_ord[[2]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[2]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[2]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

dev.off()

tiff(filename="C:/Users/joshk/Documents/BCCH_Dominance/BCCH-Dominance-Analyses/Figures/SE_Hier.tiff", width = 5, height = 7, units = "in", res = 2000)

plot(Hier_Plot[[3]], layout = att_ord[[3]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[3]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[3]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

dev.off()

tiff(filename="C:/Users/joshk/Documents/BCCH_Dominance/BCCH-Dominance-Analyses/Figures/SW_Hier.tiff", width = 5, height = 7, units = "in", res = 2000)

plot(Hier_Plot[[4]], layout = att_ord[[4]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[4]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[4]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

dev.off()

# Combined

pdf("/home/joshk/Desktop/Joined_Plot.pdf", width = 800, height = 1000)

op = par(mfrow = c(2, 2), mai = c(0.05, 0.05, 0.05, 0.05))

plot(Hier_Plot[[1]], layout = att_ord[[1]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[1]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[1]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

plot(Hier_Plot[[2]], layout = att_ord[[2]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[2]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[2]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

plot(Hier_Plot[[3]], layout = att_ord[[3]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[3]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[3]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

plot(Hier_Plot[[4]], layout = att_ord[[4]], vertex.size = nodesize,
     vertex.label = NA,
     vertex.color = c(viridis::viridis(n = 2)[1],viridis::viridis(n = 2)[2],
                      viridis::viridis(n = 2)[2],viridis::viridis(n = 2)[1],
                      viridis::viridis(n = 2)[2]),
     edge.width = tiewidth[[4]],
     edge.arrow.size = 0.01,
     edge.curved = TRUE,
     edge.loop.angle = 1,
     edge.color = E(Hier_Plot[[4]])$color, asp = 0,
     axes=FALSE, xlim = c(-1.5,-0.5))

par(op)
dev.off()

# In ggplot 

extractLegend <- function(xplot) {
  grobs <- ggplot_gtable(ggplot_build(xplot))
  g_title <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[g_title]]
}

{
weights_df = vector("list", 4)
layout_df = vector("list", 4)
out_plots = vector("list", 4)

for (i in 1:4){
  layout_df[[i]] = as.data.frame(att_ord[[i]])
  layout_df[[i]]$IDs = At_list[[i]]$ID  
  weights_df[[i]] = get.data.frame(Hier_Plot[[i]]) 
  weights_df[[i]] = weights_df[[i]] %>% 
    mutate(from.x = layout_df[[i]]$x[match(weights_df[[i]]$from, layout_df[[i]]$IDs)],
    from.y = layout_df[[i]]$y[match(weights_df[[i]]$from, layout_df[[i]]$IDs)],
    to.x = layout_df[[i]]$x[match(weights_df[[i]]$to, layout_df[[i]]$IDs)],
    to.y = layout_df[[i]]$y[match(weights_df[[i]]$to, layout_df[[i]]$IDs)])
  weights_df[[i]]$Direction = ifelse(weights_df[[i]]$to.y > weights_df[[i]]$from.y, "Up", "Down")
  weights_df[[i]]$Direction = factor(weights_df[[i]]$Direction, levels = c("Up", "Down"))

  layout_df[[i]]$Status = ifelse(layout_df[[i]]$y %in% c(-1,-2), "Dominant", "Subordinate")

  out_plots[[i]] = ggplot() +
    geom_curve(data = weights_df[[i]], aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour = Direction), size = weights_df[[i]]$Wins/8) +
    geom_point(data = layout_df[[i]], aes(x = x, y = y, fill = Status), pch = 21, size = 8, colour = "black") +
    #geom_text(data=fr.all.df,aes(x=V1,y=V2,label=species)) + # add the node labels
    scale_x_continuous(expand=c(0,1))+  # expand the x limits 
    scale_y_continuous(expand=c(0,1))+ # expand the y limits
    scale_fill_viridis_d() + scale_colour_viridis_d() +  
    theme_void() + 
    theme(legend.position = "none", 
    plot.margin = rep(unit(0,"null"),4),
    panel.margin = unit(0,"null"))
}
    legend_plot = ggplot() +
    geom_point(data = layout_df[[1]], aes(x = x, y = y, fill = Status), pch = 21, size = 8, colour = "black") +
    scale_fill_viridis_d() +
    theme(legend.key = element_rect(fill = "white"))

    Leg = extractLegend(legend_plot)
    
    b_plot = ggplot(layout_df[[1]]) +
    theme_void() + 
    theme(legend.position = "none", 
    plot.margin = rep(unit(0,"null"),4),
    panel.margin = unit(0,"null"), panel.background = element_rect(fill = "white"))
}

lay = rbind(c(1,1,1,2,2,2,5,5),
            c(1,1,1,2,2,2,5,5),
            c(1,1,1,2,2,2,5,5),
            c(1,1,1,2,2,2,5,5),
            c(1,1,1,2,2,2,6,6),  
            c(3,3,3,4,4,4,6,6), 
            c(3,3,3,4,4,4,5,5),  
            c(3,3,3,4,4,4,5,5),  
            c(3,3,3,4,4,4,5,5),  
            c(3,3,3,4,4,4,5,5))

Grouped_Plot = arrangeGrob(grobs = list(out_plots[[1]], out_plots[[2]], out_plots[[3]], out_plots[[4]], b_plot, Leg), layout_matrix = lay)

Grouped_Plot = gtable_add_grob(Grouped_Plot, grobs = rectGrob(gp = gpar(lwd = 10, col = "black", fill = NA)), t = 1, b = 10, l = 1, r = 6)

ggsave("/home/joshk/Desktop/Text.pdf", width = 8, height = 10, Grouped_Plot, dpi = 1200)


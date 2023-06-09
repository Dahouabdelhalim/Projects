library(ggplot2)
library(pander)
library(dplyr)
library(qdapDictionaries)

# Retrieve data
data <- read.csv("revised_data.csv",
                 stringsAsFactors = F) # 22339
data$text <- paste(data$TI,
                   data$AB,
                   sep = ". ")
write.csv(data, 
          "revised_data.csv",
          row.names = F)
rm(list = ls())

# Subject
data <- read.csv("revised_data.csv",
                 stringsAsFactors = F)
subject.table <- data.frame(id = character(0),
                            subject = character(0))
for (i in 1:nrow(data)) {
  subject <- strsplit(data$WC[i], "; ")[[1]]
  subject.df <- data.frame(id = rep(data$UT[i], length(subject)),
                           subject = subject,
                           stringsAsFactors = F)
  subject.table <- rbind(subject.table, subject.df)
}
write.csv(subject.table, "subject.table.csv", row.names = F) # 46556
rm(list = ls())

# Abstract
abstract <- read.csv("abstract.sample.new.csv",
                     stringsAsFactors = F) # 22896
data <- read.csv("revised_data.csv",
                 stringsAsFactors = F)
data <- data[data$UT %in% abstract$id,] # 17184
write.csv(data, "data_final.csv", row.names = F)
rm(list = ls())

# Term table
data <- read.csv("data_final.csv",
                 stringsAsFactors = F)
term <- read.csv("nlp_result_final_dedup.csv",
                 stringsAsFactors = F) # 44485
for (i in 1:nrow(term)) {
  term$relation.1[i] <- strsplit(term$relation[i], ":")[[1]][1]
  term$POS.1[i] <- substr(term$POS[i], 1, 2)
}

term.1 <- data.frame(id = character(0),
                     sentence = numeric(0),
                     term = character(0),
                     lemma = character(0),
                     governor = character(0),
                     relation = character(0),
                     target = character(0),
                     POS = character(0),
                     NER = character(0),
                     direct = numeric(0),
                     relation.1 = character(0),
                     POS.1 = character(0))
for (i in 1:nrow(term)) {
  if (grepl("\\\\/", term$target[i])) {
    target = strsplit(term$target[i], "/")[[1]]
    target = target[target != "ArtificialIntelligence"]
    term$target[i] <- "ArtificialIntelligence"
    term.df <- rbind(term[i,],
                     data.frame(id = rep(term$id[i], length(target)),
                                sentence = rep(term$sentence[i], length(target)),
                                term = target,
                                lemma = target,
                                governor = rep("governor", length(target)),
                                relation = rep("conj:and", length(target)),
                                target = rep("ArtificialIntelligence", length(target)),
                                POS = rep("NN", length(target)),
                                NER = rep("O", length(target)),
                                direct = rep(term$direct[i], length(target)),
                                relation.1 = rep("conj", length(target)),
                                POS.1 =rep("NN", length(target))))
    term.1 <- rbind(term.1, term.df)
  } else {
    term.1 <- rbind(term.1, term[i,])
  }
}
row.names(term.1) <- 1:nrow(term.1)

term.2 <- data.frame(id = character(0),
                     sentence = numeric(0),
                     term = character(0),
                     lemma = character(0),
                     governor = character(0),
                     relation = character(0),
                     target = character(0),
                     POS = character(0),
                     NER = character(0),
                     relation.1 = character(0),
                     POS.1 = character(0))
for (i in 1:nrow(term.1)) {
  if (grepl("\\\\/", term.1$lemma[i])) {
    target = strsplit(term.1$lemma[i], "/")[[1]]
    if (sum(nchar(target) > 2) ==  length(target)) {
      term.1$term[i] <- target[1]
      term.1$lemma[i] <- target[1]
      target <- target[-1]
      term.df <- rbind(term.1[i,],
                       data.frame(id = rep(term.1$id[i], length(target)),
                                  sentence = rep(term.1$sentence[i], length(target)),
                                  term = target,
                                  lemma = target,
                                  governor = rep("governor", length(target)),
                                  relation = rep(term.1$relation[i], length(target)),
                                  target = rep(term.1$target[i], length(target)),
                                  POS = rep(term.1$POS[i], length(target)),
                                  NER = rep("O", length(target)),
                                  relation.1 = rep(term.1$relation.1[i], length(target)),
                                  POS.1 =rep(term.1$POS.1[i], length(target))))
      term.2 <- rbind(term.2, term.df)
    } else {
      term.2 <- rbind(term.2, term.1[i,])
    }
  } else {
    term.2 <- rbind(term.2, term.1[i,])
  }
}
row.names(term.2) <- 1:nrow(term.2)

term.3 <- term.2[,-2]
term.3 <- term.3[duplicated(term.3) == F,]
term.4 <- term.3[term.3$POS.1 %in% c("JJ", "RB", "NN", "VB"),] # 31208 words
term.4 <- term.4[term.4$id %in% data$UT,] # 31208
for (i in 1:nrow(term.4)) {
  if (term.4$lemma[i] == "computerscus") {
    term.4$lemma[i] <- "computersci"
  }
  if (term.4$lemma[i] == "cognitivescus") {
    term.4$lemma[i] <- "cognitivesci"
  }
}

write.csv(term.4, "nlp_result_final_dedup_deselect.csv", row.names = F)
rm(list = ls())

# Term
term <- read.csv("nlp_result_final_dedup_deselect.csv",
                 stringsAsFactors = F) # 31081
subject <- read.csv("subject.table.csv",
                    stringsAsFactors = F) # 46556
data <- read.csv("data_final.csv",
                 stringsAsFactors = F)
subject <- subject[subject$id %in% data$UT,] # 15,524
write.csv(subject, "subject.table.csv", row.names = F) # 35753
subject.table <- data.frame(table(subject$subject))
subject.table <- subject.table[order(subject.table$Freq,
                                     decreasing = T),]
colnames(subject.table) <- c("Subject", "Count")
row.names(subject.table) <- 1:nrow(subject.table)
write.csv(subject.table, "subject_table.csv", row.names = F)
rm(list = ls())

# Total frequency of words

data <- read.csv("data_final.csv",
                 stringsAsFactors = F)
# x = "a" (all), or "conj" (conjucture), or "verb"
get.term <- function(x) {
  term <- read.csv("nlp_result_final_dedup_deselect.csv",
                   stringsAsFactors = F) # 18,733
  if (x == "a") {
    term <- term
  } else if (x == "conj") {
    term <- term[term$relation.1 == "conj",]
  } else if (x == "verb") {
    term <- term[term$POS.1 == "VB",]
  } else if (x == "noun") {
    term <- term[term$POS.1 == "NN",]
  }
  
  return(term)
}

# x = "a" (all), or "conj" (conjucture), or "verb"
get.table <- function(x) {
  term <- get.term(x)
  term.table <- data.frame(table(term$lemma), stringsAsFactors = F)
  term.table <- term.table[order(term.table$Freq, decreasing = T),]
  colnames(term.table) <- c("Term", "Count")
  row.names(term.table) <- 1:nrow(term.table)
  term.table$Percentage <- term.table$Count / sum(term.table$Count)
  term.table$Term <- as.character(term.table$Term)
  
  return(term.table)
}

# x = "a" (all), or "conj" (conjucture), or "verb"
get.slope <- function(x) {
  term <- get.term(x)
  term.table <- get.table(x)

  term.table.1 <- term.table[,1:2]
  min.list <- c(1995, 2000, 2005, 2010, 2015)
  max.list <- c(1999, 2004, 2009, 2014, 2017)
  time.table <- data.frame(time = numeric(0),
                          count = numeric(0),
                          stringsAsFactors = F)

  for (i in 1:length(min.list)) {
    data.sub <- data$UT[data$PY > min.list[i]-1 & data$PY < max.list[i] +1]
    term.sub <- term[term$id %in% data.sub,]
    time.table[i, 1] <- i
    time.table[i, 2] <- nrow(term.sub)
    term.table.sub <- data.frame(table(term.sub$lemma), stringsAsFactors = F)
    term.table.sub <- term.table.sub[order(term.table.sub$Freq, decreasing = T),]
    colnames(term.table.sub) <- c("Term", "Frequency")
    rownames(term.table.sub) <- 1:nrow(term.table.sub)
    term.table.1 <- merge(x = term.table.1,
                          y = term.table.sub,
                          by.x = "Term",
                          by.y = "Term",
                          all.x = T,
                          all.y = T)
  }
  time.table[,3] <- time.table[,2] / time.table[4,2]

  for (i in 1:nrow(term.table.1)) {
    for (j in 1:ncol(term.table.1)) {
      if (is.na(term.table.1[i, j]) == T) {
        term.table.1[i, j] <- 0
      }
    }
    term.table.1[i, 2] <- sum(term.table.1[i, 3:7])
  }
  term.table.1 <- term.table.1[term.table.1[,2] >= 18,]

  for (i in 1:nrow(term.table.1)) {
    if (term.table.1[i, 2] > 0) {
      for (j in 1:nrow(time.table)) {
        term.table.1[i, j + 2] <- term.table.1[i, j + 2] / time.table[j, 3]
      }
      new.total = sum(term.table.1[i, 3:(nrow(time.table) + 2)])
      for (j in 1:nrow(time.table)) {
        term.table.1[i, j + 2] <- term.table.1[i, j + 2] / new.total
      }
      model <- lm(as.numeric(term.table.1[i, 3:(nrow(time.table) + 2)]) ~ c(1:nrow(time.table)))
      term.table.1[i, "slope"] <- model$coefficients[2]
    }
  }

  term.table.1[,"slope"] <- round(term.table.1[,"slope"], digits = 4)
  colnames(term.table.1)[2:8] <- c("total",
                                   "95-99",
                                   "00-04",
                                   "05-09",
                                   "10-14",
                                   "15-17")
  return(term.table.1)
}

term.table <- get.slope("a")
colnames(term.table)[2:8] <- c("total", "95-99", "00-04", "05-09", "10-14", "15-17", "Slope")
write.csv(term.table, "term_slope_all.csv", row.names = F)
term.table.conj <- get.slope("conj")
colnames(term.table.conj)[2:8] <- c("total", "95-99", "00-04", "05-09", "10-14", "15-17", "Slope")
write.csv(term.table.conj, "term_slope_conj.csv", row.names = F)
# term.table.noun <- get.slope("noun")
# colnames(term.table.noun)[2:8] <- c("total", "95-99", "00-04", "05-09", "10-14", "15-17", "Slope")
# write.csv(term.table.noun, "term_slope_noun.csv", row.names = F)
# rm(list = ls())

# Domain
domain.1 <- read.csv("domain_mapping.csv",
                     stringsAsFactors = F)
domain.2 <- read.csv("domain.csv",
                   stringsAsFactors = F)
domain <- merge(x = domain.1[, 3:4],
                y = domain.2,
                by.x = "field",
                by.y = "domain",
                all = T)
for (i in 1:nrow(domain)) {
  if (is.na(domain$new_domain[i]) == T) {
    if (domain$field[i] == "Arts & Humanities, general") {
      domain$new_domain[i] <- "Social science and humanities"
    } else {
      domain$new_domain[i] <- "Science"
    }
  }
}

domain <- read.csv("domain_all.csv",
                   stringsAsFactors = F)
domain$term <- tolower(domain$term)
subject <- read.csv("subject.table.csv",
                    stringsAsFactors = F)
subject$subject <- tolower(subject$subject)
subject <- merge(x = subject[,1:2],
                 y = domain[,2:3],
                 by.x = "subject",
                 by.y = "term",
                 all.x = T,
                 all.y = F)
for (i in 1:nrow(subject)) {
  if (is.na(subject$new_domain[i]) == T) {
    subject$new_domain[i] <- "Social science and humanities"
  }
}

write.csv(subject, "subject.table.csv",
          row.names = F)
rm(list = ls())

# Summarize subject
subject <- read.csv("subject.table.csv",
                    stringsAsFactors = F)
subject.1 <- subject[,2:3]
subject.1 <- subject.1[duplicated(subject.1) == 0,] # 11,308
subject.table.1 <- data.frame(table(subject.1$new_domain))
subject.table.1 <- subject.table.1[order(subject.table.1$Freq, decreasing = T),]
row.names(subject.table.1) <- 1:nrow(subject.table.1)
colnames(subject.table.1) <- c("Domain", "Count")
write.csv(subject.table.1, "subject.summary.table.csv",
          row.names = F)
rm(list = ls())

# Calculate publication slopes

wos.term <- read.csv("term_wos_publication.csv",
                     stringsAsFactors = F)
wos.pub <- read.csv("wos_publication.csv",
                    stringsAsFactors = F)
min.list <- c(1995, 2000, 2005, 2010, 2015)
max.list <- c(1999, 2004, 2009, 2014, 2017)
term.list <- as.character(unique(wos.term$Term))
term.publication.freq <- data.frame(term = character(0),
                                    frequency.1 = numeric(0),
                                    frequency.2 = numeric(0),
                                    frequency.3 = numeric(0),
                                    frequency.4 = numeric(0),
                                    frequency.5 = numeric(0),
                                    stringsAsFactors = F)

for (i in 1:8) {
  term.publication.freq[i, 1] <- term.list[i]
  for (j in 1:5) {
    year.min <- min.list[j]
    year.max <- max.list[j]
    term.publication.freq[i, j+1] <- sum(wos.term$Count[which(wos.term$Term == term.list[i] &
                                                                wos.term$Year > year.min -1 &
                                                                wos.term$Year < year.max + 1)])
  }
  term.publication.freq[i, 7] <- sum(term.publication.freq[i, 2:6])
}

year.publication.freq <- data.frame(Year = numeric(0),
                                    count = numeric(0),
                                    stringsAsFactors = F)
for (i in 1:5) {
  year.publication.freq[i, 1] <- paste("period.", i, sep = "")
  year.min <- min.list[i]
  year.max <- max.list[i]
  year.publication.freq[i, 2] <- sum(wos.pub$Publication[wos.pub$Year > year.min -1 & 
                                                           wos.pub$Year < year.max +1])
}

year.publication.freq[,3] <- year.publication.freq[,2] / year.publication.freq[3, 2]

for (i in 1:5) {
  term.publication.freq[, i + 1] <- term.publication.freq[, i + 1] / year.publication.freq[i, 3]
}

for (i in 1:8) {
  term.publication.freq[i, 7] <- sum(term.publication.freq[i, 2:6])
  for (j in 1:5) {
    term.publication.freq[i, j + 1] <- term.publication.freq[i, j + 1] / term.publication.freq[i, 7]
  }
  term.publication.freq[i, 7] <- sum(term.publication.freq[i, 2:6])
  model <- lm(as.numeric(term.publication.freq[i, 2:6]) ~ c(1:5))
  term.publication.freq[i, "slope"] <- model$coefficients[2]
}
write.csv(term.publication.freq, "term_pub_freq_all.csv", row.names = F)
rm(list = ls())

# Form a table

term.publication.freq <- read.csv("term_pub_freq_all.csv", stringsAsFactors = F)
# term.publication.freq <- term.publication.freq[! term.publication.freq$term %in% c("cybernetics", "machine learning"),]
term.table.conj <- read.csv("term_slope_conj.csv", stringsAsFactors = F)

term.frequency.visualization <- data.frame(term = character(0),
                                           group = character(0),
                                           period = numeric(0),
                                           ratio = numeric(0),
                                           stringsAsFactors = F)
term.1 <- c("Human-computer interaction", "Robot", "Natural language processing", "Artificial life", "Expert system",
            "Cognitive psychology", "cybernetics", "machine learning")
term.2 <- c("humancomputerinteract", "robot", "naturallanguageprocess", "artificiallif", "expertsystem", "cognitivepsychologi", "cybernetics", "machinelearn")

for (i in 1:8) {
  term.a <- term.1[i]
  term.b <- term.2[i]
  for (j in 1:5) {
    period = j
    ratio.1 <- term.publication.freq[which(term.publication.freq$term == term.a), j + 1]
    ratio.2 <- term.table.conj[which(term.table.conj$Term == term.b), j + 2]
    term.frequency.visualization.df.1 <- data.frame(term = term.a,
                                                    group = "P",
                                                    period = period,
                                                    ratio = ratio.1,
                                                    stringsAsFactors = F)
    term.frequency.visualization.df.2 <- data.frame(term = term.a,
                                                    group = "T",
                                                    period = period,
                                                    ratio = ratio.2,
                                                    stringsAsFactors = F)
    term.frequency.visualization <- rbind(term.frequency.visualization,
                                          term.frequency.visualization.df.1)
    term.frequency.visualization <- rbind(term.frequency.visualization,
                                          term.frequency.visualization.df.2)
  }
}

write.csv(term.frequency.visualization, "term.visualization.csv", row.names = F)


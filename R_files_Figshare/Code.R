load("dataset.RData")
library(extrafont)
library(scales)
library(MASS)
library(texreg)
library(ggplot2)
library(tidyverse)
library(tikzDevice)
library(extrafont)
library(kableExtra)
library(Hmisc)
options(scipen=999)

### Figure 1
dataset %>%
  ggplot(aes(hour, fill = time_cat)) +
  geom_histogram(bins = 50) + scale_x_datetime(breaks = "6 hours",
                                               labels = date_format("%d-%m\\n%H hrs.")) +
  xlab("") + ylab("N of Tweets") + scale_fill_manual(name = "Time",
                                                     labels = c("Before/After the event",
                                                                "During the event"),
                                                     values=c("#C0C0C0", "#808080")) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE)) +
  theme(text = element_text(family="Calibri", face="bold", size=16),
        legend.position="bottom")

## Figure 2

dataset %>%
  ggplot(aes(rt_count)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE)) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE)) +
  theme(text = element_text(family="Calibri", face="bold", size=16)) +
  xlab("N of retweets") + ylab("N of tweets")

## Figure 3

dataset %>% dplyr::select(EmoWord, MoralWord, EmoMoralWord) %>%
  gather("word_type", "word_count") %>% 
  mutate_at(vars("word_type"), as.character) %>%
  mutate(word_type = ifelse(str_detect(word_type, "EmoMoralWord"),
                            "Moral-emotional words",
                            ifelse(str_detect(word_type, "MoralWord"),
                                   "Moral words", "Emotional words"))) %>%
  ggplot(aes(word_count)) +
  geom_histogram(bins = 4) + facet_grid(.~word_type) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE)) +
  theme(text = element_text(family="Calibri", face="bold", size=16)) +
  xlab("Q of words") + ylab("N of tweets")

## Models

m1 <- glm.nb(rt_count ~ time_cat*OnlyMoralWord +
               time_cat*OnlyEmoWord + time_cat*EmoMoralWord +
               url + media +
               foll_div10,
             dataset, control = glm.control(maxit = 200))

m1.pos <- glm(rt_count ~ time_cat*OnlyMoralWord +
                time_cat*OnlyEmoWord + time_cat*EmoMoralWord +
                url + media +
                foll_div10,
              dataset, family = "poisson")

## Table 6: Likelihood ratio test

chi <- 1 - pchisq(-2*(logLik(m1.pos)[1]-logLik(m1)[1]), df = 1)

data.frame(row.names = c("Poisson", "Negative Binomial"),
           LogLik = c(2*logLik(m1.pos)[1], 2*logLik(m1)[1]),
           df = c("", 1),
           Pr = c("", chi))

## Table 7: Correlation Matrix

dataset %>%
  dplyr::select(rt_count, OnlyMoralWord,
                OnlyEmoWord, EmoMoralWord,
                url, media, foll_div10) -> sel_vars

rcorr(as.matrix(sel_vars))[[1]]

## Effects

### During the event
100*(exp(m1$coefficients[[2]])-1)
### Moral Words
100*(exp(m1$coefficients[[3]])-1)
### Emotional Words
100*(exp(m1$coefficients[[4]])-1)
### Moral-emotional words
100*(exp(m1$coefficients[[5]])-1)
### Moral words during the event
100*(exp(m1$coefficients[[3]]+m1$coefficients[[9]])-1)
### Emotional words during the event
100*(exp(m1$coefficients[[4]]+m1$coefficients[[10]])-1)
### Moral-emotional words during the event
100*(exp(m1$coefficients[[5]]+m1$coefficients[[11]])-1)

## Figure 4

data.frame(time_cat = rep(factor(c("no_cuenta", "cuenta")),
                           each = 12),
           OnlyMoralWord = rep(c(rep(0, 8), 0:3),
                                     times = 2),
           OnlyEmoWord = rep(c(rep(0, 4), 0:3, rep(0, 4)),
                                   times = 2),
           EmoMoralWord = rep(c(0:3, rep(0, 8))),
           url = FALSE,
           media = FALSE,
           foll_div10 = mean(dataset$foll_div10)) -> pred_data

cbind(pred_data,
      predict(m1,
              pred_data,
              type = "link",
              se.fit = TRUE)) %>%
  mutate(LL = fit - 1.96*se.fit,
         UL = fit + 1.96*se.fit) %>%
  mutate(word_cat =  rep(c("Moral-emotional words", "Emotional words", "Moral words"),
                         each = 4, times = 2),
         word_freq = rep(0:3, times = 6)) %>%
  mutate_at(vars(c("fit", "LL", "UL")),
            .funs = list(exp = ~exp(.))) %>%
  mutate(time_cat = ifelse(time_cat == "cuenta",
                            "During the event",
                            "Before/after the event")) %>%
  ggplot(aes(word_freq, fit_exp)) +
  geom_line(aes(colour = time_cat), size = 1) +
  facet_grid(~word_cat, scales="free_y") +
  geom_ribbon(aes(ymin = LL_exp, ymax = UL_exp, fill = time_cat), alpha = .25) +
  labs(x = "N of words", y = "Predicted retweets", color = "", fill = "") +
  scale_fill_manual(values=c("#A9A9A9", "#696969")) +
  scale_color_manual(values=c("#A9A9A9", "#696969")) +
  theme(text = element_text(family="Calibri", face="bold", size=16),
        legend.position="bottom")
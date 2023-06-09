## Install packages ###################################
install.packages("install.load")
library("install.load")
install_load("readr",
             "MASS",
             "dplyr",
             "forcats",
             "magrittr",
             "ggplot2",
             "cowplot",
             "viridis",
             "tidyr",
             "lme4")


## ggplot theme  ###################################
theme_tom <- function(base_size = 8, base_family = "") {
          # Starts with theme_grey and then modify some parts
          theme_grey(base_size = base_size, base_family = base_family) %+replace%
                    theme(
                              # theme_bw() stuff:
                              legend.key        = element_rect(colour = "grey80"),
                              panel.background  = element_rect(fill = "white", colour = NA),
                              panel.border      = element_rect(fill = NA, colour = "grey50"),
                              panel.grid.major  = element_line(colour = "grey90", size = 0.2),
                              panel.grid.minor  = element_line(colour = "grey98", size = 0.5),
                              strip.background  = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
                              # Custom stuff
                              axis.text         = element_text(size = 8),
                              axis.ticks        = element_line(colour = "black", size = 0.2),
                              plot.title        = element_text(size = 12, hjust = 0, face = "bold", margin=margin(b = 3, unit = "pt")),
                              legend.text       = element_text(size = 7),
                              legend.key.size   = unit(3, "mm"),
                              strip.text        = element_text(size = 8),
                              legend.background = element_rect(colour = 1, fill = scales::alpha("white", 0.9), size = 0.25)
                    )
}


## Load & tidy data ###################################
# Load data - 10 ornithology journals
data <- read_csv("data_ornithology.csv") %>% 
          mutate(journal = fct_reorder(journal, jif))

# Load data - 2014 data augmented with generalist journals
data2 <- read_csv("data_cites.csv") 

# Get journal impact factors
jifs <- unique(select(data, journal, jif))


# Summarise AAS by year and source
year_summary <- data %>% group_by(year) %>% summarise(score = mean(score))
src_summary  <- data %>% group_by(journal) %>% summarise(score = mean(score), jif = mean(jif))

# Multiply mentions by weights, sum to get 'estimated' score
data_weighted <- data %>% select(journal, jif, year, score, news, bloggers, wikipedia, tweeters, f1000, google, facebook, reddit) %>% 
          mutate(news = news * 8,
                 bloggers = bloggers * 5,
                 wikipedia = wikipedia * 3,
                 facebook = facebook * 0.25,
                 reddit =  reddit * 0.25) %>% 
          mutate(score_estimate = news + bloggers + wikipedia + tweeters + f1000 + google + facebook + reddit)
cor.test(~score + score_estimate, data = data_weighted) # Convincing correlation between the two

# Overall contribution per platform
contr_tot <- data_weighted %>% select(-journal, -jif, -score, -year) %>% summarise_each(funs = 'mean') %>% 
          gather(platform, score, -score_estimate) %>% 
          mutate(score = score/score_estimate) %>% select(-score_estimate) %>%  
          mutate(platform = fct_reorder(platform, score)) %>%
          filter(! platform %in% c("reddit", "f1000")) %>% arrange(-score)


# Contribution by year/source
contr_yr <- data_weighted %>% select(-journal, -jif, -score, -score_estimate) %>% 
          gather(platform, score, -year) %>% group_by(year, platform) %>% 
          summarise(score = mean(score)) %>% 
          mutate(platform = fct_relevel(platform, "tweeters", "news", "bloggers")) %>% 
          mutate(platform = fct_recode(platform, Twitter = "tweeters", News = "news", Blogging = "bloggers", Facebook = "facebook"))
          
contr_src <- data_weighted %>% select(-year, -jif, -score, -score_estimate) %>% 
          gather(platform, score, -journal) %>% group_by(journal, platform) %>% 
          summarise(score = mean(score)) %>% 
          mutate(platform = fct_relevel(platform, "tweeters", "news", "bloggers")) %>% 
          mutate(platform = fct_recode(platform, Twitter = "tweeters", News = "news", Blogging = "bloggers", Facebook = "facebook"))


# Summary stats
nrow(data)  # 2667 articles
nrow(data2) # 878  articles
table(data2$subject)
data  %>% group_by(journal) %>% summarise(n = length(score)) %>% arrange(n) 
data2 %>% group_by(journal) %>% summarise(n = length(score)) %>% arrange(-n) 

## SECTION 1: Patterns in Altmetric Attention Scores ###################################
## Simple linear models
lmod1 <- lm(log(score) ~ year, data = data)
lmod2 <- lm(log(score) ~ journal, data = data)
lmod3 <- lm(log(score) ~ jif, data = data)
summary(lmod1)
summary(lmod2)
summary(lmod3)

## Summary stats
data %>% group_by(year)    %>% summarise(mean = mean(score)) %>% arrange(mean) # 2.8 to 18.3 (2012 to 2016)
data %>% group_by(journal) %>% summarise(mean = mean(score)) %>% arrange(mean) # 4.0 to 14.7 (Emu to Auk)

## Plot
fig1a <- ggplot(data, aes(x = year, y = score, group = year)) +
          geom_point(aes(col = log(score)), position = position_jitter(0.12, 0), alpha = 0.5) +
          geom_point(data = year_summary, fill = "white", col = 1, shape = 21, size = 2) +
          geom_path(data = year_summary, aes(group = 1), linetype = 2) +
          theme_tom() + labs(x = "Year of publication", y = "Altmetric Attention Score") +
          scale_color_viridis(direction = -1, option = "D") +
          guides(fill = F, col = F) +
          scale_y_log10()

fig1b <- ggplot(data, aes(x = journal, y = score, group = journal)) +
          geom_point(aes(col = log(score)), position = position_jitter(0.25, 0), alpha = 0.5) +
          geom_point(data = src_summary, fill = "white", col = 1, shape = 21, size = 2) +
          geom_path(data = src_summary, aes(group = 1), linetype = 2) +
          theme_tom() + labs(x = NULL, y = "Altmetric Attention Score") +
          scale_color_viridis(direction = -1, option = "D") +
          guides(fill = F, col = F) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_y_log10()

plot_grid(fig1a, fig1b, nrow = 2, labels = c("a", "b"), rel_heights = c(1, 1.3))
# ggsave("fig1.png", width = 5.5, height = 6, dpi = 600)


## SECTION 2: Which online platforms drive variation in Altmetric Attention Score? ###################################
## Summary stats
sum(contr_tot$score) # >99.99% by top 6 platforms

## Neg bin regression model
nbmod1 <- glm.nb(tweeters ~ year + jif, data = data_weighted)
nbmod2 <- glm.nb(facebook ~ year + jif, data = data_weighted)
nbmod3 <- glm.nb(news ~ year + jif, data = data_weighted)
nbmod4 <- glm.nb(bloggers ~ year + jif, data = data_weighted)

confint(nbmod1)
confint(nbmod2)
confint(nbmod3)
confint(nbmod4)

## Plot overall contribution
ggplot(contr_tot, aes(x = NA, y = score, fill = platform)) + 
          geom_col(col = "white") + 
          coord_flip() + scale_y_continuous(labels = scales:::percent, breaks = seq(0, 1, by = 0.1)) +
          scale_fill_viridis(direction = -1, option = "D", discrete = T, 
                             labels = c("Google+", "Wikipedia", "Facebook", "Blogging", "News", "Twitter")) +
          labs(y = "Total contribution to Altmetric Attention Score", x = NULL, fill = NULL) +
          theme_tom() + theme(legend.position = "bottom", legend.direction = "horizontal",
                              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                              panel.grid = element_blank()) +
          guides(fill = guide_legend(nrow = 1, reverse = T))
# ggsave("fig2.png", width = 5.5, height = 2, dpi = 600)

## Plot contribution by year/journal
fig3a <- ggplot(filter(contr_yr, platform %in% c("Twitter", "Blogging", "News", "Facebook")), aes(x = year, y = score, group = platform)) + 
          geom_line(aes(colour = platform)) +  
          geom_point(shape = 21, col = 1, aes(fill = platform), size = 2) + 
          theme_tom() + labs(x = "Year of publication", y = "Average weighted mentions") +
          scale_color_viridis(direction = 1, option = "D", discrete = T, name = NULL, begin = 0, end = 0.6) +
          scale_fill_viridis(direction = 1, option = "D", discrete = T, name = NULL, begin = 0, end = 0.6) +
          theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
          guides(col = F, fill = F) + facet_wrap(~platform, nrow = 1, scales = "free_y")

fig3b <- ggplot(filter(contr_src, platform %in% c("Twitter", "Blogging", "News", "Facebook")), aes(x = journal, y = score, group = platform)) + 
          geom_line(aes(colour = platform)) +  
          geom_point(shape = 21, col = 1, aes(fill = platform), size = 2) + 
          theme_tom() + labs(x = NULL, y = "Average weighted mentions") +
          scale_color_viridis(direction = 1, option = "D", discrete = T, begin = 0, end = 0.6) +
          scale_fill_viridis(direction = 1, option = "D", discrete = T, begin = 0, end = 0.6) +
          theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.5)) +
          guides(col = F, fill = F) + facet_wrap(~platform, nrow = 1, scales = "free_y")

plot_grid(fig3a, fig3b, nrow = 2, labels = c("a", "b"), rel_heights = c(1, 1.2))
# ggsave("fig3.png", width = 5.5, height = 5, dpi = 600)


## SECTION 3: Is there an association between Altmetric Attention Score and citation count? ###################################
# Log transform variables
data2 <- data2 %>% mutate(logcit = log(citations + 0.1),
                          logscore = log(score + 0.1),
                          logjif = log(jif))

## Linear mixed model
# Fit model
lmmod1 <- lmer(logcit ~ scale(logscore) * scale(logjif) + (1|journal), 
               data = data2, REML = F, na.action = 'na.fail')
confint(lmmod1)

# New data
newx1 <- expand.grid('logscore' = seq(log(min(data2$score)+0.1), log(max(data2$score)+0.1), length.out = 50),
                     'logjif' = log(c(min(data2$jif), median(data2$jif), max(data2$jif))))

# Get fitted values with bootstrapped sd
predfun1 <- function(x){predict(x, newdata = newx1, re.form = NA)}
booted1 <- bootMer(lmmod1, predfun1, nsim = 100)
fit1 <- data.frame(newx1, citations = booted1$t0, 
                   low = apply(booted1$t, 2, quantile, probs = 0.025),
                   upp = apply(booted1$t, 2, quantile, probs = 0.975)) %>% 
          as.tbl %>% mutate(jif = factor(exp(logjif)))

# Plot predictions
fig4a <- ggplot(data = fit1) +
          geom_line(aes(x = exp(logscore), y = exp(citations), group = jif, col = jif)) +
          geom_ribbon(aes(x = exp(logscore), ymin = exp(low), ymax = exp(upp), group = jif, fill = jif), alpha = 0.3) +
          labs(x = "Altmetric Attention Score", y = "Citation count", 
               fill = "Journal Impact Factor", col = "Journal Impact Factor") +
          theme_tom() + guides(fill = F, col = F) +
          # scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000)) +
          # scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
          scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T)
fig4a

## Generalised linear mixed model (binomial)
# Fit model
glmmod1 <- glmer(cited ~ scale(logscore) * scale(logjif) + (1|journal), 
                 data = data2, na.action = 'na.fail', family = 'binomial')
confint(glmmod1)

# New data
newx2 <- expand.grid('logscore' = seq(log(min(data2$score)+0.1), log(max(data2$score)+0.1), length.out = 50),
                     'logjif' = log(c(min(data2$jif), median(data2$jif), max(data2$jif))))

# Get fitted values with bootstrapped sd
predfun2 <- function(x){predict(x, newdata = newx2, type = 'response', re.form = NA)}
booted2 <- bootMer(glmmod1, predfun2, nsim = 100)
fit2 <- data.frame(newx2, citations = booted2$t0,
                   low = apply(booted2$t, 2, quantile, probs = 0.025),
                   upp = apply(booted2$t, 2, quantile, probs = 0.975)) %>% 
          as.tbl %>% mutate(jif = factor(exp(logjif)))

# Plot predictions
fig4b <- ggplot(data = fit2) +
          geom_line(aes(x = exp(logscore), y = (citations), group = jif, col = jif)) +
          geom_ribbon(aes(x = exp(logscore), ymin = (low), ymax = (upp), group = jif, fill = jif), alpha = 0.3) +
          labs(x = "Altmetric Attention Score", y = "Probability of being cited", 
               fill = "Journal Impact Factor", col = "Journal Impact Factor") +
          theme_tom() + gglegend.lowright() +
          scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000)) +
          scale_fill_viridis(discrete = T) + scale_color_viridis(discrete = T)

## Merge plots and save
plot_grid(fig4a, fig4b, nrow = 2, labels = c("a", "b"))
# ggsave("fig4.png", width = 5.5, height = 6, dpi = 600)


## MISC ##########################
# Correlation between twitter followers and twitter mentions among journals
contr_src %>% filter(platform == "Twitter") %>% 
          mutate(followers = c(0, 4445, 210, 10800, 6670, 2360, 0, 255, 4910, 3660)) %>% 
          # filter(journal != "Bird Study") %>% 
          cor.test(formula = ~followers+absolute, data = .)



###################################################################
#
# ANALYSIS A: INJURIES
#
###################################################################
###################################################################

# Load data:
dat.injuries <- read.csv2(file = "AWD_injuries.csv", stringsAsFactors = FALSE)

# data structure:
dat.injuries$sex <- as.factor(dat.injuries$sex)
dat.injuries$PFQ <- as.factor(dat.injuries$PFQ)
str(dat.injuries)

# data summary:
summary(dat.injuries)

# how many individuals in each PFQ?
aggregate(idDog ~ PFQ, data = dat.injuries, FUN = function(x) length(unique(x)))

# Test if individuals with higher PFQ are more likely to get injuries compared to the reference category PFQ=2
# --> We run a binomial GLMM with random effect individual to account for repeated observations of the same dog
# --> We add sex as fixed effect to control for potential confounding effects of sex
library(lme4)
M.injuries <- glmer(formula = injured ~ PFQ + sex + (1|idDog), data = dat.injuries, family = "binomial")
summary(M.injuries)

# PLOT outpot of covariates affecting injuries:
library(ggplot2)
library(cowplot)

df.injured <- as.data.frame(summary(M.injuries)$coefficients)
str(df.injured)

# add bounds of 95% confidence interval:
df.injured$CI.lower <- df.injured$Estimate-1.96*df.injured$`Std. Error`
df.injured$CI.upper <- df.injured$Estimate+1.96*df.injured$`Std. Error`

# transform coefficient estimate into odds-ratio:
df.injured$odds <- exp(df.injured$Estimate)
df.injured$odds.lower <- exp(df.injured$CI.lower)
df.injured$odds.upper <- exp(df.injured$CI.upper)

# add alternative representation of p-value:
df.injured$p.value <- round(df.injured$`Pr(>|z|)`,2)
df.injured$p.label <- ifelse(test = df.injured$p.value>0, yes = paste0("p = ",sprintf("%0.2f",df.injured$p.value)), no = "p < 0.01")

# format names of covariates:
df.injured$term <- rownames(df.injured)
df.injured$term <- gsub(pattern = "\\\\(", replacement = "", x = df.injured$term)
df.injured$term <- gsub(pattern = "\\\\)", replacement = "", x = df.injured$term)
df.injured$term<- gsub(pattern = "PFQ", replacement = "PFQ = ", x = df.injured$term)
df.injured$term <- gsub(pattern = "\\\\*", replacement = "", x = df.injured$term)
df.injured$term <- gsub(pattern = "sexM", replacement = "Sex (M)", x = df.injured$term)

df.injured

# we remove intercept and define the order in which covariates appear in the plot
df.injured <- subset(df.injured, term!="Intercept")
df.injured$term.sorted <- factor(x = df.injured$term, levels = c("Sex (M)", "PFQ = 7", "PFQ = 6", "PFQ = 5", "PFQ = 4","PFQ = 3"), ordered = TRUE)

theme_set(theme_minimal())
fontSize <- 20
font <- "Helvetica"

p.injured <- ggplot(data = df.injured, aes(x = term.sorted, label = p.label))+
  geom_pointrange(
    aes(y = odds, ymin = odds.lower, ymax=odds.upper), 
    size = 1.1, col="black")+
  xlab("")+
  scale_y_log10(breaks = c(0,0.1,0.5,1,2,5), labels = c(0,0.1,0.5,1,2,5), limits = c(0.1, 20), name="Odds ratio      ")+
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.4)+
  geom_text(y=1.2, size = 5)+
  coord_flip()+
  theme(axis.text = element_text(colour = "black", size = fontSize, family = font),
        title = element_text(colour = "black", face = "bold", size = fontSize, family = font, hjust = 0, vjust = 0),
        axis.title.x = element_text(colour = "black", size = fontSize, family = font, hjust = 0.5, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = "black", face = "bold", size = fontSize, family = font),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = fontSize, family = font),
        legend.title = element_text(colour = "black", size = fontSize, family = font),
        legend.justification = 0.25,
        plot.margin = margin(5.5,100,5.5,5.5)
  ) 

print(p.injured)



###################################################################
###################################################################
#
# ANALYSIS B: BLOOD
#
###################################################################
###################################################################

# Load data:
dat.blood <- read.csv2(file = "AWD_blood.csv", stringsAsFactors = FALSE)

# data structure:
dat.blood$sex <- as.factor(dat.blood$sex)
dat.blood$PFQ <- as.factor(dat.blood$PFQ)
str(dat.blood)

# data summary:
summary(dat.blood)

# how many individuals in each PFQ?
aggregate(idDog ~ PFQ, data = dat.blood, FUN = function(x) length(unique(x)))


# Test if individuals with higher PFQ are more likely to get blood marks on legs compared to the reference category PFQ=2
# --> We run a binomial GLMM with random effect individual to account for repeated observations of the same dog
# --> We add sex as fixed effect to control for potential confounding effects of sex
M.blood <- glmer(formula = Blood_Legs ~ PFQ + sex + (1|idDog), data = dat.blood, family = "binomial")
summary(M.blood)


# PLOT outpot of covariates affecting blood marks on legs:
df.blood <- as.data.frame(summary(M.blood)$coefficients)
str(df.blood)

# add bounds of 95% confidence interval:
df.blood$CI.lower <- df.blood$Estimate-1.96*df.blood$`Std. Error`
df.blood$CI.upper <- df.blood$Estimate+1.96*df.blood$`Std. Error`

# transform coefficient estimate into odds-ratio:
df.blood$odds <- exp(df.blood$Estimate)
df.blood$odds.lower <- exp(df.blood$CI.lower)
df.blood$odds.upper <- exp(df.blood$CI.upper)

# add alternative representation of p-value:
df.blood$p.value <- round(df.blood$`Pr(>|z|)`,2)
df.blood$p.label <- ifelse(test = df.blood$p.value>0, yes = paste0("p = ",sprintf("%0.2f",df.blood$p.value)), no = "p < 0.01")

# format names of covariates:
df.blood$term <- rownames(df.blood)
df.blood$term <- gsub(pattern = "\\\\(", replacement = "", x = df.blood$term)
df.blood$term <- gsub(pattern = "\\\\)", replacement = "", x = df.blood$term)
df.blood$term<- gsub(pattern = "PFQ", replacement = "PFQ = ", x = df.blood$term)
df.blood$term <- gsub(pattern = "\\\\*", replacement = "", x = df.blood$term)
df.blood$term <- gsub(pattern = "sexM", replacement = "Sex (M)", x = df.blood$term)

df.blood

# we remove intercept and define the order in which covariates appear in the plot
df.blood <- subset(df.blood, term!="Intercept")
df.blood$term.sorted <- factor(x = df.blood$term, levels = c("Sex (M)", "PFQ = 7", "PFQ = 6", "PFQ = 5", "PFQ = 4","PFQ = 3"), ordered = TRUE)

theme_set(theme_minimal())
fontSize <- 20
font <- "Helvetica"

p.blood <- ggplot(data = df.blood, aes(x = term.sorted, label = p.label))+
  geom_pointrange(
    aes(y = odds, ymin = odds.lower, ymax=odds.upper), 
    size = 1.1, col="black")+
  xlab("")+
  scale_y_log10(breaks = c(0,0.1,0.5,1,2,5), labels = c(0,0.1,0.5,1,2,5), limits = c(NA, 20), name="Odds ratio")+
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.4)+
  geom_text(y=1.2, size = 5)+
  coord_flip()+
  theme(axis.text = element_text(colour = "black", size = fontSize, family = font),
        title = element_text(colour = "black", face = "bold", size = fontSize, family = font, hjust = 0, vjust = 0),
        axis.title.x = element_text(colour = "black", size = fontSize, family = font, hjust = 0.5, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = "black", face = "bold", size = fontSize, family = font),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = fontSize, family = font),
        legend.title = element_text(colour = "black", size = fontSize, family = font),
        legend.justification = 0.25,
        plot.margin = margin(5.5,100,5.5,5.5)
  ) 

print(p.blood)


###################################################################
# Combine both plots:

p.combinded <- plot_grid(p.injured, NULL, p.blood, nrow = 1, ncol = 3, labels = c("(a)","","(b)"), label_size = fontSize, 
                         label_x = c(-0.01,0,-0.05), label_y = rep(1.0,3), rel_widths = c(1,0.07,1), rel_heights = rep(1,3))

p.combinded

# Export as PDF:
ggsave(filename = "Fig_injuries_blood.pdf", plot = p.combinded, width = 16, height = 8, family = "Helvetica", 
       device = cairo_pdf)


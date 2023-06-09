## Analyses for paper "Children's Reliance on Non-verbal Cues of a Robot Versus a Human", by Josje Verhagen*, Rianne van den Berghe*, Ora Oudgenoeg, Aylin Kuntay & Paul Leseman (*co-first author)

## Study 1

## RQ 1 (Effects of Speaker (Agent) and Label familiarity (Label) on pointing)

data <- read.table ("XXX[File location]\\\\Data_Robot Human Nonverbal Cues_PLoS ONE_2019.txt", sep = "\\t", header = TRUE)
data2 <- subset(data,data$NrItemsSelected == 1 & SuccessfulTrial == 1 & NonverbalCue == "Point" & Perception >= 0)
data2$subject <- as.factor(data2$Subject)
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Human+Robot")
contrasts (data2$Agent) <- contrast
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Familiar+Novel")
contrasts (data2$Label) <- contrast

data2.model = glmer(PointGazeFollowing ~ (Agent|Subject) + Agent*Label, data = data2, family = "binomial", glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(data2.model)

## RQ2

data <- read.table ("XXX[File location]\\\\Data_Robot Human Nonverbal Cues_PLoS ONE_2019.txt", sep = "\\t", header = TRUE)
data2 <- subset(data,data$NrItemsSelected == 1 & SuccessfulTrial == 1 & NonverbalCue == "Point" & Perception >= 0)

Perception <- aggregate(Perception ~ Subject, data2, mean)
data2$Perception <- data2$Perception - mean(Perception$Perception)

data2$subject <- as.factor(data2$Subject)
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Human+Robot")
contrasts (data2$Agent) <- contrast
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Familiar+Novel")
contrasts (data2$Label) <- contrast

data2.model = glmer(PointGazeFollowing ~ (Agent|Subject) + Agent*Label*Perception, data = data2, family = "binomial", glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(data2.model)


## Study 2

## RQ 1 (Effects of Speaker (Agent) and Label familiarity (Label) on gaze)

data <- read.table ("XXX[File location]\\\\Data_Robot Human Nonverbal Cues_PLoS ONE_2019.txt", sep = "\\t", header = TRUE)
data2 <- subset(data,data$NrItemsSelected == 1 & SuccessfulTrial == 1 & NonverbalCue == "Gaze" & Perception >= 0)
data2$subject <- as.factor(data2$Subject)
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Human+Robot")
contrasts (data2$Agent) <- contrast
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Familiar+Novel")
contrasts (data2$Label) <- contrast

data2.model = glmer(PointGazeFollowing ~ (Agent|Subject) + Agent*Label, data = data2, family = "binomial", glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(data2.model)

## RQ2

data <- read.table ("XXX[File location]\\\\Data_Robot Human Nonverbal Cues_PLoS ONE_2019.txt", sep = "\\t", header = TRUE)
data2 <- subset(data,data$NrItemsSelected == 1 & SuccessfulTrial == 1 & NonverbalCue == "Gaze" & Perception >= 0)

Perception <- aggregate(Perception ~ Subject, data2, mean)
data2$Perception <- data2$Perception - mean(Perception$Perception)

data2$subject <- as.factor(data2$Subject)
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Human+Robot")
contrasts (data2$Agent) <- contrast
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Familiar+Novel")
contrasts (data2$Label) <- contrast

data2.model = glmer(PointGazeFollowing ~ (Agent|Subject) + Agent*Label*Perception, data = data2, family = "binomial", glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(data2.model)


# Comparison non-verbal cue following between studies (Pointing vs. Eye gaze)

data <- read.table ("XXX[File location]\\\\Data_Robot Human Nonverbal Cues_PLoS ONE_2019.txt", sep = "\\t", header = TRUE)
data2 <- subset(data,data$NrItemsSelected == 1 & SuccessfulTrial == 1 & Perception >= 0)
data2$subject <- as.factor(data2$Subject)
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Human+Robot")
contrasts (data2$Agent) <- contrast
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Familiar+Novel")
contrasts (data2$Label) <- contrast
contrast <- cbind(c(-1/2,+1/2))
colnames (contrast) <- c("-Gaze+Point")
contrasts (data2$NonverbalCue) <- contrast

data2.model = glmer(PointGazeFollowing ~ (Agent|Subject) + Agent*Label*NonverbalCue, data = data2, family = "binomial", glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(data2.model)

## Interaction plot

data2$Perception <- factor(data2$Perception, levels=c("1","2") )
interaction.plot(data2$Perception, data2$Label, data2$PointGazeFollowing, xlab="Perception", ylab="Point Following", trace.label="Label")





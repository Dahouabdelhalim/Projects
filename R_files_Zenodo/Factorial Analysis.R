# Screening Study Factorial Analysis
# Tyler Hudson
# 11/15/17

studydata <- read.csv(file.choose())
Diameter <- factor(studydata$Diameter)
Heat.intensity <- factor(studydata$Heat.Intensity)
Crossflow <- factor(studydata$Crossflow)
Species <- factor(studydata$Species)
Moisture.content <- factor(studydata$Moisture.content)
Condition <- factor(studydata$Condition)

# use this one for natural vs dowel analysis
anova <- aov(Average.Frames ~ Diameter * Heat.intensity * Crossflow * Condition * Moisture.content, data = studydata)

# use this one for species analysis
anova <- aov(Average.Frames ~ Diameter * Heat.intensity * Crossflow * Moisture.content * Species, data = studydata)

summarytable <- summary(anova)
model.tables(anova,"means")
results <- TukeyHSD(anova,ordered=TRUE)


#Save data to txt file, make sure you rename
capture.output(summarytable,file="3reps.txt")

#Save data to txt file, make sure you rename
capture.output(summarytable,file="NvsD.txt")

# species only results


interaction.plot(Diameter,Moisture.content,studydata$Average.Frames, ylab = "Average Frames")
#gp = ggplot(data=studydata, aes(x=Species,y=Average.Frames,colour=Diameter, group=Diameter))
#gp + geom_line(size=.6) + 
     geom_point(size=3) 
  
png("file.png",1000,1000)   
plotmeans(studydata$Average.Frames[Diameter==6]~Species[Diameter==6], ylim=c(100,1500),ylab="Average Frames")
plotmeans(studydata$Average.Frames[Diameter==2]~Species[Diameter==2],add=TRUE)
dev.off()

df<-with(studydata, aggregate(studydata$Average.Frame, list(Species=Species, Diameter=Diameter), mean))
df$se<-with(studydata, aggregate(studydata$Average.Frame, list(Species=Species, Diameter=Diameter), function(x) sd(x)/sqrt(10)))[,3]
gp <- ggplot(df, aes(x=Species, y=x, colour=Diameter, group=Diameter))
gp + geom_line(aes(linetype=factor(Diameter))) + 
     geom_point(aes(shape=factor(Diameter))) + 
     geom_errorbar(aes(ymax=x+se, ymin=x-se))

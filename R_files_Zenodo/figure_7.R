rm(list = ls())

setwd(".../Data_and_Code")
df <- read.csv("plant_age_height_circumference.csv")
head(df)
dim(df)
df



pdf("Figure_7_age-height-circumference-combined.pdf", width=7.5, height=4)

par(mfrow = c(1, 2))
# par(oma = c(0, 0, 0, 0)) # make room  for the overall x and y axis titles
par(mar = c(0, 1, 0, 1)) # make the plots be closer together
par(ps = 8.5, cex = 1, cex.main = 1)


    head(df)

    par(mar=c(4.2,4.2,0.2,0.2))
    sub1 <- subset(df, Site == "1")
    with (newdata <- subset(df, Site == "1"), plot(Height, Age, xlim = c(110,375), ylim=c(20,120), pch=19, col="gray80", xlab="Height (cm)", ylab="Age (yrs)", axes=F))    # remove axes=F to print all four borders
    axis(side=1)
    axis(side=2)  # to add sides to all axes, below is box() function
    reg1 = lm(Age ~ Height, data = sub1)
    clipplot(abline(reg1, lwd=2, col="gray80", ), xlim=c(min(sub1$Height), max(sub1$Height))) 
    summary(reg1)
    anova(reg1)
    
    sub2 <- subset(df, Site== "2")
    points(sub2$Height, sub2$Age, pch=19, col="blue")
    reg2 = lm(Age ~ Height, data = sub2)
    clipplot(abline(reg2, lwd=2, col="blue", ), xlim=c(min(sub2$Height), max(sub2$Height))) 
    summary(reg2)
    anova(reg2)
    
    sub3 <- subset(df, Site== "3")
    points(sub3$Height, sub3$Age, pch=19, col="orange")
    reg3 = lm(Age ~ Height, data = sub3)
    clipplot(abline(reg3, lwd=2, col="orange", ), xlim=c(min(sub3$Height), max(sub3$Height))) 
    summary(reg3)
    anova(reg3)
    
    legend_texts = expression("Site 1, p = 0.019, R"^2 == 0.57, "Site 2a, p = 0.0025, R"^2 == 0.71, "Site 2b, p = 0.002, R"^2 == 0.54)
    legend(x=105, y=120, col=c("gray80","blue", "orange"), lty=1,lwd=2, cex=0.9, legend=legend_texts, bty='n')
    
    
    par(xpd=T, mar=par()$mar+c(0,10,0,0))
    # Restore default clipping rect
    
    text(x=42, y=120, "(a)", font=2)
    
    

    
    
    par(mar=c(4.2,4.2,0.2,0.2))
    sub1 <- subset(df, Site == "1")
    with (newdata <- subset(df, Site == "1"), plot(Height, Perimeter_cm, xlim = c(110,375), ylim=c(10,30), pch=19, col="gray80", 
                                                        xlab="Height (cm)", ylab="Circumference (cm)", axes=F))	# remove axes=F to print all four borders
    axis(side=1)
    axis(side=2)  # to add sides to all axes, below is box() function
    reg1 = lm(Perimeter_cm ~ Height, data = sub1)
    clipplot(abline(reg1, lwd=2, col="gray80", ), xlim=c(min(sub1$Height), max(sub1$Height))) 
    summary(reg1)
    anova(reg1)
    
    sub2 <- subset(df, Site== "2")
    points(sub2$Height, sub2$Perimeter_cm, pch=19, col="blue")
    reg2 = lm(Perimeter_cm ~ Height, data = sub2)
    clipplot(abline(reg2, lwd=2, col="blue", ), xlim=c(min(sub2$Height), max(sub2$Height))) 
    summary(reg2)
    anova(reg2)
    
    sub3 <- subset(df, Site== "3")
    sub3a <- subset(sub3, Perimeter_cm < 34)
    rm(sub3)
    points(sub3a$Height, sub3a$Perimeter_cm, pch=19, col="orange")
    reg3 = lm(Perimeter_cm ~ Height, data = sub3a)
    clipplot(abline(reg3, lwd=2, col="orange", ), xlim=c(min(sub3a$Height), max(sub3a$Height))) 
    summary(reg3)
    anova(reg3)
    
    legend_texts = expression("Site 1, p = 0.14, R"^2 == 0.28, "Site 2a, p = 0.10, R"^2 == 0.34, "Site 2b, p = 0.057, R"^2 == 0.27)
    legend(x=105, y=30, col=c("gray80","blue", "orange"), lty=1,lwd=2, cex=0.9, legend=legend_texts, bty='n')
    
    par(xpd=T, mar=par()$mar+c(0,6,0,0))
    text(x=42, y=30, "(b)", font=2)


dev.off()


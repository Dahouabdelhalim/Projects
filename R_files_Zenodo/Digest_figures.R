### Last modified 10/01/2019
### zachmiller@uchicago.edu

### Code to generate simulated data and figures for:
### ZR Miller (2019), Evolution, Digest: Does sexual conflict complicate a trade-off between fecundity and survival?
### https://doi.org/10.1111/evo.13855

library(tidyverse)
library(gridExtra)

set.seed(123)


### Fig 1a ###

# Set figure parameters (MODIFY HERE)
num.pts <- 30 # number of points to plot for 1a; approximately twice as many points will be plotted in 1b (to maintain similar density of points)
eps <- 0.1 # set the spread of points around the trade-off
margins <- 0.05 # set the amount of whitespace between constraints and points 
# END MODIFY

# Simulate datapoints
x <- seq(0, 1, length.out = num.pts) # generate points along the line y = 1 - x
y <- 1 - x

df <- data.frame(x = x + eps * runif(num.pts, -1, 1), y = y + eps * runif(num.pts, -1, 1)) # add some random noise to datapoints
constraint <- max(rowSums(df)) # this sets the location of the constraint (blue dashed) line
bottom.cutoff <- min(rowSums(df)) # this sets the location of the "excluded region" (red triangle)

# Build plot
p1 <- ggplot(df) + 
  aes(x = x, y = y) + 
  geom_point(size = 5, alpha = 0.5) + 
  theme_classic() + 
  xlab("Reproductive effort") + 
  ylab("Survival") + 
  scale_x_continuous(limits = c(0,1), breaks = c(0,1), labels = c("Low", "High")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("Low", "High")) + 
  geom_abline(slope = -1, intercept = constraint + margins, # add constraint (blue dashed) line just above points
              linetype = "dashed", size = 1.5, color = "blue") + 
  geom_polygon(data = data.frame(x = c(0, 0, bottom.cutoff - margins), # add excluded region triangle just below points
                                 y = c(0, bottom.cutoff - margins, 0)),
              aes(x = x, y = y), 
              alpha = 0.1, fill = "red", colour = "red", size = 1.2) + 
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

#show(p1)
               

### Fig 1b ###

# requires "constraint" from 1a (above)

# Simulate datapoints 
x <- runif(4 * num.pts, 0, 1) # generate datapoints uniformly distributed on 1 x 1 square 
y <- runif(4 * num.pts, 0, 1)
permitted <- which(y < constraint - x) # reject points below the constraint line

df <- data.frame(x = x, y = y)
df <- df[permitted, ] # "bad" points are excluded here

# Build plot
p2 <- ggplot(df) + aes(x = x, y = y) + 
  geom_point(size = 5, alpha = 0.5) + 
  theme_classic() + 
  xlab("Reproductive effort") + 
  ylab("Survival") + 
  scale_x_continuous(limits = c(0,1), breaks = c(0,1), labels = c("Low", "High")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("Low", "High")) + 
  geom_abline(slope = -1, intercept = constraint + margins,  # add constraint (blue dashed) line just above points
              linetype = "dashed", size = 1.5, color = "blue") + 
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size=14))

#show(p2)

### Plot together and print to pdf ###

pdf("Digest_figure.pdf", height = 4, width = 9)
grid.arrange(p1, p2, ncol = 2)
dev.off()

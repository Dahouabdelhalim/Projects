par(mfrow = c(2,3))
par(mar = c(2, 2, 1, 0))

#https://stackoverflow.com/questions/37071001/barplot-greek-letters-on-y-axis-in-r


x <- seq(-0.6,.6, length=10)
color.names = new.pal[c(1,2,3)]



# 1. Autumn, Cost per thread####
data <- matrix(c(.37,-.37,.07,-.07,0,0,0,0), ncol = 4)

head(data)
rownames(data) <- c('+ uncertainty','- uncertainty')
colnames(data) <- c('SE','Energetic \\n conversion
                    coefficient (C.F.)', 
                    'Metabolic \\n cost
                    (b)',
                    'Intake \\nexponent (d)')
names <- c('SE','C.F.', 
           'b',
           'd')
colnames(data) <- names

barplot(data[1,], horiz = T, las=1, xlim = c(-0.5,.5), xaxt='n', ylab='',
        beside=T, col=c('black'))
barplot(data[2,], horiz = T, las=1, xlim = c(-0.5,.5), xaxt='n',
        yaxt='n',                 #To prevent double printing of y-labels.
        beside=T, col=c('white'), add = TRUE)
axis(1, at=pretty(x),  lab=paste0(pretty(x)+1), las=TRUE)
#title(main = expression("Cost per thread ("*italic(h)*", J/thread)"))

# 2. Autumn, Proportion of energy budget####

head(data)
rownames(data) <- c('+ uncertainty','- uncertainty')
colnames(data) <- c('SE','Energetic \\n conversion
                    coefficient (C.F.)', 
                    'Metabolic \\n cost
                    (b)',
                    'Intake \\nexponent (d)')

colnames(data) <- c('SE','C.F.', 
                    'b',
                    'd')

new.pal
#Order SE -> d
#never -> daily
#high val
data_pos <- matrix(c(.18,.148,.081,.02,.02,.02,-.062,-.062,-.062,.001,.002,.003), ncol = 4)
#low val
data_neg <- matrix(c(-.18,-.148,-.081,-.022,-.022,-.022,.071,.071,.071,-.001,-0.002,-0.003), ncol = 4)
#color.names = c("dark blue","light blue","forest green") #Columns are the 3-4 bars in one grouping
colnames(data_pos) <- color.names
colnames(data_neg) <- color.names

barplot(data_pos,beside=T, horiz = T, las=1, xlim = c(-0.5,.5),xaxt='n',
        yaxt = 'n', # repetitive in panel graph so not including y-axis
        ylab='', col=color.names, axis.lty="solid")
barplot(data_neg,beside=T, horiz = T, xlim = c(-0.5,.5),xaxt='n',
        yaxt='n', col="white", 
        #axis.lty="solid", 
        add = TRUE)
axis(1, at=pretty(x),  lab=paste0(pretty(x)+1), las=TRUE)
#title(main = expression("Proportion \\nof energy budget"))

# 3. Autumn, Proportion of cost ####

head(data)
rownames(data) <- c('+ uncertainty','- uncertainty')
colnames(data) <- c('SE','Energetic \\n conversion
                    coefficient (C.F.)', 
                    'Metabolic \\n cost
                    (b)',
                    'Intake \\nexponent (d)')

colnames(data) <- c('SE','C.F.', 
                    'b',
                    'd')

new.pal
#Order SE -> d
#never -> daily
#high val
data_pos <- matrix(c(.132,.080,.034,.051,.038,.023,-.142,-.109,-.071,0,0,0), ncol = 4)
#low val
data_neg <- matrix(c(-.132,-.08,-.034,-.054,-.040,-.025,.200,.142,.083,0,0,0), ncol = 4)
#color.names = c("dark blue","light blue","forest green") #Columns are the 3-4 bars in one grouping
colnames(data_pos) <- color.names
colnames(data_neg) <- color.names

barplot(data_pos,beside=T, horiz = T, las=1, xlim = c(-0.5,.5),xaxt='n',
        yaxt = 'n', # repetitive in panel graph so not including y-axis
        ylab='', col=color.names, axis.lty="solid")
barplot(data_neg,beside=T, horiz = T, xlim = c(-0.5,.5),xaxt='n',
        yaxt='n', col="white", 
        #axis.lty="solid", 
        add = TRUE)
axis(1, at=pretty(x),  lab=paste0(pretty(x)+1), las=TRUE)
#title(main = expression("Proportion \\nof cost"))

# 4. Spring, Cost per thread####
data <- matrix(c(.34,-.34,.07,-.14,0,0,0,0), ncol = 4)

head(data)
rownames(data) <- c('+ uncertainty','- uncertainty')
colnames(data) <- c('SE','Energetic \\n conversion
                    coefficient (C.F.)', 
                    'Metabolic \\n cost
                    (b)',
                    'Intake \\nexponent (d)')
names <- c('SE','C.F.', 
           'b',
           'd')
colnames(data) <- names

barplot(data[1,], horiz = T, las=1, xlim = c(-0.5,.5), xaxt='n', ylab='',
        beside=T, col=c('black'))
barplot(data[2,], horiz = T, las=1, xlim = c(-0.5,.5), xaxt='n',
        yaxt='n',                 #To prevent double printing of y-labels.
        beside=T, col=c('white'), add = TRUE)
axis(1, at=pretty(x),  lab=paste0(pretty(x)+1), las=TRUE)
#title(main = expression("Cost per thread ("*italic(h)*", J/thread)"))

# 5. Spring, Proportion of energy budget####

head(data)
rownames(data) <- c('+ uncertainty','- uncertainty')
colnames(data) <- c('SE','Energetic \\n conversion
                    coefficient (C.F.)', 
                    'Metabolic \\n cost
                    (b)',
                    'Intake \\nexponent (d)')

colnames(data) <- c('SE','C.F.', 
                    'b',
                    'd')

new.pal
#Order SE -> d
#never -> daily
#high val
data_pos <- matrix(c(0.031,0.144,0.083,0.028,0.028,.028,-0.097,-0.097,-0.097,.0005,0.0013,-.0001), ncol = 4)
#low val
data_neg <- matrix(c(-0.031,-0.144,-0.083,-0.031,-0.031,-0.031,0.120,0.120,0.120,-.0005,-0.0013,+0.0001), ncol = 4)
#color.names = c("dark blue","light blue","forest green") #Columns are the 3-4 bars in one grouping
colnames(data_pos) <- color.names
colnames(data_neg) <- color.names

barplot(data_pos,beside=T, horiz = T, las=1, xlim = c(-0.5,.5),xaxt='n',
        yaxt = 'n', # repetitive in panel graph so not including y-axis
        ylab='', col=color.names, axis.lty="solid")
barplot(data_neg,beside=T, horiz = T, xlim = c(-0.5,.5),xaxt='n',
        yaxt='n', col="white", 
        #axis.lty="solid", 
        add = TRUE)
axis(1, at=pretty(x),  lab=paste0(pretty(x)+1), las=TRUE)
#title(main = expression("Proportion \\nof energy budget"))

# 6. Spring, Proportion of cost ####

head(data)
rownames(data) <- c('+ uncertainty','- uncertainty')
colnames(data) <- c('SE','Energetic \\n conversion
                    coefficient (C.F.)', 
                    'Metabolic \\n cost
                    (b)',
                    'Intake \\nexponent (d)')

colnames(data) <- c('SE','C.F.', 
                    'b',
                    'd')

new.pal
#Order SE -> d
#never -> daily
#high val
data_pos <- matrix(c(0.049,0.126,0.049,0.069,0.058,0.042,-0.203,-0.178,-0.136,0,0,0), ncol = 4)
#low val
data_neg <- matrix(c(-0.049,-0.126,-0.049,-0.070,-0.060,-0.044,0.342,0.278,0.188,0,0,0), ncol = 4)
#color.names = c("dark blue","light blue","forest green") #Columns are the 3-4 bars in one grouping
colnames(data_pos) <- color.names
colnames(data_neg) <- color.names

barplot(data_pos,beside=T, horiz = T, las=1, xlim = c(-0.5,.5),xaxt='n',
        yaxt = 'n', # repetitive in panel graph so not including y-axis
        ylab='', col=color.names, axis.lty="solid")
barplot(data_neg,beside=T, horiz = T, xlim = c(-0.5,.5),xaxt='n',
        yaxt='n', col="white", 
        #axis.lty="solid", 
        add = TRUE)
axis(1, at=pretty(x),  lab=paste0(pretty(x)+1), las=TRUE)
#title(main = expression("Proportion \\nof cost"))




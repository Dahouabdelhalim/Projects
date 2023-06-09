

# Some replication code
# make sure you are using the right path
load("EditorialData.Rda")

dta <- ReplicationData   # use smaller dataset name, it is more user-friendly but not necessary

# Figure not shown in the paper: Number of papers published per volume
barplot(table(dta$Vol), xlab="Volume", ylab="Number of papers published per volume")
# Comment: Notice that some data are missing because some papers have not yet been assigned to a volume/issue (so there are missing data)

# An example to generate Figure 1 (the other figures are generated using similar code)
d <- dta 
d[is.na(d$Vol),]$Vol <- 26 # I assign all papers with a missing volume, to an imaginary volume '26' (this is numeric and will later be recoded)
t <- table(d$Vol, d$Paradigm); t <- prop.table(t, 1)  * 100 ;  t <- as.data.frame(t) # I prepare a table with data to be visualized


t$Var1 <- recode_factor(t$Var1, "26" = "25+") # I change "26" to "25+" because it looks better (these are papers which have been accepted but not assigned to a volume yet)


level_order <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '25+') # this is necessary, otherwise some R versions will not plot the '25+' category in the right order


tiff("Editorial_Figure1.tiff", units="in", width=8, height=4, res=300)  # This is to save the Figure

ggplot(t %>% filter(Var2 %in% c("Qualitative","Quantitative")), aes(x = factor(Var1, level = level_order), y = Freq , linetype = Var2, group=Var2)) +   geom_line() + xlab("Time (Volumes)") + ylab("% of items") + scale_fill_grey() + theme_bw() + theme(legend.position="none") +  geom_text(x=18, y=75, label="Qualitative", color="black") +  geom_text(x=18, y=6.5, label="Quantitative" , color="black")  # This is to plot the figure

dev.off() # close the file
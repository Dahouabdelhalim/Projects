
#### Script 3 - Proportions of journals that fulfill the criteria based on compliance language 
#### Get proportions of journals that fulfill a given criterion and plot based on compliance language (i.e., requirement(3), recommendation(2), suggestions(1) or absent (0))

# Import data from existing script (Relative Path)
source("Script1_Data-preparation.R")

## Include the line of code below if interested only in journals that have animal care policies (i.e., excluding journals that have no policies)
# regression_data <- regression_data[!(regression_data$final_any_ac_state_YN=="0"),]

## Proportion tables, where: 0 = None; 1 = May; 2 = Should; 3 = Must. 

# Does the journal have any statement related to animal care?

options(digits=2)
transform(as.data.frame(table(regression_data$final_any_ac_state_YN)), percentage_column=Freq/nrow(regression_data)*100)

# Does the journal have a statement related to best practices for field work?
options(digits=2)
transform(as.data.frame(table(regression_data$final_field_YN)), percentage_column=Freq/nrow(regression_data)*100)

# Does the journal require that authors state that an animal care permit or approval was granted?
options(digits=2)
transform(as.data.frame(table(regression_data$final_ac_approve)), percentage_column=Freq/nrow(regression_data)*100)

# Does the journal require that authors specify the institutional authority that granted animal care licenses or permits?
options(digits=2)
transform(as.data.frame(table(regression_data$final_ac_inst)), percentage_column=Freq/nrow(regression_data)*100)

# Does the journal require that authors provide animal care license or permit numbers or documentation?
options(digits=2)
transform(as.data.frame(table(regression_data$final_ac_num)), percentage_column=Freq/nrow(regression_data)*100)

# Does the journal have a statement encouraging authors to adopt the 3Rs?
options(digits=2)
transform(as.data.frame(table(regression_data$final_three_r)), percentage_column=Freq/nrow(regression_data)*100)

# Does the journal have a statement informing authors that adherence to relevant animal care journal policies is a condition of publication? 
options(digits=2)
transform(as.data.frame(table(regression_data$condition)), percentage_column=Freq/nrow(regression_data)*100)


regression_data<-regression_data%>%mutate(final_ac_approve=as.numeric(final_ac_approve),
                                          final_ac_inst=as.numeric(final_ac_inst),
                                          final_ac_num=as.numeric(final_ac_num),
                                          final_three_r=as.numeric(final_three_r),
                                          condition=as.numeric(condition))

# Combine above
plot_data<-transform(as.data.frame(table(regression_data$final_ac_approve)), percentage_column=Freq/nrow(regression_data)*100, criterion="final_ac_approve")%>%
  bind_rows(transform(as.data.frame(table(regression_data$final_ac_inst)), percentage_column=Freq/nrow(regression_data)*100, criterion="final_ac_inst"))%>%
  bind_rows(transform(as.data.frame(table(regression_data$final_ac_num)), percentage_column=Freq/nrow(regression_data)*100, criterion="final_ac_num"))%>%
  bind_rows(transform(as.data.frame(table(regression_data$final_three_r)), percentage_column=Freq/nrow(regression_data)*100, criterion="final_three_r"))%>%
  bind_rows(transform(as.data.frame(table(regression_data$condition), stringsAsFactors=FALSE), percentage_column=Freq/nrow(regression_data)*100, criterion="condition"))


## Plot Proportions 
plot_data<-plot_data%>%
  mutate(position=as.numeric(Var1),
         linecolour=Var1,
         linecolour=ifelse(linecolour=="0", "grey70", ifelse(linecolour=="1", "grey60", ifelse(linecolour=="2", "grey45", "black"))),
         mustmayshould=ifelse(Var1=="0", "None", ifelse(Var1=="1", "May", ifelse(Var1=="2", "Should", "Must"))),
         position=ifelse(criterion=="final_three_r", position+6.8, ifelse(criterion=="final_ac_num", position+12.6, ifelse(criterion=="final_ac_inst", position+18.7, ifelse(criterion=="final_ac_approve", position+24.9, position+0.9)))))

plot_criteria<-function(criterion_sub){
  plot(1, 1, xlim = c(0, 100) , ylim = c(0, 1+max(criterion_sub$position)), xlab = "", ylab = "", type = "n",
       axes = FALSE)
  lines_in_plot<-group_by(criterion_sub, criterion)%>%summarize(min=min(position), max=max(position))
  for (i in c(20, 40, 60, 80)) with(lines_in_plot, segments(i, min-0.5, i, 2+max, lty=4, lwd=0.6))
  with(criterion_sub, segments(0, position, percentage_column, position, lwd = 1,  col=criterion_sub$linecolour))
  with(criterion_sub, points(percentage_column, position, pch=20, lwd=0.5,  col=criterion_sub$linecolour))
  axis(2, at=(criterion_sub$position), labels=criterion_sub$mustmayshould, line=-1, las=1, cex.axis=0.65, tick=FALSE, lwd=3)
}


tiff("Proportions plot.tif", height=5, width=7, units="in", res=300)
par(mar=c(0.02,14,0.02,0), oma=c(1,1,1,1))
plot_criteria(plot_data)
axis(1, line=-1)
mtext(side=2, "Animal care compliance a condition of publication", at=1.4+max(subset(plot_data,criterion=="condition")$position) , las=1.5, font=1.5, cex = 0.8)
mtext(side=2, "Adopt the 3Rs", at=1+max(subset(plot_data,criterion=="final_three_r")$position) , las=1.5, font=1.5, cex = 0.8)
mtext(side=2, "Provide animal care approval documentation", at=1.4+max(subset(plot_data,criterion=="final_ac_num")$position) , las=1.5, font=1.5, cex = 0.8)
mtext(side=2, "Identify animal care approval granting authority", at=1.5+max(subset(plot_data,criterion=="final_ac_inst")$position) , las=1.5, font=1.5, cex = 0.8)
mtext(side=2, "Declare whether animal care approval granted", at=1.65 +max(subset(plot_data,criterion=="final_ac_approve")$position) , las=1.5, font=1.5, cex = 0.8)


# title("Percent of journals that fulfill criteria based on compliance language", outer=TRUE, cex.main=1, adj=0)
dev.off() #close image

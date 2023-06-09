
scaler=62.7 #pixel/cm conversion

reps=1:10
days=4:30

pdf(file='/PATH/TO/OUTPUT/DIRECTORY/acute_qc_heights.pdf', width=5.5, height=16.5)

layout(matrix(1:10, 5,2, byrow=TRUE))

outliers<-c()

for (line in c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100')){ #

	acute_heights<-c()

	qc_data<-read.csv('/PATH/TO/DATA/FILE/ground_truthing_heights.csv')

	anno_data<-read.csv(paste('/PATH/TO/DATA/FILE/',line,'_homology_groups_annotated.csv',sep=''))

	anno_base<-colMeans(anno_data[anno_data[,1]=='base',4:5])

	if(line=='A10'){reps=1:10; bounds=c(0,12)}
	if(line=='RIL39'){reps=c(1:10); bounds=c(0,25)}
	if(line=='RIL159'){reps=c(1:10); bounds=c(0,25)}
	if(line=='RIL159'){reps=c(1:6,9:10); bounds=c(0,25)}
	if(line=='B100'){reps=c(1:5); bounds=c(0,25)}

	for (rep in reps){

		rep_data<-anno_data[grepl(paste('_rep',rep,'_',sep=''), anno_data[,2]),]
		rep_base<-colMeans(rep_data[rep_data[,1]=='base',4:5])

		if(sum(is.na(rep_base))>1){
			rep_base<-anno_base
		}

		for (day in days){

			day_data<-rep_data[grepl(paste('_d',day,'_',sep=''), rep_data[,2]),]
			day_ligs<-day_data[grepl('lig', day_data[,1]),]

			qc_row<-qc_data[qc_data[,1]==line & qc_data[,2]==rep & qc_data[,3]==day,]
			
			if(nrow(day_ligs)>0 & nrow(qc_row)>0){

				max_lig_hgt_pos<-day_ligs[which.min(day_ligs[,5]),4:5]
				
				lfno<-as.numeric(strsplit(as.character(sort(day_ligs[,1])[nrow(day_ligs)]),'lig')[[1]][2])

				max_lig_hgt<-sqrt((max_lig_hgt_pos[1]-rep_base[1])^2+(max_lig_hgt_pos[2]-rep_base[2])^2)/scaler

				acute_heights<-rbind(acute_heights, c(line, rep, day, max_lig_hgt, qc_row[,4]/scaler, lfno, qc_row[,5]+1))
			}
		}
	}

	colnames(acute_heights)<-c('line', 'rep', 'day', 'lig_hgt', 'QC_lig_hgt', 'lfno', 'QC_lfno')

	plot(acute_heights[,c(4,5)], xlim=bounds, ylim=bounds, main=paste(line, 'Maximum Ligular Height'), xlab='Estimated height (cm)', ylab='Groundtruthed height (cm)', pch=19, cex=0.45)
	legend('topleft', lty=c(3,2,0), pch=c(-1,-1,19), col=c('black', 'blue', 'red'), legend=c('Ideal Fit (Perfect Correlation)', 'Realized Standard ~ Estimate Correlation', '>0.5 Fold Change Outlier'), cex=0.6)

	max_hgt_lm<-lm(as.numeric(acute_heights[,5])~as.numeric(acute_heights[,4]))

	abline(0, 1, lty=3); abline(max_hgt_lm, lty=2, col='blue')
	#lines(seq(bounds[1], bounds[2], by=0.1), seq(bounds[1], bounds[2], by=0.1)*0.5, col='red', lty=2)
	#lines(seq(bounds[1], bounds[2], by=0.1), seq(bounds[1], bounds[2], by=0.1)*1.5, col='red', lty=2)

	text(bounds[2]*0.2, bounds[2]*0.7, label=paste("R-squared", round(summary(max_hgt_lm)$r.squared,3)), col='blue', cex=0.9)
	text(bounds[2]*0.2, bounds[2]*0.6, label=paste("Slope", round(coef(max_hgt_lm)[2],3)), col='blue', cex=0.9)

	fold<-as.numeric(acute_heights[,5])/as.numeric(acute_heights[,4])

	points(acute_heights[fold<0.5 | fold>1.5,4], acute_heights[fold<0.5 | fold>1.5,5], pch=19, cex=0.25, col='red')

	hist(fold, xlim=c(0,4), breaks=seq(0,4,0.05), xlab='Dividend Ratio', main='Standard/Estimate Fold Change')
	abline(v=c(0.5,1.5), lty=2, col='red')

	outliers<-rbind(outliers, acute_heights[fold<0.5 | fold>1.5,])

}

dev.off()

print(outliers)




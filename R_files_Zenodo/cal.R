#performs a linear calibration on a subset of data, then validation on a subset of witheld data.

#x is a vector or matrix of dependent variables
#y is a vector of the independent variable
#prct.cal is the fraction of data to use in the calibration (between 0 and 1)
#n.it is the number of random assignments to calibration/validation datasets

cal.val<-function(x,y,prct.cal,n.it) {
	cal.valRMSE<-rep(0,n.it)
	cal.valCE<-rep(0,n.it)
	cal.valr2<-rep(0,n.it)
	cal.valBIC<-rep(0,n.it)
	cal.valint<-rep(0,n.it)
	if(is.vector(x)==TRUE) {
		cal.valcoeff<-rep(0,n.it)
		} else {
		cal.valcoeff<-matrix(,nrow=n.it,ncol=ncol(x))
		}
	
	for (j in 1:n.it) {
	n.cal<-round(length(y)*prct.cal)															#number of records in calibration period
	n.val<-length(y)-n.cal
	rndm<-sample(seq(1:length(y)))																#randomly order data
	cal.i<-y[rndm[1:n.cal]]																		#select the first n.cal for calibration
	val.i<-y[rndm[(n.cal+1):length(y)]]
	if(is.vector(x)==TRUE) {
		cal.d<-x[rndm[1:n.cal]]
		val.d<-x[rndm[(n.cal+1):length(y)]]	
	} else {
		cal.d<-x[rndm[1:n.cal],]
		val.d<-x[rndm[(n.cal+1):length(y)],]
	}
	pred<-rep(0,length(val.i))																	#empty matrix to be filled with predicted values
	
	if(is.vector(x)==TRUE) {
	fit<-lm(cal.i~cal.d)
	pred<-summary(fit)$coeff[1,1]+summary(fit)$coeff[2,1]*val.d
	} else {
		if(ncol(x)==2) {
		fit<-lm(cal.i~cal.d[,1]+cal.d[,2])
		pred<-summary(fit)$coeff[1,1]+summary(fit)$coeff[2,1]*val.d[,1]+summary(fit)$coeff[3,1]*val.d[,2]
		} else {
			if(ncol(x)==3) {
				fit<-lm(cal.i~cal.d[,1]+cal.d[,2]+cal.d[,3])
				pred<-summary(fit)$coeff[1,1]+summary(fit)$coeff[2,1]*val.d[,1]+summary(fit)$coeff[3,1]*val.d[,2]+summary(fit)$coeff[4,1]*val.d[,3]	
				} else {
					if(ncol(x)==4) {
						fit<-lm(cal.i~cal.d[,1]+cal.d[,2]+cal.d[,3]+cal.d[,4])
						pred<-summary(fit)$coeff[1,1]+summary(fit)$coeff[2,1]*val.d[,1]+summary(fit)$coeff[3,1]*val.d[,2]+summary(fit)$coeff[4,1]*val.d[,3]+summary(fit)$coeff[5,1]*val.d[,4]
						} else {
							print("You picked more than 4 variables and need to edit the code")
						}
				}
		}
}

cal.valr2[j]<-summary(fit)$r.squared
cal.valBIC[j]<-BIC(fit)
cal.valint[j]<-summary(fit)$coeff[1,1]
if(is.vector(x)==TRUE) {
	cal.valcoeff[j]<-summary(fit)$coeff[2,1]
	} else {
	cal.valcoeff[j,]<-summary(fit)$coeff[2:(ncol(x)+1),1]
	}
cal.valRMSE[j]<-sqrt(sum((pred-val.i)^2,na.rm=T)/n.val)
cal.valCE[j]<-1-(sum((pred-val.i)^2,na.rm=T)/sum((mean(val.i)-val.i)^2,na.rm=T))

}
cal.valsummary<-matrix(,nrow=9,ncol=2)
cal.valsummary[1,1]<-mean(cal.valr2)
cal.valsummary[1,2]<-sd(cal.valr2)
cal.valsummary[2,1]<-mean(cal.valBIC)
cal.valsummary[2,2]<-sd(cal.valBIC)
cal.valsummary[3,1]<-mean(cal.valint)
cal.valsummary[3,2]<-sd(cal.valint)
if(is.vector(x)==TRUE) {
	cal.valsummary[4,1]<-mean(cal.valcoeff)
	cal.valsummary[4,2]<-sd(cal.valcoeff)
	} else {
	cal.valsummary[4:(ncol(x)+3),1]<-colMeans(cal.valcoeff)
	cal.valsummary[4:(ncol(x)+3),2]<-apply(cal.valcoeff,2,sd)
	}
cal.valsummary[8,1]<-mean(cal.valRMSE)
cal.valsummary[8,2]<-sd(cal.valRMSE)
cal.valsummary[9,1]<-mean(cal.valCE)
cal.valsummary[9,2]<-sd(cal.valCE)

cal.valsummary<<-cal.valsummary

}
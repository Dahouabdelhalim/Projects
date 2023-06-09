
#### Load the data set from forecast_data.rds and get some of its basic features

LIST_forecast=readRDS("forecast_data.rds")

## Number of geographic entities available across time

nb=simplify2array(lapply(LIST_forecast,length))
plot(as.Date(names(LIST_forecast)),nb,ylim=c(0,max(nb)),xlab="Date",ylab="Number of geographic entities",axes=FALSE)
xdates=c("2020-04-01","2020-10-01","2021-04-01","2021-10-01")
axis(1,at=as.Date(xdates),xdates)
axis(2)

## Available dates at which the forecasts were run

names(LIST_forecast)

## Names of available geographic entities at a specific date at which forecasts were run

names(LIST_forecast[['2020-10-01']])


#### Extract forecasts corresponding to specific dates at which the forecasts were run and to specific geographic entities

## Extraction function
## Arguments of the extract() function
## LIST: object with the same structure than LIST_forecast (i.e., a list of lists of items)
## geographic_entity: vector of strings corresponding to the geographic entites for which forecasts are requested (any geographic entity for which the forecast is not available at the specified date(s) is ignored)
## date: vector of strings corresponding to dates at which the forecasts are requested; default is NULL, which means that all available dates are included (any date not available in the series of forecasts is ignored)
## Value of the extract function: a list of lists with the same structure than LIST_forecast (i.e., a list of lists of items) 

extract=function(LIST,geographic_entity,date=NULL){
	res=NULL
	if(is.null(date)){
		date=names(LIST)
	}
	for(i in 1:length(LIST)){
		if(names(LIST)[i]%in%date){
			res.i=NULL
			names.g=names(LIST[[i]])
			k=which(names.g%in%geographic_entity)
			if(length(k)>0){
				for(m in 1:length(k)){
					res.i=c(res.i,list(LIST[[i]][[k[m]]]))
					names(res.i)[length(res.i)]=names.g[k[m]]
				}
			}
			res=c(res,list(res.i))
			names(res)[length(res)]=names(LIST)[i]
		}
	}
	return(res)
}

## Examples of forecast extraction and plot

forecast1=extract(LIST_forecast,"France",date="2020-12-01")

forecast2=extract(LIST_forecast,c("France","US"),date="2020-12-01")

forecast3=extract(LIST_forecast,"France",date=c("2020-12-01","2021-02-01"))

dates=c("2020-12-01","2021-02-01")
GEs=c("France","US")

forecast4=extract(LIST_forecast,GEs,dates)

par(mfrow=c(2,2))
for(fc.date in dates){
	for(fc.GE in GEs){
		fc=as.data.frame(forecast4[[fc.date]][[fc.GE]])
		time=as.Date(fc$time,origin=as.Date("2020-01-01"))
	plot(time,fc$death,xlab="Date",ylab="Cumulative mortality",main=paste("Forecast drawn on",fc.date,"for",fc.GE))
	lines(time,fc$M0.forecast)
	lines(time,fc$M0.Qinf,lty="dashed")
	lines(time,fc$M0.Qsup,lty="dashed")
	}
}

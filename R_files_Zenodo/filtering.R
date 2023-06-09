#Select only entries that have label as a value
filter_field_label <- function(data,column,label){
	newdata <- data[ which(data[, substr(names(data), 1, 4) == column]==label), ]

	return(newdata)
}


filter_finished <- function(data){
  newdata <- data[ which(data$FINISHED==TRUE), ]
  
  return(newdata)
}

one.bootstrap.replicate <- function(to.bootstrap.long = droplevels(attplot), to.bootstrap.wide = droplevels(attplot.wide),
	studies = unique(attplot.wide$UI)){
	
	n.studies <- length(studies)

	r.studies <- sample(studies, length(studies), replace = TRUE)
	r.studyXdriver <- data.frame(UI = rep(r.studies, each=nlevels(to.bootstrap.long[,3])), IPBES.direct.driver = rep(levels(to.bootstrap.long[,3]), n.studies))
	long.replicate <- merge(r.studyXdriver, to.bootstrap.long, by.x = c("UI", "IPBES.direct.driver"), by.y = c("UI", "IPBES.direct.driver"), all.x=TRUE, all.y=FALSE)
	
	matching.rows <- match(r.studies, to.bootstrap.wide$UI)
	wide.replicate <- to.bootstrap.wide[matching.rows,]

	return(list(long = long.replicate, wide = wide.replicate))
	}



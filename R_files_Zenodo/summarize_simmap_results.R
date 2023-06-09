ntree <- 250  	## 250 BEAST trees evaluated
nsim <- 10     	## 10 character simulations per tree
breakpoint <- 3 	## breakpoint age between "old" & "young" splits
min_age <- 0 	## don't consider splits younger than this

simmaps.summary <- matrix(data=rep(0,ntree*nsim*4),nrow=ntree*nsim,ncol=4,dimnames=list(1:(ntree*nsim),c("h_old","h_young","r_old","r_young")))

for (tr in 1:(ntree*nsim)) {

	## HOST get a list of speciation events associated with host shift or not (from single host-mapped tree)
	root.depth <- max(node.depth.edgelength(host.simmaps[[tr]]))
	node.ages <- root.depth - node.depth.edgelength(host.simmaps[[tr]])
	nodes <- (Ntip(host.simmaps[[tr]])+1):(Ntip(host.simmaps[[tr]])+Nnode(host.simmaps[[tr]]))
	age <- node.ages[nodes]
	shift <- vector(mode="logical",length=Nnode(host.simmaps[[tr]]))
	spec.list.host <- data.frame(age,shift,row.names=nodes)

	for (node in nodes) {
		desc.edges <- which(host.simmaps[[tr]]$edge[,1]==node)
		desc.host1 <- names(host.simmaps[[tr]]$maps[[desc.edges[1]]])[length(names(host.simmaps[[tr]]$maps[[desc.edges[1]]]))] # get the name of the last host inferred along the branch
		desc.host2 <- names(host.simmaps[[tr]]$maps[[desc.edges[2]]])[length(names(host.simmaps[[tr]]$maps[[desc.edges[2]]]))]
		desc.tissue1 <- names(tissue.simmaps[[tr]]$maps[[desc.edges[1]]])[length(names(tissue.simmaps[[tr]]$maps[[desc.edges[1]]]))]
		desc.tissue2 <- names(tissue.simmaps[[tr]]$maps[[desc.edges[2]]])[length(names(tissue.simmaps[[tr]]$maps[[desc.edges[2]]]))]
		## if((desc.host1 != desc.host2)|(desc.tissue1 != desc.tissue2)) {spec.list.host[node-Ntip(host.simmaps[[tr]]),2] <- TRUE} ## use this to count shift in host species OR host tissue as a host shift (probably subject to bias due to polymorphic states)
		if(desc.host1 != desc.host2) {spec.list.host[node-Ntip(host.simmaps[[tr]]),2] <- TRUE} ## use this to count ounly shift in host species as a host shift
		}

	## GEOGRAPHY get a list of speciation events associated with regional dispersal or not (from single region-mapped tree)
	spec.list.region <- data.frame(age,shift,row.names=nodes)

	for (node in nodes) {
		desc.edges <- which(region.simmaps[[tr]]$edge[,1]==node)
		desc.region1 <- names(region.simmaps[[tr]]$maps[[desc.edges[1]]])[length(names(region.simmaps[[tr]]$maps[[desc.edges[1]]]))]
		desc.region2 <- names(region.simmaps[[tr]]$maps[[desc.edges[2]]])[length(names(region.simmaps[[tr]]$maps[[desc.edges[2]]]))]
		if(desc.region1 != desc.region2) {spec.list.region[node-Ntip(region.simmaps[[tr]]),2] <- TRUE}
		}

	### getting barplot stats for single mapped tree 
	h.shift.old <- spec.list.host[which(spec.list.host[,2] & spec.list.host[,1] > breakpoint),1]
	h.shift.young <- spec.list.host[which(spec.list.host[,2] & spec.list.host[,1] < breakpoint & spec.list.host[,1] > min_age),1]
	h.noshift.old <- spec.list.host[which(!spec.list.host[,2] & spec.list.host[,1] > breakpoint),1]
	h.noshift.young <- spec.list.host[which(!spec.list.host[,2] & spec.list.host[,1] < breakpoint & spec.list.host[,1] > min_age),1]
	freq.h.shift.old <- length(h.shift.old) / ((length(h.shift.old) + length(h.noshift.old)))
	freq.h.shift.young <- length(h.shift.young) / ((length(h.shift.young) + length(h.noshift.young)))

	r.shift.old <- spec.list.region[which(spec.list.region[,2] & spec.list.region[,1] > breakpoint),1]
	r.shift.young <- spec.list.region[which(spec.list.region[,2] & spec.list.region[,1] < breakpoint & spec.list.region[,1] > min_age),1]
	r.noshift.old <- spec.list.region[which(!spec.list.region[,2] & spec.list.region[,1] > breakpoint),1]
	r.noshift.young <- spec.list.region[which(!spec.list.region[,2] & spec.list.region[,1] < breakpoint & spec.list.region[,1] > min_age),1]
	freq.r.shift.old <- length(r.shift.old) / ((length(r.shift.old) + length(r.noshift.old)))
	freq.r.shift.young <- length(r.shift.young) / ((length(r.shift.young) + length(r.noshift.young)))

	simmaps.summary[tr,] <- c(freq.h.shift.old,freq.h.shift.young,freq.r.shift.old,freq.r.shift.young)
	}

boxplot(simmaps.summary,col="gray",boxwex=0.5,staplewex=0.8,ylab="proportion of splits",main="Proportion of splits, by age,\\nwith differences in host and region")

rm(ntree,tr,nsim,nodes,node,min_age,node.ages,r.shift.old,r.shift.young,r.noshift.old,r.noshift.young,freq.r.shift.old,freq.r.shift.young,freq.h.shift.old,freq.h.shift.young)
rm(desc.edges,desc.host1,desc.host2,desc.region1,desc.region2,desc.tissue1,desc.tissue2,shift,age,breakpoint,h.noshift.old,h.noshift.young,h.shift.old,h.shift.young)
rm(root.depth,spec.list.host,spec.list.region)


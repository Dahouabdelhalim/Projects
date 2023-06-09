## do an MPR character reconstruction

host.MPR <- MPR(bleph.char[bleph.out.tree$tip.label,"host"],bleph.out.tree,"Ble141")
#plot(bleph.out.tree,cex=0.5,no.margin=TRUE)
#nodelabels(host.MPR[,"upper"],cex=0.5)

region.MPR <- MPR(bleph.char[bleph.out.tree$tip.label,"region"],bleph.out.tree,"Ble141")
#plot(bleph.out.tree,cex=0.5,no.margin=TRUE)
#nodelabels(region.MPR[,"upper"],cex=0.5)

tissue.MPR <- MPR(bleph.char[bleph.out.tree$tip.label,"tissue"],bleph.out.tree,"Ble141")
#plot(bleph.out.tree,cex=0.5,no.margin=TRUE)
#nodelabels(tissue.MPR[,"upper"],cex=0.5)

## get numbers from MPR reconstruction of actual character states
breakpoint <- 3

min_age <- 0
phy <- bleph.out.tree

root.depth <- max(node.depth.edgelength(phy))
node.ages <- root.depth - node.depth.edgelength(phy)
nodes <- (Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
age <- node.ages[nodes]
shift <- vector(mode="logical",length=Nnode(phy))
spec.list.host <- data.frame(age,shift,row.names=nodes)
spec.list.region <- data.frame(age,shift,row.names=nodes)
for (node in nodes) {
	desc.edges <- phy$edge[which(phy$edge[,1]==node),2]
	if (desc.edges[1] > Ntip(phy)) {
		desc.host1 <- host.MPR[which(rownames(host.MPR)==desc.edges[1]),"lower"]
		desc.region1 <- region.MPR[which(rownames(region.MPR)==desc.edges[1]),"lower"]
		desc.tissue1 <- tissue.MPR[which(rownames(tissue.MPR)==desc.edges[1]),"lower"]
		}
	else {
		desc.host1 <- bleph.char[phy$tip.label[desc.edges[1]],"host"]
		desc.region1 <- bleph.char[phy$tip.label[desc.edges[1]],"region"]
		desc.tissue1 <- bleph.char[phy$tip.label[desc.edges[1]],"tissue"]
		}
	if (desc.edges[2] > Ntip(phy)) {
		desc.host2 <- host.MPR[which(rownames(host.MPR)==desc.edges[2]),"lower"]
		desc.region2 <- region.MPR[which(rownames(region.MPR)==desc.edges[2]),"lower"]
		desc.tissue2 <- tissue.MPR[which(rownames(tissue.MPR)==desc.edges[2]),"lower"]
		}
	else { 
		desc.host2 <- bleph.char[phy$tip.label[desc.edges[2]],"host"]
		desc.region2 <- bleph.char[phy$tip.label[desc.edges[2]],"region"]
		desc.tissue2 <- bleph.char[phy$tip.label[desc.edges[2]],"tissue"]
		}

	if((desc.host1 != desc.host2)|(desc.tissue1 != desc.tissue2)) {spec.list.host[node-Ntip(phy),2] <- TRUE} ## use this to count shift in host species OR host tissue as a host shift (probably subject to bias due to polymorphic states)
	##if(desc.host1 != desc.host2) {spec.list.host[node-Ntip(host.simmaps[[tr]]),2] <- TRUE} ## use this to count ounly shift in host species as a host shift
	if(desc.region1 != desc.region2) {spec.list.region[node-Ntip(phy),2] <- TRUE}
}

### getting barplot stats for single tree with actual character traits 
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

MPR.shifts <- c(freq.h.shift.old,freq.h.shift.young,freq.r.shift.old,freq.r.shift.young)

## node shift counts for shuffled data
nsim <- 1000
MPR.summary <- matrix(data=rep(0,nsim),nrow=nsim,ncol=4,dimnames=list(1:nsim,c("h_old","h_young","r_old","r_young")))
names(MPR.shifts) <- c("h_old","h_young","r_old","r_young")


for (i in 1:nsim) {
	root.depth <- max(node.depth.edgelength(phy))
	node.ages <- root.depth - node.depth.edgelength(phy)
	nodes <- (Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
	age <- node.ages[nodes]
	shift <- vector(mode="logical",length=Nnode(phy))
	spec.list.host <- data.frame(age,shift,row.names=nodes)
	spec.list.region <- data.frame(age,shift,row.names=nodes)
	shuffle.order <- c(sample(1:(nrow(bleph.char)-1)),116)
	char.shuffle <- bleph.char
	char.shuffle[,"host"] <- bleph.char[shuffle.order,"host"]
	char.shuffle[,"region"] <- bleph.char[shuffle.order,"region"]
	char.shuffle[,"tissue"] <- bleph.char[shuffle.order,"tissue"]
	host.shuffle.MPR <- MPR(char.shuffle[phy$tip.label,"host"],phy,"Ble141")
	region.shuffle.MPR <- MPR(char.shuffle[phy$tip.label,"region"],phy,"Ble141")
	tissue.shuffle.MPR <- MPR(char.shuffle[phy$tip.label,"tissue"],phy,"Ble141")

	for (node in nodes) {
		desc.edges <- phy$edge[which(phy$edge[,1]==node),2]
		if (desc.edges[1] > Ntip(phy)) {
			desc.host1 <- host.shuffle.MPR[which(rownames(host.shuffle.MPR)==desc.edges[1]),"lower"]
			desc.region1 <- region.shuffle.MPR[which(rownames(region.shuffle.MPR)==desc.edges[1]),"lower"]
			desc.tissue1 <- tissue.shuffle.MPR[which(rownames(tissue.shuffle.MPR)==desc.edges[1]),"lower"]
			}
		else {
			desc.host1 <- char.shuffle[phy$tip.label[desc.edges[1]],"host"]
			desc.region1 <- char.shuffle[phy$tip.label[desc.edges[1]],"region"]
			desc.tissue1 <- char.shuffle[phy$tip.label[desc.edges[1]],"tissue"]
			}
		if (desc.edges[2] > Ntip(phy)) {
			desc.host2 <- host.shuffle.MPR[which(rownames(host.shuffle.MPR)==desc.edges[2]),"lower"]
			desc.region2 <- region.shuffle.MPR[which(rownames(region.shuffle.MPR)==desc.edges[2]),"lower"]
			desc.tissue2 <- tissue.shuffle.MPR[which(rownames(tissue.shuffle.MPR)==desc.edges[2]),"lower"]
			}
		else { 
			desc.host2 <- char.shuffle[phy$tip.label[desc.edges[2]],"host"]
			desc.region2 <- char.shuffle[phy$tip.label[desc.edges[2]],"region"]
			desc.tissue2 <- char.shuffle[phy$tip.label[desc.edges[2]],"tissue"]
			}

		if((desc.host1 != desc.host2)|(desc.tissue1 != desc.tissue2)) {spec.list.host[node-Ntip(phy),2] <- TRUE} ## use this to count shift in host species OR host tissue as a host shift (probably subject to bias due to polymorphic states)
		##if(desc.host1 != desc.host2) {spec.list.host[node-Ntip(host.simmaps[[tr]]),2] <- TRUE} ## use this to count ounly shift in host species as a host shift
		if(desc.region1 != desc.region2) {spec.list.region[node-Ntip(phy),2] <- TRUE}
	}


	### getting barplot stats for single tree 
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

	MPR.summary[i,] <- c(freq.h.shift.old,freq.h.shift.young,freq.r.shift.old,freq.r.shift.young)
	}

diffs <- data.frame(host=(MPR.summary[,"h_old"]-MPR.summary[,"h_young"]),region=(MPR.summary[,"r_old"]-MPR.summary[,"r_young"]))
boxplot(diffs,lab=c(5,5,7),ylim=c(-.7,.4),ylab="difference between old and young split proportions",boxwex=0.5,main="Difference in old and young splits in shift proportions, based \\non MPR reconstruction of actual and shuffled character states")
points(x=c(1,2),y=c(MPR.shifts["h_old"]-MPR.shifts["h_young"],MPR.shifts["r_old"]-MPR.shifts["r_young"]),pch=23,cex=2,col="blue")
points(x=c(1,2),y=c(MPR.shifts["h_old"]-MPR.shifts["h_young"],MPR.shifts["r_old"]-MPR.shifts["r_young"]),pch=19,col="blue")
text(1.5,0.1,labels="actual",col="blue")
text(1.5,-0.3,labels="shuffled")

rm(age,breakpoint,char.shuffle,desc.edges,desc.host1,desc.host2,desc.region1,desc.region2,desc.tissue1,desc.tissue2,diffs,freq.h.shift.old)
rm(freq.h.shift.young,freq.r.shift.old,freq.r.shift.young,h.noshift.old,h.noshift.young,h.shift.old,h.shift.young,host.MPR,host.shuffle.MPR,i)
rm(min_age,MPR.shifts,MPR.summary,node,node.ages,nodes,nsim,phy,r.noshift.old,r.noshift.young,r.shift.old,r.shift.young,region.MPR,region.shuffle.MPR)
rm(root.depth,shift,shuffle.order,spec.list.host,spec.list.region,tissue.MPR,tissue.shuffle.MPR)




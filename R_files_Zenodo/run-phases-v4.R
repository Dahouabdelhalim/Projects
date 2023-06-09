# run phases Ward 
# version 20-1-29
# new:  if a better ram is reached in a phase it will be in the flock of the next phase


# we have the adapted seed
# for each cluster in Ward-15 for each resolution level and for several  phases and we run 16 experiments
# in more details:
# if at least eight different rams are determined the eight best ones are the flock in the next phase
# if there are less than eight different rams we make the flock form mutants of the best ram

# we process memetics for each chosen resolution level until less than eight different rams are produced
# R-levels 1/20, 1/10, 1/5, 1/4, 1/3 (set in adapt-seed)
# if in a phase less than eight different rams occur and the best ram has more than two duplicates we stop


#"R173"

print("seed.nr:", quot=FALSE)
print(seed.nr)

for(rel.R in 1/R.levels)
{
	print("1/rel.R:", quot=FALSE)
	print(1/rel.R)
	go.on <- TRUE
	phase <- 0
# if in a phase we have less than eight different rams 
# we proceed with a flock of mutants of the best ram
	flock.of.mutants <- TRUE
# we collect ram data
	data.results <- data.frame(seed.no=seed.nr, ph=1, mv =mv, rv=rv, f=TRUE, exp=1,  n.links=0, Psi=0, valid=TRUE, nbc=0, dist=0) 
	data.results <- data.results[0,]
	while(go.on)
	{
		print("phase:", quot=FALSE)
		print(phase <- phase + 1)
		parameters <- data.frame(mut.var=mv, ren.var=rv, seed.source="Ls", seed=seed.nr, phase=phase, n.links = m)
# further parameters with:
		source("scripts/parameters-v2.R")
		new.results <- data.frame(seed.no= rep(seed.nr, 16), ph=phase, mv =mv, rv=rv, f=rep(TRUE, 16), exp=1:16,  n.links=0, Psi=0, valid=TRUE, nbc=0, dist=0) 
		for(experiment in 1:16)
		{
			print("experiment:", quot=FALSE)
			print(experiment)
			source("scripts/seed-v3.R")
			print(generation)
			print(nr.rams)
			print(vecPsiL(ram))
			print(vecKinL(ram)/2)
			print(ram.data[nr.rams,1])
			last <- nrow(ram.data)
			new.results$n.links[experiment] <- vecKinL(ram)/2
			new.results$Psi[experiment] <- vecPsiL(ram)
		}
		data.results <- rbind(data.results, new.results)
		save(data.results, file=paste(path, "/data.results-inv-R-", 1/rel.R, "-seed-", seed.nr, "-mut.var-", mut.var, ".RObj", sep=""))
# have we enough different rams for a new flock ?
# if Psi or n.links are different L should be different = not duplicated = nd
		nd <- !duplicated(round(new.results[,7:8], 14))
# if in a phase we have less than eight different rams 
# we proceed with a flock of mutants of the best ram
# if furthermore we have more than three duplicates of the best ram we stop
# data of all 16 results:
		rams <- which(nd)[order(new.results$Psi[nd])]
		if (sum(nd) < pop.size && sum(round(new.results$Psi, 14)==min(round(new.results$Psi, 14))) > 3) go.on <- FALSE
# are there enough different rams:
		if (sum(nd) < pop.size) 
		{
			flock.of.mutants <- TRUE
# genes of best ram of current phase
			load(paste(path, "/final.ram.genes/inv-R-", 1/rel.R, "/final.ram.genes-seed-", seed.nr, "-phase-", phase, "-mut.var-", mut.var, "-ren.var-", ren.var, "-exp-", rams[1], ".RObj", sep="" ))
			best.ram.genes <- ram.genes
		} else {
			flock.of.mutants <- FALSE
			list.flock <- vector("list", length=pop.size)
# genes of eight best rams
			for(i in 1:pop.size)
			{
				load(paste(path, "/final.ram.genes/inv-R-", 1/rel.R, "/final.ram.genes-seed-", seed.nr, "-phase-", phase, "-mut.var-", mut.var, "-ren.var-", ren.var, "-exp-", rams[i], ".RObj", sep="" ))
				list.flock[[i]] <- ram.genes
			}
			
		}
# in seed.R we then make flock of mutants of the best ram 
	}
}

#Archive of analyses for PeerJ submission, Puget Sound Urbanization eDNA, July 2016


##########LOAD LIBRARIES and FUNCTIONS#############
libs=c("vegan","plotrix","MASS","DESeq2","data.table", "gdata","lattice","plyr","dplyr", "lme4", "arm", "gridExtra", "ggplot2", "eeptools", "taxize", "lmtest")
	lapply(libs, require, character.only = TRUE)  #"gsubfn"

bit=function(x){print(x[c(1:10),c(1:10)]); print(dim(x))}  #shows top-left corner of a data.frame/matrix; useful for large datasets

source("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/AnalysisFinal/plot.stacked.2.R")  #for stacked area plots
	#create color palette
	gg_colors <- function(n) {
		  hues = seq(15, 375, length=n+1)
		  hcl(h=hues, l=65, c=100)[1:n]
		}
		
REsim <- function(x, whichel=NULL, nsims){
  require(plyr)
  mysim <- sim(x, n.sims = nsims)
  if(missing(whichel)){
    dat <- plyr::adply(mysim@ranef[[1]], c(2, 3), plyr::each(c(mean, median, sd)))
    warning("Only returning 1st random effect because whichel not specified")
  } else{
    dat <- plyr::adply(mysim@ranef[[whichel]], c(2, 3), plyr::each(c(mean, median, sd)))
  }
  return(dat)
}  #function for random effects simulation

###################################################


#Load data
#########
	load("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Analysis_PeerJ/PS_urban_eDNA_data.rData")
		row.names(presAbs_transects)<-transect_means$Sample.ID
#########		
		
		#checking co-variates; marine variables do not strongly vary with urbanization.
		names(URBAN) # data.frame with environmental covariates is called "URBAN"
		res=NA; for (i in c(7,9:63)){res[i]<-summary(lm(URBAN[,8]~URBAN[,i]))$r.squared}
		res=data.frame(names(URBAN), res); res<-res[order(res[,2]),]
		
		
########################################################################
########\\subsection*{Sequence Data}  #stats for how many reads, etc
########################################################################
#	TaxonSummary=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/AnalysisFinal/TaxonCountsSummaryTable.csv", row.names=1)
		
	sum(DECONTAM) #read counts before rarefaction
		
	sum(master$OTU_abundance)  #total reads passing QC
	sum(master$N_DUPs)  #total number of unique sequences (DUPs) passing QC
	length(master$OTU_name_swarm)  #number of OTUs passing QC
		
	sum(master[!is.na(master$OTU_taxon_phylum),3])  #total reads annotated at some level (i.e., at least phylum)
	sum(master[!is.na(master$OTU_taxon_phylum),1])  #total unique sequences (DUPs) annotated at some level (i.e., at least phylum)
	length(master[!is.na(master$OTU_taxon_phylum),2])  #number of OTUs with some taxonomic assignment

	sum(master[!is.na(master$OTU_taxon_phylum),3])/sum(master$OTU_abundance)  	#fraction of reads with taxonomic assignment
	sum(master[!is.na(master$OTU_taxon_phylum),1])/sum(master$N_DUPs)  #fraction of unique sequences (DUPs) with taxonomic assignment
	length(master[!is.na(master$OTU_taxon_phylum),2])/length(master$OTU_name_swarm)  #fraction of OTUs with some taxonomic assignment
	
	sum(master[!is.na(master$OTU_taxon_family),3])/sum(master$OTU_abundance) #fraction of reads with taxonomy resolved at Family level

	colSums(TaxonSummary[,2:4], na.rm=T) #summary of taxonomy (reads)
		length(unique(master$OTU_taxon_MEGAN))  #total number of unique taxa 

	table(master$OTU_taxon_phylum)  #unique OTUs by phylum
		table(master$OTU_taxon_phylum)/sum(table(master$OTU_taxon_phylum))  #OTUs by phylum, proportions of annotated reads

	sum(master[master$eVal<1e-30,3], na.rm=T)/sum(master$OTU_abundance) #high-confidence annotations (reads)
	length(master[master$eVal<1e-30,3])/nrow(master)  #high-confidence annotations (OTUs, as fraction of all OTUs)
	length(master[master$eVal<1e-30,3])/sum(!is.na(master$OTU_taxon_MEGAN))  #high-confidence annotations (OTUs, as fraction of all annotated OTUs)

		#table of blast hits by evalue.  maybe a set of mini-histograms (eval freq) for OTUs by phylum?
		uniquePhyla=unique(master$OTU_taxon_phylum) ; uniquePhyla= uniquePhyla[!is.na(uniquePhyla)]
		eVal_summary<-list(NA)
		for (i in 1:length(uniquePhyla)){
			eVal_summary[[i]]<- table(master[master$OTU_taxon_phylum== uniquePhyla[i],9])
		}
				q=matrix(NA, nrow= length(uniquePhyla), ncol=11) ; colnames(q)=names(eVal_summary[[1]])
					for (i in 1:10){
						q[i,]=eVal_summary[[i]][match(colnames(q), names(eVal_summary[[i]]))]
					}
					q[is.na(q)]<-0 ; row.names(q)=uniquePhyla; q=q[,order(colnames(q), decreasing=T)]
						#write.csv(q, "/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/OTU_eVals_byPhylum.csv")

		#generate table of all family-level annotations
			a=aggregate(master[,c(3)], by=list(master$OTU_taxon_family), length) ; names(a)=c("Taxon", "OTUs")
			b=aggregate(master[,c(3)], by=list(master$OTU_taxon_family), sum) ; names(b)=c("Taxon", "Reads")
			familySummary= merge(a, b) ; familySummary=familySummary[order(familySummary[,3], decreasing=T),]

		#add natural history to family data
		Fams_w_NatHist=merge(familySummary, eDNANatHist, by="Taxon", all=F)
		Fams_w_NatHist= Fams_w_NatHist[,-c(4,5)]; names(Fams_w_NatHist)[2:3]<-c("OTUs", "Reads")
	#		write.csv(Fams_w_NatHist, "/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/eDNAFamilies_w_NatHist.csv",row.names=F)

		#are more common OTUs more likely to be annotated?
		isAnnotated=rep(1, times=nrow(master))
		isAnnotated[which(is.na(master$OTU_taxon_MEGAN))]<-0
		
		wilcox.test(master[isAnnotated==1,3], master[isAnnotated==0,3])




		#is there simply more DNA in the water at more urbanized sites?  #no
			a=merge(meta, completeMetadata) ; a=a[a$PSeelgrassSite=="Y",] ; names(a)[9]="Site.code"
				b=data.frame(URBAN$Site.code, URBAN$Area.Weighted.Mean.Percent.Imperviousness) ; names(b)=c("Site.code", "Imperv")
				c=merge(a, b, "Site.code")
			summary(lm(mass_DNA_total_ng~ Imperv, data=c))




########################################################################
##Alpha diversity mean number of OTUs/families per site, using transects as replicates
########################################################################
	DNA_site_MeanRichness = aggregate(specnumber(presAbs_transects), by=list(substr(row.names(presAbs_transects),1,3)), mean)[,2]
	mean_eDNA_families=aggregate(specnumber(families), list(substr(row.names(presAbs_transects), 1,3)), mean)[,2]


#also calculate total richness at each site, in case it's useful later
	DNA_site_TotalRichness= aggregate(presAbs_transects, by=list(substr(row.names(presAbs_transects),1,3)), sum)
		DNA_site_TotalRichness= DNA_site_TotalRichness[,-1] ; DNA_site_TotalRichness<-specnumber(DNA_site_TotalRichness)
		names(DNA_site_TotalRichness)<-unique(substr(row.names(presAbs_transects),1,3))
	DNA_TotalFamilies=specnumber(aggregate(families, list(substr(row.names(families), 1,3)), sum))





########################################################################
####BETA Diversity
#note beta diversity means many things to many people (see Anderson 2011 Ecol Letters, Koleff et al. 2003 J. Animal Ecol).  Here, looking at how different transects are within sites
########################################################################
	betadiv= temp[,dna_index] ; betadiv[betadiv>0]<-1 #create pres-abs df for beta diversity analysis	
	distance=NA; betaList=list(NA)  #create empty objects for results .  betaList will keep all pairwise interactions ; "distance" will use a parallel calculation for the same idea
		index=1
		for (i in unique(temp$Site.code)){
			site.df=betadiv[temp$Site.code%in%i,]
				betaList[[index]]<-as.vector(betadiver(site.df, method=1))
				distance[index]=mean(as.vector(vegdist(site.df, method="jaccard")))
				index=index+1
			}

	beta.df<-as.data.frame(betaList) ; names(beta.df)<-unique(temp$Site.code) ;beta.df[2:3,3]<-NA  #save as data frame; replace repeated values as NA to reflect fact that a two-sample comparison has only a single beta value (that sample didn't have all three transects amplify)
	
	beta<-colMeans(beta.df, na.rm=T)  #keep means as "beta"

	#plot(beta~unique(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness), pch=19, col='black', ylim=c(0.4,.95))
	#for (i in 1:length(unique(temp$Site.code))){
	#	points(beta.df[,i]~rep(unique(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness)[i], times=length(!is.na(beta.df[,i]))))
	#}

	##between sites within categories; could use mean of all cross-site pairwise transect comparisons.  here, using site-wide richness to compare one site to another.
	siteOTUvectors<-aggregate(presAbs_transects, by=list(substr(row.names(presAbs_transects),1,3)), sum) ; row.names(siteOTUvectors)<-siteOTUvectors[,1]
		siteOTUvectors<-siteOTUvectors[,-1] ; siteOTUvectors[siteOTUvectors>0]<-1 #convert to PresAbs
			urbanCats<-cofactor.df$category.x[match(row.names(siteOTUvectors), cofactor.df$Site.code)]
 		moreUrbanSiteOTUs<-siteOTUvectors[urbanCats=="more urban",]
 		lessUrbanSiteOTUs<-siteOTUvectors[urbanCats=="less urban",] 			
				#calculate betas
					MoreUrbanSiteBetas<-as.vector(vegdist(moreUrbanSiteOTUs, method="jaccard"))
					LessUrbanSiteBetas<-as.vector(vegdist(lessUrbanSiteOTUs, method="jaccard"))
					AllSiteBetas<-as.vector(vegdist(siteOTUvectors, method="jaccard"))

				BetweenSiteBetas<-data.frame(c(MoreUrbanSiteBetas, LessUrbanSiteBetas, AllSiteBetas),c(rep("More Urban", length(MoreUrbanSiteBetas)), rep("Less Urban", length(LessUrbanSiteBetas)), rep("All Sites", length(AllSiteBetas)))) ; names(BetweenSiteBetas)<-c("WhittakersBeta", "SiteCategory")
				boxplot(WhittakersBeta~SiteCategory, data=BetweenSiteBetas)

				wilcox.test(MoreUrbanSiteBetas, LessUrbanSiteBetas) ; summary(MoreUrbanSiteBetas);  summary(LessUrbanSiteBetas)


	
	#RaupCrick accounts for the fact that an increase in richness can, itself, explain a change in beta diversity.  This metric estimates beta diversity by permuting the underlying matrix of occurrences of taxa across sites and across taxa, and using this permuted matrix as a null model.  In effect, it asks, "given the arrangement of alpha diversity that we see, what beta diversity would we see by chance alone?"  Here, I did this to make sure that our beta diversity trend is real.  It seems that the trend is real, because beta diversity still declines with imperviousness, independent of changes in richness.
	raupCrick_dna=NA
	raupCrick_matrix= as.matrix(raupcrick(temp[,dna_index], null="swap"))
		for (i in 1:length(unique(temp$Site.code))){
			colindex=which(substr(colnames(raupCrick_matrix),1,3)==unique(temp$Site.code)[i])
				site.df= raupCrick_matrix[colindex, colindex]
			raupCrick_dna[i]=mean(site.df[lower.tri(site.df)])
		}

	plot(raupCrick_dna~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #raup-Crick essentially maintains beta diversity relationship, although this depends upon the null model used. Models r0 and r1 generate nonsensical results.  Here, using "swap", the algorithm permutes the original binary data matrix while leaving row and col sums intact, preserving the probabilities of picking a given taxon at random.

########################################################################
##GAMMA diversity
########################################################################
			gamma.df=data.frame(DNA_site_TotalRichness,unique(cofactor.df[,c(1,7)])[,2]) ; names(gamma.df)<-c("GammaDiversity","UrbanCategory") ; gamma.df=data.frame(gamma.df, unique(cofactor.df[,c(1,6,1737)]))
				#plot(GammaDiversity~ UrbanCategory, data=gamma.df, col=gg_colors(2)[as.numeric(as.factor(gamma.df$UrbanCategory))], ylab="Gamma Diversity (OTUs/site)", xlab="Urban Category")  #if we're looking at total OTUs per site, rather than regional diversity

	#how many OTUs/Families total in more/less urban sites?
	moreurban=presAbs_transects[cofactor.df$category.x=="more urban",]
	lessurban=presAbs_transects[cofactor.df$category.x=="less urban",]
		sum(colSums(lessurban)>0) ; sum(colSums(moreurban)>0)  #total OTUs in each category of sites
		sum(colSums(families[cofactor.df$category.x=="more urban",])>0) ; sum(colSums(families[cofactor.df$category.x=="less urban",])>0)


########################################################################
#########site pairs
########################################################################


	OTUrichnessvec=NA  #here, counting richness of the site as the mean number of unique OTUs that occur transects in that site (i.e., alpha diversity)
	Manualrichnessvec=NA  #here, counting richness of the site as the number of unique OTUs that occur at any transect in that site (i.e., overall occurrence)
	sitenamevec=NA
		for (i in unique(substr(row.names(temp),1,3))){
			site.df= temp[substr(row.names(temp), 1, 3)%in%i,]
				Manualrichnessvec <-c(Manualrichnessvec, specnumber(colSums(site.df[,manual_index])))
				OTUrichnessvec <-c(OTUrichnessvec, mean(specnumber(site.df[,dna_index])))
			sitenamevec<-c(sitenamevec, i)
			}
	Manualrichnessvec = Manualrichnessvec[-1] ; sitenamevec = sitenamevec[-1] ; OTUrichnessvec = OTUrichnessvec[-1]
		deltaRichness.df=data.frame(sitenamevec, OTUrichnessvec, Manualrichnessvec) ; names(deltaRichness.df)=c("Site.code", "OTUrichness", "Manualrichness")
			urbanizationdata=envir[,c(3,6,7,10, 13, 14)]
			deltaRichness.df = merge(deltaRichness.df, urbanizationdata) ; deltaRichness.df = merge(deltaRichness.df, URBAN[,c(2,13,29,31)]) ; deltaRichness.df =unique(deltaRichness.df)


########################################################################
#PLOTTING DIVERSITY
####Combined Figure2 Plot. alpha, beta, gamma diversity
########################################################################

	#necessary calcs for the following figure
		mean_eDNA_families=aggregate(specnumber(families), list(substr(row.names(families), 1,3)), mean)[,2]
		deltaRichness.df=data.frame(deltaRichness.df, unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)); names(deltaRichness.df)[12]="Imperv"
		betaSitePairs=data.frame(unique(cofactor.df[,c(1,6,1737)]), beta)
#			gamma.df=data.frame(DNA_site_TotalRichness,unique(cofactor.df[,c(1,7)])[,2]) ; names(gamma.df)<-c("GammaDiversity","UrbanCategory") ; gamma.df=data.frame(gamma.df, unique(cofactor.df[,c(1,6,1737)]))



def.par <- par(no.readonly = TRUE) # save default, for resetting...


nf <- layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), widths = rep.int(1.5, 2), respect = TRUE)
#layout.show(nf)

	par(mgp=c(2,.5,0)) #distance between axis titles and the axes themselves

#SUB1 #eDNA OTU richness and Family Richness
	
	par(mar=c(1,4,1,4))
	plot(specnumber(cofactor.df[,dna_index])~ cofactor.df $Area.Weighted.Mean.Percent.Imperviousness, pch=19, ylab="Alpha Diversity (OTUs/Transect)", xlab="", col=gg_colors(2)[2], cex=0.5, xlim=c(1,45), axes=F) ; axis(1, labels=F) ; axis(2, las=2) ; box(bty="o")
		points(unique(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness),DNA_site_MeanRichness, pch=19, col=gg_colors(2)[2], cex=1.2)
		abline(lm(DNA_site_MeanRichness ~ unique(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness)), col=gg_colors(1), lwd=2)

		#scaling factor to plot N families on same plot : min(specnumber(cofactor.df[,dna_index]))/min(specnumber(families))  #5.3
		points(specnumber(families)*5.3~ I(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness+2), col='lightgrey', pch=19, cex=0.5)  #eDNA
			points(I(unique(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness)+2), mean_eDNA_families*5.3, pch=19, col='lightgrey', cex=1.2)
			abline(lm(mean_eDNA_families*5.3~I(unique(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness)+2)), lty=2, col=gg_colors(1))
		axis(4, seq(50,350,50), labels=round(seq(50,350,50)/5.3))
		mtext("Taxonomic Family Richness", side=4, line=2, cex=0.7)
		legend("topleft", c("OTUs", "Taxonomic Families"), fill=c(gg_colors(2)[2], "lightgrey"), cex=.6)

#SUB2  #richness dot-and-line
#NOTE: needs original data points added (not just means)
	par(mar=c(1,4,1,4))
			plot(deltaRichness.df $OTUrichness ~ deltaRichness.df$Imperv, pch=19, type="n", ylim=c(0.9*min(OTUrichnessvec),1.1*max(OTUrichnessvec)), xlab="Watershed Imperviousness (Percent)", ylab="Mean Alpha Diversity (OTUs/Transect)", xlim=c(1,45), axes=F) ; axis(1, labels=F) ; axis(2, las=2) ; box(bty="o")
								#abline(lm(deltaRichness.df$OTUrichness~deltaRichness.df$Imperv), col="black", lty=2)
			for(i in 1:length(unique(deltaRichness.df$pair.number))){
			site.df= deltaRichness.df[deltaRichness.df$pair.number==unique(deltaRichness.df$pair.number)[i],]
				points(site.df$OTUrichness ~site.df$Imperv, pch=19, col=gg_colors(4)[i])
					lines(site.df$Imperv,site.df$OTUrichness, col=gg_colors(4)[i])
					#text(site.df$PC1.74,site.df$OTUrichness +0.5, i)
			}
				legend("topleft", c("RB-SM","BG-CW","PC-CC","SI-MA"), fill=gg_colors(4), cex=0.6)


#SUB3  Beta Div by Imperviousness
	par(mar=c(3,4,1,4))
	g=data.frame(betaMeans.df[,1:1000],unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)) ; names(g)=c(paste0('sample', 1:1000), "Imperv")
	g=reshape(g, varying=list(1:1000), direction="long")
	#it is a huge pain to do a boxplot with an x-axis proportional to the numeric values of the grouping factors...
	xvals=round(g$Imperv,1)

	boxplot(g$sample1 ~ 
		factor(xvals, levels = seq(min(xvals), max(xvals),0.1)),
		data = g, 
		xaxt = 'n',
		boxwex=10,
		cex=.5,
		border=gg_colors(2)[2],
		xlab="Watershed Imperviousness (Percent)",
		ylab="Transect-Level Beta Diversity"
		) 
	axis(1, at = c(-5,84,184,284,384), labels=c(0,10,20,30,40)) 
		lin.mod=lm(sort(rowMeans(betaMeans.df), decreasing=T)~c(1,8,25,56,238,321,336,406)) #note I needed to mess with x-values for plotting here
	abline(lin.mod, col=gg_colors(1), lwd=2)


#SUB4	Pairs, beta dot-and-line (needs original data points added, not just means)
	par(mar=c(3,4,1,4))
			plot(betaSitePairs $beta ~ betaSitePairs $Area.Weighted.Mean.Percent.Imperviousness, pch=19, type="n", ylim=c(0.9*min(betaSitePairs$beta),1.1*max(betaSitePairs$beta)), xlab="Watershed Imperviousness (Percent)", ylab="Mean Transect-Level Beta Diversity", xlim=c(1,45), axes=F) ; axis(1, labels=T) ; axis(2, las=2) ; box(bty="o")
								#abline(lm(betaSitePairs $beta~ betaSitePairs $Area.Weighted.Mean.Percent.Imperviousness), col="black", lty=2)
			for(i in 1:length(unique(betaSitePairs $pair.number.x))){
			site.df= betaSitePairs[betaSitePairs $pair.number.x==unique(betaSitePairs $pair.number.x)[i],]
				points(site.df$beta ~site.df$Area.Weighted.Mean.Percent.Imperviousness, pch=19, col=gg_colors(4)[i])
					lines(site.df$Area.Weighted.Mean.Percent.Imperviousness,site.df$beta, col=gg_colors(4)[i])
					#text(site.df$PC1.74,site.df$OTUrichness +0.5, i)
			}
				legend("topright", c("RB-SM","BG-CW","PC-CC","SI-MA"), fill=gg_colors(4), cex=0.6)





#SUB5-6 ; gamma accumulation; depends on /Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Analysis_UrbanizationSecondSubmission/GammaDiversity_accumulationCurve.R
	boxplot(gammaResultsAll ~ c(rep(1:8, rep(NsamplesPerSite,8))), ylim=c(100,1800), border =gg_colors(3)[1], xlab="N Sites", ylab="Unique OTUs")
			ydata= gammaResultsAll; xdata= c(rep(1:8, rep(NsamplesPerSite,8)))
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col=gg_colors(3)[1], lwd=2)
	
	boxplot(gammaResultsMoreUrban ~ c(rep(1:4, rep(NsamplesPerSite,4))), add=T, border =gg_colors(3)[2])
		ydata= gammaResultsMoreUrban; xdata= c(rep(1:4, rep(NsamplesPerSite,4)))
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,4), add=T, col=gg_colors(3)[2], lwd=2)

	boxplot(gammaResultsLessUrban ~ c(rep(1:4, rep(NsamplesPerSite,4))), add=T, border =gg_colors(3)[3])
		ydata= gammaResultsLessUrban; xdata= c(rep(1:4, rep(NsamplesPerSite,4)))
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,4), add=T, col=gg_colors(3)[3], lwd=2)

	legend("bottomright",legend=c("All Sites", "More Urban Sites", "Less Urban Sites"), fill=gg_colors(3))


par(def.par)  #- reset to default



########################################################################
##Regressions with Imperviousness
########################################################################
###Imperviousness Correlations w OTU richness and Family Richness
	summary(lm(DNA_site_MeanRichness~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))) #main correlation, using site means.
	summary(glm(round(mean_eDNA_families)~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness), family='poisson')) #number of taxonomic Families w imperviousness (using Poisson with site means, because of small numbers of counts)


#rare vs. common OTUs
#rare OTUs actually increase in relative frequency at high imperviousness
		common=temp[,dna_index][,colSums(temp[,dna_index])>100]
		rare=temp[,dna_index][,colSums(temp[,dna_index])<11]
			summary(lm(I(specnumber(rare)/specnumber(common))~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))


		
########################################################################
###########mixed-effects models
########################################################################
	#Set up dataframe with OTU richness, Imperv surface cover, covariates, and Family richness
		rareRichness = specnumber(temp[,dna_index], MARGIN=1)  #note: transect-level
		Imperv=cofactor.df$Area.Weighted.Mean.Percent.Imperviousness
			rareRichness=data.frame(names(rareRichness), rareRichness) ; names(rareRichness)[1]="Sample.ID"
	LMEtemp=merge(temp[,1:7], rareRichness, "Sample.ID", all.x=F)  #add to analysis df
		LMEtemp =data.frame(LMEtemp, specnumber(families)) ; names(LMEtemp)[9]="FamilyRichness"
						xyplot(rareRichness ~ category | pair.number, data= LMEtemp)  #vizualize the problem
						xyplot(FamilyRichness ~ category | pair.number, data= LMEtemp)  #vizualize the problem
					

	#LME for richness by imperviousness, 
		#syntax:  responseVar~fixed+(random) ; the number 1 indicates the intercept as a variable, such that random var (1|pair.number) means that each pair has its own intercept (but the same slope), and random var (Imperv|pair.number) estimates a separate error term the Imperv slope for each pair.
		
		
		(RandInterceptFixedSlope=lmer(rareRichness ~ Imperv + (1|pair.number) + (1 | Site.code), data= LMEtemp, REML=T))
	
		ranef(RandInterceptFixedSlope, condVar=TRUE, whichel = "pair.number")  #show the random effects
		fixef(RandInterceptFixedSlope, condVar=TRUE, whichel = "pair.number")  #show the fixed effect
		REsim(RandInterceptFixedSlope, whichel = "pair.number", nsims = 1000)  #simulate the random effects

		1-var(residuals(RandInterceptFixedSlope))/(var(model.response(model.frame(RandInterceptFixedSlope))))  #like an r-squared; see http://stats.stackexchange.com/questions/95054/how-to-get-the-overall-effect-for-linear-mixed-model-in-lme4-in-r


########################################################################
### INDIVIDUAL OTUs, pres/abs correlations with Imperviousness  [LOGISTIC REGRESSION]  ; see /Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Analysis/DiversityAnalysis/DESeq2_working.R for deseq analysis of abundances, if interested
########################################################################
		namevec=NA ; pvalvec=NA ; logitvec=NA ; AICvec=NA ; ORvec=NA
			for (i in 1:ncol(presAbs_transects)){
				namevec[i]=names(presAbs_transects)[i]
					mod=glm(presAbs_transects[,i]~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, family = "binomial"(link='logit'))
				pvalvec[i]=summary(mod)$coefficients[8]
				logitvec[i]=summary(mod)$coefficients[2]
				AICvec[i]=AIC(mod)
				ORvec[i]=exp(summary(mod)$coefficients[2])
			}		
		results=data.frame(namevec, pvalvec, logitvec, AICvec, ORvec)  ;  results=results[order(results[,5], decreasing=T),]
		
		#OTUs most positively correlated with Imperv :
		increasing=subset(results, results$logitvec>0&results$pvalvec<0.05)  ; nrow(increasing)
		a= masterDECONTAM[match(increasing$namevec, masterDECONTAM $OTU_name_swarm),]  #bivalves and echinoids
		table(a$OTU_taxon_family) ; table(a$OTU_taxon_class)

				#null hypothesis: are these taxa overrepresented, relative to their proportion in the overall dataset?
				table(master$OTU_taxon_family)
				#e.g., Mytilus has 14% of OTUs in cleaned database:  217/sum(table(littleMaster$OTU_taxon_family))

		#OTUs most negatively correlated with Imperv :
		decreasing=subset(results, results$logitvec<0&results$pvalvec<0.05)
		b= masterDECONTAM[match(decreasing$namevec, masterDECONTAM $OTU_name_swarm),]  #barnacles, but nearly nothing
		table(b$OTU_taxon_family)



########################################################################
###number of families per class or order
########################################################################
		fams_temp=merge(cofactor.df, families1, "Sample.ID", all.x=F)  #use 'phyla1', 'classes1', etc... from above ; This shows how many unique OTUs in Mytilidae, e.g., occurred in each transect
			row.names(fams_temp)= cofactor.df $Sample.ID
		fams2classes=lineages[match(names(fams_temp[,1793:ncol(fams_temp)]), lineages$OTU_taxon_family),7]		
		fams2orders=lineages[match(names(fams_temp[,1793:ncol(fams_temp)]), lineages$OTU_taxon_family),6]

		#Example analysis:
		clams=subset(masterDECONTAM, masterDECONTAM$OTU_taxon_class=="Bivalvia") 

		df=data.frame(cofactor.df $Area.Weighted.Mean.Percent.Imperviousness,specnumber(cofactor.df[,(names(cofactor.df)%in%clams$OTU_name_swarm)])) ; names(df)=c("Imperv","BivalveOTUrichness")
			mean(df[df[,1]<10, 2]) #mean number bivalve OTUs at imperv less than 10
			mean(df[df[,1]>25, 2]) #mean number bivalve OTUs at imperv greater than 25
			summary(glm(BivalveOTUrichness ~Imperv, family="poisson", data=df))
	
		df=data.frame(specnumber(fams_temp[, which(fams2classes=="Bivalvia")+1792]), cofactor.df$Area.Weighted.Mean.Percent.Imperviousness) ; names(df)=c("Richness", "Imperv")
			mean(df[df[,2]<10, 1]) #mean number bivalve families at imperv less than 10
			mean(df[df[,2]>25, 1]) #mean number bivalve families at imperv greater than 25
			summary(glm(Richness ~Imperv, family="poisson", data=df))



########################################################################
########\\subsection*{eDNA Variance Within and Among Sites}
########################################################################
		#PERMANOVA   #using presence/absence version of normalized counts, with no transect-level averaging
			distmat=vegdist(presAbs, method="jaccard") #bray and jaccard give very similar answers.  euclidian gives more variance to the residuals.
				transects=substr(row.names(presAbs), 1,5)
				sites=substr(row.names(presAbs), 1,3)
				adonis(distmat~sites+transects)
				
				#with count data rather than presence/Absence
				distmat=vegdist(t(rareData), method="bray") #bray and jaccard give very similar answers.  euclidian gives more variance to the residuals.
				transects=substr(row.names(t(rareData)), 1,5)
				sites=substr(row.names(t(rareData)), 1,3)
				adonis(distmat~sites+transects)

########################################################################
#constrained correspondence analysis; NOTE: could constrain by site pair number or by urban category; both very strong results
########################################################################
		#transect-level ordination
		dna.cca=cca(presAbs_transects~cofactor.df$category.x)
			ordiplot(dna.cca, display="sites")
			samplingsites=substr(row.names(presAbs_transects),1,3)
				ordihull(dna.cca,groups= samplingsites,draw="polygon",col="grey90",label=F)
				orditorp(dna.cca,display="sites")

		dna.pca=prcomp(presAbs_transects)
#			pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/SecondSubmission/SupplFig_DNA_PCA.pdf", width=7, height=7)
			ordiplot(dna.pca, display="sites", main="PCA on eDNA OTU Data\\nPresence-Absence")
			samplingsites=substr(row.names(presAbs_transects),1,3)
				# ordihull(dna.pca,groups= samplingsites,draw="polygon",col="grey90",label=F)
				# orditorp(dna.pca,display="sites")
			ordihull(dna.pca,groups= samplingsites,draw="polygon",col= gg_colors(2)[1],label=F, show.groups=unique(samplingsites)[c(1,5,6,7)])
			ordihull(dna.pca,groups= samplingsites,draw="polygon",col= gg_colors(2)[2],label=F, show.groups=unique(samplingsites)[c(2:4,8)])
				orditorp(dna.pca,display="sites")
										legend("topright", c("More Urban", "Less Urban"), fill=gg_colors(2), cex=0.8)
#			dev.off()

		###Site-level ordination

		a=aggregate(presAbs_transects, by=list(substr(row.names(presAbs_transects),1,3)), sum) ; a=a[,-1]
			a[a>0]<-1
		dna.cca=cca(a~unique(cofactor.df[,c(1,7)])[,2])



########################################################################
####Functional Diversity Analysis
########################################################################
natHist_index=c(5:14)
	functionalDiv_Families=matrix(NA, nrow(families), length(natHist_index)); index=1
		for (i in unique(row.names(families))){
			transect.df=families[i, families[i,]>0]
			functionalDiv_Families[index,]<-colSums(eDNANatHist[match(names(transect.df), eDNANatHist[,1]), natHist_index])  #just counts each family once, no matter how many OTUs present at a transect
				index=index+1
			}
		functionalDiv_Families =as.data.frame(functionalDiv_Families) ; row.names(functionalDiv_Families)<-row.names(families) ; colnames(functionalDiv_Families)<-names(eDNANatHist)[natHist_index]    #Note that each family can/will fit into more than one functional category (e.g., epifauna and mobile)


##COUNT data to answer question, how common are non-marine species in the data? 

sum(eDNANatHist[eDNANatHist$Habitat.terrestrial==1| eDNANatHist$Habitat.freshwater==1,3]) #total freshwater and terrestrial reads, including human; use [2:37] to exclude human
sum(eDNANatHist[eDNANatHist$Habitat.intertidal==1| eDNANatHist$Habitat.subtidal==1,3]) #total marine reads

	#Plotting functional groups
	plot(rowSums(functionalDiv_Families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)  #crude measure of functional diversity; summing rows
	
		par(mfrow=c(5,2))
		par(mar=c(1,2,1,1))
		par(mgp=c(1,.5,0))
		for (i in 1:ncol(functionalDiv_Families)){
		plot(functionalDiv_Families[,i]~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, main="", pch=19, ylab="") 
			title(names(functionalDiv_Families)[i], line=-7)
		mod=lm(functionalDiv_Families[,i]~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness) #; if(summary(mod)$coefficients[8]<0.05) abline(mod, col='red') else abline(mod, col='lightgrey')
		
		
		ydata=functionalDiv_Families[,i]
		xdata=cofactor.df$Area.Weighted.Mean.Percent.Imperviousness
		nonlin.mod=nls(ydata~a*(xdata^2)+b*xdata+c, start=list(a=1, b=1, c=1))
			#lines(sort(xdata), predict(nonlin.mod)[order(xdata)])
		
		if(lrtest(mod, nonlin.mod)$Pr[2]<0.05) lines(sort(xdata), predict(nonlin.mod)[order(xdata)]) else if(summary(mod)$coefficients[8]<0.05) abline(mod, col='red') else abline(mod, col='lightgrey')
		}


###Assigning Niche Categories
		NatHistAttributes=NA; for(i in 1:10){NatHistAttributes[i]<-strsplit(names(eDNANatHist)[natHist_index], "\\\\.")[[i]][2]}
		
		TheoreticalFunctionalNiches=expand.grid(NatHistAttributes[1:4], NatHistAttributes[5:8], NatHistAttributes[9:10])
			TheoreticalFunctionalNiches <- TheoreticalFunctionalNiches[-c(3,4,19,20,24,28,32),] #remove logically impossible niches
				row.names(TheoreticalFunctionalNiches)<-1:nrow(TheoreticalFunctionalNiches)
		
		bestmatch=NA
		listMatches=list(NA)
		for (i in 1:nrow(eDNANatHist)){
			test=NatHistAttributes[eDNANatHist[i, natHist_index]==1]
			for (j in 1:nrow(TheoreticalFunctionalNiches)){bestmatch[j]<-sum(test%in%unlist(TheoreticalFunctionalNiches[j,]))}	
			listMatches[[i]]<-TheoreticalFunctionalNiches[bestmatch==3, 1:3]
		}
		
		niche=list(NA)
		for (i in 1:nrow(eDNANatHist)){
			for (j in 1:nrow(listMatches[[i]])){
				if(nrow(listMatches[[i]])==1) niche[[i]]<-paste(as.character(unlist(listMatches[[i]])), collapse="_")
				if(nrow(listMatches[[i]])>1) if(j==1) niche[[i]]<-paste(as.character(unlist(listMatches[[i]][j,])), collapse="_")
					if(j>1) niche[[i]][j]<-paste(as.character(unlist(listMatches[[i]][j,])), collapse="_")
			}
		}
		
		observedNicheNames=unique(unlist(niche))
		observedNiches=as.data.frame(matrix(0, nrow=nrow(eDNANatHist), ncol=length(observedNicheNames)))
		for (i in 1:nrow(eDNANatHist)){
			observedNiches[i,match(niche[[i]],observedNicheNames)]<-1
		}
		row.names(observedNiches)<-eDNANatHist$Taxon ; colnames(observedNiches)<-observedNicheNames

	#assign taxonomic families to observed niches and map back to transect-level data
	Niche_Families=as.data.frame(matrix(NA, nrow(families), ncol(observedNiches))); index=1
		for (i in unique(row.names(families))){
			transect.df=families[i, families[i,]>0]
			Niche_Families[index,]<-colSums(observedNiches[match(names(transect.df), row.names(observedNiches)),])  #just counts each family once, no matter how many OTUs present at a transect
				index=index+1
			}
			
		row.names(Niche_Families)<-row.names(families) ; colnames(Niche_Families)<-names(observedNiches)    #Note that each family can be in more than one niche, because of variation in taxa below the family level.


	#Alpha functional div -- number of unique niches per transect. 

	#number of functional groups goes up with imperviousness...
plot(specnumber(Niche_Families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)  
	summary(lm(specnumber(Niche_Families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))

	#... but this is a function of the overall increase in OTU and taxon richness. Normalizing by OTU or taxon richness reveals a decline in functional diversity
plot(specnumber(Niche_Families)/specnumber(temp[,dna_index])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)  #N niches per OTU
	summary(lm(specnumber(Niche_Families)/specnumber(temp[,dna_index])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #strong decrease with imperviousness
plot(specnumber(Niche_Families)/specnumber(families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)  #N niches per family
	summary(lm(specnumber(Niche_Families)/specnumber(families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #strong decrease with imperviousness


		#dealing with pseudo-replication by using site means in calculations
		NicheFamiliesMeans<-aggregate(specnumber(Niche_Families), list(substr(names(specnumber(Niche_Families)),1,3)), mean)
		NicheFamiliesSums<-specnumber(aggregate(Niche_Families, list(substr(names(specnumber(Niche_Families)),1,3)), mean))
		
		FamilyNumberMeans<-	aggregate(specnumber(families), list(substr(names(specnumber(families)),1,3)), mean)
		FamilyNumberSums<-	specnumber(aggregate(families, list(substr(names(specnumber(families)),1,3)), sum))

			# plot(I(NicheFamiliesMeans[,2]/FamilyNumberMeans[,2])~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
				summary(lm(I(NicheFamiliesMeans[,2]/FamilyNumberMeans[,2])~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)))
				summary(lm(NicheFamiliesMeans[,2]~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)))


			##PLOTTING
#			pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/SecondSubmission/Niches_by_Imperviousness1x2.pdf", width=8, height=3.3)
			par(mfrow=c(1,2))
			par(mar=c(4,5,1,2))
			par(mgp=c(2,1,0))
				plot(NicheFamiliesMeans[,2]~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness), pch=19, ylab="Life-History Richness \\n(Life-Histories/Site)", xlab="Watershed Imperviousness (Percent)", ylim=c(10,18), col=gg_colors(2)[2])
					abline(lm(NicheFamiliesMeans[,2]~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)), col=gg_colors(1), lwd=2)
					points(specnumber(Niche_Families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, cex=0.5, pch=19, col=gg_colors(2)[2])
					
				
				plot(I(NicheFamiliesMeans[,2]/FamilyNumberMeans[,2])~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness), pch=19, ylab="Per-Taxon Life-History Richness \\n(Life-Histories/Site/Family)", xlab="Watershed Imperviousness (Percent)", ylim=c(0.2,1), col=gg_colors(2)[2])
					abline(lm(I(NicheFamiliesMeans[,2]/FamilyNumberMeans[,2])~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)), col='grey', lwd=2)
						meanFamsPerSite=NA
						for(i in 1:8){meanFamsPerSite=c(meanFamsPerSite, rep(FamilyNumberMeans[i,2], length(grep(NicheFamiliesMeans[i,1], row.names(Niche_Families)))))} ; meanFamsPerSite= meanFamsPerSite[-1]
					points(I(specnumber(Niche_Families)/meanFamsPerSite)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, cex=0.5, pch=19, col=gg_colors(2)[2])
#			dev.off()



				#splitting into more/less urban categories for wilcox.test or t.test, etc
				more1<-NicheFamiliesMeans[unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)>10, 2] #moreUrban
				less1<-NicheFamiliesMeans[unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)<10, 2] #lessUrban
					wilcox.test(more1, less1)

				more2<- (NicheFamiliesMeans[unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)>10, 2])/(FamilyNumberMeans[unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)>10, 2])
				less2<- (NicheFamiliesMeans[unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)<10, 2])/(FamilyNumberMeans[unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)<10, 2])




	Niche_Families <- Niche_Families[,order(colSums(Niche_Families), decreasing=T)]  #order cols by frequency of niche	
	orderedNicheFamilies <- Niche_Families[order(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, decreasing=F),]  #order rows by imperviousness, for plotting
	#barplot(t(as.matrix(orderedNicheFamilies[,7:10])), las=2, col=gg_colors(6), legend=T, args.legend=list(x="topleft"), ylab="N Families")
	#barplot(t(as.matrix(orderedNicheFamilies/rowSums(orderedNicheFamilies))), las=2, col=gg_colors(6), legend=F, args.legend=list(x="topleft"), ylab="N Families")

	#plotting these niches by imperviousness, after normalizing by N families present at the site. 
		par(mfrow=c(5,2))
		par(mar=c(1,2,1,1))
		par(mgp=c(1,.5,0))
		for (i in 10:19){
		plot(Niche_Families[,i]/specnumber(families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, main="", pch=19, ylab="") 
			title(names(Niche_Families)[i], line=-7)
		mod=lm(Niche_Families[,i]/specnumber(families)~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness) #; if(summary(mod)$coefficients[8]<0.05) abline(mod, col='red') else abline(mod, col='lightgrey')
		
		
		ydata= Niche_Families[,i]/specnumber(families)
		xdata=cofactor.df$Area.Weighted.Mean.Percent.Imperviousness
		nonlin.mod=nls(ydata~a*(xdata^2)+b*xdata+c, start=list(a=1, b=1, c=1))
			#lines(sort(xdata), predict(nonlin.mod)[order(xdata)])
		
		if(lrtest(mod, nonlin.mod)$Pr[2]<0.05) lines(sort(xdata), predict(nonlin.mod)[order(xdata)]) else if(summary(mod)$coefficients[8]<0.05) abline(mod, col='red') else abline(mod, col='lightgrey')
		}
		
		


	#Beta functional div -- how different are transects within a site, with respect to functional diversity?
		#use presAbs version of Niche_Families, so as to count each unique niche only once?
		NichePresAbs=Niche_Families ; NichePresAbs[NichePresAbs>0]<-1


	NicheDistance=NA; NicheBetaList=list(NA)  #create empty objects for results .  betaList will keep all pairwise interactions ; "distance" will use a parallel calculation for the same idea
		index=1
		for (i in unique(temp$Site.code)){
			site.df= NichePresAbs[temp$Site.code%in%i,]
				NicheBetaList[[index]]<-as.vector(betadiver(site.df, method=1))
				NicheDistance[index]=mean(as.vector(vegdist(site.df, method="jaccard")))
				index=index+1
			}

	NicheBeta.df<-as.data.frame(NicheBetaList) ; names(NicheBeta.df)<-unique(temp$Site.code) ; NicheBeta.df[2:3,3]<-NA  #save as data frame; replace repeated values as NA to reflect fact that a two-sample comparison has only a single beta value (that sample didn't have all three transects amplify)
	
	NicheBeta<-colMeans(NicheBeta.df, na.rm=T)  #keep means as "beta"

	plot(NicheBeta~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
		summary(lm(NicheBeta~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness)))


	#Gamma diversity (site-level)
		#nearly all functional groups occur at all sites at least once; real differences are at the transect level
			specnumber(aggregate(Niche_Families, by=list(substr(row.names(Niche_Families),1,3)), sum))


#more- and less-urban sites sort out well by functional distances as well
pca.niche<-prcomp(Niche_Families)
#			pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/SecondSubmission/PCA_LifeHistory.pdf", width=7, height=7)
			ordiplot(pca.niche, display="sites", main="Family-Level Life-History Categories")
			samplingsites=substr(row.names(Niche_Families),1,3)
				ordihull(pca.niche,groups= samplingsites,draw="polygon",col=gg_colors(8),label=F)
						ordihull(pca.niche,groups= samplingsites,draw="polygon",col= gg_colors(2)[1],label=F, show.groups=unique(samplingsites)[c(1,5,6,7)])
						ordihull(pca.niche,groups= samplingsites,draw="polygon",col= gg_colors(2)[2],label=F, show.groups=unique(samplingsites)[c(2:4,8)])
				orditorp(pca.niche,display="sites")
										legend("topright", c("More Urban", "Less Urban"), fill=gg_colors(2), cex=0.8)
#			dev.off()

			#niche loadings:
			pca.niche[[2]][,c(1:2)]



			# ordiplot(pca.niche, display="sites")
			# urbancategory=cofactor.df$category.x
						# ordihull(pca.niche,groups= urbancategory,draw="polygon",col= gg_colors(2)[1],label=F, show.groups="more urban")
						# ordihull(pca.niche,groups= urbancategory,draw="polygon",col= gg_colors(2)[2],label=F, show.groups="less urban")
						# orditorp(pca.niche,display="sites")
						# legend("topright", c("More Urban", "Less Urban"), fill=gg_colors(2), cex=0.8)


########################################################################
#####Human footprint
########################################################################
	#genera=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/AnalysisFinal/allGenus_mod.csv", header=F)
		Imperv=cofactor.df$Area.Weighted.Mean.Percent.Imperviousness
		#look for cultivated species:
		#Modiolus richness is signif correlated w Imperv
		#summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Modiolus",1]])~Imperv, family='quasipoisson'))  #cultivated?
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Panopea",1]])~Imperv, family='poisson'))   #cultivated, but obv also naturally occurring
				
		#Human richness same
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Homo",1]])~ Imperv, family='poisson'))
		
		#putative invasives/introduced/reintroduced
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Mercenaria",1]])~Imperv, family='quasipoisson'))  #marg
		summary(glm(abs(cofactor.df[,genera[genera[,3]=="Nuttallia",1]])~Imperv, family='quasipoisson'))  #marg
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Anguinella",1]])~Imperv, family='quasipoisson'))  #nonsig 
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Mya",1]])~Imperv, family='poisson'))  #highly sig 
		
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Bos",1]])~Imperv, family='quasipoisson'))  #sig 
		summary(glm(specnumber(cofactor.df[,genera[genera[,3]=="Gallus",1]])~Imperv, family='quasipoisson'))  # not sig, but tiny numbers of occurrences.  Canidae, similarly, has only one occur


########################################################################
####NATURAL HISTORY COMPARISON
#using Fams_w_NatHist and manual_FamilyCounts, from above
########################################################################
				
				#note: keeping this for possible future use to compare more/less urban sites, rather than comparing manual and eDNA counts, as here
				table(eDNANatHist$Category.infauna) ; table(ManualNatHist$Category.infauna)
					binom.test(c(39,0), p=23/135)  #are the manual and eDNA taxa drawn from the same distribution, in terms of frequency of infaunal taxa?  no.  use binomial test, with expected probability set by the observed prob in the eDNA taxa


				eDNA_infauna=eDNANatHist$Taxon[eDNANatHist$Category.infauna==1] #infaunal taxa
				eDNA_epifauna=eDNANatHist$Taxon[eDNANatHist$Category.epifauna==1] #epifaunal taxa
				eDNA_motile=eDNANatHist$Taxon[eDNANatHist$Mobility.motile ==1] #motile
				eDNA_sessile=eDNANatHist$Taxon[eDNANatHist$Mobility.sessile ==1] #sessile
				eDNA_fish=eDNANatHist$Taxon[grep("Fish", eDNANatHist $Description)]
				eDNA_inverts=eDNANatHist$Taxon[grep("Invert", eDNANatHist $Description)]
				eDNA_intertidal=eDNANatHist$Taxon[eDNANatHist$Habitat.intertidal ==1] #intertidal taxa
				eDNA_subtidal=eDNANatHist$Taxon[eDNANatHist$Habitat.subtidal ==1] #subtidal taxa
				eDNA_terrestrial=eDNANatHist$Taxon[eDNANatHist$Habitat.terrestrial ==1] #subtidal taxa
				
		summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_intertidal]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  # correlation bw imperviousness and intertidal OTU richness
		summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_subtidal]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  # correlation bw imperviousness and subtidal OTU richness
		summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_terrestrial]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  # no correlation bw imperviousness and terrestrial OTU richness
		
		
		summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in%eDNA_fish]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #no correlation bw imperviousness and fish OTU richness
		summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_inverts]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #pos correlation bw imperviousness and invert OTU richness
		
				
				summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in%eDNA_infauna]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #correlation bw imperviousness and infaunal OTU richness
				summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in%eDNA_epifauna]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #correlation bw imperviousness and epifaunal OTU richness
				summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_motile]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #correlation bw imperviousness and motile OTU richness
				summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_sessile]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #correlation bw imperviousness and sessile OTU richness
	
	#sessile OTUs triple
	nh.temp<-data.frame(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in% eDNA_sessile]]), cofactor.df$Area.Weighted.Mean.Percent.Imperviousness); names(nh.temp)=c("sessileOTUs", "Imperv")
	mean(nh.temp$sessileOTUs[nh.temp$Imperv<25])
	mean(nh.temp$sessileOTUs[nh.temp$Imperv>25])
		
		
				summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in%eDNANatHist$Taxon[eDNANatHist$Category.demersal ==1]]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #correlation bw imperviousness and demersal OTU richness
				summary(lm(specnumber(cofactor.df[,master$OTU_name_swarm[master$OTU_taxon_family%in%eDNANatHist$Taxon[eDNANatHist$Category.demersal ==0]]])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))  #correlation bw imperviousness and non-demersal OTU richness






########################################################################
#####FOR SUPPLEMENT############
########################################################################
	#with only top 100 OTUs
	summary(lm(specnumber(temp[,names(head(sort(colSums(temp[,dna_index]), decreasing=T), 100))])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
	
	#with only bottom 500 OTUs
	summary(lm(specnumber(temp[,names(head(sort(colSums(temp[,dna_index]), decreasing=F), 500))])~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
	
	#with raw data (no rarefaction or decontamination)
	z=t(OTUs[, PSeelgrassSites_tags]) ; row.names(z)=PSeelgrassSites$sample_name
					row.names(z)=gsub("1A","_1", row.names(z))
					row.names(z)=gsub("2A","_2", row.names(z))
					row.names(z)=gsub("3A","_3", row.names(z))
					row.names(z)=gsub("3B","_3", row.names(z))
	z1=aggregate(z, list(row.names(z)), mean) ; names(z1)[1]="Sample.ID"
	z2=merge(z1, cofactor.df[,c(2,1737)], "Sample.ID")
	summary(lm(specnumber(z2[,2:6198])~z2[,6199]))

	##plot for supplement; data subsets and superset
#pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/Scatter3x1_Supplement_dataSubsets.pdf", height=7 , width=5)
	par(mfrow=c(3,1))
	
		head.df=data.frame(specnumber(temp[,names(head(sort(colSums(temp[,dna_index]), decreasing=T), 100))])); names(head.df)[1]="richness"
		head.means=aggregate(head.df$richness, list(substr(row.names(head.df),1,3)), mean)[,2]
	plot(head.df$richness~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, pch=19, xlab="Watershed Imperviousness (Percent)", ylab="OTU Richness", main="100 Most Common OTUs" , col=gg_colors(2)[2], cex=0.6)
	points(head.means ~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness), pch=19, cex=1.5, col=gg_colors(2)[2])
		mod1=lm(head.means ~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
		abline(mod1, col=gg_colors(1))
		text(40, 20, paste0("R^2=",round(summary(mod1)$r.squared,3)," \\np=",round(summary(mod1)$coefficients[8], 3)))
	
		tail.df=data.frame(specnumber(temp[,names(head(sort(colSums(temp[,dna_index]), decreasing=F), 500))])) ; names(tail.df)[1]="richness"
		tail.means=aggregate(tail.df$richness, list(substr(row.names(tail.df),1,3)), mean)[,2]
	plot(tail.df$richness~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, pch=19, xlab="Watershed Imperviousness (Percent)", ylab="OTU Richness", main="500 Least Common OTUs", col=gg_colors(2)[2], cex=0.6)
	points(tail.means ~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness), pch=19, cex=1.5, col=gg_colors(2)[2])
		mod2=lm(tail.means ~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
		abline(mod2, col=gg_colors(1))
		text(40, 20, paste0("R^2=",round(summary(mod2)$r.squared,3)," \\np=",round(summary(mod2)$coefficients[8], 4)))
	
		all.df=data.frame(specnumber(z1[,2:6198])) ; names(all.df)[1]="richness"
		all.means=aggregate(all.df $richness, list(substr(row.names(tail.df),1,3)), mean)[,2]
	plot(all.df$richness~cofactor.df$Area.Weighted.Mean.Percent.Imperviousness, pch=19, xlab="Watershed Imperviousness (Percent)", ylab="OTU Richness", main="Raw OTU Data (N = 6198) \\nWithout Normalization or Decontamination", col=gg_colors(2)[2], cex=0.6)
	points(all.means~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness), pch=19, cex=1.5, col=gg_colors(2)[2])
		mod3=lm(all.means~unique(cofactor.df$Area.Weighted.Mean.Percent.Imperviousness))
		abline(mod3, col=gg_colors(1))
		text(40, 150, paste0("R^2=",round(summary(mod3)$r.squared,3)," \\np=",round(summary(mod3)$coefficients[8], 4)))
#	dev.off()

	


#####taxon stacked area chart 
		
				#create df for plotting with a subset of data, picking taxon of interest
			forplot=merge(cofactor.df[,c(2, 1742, 1737)], classes1, "Sample.ID") ; row.names(forplot)=row.names(classes1) ; forplot=forplot[,-1]
					forplotMax=aggregate(forplot, list(substr(row.names(forplot),1,3)), max)
						row.names(forplotMax)= forplotMax[,1] ; forplotMax = forplotMax[,-1]
					forplotMax = forplotMax[order(forplotMax$Area.Weighted.Mean.Percent.Imperviousness),]
		#####FIGURES####
#			pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/SecondSubmission/SupplFig_AreaStackedImperviousness.pdf" , width=9.6, height=5.9)
			plot.stacked(forplotMax $Area.Weighted.Mean.Percent.Imperviousness, forplotMax[,3:ncol(forplotMax)][,which(colSums(forplotMax[,3:ncol(forplotMax)])>20)], col= gg_colors(length(which(colSums(forplotMax[,3:ncol(forplotMax)])>20))), ylab="OTU Richness by Class", xlab="Imperviousness (Percent)", main="eDNA OTU Richness by Class")
				abline(v=forplotMax $Area.Weighted.Mean.Percent.Imperviousness, col="lightgrey")
				legend(min(forplotMax $Area.Weighted.Mean.Percent.Imperviousness),350, names(which(colSums(forplotMax[,3:ncol(forplotMax)])>20)), fill= gg_colors(length(which(colSums(forplotMax[,3:ncol(forplotMax)])>20))), cex=0.5, ncol=2, bg="white")
#			dev.off()
			
			forplotMax = forplotMax[order(forplotMax$Population.density),]
#			pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/StackedArea_humanPop_w_vLines.pdf" , width=9.6, height=5.9)
			plot.stacked(forplotMax $Population.density, forplotMax[,3:ncol(forplotMax)][,which(colSums(forplotMax[,3:ncol(forplotMax)])>20)], col= gg_colors(length(which(colSums(forplotMax[,3:ncol(forplotMax)])>20))), ylab="OTU Richness by Class", xlab="Human Population Density", main="eDNA OTU Richness by Class")
				abline(v=forplotMax$Population.density, col='lightgrey')
				legend(min(forplotMax $Population.density),350, names(which(colSums(forplotMax[,3:ncol(forplotMax)])>20)), fill= gg_colors(length(which(colSums(forplotMax[,3:ncol(forplotMax)])>20))), cex=0.5, ncol=2, bg="white")
#			dev.off()



####FIGURE -- multipanel trends; N families per class; POISSON regression
			focalClasses=names(table(fams2classes)[table(fams2classes)>5])
	#pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Figures/FIGURE_NFamilies_Class_20160204.pdf", width=6, height=7)			
			par(mfrow=c(4,2))
			par(mar=c(3,4,2,2))
			par(oma=c(2,2,1,2))
			par(mgp=c(2,1,0))
			#mean richness by imperv
			for (i in 1:(length(focalClasses)-2)){
					df=data.frame(specnumber(fams_temp[, which(fams2classes==focalClasses[i])+1792]), cofactor.df $Area.Weighted.Mean.Percent.Imperviousness) ; names(df)=c("Richness", "Imperv")
					meanRichness=aggregate(df$Richness, list(df$Imperv), mean) ; meanRichness= meanRichness[match(unique(df$Imperv), meanRichness[,1]),2]
					plot(Richness~Imperv, data=df, main= focalClasses[i], col=gg_colors(2)[2], pch=19, ylab="Number of Familes", xlab="", cex=0.6)
					points(meanRichness~unique(df$Imperv), col=gg_colors(2)[2], pch=19, cex=1.5 )
						mod=glm(round(meanRichness)~unique(df$Imperv), data=df, family='quasipoisson')
						if(summary(mod)$coefficients[8]<0.05) 
							lines(sort(data.frame(unique(df$Imperv),predict.glm(mod, type='response'))), col=gg_colors(1))
						if(summary(mod)$coefficients[8]>0.05) 
							lines(sort(data.frame(unique(df$Imperv),predict.glm(mod, type='response'))), col='lightgrey', lty=2)
				}
			for (i in 6:7){
					df=data.frame(specnumber(fams_temp[, which(fams2classes==focalClasses[i])+ 1792]), cofactor.df $Area.Weighted.Mean.Percent.Imperviousness) ; names(df)=c("Richness", "Imperv")
					meanRichness=aggregate(df$Richness, list(df$Imperv), mean) ; meanRichness= meanRichness[match(unique(df$Imperv), meanRichness[,1]),2]
					plot(Richness~Imperv, data=df, main= focalClasses[i], col=gg_colors(2)[2], pch=19, ylab="Number of Familes", xlab="Impervious Surface Cover", cex=0.6)
					points(meanRichness~unique(df$Imperv), col=gg_colors(2)[2], pch=19, cex=1.5)
						mod=glm(round(meanRichness)~unique(df$Imperv), data=df, family='quasipoisson')
						if(summary(mod)$coefficients[8]<0.05) 
							lines(sort(data.frame(unique(df$Imperv),predict.glm(mod, type='response'))), col=gg_colors(1))
						if(summary(mod)$coefficients[8]>0.05) 
							lines(sort(data.frame(unique(df$Imperv),predict.glm(mod, type='response'))), col='lightgrey', lty=2)
				}
	#dev.off()




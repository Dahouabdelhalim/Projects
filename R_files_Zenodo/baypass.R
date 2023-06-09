#!/usr/bin/Rscript
# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of pyRona.
# pyRona is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyRona is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyRona. If not, see <http://www.gnu.org/licenses/>.

require(corrplot)
require(ape)
require(mvtnorm)
require(geigen)


## Source the baypass R functions

# Full path to "baypass_utils.R" provided with baypass.
# Type: str
# Example: "~/Software/Science/baypass_2.1/utils/baypass_utils.R"
source("./R/baypass_utils.R")

## Define some variables. This is where we define local files & paths.

# Path to baypass binary.
# Type: str
# Example: "~/Software/Science/baypass_2.1/sources/g_baypass"
baypass_executable = "/usr/local/bin/baypass.2.1/sources/g_baypass"

# Path to file with population names.
# Type: str
# Example: "~/baypass_analysis/popnames_.txt"
popname_file = "popname"

# Path to ENVFILE, as described in the manual.
# Type: str
# Example: "~/baypass_analysis/ENVFILE"
envfile = "bio02"

# Path to geno file, as described in the manual.
# Type: str
# Example: "~/baypass_analysis/data.baypass"
geno_file = "./cs.north.baypass.txt"

# Results perfix, like the name of your analysis.
# Type: str
# Example:"Qsuber"
prefix = "cs"

# Number of populations.
# Type: int
# Example: 16
num_pops = 7

# Number os SNPs.
# Type: int
# Example: 325
num_SNPs = 117626
snp.info.file = "cs.snp.id"
# Number of threads to use.
# Type: int
# Example: 8
num_threads = 50

# Should we normalize the data on ENVFILE?
# Type: BOOL (TRUE or FALSE)
# Set to false if your ENVFILE is already normalized, TRUE if otherwise
scale_cov = TRUE


## Everything below this point should be fully automated.

scalecov <- if (scale_cov) {" -scalecov "} else {""}

basepath = dirname(geno_file)
coredir = paste(basepath, "/core/", sep="")
mcmc_coredir = paste(basepath, "/bio02_mcmc_core/", sep="")
mcmc_stddir = paste(basepath, "/bio02_mcmc_std/", sep="")
mcmc_auxdir = paste(basepath, "/bio02_mcmc_aux/", sep="")

dir.create(coredir)
dir.create(mcmc_coredir)
dir.create(mcmc_stddir)
dir.create(mcmc_auxdir)

core_omega_file = paste(coredir, prefix, "_core_mat_omega.out", sep="")
core_pi_xtx_file = paste(coredir, prefix, "_core_summary_pi_xtx.out", sep="")
core_summary_beta_params = paste(coredir, prefix, "_core_summary_beta_params.out", sep="")
pod_mat_omega = paste(coredir, prefix, "_core_POD_mat_omega.out", sep="")
pod_summary_beta_params = paste(coredir, prefix, "_core_POD_summary_beta_params.out", sep="")
pod_summary_pi_xtx = paste(coredir, prefix, "_core_POD_summary_pi_xtx.out", sep="")
covis_summary_betai_reg = paste(mcmc_coredir, prefix, "_mcmc_core_summary_betai_reg.out", sep="")
covis2_summary_betai_reg = paste(mcmc_coredir, prefix, "_mcmc_core2_summary_betai_reg.out", sep="")
covmcmc_summary_betai = paste(mcmc_stddir, prefix, "_mcmc_std_summary_betai.out", sep="")
covmcmc_summary_pi_xtx = paste(mcmc_stddir, prefix, "_mcmc_std_summary_pi_xtx.out", sep="")
covaux_summary_betai = paste(mcmc_auxdir, prefix, "_mcmc_aux_summary_betai.out", sep="")
covaux_summary_pi_xtx = paste(mcmc_auxdir, prefix, "_mcmc_aux_summary_pi_xtx.out", sep="")

#### Run the first command:
command1 = paste(baypass_executable, " -npop ", num_pops, " -gfile ",
                 geno_file, " -outprefix ", coredir, prefix, "_core",
                 " -nthreads ", num_threads, sep="")
system(command=command1)


## upload estimate of omega
omega=as.matrix(read.table(core_omega_file))
pop.names = c(as.matrix(read.table(popname_file)))
#
dimnames(omega)=list(pop.names,pop.names)
## Compute and visualize the correlation matrix
cor.mat=cov2cor(omega)
pdf(file=paste0(coredir,"omega_corr.pdf"),width = 16, height = 16)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))
dev.off()

## Visualize the correlation matrix as hierarchical clustering tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
pdf(file=paste0(coredir, "Hier_clust_tree.pdf"),width = 16, height = 8)
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

## Estimates of the XtX differentiation measures
anacore.snp.res=read.table(core_pi_xtx_file,h=T)
pdf(file=paste0(coredir, "XtX_diff.pdf"),width = 16, height = 8)
plot(anacore.snp.res$M_XtX)
dev.off()

# Get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution
pi.beta.coef=read.table(core_summary_beta_params,h=T)$Mean
# Upload the original data to obtain total allele count
current.data<-geno2YN(geno_file)
# Create the POD
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=num_SNPs,
                           sample.size=current.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0,suffix="btapods")


file.rename("G.btapods", paste(coredir, "G.btapods", sep=""))

###
command2 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", coredir,
                 "G.btapods", " -outprefix ", coredir, prefix, "_core_POD",
                 " -nthreads ", num_threads, sep="")
system(command=command2)


#######################################################
# Sanity Check: Compare POD and original data estimates
#######################################################
# Get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table(pod_mat_omega))
pdf(file=paste0(coredir, "Omega_estimate_from_POD.pdf"),width = 16, height = 8)
plot(pod.omega,omega) ; abline(a=0,b=1)
dev.off()
fmd.dist(pod.omega,omega)

# Get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table(pod_summary_beta_params,h=T)$Mean
pdf(file=paste0(coredir, "POD_estimates.pdf"),width = 16, height = 8)
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)
dev.off
#
#######################################################
# XtX calibration
#######################################################
# Get the pod XtX
pod.xtx=read.table(pod_summary_pi_xtx,h=T)$M_XtX
# Compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.99)
# Add the thresh to the actual XtX plot
snp.info=read.table(snp.info.file,header = T)
chr.id=as.factor(snp.info$CHROM)
mycolors=rep(rainbow(2),20)

pdf(file=paste0(coredir, "XtX_POD_diff.pdf"),width = 16, height = 8)
plot(anacore.snp.res$M_XtX,col=mycolors[chr.id])
abline(h=pod.thresh,lty=2)
dev.off()

anacore.snp.res.1 <- cbind(snp.info,anacore.snp.res)
write.table(anacore.snp.res.1,"anacore.snp.restults",row.names=T,sep=" ",quote=F)
anacore.snp.res.row <- read.table("anacore.snp.restults",sep = " ",header = T)

outliers <- anacore.snp.res.1[,c("ID")][which(anacore.snp.res.1$M_XtX>pod.thresh)]
anacore.snp.res.1[,"ID"][which(anacore.snp.res.1$M_XtX>pod.thresh)]
xtx.outliers <- anacore.snp.res.1[,"ID"][which(anacore.snp.res.1[,"M_XtX"]>pod.thresh)]
write.table(xtx.outliers,file=paste0(coredir,"bio02.xtx.ouliers.snp.txt"),sep = " ",row.names = F,col.names = F,quote = F)
###

###
command3 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_coredir, prefix, "_mcmc_core",
                 " -nthreads ", num_threads, " -efile ", envfile, scalecov,
                 sep="")
system(command=command3)

covis.snp.res=read.table(covis_summary_betai_reg,h=T)
graphics.off()

pdf(file=paste0(mcmc_coredir, "BFs_layout.pdf"),width = 16, height = 8)
layout(matrix(1:3,3,1))
plot(covis.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()


###
command4 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_coredir, prefix, "_mcmc_core2",
                 " -nthreads ", num_threads, " -efile ", envfile, scalecov,
                 " -omegafile ", core_omega_file, sep="")
system(command=command4)


covis2.snp.res=read.table(covis2_summary_betai_reg,h=T)
graphics.off()

pdf(file=paste0(mcmc_coredir, "BFs_layout_pass2.pdf"),width = 16, height = 8)
layout(matrix(1:3,3,1))
plot(covis2.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis2.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis2.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

###
command5 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_stddir, prefix, "_mcmc_std",
                 " -nthreads ", num_threads, " -efile ", envfile, scalecov,
                 " -omegafile ", core_omega_file, " -covmcmc", sep="")
system(command=command5)

covmcmc.snp.res=read.table(covmcmc_summary_betai,h=T)
covmcmc.snp.xtx=read.table(covmcmc_summary_pi_xtx,h=T)$M_XtX
graphics.off()

pdf(file=paste0(mcmc_stddir, "BFs_layout.pdf"),width = 16, height = 8)
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res$eBPmc,xlab="SNP",ylab="eBPmc")
plot(covmcmc.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covmcmc.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")
dev.off()


###
command6 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_auxdir, prefix, "_mcmc_aux",
                 " -nthreads ", num_threads, " -efile ", envfile, scalecov,
                 " -omegafile ", core_omega_file, " -auxmodel", sep="")
system(command=command6)

covaux.snp.res.raw=read.table(covaux_summary_betai,h=T)
covaux.snp.xtx.raw=read.table(covaux_summary_pi_xtx,h=T)$M_XtX
snp.info=read.table(snp.info.file,header = T)
covaux.snp.res=cbind(snp.info,covaux.snp.res.raw)
covaux.snp.xtx=cbind(snp.info,covaux.snp.xtx.raw)

mycolors=rep(rainbow(2),20)


pdf(file=paste0(mcmc_auxdir, "BFs_layout.pdf"),width = 16, height = 16)
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB.,xlab="SNP",ylab="BFmc (in dB)",col=mycolors[chr.id])
plot(covaux.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"),col=mycolors[chr.id])
plot(covaux.snp.xtx$covaux.snp.xtx.raw,xlab="SNP",ylab="XtX corrected for SMS",col=mycolors[chr.id])
dev.off()

covaus.outliers=covaux.snp.res[,"ID"][which(covaux.snp.res[,"BF.dB."]>20)]
write.table(covaus.outliers,file=paste0(mcmc_auxdir,"bio02.covaus.ouliers.BF20.snp.txt"),sep = " ",row.names = F,col.names = F,quote = F)

library(deSolve)
library(ggplot2)
library(data.table)
library(patchwork)


multispecies_multipatch_SIZ = function(t, x, params){
	# Multispecies, multipatch model: 2 species by 2 patches
	#
	# Parameters
	# ----------
	# t : array
	# 	time array
	# x : array
	# 	State variables
	# params : list
	# 	List of parameters
	#
	# Returns
	# -------
	# : ODE solution

	# Species by patch indexing. So S21 is susceptible species 2 in patch 1.
	S11 = x['S11']
	S21 = x['S21']
	I11 = x['I11']
	I21 = x['I21']
	Z1 = x['Z1']
	S12 = x['S12']
	S22 = x['S22']
	I12 = x['I12']
	I22 = x['I22']
	Z2 = x['Z2']

	beta = params[['beta']] # Transmission rate
	d = params[['d']] # Host death rate
	nu = params[['nu']] # recovery rate
	alpha = params[['alpha']] # death rate
	lambda = params[['lambda']] # shedding rate
	b = params[['b']] # Total birth rate
	gamma = params[['gamma']] # Pathogen decay
	phi = params[['phi']] # Dispersal rate
	A = params[['A']] # Relative area matrix

	# Patch 1
	# Species by patch indexing
	dS11 = b[1, 1] - d[1, 1]*S11 - beta[1, 1]*Z1*S11 + nu[1, 1]*I11 - phi[1, 1]*S11 + phi[1, 2]*S12*A[2, 1]
	dS21 = b[2, 1] - d[2, 1]*S21 - beta[2, 1]*Z1*S21 + nu[2, 1]*I21 - phi[2, 1]*S21 + phi[2, 2]*S22*A[2, 1]
	dI11 = beta[1, 1]*Z1*S11 - I11*(nu[1, 1] + alpha[1, 1] + d[1, 1]) - phi[1, 1]*I11 + phi[1, 2]*I12*A[2, 1]
	dI21 = beta[2, 1]*Z1*S21 - I21*(nu[2, 1] + alpha[2, 1] + d[2, 1]) - phi[2, 1]*I21 + phi[2, 2]*I22*A[2, 1]
	dZ1 = lambda[1, 1]*I11 + lambda[2, 1]*I21 - gamma[1]*Z1

	# Patch 2
	dS12 = b[1, 2] - d[1, 2]*S12 - beta[1, 2]*Z2*S12 + nu[1, 2]*I12 - phi[1, 2]*S12 + phi[1, 1]*S11*A[1, 2]
	dS22 = b[2, 2] - d[2, 2]*S22 - beta[2, 2]*Z2*S22 + nu[2, 2]*I22 - phi[2, 2]*S22 + phi[2, 1]*S21*A[1, 2]
	dI12 = beta[1, 2]*Z2*S12 - I12*(nu[1, 2] + alpha[1, 2] + d[1, 2]) - phi[1, 2]*I12 + phi[1, 1]*I11*A[1, 2]
	dI22 = beta[2, 2]*Z2*S22 - I22*(nu[2, 2] + alpha[2, 2] + d[2, 2]) - phi[2, 2]*I22 + phi[2, 1]*I21*A[1, 2]
	dZ2 = lambda[1, 2]*I12 + lambda[2, 2]*I22 - gamma[2]*Z2

	return(list(c(dS11, dS21, dI11, dI21, dZ1, dS12, dS22, dI12, dI22, dZ2)))

}


get_params = function(type){
	# Rows are species, columns are patches

	if(type == "reversed"){

		# Set parameters such that observed maintenance potential without
		# correcting for spatial context is opposite of the true maintenance
		# potential

		# Transmission rate
		beta = matrix(c(0.4, 2, 
						1, 0.4), nrow=2, ncol=2, byrow=TRUE)

		# Host death rate
		d = matrix(c(1, 1,
					 1, 1), nrow=2, ncol=2, byrow=TRUE)

		# Host recovery rate
		nu = matrix(c(0.5, 0.5,
					  0.5, 0.5), nrow=2, ncol=2, byrow=TRUE)

		# Host disease-induced mortality rate
		alpha = matrix(c(0, 0,
					     0, 0), nrow=2, ncol=2, byrow=TRUE)

		# Pathogen shedding rate
		lambda = matrix(c(1, 1,
					      1, 1), nrow=2, ncol=2, byrow=TRUE)

		# Host total birth rate
		b = matrix(c(3, 3,
					 3, 3), nrow=2, ncol=2, byrow=TRUE)

		# Pathogen decay rate
		gamma = c(1, 1)

		# Host dispersal rate
		phi = matrix(c(0.01, 10,
					   10, 0.01), nrow=2, ncol=2, byrow=TRUE)

		# Relative areas
		A = matrix(c(1, 1,
					 1, 1), nrow=2, ncol=2, byrow=TRUE)

		params = list(beta=beta, d=d, nu=nu, alpha=alpha, 
					  lambda=lambda, b=b, gamma=gamma, 
					  phi=phi, A=A)
	} else if(type == "matches"){

		# Set parameters such that observed prevalence matches the ordering
		# of R0

		# Transmission rate
		beta = matrix(c(0.3, 2, 
						1, 0.4), nrow=2, ncol=2, byrow=TRUE)

		# Host death rate
		d = matrix(c(1, 1,
					 1, 1), nrow=2, ncol=2, byrow=TRUE)

		# Host recovery rate
		nu = matrix(c(0.5, 0.5,
					  0.5, 0.5), nrow=2, ncol=2, byrow=TRUE)

		# Host disease-induced mortality rate
		alpha = matrix(c(0, 0,
					     0, 0), nrow=2, ncol=2, byrow=TRUE)

		# Pathogen shedding rate
		lambda = matrix(c(1, 1,
					      1, 1), nrow=2, ncol=2, byrow=TRUE)

		# Host total birth rate
		b = matrix(c(3, 3,
					 3, 3), nrow=2, ncol=2, byrow=TRUE)

		# Pathogen decay rate
		gamma = c(1, 1)

		# Host dispersal rate
		phi = matrix(c(1, 1,
					   1, 1), nrow=2, ncol=2, byrow=TRUE)

		# Relative areas
		A = matrix(c(1, 1,
					 1, 1), nrow=2, ncol=2, byrow=TRUE)

		params = list(beta=beta, d=d, nu=nu, alpha=alpha, 
					  lambda=lambda, b=b, gamma=gamma, 
					  phi=phi, A=A)

	} else{

		params = list()
	}

	return(params)
}


types = c("matches", "reversed")

# Loop over options
for(t in 1:2){

	type = types[t]

	# Get parameters and simulate
	params = get_params(type)
	init = c(S11=2.9, S21=2.9, I11=0.1, I21=0, Z1=0, 
			 S12=2.9, S22=2.9, I12=0.1, I22=0, Z2=0) 

	times = seq(0, 100, len=500)
	res = as.data.table(data.frame(ode(init, times, multispecies_multipatch_SIZ, params)))


	# Format for plotting
	res_melt = melt(res, id.var="time", variable.name="svar", value.name="val")
	res_S = res_melt[svar %like% "S"]
	res_I = res_melt[svar %like% "I"]
	res_S$prev = res_I$val / (res_I$val + res_S$val) # Compute prevalence
	res_S$density = res_I$val + res_S$val

	# Plot for Patch 1 
	patch1_dat = droplevels(res_S[svar %like% "S.1"])
	levels(patch1_dat$svar) = c("Species 1", "Species 2")

	prev_plot1 = ggplot(patch1_dat) + geom_line(aes(x=time, y=prev, color=svar)) + 
							   ylab("Prevalence") + xlab("Time (arbitrary units)") +
							   theme_classic() + theme(legend.position=c(0.5, 0.2), 
							   						   legend.title = element_blank())
	prev_plot1
	print(tail(patch1_dat[time == 100]))
	# ggsave("../results/patch1_plot_opposite.pdf", width=4, height=3)

	# Plot for Patch 2
	patch2_dat = droplevels(res_S[svar %like% "S.2"])
	levels(patch2_dat$svar) = c("Species 1", "Species 2")

	prev_plot2 = ggplot(patch2_dat) + geom_line(aes(x=time, y=prev, color=svar)) + 
							   ylab("Prevalence") + xlab("Time (arbitrary units)") +
							   theme_classic() + theme(legend.position=c(0.5, 0.2), 
							   						   legend.title = element_blank())
	prev_plot2

	# True species by patch level R0
	R0_true_part = with(params, (b / d)*beta*lambda / (alpha + nu + d))
	R0_true = with(params, t(t(R0_true_part) * (1 / gamma)))
	R0_true

	# Community-corrected R0
	R0_corrected = array(NA, dim=c(2, 2))
	for(s in 1:2){
		for(p in 1:2){

			id_focal = paste0("S", s, p)
			t_focal = res_S[(svar == id_focal) & time == 100]
			focal_dens = t_focal$density
			focal_prev = t_focal$prev
			focal_lambda = params$lambda[s, p]

			# Only correct for community bias
			correct_fact = 0
			for(i in 1:2){

				id_other = paste0("S", i, p)
				t_other = res_S[(svar == id_other) & time == 100]
				other_dens = t_other$density
				other_prev = t_other$prev
				other_lambda = params$lambda[i, p]

				correct_fact = correct_fact + (other_dens / focal_dens)*(other_prev / focal_prev)*(other_lambda / focal_lambda)
			}

			R0_corrected[s, p] = 1 / ((1 - focal_prev)*correct_fact)
		}
	}
	R0_corrected

	R0_dat = data.frame(R0true=R0_true[, 1],
						R0false=R0_corrected[, 1],
						species=c("Species 1", "Species 2"))

	# Plot true and observed maintenance potential
	r0plot1 = ggplot(R0_dat) + geom_point(aes(x=R0false, y=R0true, color=species), size=4) + 
							   ylab("R0,s,p\\n(True maintenance potential\\nwith spatial context)") +
							   xlab("R0,s\\n(Biased maintenance potential\\nwithout spatial context)") +
							   geom_abline(aes(slope=0, intercept=1), linetype="dashed") +
							   geom_vline(aes(xintercept=1), linetype="dashed") +
							   theme_classic() + theme(legend.title = element_blank(),
											 			 	   legend.position="none")

	p1_plot = prev_plot1 + r0plot1
	ggsave(paste0("../results/patch", t, "_both_plots_", type, ".pdf"), width=6, height=3)
}














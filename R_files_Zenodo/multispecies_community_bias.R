library(deSolve)
library(ggplot2)
library(data.table)
library(patchwork)

multispecies_SIZ = function(t, x, params){
	# Multispecies SIZ model to test community bias
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

	# Number of species in community
	spp_num = params[['spp_num']]

	d_vect = params[['d']] # Host death rate
	beta_vect = params[['beta']] # Transmission rate
	nu_vect = params[['nu']] # recovery rate
	alpha_vect = params[['alpha']] # death rate
	lambda_vect = params[['lambda']] # shedding rate
	b_vect = params[['b']] # Total birth rate
	gamma = params[['gamma']] # Parasite death rate

	# Extract state variables
	Svals = x[1:spp_num] # Susceptible species
	Ivals = x[(spp_num + 1):(2*spp_num)] # Infected species
	Zvals = x[length(x)] # Pathogen

	# Define the ODEs
	dS_vect = b_vect - d_vect*Svals - beta_vect*Zvals*Svals + nu_vect*Ivals
	dI_vect = beta_vect*Zvals*Svals - (nu_vect + d_vect + alpha_vect)*Ivals
	dZ = sum(lambda_vect*Ivals) - gamma*Zvals

	return(list(c(dS_vect, dI_vect, dZ)))

}



get_params = function(p){
	# Extract different parameter sets for multi-species model with a single patch.
	# 
	# Parameters
	# ----------
	# p : character
	# 	Either 
	#	"equal": No relationship between R0 and prev
	#   "positive": A positive relationship between R0 and prev
	#   "negative": A negative relationships between R0 and prev
	#
	# Returns
	# -------
	# : list of parameters

	if(p == 'equal'){

		# No relationship between R0 and prev
		spp_num = 4
		params = list()
		params[['spp_num']] = spp_num
		params[['b']] = rep(2, spp_num)
		params[['d']] = 1:spp_num
		params[['alpha']] = rep(0, spp_num)
		params[['lambda']] = rep(1, spp_num)
		params[['gamma']] = 1.5
		params[['nu']] = rep(0, spp_num)
		params[['beta']] = 1:4 #rep(1, spp_num)

	} else if(p == "positive"){

		# Positive relationship between R0 and prev
		spp_num = 4
		params = list()
		params[['spp_num']] = spp_num
		params[['b']] = rep(2, spp_num)
		params[['d']] = 1:spp_num
		params[['alpha']] = rep(0, spp_num)
		params[['lambda']] = rep(1, spp_num)
		params[['gamma']] = 1
		params[['nu']] = rep(0.1, spp_num)
		params[['beta']] = seq(2, 0.5, len=spp_num)

	} else {

		# Negative relationship between R0 and prev
		spp_num = 4
		params = list()
		params[['spp_num']] = spp_num
		params[['b']] = rep(2, spp_num)
		params[['d']] = rep(1, spp_num)
		params[['alpha']] = rep(0, spp_num)
		params[['lambda']] = c(4 / 1.5, 1.75 / 3, 1 / 6, 0.01)
		params[['gamma']] = 4
		params[['nu']] = rep(0.1, spp_num)
		params[['beta']] = 1:4
	}

	return(params)
}


make_plot = function(opt, i){
	# Make side by side plots for prevalence dynamics and enzootic prev vs. R0
	#
	# Parameters
	# -----------
	# opt : string
	# 		equal, positive, or negative
	# i : int
	#  Specifies plot formatting
	#
	# Returns
	# -------
	# : ggplot

	params = get_params(opt)
	R0s = ((params$b / params$d)*params$beta * params$lambda) / ((params$alpha + params$d + params$nu)*params$gamma)
	print(R0s)

	spp_num = 4
	init = c(rep(params[['b']] / params[['d']]), rep(0.001, spp_num), 0)
	times = seq(0, 50, len=100)
	res = as.data.table(data.frame(ode(init, times, multispecies_SIZ, params)))
	names = c("time", paste0("S", 1:4), paste0("I", 1:4), "Z")
	colnames(res) = names

	# Calculate prevalence
	for(s in 1:spp_num){
		prev = res[[paste0("I", s)]] / (res[[paste0("I", s)]] + res[[paste0("S", s)]])
		res[, paste0("Prev., species ", s):=prev]
	}

	res_melt = melt(res, id.var="time", variable.name="svar", value.name="val")
	res_prev = res_melt[svar %like% "Prev"]

	prev_plot = ggplot(res_prev) + geom_line(aes(x=time, y=val, color=svar)) + 
							   ylab("Prevalence") + xlab("Time (arbitrary units)") +
							   theme_classic() + 
							   theme(legend.title = element_blank(),
							   		 legend.position="none")


	equil_prev = res_prev[, .(max_prev=max(val)),  by=.(svar)]
	equil_prev$svar = paste("Species", 1:4)

	if(i == 1){
		prev_plot = prev_plot + xlab("")
		rank_plot = ggplot(equil_prev) + geom_point(aes(x=R0s, y=max_prev, color=svar), size=3) + 
						 			 ylab("Enzootic prevalence") + xlab("") +
						 			 theme_classic() + 
						 			 theme(legend.title = element_blank(),
						 			 	   legend.position=c(0.8, 0.4))
	} else if(i == 2) {
		prev_plot = prev_plot + xlab("")
		rank_plot = ggplot(equil_prev) + geom_point(aes(x=R0s, y=max_prev, color=svar), size=3) + 
						 			 ylab("Enzootic prevalence") + xlab("") +
						 			 theme_classic() + 
						 			 theme(legend.title = element_blank(),
						 			 	   legend.position="none")
	} else{
		rank_plot = ggplot(equil_prev) + geom_point(aes(x=R0s, y=max_prev, color=svar), size=3) + 
						 			 ylab("Enzootic prevalence") + xlab("Species R0\\n(i.e., maintenance potential)") +
						 			 theme_classic() + 
						 			 theme(legend.title = element_blank(),
						 			 	   legend.position="none")
	}

	tplot = prev_plot + rank_plot
	return(list(plot=tplot, res=res))

}

# Make community plots
tplot1 = make_plot("positive", 1)[['plot']]
res1 = make_plot("positive", 1)[['res']]
tplot2 = make_plot("equal", 2)[['plot']]
res2 = make_plot("equal", 2)[['res']]
tplot3 = make_plot("negative", 3)[['plot']]

full_plot = (tplot1 / tplot2 / tplot3) + plot_annotation(tag_levels = 'A', tag_suffix=".")
ggsave("../results/community_context_example.pdf", width=7, height=8)





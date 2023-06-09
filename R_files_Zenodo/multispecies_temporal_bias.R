library(deSolve)
library(ggplot2)
library(data.table)
library(patchwork)


seasonal_fxn = function(t, val_max, val_min, period, shift){
	# Sinusoidal transmission function for environmentally forced transmission
	# coefficient
	#
	# Parameters
	# ----------
	# t : numeric
	# 	Time in the season
	# val_max : numeric
	#	Max value of parameter
	# val_min : numeric
	# 	Min value of parameter
	# period : numeric
	# 	Length of cycle
	# shift : numeric
	# 	shift when the cycle peaks
	#
	# Returns
	# -------
	# : seasonal parameter value

	val = ((val_max - val_min) / 2)*(1 - cos(2*pi*(t - shift) / period)) + val_min
	return(val)

}

multispecies_SIZ_temporal = function(t, x, params){
	# Model with temporally varying vital rates
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

	# Time varying transmission
	beta_min = params[['beta_min']]
	beta_max = params[['beta_max']] 
	beta_period = params[['beta_period']] 
	beta_shift = params[['beta_shift']]
	beta_vect = seasonal_fxn(t, beta_max, beta_min, beta_period, beta_shift)

	# Time-varying pathogen decay	
	gamma_min = params[['gamma_min']]
	gamma_max = params[['gamma_max']]
	gamma_period = params[['gamma_period']]
	gamma_shift = params[['gamma_shift']]
	gamma = seasonal_fxn(t, gamma_max, gamma_min, gamma_period, gamma_shift)

	# Time varying shedding
	lambda_min = params[['lambda_min']]
	lambda_max = params[['lambda_max']] 
	lambda_period = params[['lambda_period']] 
	lambda_shift = params[['lambda_shift']]
	lambda_vect = seasonal_fxn(t, lambda_max, lambda_min, lambda_period, lambda_shift)

	nu_vect = params[['nu']] # recovery rate
	alpha_vect = params[['alpha']] # death rate
	b_vect = params[['b']] # Total birth rate

	# Extract state variables
	Svals = x[1:spp_num]
	Ivals = x[(spp_num + 1):(2*spp_num)]
	Zvals = x[length(x)]

	# Define the ODEs
	dS_vect = b_vect - d_vect*Svals - beta_vect*Zvals*Svals + nu_vect*Ivals
	dI_vect = beta_vect*Zvals*Svals - (nu_vect + d_vect + alpha_vect)*Ivals
	dZ = sum(lambda_vect*Ivals) - gamma*Zvals

	return(list(c(dS_vect, dI_vect, dZ)))

}

get_params = function(param_type){
	# Assign parameters for different scenarios of temporal bias
	# 
	# Parameters
	# ----------
	# param_type : character
	# 	Either "matching" or "other". "other" is what is used in the manuscript

	if(param_type == "matching"){

		spp_num = 3
		params = list()
		params[['spp_num']] = spp_num
		params[['b']] = rep(2, spp_num) # host birth rate
		params[['d']] = rep(2, spp_num) # Host death rate
		params[['alpha']] = rep(0, spp_num) # parasite-induced mortality rate
		params[['lambda']] = rep(1, spp_num) # Shedding rate
		params[['nu']] = rep(0, spp_num) # Pathogen decay

		# Pathogen decay
		params[['gamma_min']] = 2.5
		params[['gamma_max']] = 4.5
		params[['gamma_period']] = 50
		params[['gamma_shift']] = 0

		# Transmission rate
		params[['beta_min']] = rep(1, spp_num)  #rep(1, spp_num)
		params[['beta_max']] = rep(1, spp_num)  #rep(1, spp_num)
		params[['beta_period']] = rep(50, spp_num)  #rep(1, spp_num)
		params[['beta_shift']] = c(0, 0, 0)

		# Shedding rate
		params[['lambda_min']] = rep(0.1, spp_num)  #rep(1, spp_num)
		params[['lambda_max']] = rep(10, spp_num)  #rep(1, spp_num)
		params[['lambda_period']] = rep(50, spp_num)  #rep(1, spp_num)
		params[['lambda_shift']] = c(0, 15, 30)
	}

	if(param_type == "other"){

		spp_num = 2
		params = list()
		params[['spp_num']] = spp_num
		params[['b']] = rep(2, spp_num)
		params[['d']] = rep(2, spp_num)
		params[['alpha']] = rep(0, spp_num)
		params[['lambda']] = rep(1, spp_num)
		params[['nu']] = rep(0, spp_num)

		# Pathogen decay. Assuming no fluctuation here
		params[['gamma_min']] = 4.
		params[['gamma_max']] = 4.
		params[['gamma_period']] = 50
		params[['gamma_shift']] = 10

		# Species 2 has higher transmission rate and no fluctuation
		params[['beta_min']] = c(1, 3)
		params[['beta_max']] = c(1, 3) 
		params[['beta_period']] = rep(50, spp_num)  
		params[['beta_shift']] = c(0, 0)

		# Species 1 shedding fluctuates through the year
		# while species 2 shedding stays constant.
		params[['lambda_min']] = c(0.01, 2)
		params[['lambda_max']] = c(12, 2)  
		params[['lambda_period']] = rep(50, spp_num)  
		params[['lambda_shift']] = c(0, 15)
	}

	return(params)

}

params = get_params("other")
spp_num = params[['spp_num']]
init =  c(rep(params[['b']] / params[['d']]), rep(0.001, spp_num), 0)
times = seq(0, 200, len=400)
res = data.table(data.frame(ode(init, times, multispecies_SIZ_temporal, params)))
names = c("time", paste0("S", 1:spp_num), paste0("I", 1:spp_num), "Z")
colnames(res) = names

# Calculate prevalence
for(s in 1:spp_num){
	prev = res[[paste0("I", s)]] / (res[[paste0("I", s)]] + res[[paste0("S", s)]])
	N = (res[[paste0("I", s)]] + res[[paste0("S", s)]])
	res[, paste0("Prev., species ", s):=prev]
	res[, paste0("N", s):=N]
}

res_melt = melt(res, id.var="time", variable.name="svar", value.name="val")
res_prev = res_melt[svar %like% "Prev"]
res_N = res_melt[svar %like% "N"]
res_Z = res_melt[svar %like% "Z"]

# Calculate instantaneous R0 through time.
R0_int = array(NA, dim=c(length(times), spp_num))
for(i in 1:length(times)){

	# Extract pop size
	t = times[i]
	N = res_N[time == t][order(svar)]$val

	betas = seasonal_fxn(t, params[['beta_max']], 
							params[['beta_min']], 
							params[['beta_period']],
							params[['beta_shift']])
	lambdas = seasonal_fxn(t, params[['lambda_max']], 
							  params[['lambda_min']], 
							  params[['lambda_period']],
							  params[['lambda_shift']])
	gamma = seasonal_fxn(t, params[['gamma_max']], 
							  params[['gamma_min']],
							  params[['gamma_period']], 
							  params[['gamma_shift']])

	R0 = betas*N*lambdas / ((params[['alpha']] + params[['nu']] + params[['d']])*gamma)
	R0_int[i, ] = R0

}


# Calculate uncorrected community R0 through time
R0_uncorrected = array(NA, dim=c(length(times), spp_num))
for(t in 1:length(times))  {

	# Extract shedding
	lambdas = seasonal_fxn(times[t], params[['lambda_max']], 
							  params[['lambda_min']], 
							  params[['lambda_period']],
							  params[['lambda_shift']])

	# Extract prevalence
	tprev = res_prev[time == times[t]]$val
	tN = res_N[time == times[t]]$val
	tR0s = array(NA, dim=spp_num)

	for(i in 1:spp_num){

		ratio_sum = 0
		for(j in 1:spp_num){
			tp_ratio = (tprev[j] / tprev[i]) * (tN[j] / tN[i]) * (lambdas[j] / lambdas[i])
			ratio_sum = ratio_sum + tp_ratio
		}
		tR0s[i] = 1 / ((1 - tprev[i])*ratio_sum)
	}

	R0_uncorrected[t, ] = tR0s

}

# Plot the results and compare R0 estimates

R0_df = data.table(R0_int)
R0_df_uncorrected = data.table(R0_uncorrected)
colnames(R0_df) = paste0("Species ", 1:spp_num)
colnames(R0_df_uncorrected) = paste0("Species ", 1:spp_num)
R0_df[, time:=times]
R0_df[, kind:="True\\n(with temporal context)"]
R0_df_uncorrected[, time:=times]
R0_df_uncorrected[, kind:="Biased\\n(w/out temporal context)"]

R0_melted = melt(R0_df, id.vars=c("time", "kind"))
R0_melted_uncorrected = melt(R0_df_uncorrected, id.vars=c("time", "kind"))
R0_all = data.table(rbind(R0_melted, R0_melted_uncorrected))
R0_all$kind = factor(R0_all$kind, levels=c("True\\n(with temporal context)", "Biased\\n(w/out temporal context)"))

tval = 50
prev_plot = ggplot(res_prev[time >= tval]) + geom_line(aes(x=time - tval, y=val, color=svar)) + 
						   ylab("Prevalence") + xlab("Time (arbitrary units)") +
						   theme_classic() + 
						   theme(legend.title = element_blank(),
						   		 legend.position="none")
prev_plot

R0_biased = R0_all[time >= tval][kind %like% "Biased"]
sp1 = R0_biased[variable == "Species 1"]
sp2 = R0_biased[variable == "Species 2"]
bias_values = sp1$value / sp2$value

R0_true = R0_all[time >= tval][kind %like% "True"]
sp1 = R0_true[variable == "Species 1"]
sp2 = R0_true[variable == "Species 2"]
true_values = sp1$value / sp2$value


R0_plot = ggplot(R0_all[time >= tval]) + geom_line(aes(x=time - tval, y=value, color=variable, linetype=kind)) + 
							  theme_classic() + ylab("Maintenace potential\\n(Instantaneous R0)") + 
							  xlab("Time (arbitrary units)") + 
							  scale_color_discrete(name="") +
							  scale_linetype_discrete(name="R0 estimate") 

R0_plot

full_plot = prev_plot / R0_plot + plot_annotation(tag_levels = 'A', tag_suffix=".")
full_plot

ggsave("../results/temporal_context_example.pdf", width=7, height=7)






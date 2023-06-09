###############################################################
#                                                             #
# Possibilistic Network.                                      #
# Social Polarization in the Metropolitan Area of Marseille.  #
# An R script to build and execute a possibilistic network.   #
#                                                             #
# AuthorS: Andrea G. B. Tettamanzi                            #
#          Université Côte d'Azur, CNRS, I3S, UMR7271         #
#          06903 Sophia Antipolis, France                     #
#          E-mail: andrea.tettamanzi@unice.fr                 #
#                                                             #
#          Giovanni Fusco                                     #
#          Université Côte d'Azur, CNRS, ESPACE, UMR7300      #
#          98 Bd Herriot, BP 3209, 06204 Nice, France         #
#          E-mail: giovanni.fusco@unice.fr                    #
#                                                             #
###############################################################

# This script builds and executes a possibilistic network which has the same structure as the Baysesian network 
# proposed by F. Scarella [2] in order to model social polarization in the metropolitan area of Marseille (France).
# The model can be used to infer a trend scenario of social polarization of the 439 municipalities in the metropolitan
# area of Marseille (France) in a 10 years time.
# Within the model, a valorized municipality is defined as a municipality where executives and professionals are
# overrepresented within its resident population; a devalorized municipality is defined as a municipality where the
# unemployed are overrepresented within its resident population.
# The initial data for the study area were elaborated for 2009 and are in the Data_Marseille_2009.txt file.
# Other auxiliary files are:
#   - DependentVariables.txt containing the dependency structure of the possibilistic network
#   - VariableModalities.txt containing the values of each variable of the network
#   - VariableModalities_PlainEnglish.txt gives plain English names to variables and modalities, but is not used by 
#     the R-script.
#   - ModelStructure.png visualizes the DAG structure of the possibilistic network
#
# The network is build using uncertain logical gates ([3], [4]), which are the possibilistic counterpart of
# the noisy logical gates used in Bayesian networks [2]. More thorough presentations of the model and of the model 
# results are available in [5] and in [6].
# Model results and comparison with Bayesian network results can be explored through an interactive data-visualization 
# at the following address: https://public.tableau.com/profile/fusco#!/vizhome/RepresentingUncertainFutures/Story1

# Acknowledgement:
# This script was produced within the Géo-Incertitude project (2014-2015, CNRS grant of the PEPS HuMaIn program). 

# References:
# [1] Francisco Díez and Marek Druzdzel. "Canonical Probabilistic Models for Knowledge Engineering",
#     Tech. Rep. CISIAN-06-01, version 0.9, April 28, 2007.
# [2] Floriane Scarella, La ségrégation résidentielle dans l'espace-temps métropolitain: 
#     analyse spatiale et géo-prospective des dynamiques résidentielles de la métropole azuréenne, 
#     PhD dissertation, University of Nice Sophia Antipolis, 2014.
# [3] Matteo Caglioni, Didier Dubois, Giovanni Fusco, Diego Moreno, Henri Prade, Floriane Scarella,
#     and Andrea Tettamanzi. "Mise en oeuvre pratique de réseaux possibilistes pour modéliser
#     la spécialisation sociale dans les espaces métropolisés", LFA 2014 - Cargèse 22-24 novembre 2014, 
#     Cépaduès, Toulouse, ISBN : 9782364931565,  pp. 267-274
# [4] Didier Dubois, Giovanni Fusco, Henri Prade, and Andrea Tettamanzi, "Uncertain Logical Gates in Possibilistic Networks.
#     An Application to Human Geography". In Ch. Beierle and A. Dekhtyar (Eds.). Scalable Uncertainty Management - 
#     9th International Conference, SUM 2015, Québec City, QC, Canada, September 16-18, 2015. Proceedings (ISBN: 978-3-319-23539-4), 
#     Lecture Notes in Artificial Intelligence, vol. 9310, Springer, pp. 249-263.
# [5] Didier Dubois, Giovanni Fusco, Henri Prade, and Andrea Tettamanzi, "Uncertain Logical Gates in Possibilistic Networks:
#     Theory and application to human geography", International Journal of Approximate Reasoning, 2016 (in progress)
# [6] Giovanni Fusco, Cristina Cao, Didier Dubois, Henri Prade, Floriane Scarella, and Andrea Tettamanzi, 2015, 
#     Social polarization in the metropolitan area of Marseille. Modelling uncertain knowledge with probabilistic and 
#     possibilistic networks, ECTQG 2015 - XIX European Colloquium on Theoretical and Quantitative Geography, Bari (Italy),
#     September 3rd-7th 2015, Proceedings, Plurimondi. An International Forum for Research and Debate on Human Settlements, 8 p.

# To install the graph, RBGL and Rgraphviz libraries if not previously installed :
# source("http://bioconductor.org/biocLite.R")
# biocLite("graph")
# biocLite("RBGL")
# biocLite("Rgraphviz")

library(graph)
library(gRbase)

############################################################
# Part 1. General Functions and Utilities                  #
############################################################

# Set the CPT for variable v with conditional possibilities cp.
# Arguments:
# v - a string, the name of the variable.
# cp - an array of conditional possibilities, which will be used
#      to fill the table (possibly by recycling up to the desired length).
#      The leftmost indices (modalities of the parent variables) move fastest.
#      By convention, v is the leftmost (first) variable.
#
set.cpt <- function(v, cp)
{
	p <- parents(v, g)
	d <- length(poss.distr[[v]])
	n <- list()
	n[[v]] <- names(poss.distr[[v]])
	for(i in p)
	{
		d <- c(d, length(poss.distr[[i]]))
		n[[i]] <- names(poss.distr[[i]])
	}
	cpt[[v]] <<- array(rep_len(cp, prod(d))) # recycle cp to obtain the right length
	dim(cpt[[v]]) <<- d
	# print(c(v, p))
	dimnames(cpt[[v]]) <<- n
}

# Retrieve Possibility(v = modality | combination of values of parent variables)
# Arguments:
# v           - a string, the name of the variable
# modality    - a string, the name of one of the modalities of v
# combination - an array of strings, the names of the modalities of the parent variables
#
get.cp <- function(v, modality, combination)
{
	cmd <- sprintf("cpt$%s[\\"%s\\"", v, modality)
	for(parent.modality in combination)
		cmd <- sprintf("%s, \\"%s\\"", cmd, parent.modality)
	cmd <- paste(cmd, "]", sep = "")
	eval(parse(text = cmd))
}

# Set Possibility(v = modality | combination of values of parent variables)
# Arguments:
# v           - a string, the name of the variable
# modality    - a string, the name of one of the modalities of v
# combination - an array of strings, the names of the modalities of the parent variables
# poss        - the new value to be set
#
set.cp <- function(v, modality, combination, poss)
{
	cmd <- sprintf("cpt$%s[\\"%s\\"", v, modality)
	for(parent.modality in combination)
		cmd <- sprintf("%s, \\"%s\\"", cmd, parent.modality)
	cmd <- paste(cmd, "] <<- ", poss, sep = "")
	eval(parse(text = cmd))
}

# Generate a conditional possibility table with an uncertain MAX 
# Arguments:
# v   - a string, the name of the variable.
# prm - the parameters, organized as a list of records where each record is composed of:
#       prm[[i]]$pv    - the name of a parent variable (NULL corresponds to TRUE condition, i.e., leak)
#       prm[[i]]$mod   - a modality of the parent variables: if pv == mod, this record is used;
#       prm[[i]]$k     - an array of normalized possibilities, one for each modality of v.
# N.B.: the leak factors are expressed by a record with NULL pv (i.e., default condition)
#

uncertain.MAX <- function(v, prm)
{
	# Initialize the table of conditional possibilities to zero:
	set.cpt(v, 0)
	# Get the modalities of v and its parent variables:
	modalities <- names(poss.distr[[v]])
	p <- parents(v, g)
	cat("Constructing the possibility distribution of", v, "based on the uncertain MAX...\\n")
	# Let l be a list of vectors, each vector containing the modalities of a parent variable:
	l <- list()
	for(i in p)
		l[[i]] <- names(poss.distr[[i]])
	# We use expand.grid() to create a dataframe of all possible combinations of the
	# modalities of the parent variables:
	combinations <- expand.grid(l)
	for(c in 1:dim(combinations)[1])
	{
		# x will be a table holding the values of the parent variables for this combination
		x <- combinations[c, ]
		# Extract from prm the normalized possibilities whose conditions are satisfied by x...
		paramlist <- list() # This will be the list of "active" rows of prm
		for(row in prm)
			if(is.null(row$pv) || Reduce("&&", x[row$pv]==row$mod))
				paramlist[[length(paramlist) + 1]] <- row$k
		# For each combination of the parameters in the extracted rows...
		# (we use expand.grid() to create a dataframe of combinations
		# of indices in {1, ..., length(modalities)})
		args = ""
		for(r in 1:length(paramlist))
		{
			if(nchar(args) > 0)
				args <- paste(args, ", ", sep = "")
			args <- paste(args, "seq(1, ", length(modalities), ")", sep = "")
		}
		ind <- eval(parse(text = paste("expand.grid(", args, ")", sep = "")))
		# Finally, compute the conditional possibilities for v given x:
		for(i in 1:dim(ind)[1])
		{
			# Compute the minimum of the parameters, which we will call "beta":
			beta <- 1
			for(r in 1:length(paramlist))
				beta <- min(beta, paramlist[[r]][ind[i, r]])
			# Determine which modality this beta belongs to:
			# (here, we make the assumption that v is ordinal, and
			# the modalities of v are in increasing order of "priority",
			# i.e., the lowest-priority first and the highest-priority last)
			j <- max(ind[i, ])
			if(get.cp(v, modalities[j], x) < beta)
				set.cp(v, modalities[j], x, beta)
		}
	}
}

# Generate a conditional possibility table with an uncertain MAX-threshold 
# Arguments:
# v   - a string, the name of the variable.
# prm - the parameters, organized as a list of records where each record is composed of:
#       prm[[i]]$pv    - the name of a parent variable (NULL corresponds to TRUE condition, i.e., leak)
#       prm[[i]]$mod   - a modality of the parent variables: if pv == mod, this record is used;
#       prm[[i]]$k     - an array of normalized possibilities, one for each modality of v.
# thr   a vector of thresholds, one for each modality of v, which represent
#       the minimal number of conditions for which k > 0 (i.e., that modality of v is not impossible).
# N.B.: the leak factors are expressed by a record with NULL pv (i.e., default condition)
#

uncertain.MAX.threshold <- function(v, prm, thr)
{
	# Initialize the table of conditional possibilities to zero:
	set.cpt(v, 0)
	# Get the modalities of v and its parent variables:
	modalities <- names(poss.distr[[v]])
	p <- parents(v, g)
	cat("Constructing the possibility distribution of", v, "based on the uncertain MAX-threshold...\\n")
	# Let l be a list of vectors, each vector containing the modalities of a parent variable:
	l <- list()
	for(i in p)
		l[[i]] <- names(poss.distr[[i]])
	# Initialize a vector of leak coefficients, to use to count the number of possible combinations:
	leak <- rep(0, length(modalities))
	for(row in prm)
		if(is.null(row$pv))
			leak <- pmax(leak, row$k)
	# We use expand.grid() to create a dataframe of all possible combinations of the
	# modalities of the parent variables:
	combinations <- expand.grid(l)
	for(c in 1:dim(combinations)[1])
	{
		# x will be a table holding the values of the parent variables for this combination
		x <- combinations[c, ]
		count <- rep(0, length(modalities)) # number of combinations which are not impossible for each modality 
		# Extract from prm the normalized possibilities whose conditions are satisfied by x...
		paramlist <- list() # This will be the list of "active" rows of prm
		for(row in prm)
			if(is.null(row$pv) || Reduce("&&", x[row$pv]==row$mod))
				paramlist[[length(paramlist) + 1]] <- row$k
		# For each combination of the parameters in the extracted rows...
		# (we use expand.grid() to create a dataframe of combinations
		# of indices in {1, ..., length(modalities)})
		args = ""
		for(r in 1:length(paramlist))
		{
			if(nchar(args) > 0)
				args <- paste(args, ", ", sep = "")
			args <- paste(args, "seq(1, ", length(modalities), ")", sep = "")
		}
		ind <- eval(parse(text = paste("expand.grid(", args, ")", sep = "")))
		# Finally, compute the conditional possibilities for v given x:
		for(i in 1:dim(ind)[1])
		{
			# Compute the minimum of the parameters, which we will call "beta":
			beta <- 1
			for(r in 1:length(paramlist))
				beta <- min(beta, paramlist[[r]][ind[i, r]])
			# Determine which modality this beta belongs to:
			# (here, we make the assumption that v is ordinal, and
			# the modalities of v are in increasing order of "priority",
			# i.e., the lowest-priority first and the highest-priority last)
			j <- max(ind[i, ])
			# Update the count for checking the threshold:
			if(beta > leak[j])
				count[j] <- count[j] + 1
			if(count[j] >= thr[j])
				beta <- 1
			# Modify the CPT:
			if(get.cp(v, modalities[j], x) < beta)
				set.cp(v, modalities[j], x, beta)
		}
	}
}


# Generate a conditional possibility table with a uncertain OR.
# The uncertain OR is a special case of a uncertain MAX, where all the variables involved
# are binary (i.e., they have two modalities).
# Arguments:
# v    - a string, the name of the variable.
# prm - the parameters, organized as a list of records where each record is composed of:
#       prm[[i]]$pv    - the name of a parent variable (NULL corresponds to TRUE condition, i.e., leak)
#       prm[[i]]$mod   - a modality of the parent variables: if pv == mod, this record is used;
#       prm[[i]]$k     - the possibility that the parent variable is inhibited (or the leak if pv==NULL).
# N.B.: the leak factors are expressed by a record with NULL pv (i.e., default condition)
#
uncertain.OR <- function(v, prm)
{
	p <- parents(v, g)
	# translate the parameters into the format accepted by uncertain.MAX():
	params <- list()
	for(i in 1:(length(prm)-1))
	{
    params[[i]] <- list(
      pv = prm[[i]]$pv,
      mod = prm[[i]]$mod,
      k = c(prm[[i]]$k, 1.0)
      )
	}
    
   	params[[length(prm)]] <- list(
			pv = prm[[length(prm)]]$pv,
			mod = prm[[length(prm)]]$mod,
			k = c(1.0, prm[[length(prm)]]$k)
		)
	
	# call the uncertain MAX with the appropriate parameters
	  uncertain.MAX(v, params)
}


# Generate a conditional possibility table with a uncertain AND
# Arguments:
# v    - a string, the name of the variable.
# prm - the parameters, organized as a list of records where each record is composed of:
#       prm[[i]]$pv    - the name of a parent variable (NULL corresponds to TRUE condition, i.e., leak)
#       prm[[i]]$mod   - a modality of the parent variables: if pv == mod, this record is used;
#       prm[[i]]$k     - the possibility that the parent variable is inhibited (or the leak if pv==NULL).
# N.B.: the leak factors are expressed by a record with NULL pv (i.e., default condition)
#
uncertain.AND <- function(v, prm)
{
	# Initialize the table of conditional possibilities to one:
	set.cpt(v, 1)
	# Get the modalities of v and its parent variables:
	modalities <- names(poss.distr[[v]])
	p <- parents(v, g)
	cat("Constructing the possibility distribution of", v, "based on the uncertain AND...\\n")
	# Let l be a list of vectors, each vector containing the modalities of a parent variable:
	l <- list()
	for(i in p)
		l[[i]] <- names(poss.distr[[i]])
	# We use expand.grid() to create a dataframe of all possible combinations of the
	# modalities of the parent variables:
	combinations <- expand.grid(l)
	for(c in 1:dim(combinations)[1])
	{
		# x will be a vector holding a single combination of modalities
		x <- combinations[c, ]
		beta <- 0
		leak <- 1
		default <- FALSE
		
    # Compute from prm the possibilities based on the conditions satisfied by x...

    for(row in prm)
		{
		  if(is.null(row$pv))
		    leak <- row$k
		  else if(Reduce("&&", x[row$pv] == row$mod))
        beta <- max(beta, row$k)
		  else default <- TRUE # at least one of the conditions is not satisfied!
		}
		
    
    # (here, we make the assumption that v is binary, and
		# the modalities of v are in increasing order of "severity",
		# i.e., the "absence" first and the "presence" last)
		if(default)
		{
			set.cp(v, modalities[1], x, 1)
			set.cp(v, modalities[2], x, leak)
		}
		else
		{
			set.cp(v, modalities[1], x, beta)
			set.cp(v, modalities[2], x, 1)
		}
	}
}



# Show the CPT for variable v in a human-friendly format
#
show.cpt <- function(v)
{
	p <- parents(v, g)
	# cat("Computing the possibility distribution of", v, "...\\n")
	# Let l be a list of vectors, each vector containing the modalities of a parent variable:
	l <- list()
	for(i in p)
		l[[i]] <- names(poss.distr[[i]])
	# We use expand.grid() to create a dataframe of all possible combinations of the
	# modalities of the parent variables:
	combinations <- expand.grid(l)
	for(modality in names(poss.distr[[v]]))
	{
		for(i in 1:dim(combinations)[1])
		{
			cat(sprintf("Poss(%s = %s | ", v, modality))
			# Retrieve Possibility(v = modality | i)
			cmd <- sprintf("cpt$%s[\\"%s\\"", v, modality)
			for(j in 1:dim(combinations)[2])
			{
				if(j>1) cat(", ")
				cat(sprintf("%s = %s", names(combinations)[j], combinations[i, j]))
				cmd <- sprintf("%s, \\"%s\\"", cmd, combinations[i, j])
			}
			cat(") = ")
			cmd <- paste(cmd, "]")
			cond.poss <- eval(parse(text = cmd))
			cat(cond.poss, "\\n")
		}
	}
}

# Compute the possibility distribution for variable v
#
compute.poss <- function(v)
{
	p <- parents(v, g)
	# cat("Computing the possibility distribution of", v, "...\\n")
	# Let l be a list of vectors, each vector containing the modalities of a parent variable:
	l <- list()
	for(i in p)
		l[[i]] <- names(poss.distr[[i]])
	# We use expand.grid() to create a dataframe of all possible combinations of the
	# modalities of the parent variables:
	combinations <- expand.grid(l)
	for(modality in names(poss.distr[[v]]))
	{
		max.poss <- 0.0
		# for each combination of the modalities of the parent variables,
		# compute the minimum of the possibility of that combination of modalities and
		# of the conditional possibility of v = modality given that combination;
		# take the maximum of those minima.
		for(i in 1:dim(combinations)[1])
		{
			min.poss <- 1.0
			for(x in p)
				min.poss <- min(min.poss, poss.distr[[x]][combinations[i, x]])
			min.poss <- min(min.poss, get.cp(v, modality, combinations[i, ]))
			max.poss <- max(max.poss, min.poss)
		}
		poss.distr[[v]][modality] <<- max.poss
	}
}

# Propagate possibilities forward
#
propagate <- function()
{
	# Get a topologically sorted list of the variables:
	sorted.vars <- topoSort(g)
	for(v in sorted.vars)
	{
		p <- parents(v, g)
		if(length(p) > 0)
			compute.poss(v)
	}
}

############################################################
# Part 2. Definition of the Possibilistic Network          #
############################################################

# creates a dag from a txt file where each columnt lists a dependent variable followed by its parents
var.df <- read.table ("DependentVariables.txt", sep="\\t", na.strings = "", header = FALSE)

# the text calling the dag function (dag.txt) has to be written iteratively by parsing the var.df dataframe
dag.txt <- "dag("

for (i in 1:dim(var.df)[2])
  {
    dag.txt <- paste (dag.txt, "c(", sep = "")
    #show (var.df[[i]])
    for (j  in 1:dim(var.df)[1]) 
     {
       if( is.na(var.df[[i]][j]) == FALSE)
          dag.txt <- paste(dag.txt, "\\"", var.df[[i]][j],"\\"",",",  sep = "")
     }
    dag.txt <- paste(dag.txt, "),",  sep = "")
  }

dag.txt <- paste (dag.txt, ")", sep ="")
dag.txt <- gsub (",)",")", dag.txt)

# the dag function is called by evaluating dag.txt
g <- eval(parse(text=dag.txt))
show (g)

# Prepare a list of possibility distributions, one per variable:
poss.distr <- list()
length(poss.distr) <- length(nodes(g))
names(poss.distr) <- nodes(g)

# Prepare a list of conditional probability matrices, one per variable:
cpt <- list()
length(cpt) <- length(nodes(g))
names(cpt) <- nodes(g)

############################################################
# Part 3. Initialization of the distributions and CPTs     #
#                                                          #
# N.B.: modalities of all variables MUST be ordered        #
#       by increasing "severity" / "absence" < "presence"  #
#                                                          #
############################################################

# Reads text file listing the modalities of each variable
var.df <- read.table("VariableModalities.txt", sep="\\t", colClasses = "character", na.strings = "", header = TRUE)


# Checks that the names of Variables in the file corresponds to the names of variables in the dag
a <- sort(names(poss.distr))
b <- sort(names(var.df))
if (identical(a,b))
  {
    print("Variable names consistent with variable names in graph")
  } else  {
    print("WARNING: Variable names not consistent with variable names in graph!")
  }
rm (a, b)  

# Initializes complete-ignorance possibility distributions for all the variables:
for (v in names(poss.distr))
{
  poss.distr[[v]] <- 1.0
  b <-  var.df[[v]][1]
 
  for (j  in 2:dim(var.df)[1]) 
  {
    if( is.na(var.df[[v]][j]) == FALSE)
      {
        poss.distr[[v]] <- c(poss.distr[[v]], 1.0)
        b <- c (b, var.df[[v]][j]) 
      }
  }

names(poss.distr[[v]]) <- b
rm (b)
#show(poss.distr[[v]])
}
rm (var.df)

# Set up the Conditional Possibility Tables:

set.cpt("SituationT2", c(
# SituationT2
# valorisé, autre, dévalorisé| SituationT1,   Évolution
	1,      0,     0,       # |  valorisé,         aucune
	0,      1,     1,       # |  valorisé, dévalorisation
  1,      0,     0,       # |  valorisé,   valorisation
	0,      1,     0,       # |     autre,         aucune
	0,      0,     1,       # |     autre, dévalorisation
	1,      0,     0,       # |     autre,   valorisation
	0,      0,     1,       # |dévalorisé,         aucune
	0,      0,     1,       # |dévalorisé, dévalorisation
	1,      1,     0        # |dévalorisé,   valorisation
))

# parc_attractif is a intermediate deterministic variable in order to calculate possibilities for atout_gentrification
set.cpt("parc_attractif", c(
  # parc_attractif
  # no,    yes      | bati_maitrisé, diversification_logements, espace_temps
  0,      1,     # |           yes,                       yes,       centre
  0,      1,      # |           yes,                        no,       centre
  0,      1,      # |           yes,                       yes,       couronne1
  0,      1,      # |           yes,                        no,       couronne1
  0,      1,      # |           yes,                       yes,       couronne2
  0,      1,      # |           yes,                        no,       couronne2
  0,      1,      # |            no,                       yes,       centre
  1,      0,      # |            no,                        no,       centre
  0,      1,      # |            no,                       yes,       couronne1
  1,      0,      # |            no,                        no,       couronne1
  1,      0,      # |            no,                       yes,       couronne2
  1,      0       # |            no,                        no,       couronne2
))

# prm - the parameters, organized as a list of records where each record is composed of:
#       prm[[i]]$pv    - the name of a parent variable (NULL corresponds to TRUE condition, i.e., leak)
#       prm[[i]]$mod   - a modality of the parent variables: if pv == mod, this record is used;
#       prm[[i]]$k     - an array of normalized possibilities, one for each modality of v.
# N.B.: the leak factors are expressed by a record with NULL pv (i.e., default condition)
uncertain.MAX.threshold ("Évolution", list(
  list(pv = "Stationnarité", mod = "no",      k = c(1, 0.3, 0.2)),
  list(pv = "PosVoisinage",  mod = "marge",   k = c(1, 0.3, 0.2)),
  list(pv = "PosVoisinage",  mod = "enclave", k = c(1, 0.3, 0.2)),
  list(pv = "PercVal",       mod = "yes",     k = c(1, 0, 0.4)),
  list(pv = "PercDéval",     mod = "yes",     k = c(1, 0.4, 0)),
  list(pv = "atout_gentrification", mod = "yes",     k = c(1, 0, 0.4)),
  list(pv = "frein_gentrification",mod = "yes",     k = c(1, 0.4, 0)),
  list(pv = NULL,            mod = NULL,      k = c(1, 0.1, 0.1))  # <---- LEAK
), c(6, 3, 3)
)

# diversification_logements needs an uncertain.MAX function because of the very weak possibilistic force of type_constr_neuve ="collectif_sup50"
uncertain.MAX("diversification_logements", list(
  list(pv = c("type_constr_neuve","individuel_sup80"),     mod = c("collectif_sup50","no"),   k = c(.2, 1)),
  list(pv = c("type_constr_neuve","individuel_sup80"),     mod = c("ind_majoritaire","no"), k = c(.7, 1)),
  list(pv = "type_constr_neuve",     mod = "collectif_sup50", k = c(1, .4)),
  list(pv = NULL,           mod = NULL,       k = c(1, .2))  # <---- LEAK
))


uncertain.AND("PercVal", list(
	list(pv = "ImpactFD",     mod = "cadres",   k = .6),
	list(pv = "SitDomVoisinage", mod = "valorisé", k = .2),
	list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))

uncertain.AND("PercDéval", list(
  list(pv = "ImpactFD",     mod = "chômeurs",   k = .2),
  list(pv = "SitDomVoisinage", mod = "dévalorisé", k = .6),
  list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))

uncertain.AND("atout_gentrification", list(
  list(pv = "aménités_env",     mod = "yes",   k = .6),
  list(pv = "parc_attractif",     mod = "yes", k = .6),
  list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))

uncertain.OR("aménités_env", list(
  list(pv = "espaces_naturels",     mod = "yes",   k = .4),
  list(pv = "agricole_valorisé",     mod = "yes", k = .2),
  list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))

uncertain.AND("bati_maitrisé", list(
  list(pv = "bati_diffus",     mod = "no",   k = .7),
  list(pv = "pression_constr_neuve",     mod = "no", k = .4),
  list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))


uncertain.OR("frein_gentrification", list(
  list(pv = c("diversification_logements","parc_peu_attract"),     mod = c("no","yes"),   k = .2),
  list(pv = c("parc_peu_attract","chomage_problématique"),     mod = c("yes","yes"), k = .4),
  list(pv = c("diversification_logements","chomage_problématique"),     mod = c("no","yes"), k = .6),
  list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))

uncertain.OR("parc_peu_attract", list(
  list(pv = "espace_temps",     mod = "couronne2",   k = .6),
  list(pv = "parc_ancien_préponderant",     mod = "yes",   k = .4),
  list(pv = "RS_LV_préponderants",     mod = "yes",   k = .4),
  list(pv = "HLM_important",     mod = "yes",   k = .4),
  list(pv = NULL,           mod = NULL,       k = .2)  # <---- LEAK
))

############################################################
# Part 4. Batch Execution                                  #
############################################################

# Execute the possibilistic network in batch mode
# on the dataset contained in the given file
#
batch <- function(filename = "Data_Marseille_2009.txt")
{
	res.df <<- read.table(filename, header = TRUE, na.strings = "*")
	Évolution.list <- list()
	SitutationT2.list <- list()
	for(modality in names(poss.distr$Évolution))
	  Évolution.list[[paste("Évolution", modality, sep = ".")]] <- numeric()
	for(modality in names(poss.distr$SituationT2))
	  SitutationT2.list[[paste("SituationT2", modality, sep = ".")]] <- numeric()
	for(row in 1:dim(res.df)[1])
	{
		# 1. Set the possibility distributions for input variables:
		for(v in nodes(g))
		{
			p <- parents(v, g)
			if(length(p) == 0)
			{
			  if (is.na(res.df[row,v]))
			  {
			    #sprintf ("NA:", res.df[[row]][1])
			    
			    #show (paste (row, v, "NA",res.df[row,v], sep = " : "))
			    for(modality in names(poss.distr[[v]]))
			      poss.distr[[v]][modality] <<- 1.0
			    #show (poss.distr[[v]])
			  } else  {
				for(modality in names(poss.distr[[v]]))
					poss.distr[[v]][modality] <<- 0.0
				poss.distr[[v]][as.character(res.df[row, v])] <<- 1.0
				#how (poss.distr[[v]])
			  }
			}
		}

		# 2. Execute the network:
		propagate()

		# 3. Write the result (for the time being, on the console):
		cat("IDnum", res.df$IDnum[row], ":\\n")
		print(poss.distr$SituationT2)
		for(modality in names(poss.distr$Évolution))
		  Évolution.list[[paste("Évolution",modality, sep = ".")]][row] <- poss.distr$Évolution[modality]
    for(modality in names(poss.distr$SituationT2))
			SitutationT2.list[[paste("SituationT2", modality, sep = ".")]][row] <- poss.distr$SituationT2[modality]
	}
	res.df <<- cbind(res.df, Évolution.list, SitutationT2.list)
  write.csv(res.df, "results10-01th3.csv")
}


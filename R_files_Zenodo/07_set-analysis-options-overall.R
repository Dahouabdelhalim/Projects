# This block of code allows the user to control the way the data are analysed

# Specify whether Other is to be included or excluded
exclude.other.drivers <- "Yes" # options: "No", "Yes"
# --> "No": Other drivers are included in the comparisons and weights
# --> "Yes": Comparisons involving Other drivers are removed so that only the "Big 5" are considered

# Specify the hypothesis made to assign a relative impact to drivers that are not analysed in a study 
include.non.assessed.drivers <- "Less" # options: "No", "Less", "Equal", "Unknown" 
# --> "No": non-assessed drivers have NO impact (hypothesis 1)
# --> "Less": non-assessed drivers have LESS impact than lowest-ranked assess driver (hypothesis 2)
# --> "Equal": non-assessed drivers have SAME impact as lowest-ranked assessed driver (hypothesis 2bis)
# --> "Unknown": non-assessed drivers have UNKNOW impact - i.e. average impact of assessed drivers (hypothesis 3)

# Specify whether magnitude of impact should be used instead of ranks to estimate the relative impact of drivers for the studies where this information is available
combine.rank.magnitude <- "No" # options:  "Yes", "No"
# --> "Yes": use magnitude instead of rank information when magnitude is reliably estimated
# --> "No": use only rank information for all studies even if magnitude is reliably estimated (we decided to go for this option)

# Specify the level at which the relative impacts of drivers are to be assessed
analysis.level <- "overall" # options: "indicators", "ebv", "indicators.x.ebv", "overall"
# --> "indicators": assessment of relative impact of drivers at the level of each individual indicators (i.e. aggregated across studies within each indicator)
# --> "ebv": assessment of relative impact of drivers at the level of each ebv class (i.e. aggregated across the individual indicators within each ebv class)
# --> "indicators.x.ebv": assessment of relative impact of drivers at the level of each individual indicators  for the ebv level separately
# --> "overall": assessment of relative impact of drivers at the overall level (i.e. aggregated across all indicators)

# Specify whether the analysis should be restricted to a selected list of indicators 
indicator.selection <- "No" # options: "Yes", "No"
# --> "No": analysis based on all indicators included in the data files
# --> "Yes": analysis based on the indicators selected below
selected.indicator <- c("FSS","IFL","IUC","LCC","LPI","LSR","MFC","MLF","MSA","NPP","RLI") # acronyms of the indicators selected for the analysis

# Specify whether the analysis should be excluded some indicators 
indicator.exclusion <- "No" # options: "Yes", "No"
# --> "No": analysis based on all indicators included in the data files
# --> "Yes": analysis based on all indicators except those selected below
excluded.indicator <- c("") # acronyms of the indicators excluded from the analysis

# Specify whether the analysis should be restricted to a particular realm 
realm.selection <- "No" # options: "Yes", "No"
# --> "No": analysis based on all realms
# --> "Yes": analysis based on the realm selected below
selected.realm <- "Marine" # options: "Terrestrial", "Freshwater", "Marine", "All realms" ("All realms" = studies with assessments across the three realms)
include.all.realms <- "No" # options: "Yes", "No"
# --> "No": exclude results from studies with cross-realm assessments ("All realms")
# --> "Yes": include results from studies with cross-realm assessments ("All realms") along with studies with assessments on the selected realm

# Specify whether the analysis should be restricted to a particular IPBES region 
region.selection <- "No" # options: "Yes", "No"
# --> "No": analysis based on all IPBES regions
# --> "Yes": analysis based on the IPBES region selected below
selected.region <- "Europe and Central Asia" # options: "Americas", "Africa", "Europe and Central Asia", "Asia-Pacific", "All regions" ("All regions" = studies with assessments across the four regions)
include.all.regions <- "No" # options: "Yes", "No"
# --> "No": exclude results from studies with cross-region assessments ("All regions")
# --> "Yes": include results from studies with cross-region assessments ("All regions") along with studies with assessments on the selected region

# Specify whether the analysis should be restricted to a particular climatic domain 
domain.selection <- "No" # options: "Yes", "No"
# --> "No": analysis based on all domains
# --> "Yes": analysis based on the domain selected below
selected.domain <- "Polar" # options: "Polar", "Boreal", "Temperate", "Subtropical", "Tropical", "All domains" ("All domains" = studies with assessments across the five domains)
include.all.domains <- "No" # options: "Yes", "No"
# --> "Yes": exclude results from studies with cross-domain assessments ("All domains")
# --> "Yes": include results from studies with cross-domain assessments ("All domains") along with studies with assessments on the selected domain

# Specify whether the analysis should be restricted to studies covering (a) particular spatial scale(s)  
scale.selection <- "No" # options: "Yes", "No"
# --> "No": analysis based on all spatial scales*
# --> "Yes": analysis based on the spatial scale(s) selected below
# * even if all spatial scales are included in the analysis, a weighting procedure is available below so that the weight of studies in the analysis
#   may vary according to their spatial coverage 
selected.scale <- c("local","regional","continental","global") # options: "local", "regional", "continental", "global"

# Specify whether the analysis should be restricted to studies analysing the impact of a specific number of drivers
n.drivers.selection <- "No" # options: "Yes", "No"
# --> "No": analysis based on studies assessing the relative impacts of 2 to 6 drivers
# --> "Yes": analysis based on studies assessing the relative impacts of the number of drivers specified below
# * even if all studies assessing the relative impacts of varying drivers are included in the analysis, a weighting procedure is available below so that 
#   the weight of studies in the analysis may vary according to the number of drivers they assessed 
selected.n.drivers <- c(2,3,4,5,6) # options: number of drivers whose impacts are assessed for a study to be included in the analysis

# Specify whether the analysis should exclude information as define in a specific column of the dataset to avoid redundant information among/within studies
# (recommended option)
manual.selection <- "Yes" # options: "Yes", "No"
# --> "Yes": selection of studies or levels of analysis based on the column "Include"*
# --> "No": no selection of studies or levels of analysis based on the column "Include"*
# * "Include": column from the data files set to avoid including redundant pieces of information in the analysis (depending on the level of analysis)

# Specify whether the analysis should weight the importance of studies according to their spatial coverage and the number of drivers they assessed
weighted.mean.procedure <- "Yes" # options: "Yes", "No"
# --> "Yes": number of drivers assessed in the study and spatial coverage used as weights when averaging the importace of different drivers across studies
# --> "No": arithmetic mean when averaging the importace of different drivers across studies

# Specify how studies should be weighted according to their spatial scale; set all to 1 for equal weights (applies to analyses of both long and wide data)
scale.weights <- list("local"=1, "regional"=3, "continental"=7, "global"=15)

# Specify how studies should be weighted according to the number of drivers they study; set all to 1 for equal weights (applies only to long data)
ndriver.weights <- c(1, 3, 6, 10, 15, 21) #If set to triangular numbers, weights = number of pairwise comparisons made, as in analysis of wide data

# Specify whether indicator weights should be used to equalise the contribution of different indicators
indicator.reweights <- "sumscale" #If set to "sumscale", they are inversely proportional to the summed scale weight for each indicator; 
#                               if set to "none", evidence is not weighted by indicator rarity;
#                               if set to "comparisons", then weights are inversely proportional to the number of comparisons for each indicator;
#                               if set to "equalising", weights are inversely proportional to the sum of (scale weight * (Ndrivers - 1)) for each indicator.
#                               if set to "old_fix", they are inversely proportional to summed scale weight in original data AND ARE CONSTANT IN BOOTSTRAP REPLICATES

# Specify the 'score' given to each driver in the event of a draw (as a proportion of the score for a win)
draw.score <- 0.5 # options: 0, 0.5

# Specify whether to do a bootstrap analysis.level
do.bootstrap <- "Yes" # options: "Yes", "No"

# Specify number of bootstrap replicates (not used if bootstrap not done) and randomisations
nreps <- 10000 #options: Any positive integer; 10000 is sensible number when code is working properly, as it allows reasonably accurate estimation of p-values near 0.05
nrands <- 10000 #options: Any positive integer; 10000 is sensible number when code is working properly, as it allows reasonably accurate estimation of p-values near 0.05

# Specify whether the script should export results with estimation of relative impact of drivers 
export.tables <- "Yes" # options: "Yes", "No"
# --> "Yes": write csv tables with the results of aggregation and bootstrap
# --> "No": not implemented

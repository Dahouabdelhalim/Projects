# Appendix_AnalyzeData.R
# ======================
# Paper: A general approach for quantifying microbial effects on plant competition
# Authors: Po-Ju Ke and Joe Wan
# 
# Given simulated plant performance data, fit competitive responses and 
# predict the competitive outcome using the method outlined in the main text.

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(broom)

##### 1. Functions for performing invasion analysis #####
 
get.resident.eq <- function(c, b.con, a.con) {
  # Calculates the resident equilibrium from the resident's conspecific 
  # interactions.
  # 
  # Arguments:
  #   c: A numeric value or vector with the model intercept
  #   b.con: A numeric value or vector with the model's conspecific linear 
  #     coefficient (directly from the model)
  #   a.con: A numeric value or vector with the model's conspecific quadratic 
  #     coefficient (directly from the model; use 0 for a linear model)
  
  result <- case_when(
    # The formulas do not apply when the intercept is negative
    c <= 0 ~ as.numeric(NA),
    # Linear case: the equilibrium is simply -c / b, provided b is negative.
    a.con == 0 ~ if_else(b.con<0, -c/b.con, as.numeric(NA)),
    # Quadratic case: we take the smallest positive root from the quadratic 
    # formula; the conditions ensure that there is a real and positive root
    a.con < 0 | (b.con<0 & b.con^2>4*a.con*c) ~ 
      (-b.con - sqrt(b.con^2-4*a.con*c))/(2*a.con))
  return(result)
}

get.invasion.performance <- function(c, b.het, a.het, res.eq) {
  # Calculates the invasion performance from the invader's heterospecific 
  # interactions and the previously-calculated resident equilbrium.
  # 
  # Arguments:
  #   c: A numeric value or vector with the model intercept
  #   b.het: A numeric value or vector with the model's heterospecific linear 
  #     coefficient (directly from the model)
  #   a.het: A numeric value or vector with the model's heterospecific quadratic 
  #     coefficient (directly from the model; use 0 for a linear model)
  #   res.eq: A numeric value or vector with the previously-calculated resident 
  #     equilibrium.
  
  # We simply plug in the resident equilibrium into the quadratic predicting 
  # performance under heterospecific competition
  result <- c + b.het*res.eq + a.het*res.eq^2
  return(result)
}

get.maximum <- function(c, b, a) {
  # Calculates the maximum performance for a species' con- or heterospecific 
  # response.
  # 
  # Arguments:
  #   c: A numeric value or vector with the model intercept
  #   b: A numeric value or vector with the model's con- or heterospecific  
  #     linear coefficient (directly from the model)
  #   a: A numeric value or vector with the model's con- or heterospecific 
  #     quadratic coefficient (directly from the model; use 0 for a linear 
  #     model)
  result <- case_when(
    # The formulas do not apply when the intercept is negative
    c <= 0 ~ as.numeric(NA),
    # Linear case: if b is negative, the maximum is the intercept; otherwise
    # it is undefined
    a == 0 ~ if_else(b<=0, c, as.numeric(NA)),
    # Quadratic case: if response is concave down, the maximum is the vertex
    # of the parabola, or the intercept if the vertex is at a negative x-value
    a < 0 ~ if_else(b<=0, c, c - b^2/(4*a)),
    # If the response is concave up and there is a positive root, the maximum
    # is at the intercept; otherwise it is undefined
    a > 0 ~  if_else(b<0 & b^2>4*a*c, c, as.numeric(NA)))
  return(result)
}


##### 2. Process and fit simulated data #####

# Load data
# Ensure that the dataset 'Appendix_Data.csv' is in the working directory
appendix.data <- read_csv('Appendix_Data.csv')
appendix.data <- appendix.data %>%
  # Create separate predictors for conspecific and heterospecific density
  mutate(ConDensity=if_else(FocalSpecies==CompetitorSpecies, CompetitorDensity, 0),
         HetDensity=if_else(FocalSpecies==CompetitorSpecies, 0, CompetitorDensity)) %>%
  # Classify con- and heterospecific competition treatments
  mutate(CompetitionTreatment=if_else(
    CompetitorSpecies==FocalSpecies, 'Conspecific', 'Heterospecific'))

# Fit model for biomass of each species. Here, for generality, we fit a 
# quadratic model to each response, though in practice it may be best to examine
# the fits and only use a nonlinear fit when necessary.
{
  # Initialize results: a data frame with predicted values, a data frame for the
  # model summaries, and a list of the models
  predictions <- data.frame()
  summaries <- data.frame()
  models <- list()
  # Loop through species
  for (cur.species in c('Plant_i', 'Plant_j')) {
    # Subset the data to the current focal species
    cur.data <- filter(appendix.data, FocalSpecies==cur.species)
    
    # Fit a model to the focal species' performance
    cur.model <- lm(
      Performance ~ -1 + SoilTreatment +
        SoilTreatment:(ConDensity + I(ConDensity^2) + HetDensity + I(HetDensity^2)), 
      cur.data)
    # Add model to outputs
    models[[cur.species]] <- cur.model
    cur.summary <- summary(cur.model) %>% 
      tidy %>% mutate(FocalSpecies=cur.species) %>% select(FocalSpecies, everything())
    summaries <- bind_rows(summaries, cur.summary)
    
    # Generate input data for predictions
    cur.predictions <- with(cur.data, 
                            expand_grid(FocalSpecies=cur.species, 
                                        SoilTreatment=unique(SoilTreatment),
                                        CompetitorSpecies=unique(CompetitorSpecies), 
                                        CompetitorDensity=seq(0, 15, by=0.01))) %>% 
      mutate(ConDensity=if_else(FocalSpecies==CompetitorSpecies, CompetitorDensity, 0),
             HetDensity=if_else(FocalSpecies==CompetitorSpecies, 0, CompetitorDensity),
             CompetitionTreatment=if_else(
               CompetitorSpecies==FocalSpecies, 'Conspecific', 'Heterospecific')) %>%
      filter(SoilTreatment=='- Microbes' | str_sub(SoilTreatment, -1)==str_sub(CompetitorSpecies, -1))
    
    # Predict and add to combined output
    cur.predictions$Performance <- predict(cur.model, cur.predictions)
    predictions <- bind_rows(predictions, cur.predictions)
  }
}

# Process the combined summary tables to extract model terms
coefficients <- summaries %>%
  # Process the 'term' column to extract treatment and predictor
  mutate(
    SoilTreatment=str_extract(summaries$term, '(?<=SoilTreatment)[^:]+'),
    Term=str_extract(summaries$term, '(?<=:).+')) %>%
  # Find which coefficient the term corresponds to
  mutate(Coefficient=case_when(
    is.na(Term) ~ 'c', 
    Term=='ConDensity' ~ 'b.con', Term=='I(ConDensity^2)' ~ 'a.con',
    Term=='HetDensity' ~ 'b.het', Term=='I(HetDensity^2)' ~ 'a.het')) %>%
  # Get each predictor in its own column
  transmute(FocalSpecies, SoilTreatment, Coefficient, Value=estimate) %>%
  spread(Coefficient, Value)


##### 3. Perform invasion analysis #####

# Step 1. Calculate resident equilibrium for each plant/soil combination
resident.equilibria <- coefficients %>%
  mutate(ResidentEquilibrium=get.resident.eq(c, b.con, a.con)) %>%
  select(FocalSpecies, SoilTreatment, ResidentEquilibrium)
resident.equilibria

# Step 2. Calculate invasion performance
invasion.analysis <- coefficients %>% 
  # The focal species is the invader, and the other species is the resident
  mutate(InvaderSpecies=FocalSpecies,
         ResidentSpecies=if_else(FocalSpecies=='Plant_i', 'Plant_j', 'Plant_i')) %>%
  # Add a column with the resident equilibrium from Step 1
  left_join(resident.equilibria, 
            by=c('SoilTreatment', 'ResidentSpecies'='FocalSpecies')) %>%
  # Calculate the invasion performance
  mutate(InvasionPerformance=get.invasion.performance(c, b.het, a.het, ResidentEquilibrium)) %>%
  select(FocalSpecies, ResidentSpecies, SoilTreatment, ResidentEquilibrium, InvasionPerformance) %>%
  filter(!is.na(InvasionPerformance))
invasion.analysis

# Step 3. For MCT components, calculate reference performance
# Recommended general approach: *maximum* performance across all treatments
references <- coefficients %>%
  # Apply the get.maximum function to con- and heterospecific responses, then 
  # take the maximum of those (when defined)
  mutate(max.con=get.maximum(c, b.con, a.con), 
         max.het=get.maximum(c, b.het, a.het),
         max.all=pmax(max.con, max.het, na.rm=T),
         # For diagnostics, indicate whether maximum was from con- or hetero-
         max.source=if_else(replace_na(max.con,0)>replace_na(max.het,0), 
                            'Conspecific', 'Heterospecific')) %>%
  # Take the maximum value for the focal species across all treatments. In
  # more complex experimental designs, this must be done for each "plant--soil
  # system"
  group_by(FocalSpecies) %>% arrange(desc(max.all)) %>% slice(1) %>% ungroup %>%
  transmute(FocalSpecies, ReferencePerformance=max.all,
            Reference='Maximum',
            ReferenceSource=paste(SoilTreatment, max.source))
# Approach from Ke & Wan (2020): *intercept* from the *sterile* treatment.
# This is identical under the assumptions of Ke & Wan of primarily pathogenic 
# soils and no nonlinear responses, but may give undefined niche and fitness in
# cases where soils/competitors have positive effects
references <- coefficients %>%
  filter(SoilTreatment=='- Microbes') %>%
  transmute(FocalSpecies, ReferencePerformance=c,
            Reference='Sterile Intercept', ReferenceSource=SoilTreatment) %>%
  bind_rows(references, .)
references

# Step 4. Calculate relative reduction in performance, Carroll (2011)'s S, using
# invasion performance (IP) and reference performance (RP): S = (RP - IP)/RP
S.calcs <- invasion.analysis %>%
  left_join(references, by='FocalSpecies') %>%
  transmute(FocalSpecies,
            TreatmentType=if_else(SoilTreatment=='- Microbes', 'Sterile', 'Live'),
            IP=InvasionPerformance, Reference, RP=ReferencePerformance,
            S=(RP-IP)/RP, Reference, Persistence=S<1)
S.calcs

# Step 5. Calculate MCT components from S, and make predictions of coexistence
results <- 
  # Join calculations for plant i with those for plant j
  left_join(
    filter(S.calcs, FocalSpecies=='Plant_i'),
    filter(S.calcs, FocalSpecies=='Plant_j'),
    by=c('TreatmentType','Reference'), suffix=c('.i', '.j')) %>%
  # Calculate outcome of competition using invasion performance: positive 
  # invasion performance indicates that a plant persists
  mutate(Outcome=case_when(
    IP.i>0 & IP.j>0 ~ 'Coexistence',
    IP.i>0 & IP.j<=0 ~ 'Plant i wins',
    IP.i<=0 & IP.j>0 ~ 'Plant j wins',
    IP.i<=0 & IP.j<=0 ~ 'Priority effects')) %>%
  # Calculate MCT components: rho is the geometric mean of relative reductions,
  # while fitness ratio is the ratio of relative reduction and rho
  mutate(rho=sqrt(S.i*S.j), fi.fj=S.j/rho, fj.fi=S.i/rho) %>%
  # Select columns to keep
  select(TreatmentType, Reference, 
         IP.i, IP.j, RP.i, RP.j, Persistence.i, Persistence.j, 
         S.i, S.j, rho, fi.fj, fj.fi, Outcome) %>%
  arrange(Reference)
results

# Predictions: in sterile soil, plant i wins; in live soil, the species coexist.
# Predictions using invasion performance and using the niche/fitness components
# give the same answer; the two methods for calculating the S values give 
# slightly different value for niche/fitness but equivalent predictions.

##### 4. Visualize model fits and invasion calculations #####

# Visualize model fits: here, we plot each focal plant's responses in each soil,
# since each focal plant/soil combination shares an intercept.
# Note that not all combinations occur: each live treatment only uses the soil 
# of the competitor species.
ggplot(mapping=aes(x=CompetitorDensity, y=Performance, color=CompetitorSpecies)) +
  # Add solid lines for x- and y-axes
  geom_hline(yintercept=0, size=0.25) + geom_vline(xintercept=0, size=0.25) +
  # Add density/performance measurements
  geom_point(data=appendix.data) + 
  # Add model fits
  geom_line(data=filter(predictions, Performance>0)) + 
  # Use facets to show focal plant/soil combinations separately
  facet_grid(paste('Focal:', FocalSpecies) ~ SoilTreatment, scales='free') +
  theme_minimal()

# Visualize invasion analysis: here, we plot each combination of focal 
# plant/competition treatment, visualizing the effect of soil on the competitive 
# response.
# Note that each column of panels shares the same *competitor* species; thus,
# the bottom left panel shows the invasion performance of plant j, and the 
# bottom right panel shows the invasion performance of plant i.
ggplot(mapping=aes(x=CompetitorDensity, y=Performance, 
                   color=SoilTreatment)) +
  # Add solid lines for x- and y-axes
  geom_hline(yintercept=0, size=0.25) + geom_vline(xintercept=0, size=0.25) +
  # Add density/performance measurements
  geom_point(size=1, data=appendix.data) +
  # Add model fits
  geom_line(data=filter(predictions, Performance>=-0.2, CompetitorDensity<=8.5)) + 
  # Add resident equilibrium calculation
  geom_vline(aes(xintercept=ResidentEquilibrium, color=SoilTreatment), 
             linetype='dashed', 
             # Show resident equilibrium for the competitor species in each panel
             data=mutate(
               filter(resident.equilibria, !is.na(ResidentEquilibrium)), 
               CompetitorSpecies=FocalSpecies)) +
  # Add invasion performance calculation
  geom_point(aes(x=ResidentEquilibrium, y=InvasionPerformance), 
             shape=17, size=2,
             # Show invasion performance in the heterospecific panel
             data=mutate(invasion.analysis, 
                         CompetitionTreatment='Heterospecific', 
                         CompetitorSpecies=ResidentSpecies)) +
  # Add a horizontal line segment showing invasion performance
  geom_segment(aes(x=0, y=InvasionPerformance,
                   xend=ResidentEquilibrium, yend=InvasionPerformance),
               linetype='dashed',
               data=mutate(invasion.analysis, 
                           CompetitionTreatment='Heterospecific', 
                           CompetitorSpecies=ResidentSpecies)) +
  # Use facets to show competitor identity (column) and hetero-/conspecific 
  # (row) combinations separately
  facet_grid(CompetitionTreatment ~ paste('Competitor: ', CompetitorSpecies),
             scales="free") +
  theme_minimal()

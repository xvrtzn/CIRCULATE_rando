
#  From CIRCULATE protocol
# Before randomization several points have to be checked:
#   ???	Signature of ICF,
# ???	Eligibility criteria,
# ???	ctDNA assessment performed and filled in CRF.
# -	Patients who are positive for postoperative ctDNA (ctDNA+) will be randomized (2:1) 
# to adjuvant chemotherapy (n=132) or to follow up (n=66).
# The stratification factors are: Center / T3 vs T4a / emergency resection.

# -	Patients who are negative for postoperative ctDNA (ctDNA-) will be randomized (1:4) 
# to follow-up within CIRCULATE or to routine follow up outside the trial  
# The stratification factors are: Center / T3 vs T4a / emergency resection.

# Stratification factors 
# nb of centers (100 levels) : approx 100 => Center001 -> Center100
# tumor stage (2 levels) : T3, T4a
# emergency resection (2 levels) : Yes, No

# Method as described in https://www.lexjansen.com/nesug/nesug04/ap/ap07.pdf

library(tidyverse)
library(pbapply) # Apply with progress bar

# Global parameters 
theme_set(theme_bw())

# create seed for program to be rerun with same results
set.seed(42)

# Reinitialize random seed
# set.seed(NULL)

# Load utility functions
source("randoUtils.r")


## Liste des param�tres pour l'�tude CIRCULATE ------

numberOfCenters <- 100
allTreatments <- c("A", "B")
nbTreatments <- length(allTreatments)
allCenters <- sprintf("Center%03d", as.numeric(1:numberOfCenters))
allStages <- c("T3", "T4a")
allEmergencyResections <- c("Yes", "No")
treatmentAllocationOdds <- c(2, 1)
treatmentAllocationProbs_vec <- treatmentAllocationOdds / sum(treatmentAllocationOdds)

# Create complete DF of all possible combinations of stratification factors + treatments 
allCases <- 
  crossing(
    Center = allCenters, 
    Stage = allStages, 
    EmergencyResection = allEmergencyResections,
    Treatment = allTreatments
    )


# define study-level parameters for randomization
# Treatment allocation probabilities
treatmentAllocationProbs <- 
  tribble(
    ~Treatment, ~Weight, ~Odds,
    allTreatments, treatmentAllocationProbs_vec, treatmentAllocationOdds
  ) %>%  
  unnest()

## Build a virtual trial by incrementally randomizing randomly sampled patients ----

performOneTrial <- function(nPatients, dist_method, choice_method) {
  
  # Create first patient and allocate treatment randomly (using the predefined ratio ebtween treatment groups)
  firstPatient_fullTrial <- samplePatients(1, allCases) 
  firstPatient_fullTrial$Treatment <- 
    sample( allTreatments, 1, replace = TRUE, prob = treatmentAllocationProbs_vec )
  
  patientList_fullTrial <- firstPatient_fullTrial
  
  for (idx in 1:(nPatients-1)) {
    
    # Create new patient
    newPatient_fullTrial <- samplePatients(1, allCases) %>% 
      mutate( patientID = max(patientList_fullTrial$patientID) + 1 )
    
    # Compute imabalance scores for new patient
    ImbalanceScores <- patientList_fullTrial %>% 
      ComputePatientListBreakdown(allCases) %>% 
      ComputeExpectedValues(treatmentAllocationProbs, newPatient_fullTrial) %>% 
      ComputeImbalanceScores(allTreatments )
    
    # Allocate treatment for new patient
    newPatient_fullTrial$Treatment <- 
      selectTreatment(ImbalanceScores, treatmentAllocationProbs,  dist_method, choice_method)
    
    patientList_fullTrial <- rbind(patientList_fullTrial, newPatient_fullTrial) 
    
  }
  
  return( patientList_fullTrial )
  
}

oneTrial <- function( Params ) {
  performOneTrial(200, Params[1], Params[2]) %>% 
    mutate( nb_A_before = cumsum(Treatment=="A"), 
            nb_B_before = cumsum(Treatment=="B"), 
            prop_A_B_before = nb_A_before / nb_B_before, 
            distParam = Params[1], 
            choiceParam = Params[2] ) 
    
}

test_oneTrial <- oneTrial( c("MAX", "PROP") )

# Now visualize allocation proportion as a function of n
# Visualize 10 replications of one given set of parameters
pbreplicate(3, oneTrial( c("MAX", "PROP") ), simplify = FALSE) %>% 
  bind_rows( .id = "TrialIDX" ) %>% 
  ggplot( aes( x = patientID, y = prop_A_B_before, group = TrialIDX, color = Stage ) ) +
  geom_line() +
  facet_grid(distParam+Stage~choiceParam+EmergencyResection) +
  ylim(0, 5)

# Create complete simulation with all combinations of parameters
# Create list of all possible randomization parameters
allParams <- 
  crossing( distParam = c("RANGE", "VAR", "MAX"), 
            choiceParam = c("BEST", "PROB", "PROP")) %>% 
  t() %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  as.list()

df_all <- pbreplicate(5, map_df(allParams, oneTrial), simplify = FALSE) %>% 
  bind_rows( .id = "TrialIDX" ) 

# View all data - facet by Minimization methods + stratification factors
df_all %>%
  ggplot( aes( x = patientID, 
               y = prop_A_B_before, 
               group = TrialIDX ) ) +
  geom_line( alpha = I(0.3), size=1) +
  facet_grid( distParam + Stage ~ choiceParam + EmergencyResection ) +
  # facet_grid( distParam  ~ choiceParam  ) +
  ylim(0, 5)

# Try summarising to enhance visualization - not great yet
df_all %>%
  # mutate( TrialIDX = paste0( "0", as.numeric(TrialIDX) ) )%>% 
  # bind_rows(df_all1, .id = "Weight" ) %>% 
  ggplot( aes( x = patientID, y = prop_A_B_before ) ) +
  stat_summary(fun.data="mean_sdl", geom="ribbon", alpha = 0.2 ) +
  stat_summary(fun.y=median, geom="line", color = "blue" ) +
  facet_grid(distParam~choiceParam) +
  ylim(0, 5)

# Code to load data on laptop
# saveRDS(df_all, "rando_simResults.Rds")
df_all = readRDS("rando_simResults.Rds")


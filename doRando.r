
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

# create seed for program to be rerun with same results
# TODO :  must be returned as well
set.seed(42)

set.seed(NULL)


###########################

numberOfCenters <- 100
allTreatments <- c("A", "B", "C")
allCenters <- sprintf("Center%03d", as.numeric(1:numberOfCenters))
allStages <- c("T3", "T4a", "T5")
allEmergencyResections <- c("Yes", "No", "Maybe", "Whatever")
samplingScheme <- c(2, 1, 1)
allTreatmentsProbs <- samplingScheme / sum(samplingScheme)

# Create complete DF of all possible combinations of stratification factors + treatments 
allCases <- 
  crossing(
    Center = allCenters, 
    Stage = allStages, 
    EmergencyResection = allEmergencyResections,
    Treatment = allTreatments
    )

# Generate a random list of n patients already allocated with a treatment 
# Treatment are randomly allocated with unbalanced probs (not proper randomization but this is a start)

samplePatients <- function( n = 1 ) {
  tibble( 
    Center = sample( allCenters, n, replace = TRUE ), 
    Stage = sample( allStages, n, replace = TRUE  ), 
    EmergencyResection = sample( allEmergencyResections, n, replace = TRUE  ),
    Treatment = sample( allTreatments, n, replace = TRUE, prob = allTreatmentsProbs )
  ) 
}

patientList <- 
  samplePatients(100) %>% 
  mutate(patientID = row_number()) 

################################################################
# Minimization algorithm - 2 treatments
# TODO : generalize algorithm to more than 2 trts

# define study-level parameters for randomization
# Treatment allocation probabilities
treatmentAllocationProbs <- 
  tribble(
    ~Treatment, ~Weight, 
    allTreatments, allTreatmentsProbs
  ) %>%  
  unnest()

# Get weights for stratification factors into a vector for future use
# TODO

# Create one new patient to be randomized
newPatient <- samplePatients(1) %>% 
  select(-Treatment) %>% 
  mutate(patientID = nrow(patientList) + 1 )

# Create table summarizing number of cases within each stratum of the 
# stratification factors in the current patient list
# DONE - the count should return zero for centers that are not used
# DONE - this needs to be generalized to handle more stratif factors (using gather/spread?)

# Create Table 3 - Breakdown of patient list before randomization 
# of new patient (complete with cases not seen before)
ComputePatientListBreakdown <- function(allCases, patientList) {
  allCases %>% 
    left_join(patientList, by = names(allCases)) %>% 
    gather( key = "key", value = "value", -patientID, -Treatment ) %>% 
    group_by( Treatment, key, value ) %>% 
    summarise( n = sum(!is.na(patientID) ) )
}

patientListBreakdown <- ComputePatientListBreakdown(allCases, patientList)

# Create Tables 7 and 8 - expected values only
ComputeExpectedValues <- function(treatmentAllocationProbs, newPatient, patientListBreakdown) {
  
  newPatient %>% 
    select( -patientID ) %>% 
    gather( key = "key", value = "value" ) %>% 
    left_join( patientListBreakdown, by = c( "key", "value" ) ) %>% 
    group_by( key ) %>% 
    mutate( tot_byStratFactor = sum( n ) ) %>% 
    left_join( treatmentAllocationProbs, by = "Treatment") %>% 
    mutate( expected = ( tot_byStratFactor + 1 ) * Weight )
  
}

ExpectedValues <- ComputeExpectedValues(treatmentAllocationProbs, newPatient, patientListBreakdown)

# Create Table 10 - differences observed-expected, distances (3 methods) and 
# imbalance scores (3 methods)
ComputeImbalanceScores <- function(ExpectedValues, allTreatments) {
  
  ExpectedValues  %>% 
    crossing( newTreatment = allTreatments ) %>% 
    mutate( n_new = ifelse(Treatment == newTreatment, n + 1, n) ) %>% 
    mutate( diff = n_new - expected ) %>%
    group_by( newTreatment, key ) %>% 
    summarise( dist_range = max(diff)-min(diff), 
               dist_var = var(diff) * ( n() - 1 ) / n(), 
               dist_max = max(diff) )  %>% 
    summarise( imbalanceScore_range = sum(dist_range), 
               imbalanceScore_var = sum(dist_var), 
               imbalanceScore_max = sum(dist_max)) 
  
}

ImbalanceScores <- ComputeImbalanceScores(ExpectedValues, allTreatments )

# Select treatment 
# Method BEST : assign treatment with lowest imbalance score
# Method PROB : assign treatment with lowest imbalance score, using fixed probabilities
# Method PROP : assign treatment with lowest imbalance score, using probabilities 
#               based on imbalance score values
selectTreatment <- function(ImbalanceScores, dist_method, choice_method) {
  
  dist_colName <- 
    case_when(
      dist_method == "RANGE"  ~ "imbalanceScore_range",
      dist_method == "VAR"  ~ "imbalanceScore_var",
      dist_method == "MAX"  ~ "imbalanceScore_max",
      TRUE ~ "UNKNOWN_COLNAME"
    )

  dist_colName <- enquo(dist_colName)
    
  reducedImbalanceScores <- ImbalanceScores %>% 
    select( newTreatment, imbalanceScore = !!dist_colName )
  
  assignedTreatment <- "UNKNOWN"
  
  if (choice_method == "BEST") {
    
    # print("BEST method chosen")

    assignedTreatment <- reducedImbalanceScores %>% 
      filter( imbalanceScore == min(imbalanceScore) ) %>% 
      dplyr::pull(newTreatment)

  } else if (choice_method == "PROB") {
    
    # print("PROB method chosen")
    
    sortedTreatments <- ImbalanceScores %>% 
      select(newTreatment, imbalanceScore = imbalanceScore_range) %>% 
      arrange(imbalanceScore)  %>% 
      dplyr::pull(newTreatment)
    
    assignedTreatment <- sample( sortedTreatments, 1, prob =  c(0.75, 0.25) )

  } else if (choice_method == "PROP") {

    # print("PROP method chosen")
    
    sortedImbalanceScores <- ImbalanceScores %>% 
      select(newTreatment, imbalanceScore = imbalanceScore_range) %>% 
      arrange(imbalanceScore)  
      
    sortedTreatments <- dplyr::pull(sortedImbalanceScores, newTreatment)
    sortedScores <- dplyr::pull(sortedImbalanceScores, imbalanceScore)
    
    assignedTreatment <- sample( sortedTreatments, 
                                 1, 
                                 prob =  1 - sortedScores / sum(sortedScores) 
                                 )

  } else {
    
    print("Unknown method to choose treatment")
    
  }
  
  return(assignedTreatment)

}

# Test once
allocatedNewPatient <- newPatient 
allocatedNewPatient$Treatment <- selectTreatment(ImbalanceScores, "RANGE", "BEST")
allocatedNewPatient$status <- "NEW"

# Test 1000 times and plot histogram of treatment allocation
replicate( 1000, selectTreatment(ImbalanceScores, "VAR", "PROP") ) %>% 
  as.tibble() %>% 
  ggplot( aes( x = value ) ) +
  geom_bar()

patientList %>% 
  mutate( status = "OLD" ) %>% 
  bind_rows(allocatedNewPatient) %>% 
  ggplot( aes( x = Treatment, fill = Stage ) ) +
  geom_bar( position = position_dodge() ) + 
  facet_grid(status ~ EmergencyResection )

# It looks like the OLD patient list is not properly stratified. This is normal, as it has been 
# randomly created without stratification!!! To create a properly randomized patient list, randomization
# would need to be used incrementally to built the list



# Generate a random list of n patients already allocated with a treatment 
# Treatment are randomly allocated with unbalanced probs (not proper randomization but this is a start)

samplePatients <- function( n = 1, allCases ) {
  allCases %>% 
    dplyr::select(-Treatment) %>% 
    sample_n( n, replace = TRUE ) %>% 
    mutate( patientID = row_number() )
}


# Create table summarizing number of cases within each stratum of the 
# stratification factors in the current patient list
# DONE - the count should return zero for centers that are not used
# DONE - this needs to be generalized to handle more stratif factors (using gather/spread?)

# Create Table 3 - Breakdown of patient list before randomization 
# of new patient (complete with cases not seen before)
ComputePatientListBreakdown <- function(patientList, allCases) {
  allCases %>% 
    left_join(patientList, by = names(allCases)) %>% 
    gather( key = "key", value = "value", -patientID, -Treatment ) %>% 
    group_by( Treatment, key, value ) %>% 
    summarise( n = sum(!is.na(patientID) ) )
}


# Create Tables 7 and 8 - expected values only
ComputeExpectedValues <- function(patientListBreakdown, treatmentAllocationProbs, newPatient) {
  
  newPatient %>% 
    select( -patientID ) %>% 
    gather( key = "key", value = "value" ) %>% 
    left_join( patientListBreakdown, by = c( "key", "value" ) ) %>% 
    group_by( key ) %>% 
    mutate( tot_byStratFactor = sum( n ) ) %>% 
    left_join( treatmentAllocationProbs, by = "Treatment") %>% 
    mutate( expected = ( tot_byStratFactor + 1 ) * Weight )
  
}



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


# Select treatment 
# Method BEST : assign treatment with lowest imbalance score
# Method PROB : assign treatment with lowest imbalance score, using fixed probabilities
# Method PROP : assign treatment with lowest imbalance score, using probabilities 
#               based on imbalance score values
selectTreatment <- function(ImbalanceScores, treatmentAllocationProbs, dist_method, choice_method) {
  
  dist_colName <- 
    case_when(
      dist_method == "RANGE"  ~ "imbalanceScore_range",
      dist_method == "VAR"  ~ "imbalanceScore_var",
      dist_method == "MAX"  ~ "imbalanceScore_max",
      TRUE ~ "UNKNOWN_COLNAME"
    )
  
  dist_colName <- enquo(dist_colName)
  
  reducedImbalanceScores <- ImbalanceScores %>% 
    select( newTreatment, imbalanceScore = !!dist_colName ) %>% 
    left_join( treatmentAllocationProbs, by = c("newTreatment" = "Treatment"))
  
  assignedTreatment <- "UNKNOWN"
  
  if (choice_method == "BEST") {
    
    # print("BEST method chosen")
    
    assignedTreatment <- reducedImbalanceScores %>% 
      filter( imbalanceScore == min(imbalanceScore) ) %>% 
      dplyr::pull(newTreatment)
    
  } else if (choice_method == "PROB") {
    
    # print("PROB method chosen")

    sortedTreatments <- reducedImbalanceScores %>% 
      arrange(imbalanceScore)  %>% 
      dplyr::pull(newTreatment)
    
    # Choose the treatment with lowest imbalanced score, with fixed probablity 75%
    assignedTreatment <- sample( sortedTreatments, 1, prob = c(0.75,0.25)  )
    
  } else if (choice_method == "PROP") {
    
    # print("PROP method chosen")
    
    sortedImbalanceScores <- reducedImbalanceScores %>% 
      mutate( Weight_inverseImbalanceScores =  1 - imbalanceScore / sum(imbalanceScore), 
              Weight_Odds = Weight_inverseImbalanceScores * Odds,  
              TrtAllocProb = Weight_Odds / sum(Weight_Odds) ) %>% # Normalize to transform into probability
      arrange(imbalanceScore)  
      
    
    sortedTreatments <- dplyr::pull(sortedImbalanceScores, newTreatment)
    TrtAllocProb <- dplyr::pull(sortedImbalanceScores, TrtAllocProb)
    
    # Choose the treatment with lowest imbalanced score, with weighted probablity
    assignedTreatment <- sample( sortedTreatments, 1, prob =  TrtAllocProb )

  } else {
    
    print("Unknown method to choose treatment")
    
  }
  
  return(assignedTreatment)
  
}



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

###########################

numberOfCenters <- 100
allTreatments <- c('A', 'B')
allCenters <- sprintf("Center%03d", as.numeric(1:numberOfCenters))
allStages <- c('T3', 'T4a')
allEmergencyResections <- c('Yes', 'No')
treatmentAllcocationProbs <- c(0.6, 0.4)


# Generate a random list of n patients already allocated with a treatment 
# Treatment are randomly allocated with unbalanced probs (not proper randomization but this is a start)

samplePatients <- function( n = 1 ) {
  tibble( 
    Center = sample( allCenters, n, replace = TRUE ), 
    Stage = sample( allStages, n, replace = TRUE  ), 
    EmergencyResection = sample( allEmergencyResections, n, replace = TRUE  ),
    Treatment = sample( allTreatments, n, replace = TRUE, prob = treatmentAllcocationProbs )
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
treatmentAllocationProbs <- tribble(
  ~Treatment, ~Weight, 
  allTreatments, c(0.6, 0.4)
) %>%  unnest()

# Get weights for stratification factors into a vector for future use
# TODO


# Create one new patient to be randomized
newPatient <- samplePatients(1)

# Create table summarizing number of cases within each stratum of the 
# stratification factors in the current patient list
# DONE - the count should return zero for centers that are not used
# TODO - this needs to be generalized to handle more stratif factors (using gather/spread?)

# Create *complete* list of combinations of stratification variables, together with count of occurrences in actual patientList
# The first mutate is used to complete the data frame with all possible combinations (even the ones not seen yet). Makes it tough to generalize to any number of stratification factor though

patientList_table <- patientList %>% 
  mutate( Center = factor(Center, levels = allCenters),
          Stage = factor(Stage, levels = allStages), 
          EmergencyResection = factor(EmergencyResection, levels = allEmergencyResections),
          Treatment = factor(Treatment, levels = allTreatments) ) %>% 
  complete(Treatment, Center, Stage, EmergencyResection) %>% 
  group_by(Treatment, Center, Stage, EmergencyResection) %>%
  summarise(non_na_count = sum(!is.na(patientID))) %>% 
  mutate_if( is.factor, as.character)

# Summary of counts for each stratification variable
patientList_table_summ <- patientList_table %>% 
  gather(key = "key", value = "value", -non_na_count, -Treatment) %>% 
  group_by(Treatment, key, value) %>% 
  summarise( n = sum(non_na_count) ) 

# Table 4 in Levy and Blood paper 
# This works for any number of stratification factors
TAB4 <- newPatient %>% 
  gather(key = "key", value = "value") %>% 
  left_join(patientList_table_summ, by=c('key', 'value'))

# Show as a standard table
TAB4 %>% ungroup() %>% select(-key) %>% spread( key = value, value = n)

# TODO - build the table of counts for the different allocation strategies
# If possible, make it work with more than two treatments

# TODO - table of expected values for each cell of TAB4 (tables 7 and 8 in paper)
TAB78 <- TAB4 %>% 
  group_by(key) %>% 
  mutate( tot = sum(n)) %>% 
  left_join(treatmentAllocationProbs, by = "Treatment") %>% 
  mutate( expected = tot * Weight ) 
  

# TODO - Compute Imbalance score (use RANGE method)
# TODO (later) - implement VARIANCE and MAX methods
TAB78_test <- tribble(
  ~Treatment, ~key, ~value, ~n_plusone, ~expected,
  c('A', 'B', 'A', 'B', 'A', 'B', 'A', 'B'), 
  c('F1', 'F1', 'F1', 'F1', 'F2', 'F2', 'F2', 'F2'), 
  c('five', 'five', 'five', 'five', 'three', 'three', 'three', 'three'), 
  c(7,5,6,6, 5,4,4,5), 
  c(6,6,6,6,4.5,4.5,4.5,4.5)
) %>% unnest()


Imbalance_table <- TAB78_test %>% 
  mutate( diff = n_plusone-expected )
  


# TODO - implement treatment selection method



################################ SANDBOX
# Create data frame of stratification information
# TODO : how to store weights in this DF?
stratification_df <- 
  tribble(
    ~Factor, ~Value, ~Weight,
    "Center", allCenters, 1/3,
    "Stage", allStages, 1/3,
    "EmergencyResection", allEmergencyResections , 1/3
  ) %>% unnest() 



# Create all possible unique patient characterstics - NOT USED
allPossiblePatients <- expand.grid( Center = allCenters, 
                                    Stage = allStages, 
                                    EmergencyResection = allEmergencyResections,
                                    Treatment = allTreatments
)

# TODO : count number of times a row of allPossiblePatients is present in patientList - ????
test <- intersect(select(patientList, -patientID), allPossiblePatients)

############################################################
# First try a simple random allocation
simpleRando <- function(df, trt) {
  
  df$Treatment <- sample( trt, nrow(df), replace = TRUE, prob = c(0.6, 0.4) )
  
  return(df)
  
}

# Is the sampling balanced? Looks like it
simpleRando(patientList, allTreatments) %>% 
  group_by(Treatment) %>% 
  count()

# Example illustration of an unbalanced randomization
# Sample A,B 100 times unsing unbalanced probs, repeat 1000 times and check histogram of number of As
replicate(1000, 
          sample(treatments, 100, replace=T, prob = c(0.6, 0.4))  %>% 
            as.vector() %>% 
            table() %>% 
            pluck(1) ) %>% 
  as_tibble() %>% 
  ggplot( aes(x=value) ) +
  geom_histogram( binwidth = 1 )

################ TEST DATA (from Levy and Blood article)

# Tabel of already included patients
patientList_test <- tribble(
  ~Treatment, ~patientID, ~stratFactor1, ~stratFactor2,
  "B", "11001", "6", "4",
  "A", "11002", "5", "3",
  "A", "11003", "5", "3",
  "B", "11004", "6", "3",
  "A", "11005", "5", "4",
  "B", "11006", "5", "4",
  "A", "11007", "6", "4",
  "B", "11008", "5", "3",
  "A", "11009", "6", "4",
  "A", "11010", "5", "4",
  "B", "11011", "6", "3",
  "A", "11012", "5", "3",
  "A", "11013", "5", "3",
  "B", "11014", "5", "3",
  "B", "11015", "5", "4",
  "B", "11016", "5", "4"
)

# Function to get all unique values of a given column in a data frame
allUniqueValues <- function(df, col) {
  col <- enquo(col)

  df %>%
    dplyr::pull(!!col) %>% 
    unique() 
}

allUniqueValues(patientList_test, Treatment)
allUniqueValues(patientList_test, stratFactor2)

test <- patientList_test %>% 
  complete( Treatment, stratFactor1, stratFactor2 ) %>% 
  group_by( Treatment, stratFactor1, stratFactor2 ) %>%
  summarise(non_na_count = sum(!is.na(patientID))) 

# New patient for which to determine best treatment allocation
newPatient_test <- data_frame( patientID = "11017", stratFactor1 = "5",  stratFactor2 = "3")



# Table of treatment allocation probabilities
treatmentAllocationProbs_test <- tribble(
  ~Treatment, ~Weight, 
  c("A", "B"), c(0.5, 0.5)
) %>%  unnest()


# Create Table 3
TAB3_test <- patientList_test %>% 
  gather(key, value, -patientID, -Treatment) %>% 
  group_by(Treatment, key, value) %>% 
  count()  
#  unite( key, value, col = "stratFactor", sep="_") 

# Create Table 4
TAB4_test <- newPatient_test %>% 
  select(-patientID) %>% 
  gather(key = "key", value = "value") %>% 
  left_join(TAB3_test, by=c('key', 'value'))

# Create Table 10
TAB_expected_test  %>% 
  crossing(newTreatment=c("A", "B")) %>% 
  mutate( n_new = ifelse(Treatment == newTreatment, n+1, n) ) %>% 
  mutate( diff = n_new - expected ) %>%
  arrange(newTreatment, key) %>% 
  group_by(newTreatment, key) %>% 
  summarise( dist_range = abs(max(diff)-min(diff)), dist_var = var(diff), dist_max = max(diff) )  %>% 
  summarise( imbalanceScore_range = sum(dist_range), imbalanceScore_var = sum(dist_var), imbalanceScore_max = sum(dist_max),)


crossing(stratFactor1 = c("5", "6"),stratFactor2 = c("3", "4"), Treatment=c("A", "B"))


# Create Tables 7 and 8
TAB_expected_test <- TAB4_test %>% 
  group_by(key) %>% 
  mutate( tot = sum(n) ) %>% 
  left_join(treatmentAllocationProbs_test, by = "Treatment") %>% 
  mutate( expected = (tot+1) * Weight )



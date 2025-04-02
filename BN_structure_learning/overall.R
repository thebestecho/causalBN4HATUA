library(readr)
library(dplyr)
library(bnlearn)
library(boot)


load("20240328_data_subset.RData")
load("20240528_bl_no_hos_no_site.RData")

bl <- as.data.frame(bl)

## Define the bootstrapping function to generate bootstrapped datasets
## selecting rows from the original dataset based on the randomly generated indices with replacement
bootstrap_data <- function(data, indices) {
  
  boot_sample <- data[indices, ]
  return(boot_sample)

}


## Define a function for structure learning using structural.em
learn_structure <- function(boot_data) {

  structural.em(boot_data, maximize = "hc", 
                maximize.args = list(score = "bde", 
                                     iss = 18,
                                     maxp = 3,
                                     max.iter = 1000,
                                     blacklist = bl),
                max.iter = 1000,
                return.all = TRUE)

}


## Function to perform bootstrapping and structure learning
bootstrap_and_learn <- function(seed,index) {
  
  seed_value <- seed + index

  n <- nrow(data_subset)
  
  ## Generate random indices for bootstrapping
  set.seed(seed_value)
  boot_indices <- sample(n, replace = TRUE)
  
  ## Use the indices to create a bootstrapped dataset
  boot_data <- bootstrap_data(data_subset, boot_indices)
  
  ## Learn the structure using structural.em
  learned <- learn_structure(boot_data)
  
  return(learned)
  
}


## ---------------------------------- Get 1000 average networks ----------------------------------

## Set the number of bootstrap samples
args <- commandArgs(trailingOnly = TRUE)
array_ID <- as.numeric(args[1])  # Convert the argument to numeric

print(array_ID)

## Perform bootstrap resampling and collect learned structures
results <- bootstrap_and_learn(seed = 121, index = array_ID)

## Obtain the structure
learned_structure <- results$dag

## Obtain the imputed data set
imputed_data <- results$imputed

# Directories for saving the results
structure_dir <- "./overall_structures"
imputed_dir <- "./overall_imputed"

# Create the directories if they don't exist
if (!dir.exists(structure_dir)) dir.create(structure_dir, recursive = TRUE)
if (!dir.exists(imputed_dir)) dir.create(imputed_dir, recursive = TRUE)

# Save the structure and imputed data
save(learned_structure, file = paste0(structure_dir, "/structure_", array_ID, ".RData"))
write.csv(imputed_data, file = paste0(imputed_dir, "/imputed_", array_ID, ".csv"), row.names = FALSE)


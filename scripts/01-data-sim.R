# Load the required libraries
library(MASS)
library(data.table)
library(furrr)

# Load the helper functions
source("scripts/00-functions.R")

# Set the Parameters for the Simulation
alpha = -4.50               # Intercept
beta = 0.00                 # Treatment Effect
n = 50000                   # Sample Size
sigma = c(1.0, 1.0)         # Standard Deviation of the Potential Outcomes
rho = 0.00                  # Correlation between the Potential Outcomes

# Set the Parallel Processing Plan
plan(multisession, workers = 4)

# Simulate 10,000 datasets
aa_datasets <- future_map(
    .x = 1:10000,
    ~ sim_data(alpha, beta, n, sigma, rho, 0.5),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

# Close the Parallel Processing Plan
plan(sequential)

# Append the Results to a Single Data Frame
aa_datasets <- rbindlist(aa_datasets, idcol = "dataset_id")
aa_datasets[, treat := factor(treat_id, labels = c("Control", "Treatment"))]

# Write the data to a CSV file
fwrite(aa_datasets, "data/aa_datasets.csv")

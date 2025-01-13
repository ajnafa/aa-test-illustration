# Load the libraries
library(data.table)
library(marginaleffects)

# Load the helper functions
source("scripts/00-functions.R")

# Read in the simulated datasets
aa_datasets <- fread("data/aa_datasets.csv")

# Sort the data by dataset and treatment id
aa_datasets <- aa_datasets[order(dataset_id, treat_id)]

# Parallel computation via furrr
plan(multisession(workers = 3))

# Fit models and add the draws to the simulation data frame
sim_estimates <- future_map(
    .x = unique(aa_datasets$dataset_id),
    .f = ~ freq_model(data = aa_datasets[dataset_id == .x]),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
)

# Combine the results
out <- rbindlist(sim_estimates, idcol = "dataset_id")

# Write the data to a CSV file
fwrite(out, "data/aa_frequentist_estimates.csv")



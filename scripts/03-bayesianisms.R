# Load the libraries
library(data.table)
library(cmdstanr)
library(furrr)

# Load the helper functions
source("scripts/00-functions.R")

# Read in the simulated datasets
aa_datasets <- fread("data/aa_datasets.csv")

# Sort the data by dataset and treatment id
aa_datasets <- aa_datasets[order(dataset_id, treat_id)]

# Nest the data by dataset
aa_data <- aa_datasets[, 
    list(
        dataset=list(data.table(treat_id, trials, Y_obs, Y_mis, treat))
    ), by=dataset_id
]

# Install cmdstan if not alreafy installed, requires 2.34 or later
#cmdstanr::install_cmdstan(cores = 4)

# Compile the Stan model
model <- cmdstan_model("models/binomial-logit.stan")

# Make the Stan data for the models
stan_data_list <- purrr::map(
    .x = seq_along(aa_data$dataset_id),
    .f = ~ make_stan_data(
        data = aa_data$dataset[[.x]],
        mu_alpha = -3.5,
        sigma_alpha = 1,
        mu_beta = 0,
        sigma_beta = 0.5
    )
)

# Parallel computation via furrr
plan(tweak(multisession, workers = 3))

# Fit models and add the draws to the simulation data frame
sim_draws <- future_map(
    .x = stan_data_list,
    .f = ~ bayes_models(
      stan_data = .x,
      stan_model = model,
      sampling = 1e3,
      warmup = 1e3,
      chains = 4
    ),
    .options = furrr_options(
      scheduling = 1,
      seed = TRUE,
      prefix = "prefix"
    ),
    .progress = TRUE
)

# Close the parallel processing plan
plan(sequential)

# Append the estimates to a single data frame
draws <- rbindlist(sim_draws, idcol = "dataset_id")

# Write the estimates to a CSV file
fwrite(draws, "data/aa_bayesian_estimates.csv")

# Function for simulating the data
sim_data <- function(alpha, beta, n, sigma, rho, prob) {

    # Generate the Treatment Indicator
    W = rbinom(n, 1, prob)

    # Linear Predictor for the Potential Outcomes
    mu = c(alpha, alpha + beta)

    # Covariance Matrix for the Potential Outcomes
    Sigma = diag(length(mu))

    # Potential Outcomes on the Latent Scale
    Mu = MASS::mvrnorm(n, mu, Sigma)

    # Observed and Missing Outcomes on the Response Scale
    Y_obs = rbinom(n, 1, plogis(Mu * W + (1 - W) * Mu))
    Y_mis = rbinom(n, 1, plogis(Mu * (1 - W) + W * Mu))

    # Create the Data Frame
    data = data.table(
        treat_id = W + 1,
        Y_obs = Y_obs,
        Y_mis = Y_mis
    )

    data = data[,
        .(
            trials = .N,
            Y_obs = sum(Y_obs),
            Y_mis = sum(Y_mis)
        ), by = treat_id
    ]

    return(data)
}

# Function for Building the Data to Pass to the Stan Model
make_stan_data <- function(data, mu_alpha, sigma_alpha, mu_beta, sigma_beta, ...) {
    
    # Prepare the data for the Binomial Model
    X <- model.matrix(~ treat, data = data)
    stan_data <- list(
        N = nrow(data),
        D = ncol(X),
        Y = data$Y_obs,
        K = data$trials,
        P = X,
        mu_alpha = mu_alpha,
        sigma_alpha = sigma_alpha,
        mu_beta = mu_beta,
        sigma_beta = sigma_beta
    )

    return(stan_data)
}


# Function for fitting the Stan Models in Parallel
bayes_models <- function(stan_model, 
                         stan_data,
                         sampling, 
                         warmup, 
                         chains,
                        ...) {
  
  # Set the initial number of draws to 0
  min_draws <- 0
  
  # We need a while loop here because sometimes one of the chains will fail
  # unexpectedly, so we have to re-run the model. It's not clear why this happens 
  # and it only seems to be an issue when running models that finish quickly in 
  # parallel.
  while (min_draws < (sampling * chains)) {
    
    # Fit the Stan Model
    fit <- stan_model$sample(
      data = stan_data,
      sig_figs = 5,
      parallel_chains = chains,
      iter_warmup = warmup,
      iter_sampling = sampling,
      refresh = 0,
      show_messages = FALSE,
      ...
    )
    
    # Update the check
    min_draws <- posterior::ndraws(fit$draws())
  }
  
  # Extract the posterior draws
  draws <- fit$draws(c("post_prob"), format = "draws_df")

  # Store the probabilities
  estimates <- data.table(
      control_prob = mean(draws$`post_prob[1]`),
      treat_prob = mean(draws$`post_prob[2]`)
  )
  
  # Return the data frame of draws
  return(estimates)
}


freq_model <- function(data, ...) {
    
    # Fit the model
    fit <- glm(
        cbind(Y_obs, trials - Y_obs) ~ treat, 
        family = binomial(),
        data = data
    )

    # Compute the average marginal effects
    out <- avg_comparisons(
        fit,
        at = "treat",
        type = "response"
    ) |> tidy()

    return(out)
}
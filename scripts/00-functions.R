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

# Custom theme for data visualizations
plot_theme <- function(title_size = NULL, 
                       xaxis_size = NULL, 
                       yaxis_size = NULL, 
                       strip_size = NULL, 
                       strip_face = NULL, 
                       caption.hjust = 1, 
                       caption.vjust = 0, 
                       x_axis_face = NULL, 
                       y_axis_face = NULL, 
                       transparent = FALSE, 
                       axis_text_size = NULL, 
                       legend_text_size = NULL,
                       subtitle_size = NULL,
                       caption_size = NULL,
                       base_size = 12,
                       ...) {
  .theme <- theme_minimal(base_size = base_size) + theme(
    # Specify the default settings for the plot title
    plot.title = element_text(
      size = title_size,
      face = "bold",
      family = "serif"
    ),
    # Specify the default settings for caption text
    plot.caption = element_text(
      size = caption_size,
      family = "serif",
      hjust = caption.hjust,
      vjust = caption.vjust
    ),
    # Specify the default settings for subtitle text
    plot.subtitle = element_text(
      size = subtitle_size,
      family = "serif"
    ),
    # Specify the default settings specific to the x axis title
    axis.title.y = element_text(
      size = yaxis_size, 
      face = y_axis_face, 
      family = "serif",
      margin = margin(r = 10, l = -10)
    ),
    # Specify the default settings specific to the y axis title
    axis.title.x = element_text(
      size = xaxis_size, 
      face = x_axis_face, 
      family = "serif",
      margin = margin(t = 10, b = -10)
    ),
    # Specify the default settings for x axis text
    axis.text.x = element_text(
      size = axis_text_size,
      family = "serif",
      face = x_axis_face
    ),
    # Specify the default settings for y axis text
    axis.text.y = element_text(
      size = axis_text_size,
      family = "serif",
      face = y_axis_face
    ),
    # Specify the default settings for legend titles
    legend.title = element_text(
      size = legend_text_size,
      face = "bold",
      family = "serif"
    ),
    # Specify the default settings for legend text
    legend.text = element_text(
      size = legend_text_size,
      family = "serif"
    ),
    # Set the strip background fill to blank
    strip.background = element_blank(),
    # Adjust the strip text size settings
    strip.text = element_text(
      family = "serif", 
      size = strip_size,
      face = strip_face
    ),
    # Additional Settings Passed to theme()
    ...
  )
  # Plot Transparency
  if (transparent == TRUE) {
    .theme <- .theme + theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
  }
  return(.theme)
}
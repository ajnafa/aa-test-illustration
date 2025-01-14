# Load the required libraries
library(data.table)
library(ggplot2)

# Load the helper functions
source("scripts/00-functions.R")

# Set the ggplot2 theme
theme_set(plot_theme(
    xaxis_size = 22,
    yaxis_size = 22,
    title_size = 25,
    caption_size = 18,
    axis_text_size = 20,
    strip_size = 20,
    plot.caption.position = "plot",
    plot.title.position = "plot",
    legend_text_size = 20,
    legend.position = "top",
    strip_face = "bold",
    caption.hjust = 0,
    caption.vjust = -1,
    subtitle_size = 20,
    base_size = 20,
    transparent = FALSE,
    plot.margin = margin(2, 8, 4, 5, "mm")
))

# Read in the simulated datasets
aa_freq <- fread("output/aa_frequentist_estimates.csv")
aa_bayes <- fread("output/aa_bayesian_estimates.csv")

# Plot the frequentist estimates
freq_plot <- ggplot(aa_freq, aes(x = p.value)) +
    geom_histogram(bins = 30, fill = "blue", color = "black") +
    labs(
        title = "Distribution of Frequentist p-values in 3,000 A/A Tests",
        subtitle = stringr::str_wrap(
            "Distribution of p-values for the treatment condition based on 
            3,000 simulated A/A tests where the true effect is zero. p-values 
            correspond to average marginal effects from a logistic regression model.",
            width = 105
        ),
        x = "p-value for the Treatment Condition",
        y = "Frequency"
    )

# Save the plot
ggsave(
    "figures/frequentist_pvalues.jpeg", 
    freq_plot, 
    width = 12, 
    height = 10,
    dpi = "retina"
)

# Plot the Bayesian estimates
bayes_plot <- ggplot(aa_bayes, aes(x = treat_prob)) +
    geom_histogram(bins = 30, fill = "blue", color = "black") +
    labs(
        title = "Distribution of Posterior Predictive Probabilities in 3,000 A/A Tests",
        subtitle = stringr::str_wrap(
            "Distribution of treatment probabilities for the treatment condition based on 
            3,000 simulated A/A tests where the true effect is zero. Probabilities 
            correspond to the posterior predictive probabilities from a Bayesian logistic 
            regression model.",
            width = 105
        ),
        x = "Posterior Predictive Probability",
        y = "Frequency"
    )

# Save the plot
ggsave(
    "figures/bayesian_treatment_probabilities.jpeg", 
    bayes_plot, 
    width = 12, 
    height = 12,
    dpi = "retina"
)



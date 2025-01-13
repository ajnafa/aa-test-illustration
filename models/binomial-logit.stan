/* 
    Model: Binomial Likelihood with a Logit Link, Sufficient
           Formulation of the Bernoulli Likelihood
    Author: A. Jordan Nafa
    Date: 2023-10-29
    License: MIT
*/

data {
    // Data Dimensions
    int<lower=0> N;                             // N Cells
    int<lower=1> D;                             // D Experimental Conditions

    // Input Data
    array[N] int<lower=0> Y;                    // Observed Successes
    array[N] int<lower=0> K;                    // Observed Trials
    matrix[N,D] P;                              // Design Matrix

    // Priors for the Parameters
    real mu_alpha;                              // Intercept Prior Mean
    real<lower=0> sigma_alpha;                  // Intercept Prior Std Dev
    real mu_beta;                               // Coefficients Prior Mean
    real<lower=0> sigma_beta;                   // Coefficients Prior Std Dev
}

transformed data {
    int L = D - 1;                              // Number of Coefficients
    matrix[N, L] X;                             // Design Matrix

    // Design Matrix for the Coefficients
    X = P[, 2:D];

    // Total Number of Trials
    array[D] int T = rep_array(sum(K), D);
}

parameters {
    real alpha;                                 // Intercept
    vector[L] beta;                             // Coefficients
}

model {
    // Binomial Likelihood, GLM with Logit Link
    target += binomial_logit_glm_lpmf(Y | K, X, alpha, beta);

    // Priors
    target += normal_lpdf(alpha | mu_alpha, sigma_alpha);
    target += normal_lpdf(beta | mu_beta, sigma_beta);
}

generated quantities {
    // Posterior Predictive Draws Under the Counterfactuals
    array[D] int y_rep = binomial_rng(T, inv_logit(alpha + X * beta));

    // Probabilitiy Under Y(1) and Y(0)
    vector[D] theta = to_vector(y_rep) ./ to_vector(T);

    /* One could gripe here that this isn't a proper posterior probability calculation
       and may be more appropriately called a posterior predictive probability or 
       more commonly "probability of best," but it illustrates the idea and since 
       proper posterior probabilities require bridge sampling, they're a pain to 
       calculate (and AFAIK, {bridgesampling} doesn't support models fit using cmdstanr). 

       An alternative approach would be to calculate the log likelihood and use 
       LOO to produce stacking or pseudo-BMA weights then calculate probabilities 
       based on importance sampling weights, but again, this is a pain to do when 
       you're fitting the model 3,000 times.
    */

    // Posterior Predictive Probability
    simplex[D] post_prob;
    {
        real best_arm = max(theta);
        for (d in 1:D) {
            post_prob[d] = (theta[d] >= best_arm);
        }
        // Uniform in the Case of Ties
        post_prob = post_prob / sum(post_prob);
    }
}
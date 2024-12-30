# Load necessary libraries
library(ggplot2)
library(boot)
library(readxl)
library(ks)  # For kernel density estimation

# Step 1: Load the data for Case 6 and Case 7 (adjust paths as needed)
case6 <- read_excel("Mcase6 process.xls")
case7 <- read_excel("Mcase7 process.xls")

# Step 2: Extract the relevant parameters for Case 6 and Case 7
# Columns: 15(omega_v), 16(alpha), 17 (Beta_s), 18 (Beta_a), 19 (Xi_v), 20 (Epsilon_v), 25 (Delta), 34 (p_s), 35 (p_a), 27 (R0)
params_case6 <- case6[, c(15, 16, 17, 18, 19, 20, 25, 34, 35, 27)]  # Include R0 (column 27)
params_case7 <- case7[, c(15, 16, 17, 18, 19, 20, 25, 34, 35, 27)]

# Step 3: Check the number of rows in each data frame and truncate the larger one to match the smaller one
nrow_case6 <- nrow(params_case6)
nrow_case7 <- nrow(params_case7)

# Determine the smaller size
min_rows <- min(nrow_case6, nrow_case7)

# Truncate both data frames to the size of the smaller one
params_case6_truncated <- params_case6[1:min_rows, ]
params_case7_truncated <- params_case7[1:min_rows, ]

# Step 4: Calculate the parameter differences between Case 6 and Case 7 (after truncating)
diff_params <- params_case7_truncated - params_case6_truncated

# Step 4: Bootstrap Resampling Function (using parameter differences)
bootstrap_resampling <- function(diff_params, R = 1000) {
  # Define the bootstrap function to compute the mean of the resampled differences
  bootstrap_function <- function(data, indices) {
    sample_data <- data[indices, ]
    return(apply(sample_data, 2, mean, na.rm = TRUE))  # Compute the mean of resampled data
  }
  
  # Apply bootstrap resampling for the differences
  bootstrap_results <- boot(diff_params, bootstrap_function, R = R)
  return(bootstrap_results)
}

# Apply the bootstrap resampling for Case 6 to 7
bootstrap_case6_to_7 <- bootstrap_resampling(diff_params)

# Step 5: Perform Kernel Density Estimation (KDE) on bootstrapped samples
kde_results <- lapply(1:ncol(diff_params), function(i) kde(bootstrap_case6_to_7$t[, i]))

# Step 6: Visualize the KDE results for each parameter
param_names <- c("omega_v", "alpha", "beta_s", "beta_a", "xi_v", "epsilon_v", "delta", "p_s", "p_a", "R0")

# Step 7: Plot KDE results for each parameter using ggplot2
plot_kde_results <- function(kde_results, param_names) {
  par(mfrow = c(3, 3))  # Set up a 3x3 grid for the plots
  
  for (i in 1:length(kde_results)) {
    kde_obj <- kde_results[[i]]  # Get KDE object for each parameter
    data_plot <- data.frame(Value = kde_obj$eval.points, Density = kde_obj$estimate)  # Create a data frame
    
    # Plot KDE result using ggplot2
    gg <- ggplot(data_plot, aes(x = Value, y = Density)) +
      geom_line(color = "blue", size = 1) +  # Line plot for density
      labs(title = paste("KDE for", param_names[i]), 
           x = param_names[i], y = "Density") +
      theme_minimal()
    
    print(gg)  # Print the plot
  }
}

# Step 8: Call the function to visualize KDE results
plot_kde_results(kde_results, param_names)

# Step 9: Joint likelihood calculation based on KDE (using specific parameters for each hypothesis)
joint_likelihood_estimation <- function(bootstrap_samples, observed_diffs, param_indices) {
  n <- length(param_indices)
  kde_joints <- lapply(param_indices, function(i) kde(bootstrap_samples[, i]))
  
  # Calculate the joint density at the observed differences
  joint_density <- 1
  for (i in 1:n) {
    density_at_observed <- approx(kde_joints[[i]]$eval.points, kde_joints[[i]]$estimate, xout = observed_diffs[param_indices[i]])$y
    joint_density <- joint_density * density_at_observed  # Calculate the product of densities
  }
  
  return(joint_density)
}

# Observed values of parameter differences (this is the actual difference between Case 6 and Case 7)
observed_diffs <- apply(diff_params, 2, mean)

# Step 10: Define which parameters to use for each hypothesis

# Transmission-Virulence Trade-Off: Parameters p_s, p_a, delta
param_indices_H1 <- c(9, 7)

# Vaccination-Induced Virulence: Parameters omega_v, xi_v, delta
param_indices_H3 <- c(1, 7)

# Immune Selection: Parameters alpha, xi_v, delta
param_indices_H4 <- c(2, 7)

# Transmission-Virulence Correlation: Parameters p_s, p_a, delta
param_indices_H6 <- c(8, 7)

# Short-Sighted Evolution: Parameters beta_s, beta_a, delta
param_indices_H2 <- c( 4, 7)

# Step 11: Calculate the likelihood for each hypothesis using the specific parameters

# Likelihood for Hypothesis 1: Transmission-Virulence Trade-Off
likelihood_H1 <- joint_likelihood_estimation(bootstrap_case6_to_7$t, observed_diffs, param_indices_H1)

# Likelihood for Hypothesis 2: Short-Sighted Evolution
likelihood_H2 <- joint_likelihood_estimation(bootstrap_case6_to_7$t, observed_diffs, param_indices_H2)

# Likelihood for Hypothesis 3: Vaccination-Induced Virulence
likelihood_H3 <- joint_likelihood_estimation(bootstrap_case6_to_7$t, observed_diffs, param_indices_H3)

# Likelihood for Hypothesis 4: Immune Selection
likelihood_H4 <- joint_likelihood_estimation(bootstrap_case6_to_7$t, observed_diffs, param_indices_H4)

# Likelihood for Hypothesis 6: Transmission-Virulence Correlation
#likelihood_H6 <- joint_likelihood_estimation(bootstrap_case6_to_7$t, observed_diffs, param_indices_H6)

# Step 12: Define Prior Probabilities for Each Hypothesis (Uniform prior assumed)
prior_H1 <- 1/4  # Transmission-Virulence Trade-Off
prior_H2 <- 1/4  # Short-Sighted Evolution
prior_H3 <- 1/4  # Vaccination-Induced Virulence
prior_H4 <- 1/4  # Immune Selection
#prior_H6 <- 1/5  # Transmission-Virulence Correlation

# Step 13: Compute Posterior Probabilities Using Bayes' Theorem
total_likelihood <- likelihood_H1 * prior_H1 + likelihood_H2 * prior_H2 + likelihood_H3 * prior_H3 + likelihood_H4 * prior_H4 

# Avoid division by zero by setting a small number if total likelihood is 0
if (total_likelihood == 0) {
  total_likelihood <- 1e-10  # Small value to prevent division by zero
}

posterior_H1 <- (likelihood_H1 * prior_H1) / total_likelihood
posterior_H2 <- (likelihood_H2 * prior_H2) / total_likelihood
posterior_H3 <- (likelihood_H3 * prior_H3) / total_likelihood
posterior_H4 <- (likelihood_H4 * prior_H4) / total_likelihood
#posterior_H6 <- (likelihood_H6 * prior_H6) / total_likelihood

# Step 14: Output the Posterior Probabilities
cat("Posterior probability of H1 (Transmission-Virulence Trade-Off):", posterior_H1, "\n")
cat("Posterior probability of H2 (Short-Sighted Evolution):", posterior_H2, "\n")
cat("Posterior probability of H3 (Vaccination-Induced Virulence):", posterior_H3, "\n")
cat("Posterior probability of H4 (Immune Selection):", posterior_H4, "\n")
#cat("Posterior probability of H6 (Transmission-Virulence Correlation):", posterior_H6, "\n")

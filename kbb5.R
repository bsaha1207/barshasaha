# Load necessary libraries
library(ggplot2)
library(boot)
library(readxl)
library(ks)  # For kernel density estimation

# Step 1: Load the data for Case 5 and Case 6 (adjust paths as needed)
case5 <- read_excel("Mcase5 process.xls")
case6 <- read_excel("Mcase6 process.xls")

# Step 2: Extract the relevant parameters for Case 5 and Case 6
# Columns: 15(omega_v), 16(alpha), 17 (Beta_s), 18 (Beta_a), 19 (Xi_v), 20 (Epsilon_v), 25 (Delta), 34 (p_s), 35 (p_a), 27 (R0)
params_case5 <- case5[, c(15, 16, 17, 18, 19, 20, 25, 34, 35, 27)]  # Include R0 (column 27)
params_case6 <- case6[, c(15, 16, 17, 18, 19, 20, 25, 34, 35, 27)]

# Step 3: Check the number of rows in each data frame and truncate the larger one to match the smaller one
nrow_case5 <- nrow(params_case5)
nrow_case6 <- nrow(params_case6)

# Determine the smaller size
min_rows <- min(nrow_case5, nrow_case6)

# Truncate both data frames to the size of the smaller one
params_case5_truncated <- params_case5[1:min_rows, ]
params_case6_truncated <- params_case6[1:min_rows, ]

# Step 4: Calculate the parameter differences between Case 5 and Case 6 (after truncating)
diff_params <- params_case6_truncated - params_case5_truncated

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

# Apply the bootstrap resampling for Case 5 to 6
bootstrap_case5_to_6 <- bootstrap_resampling(diff_params)

# Step 5: Perform Kernel Density Estimation (KDE) on bootstrapped samples
kde_results <- lapply(1:ncol(diff_params), function(i) kde(bootstrap_case5_to_6$t[, i]))

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

# Observed values of parameter differences (this is the actual difference between Case 5 and Case 6)
observed_diffs <- apply(diff_params, 2, mean)

# Step 10: Define which parameters to use for each hypothesis
# Vaccination-Induced Virulence: Parameters xi_v, epsilon_v, delta
param_indices_H3 <- c(5, 6, 7)

# Law of Declining Virulence: Parameters beta_s, beta_a, delta
param_indices_H5 <- c(7, 8, 9)

# Transmission-Virulence Correlation: Parameters p_s, p_a, delta
param_indices_H6 <- c(3, 4, 7)

# Step 11: Calculate the likelihood for each hypothesis using the specific parameters

# Likelihood for Hypothesis 3: Vaccination-Induced Virulence
likelihood_H3 <- joint_likelihood_estimation(bootstrap_case5_to_6$t, observed_diffs, param_indices_H3)

# Likelihood for Hypothesis 5: Law of Declining Virulence
likelihood_H5 <- joint_likelihood_estimation(bootstrap_case5_to_6$t, observed_diffs, param_indices_H5)

# Likelihood for Hypothesis 6: Transmission-Virulence Correlation
likelihood_H6 <- joint_likelihood_estimation(bootstrap_case5_to_6$t, observed_diffs, param_indices_H6)

# Step 12: Define Prior Probabilities for Each Hypothesis (Uniform prior assumed)
prior_H3 <- 1/3  # Vaccine-Induced
prior_H5 <- 1/3  # Law of Declining Virulence
prior_H6 <- 1/3  # Transmission-Virulence Correlation

# Step 13: Compute Posterior Probabilities Using Bayes' Theorem
total_likelihood <- likelihood_H3 * prior_H3 + likelihood_H5 * prior_H5+ likelihood_H6 * prior_H6

# Avoid division by zero by setting a small number if total likelihood is 0
if (total_likelihood == 0) {
  total_likelihood <- 1e-10  # Small value to prevent division by zero
}

posterior_H3 <- (likelihood_H3 * prior_H3) / total_likelihood
posterior_H5 <- (likelihood_H5 * prior_H5) / total_likelihood
posterior_H6 <- (likelihood_H6 * prior_H6) / total_likelihood

# Step 14: Output the Posterior Probabilities
cat("Posterior probability of H3 (Vaccination-Induced Virulence):", posterior_H3, "\n")
cat("Posterior probability of H5 (Law of Declining Virulence):", posterior_H5, "\n")
cat("Posterior probability of H6 (Transmission-Virulence Correlation):", posterior_H6, "\n")


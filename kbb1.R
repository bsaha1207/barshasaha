# Install necessary packages if not already installed
# install.packages("ggplot2")

# Load necessary libraries
library(ggplot2)
library(boot)
library(readxl)
library(ks)  # For kernel density estimation

# Step 1: Load the data for Case 1 and Case 2 (adjust paths as needed)
case1 <- read_excel("Mcase1 process.xls")
case2 <- read_excel("Mcase2 process.xls")

# Step 2: Extract the relevant parameters for Case 1 and Case 2
# Columns: 15(omega_v), 16(alpha), 17 (Beta_s), 18 (Beta_a), 19 (Xi_v), 20 (Epsilon_v), 25 (Delta), 34 (p_s), 35 (p_a), 27 (R0)
params_case1 <- case1[, c( 16, 17, 18, 25, 34, 35, 27)]  # Include R0 (column 27)
params_case2 <- case2[, c( 16, 17, 18,  25, 34, 35, 27)]

# Step 3: Check the number of rows in each data frame and truncate the larger one to match the smaller one
nrow_case1 <- nrow(params_case1)
nrow_case2 <- nrow(params_case2)

# Determine the smaller size
min_rows <- min(nrow_case1, nrow_case2)

# Truncate both data frames to the size of the smaller one
params_case1_truncated <- params_case1[1:min_rows, ]
params_case2_truncated <- params_case2[1:min_rows, ]

# Step 4: Calculate the parameter differences between Case 4 and Case 5 (after truncating)
diff_params <- params_case2_truncated - params_case1_truncated

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

# Apply the bootstrap resampling for Case 1 to 2
bootstrap_case1_to_2 <- bootstrap_resampling(diff_params)

# Step 5: Perform Kernel Density Estimation (KDE) on bootstrapped samples
kde_results <- lapply(1:ncol(diff_params), function(i) kde(bootstrap_case1_to_2$t[, i]))

# Step 6: Visualize the KDE results for each parameter
param_names <- c( "alpha", "beta_s", "beta_a", "delta", "p_s", "p_a", "R0")

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

#Step 9:Joint likelihood calculation based on KDE (using all parameters together, including R0)
joint_likelihood_estimation <- function(bootstrap_samples, observed_diffs) {
  n <- ncol(bootstrap_samples)
  kde_joints <- lapply(1:n, function(i) kde(bootstrap_samples[, i]))
  
  # Calculate the joint density at the observed differences
  joint_density <- 1
  for (i in 1:n) {
    density_at_observed <- approx(kde_joints[[i]]$eval.points, kde_joints[[i]]$estimate, xout = observed_diffs[i])$y
    joint_density <- joint_density * density_at_observed  # Calculate the product of densities
  }
  
  return(joint_density)
}

# Observed values of parameter differences (this is the actual difference between Case 1 and Case 2)
observed_diffs <- apply(diff_params, 2, mean)

# Step 10: Use the R0 change to differentiate hypotheses
R0_diff <- observed_diffs[7]  # R0 is in the 10th column of diff_params

# Set thresholds based on R0 change to differentiate hypotheses:
# We expect small/moderate change for H1, a larger decrease for H2, and mixed behavior for H3.

threshold_tradeoff_R0 <- -0.5  # Small decrease for H1
threshold_shortsighted_R0 <- -1.0  # Larger decrease for H2

# Step 11: Calculate the likelihood for each hypothesis

# Likelihood for Hypothesis 1: Transmission-Virulence Trade-Off
# We expect a small/moderate decrease in R0 and balanced changes in transmission and virulence.
likelihood_H1 <- ifelse(R0_diff < threshold_tradeoff_R0, 
                        joint_likelihood_estimation(bootstrap_case1_to_2$t[, c(4,5,6)], observed_diffs[c(4,5,6)]), 0)

# Likelihood for Hypothesis 2: Short-Sighted Evolution
# We expect a larger decrease in R0 (since virulence increases more aggressively, harming long-term transmission).
likelihood_H2 <- ifelse(R0_diff < threshold_shortsighted_R0, 
                        joint_likelihood_estimation(bootstrap_case1_to_2$t[, c(2,4,5,6)], observed_diffs[c(2,4,5,6)]), 0)

# Likelihood for Hypothesis 3: Vaccination-Induced Virulence
# R0 can either decrease or increase, but this depends on vaccination parameters and virulence.
#likelihood_H3 <- joint_likelihood_estimation(bootstrap_case1_to_2$t[, c(4, 5, 6)], observed_diffs[c(4, 5, 6)])

# Step 12: Define Prior Probabilities for Each Hypothesis (Uniform prior assumed)
prior_H1 <- 1/2  # Transmission-Virulence Trade-Off
prior_H2 <- 1/2  # Short-Sighted Evolution
#prior_H3 <- 1/3  # Vaccination-Induced Virulence

# Step 13: Compute Posterior Probabilities Using Bayes' Theorem
total_likelihood <- likelihood_H1 * prior_H1 + likelihood_H2 * prior_H2 

# Avoid division by zero by setting a small number if total likelihood is 0
if (total_likelihood == 0) {
  total_likelihood <- 1e-10  # Small value to prevent division by zero
}

posterior_H1 <- (likelihood_H1 * prior_H1) / total_likelihood
posterior_H2 <- (likelihood_H2 * prior_H2) / total_likelihood
#posterior_H3 <- (likelihood_H3 * prior_H3) / total_likelihood

# Step 14: Output the Posterior Probabilities
cat("Posterior probability of H1 (Transmission-Virulence Trade-Off):", posterior_H1, "\n")
cat("Posterior probability of H2 (Short-Sighted Evolution):", posterior_H2, "\n")
#cat("Posterior probability of H3 (Vaccination-Induced Virulence):", posterior_H3, "\n")

# Install plotly if not already installed
# install.packages("plotly")

library(plotly)

# Extract the bootstrapped samples for delta, p_s, and p_a
delta_samples <- bootstrap_case1_to_2$t[, 4]
p_s_samples <- bootstrap_case1_to_2$t[, 5]
p_a_samples <- bootstrap_case1_to_2$t[, 6]

# Create a 3D scatter plot of the joint density
plot_ly(x = ~delta_samples, y = ~p_s_samples, z = ~p_a_samples, type = "scatter3d", mode = "markers",
        marker = list(size = 3, color = 'blue', opacity = 0.6)) %>%
  layout(title = "3D Joint Density of Delta, p_s, and p_a",
         scene = list(xaxis = list(title = 'Delta'),
                      yaxis = list(title = 'p_s'),
                      zaxis = list(title = 'p_a')))

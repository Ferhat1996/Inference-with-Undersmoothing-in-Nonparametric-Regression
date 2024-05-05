
#**********************************Exercise_2_a*********************************#
# Load necessary libraries
#install.packages("locpol")
#install.packages("stats")

library(locpol)
library(stats)

# Define the true function m(x)
m <- function(x) { x^2 }

# Define the data-generating process
generate_data <- function(n) {
  X <- runif(n, -1, 1)
  u <- rnorm(n, 0, sqrt(0.25))
  Y <- m(X) + u
  data.frame(X = X, Y = Y)
}

# From Lecture Notes: V(x) ∼ σ^2(x)/f(x)

# triangular_kernel <- function(v) {return((1 - abs(data)) * (abs(data) <= 1)) can be used but it gives same results}
# kappa_0 = 2/3, mu_2 = 1 / 6

f_x <- 1/2 # For triangular kernel

# Define the function to compute the local linear estimator
local_linear_estimator <- function(data, x, h) {
  fit <- locpoly(data$X, data$Y, bandwidth = h, kernel = "TriweigK", degree = 1)
  approx(fit$x, fit$y, x)$y
}

# Define the function to compute the empirical coverage rate
compute_coverage_rate <- function(n, S, alpha) {
  B <- 2  # Second derivative of m(x)
  
  # Generate data once to calculate the variance
  data <- generate_data(n)
  V_LL <- var(local_linear_estimator(data, data$X, h = 0.1)) / f_x # a reasonable bandwidth (e.g., 0.1) here
  h_star <- ((V_LL)/(4 * B^2))^(1/5) * n^(-1/5)  # AMSE-optimal bandwidth
  
  
  cat(paste("Sample size:", n, "\n",
            "h_star:", h_star, "\n",
            "V(x):", V_LL, "\n"))
  
  c_alpha <- qnorm(1 - alpha/2)  # 1 - alpha/2 quantile of the normal distribution
  
  inside_interval <- replicate(S, {
    data <- generate_data(n)
    m_hat <- local_linear_estimator(data, 0, h_star)
    V_LL_hat <- var(local_linear_estimator(data, data$X, h = h_star)) / f_x 
    CI <- c(m_hat - c_alpha * sqrt(V_LL_hat/(n*h_star)), m_hat + c_alpha * sqrt(V_LL_hat/(n*h_star)))
    m(0) >= CI[1] && m(0) <= CI[2]
  })
  
  mean(inside_interval)
}

# Define the function to compute the rescaled bias
compute_rescaled_bias <- function(n, S) {
  B <- 2  # Second derivative of m(x)
  
  # Generate data once to calculate the variance
  data <- generate_data(n)
  V_LL <- var(local_linear_estimator(data, data$X, h = 0.1)) / f_x # a reasonable bandwidth (e.g., 0.1) here
  h_star <- ((V_LL)/(4 * B^2))^(1/5) * n^(-1/5)  # AMSE-optimal bandwidth
  
  biases <- replicate(S, {
    data <- generate_data(n)
    m_hat <- local_linear_estimator(data, 0, h_star)
    V_LL_hat <- var(local_linear_estimator(data, data$X, h = h_star)) / f_x 
    sqrt(n * h_star) * (m_hat - m(0))
  })
  
  mean(biases)
}

# Perform the Monte Carlo study and store results in a data frame
results <- data.frame(
  Sample_Size = numeric(),
  Empirical_Coverage_Rate = numeric(),
  Rescaled_Bias = numeric()
)

set.seed(123)  # For reproducibility
S <- 10000  # Number of simulations
alpha <- 0.05  # Significance level

for (n in c(100, 500, 5000)) {
  data <- generate_data(n)  # Ensure data is generated before using in functions
  coverage_rate <- compute_coverage_rate(n, S, alpha)
  rescaled_bias <- compute_rescaled_bias(n, S)
  
  cat(paste("Sample size:", n, "\n",
            "Empirical coverage rate of CI:", coverage_rate, "\n",
            "Average value of the rescaled bias:", rescaled_bias, "\n\n"))
  
  # Append results to the data frame
  results <- rbind(results, c(n, coverage_rate, rescaled_bias))
}

# Add column headers
colnames(results) <- c("Sample_Size", "Empirical_Coverage_Rate", "Rescaled_Bias")

# Print the results as a table
print(results)
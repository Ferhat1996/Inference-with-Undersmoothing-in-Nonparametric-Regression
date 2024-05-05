I attempted to solve this exercise using Table 1, but I was unsuccessful. Consequently, I pursued the alternative approach that I loaded.

Below, you'll find my progress on solving the problem using the information provided in Table 1:

# Load necessary libraries
library(MASS)
library(KernSmooth)
library(numDeriv)
library(np)

# Function to generate mixture distribution
generate_mixture_distribution <- function(n) {
  set.seed(42)
  s <- rbinom(n, 1, 0.5)
  X <- s * rnorm(n, mean = -1, sd = 1) + (1 - s) * rnorm(n, mean = 1.75, sd = 0.25)
  return(X)
}


# Define the true function m(x)
m <- function(x) sin(2.5 * x)

# Define the first derivative of m(x)
m_prime <- function(x) {
  return(2.5 * cos(2.5 * x))
}

# Define the second derivative of m(x)
m_dprime <- function(x) {
  return(-2.5^2 * sin(2.5 * x))
}

# Kernel function
K <- function(u) {
  return(exp(-0.5 * u^2) / sqrt(2 * pi))
}


# Function for marginal density f_X(x)
f_x <- function(X) {
  density_result <- density(X)  # Using default parameters, you may adjust them
  f_x <- approxfun(density_result$x, density_result$y, rule = 2)  # Linear interpolation
  return(f_x)
}

# Kernel function
K <- function(u) {
  return(exp(-0.5 * u^2) / sqrt(2 * pi))
}

u <- seq(-2, 2, length.out = 1000)

# Define the Nadaraya-Watson expressions as functions
B_NW <- function(x, m_dprime, m, f_x, K, h_n) {
  u <- seq(-2, 2, length.out = 1000)
  term1 <- (1/2) * m_dprime(x)
  term2 <- D(m(x), "x") * D(f_x(x), "x") / f_x(x) * integrate(u^2 * K(u), -Inf, Inf) * h_n^2
  result <- term1 + term2
  return(result)
}

V_NW <- function(x, sigma, f_x, n, h_n, K) {
  u <- seq(-2, 2, length.out = 1000)
  result <- sigma^2(x) / (f_x(x) * n * h_n) * integrate(K(u)^2, -Inf, Inf)
  return(result)
}

# Define the Local Linear Smoother expressions as functions
B_LL <- function(x, m_dprime, K, h_n) {
  u <- seq(-2, 2, length.out = 1000)
  result <- (1/2) * m_dprime(x) * integrate(u^2 * K(u), -Inf, Inf) * h_n^2
  return(result)
}

V_LL <- function(x, sigma, f_x, n, h_n, K) {
  u <- seq(-2, 2, length.out = 1000)
  result <- sigma^2(x) / (f_x(x) * n * h_n) * integrate(K(u)^2, -Inf, Inf)
  return(result)
}

# Define a sample value for x (replace this with your actual data)
x_sample <- 1

# Calculate h_LL and h_NW
h_LL <- ((V_LL)/(4 * B_LL(x_sample, m_dprime, K, h_LL)^2))^(1/5) * n^(-1/5) 
h_NW <- ((V_NW)/(4 * B_NW(x_sample, m_dprime, m, f_x, K, h_NW)^2))^(1/5) * n^(-1/5)



asymptotic_MISE_LL <- B_LL*h_LL^4 +V_LL/(n*h_LL)
asymptotic_MISE_BW <- B_NW*h_NW^4 +V_BW/(n*h_NW)
  
# Function to estimate the regression function using local linear and Nadaraya-Watson estimators
estimate_regression <- function(data, x_grid, h_LL, h_NW) {
  # Local Linear Estimator
  m_LL <- apply(data, 1, function(x) B_LL(x, m_dprime, K, h_LL))
  m_LL <- m_LL + V_LL(data$X, sigma, f_x, nrow(data), h_LL, K)
  
  # Nadaraya-Watson Estimator
  m_NW <- apply(data, 1, function(x) B_NW(x, m_dprime, m, f_x, K, h_NW))
  m_NW <- m_NW + V_NW(data$X, sigma, f_x, nrow(data), h_NW, K)
  
  return(list(m_LL = m_LL, m_NW = m_NW))
}

# Function to simulate and calculate MISE for both local linear and Nadaraya-Watson estimators
simulate_and_calculate_MISE <- function(n, B, x_grid, h_LL, h_NW) {
  MISE <- matrix(0, nrow = 2, ncol = length(n), dimnames = list(c("LL", "NW"), n))
  
  for (i in 1:length(n)) {
    for (b in 1:B) {
      data <- generate_sample(n[i])
      estimators <- estimate_regression(data, x_grid, h_LL, h_NW)
      
      ISE_ll <- mean((m_true(x_grid) - estimators$m_LL)^2)
      ISE_nw <- mean((m_true(x_grid) - estimators$m_NW)^2)
      
      MISE[1, i] <- MISE[1, i] + ISE_ll
      MISE[2, i] <- MISE[2, i] + ISE_nw
    }
    MISE[, i] <- MISE[, i] / B
  }
  
  return(MISE)
}

# Set the sample sizes for simulation
sample_sizes <- c(100, 200, 400)

# Set the number of repetitions
repetitions <- 300

# Set initial values for bandwidths
initial_h_LL <- 1
initial_h_NW <- 1

lower_bound <- -2
upper_bound <- 2
# Run the optimization to find optimal bandwidths
optimal_result_LL <- optimize(asymptotic_MISE_LL, interval = c(lower_bound, upper_bound),
                              B_LL = B_LL, f_x = f_x, sigma = sigma, X = X)
optimal_result_NW <- optimize(asymptotic_MISE_NW, interval = c(lower_bound, upper_bound),
                              B_NW = B_NW, f_x = f_x, sigma = sigma, X = X)

# Extract optimal bandwidths
h_star_LL <- optimal_result_LL$minimum
h_star_NW <- optimal_result_NW$minimum

# Run the simulation study
MISE_results <- simulate_and_calculate_MISE(sample_sizes, repetitions, x_grid, h_star_LL, h_star_NW)

# Compute the relative efficiency of the local linear estimator
RE <- (MISE_results[1, ] / MISE_results[2, ])^5/4

# Construct Table 4
table_data <- data.frame(Sample_Size = sample_sizes, MISE_Local_Linear = MISE_results[1, ], MISE_Nadaraya_Watson = MISE_results[2, ], Relative_Efficiency = RE)
print("Table 4:")
print(table_data)

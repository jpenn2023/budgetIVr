# Load libraries, register cores
library(data.table)
# 
# library(doParallel)
# cl <- makeCluster(6)
# registerDoParallel(cl)

setwd('C:/Users/k23067841/Downloads/BudgetIV_data')

library(mvtnorm)
library(matrixcalc)

# Set seed (first attempt)
set.seed(987)

# Set seed (second attempt)
set.seed(654)

random_simulation <- function(A, B, theta_true, N, sim_idx){
  
  random_parameters <- randomise_parameters_heteroskedastic(A, B, theta_true)
  
  sigma_mat <- matrix(0, 4, 4)
  
  # On-diagonal variances
  diag(sigma_mat) <- c((random_parameters$Z1)^2, (random_parameters$Z2)^2, (random_parameters$ex)^2, (random_parameters$ea)^2)
  
  # Off-diagonal cross covariances
  sigma_mat[upper.tri(sigma_mat)] <- c(random_parameters$Z1Z2 * random_parameters$Z1 * random_parameters$Z2,
                                       random_parameters$Z1ex * random_parameters$Z1 * random_parameters$ex, 
                                       random_parameters$Z2ex * random_parameters$Z2 * random_parameters$ex, 
                                       random_parameters$Z1ea * random_parameters$Z1 * random_parameters$ea,
                                       random_parameters$Z2ea * random_parameters$Z2 * random_parameters$ea,
                                       random_parameters$exea * random_parameters$ex * random_parameters$ea)
  
  sigma_mat[lower.tri(sigma_mat)] <- t(sigma_mat)[lower.tri(sigma_mat)]
  
  
  # Rejection sampling to get a valid covariance matrix (i.e., positive semidefinite)
  while (random_parameters$ea < 0 || !is.positive.semi.definite(sigma_mat)){
    
    random_parameters <- randomise_parameters_heteroskedastic(A, B, theta_true)
  
    sigma_mat <- matrix(0, 4, 4)
    
    # On-diagonal variances
    diag(sigma_mat) <- c((random_parameters$Z1)^2, (random_parameters$Z2)^2, (random_parameters$ex)^2, (random_parameters$ea)^2)
    
    # Off-diagonal cross covariances
    sigma_mat[upper.tri(sigma_mat)] <- c(random_parameters$Z1Z2 * random_parameters$Z1 * random_parameters$Z2,
                                         random_parameters$Z1ex * random_parameters$Z1 * random_parameters$ex, 
                                         random_parameters$Z2ex * random_parameters$Z2 * random_parameters$ex, 
                                         random_parameters$Z1ea * random_parameters$Z1 * random_parameters$ea,
                                         random_parameters$Z2ea * random_parameters$Z2 * random_parameters$ea,
                                         random_parameters$exea * random_parameters$ex * random_parameters$ea)
    
    sigma_mat[lower.tri(sigma_mat)] <- t(sigma_mat)[lower.tri(sigma_mat)]
    
    }
  
  dataset <- simulate_data(random_parameters, sigma_mat, N)
  fwrite(dataset, paste0('./real deal 2/', sim_idx, '.csv'))
  
  
}

randomise_parameters_heteroskedastic <- function(A, B, theta_true){
  
  latent_cov <- list("Z1gy" = (A[1] - B[1] * theta_true), "Z2gy" = (A[2] - B[2] * theta_true))
  
  confounding_coeffs <- list("Z1Z2" = NA, "Z1ex" = NA, "Z1ea" = NA, 
                             "Z2ex" = NA, "Z2ea" = NA, "exea" = NA)
  
  #variances <- list("Z1" = NA, "Z2" = NA, "ex" = NA, "ea" = NA)
  std_devs <- list("Z1" = NA, "Z2" = NA, "ex" = NA, "ea" = NA)
  
  mu_m <- list("mu_m" = NA)
  
  # Generate rho_{Z_1 e_x}
  if (B[1] > 0){ confounding_coeffs$Z1ex <- runif(1) }
  else if (B[1] < 0) { confounding_coeffs$Z1ex <- runif(1, -1, 0) }
  
  # Generate rho_{Z_2 e_x}
  if (B[2] > 0) { confounding_coeffs$Z2ex <- runif(1) }
  else if (B[2] < 0) { confounding_coeffs$Z2ex <- runif(1, -1, 0)}
  
  # Generate other confounding coefficients
  confounding_coeffs$Z1Z2 <- runif(1, -1, 1)
  confounding_coeffs$Z1ea <- runif(1, -1, 1)
  confounding_coeffs$Z2ea <- runif(1, -1, 1)
  confounding_coeffs$exea <- runif(1, -1, 1)
  
  # Randomize free standard deviation 
  std_devs$ex <- rexp(1)
  
  # Calculate standard deviations of Z_1, Z_2, now constrained
  std_devs$Z1 <- B[1]/(std_devs$ex*confounding_coeffs$Z1ex)
  std_devs$Z2 <- B[2]/(std_devs$ex*confounding_coeffs$Z2ex)
  
  # Mean of e_m (see appendix for calculation)
  mu_m$mu_m <- (1/((std_devs$Z1)^2 * (std_devs$Z2))) * ((confounding_coeffs$Z2ea * std_devs$Z2 * latent_cov$Z1gy - 
                                                      confounding_coeffs$Z1ea * std_devs$Z1 * latent_cov$Z2gy)/(
                                                        confounding_coeffs$Z2ea - confounding_coeffs$Z1ea * confounding_coeffs$Z1Z2
                                                      ))
  
  # Standard deviation of e_a (")
  std_devs$ea <- (1/(std_devs$Z1 * std_devs$Z2)) * (confounding_coeffs$Z1Z2 * std_devs$Z2 * latent_cov$Z1gy - std_devs$Z1 * latent_cov$Z2gy)
  
  return(c(confounding_coeffs, std_devs, mu_m))
  
}

simulate_data <- function(random_parameters, sigma_mat, N) {
  
  # Generate data from e_m
  
  if(random_parameters$mu_m < 0){
    em_dat <- runif(N, random_parameters$mu_m, 0)
  }
  
  else {
    em_dat <- runif(N, 0, random_parameters$mu_m)
  }
  
  my_data <- rmvnorm(N, sigma = sigma_mat)
  
  Z1_dat <- my_data[, 1]
  Z2_dat <- my_data[, 2]
  ex_dat <- my_data[, 3]
  ea_dat <- my_data[, 4]
  
  X_dat <- ex_dat
  Y_dat <- theta_true * X_dat + em_dat * Z1_dat + ea_dat
  
  return(data.table(Z1_dat, Z2_dat, X_dat, Y_dat))
  
  }

A <- c(-4, 1)
B <- c(-2, 0.9)
theta_true <- 1
N <- 1000000
unique_settings = 30

for (sim_idx in 1:unique_settings){
  
  random_simulation(A, B, theta_true, N, sim_idx)
  
}

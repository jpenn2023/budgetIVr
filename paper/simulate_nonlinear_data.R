# Code used to generate the dataset simulated_data_BudgetIV.RData in the BudgetIV package.
# 
# Here, the dataset is saved as a csv file called "my_dat R = 0.5 SNR_y = 1.csv".
# 
# The approach is outlined in appendix Appx. C.2 of Penn et al. (2025).
#
# Jordan Penn, Lee Gunderson, Gecia Bravo-Hermsdorff,
# Ricardo Silva, and David Watson. (2024). BudgetIV: Optimal Partial Identification of Causal Effects with Mostly Invalid Instruments. \emph{arXiv}
# preprint, 2411.06913.

library(MASS)
library(ggplot2)

devtools::install_github('jpenn2023/BudgetIV')

# Set seed
set.seed(123)

# Hyperparameters
d_z <- 6L
n <- 500000L
num_valid <- 3
# Simulate normal m's
m <- rnorm(d_z, mean = 1, sd=4)

# Simulate Lambda matrix
Lambda <- matrix(nrow = d_z, ncol = d_z)
for (i in 1:d_z) {
  for (j in 1:d_z) {
    if (i == j & i > num_valid) {
      # print("diag")
      # Lambda[i, j] <- rbinom(1, 1, prob = 0.9)
      Lambda[i, j] <- rnorm(1, mean=0.9, sd=0.9)
    } else if (i != j & i > num_valid & j > num_valid) {
      # print("off diag")
      # Lambda[i, j] <- rbinom(1, 1, prob = 0.3)
      Lambda[i, j] <- rnorm(1, mean=0.3, sd=0.3)
    } else {
      Lambda[i, j] <- 0L
    }
  }
}

# Z's form a Markov chain
Z <- matrix(nrow = n, ncol = d_z)
Z[, 1] <- rbinom(n, 1, prob = 0.05)
for (j in 2:3) {
  Z[, j] <- rbinom(n, 1, prob = ifelse(Z[, (j - 1)] == 1, 0.9, 0.05))
}
Z[, 4] <- rbinom(n, 1, prob = 0.05)
for (j in 5:d_z){
  Z[, j] <- rbinom(n, 1, prob = ifelse(Z[, (j - 1)] == 1, 0.9, 0.05))
}

# Z <- matrix(nrow = n, ncol = d_z)
# Z[, 1] <- rbinom(n, 1, prob = 0.05)
# Z[, 2] <- rbinom(n, 1, prob = ifelse(Z[, 1] == 1, 0.9, 0.05))
# Z[, 3] <- rbinom(n, 1, prob = 0.05)
# for (j in 4:d_z){
#   Z[, j] <- rbinom(n, 1, prob = ifelse(Z[, (j - 1)] == 1, 0.9, 0.05))
# }

# print(Z)

# Induce correlation with eps_x, eps_y
x_wts <- c(rnorm(d_z)) # Consider changing sd
y_wts <- c(0, 0, 0, rnorm(d_z - 3, mean=1)) # Consider changing sd
eps_noise <- rnorm(n)               # Consider changing sd
eps_x <- drop(Z %*% x_wts) + eps_noise
eps_y <- drop(Z %*% y_wts) + eps_noise

# Now we have the structural equation for x
X <- drop(plogis(eps_x - Z %*% m))

# phi <- sqrt(1 + (X - 0.5)^2)
phi <- (X - 0.25)^2 - 0.25^2
theta <- matrix(c(-0.25,1), nrow=1)

# print(Z[1,] %*% Lambda %*% Z[i, ])

# Varying simulation settings
# phi <- sqrt(1 + (X - 0.5)^2)
var_y <- 5
gamma <- sapply(1:n, function(i) {
  t(Z[i, ]) %*% Lambda %*% Z[i, ]
})
Z_tilde <- Z[, num_valid:d_z]
numer <- norm(cov(gamma, Z_tilde), type = '2') 
denom <- norm(cov(eps_y, Z_tilde), type = '2') 
r <- numer / denom
sim_grd <- as.data.table(
  expand.grid(snr_y = c(1/10, 1/5, 1), 
              R = c(1/2, 1, 2))
)
sim_grd[, s := sqrt((1 / var(phi)) * (var_y / (1 + 1 / snr_y)))]
sim_grd[, idx := .I]

# u, lambda_a2, lambda_a3 depend on the simulation setting
for (i in 1:nrow(sim_grd)) {
  tmp <- sim_grd[i, ]
  u <- r / tmp$R * gamma + eps_y
  lambda_a2 <- ((tmp$s * cov(phi, u)) / var(u)) * 
    (-1 + sqrt(1 + (var(u) / (tmp$s * cov(phi, u))^2) * (var_y / (1 + tmp$snr_y))))
  lambda_a3 <- tmp$R / r * lambda_a2
  
  Z_cent <- scale(Z, center = TRUE, scale = FALSE)
  Y_cent <- scale(tmp$s * phi + lambda_a3 * gamma + lambda_a2 * eps_y, center = TRUE, scale = FALSE)
  phi_cent <- scale(phi, center = TRUE, scale = FALSE)
  
  # print(nrow(Z_cent))
  # print(nrow(Y_cent))
  # print(nrow(phi_cent))
  
  beta_y <- t(t(Z_cent) %*% Y_cent) / (nrow(Z_cent) - 1)
  
  # Basis expansion in phi
  phi_basis_cent <- tmp$s * scale(cbind(X, X^2), center = TRUE, scale = FALSE)
  
  phi_basis <- expression(x, x^2)
  
  phi_basis_list <- list(expression(x), expression(x^2))
  
  phi_basis_list <- list("x", "x^2")
  
  if(is.list(phi_basis_list)){
    
    if(all(sapply(phi_basis_list, function(x) is.expression(x)))){ 
      
      phi_basis_list <- do.call(expression, phi_basis_list)
      
    }
  }
  
  
  beta_phi <- t(t(Z_cent) %*% phi_basis_cent) / (nrow(Z_cent) - 1)
  
  print(paste("s_phi:", tmp$s))
  print(paste("lambda_a2:", lambda_a2))
  print(paste("lambda_a3:", lambda_a3))
  
  var_y <- var(Y_cent)
  
  var_z <- as.vector(diag(cov(Z_cent)))
  
  SE_beta_y <- sqrt((var_y*var_z + as.vector(beta_y * beta_y))/nrow(Z_cent))
  
  alpha = 0.05
  delta_beta_y <- qnorm(1 - alpha/(2*nrow(Z_cent)))*SE_beta_y
  
  my_dat <- data.table("beta_y" = as.vector(beta_y),
                       "beta_phi_1" = t(beta_phi)[1,],
                       "beta_phi_2" = t(beta_phi)[2,],
                       "delta_beta_y" = delta_beta_y
  )
  
  print(my_dat)
  
  fwrite(my_dat, paste0("my_dat ", 'R = ', tmp$R, ' SNR_y = ', tmp$snr_y, '.csv'))
  
  #
  #
  #
  # Remaining code to run BudgetIV on each example and generate plots of results.
  #
  # If BudgetIV is not installed, install now using
  # 
  # devtools::install_github('jpenn2023/BudgetIV')
  #
  
  #
  #
  # tau_vec = c(0)
  # m_vec = c(num_valid)
  # dom_ATE <- matrix(c(0, 1), ncol=2, nrow=1)
  # 
  # x_vals <- seq(from = 0, to = 1, length.out = 500)
  # 
  # ATE_search_domain <- expand.grid("x" = x_vals)
  # 
  # X_baseline <- list("x" = c(0))
  # 
  # # print(beta_phi)
  # 
  # partial_identification_ATE <- BudgetIV(beta_y=beta_y, 
  #                                        beta_phi=beta_phi, 
  #                                        phi_basis=phi_basis, 
  #                                        tau_vec=tau_vec, 
  #                                        b_vec=m_vec, 
  #                                        ATE_search_domain=ATE_search_domain, 
  #                                        X_baseline=X_baseline, delta_beta_y = as.vector(delta_beta_y)
  # )
  # 
  # ground_truth_dgp <- data.frame("X"=X, "phi"=phi)
  # 
  # g <- ggplot(partial_identification_ATE, aes(x = x, group = curve_index)) +
  #   
  #   # Add shaded regions between lower and upper bounds
  #   geom_ribbon(aes(ymin = lower_ATE_bound, ymax = upper_ATE_bound, fill = curve_index), alpha = 0.5) +
  #   
  #   # Add titles and labels
  #   labs(#title = "ATE bounds corresponding to all feasible S",
  #     x = expression(paste('Exposure ', italic(X))),
  #     y = 'Average Treatment Effect') +
  #   theme_bw() + 
  #   theme(axis.title = element_text(size = 20),
  #         legend.position = "none",
  #         axis.text = element_text(size = 16)) +
  #   geom_line(data = ground_truth_dgp, aes(x=X, y=phi), 
  #             group="Ground truth", linewidth = 1)
  # ggsave(paste0('R = ', tmp$R, ' SNR_y = ', tmp$snr_y, '.png'), width = 8)
  
}
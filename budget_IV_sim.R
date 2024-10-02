# Load libraries
library(data.table)
library(MASS)
library(Matrix)

# Set seed
set.seed(42)

# Hyperparameters
d_z <- 7L
n <- 2000L

# Simulate normal m's
m <- rnorm(d_z, mean = 1)

# Simulate Lambda matrix
Lambda <- matrix(nrow = d_z, ncol = d_z)
for (i in 1:d_z) {
  for (j in 1:d_z) {
    if (i == j & i > 3) {
      Lambda[i, j] <- rbinom(1, 1, prob = 0.9)
    } else if (i != j & i > 3 & j > 3) {
      Lambda[i, j] <- rbinom(1, 1, prob = 0.3)
    } else {
      Lambda[i, j] <- 0L
    }
  }
}

# Z's form a Markov chain
Z <- matrix(nrow = n, ncol = d_z)
Z[, 1] <- rbinom(n, 1, prob = 0.05)
for (j in 2:d_z) {
  Z[, j] <- rbinom(n, 1, prob = ifelse(Z[, (j - 1)] == 1, 0.9, 0.05))
}

# Induce correlation with eps_x, eps_y
x_wts <- rnorm(d_z)
y_wts <- rnorm(d_z)
eps_x <- drop(Z %*% x_wts)
eps_y <- drop(Z %*% y_wts)

# Now we have the structural equation for x
x <- drop(plogis(eps_x - Z %*% m))

# Varying simulation settings
phi <- sqrt(1 + (x - 0.5)^2)
var_y <- 10
gamma <- sapply(1:n, function(i) {
  t(Z[i, ]) %*% Lambda %*% Z[i, ]
})
Z_tilde <- Z[, 4:d_z]          
numer <- norm(cov(gamma, Z_tilde), type = '2') 
denom <- norm(cov(eps_y, Z_tilde), type = '2') 
r <- numer / denom
sim_grd <- as.data.table(
  expand.grid(snr_y = c(1/3, 1, 3), 
              R = c(1/2, 1, 2))
)
sim_grd[, s := (1 / var(phi)) * (var_y / (1 + 1 / snr_y))]
sim_grd[, idx := .I]

# u, lambda_a2, lambda_a3 depend on the simulation setting
for (i in 1:nrow(sim_grd)) {
  tmp <- sim_grd[i, ]
  u <- r / tmp$R * gamma + eps_y
  lambda_a2 <- ((tmp$s * cov(phi, u)) / var(u)) * 
    (-1 + sqrt(1 + (var(u) / (tmp$s * cov(phi, u)^2)) * (var_y / (1 + tmp$snr_y))))
  lambda_a3 <- tmp$R / r * lambda_a2
}





# # Simulate the correlation matrix for Gaussian copula
# d_rho <- d_z + 2L
# R_ww <- diag(nrow = d_rho)
# for (i in 1:(d_rho - 1)) {
#   for (j in (i + 1):d_rho) {
#     R_ww[i, j] <- R_ww[j, i] <- runif(1, min = -1, max = 1)
#   }
# }
# for (k in 1:3) {
#   R_ww[k, d_rho] <- R_ww[d_rho, k] <- 0
# }
# 
# # Positive definite check
# e <- eigen(R_ww)
# pd <- all(e$values >= 0)
# if (!isTRUE(pd)) {
#   while (!isTRUE(pd)) {
#     R_ww <- nearPD(R_ww, corr = TRUE)$mat
#     for (k in 1:3) {
#       R_ww[k, d_rho] <- R_ww[d_rho, k] <- 0
#     }
#     e <- eigen(R_ww)
#     pd <- all(e$values >= 0)
#   }
# }
# 
# # Gaussian copula sim
# mvn_samples <- mvrnorm(n, mu = c(rep(-3, d_z), rep(0, 2)), Sigma = R_ww)
# pi_z <- plogis(mvn_samples[, 1:d_z])
# Z <- matrix(rbinom(n * d_z, size = 1, prob = as.numeric(pi_z)), ncol = d_z)
# eps_x <- mvn_samples[, d_z + 1L]
# eps_y <- mvn_samples[, d_z + 2L]

# W <- pnorm(mvn_samples)
# 
# # Simulate p's from a funky beta
# p <- rbeta(d_z, 1/29, 1/29)
# 
# # Simulate Z's
# Z <- sapply(1:d_z, function(j) as.numeric(W[, j] > p[j]))

# Goal is: MAFs on the order of 5%, potentially strong covariance structure

# New idea: why not just do a multivariate logistic normal distribution?
# Then we'd have some covariance structure over the logits...

# Exogenous epsilons
# eps_x <- qnorm(W[, d_z + 1L])
# eps_y <- qnorm(W[, d_z + 2L])


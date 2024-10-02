library(randcorr)
library(MASS)
library(data.table)

#set.seed(130)

setwd('C:/Users/k23067841/Downloads/BudgetIV_experiments')

generate_parameters_experiment_1 <- function(theta_true, g_true){
  
  rho <- randcorr(4) 
  
  while(sign(rho[4,1]) != sign(g_true[1]) | sign(rho[4,2]) != sign(g_true[2])){
    
    rho <- randcorr(4)
    
  }
  
  eta_e_x <- rexp(1)
  eta_e_y <- rexp(1)
  
  eta_Z_1 <- (g_true[1])/(rho[4,1] * eta_e_y)
  eta_Z_2 <- (g_true[2])/(rho[4,2] * eta_e_y)
  
  std_devs <- c(eta_Z_1, eta_Z_2, eta_e_x, eta_e_y)
  
  sigma <- rho * (std_devs %*% t(std_devs))
  
  #print(rho)
  #print(std_devs)
  #print(sigma)
  
  #evalues <- eigen(sigma)$values
  
  #print(evalues)
  
  #covariance_mat <- 
    
  return(sigma)
  
}

generate_parameters_experiment_2 <- function(theta_true, g_true, beta_x_true){
  
  rho <- matrix(0, 4, 4)
  
  evals_rho <- eigen(rho)$values
  
  num <- 0
  
  while( !all(evals_rho > 0) ){
    
    eta_e_x <- rexp(1,rate = 1)
    eta_e_y <- rexp(1,rate = 1)
    
    eta_Z_1 <- rexp(1,rate = 1)
    eta_Z_2 <- rexp(1,rate = 1)
    
    rho_Z1_ex <- beta_x_true[1]/(eta_Z_1 * eta_e_x)
    rho_Z2_ex <- beta_x_true[2]/(eta_Z_2 * eta_e_x)
    
    rho_Z1_ey <- g_true[1]/(eta_Z_1 * eta_e_y)
    rho_Z2_ey <- g_true[2]/(eta_Z_2 * eta_e_y)
    
    rho_ex_ey <- runif(1, min = -1, max = 1)
    rho_Z1_Z2 <- runif(1, min = -1, max = 1)
    
    rho <- matrix(c(1, rho_Z1_Z2, rho_Z1_ex, rho_Z1_ey, 
                    rho_Z1_Z2, 1, rho_Z2_ex, rho_Z2_ey, 
                    rho_Z1_ex, rho_Z2_ex, 1, rho_ex_ey,
                    rho_Z1_ey, rho_Z2_ey, rho_ex_ey, 1
                    ), nrow = 4, byrow = TRUE)
    
    print(rho)
    
    evals_rho <- eigen(rho)$values
    
    rho_ex_ey <- runif(1, min = -1, max = 1)
    rho_Z1_Z2 <- runif(1, min = -1, max = 1)
    
    rho <- matrix(c(1, rho_Z1_Z2, rho_Z1_ex, rho_Z1_ey, 
                    rho_Z1_Z2, 1, rho_Z2_ex, rho_Z2_ey, 
                    rho_Z1_ex, rho_Z2_ex, 1, rho_ex_ey,
                    rho_Z1_ey, rho_Z2_ey, rho_ex_ey, 1
    ), nrow = 4, byrow = TRUE)
    
    evals_rho <- eigen(rho)$values
    
    num <- num + 1
    
    print(paste0("Sad:( #", num))
    print(paste0(evals_rho))
    
  }
  
  std_devs <- c(eta_Z_1, eta_Z_2, eta_e_x, eta_e_y)
  
  sigma <- rho * (std_devs %*% t(std_devs))
  
  return(sigma)
  
  # M_x1 <- beta_x_true[1]/rho[3,1]
  # M_x2 <- beta_x_true[2]/rho[3,2]
  # 
  # M_y1 <- g_true[1]/rho[4,1]
  # M_y2 <- g_true[2]/rho[4,2]
  
  
  
}

generate_dataset_linear <- function(sigma, theta_true, N){
  
  exogenous_dataset <- mvrnorm(N, mu=c(0,0,0,0), Sigma = sigma)
  
  Z_1 <- exogenous_dataset[, 1]
  Z_2 <- exogenous_dataset[, 2]
  
  e_x <- exogenous_dataset[, 3]
  e_y <- exogenous_dataset[, 4]
  
  X <- e_x
  Y <- theta_true * X + e_y
  
  #print(X)
  
  observable_dataset <- data.table(Z_1, Z_2, X, Y)
  
  return(observable_dataset)
  
}


generate_parameters_hyperbola <- function(){
  
  
  
  
}


# a <- c(1,0,1,0)
# b <- c(0,1,1,0)
# c <- c(1,0,0,1)
# d <- c(0,1,0,1)
# 
# xyz <- rbind(a,b,c,d)
# 
# L <- c(1,1,1,1)
# 
# qr.solve(xyz, L)

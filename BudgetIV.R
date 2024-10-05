#
# Budget constraint solver
#
# Use lpSolve to get feasible region from 
#
#   (i) Cross-covariance estimates, taken as 'oracles'
#
#   (ii) Background budget constraints, parameterised by (tau_i, m_i)_{i=1}^K, also written (tau_vec, m_vec)

library(MASS)
library(lpSolve)

# install.packages("Rglpk")
library(Rglpk)

set.seed(123)

# install.packages("arrangements")
library(arrangements)

BudgetIV <- function(
    beta_y, # Cross covariance vector Cov(Y, Z_vec), or point estimate thereof. A row vector or d_Z * 1 matrix.
    beta_phi, # Cross covariance vector Cov(Phi(X), Z_vec)
    phi_basis, # The actual functions Phi(x), input as a list of expressions
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set
    m_vec, # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i 
           # where II is the indicator function 
    search_dom_ATE, # Domain of X to optimise ATE within (must be specified)
    search_res_ATE=100, # How many values of X to optimite ATE in for each component of X (set to 100 in every component)
    tol=1e-10, # Tolerance for "being along a basis vector". Default set to 1e-10
    # but choice of units for (X,Y,Z) affects the decision rule. 
    dummy_infinity=1e10, # Dummy value for tau_{K+1} := infinity. Adjust if necessary
    delta_beta_y=NA # Half-width 
) {
  
  # print(length(beta_y))
  # print(ncol(beta_phi))
  
  d_Z <- ncol(beta_y)
  
  if (all(is.na(delta_beta_y))){delta_beta_y <- numeric(d_Z)}
  
  # if (search_res_ATE)
  
  # print(beta_phi)
    
  # Warnings
  if (is.unsorted(tau_vec)) {
    stop('Please input tau constraints in increasing order.')
  }
  else if (any(duplicated(tau_vec))){ 
    stop('The same tau constraint cannot be specified twice.')
  }
  else if (is.unsorted(m_vec) || any(duplicated(m_vec))) {
    stop('The vector m must be strictly increasing, please see the definition of boldface m in the manuscript.')
  }
  else if (length(beta_y) != ncol(beta_phi)) {
    stop('Cov(Y, Z) and Cov(Phi(X), Z) must be vectors of the same length for scalar Phi(X). Please call "BudgetIV" for treatment of vector Phi(X).')
  }
  else if (length(beta_y) < nrow(beta_phi)){
    stop('BudgetIV only supports partial identification in the regime d_{Phi} <= d_{Z}')
  }
  else if (length(beta_y) != length(delta_beta_y)){
    stop('If specifying half-width errors delta_beta_y, there must be as many errors as components of beta_y. Run get_covariance with confidence_threshold set to sum numeric in (0,1) to compute these')
  }
  else if (!all(delta_beta_y >= 0)){
    stop("Please ensure hal-width errors are greater than or equal to zero")
  }
  else if (length(phi_basis) != nrow(beta_phi)){
    stop("Please ensure the number of basis features phi_i equals the number of columns (beta_phi)_i = cov(phi_i, Z)")
    }
  
  
  # Remove Z_i for which e_i \in B^{\perp}
  reduced_problem <- reduce_dZ(beta_y,beta_phi,tau_vec,m_vec,tol)
  
  proven_infeasible <- reduced_problem$proven_infeasible
  partially_identifiable <- reduced_problem$partially_identifiable
  
  if(proven_infeasible){
    
    print("Proven infeasible: too many instruments uncorrelated with exposure and strongly correlated with the outcome.")
    return(list("proven_infeasible"=TRUE, "partially_identifiable"=TRUE))
    
  }
  
  else if(!partially_identifiable){
    
    print("Result is vacuous: all theta are feasible because of too many uncorrelated instruments.")
    return(list("proven_infeasible"=FALSE, "partially_identifiable"=FALSE))
    
  }
  
  beta_y <- reduced_problem$beta_y_red
  beta_phi <- reduced_problem$beta_phi_red
  tau_vec <- reduced_problem$tau_red
  m_vec <- reduced_problem$m_red
  
  d_Z <- ncol(beta_y)
  d_Phi <- ncol(beta_phi)
  
  # If there are as many constraints as instruments, don't include tau_{K+1} = infinity
  if(m_vec[length(m_vec)] == d_Z){
    
    # (tau_1, ..., tau_K)
    taus <- tau_vec
    
    # Vector of differences (m_1, m_2 - m_1, ..., m_K - m_{K-1})
    m_deltas <- c(m_vec[1], diff(m_vec))
    
    }
  
  # Otherwise, include tau_{K+1} = infinity (problem is "under-constrained")
  else{
    
    # (tau_1, ..., tau_K, tau_{K+1})
    taus <- c(tau_vec, dummy_infinity)
    
    # Vector of differences (m_1, m_2 - m_1, ..., m_K - m_{K-1}, d_Z - m_K)
    m_deltas <- c(m_vec[1], diff(m_vec), d_Z - m_vec[length(m_vec)])
    
    }
  
  # Iterate through the values of the one-hot encoding S. The iterator maps tau_i to i, so we have to reverse this map. 
  S_perm_iter <- ipermutations(taus, freq=m_deltas) 
  
  curr_S_perm <- S_perm_iter$getnext()
  
  # For every possible perm for every unique way to satisfy the budget constraints. 
  # There are d_Z! /( m_1!(m_2 - m_1)!...(m_{K_red}-m_{K_red - 1})! ) perms. 
  while (!is.null(curr_S_perm)) {
    
    #lower_bounds <- -taus[curr_S_perm] - delta_beta_y
    # lower_bounds <- -taus[curr_S_perm] - delta_beta_y
    bounds <- taus[curr_S_perm] + delta_beta_y
    
    print(bounds)
    
    # print(upper_bounds)
    # 
    # print(lower_bounds)
    
    # lpsolve for testing intersection:
    
    # print("beta_phi")
    # print(beta_phi)
    
    f.con <- rbind(t(beta_phi), - t(beta_phi))  # Combine the upper and lower bound constraints
    f.dir <- c(rep("<=", d_Z), rep("<=", d_Z))  # Directions of the inequalities
    f.rhs <- c(beta_y + bounds, bounds - beta_y)  # Right-hand side for the inequalities
    
    f.obj <- rep(0, d_Phi)  # Objective function for theta (length m, the dimension of theta)
    
    # print("f.con")
    # print(f.con)
    # 
    # print("f.dir")
    # print(f.dir)
    # 
    # print("f.rhs")
    # print(f.rhs)
    
    solution <- lp("min", f.obj, f.con, f.dir, f.rhs, scale=TRUE)
    
    print(solution)
    
    # solution <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = TRUE)
    
    # if (solution$status == 0){
    #   
    #   # print(beta_phi)
    #   # print(beta_y)
    #   # 
    #   # # print(bounds)
    #   # print(curr_S_perm)
    #   print("Soluable")
    #   
    # }
    # 
    # else {
    #   
    #   print("Insoluble")
    #   
    # }
    
    curr_S_perm <- S_perm_iter$getnext()
    
  }
  
}

#
# Remove Z_i from cov(Z_i, g_y) calculations for all i(B_i \in B^{\perp})
# 
# 

reduce_dZ <- function(
    beta_y, # Cross covariance Cov(Y, Z_vec)
    beta_phi, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec, # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
           # where II is the indicator function 
    tol=1e-10 # Tolerance for "being along a basis vector". Default set to 1e-10, 
    # but choice of units for (X,Y,Z) affects the decision rule. 
) {
  
  changed_flag <- FALSE
  
  d_Z <- ncol(beta_y)
  K <- length(tau_vec)
  
  m_new <- m_vec
  to_remove_instruments <- rep(0, d_Z)
  to_remove_constraints <- rep(0, K)
  
  for(instrument in 1:d_Z){
    
    instrument_basis_vec <- matrix(0, ncol=1, nrow=d_Z)
    instrument_basis_vec[instrument] <- 1
    
    # print(beta_phi)
    # print(instrument_basis_vec)
    
    instrument_projection <- beta_phi %*% instrument_basis_vec
    
    if(isTRUE(all(instrument_projection < tol))){
      
      to_remove_instruments[instrument] <- 1
      
      for(k in 1:K){
        
        if(beta_y[, instrument] <= tau_vec[k]){
          m_new[k] <- m_new[k] - 1}}
      
      }
      
    collapse_pos <- 1
    
    
    for(k in 2:K){
      
      if(m_new[k] <= m_new[collapse_pos]){ to_remove_constraints[k] <- 1 }
      
      else{collapse_pos <- k}
        
      }
      
    }
  
  # print(beta_y)
  # print(beta_phi)
  
  m_red <- m_new[to_remove_constraints == 0]
  tau_red <- tau_vec[to_remove_constraints == 0]
  beta_y_red <- beta_y[, to_remove_instruments == 0 , drop=FALSE]
  beta_phi_red <- beta_phi[, to_remove_instruments == 0, drop=FALSE]
  
  # print(beta_y_red)
  # print(beta_phi_red)
  
  partially_identifiable <- isTRUE(m_red[1] > 0)
  proven_infeasible <- isTRUE(m_red[length(m_red)] > length(beta_phi_red))
  
  return(list("m_red" = m_red, "tau_red" = tau_red, "beta_y_red" = beta_y_red, "beta_phi_red" = beta_phi_red, 
              "partially_identifiable" = partially_identifiable, "proven_infeasible" = proven_infeasible))
}


# beta_phi_true <- matrix(c(1.2, -0.4, 1, 0.6, 0.3, -0.9, -0.2, 0.45), ncol=4, nrow=2)
# beta_y_true <- matrix(c(0, -3.9, 3, -2), ncol=4, nrow=1)
# delta_beta_y <- c(0.1,0.2, 0.2, 0.15)

beta_phi_true <- matrix(c(1.2, -0.4, 1, 0.6, 0.3, -0.9, -0.2, 0.45), ncol=4, nrow=2)
beta_y_true <- matrix(c(-0.3, 0.5, 3.4, 3.25), ncol=4, nrow=1)
delta_beta_y <- c(0.1,0.2, 0.2, 0.15)

tau_vec <- c(0.2,3)
m_vec <- c(2,4)
dom_ATE <- matrix(c(-4, 4), ncol=1, nrow=2)
search_res_ATE <- 500

phi_basis <- expression(x, x^2)

# print(length(phi_basis))

BudgetIV(beta_y=beta_y_true, beta_phi=beta_phi_true, phi_basis=phi_basis, tau_vec=tau_vec, m_vec=m_vec, search_dom_ATE=dom_ATE, delta_beta_y=delta_beta_y)


# deltas <- c(1,2)
# D <- 4
# instruments <- 1:D
# my_comb <- combinations(v=instruments, k = sum(deltas), freq = rep(1, length(instruments)), replace=FALSE)
# 
# print(my_comb)




# The basic idea 
# 
#   Generate a tuple of length d_Z, labelled with "bin numbers" 1 through K+1 
#   
#   Each number j can show up delta_j times, and K+1 shows up d_Z - \sum_{j=1}^K delta_j times. 
#   
#   These are permutations, used to define "multisets"? 
#   
#   Parameters: 
#               
#               v = 1:d_Z
#               
#               freq = c(detla_1, delta_2, ..., delta_{K}, delta_{K+1})
#               
#               
#   For each tuple, do a linear constraint-satisfaction test (should be O(d_Z)-complexity)
#   
#   Then, if "linear" or monotonic in each argument, optimise for theta.
#   
#   Otherwise, for each x \in dom(X), optimise ATE_{theta} (x, x_0) := theta \cdot (Phi(x) - Phi(x_0))
#   
#   
#   
#   
#   So what about the simulation study itself (i.e., the choice of Phi(x))?
#   
#   Need a simple example with nonlinear terms (want to do a grid search over x): 
#       
#       Take the hyperbola, expand in Legendre polynomials or some relevant bases(?). Let dim(Phi) = 5 or so. Bound dom(X) with compact support 
#       
#       For g_y, let's do a (nearly) linear model with interaction terms (e.g., Z_1 shows up conditional on Z_2). Let d_Z = 7
#       
#       For f_X, also do a linear model with interaction terms. 
#       
#       Let confounders with X, Y have some (weak) correlation with all Z_i, and maybe some strong correlation with one of them. 
#       
#       
#       Fix Cov(Phi(X), Z) and Cov(Y, Z) (this also fixes Cov(g_y(epsilon_y, Z), Z). Need some paper to figure this out.
#       
#       
#       Finally, randomise some of this procedure for 
# 

# deltas <- c(1,2)
# tau_vec <- c(1, 5)
# dummy_infinity <- 1e10
# taus <- c(tau_vec, dummy_infinity)
# D <- 3
# deltas <- c(deltas, D - sum(deltas))
# group_labels <- rep(taus, times = deltas)
# print(group_labels)

#permutations(group_labels, k = D, replace=FALSE)


# d_Z <- 4
# 
# tau_vec <- c(1,5)
# 
# dummy_infinity <- 1e10
# 
# tau_vec_infty <- c(tau_vec, dummy_infinity)
# 
# m_vec <- c(1,3)
# 
# deltas <- c(m_vec[1], diff(m_vec), d_Z - m_vec[length(m_vec)])
# 
# print(tau_vec_infty)
# print(deltas)
# 
# permutations(tau_vec_infty, freq=deltas)
# 



#
# Budget constraint solver
#
# Use lpSolve to get feasible region from 
#
#   (i) Cross-covariance estimates, taken as 'oracles'
#
#   (ii) Background budget constraints, parameterised by (tau_i, m_i)_{i=1}^K, also written (tau_vec, m_vec)

# install.packages("Rglpk")
# install.packages("arrangements")

set.seed(123)

library(MASS)
library(data.table)
library(ggplot2)

library(arrangements)

# library(lpSolve)
library(Rglpk)

BudgetIV <- function(
    beta_y, # Cross covariance vector Cov(Y, Z_vec), or point estimate thereof. A row vector or d_Z * 1 matrix
    beta_phi, # Cross covariance vector Cov(Phi(X), Z_vec)
    phi_basis, # The actual functions Phi(x), input as a list of expressions
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set
    m_vec, # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i 
           # where II is the indicator function 
    search_dom_ATE, # Domain of X to optimise ATE within (must be specified)
    search_res_ATE=101, # How many values of X to optimite ATE in for each component of X (set to 100 in every component)
    d_X=1, # Dimension of the exposure variable X
    X_baseline=0, # The baseline for the dose response curve y = ATE(x, x_baseline)
    tol=1e-10, # Tolerance for "being along a basis vector". Default set to 1e-10
    # but choice of units for (X,Y,Z) affects the decision rule. 
    dummy_infinity=1e10, # Dummy value for tau_{K+1} := infinity. Adjust if necessary
    delta_beta_y=NA # Half-width 
) {
  
  # print(length(beta_y))
  # print(ncol(beta_phi))
  
  d_Z <- ncol(beta_y)
  
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
  d_Phi <- nrow(beta_phi)
  
  if (all(is.na(delta_beta_y))){delta_beta_y <- numeric(d_Z)}
  
  if (length(search_res_ATE)==1){search_res_ATE <- rep(search_res_ATE, d_X)}
  
  # print(search_res_ATE)
  
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
  
  # List to contain all dose reponse ATE bounds and the corresponding decision variable S
  partial_identification_ATE <- data.table(
    "curve_index" = numeric(),
    "S" = matrix(nrow=0,ncol=d_Z),
    "x" = numeric(),
    "lower_ATE_bound" = numeric(),
    "upper_ATE_bound" = numeric(),
    "S_string" = character()
  )
  
  curve_index <- 1
  
  # Iterate through the values of the one-hot encoding S. The iterator maps tau_i to i, so we have to reverse this map. 
  S_perm_iter <- ipermutations(taus, freq=m_deltas) 
  
  curr_S_perm <- S_perm_iter$getnext()
  
  # For every possible perm for every unique way to satisfy the budget constraints. 
  # There are d_Z! /( m_1!(m_2 - m_1)!...(m_{K_red}-m_{K_red - 1})! ) perms. 
  while (!is.null(curr_S_perm)) {
    
    # print(curr_S_perm)
    
    #lower_bounds <- -taus[curr_S_perm] - delta_beta_y
    # lower_bounds <- -taus[curr_S_perm] - delta_beta_y
    bounds <- taus[curr_S_perm] + delta_beta_y
    
    # print(bounds)
    
    # print(taus[curr_S_perm])
    # print(delta_beta_y)
    # 
    # print(beta_y)
    # print(beta_phi)
    
    # print(upper_bounds)
    # 
    # print(lower_bounds)
    
    # lpsolve for testing intersection:
    
    # print("beta_phi")
    # print(beta_phi)
    
    f.con <- rbind(t(beta_phi), t(beta_phi))  # Combine the upper and lower bound constraints
    f.dir <- c(rep("<=", d_Z), rep(">=", d_Z))  # Directions of the inequalities
    f.rhs <- c(beta_y + bounds, beta_y - bounds)  # Right-hand side for the inequalities
    
    f.obj <- rep(0, d_Phi)  # Objective function for theta (length d_Phi, the dimension of theta)
    
    search_bounds <- list(lower = list(ind = c(1L, 2L), val = c(-Inf, -Inf)),
                          upper = list(ind = c(1L, 2L), val = c(Inf, Inf)))
    
    # print("f.con")
    # print(f.con)
    # 
    # print("f.dir")
    # print(f.dir)
    # 
    # print("f.rhs")
    # print(f.rhs)
    # 
    # constraint_satisfaction <- lp("min", f.obj, f.con, f.dir, f.rhs, scale=0)
    constraint_satisfaction <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = TRUE, bounds = search_bounds)
    
    # print(solution)
    
    # print(solution$optimum)
    # print(solution$status)
    
    
    
    
    # 
    # Quick and dirty solution for d_X = 1
    # 
    
    
    
    
    if (constraint_satisfaction$status == 0){
      
      # dose_response_bounds_lower <- rep(0, search_res_ATE)
      # dose_response_bounds_upper <- rep(0, search_res_ATE)
      
      print(constraint_satisfaction$solution)
      
      x_min <- dom_ATE[,1]
      x_max <- dom_ATE[,2]
      
      for (coord_index in 1:search_res_ATE) {
        
        curr_coord <- x_min + (coord_index - 1) * (x_max - x_min) / (search_res_ATE-1)
        
        # print(curr_coord)
        
        phi_curr_coord <- lapply(phi_basis, function(e) eval(e, list(x=curr_coord)))
        
        # print(c(phi_curr_coord[[1]], phi_curr_coord[[2]]))
        
        f.obj <- c(phi_curr_coord[[1]], phi_curr_coord[[2]])
        
        ATE_min_curr_coord <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = FALSE, bounds = search_bounds)
        ATE_max_curr_coord <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = TRUE, bounds = search_bounds)
        
        # dose_response_bounds_lower[coord_index] <- ATE_min_curr_coord$optimum
        # dose_response_bounds_upper[coord_index] <- ATE_max_curr_coord$optimum
       
        new_bound <- data.table(
          "curve_index" = curve_index,
          "S" = matrix(curr_S_perm, nrow = 1),
          "x" = curr_coord,
          "lower_ATE_bound" = ATE_min_curr_coord$optimum,
          "upper_ATE_bound" = ATE_max_curr_coord$optimum,
          "S_string" = paste(curr_S_perm, collapse = ", ")
        )
        
        print(new_bound)
        
        partial_identification_ATE <- rbind(partial_identification_ATE, new_bound)
         
      }
      
      curve_index <- curve_index + 1
      
    }
    
    # 
    # Grid over d_X > 1 for general BudgetIV solution
    # 
    # 
    
    
    # if (solution$status == 0){
    #   
    #   search_function <- lapply(search_res_ATE, function(n) seq_len(n))
    #   search_grid <- expand.grid(search_function)
    #   
    #   # print(search_grid)
    #   
    #   print(search_grid[[1]][5])
    #   
    #   # print(search_function)
    #   # (dom_ATE[i, 1], dom_ATE[i, 2], length.out = search_res_ATE[i])
    #   # print(search_grid[])
    #   
    #   for (search_index in 1:nrow(search_grid)) {
    #     
    #     # print(search_index)
    #     
    #     # curr_coord <-  (search_index-1)
    #     
    #     # print(curr_coord)
    # 
    #   }
    # 
    # }
    
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
  
  return(partial_identification_ATE)
  
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

beta_phi_true <- matrix(c(1.2, -0.4, 1, 0.6, 3, 9, -2, 4.5), ncol=4, nrow=2)
beta_y_true <- matrix(c(0.5, 1, 3.5, 2.9), ncol=4, nrow=1)
# beta_y_true <- 10*matrix(c(300, 50, 34, 325), ncol=4, nrow=1)
delta_beta_y <- c(0,0,0,0)

tau_vec <- c(0.2,3)
m_vec <- c(2,4)
dom_ATE <- matrix(c(-4, 4), ncol=2, nrow=1)
search_res_ATE <- 501

phi_basis <- expression((x^2 + x + 4), 0.7*(x^2 - 3*x))

# print(length(phi_basis))

partial_identification_ATE <- BudgetIV(beta_y=beta_y_true, beta_phi=beta_phi_true, phi_basis=phi_basis, tau_vec=tau_vec, m_vec=m_vec, search_dom_ATE=dom_ATE, delta_beta_y=delta_beta_y)

print(partial_identification_ATE)


my_plot <- ggplot(partial_identification_ATE, aes(x = x, group = curve_index, fill = S_string)) +
  
  # Add shaded regions between lower and upper bounds
  geom_ribbon(aes(ymin = lower_ATE_bound, ymax = upper_ATE_bound), alpha = 0.3) +
  
  # Add lines for the upper bounds
  geom_line(aes(y = lower_ATE_bound), linewidth = 0.5) +
  
  # Add lines for the lower bounds
  geom_line(aes(y = upper_ATE_bound), linewidth = 0.5) +
  
  # Add titles and labels
  labs(title = "ATE bounds corresponding to all feasible S",
       x = "x",
       y = "ATE(x, x_0 := 0)",
       fill = "Value of S") +
  theme_minimal()

print(my_plot)

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



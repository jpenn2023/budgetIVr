#
# Budget constraint solver
#
# Use lpSolve to get feasible region from 
#
#   (i) Cross-covariance estimates, taken as 'oracles'
#
#   (ii) Background budget constraints, parameterised by (tau_i, m_i)_{i=1}^K, also written (tau_vec, m_vec)

BudgetIV <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec)
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
          # where II is the indicator function 
) {
  
  
  
}

#
# Remove Z_i from cov(Z_i, g_y) calculations for all i(e_i \notin Span(B))
# 
# 

reduce_dZ <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
) {
  
  # Assuming the budget constraints are feasible, is the feasible set of theta bounded?
  bounded <- TRUE
  
  
}




#
#
#
#
# Special case of d_{Phi} = 1 yields very efficient, polytime solution.
#
#
#
#



BudgetIV_scalar_exposure_feature <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
) {
  
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
  else if (length(A) != length(B)) {
    stop('Cov(Y, Z) and Cov(Phi(X), Z) must be vectors of the same length for scalar Phi(X). Please call "BudgetIV" for treatment of vector Phi(X).')
  }
  
  
  # Run polytime tests and get rid of 'sticky' directions in Z.
  polytime_sln <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
  
  
  m_vec_red <- polytime_sln$m_new
  J_non_sticky <- polytime_sln$J
  identifiable <- polytime_sln$identifiable
  proven_infeasible <- !polytime_sln$feasible
  
  if (proven_infeasible){
    print("Infeasible!")
    return("feasible" = FALSE)
    }
  
  else if(!identifiable){
    print("Unidentifiable")
    return("identifiable" = FALSE)
    }
  
  # Position of tau_vec such that tau_vec[i] represents intervals of feasible theta rather than points
  intervals_start_flag <- 1
  
  
  #
  # Find all theta where the number of bucket constraints satisfied changes (either singularly at that point or changes as theta increases past that point) 
  #
  
  theta_points <- NA
  
  if (tau_vec[1] == 0){ 
    
    intervals_start_flag <- 2
    
    theta_points <- numeric(length(J_non_sticky))
    
    # theta_points(j) are defined when tau_1 = 0 
    for(j in J_non_sticky) {
      theta_points[j] <- A[j] / B[j]
    }}
  
  # Bounds in theta where a specific A[j] - B[j] * theta = +- tau_vec[k]
  
  tau_intervals_lower <- matrix(nrow = (length(tau_vec)+1-intervals_start_flag), ncol = length(J_non_sticky))
  
  tau_intervals_upper <- matrix(nrow = (length(tau_vec)+1-intervals_start_flag),  ncol = length(J_non_sticky))
  
  for (k in intervals_start_flag:length(tau_vec)){
    for (j_prime in 1:length(J_non_sticky)){
      j <- J_non_sticky[j_prime]
      
      tau_intervals_lower[k, j_prime] <- A[j] / B[j] - abs(tau_vec[k]) / B[j]
      
      tau_intervals_upper[k, j_prime] <- A[j] / B[j] + abs(tau_vec[k]) / B[j]
      
    }}
  
  #
  # Find the full feasible set of theta, including points (measure zero subsets) and intervals along RR. 
  #
  
  if (!is.na(theta_points)){
  
    # Points which can be tested one at a time:
    feasible_points <- rep(NA, length(theta_points))
    
    for(p in 1:length(theta_points)){
      
      if(validPoint_scalar(A, B, J_non_sticky, theta_points[p], m_vec_red, tau_vec)) {feasible_points[p] <- theta_points[p]}
      
    na.omit(feasible_points)
    sort(feasible_points)
    
    }
  }
  
  else {feasible_points <- NA}
  
  # Feasible intervals
  feasible_intervals <- matrix(nrow = 0, ncol = 2)
  
  possible_bounds <- sort(c(c(tau_intervals_lower), c(tau_intervals_upper)))
  
  #print(validPoint_scalar(A, B, J_non_sticky, -0.2, m_vec_red, tau_vec))
  
  # Search through all possible points at which a feasible interval could begin or end:
  in_feasible <- FALSE
  
  for (p in 1:(length(possible_bounds)-1)){
    
    curr_point <- (possible_bounds[p] + possible_bounds[p+1])/2
    
    # Make sure interval bounds are not biased by singularity at theta \in theta_points
    if (curr_point %in% feasible_points){while (curr_point %in% feasible_points){ curr_point <- (curr_point + possible_bounds[p])/2} }

    # Keep track of whether inside or outside feasible region, including last theta at which a feasible interval was entered  
    # Add a new feasible interval once left
    curr_point_feasible <- validPoint_scalar(A, B, J_non_sticky, curr_point, m_vec_red, tau_vec)
    
    if(curr_point_feasible && !in_feasible){
        
        last_feasible_opening = possible_bounds[p]
        in_feasible <- TRUE
        
      }
      
      else if(!curr_point_feasible && in_feasible){
        
        in_feasible <- FALSE
        feasible_intervals <- rbind(feasible_intervals, c(last_feasible_opening,  possible_bounds[p]))
        
      }
      
    }
  
  # Add final interval if it ends at the final theta
  if(in_feasible){ feasible_intervals <- rbind(feasible_intervals, c(last_feasible_opening,  possible_bounds[length(possible_bounds)])) }
  
  return(list("points" = feasible_points, "intervals" = feasible_intervals))
  
}

reduce_dZ_scalar_exposure_feature <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
) {
  
  # Redefine budget constraints once all i \in [d_Z] (B_i = 0) are put into correct groups
  m_new <- m_vec
  
  for(i in 1:length(B)){
    
    if(B[i] == 0){
      for(j in 1:length(tau_vec)){
        
        if(A[i] <= tau_vec[j]){
          m_new[j] <- m_new[j] - 1}}}}
  
  # The set {i \in [d_Z] : B_i =/= 0}, called \notin I in the manuscript.
  J_non_sticky <- which(B!=0)
  
  # max_{k \in K} (m_new[k])
  m_hardest <- max(m_new)
  
  # Return if unidentifiable (feasible set = RR)
  # follows from theorem labelled "Unidentifiability" in the manuscript.
  if(all(m_new <= 0)){return(list("J" = J_non_sticky, "m_new" = m_new, "identifiable" = FALSE, "feasible" = TRUE))}
  
  # Return if proven infeasible (feasible set empty) 
  # follows from theorem labelled "Sufficient condition for infeasibility"
  else if(m_hardest > length(J_non_sticky)){return(list("J" = J_non_sticky, "m_new" = m_new, 
                                                   "identifiable" = NA, "feasible" = FALSE))}
  
  # Otherwise, return reduced-d_Z problem
  else{return(list("J" = J_non_sticky, "m_new" = m_new, "identifiable" = TRUE, "feasible" = TRUE))}
  
}



# Validate a point in RR^{d_Z} is in the feasible set.
# Used for efficient feasible set generation when d_Phi = 1.
#
validPoint_scalar <- function(A, B, J_non_sticky, theta, m_vec_red, tau_vec){
  
  m_to_fill <- m_vec_red
  
  for(j in J_non_sticky){
    
    g_j <- A[j] - B[j] * theta
    
    for(k in 1:length(tau_vec)){
      
      if(abs(g_j) <= tau_vec[k]){ m_to_fill[k] <- m_to_fill[k] - 1 }
      
    }
  }
  
  print(m_to_fill)
  
  return(all(m_to_fill <= 0))
  
}



# A <- c(-4.31, 1.7686, 3.4342, 2.234)
# B <- c(-2.212, 0.9, 1.23343, 11)
# tau_vec <- c(0.345, 4.23)
# m_vec <- c(1, 2)

# A <- c(-4, 1)
# B <- c(-2, 0.9)
# tau_vec <- c(0.3, 10)
# m_vec <- c(1,2)

# A <- c(-1, 1.5)
# B <- c(1, 0.1)
# tau_vec <- c(1, 2)
# m_vec <- c(1, 2)

# 
# full_trial <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
# print(full_trial)

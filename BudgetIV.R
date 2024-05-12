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
# Special case of d_{Phi} = 1 yields very efficient solution.
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
  
  # Run polytime tests and get rid of 'sticky' directions in Z.
  polytime_sln <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
  
  J_non_sticky <- polytime_sln$I
  m_new <- polytime_sln$m_new
  identifiable <- polytime_sln$identifiable
  feasible <- polytime_sln$feasible
  
  print(J_non_sticky)
  print(m_new)
  print(feasible)
  print(identifiable)
  
  return(polytime_sln)
  
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
  if(all(m_new) <= 0){return(list("J" = J_non_sticky, "m_new" = m_new, "indentifiable" = FALSE, "feasible" = TRUE))}
  
  # Return if proven infeasible (feasible set empty) 
  # follows from theorem labelled "Sufficient condition for infeasibility"
  else if(m_hardest > length(J_non_sticky)){return(list("J" = J_non_sticky, "m_new" = m_new, 
                                                   "indentifiable" = NA, "feasible" = FALSE))}
  
  # Otherwise, return reduced-d_Z problem
  else{return(list("J" = J_non_sticky, "m_new" = m_new, "indentifiable" = TRUE, "feasible" = TRUE))}
  
}

A <- c(-4, 1, 3, 2)
B <- c(-2, 0.9, 0, 0)
tau_vec <- c(4, 0.3)
m_vec <- c(1, 1)

#
# Test of "reduce_dZ_scalar_exposure_feature" function, one run for each case
#
#

# Identifiable and feasible

red_trial <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(red_trial)

BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)


# get_p_vec <- function(){
#   return(list("p" = 3, "m" = list(2,1)))
# }
# 
# p_vec <- get_p_vec()
# 
# print(p_vec[2])
# m <- p_vec$m
# 
# print(m)
# 
# for(i in 1:length(m)){
#   print(m[i])
# }
# 
# my_vec <- c(0, -1, -2)
# 
# if(all(my_vec <= 0)){print("gotcha")}

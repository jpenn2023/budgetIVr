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
# Remove axes cov(Z_i, g_y) for which e_i \notin Span(B)
# 
# 

reduce_dZ <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
) {
  
  
  
}

#
# Special case of d_{Phi} = 1 yeilds very efficient solution.
#
#

BudgetIV_scalar_exposure_feature <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
) {
  
  
  
}

reduce_dZ_scalar_exposure_feature <- function(
    A, # Cross covariance Cov(Y, Z_vec)
    B, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    m_vec # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
) {
  
  
  
}
# 
# An approach to estimation with 
# 
# 
# 
# In a sense 'conservative' to consider all d_Z - n_max "valid" genes equally plausible
#   we also keep track of tau_min(n \in n_max) to get a two-factor readout.
# 
# 
# Rough outline: 
# 
# > Input \beta_x, \beta_y 
# 
# \tau_check <- argmin_{\theta \in \R} || \beta_y - \beta_x \theta ||_{\infty}
# 
# > For m \in [d_Z] (increasing){
# 
# >   \tau <- \tau_check  
# 
# >   If n_{max}(\tau_check)
# 
# >   
# 
# > }
# 
# 




budgetMR <- function(A, B){
  
  
  
  
}

find_n_max <- function(A, B, d_Z, n_max){
  
  # Return all open/closed intervals 
  
  
  
}

tau_check <- function(A, B, d_Z, n_max, tau){
  
  
  
}
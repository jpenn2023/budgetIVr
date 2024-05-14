source("BudgetIV.R")

#
# Test of "reduce_dZ_scalar_exposure_feature" function, one run for each case
#
#

# Identifiable and not provably infeasible

A <- c(-4, 1, 3, 2)
B <- c(-2, 0.9, 1, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

red_trial <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(red_trial)

# Unidentifiable

A <- c(-4, 1, 0, 0)
B <- c(-2, 0.9, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

red_trial <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(red_trial)

# Provably infeasible

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

red_trial <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(red_trial)



source("BudgetIV.R")


#
# Test checks and warnings in "BudgetIV_scalar_exposure_feature"
# 
# 
# 
# 

# Okay input

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

warning_test <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(warning_test)


# Duplicate tau 

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(0.3, 0.3)
m_vec <- c(1, 2)

warning_test <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)

# Decreasing tau_vec

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(4, 0.3)
m_vec <- c(1, 2)

warning_test <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)


# Decreasing m_vec

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(2, 1)

warning_test <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)

# Duplicate m

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(2, 2)

warning_test <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)


#
# Personal use
#
# 
# me <- TRUE
# you <- !me
# print(you)
# 
# 
# m_pty <- matrix(nrow = 0, ncol = 2)
# print(m_pty)
# 
# n_pty <- rbind(m_pty, c(-1, 2))
# print(n_pty)



# 
# 
# Test actual solution for d_Phi = 1
# 
# 

source("BudgetIV.R")

A <- c(-4, 1, 3, 2)
B <- c(-2, 0.9, 1, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

full_trial <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(full_trial)

# Unidentifiable

A <- c(-4, 1, 0, 0)
B <- c(-2, 0.9, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

full_trial <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(full_trial)

# Provably infeasible

A <- c(-4, 10, 10, 10)
B <- c(-2, 0, 0, 0)
tau_vec <- c(0.3, 4)
m_vec <- c(1, 2)

full_trial <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
print(full_trial)

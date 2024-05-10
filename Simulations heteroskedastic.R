# Simulated data with preprocessing step for Phi(X)

#preprocess_X <- function(X_dat, model){

#  dim_X <- 

#}

#
# Generate model with randomised parameters satisfying constraints:
#
# X = f_x (Z_1, Z_2, e_x)
# Y = theta_true * X + g_y(Z_1, Z_2, e_m, e_a)
#
# g_y(Z_1, Z_2, e_m, e_a) = e_m * Z_1 + e_a 
#
# e_m ~ U[0,1] (indep to everything else)
#
# Cov(Z_1, g_y) = -2 (high)
# Cov(Z_2, g_y) = 0.1 (low)
# 
#   => Cov(e_a, Z_1) = 0.1
#
# theta_true = 1
#
# Cov(Y, (Z_1, Z_2)) = (-4, 1)
#
#   => Cov(X, (Z_1, Z_2)) = (-2, 0.9)
#

#
# Free parameters to randomise
#
# f_x, which we'll pick randomly out of 3 choices 
#
#   (i) f_x(Z_1, Z_2, e_x) = Z_1 + Z_2 + e_x
#
#   (ii) f_x(Z_1, Z_2, e_x) = min(Z_1, Z_2, e_x)
#
#   (iii) f_x(Z_1, Z_2, e_x) = {0, ; a, }
#
#   Each case e_x is normally distributed.
#
#   Also need to pick Cov(e_x, (Z_1, Z_2)) within this to fit Cov(X, (Z_1, Z_2)) constraints.
#
#   
#
# 


set.seed(987)


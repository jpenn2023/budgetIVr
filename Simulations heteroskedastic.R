# Load libraries, register cores
library(data.table)

library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)


# Set seed (first attempt)
set.seed(987)

random_simulation <- function(A, B, theta_true){}

randomise_test_heteroskedastic <- function(A, B){
  
  confounding_coeffs <- list("Z1Z2" = NA, "Z1ex" = NA, "Z1ea" = NA, 
                             "Z2ex" = NA, "Z2ea" = NA, "exea" = NA)
  
  variances <- list("Z1" = NA, "Z2" = NA, "ex" = NA, "ea" = NA)
  std_devs <- list("Z1" = NA, "Z2" = NA, "ex" = NA, "ea" = NA)
  
  if (B[1] > 0){}
  
}

simulate_data <- function() {}
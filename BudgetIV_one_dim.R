#
#
#
#
# Special case of d_{Phi} = 1 yields very efficient, polytime solution.
#
#
#
#

library(data.table)
library(graphics)
setwd('C:/Users/k23067841/Downloads/BudgetIV')
source("Dataset simulations.R")
setwd('C:/Users/k23067841/Downloads/BudgetIV_experiments')

d <- c(3,2,1)
t <- c(0.3, 1, 3.14)
S <- rep(1:length(d), t)
print(d)
print(S)

#estimate_covariance_matrix <- function(dataset){
  
#  A <- c(cov(dataset$Z_1, curr_data$Y), cov(dataset$Z_2, curr_data$Y))
#  B <- c(cov(dataset$Z_1, curr_data$X), cov(dataset$Z_2, curr_data$X))
  
#}


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
  reduced_problem <- reduce_dZ_scalar_exposure_feature(A, B, tau_vec, m_vec)
  
  
  m_vec_red <- reduced_problem$m_new
  J_non_sticky <- reduced_problem$J
  identifiable <- reduced_problem$identifiable
  proven_infeasible <- !reduced_problem$feasible
  
  if (proven_infeasible){
    print("Infeasible!")
    return("feasible" = FALSE)
  }
  
  else if(!identifiable){
    print("Unidentifiable")
    return("identifiable" = FALSE)
  }
  
  # Treat tau_vec[1] = 0 and tau_vec[1] > 0 seperately, since |A_i - B_i theta| = 0 is satisfied at exactly one point in the reduced problem.
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
  
  #
  # Bounds in theta where a specific A[j] - B[j] * theta = +- tau_vec[k]
  #
  
  tau_intervals_lower <- matrix(nrow = (length(tau_vec)+1-intervals_start_flag), ncol = length(J_non_sticky))
  
  tau_intervals_upper <- matrix(nrow = (length(tau_vec)+1-intervals_start_flag),  ncol = length(J_non_sticky))
  
  for (k in 1:length(tau_vec)+1-intervals_start_flag){
    for (j_prime in 1:length(J_non_sticky)){
      j <- J_non_sticky[j_prime]
      
      tau_intervals_lower[k, j_prime] <- A[j] / B[j] - tau_vec[k+1-intervals_start_flag] / B[j]
      
      tau_intervals_upper[k, j_prime] <- A[j] / B[j] + tau_vec[k+1-intervals_start_flag] / B[j]
      
    }}
  
  #
  # Find the full feasible set of theta, including points and intervals along \mathbb{R}. 
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
  
  #print(m_to_fill)
  
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



#
# L2 ball for benchmarking
#
#
#

l2_ball <- function(A, B, tau_vec, slack){
  
  tau <- slack * norm(tau_vec, type="2")
  
  C_AA <- A %*% A
  C_AB <- A %*% B
  C_BB <- B %*% B
  
  feasible_int_l2 <- list("theta_lo"=NA, "theta_hi"=NA)
  
  if (C_AB^2 >= C_BB*(C_AA - tau^2)){
    
    theta_low <- (C_AB/C_BB)*(1 - sqrt(1 - C_BB * (C_AA - tau^2) / C_AB^2))
    
    #print(theta_low)
    
    theta_hi <- (C_AB/C_BB)*(1 + sqrt(1 - C_BB * (C_AA - tau^2) / C_AB^2))
    
    #print(theta_hi)
    
    return(list("theta_lo"=theta_low, "theta_hi"=theta_hi))
    
  }
  
  return(feasible_int_l2)
  
}

# Sanity check to plot intersection between A - B \theta and \Gamma in the special case of d_Z = 2

plot_intersection_2D <- function(A, B, tau_vec, m_vec, tau_l2){
  
  #plot(c(1, 9), 1:2, type = "n", xlab = "Time", ylab = "Distance")
  
  # 
  if(length(tau_vec) == 2){
    
    tau_1 <- tau_vec[1]
    tau_2 <- tau_vec[2]
    tau_l2 <- tau_l2
    
    #print(tau_l2)
    
    
    x <- c(-tau_1, tau_1, tau_1, tau_2, tau_2, tau_1, tau_1, -tau_1, -tau_1, -tau_2, -tau_2, -tau_1)
    y <- c(tau_2, tau_2, tau_1, tau_1, -tau_1, -tau_1, -tau_2, -tau_2, -tau_1, -tau_1, tau_1, tau_1)
    
    png(paste0("Example g plots/tau_max ", tau_2, ".png"))
    
    plot(x,y,type='n', xlim = c(-1.5*tau_2, 1.5*tau_2), ylim = c(-1.5*tau_2, 1.5*tau_2))
    grid()
    
    symbols(0,0, circles=tau_l2, bg="deepskyblue", fg="deepskyblue", add=TRUE, inches=FALSE)
    polygon(x, y, col = "green", density = 50)
    
    theta <- seq(-1000, 1000, 1)
    
    x_line <- A[1] - B[1]*theta
    y_line <- A[2] - B[2]*theta
    
    intercept <- A[2] - A[1]*B[2]/B[1]
    
    slope <- B[2]/B[1]
    
    abline(intercept, slope)
    
    dev.off()
    
  }
  
}


experiment_1 <- function(){
  
  theta_true <- 1
  g_true <- c(-2, 0.1)
  N <- 1000
  random_exeriments <- 100
  
  tau_resolution <- 1
  fixed_tau <- 0.2
  
  slack_scaling <- 1
  
  tau_vec_max <- random_exeriments/tau_resolution
  
  results <- data.table("sim_idx"=numeric(), "interval"=logical(), "feasible_start"=numeric(), "feasible_end"=numeric(), "tau_lo"=numeric(), "tau_hi"=numeric(), "var_tau"=numeric())
  results_l2 <- data.table("sim_idx"=numeric(), "theta_lo"=numeric(), "theta_hi"=numeric(), 
                           "tau_lo"=numeric(), "tau_hi"=numeric(), "var_tau"=numeric(), "slack"=numeric(), "tau_radius"=numeric())
    
  #plot(1:10, 1:10, type='n', xlim = c(-1.5*tau_vec_max, 1.5*tau_vec_max), ylim = c(-1.5*tau_vec_max, 1.5*tau_vec_max))
  
  
  for(i in 1:random_exeriments){
    
    if(i/tau_resolution < fixed_tau){
      
      tau_vec <- c(i/tau_resolution, fixed_tau)
      m_vec <- c(1,2)
      
    }
    else if(i/tau_resolution == fixed_tau){
      
      tau_vec <- c(fixed_tau)
      m_vec <- c(2)
      
    }
    else{
      
      tau_vec <- c(fixed_tau, i/tau_resolution)
      m_vec <- c(1,2)
      
    }
    
    curr_radius <- slack_scaling * norm(tau_vec, type="2")
    
    sigma_experiment <- generate_parameters_experiment_1(theta_true, g_true)
    
    dataset <- generate_dataset_experiment_1(sigma_experiment, theta_true, N)
    
    A <- c(cov(dataset$Z_1, dataset$Y), cov(dataset$Z_2, dataset$Y))
    B <- c(cov(dataset$Z_1, dataset$X), cov(dataset$Z_2, dataset$X))
    
    feasible_region <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
    
    plot_intersection_2D(A, B, tau_vec, m_vec, curr_radius)
    
    feasible_points <- feasible_region$points
    feasible_intervals <- feasible_region$intervals
  
    curr_feasible_int_l2 <- l2_ball(A, B, tau_vec, slack_scaling)
    
    #print(curr_feasible_int_l2)
    
    if (!is.na(curr_feasible_int_l2$theta_lo)){
      
      #print("hi")
      
      results_l2 <- rbind(results_l2, list("sim_idx"=i, 
                                           "theta_lo"=curr_feasible_int_l2$theta_lo, "theta_hi"=curr_feasible_int_l2$theta_hi,
                                           "tau_lo"= tau_vec[1], "tau_hi"=tau_vec[2], "var_tau"=i/tau_resolution, "slack"=slack_scaling, "tau_radius"=curr_radius))
      
    }
    
    if (!is.na(feasible_points)) {
      for (point_idx in 1:length(feasible_points)){
  
        next_res <- list("sim_idx"=i, "interval"=FALSE, "feasible_start"=feasible_points[point_idx],
                         "feasible_end"=feasible_points[point_idx], "tau_lo"=tau_vec[1], "tau_hi"=tau_vec[2], "var_tau"=i/tau_resolution)
  
        results <- rbind(results, next_res)
      }}
    
    if (length(feasible_intervals) != 0) {
      # print(feasible_intervals)
  
      # print(length(feasible_intervals))
  
      for (interval_idx in 1:(length(feasible_intervals)/2)){
  
        #print(interval_idx)
  
        next_res <- list("sim_idx"=i, "interval"=TRUE, "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2], "tau_lo"=tau_vec[1], "tau_hi"=tau_vec[2], "var_tau"=i/tau_resolution)
  
        # print(next_res)
        results <- rbind(results, next_res)
  
      }}
  
    fwrite(results, './small dphi experiments/results.csv')
    fwrite(results_l2, './small dphi experiments/results_l2.csv')
    
  }
  
  print(feasible_region)

}


experiment_2 <- function(){
  
  theta_true <- 1
  g_true <- c(-2, 0.1)
  beta_x_true <- c(2, -4)
  N <- 10000
  random_exeriments <- 100
  
  tau_resolution <- 10
  fixed_tau <- 0.2
  
  slack_scaling <- 1
  
  tau_vec_max <- random_exeriments/tau_resolution
  
  results <- data.table("sim_idx"=numeric(), "interval"=logical(), "feasible_start"=numeric(), "feasible_end"=numeric(), "tau_lo"=numeric(), "tau_hi"=numeric(), "var_tau"=numeric())
  results_l2 <- data.table("sim_idx"=numeric(), "theta_lo"=numeric(), "theta_hi"=numeric(), 
                           "tau_lo"=numeric(), "tau_hi"=numeric(), "var_tau"=numeric(), "slack"=numeric(), "tau_radius"=numeric())
  
  #plot(1:10, 1:10, type='n', xlim = c(-1.5*tau_vec_max, 1.5*tau_vec_max), ylim = c(-1.5*tau_vec_max, 1.5*tau_vec_max))
  
  
  for(i in 1:random_exeriments){
    
    if(i/tau_resolution < fixed_tau){
      
      tau_vec <- c(i/tau_resolution, fixed_tau)
      m_vec <- c(1,2)
      
    }
    else if(i/tau_resolution == fixed_tau){
      
      tau_vec <- c(fixed_tau)
      m_vec <- c(2)
      
    }
    else{
      
      tau_vec <- c(fixed_tau, i/tau_resolution)
      m_vec <- c(1,2)
      
    }
    
    curr_radius <- slack_scaling * norm(tau_vec, type="2")
    
    sigma_experiment <- generate_parameters_experiment_2(theta_true, g_true, beta_x_true)
    
    dataset <- generate_dataset(sigma_experiment, theta_true, N)
    
    A <- c(cov(dataset$Z_1, dataset$Y), cov(dataset$Z_2, dataset$Y))
    B <- c(cov(dataset$Z_1, dataset$X), cov(dataset$Z_2, dataset$X))
    
    feasible_region <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
    
    plot_intersection_2D(A, B, tau_vec, m_vec, curr_radius)
    
    feasible_points <- feasible_region$points
    feasible_intervals <- feasible_region$intervals
    
    curr_feasible_int_l2 <- l2_ball(A, B, tau_vec, slack_scaling)
    
    #print(curr_feasible_int_l2)
    
    if (!is.na(curr_feasible_int_l2$theta_lo)){
      
      #print("hi")
      
      results_l2 <- rbind(results_l2, list("sim_idx"=i, 
                                           "theta_lo"=curr_feasible_int_l2$theta_lo, "theta_hi"=curr_feasible_int_l2$theta_hi,
                                           "tau_lo"= tau_vec[1], "tau_hi"=tau_vec[2], "var_tau"=i/tau_resolution, "slack"=slack_scaling, "tau_radius"=curr_radius))
      
    }
    
    if (!is.na(feasible_points)) {
      for (point_idx in 1:length(feasible_points)){
        
        next_res <- list("sim_idx"=i, "interval"=FALSE, "feasible_start"=feasible_points[point_idx],
                         "feasible_end"=feasible_points[point_idx], "tau_lo"=tau_vec[1], "tau_hi"=tau_vec[2], "var_tau"=i/tau_resolution)
        
        results <- rbind(results, next_res)
      }}
    
    if (length(feasible_intervals) != 0) {
      # print(feasible_intervals)
      
      # print(length(feasible_intervals))
      
      for (interval_idx in 1:(length(feasible_intervals)/2)){
        
        #print(interval_idx)
        
        next_res <- list("sim_idx"=i, "interval"=TRUE, "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2], "tau_lo"=tau_vec[1], "tau_hi"=tau_vec[2], "var_tau"=i/tau_resolution)
        
        # print(next_res)
        results <- rbind(results, next_res)
        
      }}
    
    fwrite(results, './small dphi experiments 2/lazy/results.csv')
    fwrite(results_l2, './small dphi experiments 2/lazy/results_l2.csv')
    
  }
  
  print(feasible_region)
  
}

experiment_2()

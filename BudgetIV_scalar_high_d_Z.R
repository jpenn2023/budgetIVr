#library(data.table)
library(randcorr)
library(MASS)
library(ggplot2)
library(ggsci)
library(sisVIVE)
library(MendelianRandomization)
# source("Dataset simulations.R")

set.seed(42)

# Instead of Toeplitz R, sample R uniformly and set the b^* absolutely-least entries in the column R_{eps_y} to zero. 
# Then perform rejection sampling. 
# Let d_Z = 100 and b^* = 20. But run this against Gamma(b = 10, tau = 0).
# Sample fewer data for the statistic \hat{beta}_y than \hat{beta}_x and calculate \hat{SE}_y. 
# Sample 100 points for \hat{beta}_y and 10,000 for \hat{beta}_x. 

BudgetIV_scalar_exposure_feature <- function(
    beta_y, # Cross covariance Cov(Y, Z_vec)
    beta_phi, # Cross covariance Cov(Phi(X), Z_vec), now a vector
    tau_vec, # Degrees of violation of (AWE') (see manuscript). Ordered set.
    b_vec, # Ordered (increasing) set, demanding \sum_{i \in [d_Z]} { II (cov(Z_i, g_y) <= tau_i) } >= m_i, 
    # where II is the indicator function 
    delta_beta_y=NA # Box-shaped confidence set for beta_y, empty by default for point estimation
) {
  
  # Warnings
  if (is.unsorted(tau_vec)) {
    stop('Please input tau constraints in increasing order.')
  }
  else if (any(duplicated(tau_vec))){ 
    stop('The same tau constraint cannot be specified twice.')
  }
  else if (is.unsorted(b_vec) || any(duplicated(b_vec))) {
    stop('The vector m must be strictly increasing, please see the definition of boldface m in the manuscript.')
  }
  else if (length(beta_y) != length(beta_phi)) {
    stop('Cov(Y, Z) and Cov(Phi(X), Z) must be vectors of the same length for scalar Phi(X). Please call "BudgetIV" for treatment of vector Phi(X).')
  }
  
  d_Z <- length(beta_y)
  
  if(any(is.na(delta_beta_y))){delta_beta_y <- numeric(d_Z)}
  
  #
  ###    Redundant in this version since each delta_beta_y is non-zero    ### 
  #
  # Treat tau_vec[1] = 0 and tau_vec[1] > 0 seperately, since |A_i - B_i theta| = 0 is satisfied at exactly one point in the reduced problem.
  #intervals_start_flag <- 1
  #
  # Find all theta where the number of bucket constraints satisfied changes (either singularly at that point or changes as theta increases past that point) 
  #
  #theta_points <- NA
  #
  #
  # if (tau_vec[1] == 0){ 
  #   
  #   intervals_start_flag <- 2
  #   
  #   # theta_points <- numeric(length(J_non_sticky))
  #   theta_points <- numberic(d_Z)
  #   
  #   # theta_points(j) are defined when tau_1 = 0 
  #   for(j in 1:d_Z) {
  #     theta_points[j] <- beta_y[j] / beta_phi[j]
  #   }}
  #
  ###
  
  #
  # Bounds in theta where a specific beta_y[j] - beta_phi[j] * theta = +- (tau_vec[k]+delta_beta_y[j])
  #
  
  tau_intervals_lower <- matrix(nrow = (length(tau_vec)), ncol = d_Z)
  
  tau_intervals_upper <- matrix(nrow = (length(tau_vec)), ncol = d_Z)
  
  for (k in 1:length(tau_vec)){
    for (j in 1:d_Z){
      #j <- J_non_sticky[j_prime]
      
      tau_intervals_lower[k, j] <- beta_y[j] / beta_phi[j] - (tau_vec[k]+delta_beta_y[j]) / beta_phi[j]
      
      tau_intervals_upper[k, j] <- beta_y[j] / beta_phi[j] + (tau_vec[k]+delta_beta_y[j]) / beta_phi[j]
      
    }}
  
  #
  # Find the full feasible set of theta, including points and intervals along \mathbb{R}. 
  #
  
  # if (!is.na(theta_points)){
  #   
  #   # Points which can be tested one at a time:
  #   feasible_points <- rep(NA, length(theta_points))
  #   
  #   for(p in 1:length(theta_points)){
  #     
  #     if(validPoint_scalar(A, B, J_non_sticky, theta_points[p], m_vec_red, tau_vec)) {feasible_points[p] <- theta_points[p]}
  #     
  #     na.omit(feasible_points)
  #     sort(feasible_points)
  #     
  #   }
  # }
  # else {feasible_points <- NA}
  
  feasible_points <- NA
  
  # Feasible intervals
  feasible_intervals <- matrix(nrow = 0, ncol = 2)
  
  possible_bounds <- sort(c(c(tau_intervals_lower), c(tau_intervals_upper)))
  
  #print(validPoint_scalar(A, B, J_non_sticky, -0.2, m_vec_red, tau_vec))
  
  # Search through all possible points at which a feasible interval could begin or end:
  in_feasible <- FALSE
  
  #print("Just sorting and checking now:)")
  
  for (p in 1:(length(possible_bounds)-1)){
    
    #print(paste0("We're on to possible bound: ", p))
    
    curr_point <- (possible_bounds[p] + possible_bounds[p+1])/2
    
    # Make sure interval bounds are not biased by singularity at theta \in theta_points
    # if (curr_point %in% feasible_points){print("trouble here")
    #   
    #   while (curr_point %in% feasible_points){ curr_point <- (curr_point + possible_bounds[p])/2} }
    # 
    #print("no trouble here")
    
    # Keep track of whether inside or outside feasible region, including last theta at which a feasible interval was entered  
    # Add a new feasible interval once left
    curr_point_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point, b_vec, tau_vec, delta_beta_y)
    
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

validPoint_scalar <- function(beta_y, beta_phi, d_Z, theta, b_vec, tau_vec, delta_beta_y){
  
  b_to_fill <- b_vec
  
  for(j in 1:d_Z){
    
    #print("checking a point for ya")
    
    gamma_j <- beta_y[j] - beta_phi[j] * theta
    
    for(k in 1:length(tau_vec)){
      
      if(abs(gamma_j) <= tau_vec[k]+delta_beta_y[j]){ b_to_fill[k] <- b_to_fill[k] - 1 }
      
    }
  }
  
  #print(m_to_fill)
  
  return(all(b_to_fill <= 0))
  
}

generate_parameters_high_d_Z_rejection_method <- function(theta_true=1, number_valid=20, d_Z=100){
  
  R <- randcorr(d_Z+2)
  
  # Correlation between each candidate instrument and g(eps_y) = eps_y
  R_z_gy <- R[d_Z+2, 1:d_Z]
  
  smallest_indices <- order(abs(R_z_gy))[1:number_valid]
  
  R[d_Z+2, smallest_indices] <- 0
  R[smallest_indices, d_Z+2] <- 0
  
  while(any(eigen(R)$values <= 0)){
    
    R <- randcorr(d_Z+2)
    
    # Correlation between each candidate instrument and g(eps_y) = eps_y
    R_z_gy <- R[d_Z+2, 1:d_Z]
    
    smallest_indices <- order(abs(R_z_gy))[1:number_valid]
    
    R[d_Z+2, smallest_indices] <- 0
    R[smallest_indices, d_Z+2] <- 0
    
  }
  
  eta <- sqrt(rexp(d_Z+2, 1))
  
  #eta[d_Z+2] <- sqrt(rexp(1,d_Z))
  
  population_cov <- R * outer(eta, eta, "*")
  
  return(population_cov)
  
}


generate_parameters_high_d_Z_A3_method <- function(theta_true=1, number_valid=20, d_Z=100){
  
  R <- randcorr(d_Z+2)
  
  
  R[d_Z+2, ] <- 0
  R[, d_Z+2] <- 0
  R[d_Z+2, d_Z+2]
  
  R[d_Z+1, ] <- 0
  R[, d_Z+1] <- 0
  R[d_Z+1, d_Z+1] <- 1
  
  
  while(any(eigen(R)$values <= 0)){
    
    R <- randcorr(d_Z+2)
    
    R[d_Z+2, ] <- 0
    R[, d_Z+2] <- 0
    R[d_Z+2, d_Z+2] <- 1
    
    R[d_Z+1, ] <- 0
    R[, d_Z+1] <- 0
    R[d_Z+1, d_Z+1] <- 1
    
    gamma_A3 <- numeric(d_Z)
    
    gamma_A3[(number_valid+1):d_Z] <- 1
    
    delta <- sample(c(1, 2), size = d_Z, replace = TRUE)
    
  }
  
  eta <- sqrt(rexp(d_Z+2, 1))
  
  #eta[d_Z+1] <- sqrt(rexp(1, d_Z))
  
  population_cov <- R * outer(eta, eta, "*")
  
  return(list('population_cov'=population_cov, 'gamma_A3'=gamma_A3, 'delta'=delta))
  
}



generate_beta_statistics_A2 <- function(population_cov, theta_true, N_x=10000, N_y=100){
  
  d_Z = ncol(population_cov) - 2
  
  # Two sample approach, with N_y samples for \hat{beta}_y and N_x samples for \hat{beta}_x
  
  # \hat{beta}_y
  dataset_for_y <- mvrnorm(N_y, numeric(d_Z+2), Sigma = population_cov)
  
  Z_for_y <- dataset_for_y[, 1:d_Z]
  
  X_for_y <- dataset_for_y[, d_Z+1]
  Y <- theta_true * X_for_y + dataset_for_y[, d_Z+2]
  
  beta_y <- cov(Y, Z_for_y)
  
  var_y <- var(Y)
  
  var_z_for_y <- diag(cov(Z_for_y))
  
  SE_beta_y <- sqrt((var_y^2*(var_z_for_y * var_z_for_y) + beta_y * beta_y)/N_y)
  
  dataset_for_x <- mvrnorm(N_x, numeric(d_Z+2), Sigma = population_cov)
  
  Z_for_x <- dataset_for_x[, 1:d_Z]
  
  X <- dataset_for_x[, d_Z+1]
  
  beta_x <- cov(X, Z_for_x)
  
  var_x <- var(X)
  
  var_z_for_x <- diag(cov(Z_for_x))
  
  SE_beta_x <- sqrt((var_x^2*(var_z_for_x * var_z_for_x) + beta_x * beta_x)/(N_x-1))
  
  p_value_beta_x <- 2 * (1 - pnorm(abs(beta_x / SE_beta_x)))
  
  summary_stats <- data.frame('beta_y'=t(beta_y), 'SE_beta_y'=t(SE_beta_y), 'beta_x'=t(beta_x), 'p_value_beta_x'=t(p_value_beta_x))
  
  #summary_stats <- summary_stats[!is.nan(summary_stats$SE_beta_y), ]
  
  summary_stats <- summary_stats[summary_stats$p_value_beta_x <= 1e-8, ]
  
  return(summary_stats)
  
}


generate_beta_statistics_A3 <- function(population_cov, theta_true, N_x=10000, N_y=100, gamma_A3, delta){
  
  d_Z = ncol(population_cov) - 2
  
  # Two sample approach, with N_y samples for \hat{beta}_y and N_x samples for \hat{beta}_x
  
  # \hat{beta}_y
  dataset_for_y <- mvrnorm(N_y, numeric(d_Z+2), Sigma = population_cov)
  
  Z_for_y <- dataset_for_y[, 1:d_Z]
  
  X_for_y <- dataset_for_y[, d_Z+1]
  Y <- theta_true * X_for_y + dataset_for_y[, d_Z+2] + as.vector(Z_for_y %*% gamma_A3)
  
  beta_y <- cov(Y, Z_for_y)
  
  var_y <- var(Y)
  
  var_z_for_y <- diag(cov(Z_for_y))
  
  SE_beta_y <- sqrt((var_y^2*(var_z_for_y * var_z_for_y) + beta_y * beta_y)/N_y)
  
  dataset_for_x <- mvrnorm(N_x, numeric(d_Z+2), Sigma = population_cov)
  
  Z_for_x <- dataset_for_x[, 1:d_Z]
  
  X <- dataset_for_x[, d_Z+1] + as.vector(Z_for_x %*% delta)
  
  beta_x <- cov(X, Z_for_x)
  
  var_x <- var(X)
  
  var_z_for_x <- diag(cov(Z_for_x))
  
  SE_beta_x <- sqrt((var_x^2*(var_z_for_x * var_z_for_x) + beta_x * beta_x)/N_x)
  
  p_value_beta_x <- 2 * (1 - pnorm(abs(beta_x / SE_beta_x)))
  
  summary_stats <- data.frame('beta_y'=t(beta_y), 'SE_beta_y'=t(SE_beta_y), 'beta_x'=t(beta_x), 'SE_beta_x'=t(SE_beta_x), 'p_value_beta_x'=t(p_value_beta_x), 'gamma_A3'=gamma_A3, 'delta'=delta)
  
  #summary_stats <- summary_stats[!is.nan(summary_stats$SE_beta_y), ]
  
  #summary_stats <- summary_stats[summary_stats$p_value_beta_x <= 1e-8, ]
  
  return(summary_stats)
  
}


#simulated_experiment <- function(summary_stats, alpha=0.05){
#   
#   d_Z <- nrow(summary_stats)
#   
#   beta_x <- summary_stats$beta_x
#   
#   beta_y <- summary_stats$beta_y
#   
#   SE_beta_y <- summary_stats$SE_beta_y
#   
#   delta_beta_y <- qnorm(1 - alpha/(2*d_Z)) * SE_beta_y
#   
#   #
#   # The goal... (include the proposed instruments)
#   #
#   #results <- data.table("feasible_start"=numeric(), "feasible_end"=numeric(),
#   #                      "b"=numeric(), "proposed_instruments"=logical(d_Z))
#   
#   #
#   # But for now... (a lazy implementation)
#   #
#   results <- data.table("b"=numeric(), "feasible_start"=numeric(), "feasible_end"=numeric())
#   
#   tau_vec <- c(0)
#   
#   for(max_invalid in 0:d_Z){
#     
#     b_vec = c(d_Z - max_invalid)
#     
#     #print(paste0("Next b is: ", b_vec[1]))
#     
#     feasible_region <- BudgetIV_scalar_exposure_feature(beta_y, beta_x, tau_vec, b_vec, delta_beta_y)
#     
#     print(feasible_region)
#     
#     feasible_intervals <- feasible_region$intervals
#     
#     if(length(feasible_intervals) != 0){
#       
#       for (interval_idx in 1:(length(feasible_intervals)/2)){
#         
#         print(paste0("Next interval at b = ", b_vec[1], " is: ", interval_idx))
#         
#         #print(interval_idx)
#         
#         next_res <- list("b"=b_vec[1], "feasible_start"=feasible_intervals[interval_idx, 1],
#                          "feasible_end"=feasible_intervals[interval_idx, 2])
#         
#         # print(next_res)
#         results <- rbind(results, next_res)
#       }
#     }
#   }
#   
#   print(results)
#   
#   fwrite(results, './results_BudgetIV_high_d_Z_experiment.csv')
#   
# }


simulated_experiment <- function(summary_stats, alpha=0.05){

  d_Z <- nrow(summary_stats)

  beta_x <- summary_stats$beta_x

  beta_y <- summary_stats$beta_y

  SE_beta_y <- summary_stats$SE_beta_y

  delta_beta_y <- qnorm(1 - alpha/(2*d_Z)) * SE_beta_y

  #
  # The goal... (include the proposed instruments)
  #
  #results <- data.table("feasible_start"=numeric(), "feasible_end"=numeric(),
  #                      "b"=numeric(), "proposed_instruments"=logical(d_Z))

  #
  # But for now... (a lazy implementation)
  #
  results <- data.table("b"=numeric(), "feasible_start"=numeric(), "feasible_end"=numeric())

  tau_vec <- c(0)

  for(max_invalid in 0:d_Z-1){

    b_vec = c(d_Z - max_invalid)

    #print(paste0("Next b is: ", b_vec[1]))

    feasible_region <- BudgetIV_scalar_exposure_feature(beta_y, beta_x, tau_vec, b_vec, delta_beta_y)

    print(feasible_region)

    feasible_intervals <- feasible_region$intervals

    if(length(feasible_intervals) != 0){

      for (interval_idx in 1:(length(feasible_intervals)/2)){

        print(paste0("Next interval at b = ", b_vec[1], " is: ", interval_idx))

        #print(interval_idx)

        next_res <- list("b"=b_vec[1], "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2])

        # print(next_res)
        results <- rbind(results, next_res)
      }
    }
  }

  print(results)

  fwrite(results, './results_BudgetIV_high_d_Z_experiment_A3.csv')

}

benchmark_estimates <- function(summary_stats, alpha=0.05){
  
  mr_summary_stats <- mr_input(bx = summary_stats$beta_x,
                               by = summary_stats$beta_y,
                               bxse = summary_stats$SE_beta_x,
                               byse = summary_stats$SE_beta_y)
  
  mbe_solution <- mr_mbe(mr_summary_stats, alpha=alpha)
  
  median_solution <- mr_median(mr_summary_stats, alpha=alpha)
  
  mr_egger_solution <- mr_egger(mr_summary_stats, alpha=alpha)
  
  ivw_solution <- mr_ivw(mr_summary_stats, alpha=alpha)
  
  return(list('mbe'=mbe_solution, 'median'=median_solution, 'mr_egger'=mr_egger_solution, 'ivw'=ivw_solution))
  
}



d_Z <- 100
number_valid <- 30
theta_true <- 1
N_x <- 1000000
N_y <- 100000

#population_parameters <- generate_parameters_high_d_Z_A3_method(theta_true, number_valid, d_Z)

population_cov <- population_parameters$population_cov

gamma_A3 <- population_parameters$gamma_A3

delta <- population_parameters$delta

#summary_stats <- generate_beta_statistics_A3(population_cov, theta_true, N_x, N_y, gamma_A3, delta)

dataset_for_sisvive <- mvrnorm(1000000, numeric(d_Z+2), Sigma = population_cov)

Z_sisvive <- dataset_for_sisvive[, 1:d_Z]
X_sisvive <- dataset_for_sisvive[, d_Z+1]
Y_sisvive <- theta_true * X_sisvive + dataset_for_sisvive[, d_Z+2] + as.vector(Z_sisvive %*% gamma_A3)

sisVIVE_estimator <- cv.sisVIVE(Y=Y_sisvive, D=X_sisvive, Z=Z_sisvive)

sisVIVE_estimate <- sisVIVE_estimator$beta

simulated_experiment(summary_stats, 0.05)

benchmarks <- benchmark_estimates(summary_stats, 0.05)

d_Z <- nrow(summary_stats)

feasible_regions <- read.csv('./results_BudgetIV_high_d_Z_experiment_A3.csv')

g <- ggplot() +
  # geom_segment(data = feasible_intervals_l2,
  # aes(x=var_tau, xend=var_tau, y=theta_lo, yend=theta_hi,
  # color = Constraint), linewidth = 3) +
  # geom_segment(data = feasible_intervals_l1,
  #              aes(x=var_tau, xend=var_tau, y=theta_lo, yend=theta_hi,
  #                  color = Constraint), linewidth = 3) +
  geom_segment(data = feasible_regions,
               aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 2) +
  geom_hline(yintercept = 1, color = "blue", linetype = "31", linewidth=1) +
  geom_vline(xintercept = d_Z-number_valid, color = "blue", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mbe@CIUpper, color = "orange", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mbe@CILower, color = "orange", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$median@CIUpper, color = "red", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$median@CILower, color = "red", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mr_egger@CIUpper.Est, color = "gray", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mr_egger@CILower.Est, color = "gray", linetype = "31", linewidth=1) +
  geom_hline(yintercept = sisVIVE_estimate, color = "darkgreen", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$ivw@CIUpper, color = "violet", linetype = "31", linewidth=0.5) +
  geom_hline(yintercept = benchmarks$ivw@CILower, color = "violet", linetype = "31", linewidth=0.5) +
  # scale_color_manual(labels = c('Budget',
  #                               expression(italic(L)[1]),
  #                               expression(italic(L)[2])),
  #                    values = c('#FF7F0EFF', '#1F77B4FF', 'lightblue')) +
  # scale_x_continuous(limits = c(1, 150), breaks = seq(0, 150, 10)) +
  coord_cartesian(xlim=c(0,100), ylim=c(-2,2)) +
  #scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1)) +
  labs(x = paste0("Maximum number of plausibly invalid instruments, b (d_Z = ", d_Z,")"),
       y = paste("(95% coverage over) the set of feasible theta")) +
  #annotate("text", x = 1, y = 1, label = expression(paste(theta,"*")),
  #         hjust = 0.3, vjust = -0.5, color = "blue", size = 10) +
  #annotate("text", x = 66, y = 1.7, label = expression(paste(b,"*")),
  #         hjust = 0.6, vjust = -0.5, color = "blue", size = 10) +
  annotate("text", x = 1, y = -1.3, label = expression(paste("MBE")),
           hjust = 0, vjust = -0.5, color = "orange", size = 10) +
  annotate("text", x = 35, y = -1.3, label = expression(paste("Ground truth")),
           hjust = 0, vjust = -0.5, color = "blue", size = 10) +
  annotate("text", x = 1, y = -1.6, label = expression(paste("MR Median")),
           hjust = 0, vjust = -0.5, color = "red", size = 10) +
  annotate("text", x = 35, y = -1.6, label = expression(paste("sisVIVE")),
           hjust = 0, vjust = -0.5, color = "darkgreen", size = 10) +
  annotate("text", x = 1, y = -2, label = expression(paste("MR Egger")),
           hjust = 0, vjust = -0.5, color = "gray", size = 10) +
  annotate("text", x = 35, y = -1.915, label = expression(paste("IVW")),
           hjust = 0, vjust = -0.5, color = "violet", size = 10) +
  # scale_color_manual(values = c("Group 1" = "red", "Group 2" = "blue", "Group 3" = "orange", "MR Egger" = "gray", "Budget IV" = "black"),
  #                    labels = c("h", "i", "f", "MR Egger", "Budget IV")) +
  #scale_linewidth_continuous(guide = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'bottom')
# ggsave('./high_d_Z_eps_y/MANY_SAMPLES_high_d_Z_100_theta_true_1_40_valid_instruments.pdf', width = 8)
ggsave('./high_d_Z_A3_violation/MANY_SAMPLES_high_d_Z_100_theta_true_1_30_valid_instruments_zoomed_in_labelled.pdf', width = 8)

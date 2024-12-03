#
#
#
#
# Special case of d_{Phi} = 1 yields very efficient, polytime solution.
#
# Call the function experiment_in_manuscript() to reproduce the results.
#
#

library(data.table)
library(graphics)
library(stats)
library(MendelianRandomization)
library(ggplot2)
library(ggsci)
# source("Dataset simulations.R")

set.seed(42)

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

real_data_experiment_lipid_fractions <- function(alpha){
  
  MBE_data <- read.csv('MBE_data.csv')
  
  candidatesHDL = MBE_data[MBE_data$pHDL <= 1e-8, ]
  
  candidate_labels <- candidatesHDL$rsID
  d_Z <- length(candidate_labels)
  
  beta_x <- candidatesHDL$betaHDL
  
  beta_y <- candidatesHDL$betaCAD
  
  SE_beta_y <- abs(beta_y) / qnorm(1-candidatesHDL$pCAD/2)
  
  delta_beta_y <- qnorm(1 - alpha/(2*d_Z)) * SE_beta_y
  
  new_dataset <- data.frame('candidate_labels' = candidate_labels, 'beta_x' = beta_x, 
                            'beta_y' = beta_y, 'delta_beta_y' = delta_beta_y)
  
  new_dataset <- new_dataset[!is.nan(delta_beta_y), ]
  
  print(new_dataset)
  
  #print(new_dataset)
  
  # Basic idea: vary b and form a plot similar to the figure 1 in the paper
  
  beta_y <- new_dataset$beta_y
  beta_phi <- new_dataset$beta_x
  delta_beta_y <- new_dataset$delta_beta_y
  
  d_Z <- length(new_dataset$beta_y)
  
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
  
  for(max_invalid in 0:d_Z){
    
    b_vec = c(d_Z - max_invalid)
    
    #print(paste0("Next b is: ", b_vec[1]))
    
    feasible_region <- BudgetIV_scalar_exposure_feature(beta_y, beta_phi, tau_vec, b_vec, delta_beta_y)
    
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
  
  fwrite(results, './results_BudgetIV_lipid_fractions_HDL_CAD.csv')
  
}

real_data_experiment_lipid_fractions(0.05)

feasible_regions <- read.csv('./results_BudgetIV_lipid_fractions_HDL_CAD.csv')

d_Z = 88


g <- ggplot() +
  # geom_segment(data = feasible_intervals_l2,
               # aes(x=var_tau, xend=var_tau, y=theta_lo, yend=theta_hi,
                   # color = Constraint), linewidth = 3) +
  # geom_segment(data = feasible_intervals_l1,
  #              aes(x=var_tau, xend=var_tau, y=theta_lo, yend=theta_hi,
  #                  color = Constraint), linewidth = 3) +
  geom_segment(data = feasible_regions,
               aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 1) +
  # geom_hline(yintercept = 1, color = "red", linetype = "21", linewidth=1) +
  # geom_vline(xintercept = 2, color = "black", linetype = "21", linewidth=1) +
  # scale_color_manual(labels = c('Budget',
  #                               expression(italic(L)[1]),
  #                               expression(italic(L)[2])),
  #                    values = c('#FF7F0EFF', '#1F77B4FF', 'lightblue')) +
  scale_x_continuous(limits = c(1, 90), breaks = seq(0, 90, 5)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1)) +
  labs(x = expression(paste("Maximum number of plausibly invalid instruments allowed, ", b, " (d_Z = 74)")),
       y = expression(paste("(95% coverage) ATE of HDL lipid fraction of log-odds CAD risk"))) +
  # annotate("text", x = 1, y = 1, label = expression(paste(theta,"*")),
  #          hjust = 0.3, vjust = -0.5, color = "red", size = 10) +
  # annotate("text", x = 2.5, y = -2, label = expression(paste(tau,"*")),
  #          hjust = 0.6, vjust = -0.5, color = "black", size = 10) +
  #scale_linewidth_continuous(guide = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'bottom')
ggsave('./lipid_fractions_HDL_CAD.pdf', width = 8)

#plot_real_experiment()
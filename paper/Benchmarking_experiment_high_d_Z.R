
budgetIV_scalar <- function(
    beta_y,
    beta_phi,
    tau_vec,
    b_vec,
    delta_beta_y=NULL,
    bounds_only=FALSE
) {
  
  if(is.matrix(beta_y) && is.numeric(beta_y)){ 
    if(nrow(beta_y) == 1){
      beta_y <- as.vector(beta_y)
    }
    else{
      stop("Argument 'beta_y' must be a vector or a matrix with one row.")
    }
  }
  else if(!is.numeric(beta_y)){
    stop("Argument 'beta_y' must have numeric entries and be input as a matrix or vector.")
  }
  
  if(is.matrix(beta_phi)){ 
    if(nrow(beta_phi) == 1 && is.numeric(beta_phi)){
      beta_y <- as.vector(beta_phi)
    }
    else{
      stop("Argument 'beta_phi' must be a vector or a matrix with one row. 
           Please use budgetIV for non-scalar 'phi_basis'.")
    }
  }
  else if(!is.vector(beta_phi)){
    stop("Argument 'beta_phi' must have numeric entries and be input as a matrix or vector.")
  }
  
  if (is.null(delta_beta_y)){
    warning("Argument 'delta_beta_y' not specified. 
            No confidence bounds for agument 'beta_y' given: treating 'beta_y' as an oracle summary statistic.")
    delta_beta_y <- numeric(d_Z)
  }
  
  if(is.matrix(delta_beta_y) && is.numeric(delta_beta_y)){ 
    if(nrow(beta_y) == 1){
      beta_y <- as.vector(beta_y)
    }
    else{
      stop("Argument 'delta_beta_y' must be a vector or a matrix with one row.")
    }
  }
  else if(!is.numeric(delta_beta_y)){
    stop("Argument 'delta_beta_y' must have numeric entries and be input as a matrix or vector.")
  }
  
  d_Z <- ncol(beta_y)
  
  # Error messages
  if(!is.vector(beta_y)){
    stop("Argument 'beta_y' must be a vector or single-row matrix.")
  }
  else if(!is.numeric(beta_y)){
    stop("Argument 'beta_y' must have numeric entries.")
  }
  else if(!is.numeric(beta_phi)){
    stop("Argument 'beta_phi' must have numeric entries.")
  }
  else if (length(beta_y) != length(beta_phi)) {
    stop("Arguments 'beta_y' and 'beta_phi' must have the same length/number of columns.")
  }
  else if(!is.vector(tau_vec)){
    stop("Argument 'tau_vec' must be a vector. Use tau_vec = c(threshold_value) for a single budget constraint (e.g., tau_vec = c(0) for an L_0-norm constraint).")
  }
  else if(!is.numeric(tau_vec)){
    stop("Argument 'tau_vec' must have numeric entries.")
  }
  else if(!all(tau_vec >= 0)){
    stop("Argument 'tau_vec' must have positive entries.")
  }
  else if (is.unsorted(tau_vec)) {
    stop("Argument 'tau_vec' must have entries in increasing order.")
  }
  else if (any(duplicated(tau_vec))){ 
    stop("Argument 'tau_vec' must be strictly increasing, i.e., with no repeated entries.")
  }
  else if(!is.vector(b_vec)){
    stop("Argument 'b_vec' must be a vector. Use b_vec = c(budget_value) for a single budget constraint.")
  }
  else if(!is.numeric(b_vec)){
    stop("Argument 'b_vec' must have numeric entries.")
  }
  else if(!all(b_vec == as.integer(b_vec)) & is.numeric(b_vec)){
    stop("Argument 'b_vec' must have integer entries.")
  }
  else if(!all(b_vec > 0)){
    stop("Argument 'b_vec' must have entries strictly greater than zero.")
  }
  if (is.unsorted(b_vec)) {
    stop("Argument 'b_vec' must have entries in increasing order.")
  }
  else if (any(duplicated(b_vec))){ 
    stop("Argument 'b_vec' must be strictly increasing, i.e., with no repeated entries.")
  }
  else if(any(delta_beta_y < 0)){
    stop("Argument 'delta_beta_y' must have positive entries.")
  }
  else if (length(beta_y) != length(beta_phi)) {
    stop("Arguments 'beta_y' and 'beta_phi' must be vectors of the same length for scalar Phi(X). Please call 'budgetIV' for treatment of vector Phi(X).")
  }
  else if (length(delta_beta_y) != length(beta_y)){
    stop("Argument 'delta_beta_y', if given, must be of the same length as beta_y.")
  }
  else if (!is.logical(bounds_only) || length(bounds_only) != 1 || is.na(bounds_only)){
    stop("Argument 'bounds_only' must be a single TRUE or FALSE value only.")
  }
  
  print("Made it past the error checks")
  
  d_Z <- length(beta_y)
  
  tau_intervals_lower <- matrix(nrow = (length(tau_vec)), ncol = d_Z)
  
  tau_intervals_upper <- matrix(nrow = (length(tau_vec)), ncol = d_Z)
  
  for (k in 1:length(tau_vec)){
    for (j in 1:d_Z){
      
      tau_intervals_lower[k, j] <- beta_y[j] / beta_phi[j] - (tau_vec[k]+delta_beta_y[j]) / beta_phi[j]
      
      tau_intervals_upper[k, j] <- beta_y[j] / beta_phi[j] + (tau_vec[k]+delta_beta_y[j]) / beta_phi[j]
      
    }}
  
  possible_bounds <- sort(c(c(tau_intervals_lower), c(tau_intervals_upper)))
  
  in_feasible <- FALSE
  
  if (bounds_only == TRUE){
    
    causal_effect_bounds <- data.table(
      "is_point"=logical(),
      "lower_bound"=numeric(),
      "upper_bound"=numeric()
    )
    
    for (p in 1:(length(possible_bounds)-1)){
      
      curr_point <- possible_bounds[p]
      
      curr_point_interval <- (possible_bounds[p] + possible_bounds[p+1])/2
      
      curr_interval_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point_interval, b_vec, tau_vec, delta_beta_y)
      
      curr_point_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point, b_vec, tau_vec, delta_beta_y)
      
      if(curr_interval_feasible && !in_feasible){
        
        print("found a region")
        
        last_feasible_opening <- curr_point
        in_feasible <- TRUE
        
      }
      
      else if(!curr_interval_feasible && in_feasible){
        
        in_feasible <- FALSE
        
        new_interval <- data.table(
          "is_point"=FALSE,
          "lower_bound"=last_feasible_opening,
          "upper_bound"=curr_point
        )
        
        causal_effect_bounds <- rbind(causal_effect_bounds, new_interval)
        
      }
      
      else if(!in_feasible && curr_point_feasible){
        
        new_point <- data.table(
          "is_point"=TRUE,
          "lower_bound"=curr_point,
          "upper_bound"=curr_point
        )
        
        causal_effect_bounds <- rbind(causal_effect_bounds, new_point)
        
      }
      
    }
    
    curr_point <- possible_bounds[length(possible_bounds)]
    
    if(in_feasible){ 
      
      new_interval <- data.table(
        "is_point"=FALSE,
        "lower_bound"=last_feasible_opening,
        "upper_bound"=curr_point
      )
      
      causal_effect_bounds <- rbind(causal_effect_bounds, new_interval)
      
    }
    
    else if(validPoint_scalar(beta_y, beta_phi, d_Z, curr_point, b_vec, tau_vec, delta_beta_y)){
      
      new_point <- data.table(
        "is_point"=TRUE,
        "lower_bound"=curr_point,
        "upper_bound"=curr_point
      )
      
      causal_effect_bounds <- rbind(causal_effect_bounds, new_point)
      
    }
    
    return(causal_effect_bounds)
    
  }
  
  else if (bounds_only==FALSE){
    
    causal_effect_bounds <- data.table(
      "is_point"=logical(),
      "lower_bound"=numeric(),
      "upper_bound"=numeric(),
      "budget_assignment" = list()
    )
    
    # print(causal_effect_bounds)
    
    for (p in 1:(length(possible_bounds)-1)){
      
      curr_point <- possible_bounds[p]
      
      curr_point_interval <- (curr_point + possible_bounds[p+1])/2
      
      curr_interval_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point_interval, b_vec, tau_vec, delta_beta_y)
      
      curr_point_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point, b_vec, tau_vec, delta_beta_y)
      
      if(curr_interval_feasible){
        
        curr_interval_budgets <- eval_budgets(beta_y, beta_phi, d_Z, curr_point_interval, tau_vec, delta_beta_y)
        
        last_interval_feasible <- TRUE
        
        new_interval <- data.table(
          "is_point" = FALSE,
          "lower_bound" = curr_point,
          "upper_bound" = possible_bounds[p+1],
          "budget_assignment" = list(curr_interval_budgets)
        )
        
        # print(new_interval)
        
        causal_effect_bounds <- rbind(causal_effect_bounds, new_interval)
        
      }
      
      else if(curr_point_feasible && !last_interval_feasible){
        
        curr_point_budgets <- eval_budgets(beta_y, beta_phi, d_Z, curr_point, tau_vec, delta_beta_y)
        
        new_point <- data.table(
          "is_point" = TRUE,
          "lower_bound" = curr_point,
          "upper_bound" = possible_bounds[p+1],
          "budget_assignment" = list(curr_point_budgets)
        )
        
        
        causal_effect_bounds <- rbind(causal_effect_bounds, new_point)
        
      }
      
      else { 
        
        last_interval_feasible <- FALSE 
      }
      
    }
    
    curr_point <- possible_bounds[length(possible_bounds)]
    
    curr_point_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point, b_vec, tau_vec, delta_beta_y)
    
    if(curr_point_feasible && !last_interval_feasible){
      
      curr_point_budgets <- eval_budgets(beta_y, beta_phi, d_Z, curr_point, tau_vec, delta_beta_y)
      
      new_point <- data.table(
        "is_point" = TRUE,
        "lower_bound" = curr_point,
        "upper_bound" = curr_point,
        "budget_assignment" = list(curr_point_budgets)
      )
      
      causal_effect_bounds <- rbind(causal_effect_bounds, new_point)
      
    }
    
    return(causal_effect_bounds)
    
    
  }
  
}

validPoint_scalar <- function(beta_y, beta_phi, d_Z, theta, b_vec, tau_vec, delta_beta_y){
  
  b_to_fill <- b_vec
  
  for(i in 1:d_Z){
    
    beta_theta_i <- beta_phi[i] * theta
    
    tol <- .Machine$double.eps * max(beta_y[i], beta_theta_i)
    
    gamma_i <- beta_y[i] - beta_phi[i] * theta
    
    for(tau_index in 1:length(tau_vec)){
      
      if(abs(gamma_i) <= tau_vec[tau_index]+delta_beta_y[i] + tol){ b_to_fill[tau_index] <- b_to_fill[tau_index] - 1 }
      
    }
  }
  
  #print(m_to_fill)
  
  return(all(b_to_fill <= 0))
  
}



# Return the budget assignment at a single point, corresponding to a single possible 'theta'. 
# For delta_beta_y not NA, the output is the smallest (in terms of component-wise partial order)

eval_budgets <- function(beta_y, beta_phi, d_Z, theta, tau_vec, delta_beta_y){
  
  gammas <- beta_y - theta * beta_phi
  
  budget_assignments <- rep(0, d_Z)
  
  for (i in 1:d_Z){
    
    tol <- .Machine$double.eps * max(beta_y[i], beta_phi[i] * theta)
    
    tau_index = 1
    
    for(tau_index in 1:length(tau_vec) ){
      
      if (abs(gammas[i]) <= tau_vec[tau_index] + delta_beta_y[i] + tol){budget_assignments[i] <- tau_index }
      
    }
    
    if(budget_assignments[i] == 0){ budget_assignments[i] <- length(tau_vec)+1 }
    
  }
  
  return(budget_assignments)
  
}





















# library(randcorr)
# library(MASS)
# library(ggplot2)
# library(ggsci)
# library(sisVIVE)
# library(MendelianRandomization)
# library(boot)
# library(data.table)
# library(ggsci)

# # source("Dataset simulations.R")
# 

# set.seed(42)
# 
# # Instead of Toeplitz R, sample R uniformly and set the b^* absolutely-least entries in the column R_{eps_y} to zero. 
# # Then perform rejection sampling. 
# # Let d_Z = 100 and b^* = 20. But run this against Gamma(b = 10, tau = 0).
# # Sample fewer data for the statistic \hat{beta}_y than \hat{beta}_x and calculate \hat{SE}_y. 
# # Sample 100 points for \hat{beta}_y and 10,000 for \hat{beta}_x. 
# 
set.seed(42)

# setwd("C:/Users/k23067841/Downloads/budgetivr/R")
# 
# source("budgetIV_scalar.R")

setwd("C:/Users/k23067841/Downloads/Temp experiments BIV/linear high d_Z exp")

generate_parameters_high_d_Z_A3_method <- function(theta_true=1, number_valid=20, d_Z=100){
  
  R <- diag(nrow=(d_Z+2), ncol=(d_Z+2))
  
  R[d_Z+1, d_Z+2] <- runif(1,-1,1)
  
  gamma_A3 <- numeric(d_Z)
  
  gamma_A3[(number_valid+1):d_Z] <- runif(d_Z-number_valid,1,2)
  
  delta <- runif(d_Z,1,2)
  
  eta <- sqrt(rexp(d_Z+2, 1))
  
  population_cov <- R * outer(eta, eta, "*")
  
  return(list('population_cov'=population_cov, 'gamma_A3'=gamma_A3, 'delta'=delta))
  
}

generate_beta_statistics_A3 <- function(population_cov, theta_true, N_x, N_y, gamma_A3, delta){
  
  d_Z = ncol(population_cov) - 2
  
  dataset_for_y <- mvrnorm(N_y, numeric(d_Z+2), Sigma = population_cov)
  
  Z_for_y <- dataset_for_y[, 1:d_Z]
  
  X_for_y <- dataset_for_y[, d_Z+1] + as.vector(Z_for_y %*% delta)
  Y <- theta_true * X_for_y + dataset_for_y[, d_Z+2] + as.vector(Z_for_y %*% gamma_A3)
  
  beta_y <- cov(Y, Z_for_y)
  
  var_y <- var(Y)
  
  var_z_for_y <- diag(cov(Z_for_y))
  
  SE_beta_y <- sqrt((var_y*var_z_for_y + beta_y * beta_y)/N_y)
  
  dataset_for_x <- mvrnorm(N_x, numeric(d_Z+2), Sigma = population_cov)
  
  Z_for_x <- dataset_for_x[, 1:d_Z]
  
  X <- dataset_for_x[, d_Z+1] + as.vector(Z_for_x %*% delta)
  
  Y_for_x <- theta_true * X + dataset_for_x[, d_Z+2] + as.vector(Z_for_x %*% gamma_A3)
  
  beta_x <- cov(X, Z_for_x)
  
  var_x <- var(X)
  
  cov_z_for_x <- cov(Z_for_x)
  var_z_for_x <- diag(cov_z_for_x)
  
  corr_z_for_x <- cov_z_for_x * outer(1/var_z_for_x, 1/var_z_for_x, '*')
  
  SE_beta_x <- sqrt((var_x*(var_z_for_x) + beta_x * beta_x)/N_x)
  
  p_value_beta_x <- 2 * (1 - pnorm(abs(beta_x / SE_beta_x)))
  
  summary_stats <- data.frame('beta_y'=t(beta_y), 'SE_beta_y'=t(SE_beta_y), 'beta_x'=t(beta_x), 'SE_beta_x'=t(SE_beta_x), 'p_value_beta_x'=t(p_value_beta_x), 'gamma_A3'=gamma_A3, 'delta'=delta)
  
  return(list('summary_stats'=summary_stats, 
              # 'SS_massive'=SS_massive, 
              'empirical_corr_z'=corr_z_for_x))
  
}

simulated_experiment <- function(summary_stats, alpha=0.05){
  
  d_Z <- nrow(summary_stats)
  
  beta_x <- summary_stats$beta_x
  
  beta_y <- summary_stats$beta_y
  
  SE_beta_y <- summary_stats$SE_beta_y
  
  delta_beta_y <- qnorm(1 - alpha/(2*d_Z)) * SE_beta_y
  
  results <- data.table("b"=numeric(), "feasible_start"=numeric(), "feasible_end"=numeric())
  
  tau_vec <- c(0)
  
  for(max_invalid in 0:d_Z-1){
    
    b_vec = c(d_Z - max_invalid)
    
    feasible_region <- budgetIV_scalar(beta_y = beta_y,
                                       beta_phi = beta_x,
                                       tau_vec = tau_vec,
                                       b_vec = b_vec,
                                       delta_beta_y = delta_beta_y)
    
    print(feasible_region)
    
    feasible_intervals <- cbind(feasible_region$lower_bound, feasible_region$upper_bound) 
    
    if(nrow(feasible_intervals) != 0){
      
      for (interval_idx in 1:(length(feasible_intervals)/2)){
        
        print(paste0("Next interval at b = ", b_vec[1], " is: ", interval_idx))
        
        next_res <- list("b"=b_vec[1], "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2])
        
        results <- rbind(results, next_res)
      }
    }
  }
  
  # print(results)
  
  fwrite(results, 'results_BudgetIV_high_d_Z_experiment_A3.csv')
  
}

simulated_experiment_oracle <- function(summary_stats, alpha=0.05){
  
  d_Z <- nrow(summary_stats)
  
  beta_x <- summary_stats$beta_x
  
  beta_y <- summary_stats$beta_y
  
  SE_beta_y <- summary_stats$SE_beta_y
  
  delta_beta_y <- qnorm(1 - alpha/(2*d_Z)) * SE_beta_y
  
  results <- data.table("b"=numeric(), "feasible_start"=numeric(), "feasible_end"=numeric())
  
  tau_vec <- c(0.01)
  
  for(max_invalid in 0:d_Z-1){
    
    b_vec = c(d_Z - max_invalid)
    
    feasible_region <- budgetIV_scalar(beta_y = beta_y,
                                       beta_phi = beta_x,
                                       tau_vec = tau_vec,
                                       b_vec = b_vec,
                                       delta_beta_y = delta_beta_y)
    
    # print(feasible_region)
    
    feasible_intervals <- cbind(feasible_region$lower_bound, feasible_region$upper_bound) 
    
    if(nrow(feasible_intervals) != 0){
      
      for (interval_idx in 1:(length(feasible_intervals)/2)){
        
        print(paste0("Next interval at b = ", b_vec[1], " is: ", interval_idx))
        
        next_res <- list("b"=b_vec[1], "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2])
        
        results <- rbind(results, next_res)
      }
    }
  }
  
  # print(results)
  
  fwrite(results, 'oracle_results_BudgetIV_high_d_Z_experiment_A3.csv')
  
}

benchmark_estimates <- function(summary_stats, corr_z, alpha=0.05){
  
  d_Z <- nrow(summary_stats)
  
  mr_summary_stats <- mr_input(bx = summary_stats$beta_x,
                               by = summary_stats$beta_y,
                               bxse = summary_stats$SE_beta_x,
                               byse = summary_stats$SE_beta_y,
                               correlation = corr_z)
  
  mbe_solution <- mr_mbe(mr_summary_stats, alpha=alpha)
  
  median_solution <- mr_median(mr_summary_stats, alpha=alpha)
  
  mr_egger_solution <- mr_egger(mr_summary_stats, alpha=alpha)
  
  ivw_solution <- mr_ivw(mr_summary_stats, alpha=alpha)
  
  return(list('mbe'=mbe_solution, 'median'=median_solution, 'mr_egger'=mr_egger_solution, 'ivw'=ivw_solution))
  
}

# study_sizes <- c(20, 40, 60, 80, 100)
# nums_valid <- 0.3 * study_sizes
# 
# studies <- 5

d_Z <- 100
number_valid <- 30

theta_true <- 1
N_x <- 1e6
N_y <- 1e5

# d_Z <- study_sizes[study_idx]
# number_valid <- nums_valid[study_idx]

population_parameters <- generate_parameters_high_d_Z_A3_method(theta_true, number_valid, d_Z)

population_cov <- population_parameters$population_cov

gamma_A3 <- population_parameters$gamma_A3

delta <- population_parameters$delta

summary_stats_all <- generate_beta_statistics_A3(population_cov, theta_true, N_x, N_y, gamma_A3, delta)

summary_stats <- summary_stats_all$summary_stats

empirical_corr_z <- summary_stats_all$empirical_corr_z

var_z_true <- population_cov[1:d_Z,1:d_Z]

cov_eps_x_Z_true <- population_cov[(1+d_Z), 1:d_Z]

cov_eps_y_Z_true <- population_cov[(2+d_Z), 1:d_Z]

beta_x_true <- delta %*% var_z_true + cov_eps_x_Z_true

SE_beta_x_true <- numeric(length(beta_x_true))

p_value_beta_x_true <- numeric(length(beta_x_true))

beta_y_true <- theta_true * beta_x_true + gamma_A3 %*% var_z_true + cov_eps_y_Z_true

SE_beta_y_true <- numeric(length(beta_y_true))

oracle_stats <- data.frame('beta_y'=t(beta_y_true), 'SE_beta_y'=SE_beta_y_true, 'beta_x'=t(beta_x_true), 'SE_beta_x'=SE_beta_x_true, 'p_value_beta_x'=p_value_beta_x_true, 'gamma_A3'=gamma_A3, 'delta'=delta)

simulated_experiment(summary_stats, 0.05)

simulated_experiment_oracle(oracle_stats, 0.05)

benchmarks <- benchmark_estimates(summary_stats, empirical_corr_z, 0.05)

feasible_regions <- read.csv('results_BudgetIV_high_d_Z_experiment_A3.csv')

oracle_regions <- read.csv('oracle_results_BudgetIV_high_d_Z_experiment_A3.csv')

massive_results <- readRDS('posterior.rds')

massive_beta_results <- massive_results$betas

betas_posterior <- data.frame(value = massive_beta_results)

posterior_plot <- ggplot() +
  geom_histogram(data=betas_posterior, aes(x = value), bins=80) +
  theme_classic()

ggsave(paste0('posterior_hist.pdf'), width = 4)

my_colours <- pal_lancet()(7)

num_MBE <- 16

num_median <- 20

num_mr_egger <- 15

num_ivw <- 18

num_massive <- 40

benchmark_plot <- ggplot() +
  geom_hline(yintercept = theta_true, color = my_colours[1], linewidth=2) +
  geom_vline(xintercept = d_Z-number_valid, color = my_colours[1], linewidth=2) +
  geom_segment(data = feasible_regions,
               aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 5, color = "gray10") +
  geom_segment(data = oracle_regions, color = my_colours[1],
               aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 5) +
  # geom_hline(yintercept = benchmarks$mbe@CIUpper, color = my_colours[6], linetype = "31", linewidth=1) +
  # geom_hline(yintercept = benchmarks$mbe@CILower, color = my_colours[6], linetype = "31", linewidth=1) +
  # geom_hline(yintercept = benchmarks$median@CIUpper, color = my_colours[3], linetype = "27", linewidth=1.5) +
  # geom_hline(yintercept = benchmarks$median@CILower, color = my_colours[3], linetype = "27", linewidth=1.5) +
  # geom_hline(yintercept = benchmarks$mr_egger@CIUpper.Est, color = my_colours[4], linetype = "25",  linewidth=1) +
  # geom_hline(yintercept = benchmarks$mr_egger@CILower.Est, color = my_colours[4], linetype = "25", linewidth=1) +
  # geom_point(aes(x = seq(40, 100, length.out = num_mr_egger), y = rep(benchmarks$mr_egger@CIUpper.Est, num_mr_egger)), shape = 15, size = 2, fill = my_colours[4], color = my_colours[4]) + 
  # geom_point(aes(x = seq(40, 100, length.out = num_mr_egger), y = rep(benchmarks$mr_egger@CILower.Est, num_mr_egger)), shape = 15, size = 2, fill = my_colours[4], color = my_colours[4]) +
  # geom_hline(yintercept = benchmarks$ivw@CIUpper, color = my_colours[5], linetype = "31", linewidth=0.5) +
  # geom_hline(yintercept = benchmarks$ivw@CILower, color = my_colours[5], linetype = "31", linewidth=0.5) +
  # geom_point(aes(x = seq(40, 100, length.out = num_ivw), y = rep(benchmarks$ivw@CIUpper, num_ivw)), shape = 17, size = 2, fill = my_colours[5], color = my_colours[5]) + 
  # geom_point(aes(x = seq(40, 100, length.out = num_ivw), y = rep(benchmarks$ivw@CILower, num_ivw)), shape = 17, size = 2, fill = my_colours[5], color = my_colours[5]) +
  geom_hline(yintercept = 1.536057, color = my_colours[2], linetype = "31", linewidth=0.6) +
  geom_hline(yintercept = 1.551977, color = my_colours[2], linetype = "31", linewidth=0.6) +
  coord_cartesian(xlim=c(d_Z/2 -5,d_Z-1), ylim=c(0, 3)) +
  scale_x_continuous(expand = c(0, 0)) +
  # scale_x_continuous(limits = c(45, 100)) +
  labs(x = bquote(paste("Max # invalid instruments, ", italic(d_Z - b), " (", italic(d)[italic(Z)], " = ", .(d_Z),")")),
       y = expression(paste("Feasible  ", theta, " (95% confidence/credible set)"))) +
  annotate("text", x = 60, y = 1, label = expression(paste(theta,"*")),
           hjust = 0.3, vjust = -0.5, color = my_colours[1], size = 10) +
  annotate("text", x = 75, y = 2.8, label = expression(paste(b,"*")),
           hjust = 0.6, vjust = -0.5, color = my_colours[1], size = 10) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.line = element_line(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'bottom')

# benchmark_plot <- ggplot() +
#   geom_hline(yintercept = theta_true, color = "blue", linewidth=2) +
#   geom_vline(xintercept = d_Z-number_valid, color = "blue", linewidth=2) +
#   geom_segment(data = feasible_regions,
#                aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 5) +
#   geom_segment(data = oracle_regions, color = "blue",
#                aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 5) +
#   geom_hline(yintercept = benchmarks$mbe@CIUpper, color = "orange", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = benchmarks$mbe@CILower, color = "orange", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = benchmarks$median@CIUpper, color = "red", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = benchmarks$median@CILower, color = "red", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = benchmarks$mr_egger@CIUpper.Est, color = "purple", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = benchmarks$mr_egger@CILower.Est, color = "purple", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = benchmarks$ivw@CIUpper, color = "yellow", linetype = "31", linewidth=0.5) +
#   geom_hline(yintercept = benchmarks$ivw@CILower, color = "yellow", linetype = "31", linewidth=0.5) +
#   geom_hline(yintercept = 1.536057, color = "gray", linetype = "31", linewidth=1) +
#   geom_hline(yintercept = 1.551977, color = "gray", linetype = "31", linewidth=1) +
#   coord_cartesian(xlim=c(d_Z/2,d_Z-2), ylim=c(0, 3)) +
#   labs(x = bquote(paste("Max # invalid instruments, ", italic(b), " (", italic(d)[italic(Z)], " = ", .(d_Z),")")),
#        y = expression(paste("Feasible  ", theta, " (95% coverage)"))) +
#   annotate("text", x = 60, y = 1, label = expression(paste(theta,"*")),
#            hjust = 0.3, vjust = -0.5, color = "blue", size = 10) +
#   annotate("text", x = 75, y = 2.8, label = expression(paste(b,"*")),
#            hjust = 0.6, vjust = -0.5, color = "blue", size = 10) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 16),
#         axis.text = element_text(size = 16),
#         legend.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.position = 'bottom')

library(extrafont)

# Register fonts
loadfonts(device = "pdf")

benchmark_plot

ggsave(paste0('color_coded_high_d_Z_', d_Z, '_theta_true_1_', number_valid, '_valid_instruments_super_zoomed_in_labelled.pdf'), width = 6, height = 5.325,  device = cairo_pdf)

print(benchmarks)

plots <- list('posterior_plot'=posterior_plot, 'benchmark_plot'=benchmark_plot)

saveRDS(plots, 'plots.rds')


# devtools::install_github('jpenn2023/budgetivr')
# 
# library(budgetivr)
# 
# data(Do_et_al_summary_statistics)
# 
# candidatesHDL = Do_et_al_summary_statistics[Do_et_al_summary_statistics$pHDL <= 1e-8, ]
# 
# d_Z <- nrow(candidatesHDL)
# 
# SE_beta_y <- abs(candidatesHDL$betaCAD) / qnorm(1-candidatesHDL$pCAD/2)
# 
# # For 95% (asymptotic) confidence set.
# alpha <- 0.05
# 
# feasible_region <- budgetIV_scalar(
#                                    beta_y = candidatesHDL$betaCAD,
#                                    beta_phi = candidatesHDL$betaHDL,
#                                    tau_vec = c(0),
#                                    b_vec = c(30),
#                                    delta_beta_y = qnorm(1 - alpha/(2*d_Z))*SE_beta_y,
#                                    bounds_only = FALSE
#                                    )
# Experiment used to generate Fig. 5 in Penn et al. (2025).
#
# Jordan Penn, Lee Gunderson, Gecia Bravo-Hermsdorff,
# Ricardo Silva, and David Watson. (2024). BudgetIV: Optimal Partial Identification of Causal Effects with Mostly Invalid Instruments. \emph{arXiv}
# preprint, 2411.06913.

#
# Uncomment the following to install the packages used in this experiment:
# 
# devtools::install_github('jpenn2023/BudgetIV')
# 
# library(randcorr)
# library(MASS)
# library(ggplot2)
# library(ggsci)
# library(sisVIVE)
# library(MendelianRandomization)
# library(boot)
# library(data.table)
#

set.seed(42)

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
  
  return(list('summary_stats'=summary_stats, 'SS_massive'=SS_massive, 'empirical_corr_z'=corr_z_for_x))
  
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
    
    feasible_region <- BudgetIV_scalar_exposure_feature(beta_y, beta_x, tau_vec, b_vec, delta_beta_y)
    
    print(feasible_region)
    
    feasible_intervals <- feasible_region$intervals
    
    if(length(feasible_intervals) != 0){
      
      for (interval_idx in 1:(length(feasible_intervals)/2)){
        
        print(paste0("Next interval at b = ", b_vec[1], " is: ", interval_idx))
        
        next_res <- list("b"=b_vec[1], "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2])
        
        results <- rbind(results, next_res)
      }
    }
  }
  
  print(results)
  
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
    
    feasible_region <- BudgetIV_scalar_exposure_feature(beta_y, beta_x, tau_vec, b_vec, delta_beta_y)
    
    print(feasible_region)
    
    feasible_intervals <- feasible_region$intervals
    
    if(length(feasible_intervals) != 0){
      
      for (interval_idx in 1:(length(feasible_intervals)/2)){
        
        print(paste0("Next interval at b = ", b_vec[1], " is: ", interval_idx))
        
        next_res <- list("b"=b_vec[1], "feasible_start"=feasible_intervals[interval_idx, 1],
                         "feasible_end"=feasible_intervals[interval_idx, 2])
        
        results <- rbind(results, next_res)
      }
    }
  }
  
  print(results)
  
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

study_sizes <- c(20, 40, 60, 80, 100)
nums_valid <- 0.3 * study_sizes

studies <- 5

theta_true <- 1
N_x <- 1e6
N_y <- 1e5

d_Z <- study_sizes[study_idx]
number_valid <- nums_valid[study_idx]

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

benchmark_plot <- ggplot() +
  geom_hline(yintercept = theta_true, color = "blue", linewidth=2) +
  geom_vline(xintercept = d_Z-number_valid, color = "blue", linewidth=2) +
  geom_segment(data = feasible_regions,
               aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 5) +
  geom_segment(data = oracle_regions, color = "blue",
               aes(x=d_Z-b, xend=d_Z-b, y=feasible_start, yend=feasible_end), linewidth = 5) +
  geom_hline(yintercept = benchmarks$mbe@CIUpper, color = "orange", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mbe@CILower, color = "orange", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$median@CIUpper, color = "red", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$median@CILower, color = "red", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mr_egger@CIUpper.Est, color = "purple", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$mr_egger@CILower.Est, color = "purple", linetype = "31", linewidth=1) +
  geom_hline(yintercept = benchmarks$ivw@CIUpper, color = "yellow", linetype = "31", linewidth=0.5) +
  geom_hline(yintercept = benchmarks$ivw@CILower, color = "yellow", linetype = "31", linewidth=0.5) +
  geom_hline(yintercept = 1.536057, color = "gray", linetype = "31", linewidth=1) +
  geom_hline(yintercept = 1.551977, color = "gray", linetype = "31", linewidth=1) +
  coord_cartesian(xlim=c(d_Z/2,d_Z-2), ylim=c(0, 3)) +
  labs(x = bquote(paste("Max # invalid instruments, ", italic(b), " (", italic(d)[italic(Z)], " = ", .(d_Z),")")),
       y = expression(paste("Feasible  ", theta, " (95% coverage)"))) +
  annotate("text", x = 60, y = 1, label = expression(paste(theta,"*")),
           hjust = 0.3, vjust = -0.5, color = "blue", size = 10) +
  annotate("text", x = 75, y = 2.8, label = expression(paste(b,"*")),
           hjust = 0.6, vjust = -0.5, color = "blue", size = 10) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'bottom')

benchmark_plot

ggsave(paste0('color_coded_high_d_Z_', d_Z, '_theta_true_1_', number_valid, '_valid_instruments_super_zoomed_in_labelled.png'), width = 4)

print(benchmarks)

plots <- list('posterior_plot'=posterior_plot, 'benchmark_plot'=benchmark_plot)

saveRDS(plots, 'plots.rds')
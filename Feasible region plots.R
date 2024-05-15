setwd('C:/Users/k23067841/Downloads/BudgetIV')

library(data.table)

source('BudgetIV.R')

setwd('C:/Users/k23067841/Downloads/BudgetIV_data/real deal 2')

#
# Recall the ground truth: theta* = 1, A* = (-4, 1), B* = (-2, 0.9)
#
# => Cov(Z, g_y*) = (-2, 0.1)
#
#

# Ratio of the radius of the benchmark L_2 ball to the distance
# norm(tau_vec, type="2")
slack_scaling <- 1

number_runs <- 30

tau_1 <- 0.3

# p <- read.csv('1.csv')

# print(p$Z1_dat)

results <- data.table("setting"=numeric(), "sim_idx"=numeric(), "interval"=logical(), "feasible_start"=numeric(), "feasible_end"=numeric(), "tau_lo"=numeric(), "tau_hi"=numeric())
results_l2 <- data.table("setting"=numeric(), "sim_idx"=numeric(), "theta_lo"=numeric(), "theta_hi"=numeric(), 
                         "tau_lo"=numeric(), "tau_hi"=numeric(), "slack"=numeric)

for (setting in 1:number_runs){
# for (setting in 1:2){
  
  curr_data <- read.csv(paste0(setting,'.csv'))

  A <- c(cov(curr_data$Z1_dat, curr_data$Y_dat), cov(curr_data$Z2_dat, curr_data$Y_dat))
  B <- c(cov(curr_data$Z1_dat, curr_data$X_dat), cov(curr_data$Z2_dat, curr_data$X_dat))

  # Change tau_2 from 0.2 through to 6
  curr_tau_2 <- setting/5
  # 
  # print(paste0("Cov (Y, Z): ", A))
  # 
  # print(paste0("Cov (X, Z): ", B))
  # 
  # print(curr_tau_2)
  
  if (curr_tau_2 < tau_1) {
    tau_vec <- c(curr_tau_2, tau_1)
    m_vec <- c(1,2)
    }
  
  else if (curr_tau_2 > tau_1) {
    tau_vec <- c(tau_1, curr_tau_2)
    m_vec <- c(1,2)
    }
  
  else if (curr_tau_2 == tau_1) {
    tau_vec <- c(tau_1)
    m_vec <- c(2)
  }
  
  curr_feasible_int_l2 <- l2_ball(A, B, tau_vec, slack_scaling)
  
  if (length(curr_feasible_int_l2) != 0){
  
    results_l2 <- rbind(results_l2, list("setting"=setting, "sim_idx"=sim_idx, 
                                         "theta_lo"=curr_feasible_int_l2$theta_lo, "theta_hi"=curr_feasible_int_l2$theta_hi,
                                         "tau_lo"= tau_vec[1], "tau_hi"=tau_vec[2], "slack"=slack_scaling))
  
    }
  
  feasible_region <- BudgetIV_scalar_exposure_feature(A, B, tau_vec, m_vec)
  
  feasible_points <- feasible_region$points
  feasible_intervals <- feasible_region$intervals
  
  # print(feasible_intervals)


  if (!is.na(feasible_points)) {
    for (point_idx in 1:length(feasible_points)){
      
      next_res <- list("setting"=setting, "sim_idx"=sim_idx, "interval"=FALSE, "feasible_start"=feasible_points[point_idx], 
                       "feasible_end"=feasible_points[point_idx], "tau_lo"=tau_vec[1], "tau_hi"=tau_vec[2])
      
      results <- rbind(results, next_res)
  }}

  if (length(feasible_intervals) != 0) {
    # print(feasible_intervals)
    
    # print(length(feasible_intervals))
    
    for (interval_idx in 1:(length(feasible_intervals)/2)){
      
      print(interval_idx)
      
      
      next_res <- list("setting"=setting, "sim_idx"=1, "interval"=TRUE, "feasible_start"=feasible_intervals[interval_idx, 1], 
                       "feasible_end"=feasible_intervals[interval_idx, 2], "tau_lo"=tau_vec[1], "tau_hi"=tau_vec[2])
      
      # print(next_res)
      results <- rbind(results, next_res)
      
  }}
  
  # print(results)
}

# print(results)

fwrite(results, './results.csv')
fwrite(results_l2, './results_l2.csv')

# print(p)
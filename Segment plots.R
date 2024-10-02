#
# Barcode plots 1
#
#

library(data.table)
library(graphics)

#setwd('C:/Users/k23067841/Downloads/BudgetIV')

#source('Feasible region plots.R')

setwd('C:/Users/k23067841/Downloads/BudgetIV_experiments/small dphi experiments 2/lazy')

feasible_regions <- read.csv('results.csv')
feasible_intervals_l2 <- read.csv('results_l2.csv')

feasible_intervals_l2 <- na.omit(feasible_intervals_l2)

#print(feasible_regions$tau_lo)

#print(feasible_intervals_l2$theta_lo)

print(c(feasible_intervals_l2$theta_lo, rev(feasible_intervals_l2$theta_hi)))

# # Calculate min and max with a small slack
# slack <- 0.1  # Slack as a proportion of the range (10%)
# 
# x <- feasible_intervals_l2$tau_var
# y <- 
# 
# # Compute range for x and y
# x_range <- range(x)
# y_range <- range(y)
# 
# # Calculate slack values
# x_slack <- (x_range[2] - x_range[1]) * slack
# y_slack <- (y_range[2] - y_range[1]) * slack
# 
# # Set xlim and ylim with slack
# xlim <- c(x_range[1] - x_slack, x_range[2] + x_slack)
# ylim <- c(y_range[1] - y_slack, y_range[2] + y_slack)

plot(1, type="n", xlab="IV Violation parameter (tau)", ylab="Feasible regions for causal effect", xlim=c(0, 10), ylim=c(-1, 3))
#polygon(c(feasible_intervals_l2$tau_hi, rev(feasible_intervals_l2$tau_hi)), C(feasible_intervals_l2$theta_lo, rev(feasible_intervals_l2$theta_hi)))

abline(1, 0, col="red", lty =2)

segments(feasible_intervals_l2$var_tau, feasible_intervals_l2$theta_lo, y1 = feasible_intervals_l2$theta_hi, col="lightblue", lwd = 4, lend = "butt")

segments(feasible_regions$var_tau, feasible_regions$feasible_start, y1 = feasible_regions$feasible_end, lwd = 3, col="orange", lend = "butt")


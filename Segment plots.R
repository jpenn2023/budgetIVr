#
# Barcode plots 1
#
#

library(data.table)
library(graphics)

setwd('C:/Users/k23067841/Downloads/BudgetIV')

#source('Feasible region plots.R')

setwd('C:/Users/k23067841/Downloads/BudgetIV_data/real deal 2')

feasible_regions <- read.csv('results.csv')
feasible_intervals_l2 <- read.csv('results_l2.csv')

feasible_intervals_l2 <- na.omit(feasible_intervals_l2)

#print(feasible_regions$tau_lo)

#print(feasible_intervals_l2$theta_lo)

print(c(feasible_intervals_l2$theta_lo, rev(feasible_intervals_l2$theta_hi)))

plot(1, type="n", xlab="tau_2", ylab="Feasible theta", xlim=c(0, 6), ylim=c(-1.5, 5))

#polygon(c(feasible_intervals_l2$tau_hi, rev(feasible_intervals_l2$tau_hi)), C(feasible_intervals_l2$theta_lo, rev(feasible_intervals_l2$theta_hi)))

abline(1, 0, col="red", lty =2)

segments(feasible_intervals_l2$tau_hi, feasible_intervals_l2$theta_lo, y1 = feasible_regions$feasible_end, col="blue", lwd = 4, lend = "butt")

segments(feasible_regions$tau_hi, feasible_regions$feasible_start, y1 = feasible_regions$feasible_end, lwd = 3, lend = "butt")


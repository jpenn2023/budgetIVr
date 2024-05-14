setwd('C:/Users/k23067841/Downloads/BudgetIV')

source('BudgetIV.R')

setwd('C:/Users/k23067841/Downloads/BudgetIV_data/real deal')

#
# Recall the ground truth: theta* = 1, A* = (-4, 1), B* = (-2, 0.9)
#
# => Cov(Z, g_y*) = (-2, 0.1)
#
#

number_runs <- 30

p <- read.csv('1.csv')

# print(p$Z1_dat)

points <- list()
intervals <- list()


for (setting in 1:number_runs){
  
  p <- read.csv(paste0(setting,'.csv'))
  print(setting)
  
  
  
}

print(p)
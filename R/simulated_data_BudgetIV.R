#' Simulated summary statistics with invalid instruments and nonlinear treatment effect
#'
#' Example dataset from the nonlinear simulation study using 6 candidate instruments, 3 of which are invalid with violation
#' of IV assumptions (A2) and (A3). 
#' See Appx. C.2 of Penn et al. (2025) for technical details or visit the source code for reproducibility, both referenced below.
#' The ground truth causal effect is \eqn{\Phi^* (X) = (X - 0.25)^2 - 0.25^2}.
#' \eqn{\beta_{\Phi}} is taken with respect to the basis functions \eqn{\Phi = (X, X^2)}.
#' 
#' @docType data
#'
#' @usage data(simulated_data_budgetIV)
#'
#' @format A data frame with 6 rows and 4 columns.
#' 
#' @details
#' \describe{
#'  \item{\code{beta_y}}{Components of the estimator \eqn{\mathrm{Cov} (Y, Z)}.}
#'  \item{\code{beta_phi_1}}{Components of the estimator \eqn{\mathrm{Cov} ( \Phi_1 (X), Z )}.}
#'  \item{\code{beta_phi_2}}{Components of the estimator \eqn{\mathrm{Cov} ( \Phi_2 (X), Z )}.}
#'  \item{\code{delta_beta_y}}{Components of the standard error \eqn{\mathrm{Se} (\mathrm{Cov} (Y, Z))}.}
#' }
#'
#' @keywords datasets
#'
#' @references Jordan Penn, Lee Gunderson, Gecia Bravo-Hermsdorff,
#' Ricardo Silva, and David Watson. (2024). BudgetIV: Optimal Partial Identification of Causal Effects with Mostly Invalid Instruments. \emph{arXiv}
#' preprint, 2411.06913.
#' 
#' @source The code that generated this dataset was written by the authors and can be found in \url{https://github.com/jpenn2023/budgetivr/tree/main/paper/simulate_nonlinear_data.R}.
#' The dataset is saved as "my_dat R = 0.5 SNR_y = 1.csv".
#' 
#' @examples
#' data(simulated_data_budgetIV)
#'
#' beta_y <- simulated_data_budgetIV$beta_y
#' 
#' beta_phi_1 <- simulated_data_budgetIV$beta_phi_1
#' beta_phi_2 <- simulated_data_budgetIV$beta_phi_2
#' 
#' d_Z <- length(beta_phi_1)
#' 
#' beta_phi <- matrix(c(beta_phi_1, beta_phi_2), nrow = 2, byrow = TRUE)
#' 
#' delta_beta_y <- simulated_data_budgetIV$delta_beta_y
"simulated_data_budgetIV"
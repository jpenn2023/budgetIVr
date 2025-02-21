#' BudgetIV for scalar exposures
#' 
#' Partial identification and coverage of a causal effect parameter using summary statistics and budget constraint assumptions.
#' 
#' @param beta_y A \eqn{1 \times d_{Z}} matrix representing the (estimated) cross covariance \eqn{\mathrm{Cov}(Y, Z)}.
#' @param beta_phi A \eqn{d_{\Phi} \times d_{Z}} matrix representing the (estimated) cross covariance \eqn{\mathrm{Cov}(\Phi (X), Z)}.
#' @param phi_basis A 
#' @param tau_vec, A \eqn{K}-dimensional vector of strictly increasing, positive budget thresholds. \eqn{K} is the number of budget groups.
#' @param b_vec A \eqn{K}-dimensional vector of increasing positive integer budgets. 
#' Represents the constraint that at least \eqn{b_i} different values of \eqn{j}, 
#' the candidate instrument \eqn{Z_j} satisfies \eqn{\mathrm{Cov} (Y - \theta \Phi (X), Z_j) \leq \tau_j}.
#' @param delta_beta_y Either (a) a \eqn{d_{Z}}-dimensional vector of positive half-widths for box-shaped confidence bounds on \code{beta_y};
#' or (b) the empty value \code{NA} for partial identification or feasible region estimates without uncertainty quantification. 
#' @param ATE_search_domain A \eqn{d_{X} \times N} matrix consisting of all \eqn{x \in \mathbb{R}^{d_{X}}} for which we want to calculate 
#' \eqn{ATE(x, x_0)}
#' @param X_baseline The baseline treatment \eqn{x_0}. Either (a) a \eqn{d_{X}}-dimensional vector; or (b) \code{NA} (default), in which case 
#' \eqn{x_0} is set to a vector of zeros. 
#' @param dummy_infinity Dummy value for \eqn{\tau_{K+1} := + \infty}. Necessary for linear programming approach and set to \code{1e10} by default. 
#' 
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: (A1) they are associated with the treatment; 
#' (A2) they are unconfounded with the outcome; and (A3) they exclusively effect the outcome through the treatment. 
#' Assumption (A1) has a simple statistical test, whereas for many data generating processes (A2) and (A3) are 
#' unprovably false.
#' The \code{BudgetIV} and \code{BudgetIV_scalar} algorithms allow for valid causal inference when some proportion, 
#' possibly a small minority, of candidate instruments satisfy both (A2) and (A3).
#' Tuneable thresholds decided by the user also allow for bounds on the degree of invalidity for each instrument 
#' (i.e., bounds on the proportion of \eqn{\mathrm{Cov}(Y, Z)} not explained by the causal effect of \eqn{X} on \eqn{Z}).  
#' 
#' \code{BudgetIV} assumes a homogeneous treatment effect, which implies the separable structural 
#' equation \eqn{Y = \theta \Phi(X) + g_y(Z, \epsilon_x)}, where \eqn{\theta} and \eqn{\Phi(X)} are a 
#' \eqn{d_{\Phi}} dimensional vector and vector-valued function respectively. 
#' A valid basis expansion \eqn{\Phi (X)} is assumed, e.g., linear, logistic, polynomial, RBF, hazard model, neural network.
#' It is also assumed that \eqn{d_{\Phi} < d_{Z}}, which allows us to treat the basis functions as a complete linear model
#' (see Theil (1953)).
#' The parameters \eqn{\theta} captures the unknown treatment effect.
#' Violation of (A2) and/or (A3) will bias classical IV approaches through the statistical dependence
#' between \eqn{Z} and \eqn{g_y(Z, \epsilon_x)}, summarized by the covariance parameter 
#' \eqn{\gamma := \mathrm{Cov} (g_y(Z, \epsilon_x), Z)}.
#' 
#' \code{BudgetIV} constrains \eqn{\gamma} through a series of positive thresholds 
#' \eqn{0 \leq \tau_1 < \tau_2 < \ldots < \tau_K} and corresponding integer budgets \eqn{0 < b_1 < b_2 < \ldots < b_K \leq d_Z}. 
#' It is assumed for each \eqn{i \in \{ 1, \ldots, K\}} that no more than \eqn{b_i} components of \eqn{\gamma} are greater in 
#' magnitude than \eqn{\tau_i}.
#' For instance, taking \eqn{d_Z = 100}, \eqn{K = 1}, \eqn{b_1 = 5} and \eqn{\tau_1 = 0} means 
#' assuming \eqn{5\%} of the \eqn{100} candidates are valid instrumental variables (in the sense that their ratio 
#' estimates \eqn{\theta_j := \mathrm{Cov}(Y, Z_j)/\mathrm{Cov}(\Phi(X), Z_j)} are unbiased).
#' 
#' With \code{delta_beta_y = NA}, \code{BudgetIV} & \code{BudgetIV_scalar} return the identified set
#' of causal effects that agree with both the budget constraints described above and the values of
#' \eqn{\mathrm{Cov}(Y, Z)} and \eqn{\mathrm{Cov}(Y, Z)}, assumed to be exactly precise. 
#' Unlike classical partial identification methods (see Manski (1990) for a canonical example), the non-convex mixed-integer
#' budget constraints yield a possibly disconnected identified set. 
#' Each connected subset has a different interpretation as to which of the candidate instruments \eqn{Z} 
#' are valid up to each threshold.
#' \code{BudgetIV_scalar} returns these interpretations alongside the corresponding bounds on \eqn{\theta}. 
#' 
#' When \code{delta_beta_y} is not null, it is used as box-constraints to quantify uncertainty in \code{beta_y}. 
#' In the examples, \code{delta_beta_y} is calculated through a Bonferroni correction and gives an (asymptotically) 
#' valid confidence set over \code{beta_y}. 
#' Under the so-called "no measurement error" assumption (see Bowden et al. (2016)) which is commonly applied in Mendelian randomisation, it is
#' assumed that the estimate of \code{beta_y} is the dominant source of finite-sample uncertainty, with uncertainty in \code{beta_x}
#' entirely negligible. 
#' With an (asymptotically) valid confidence set for \code{delta_beta_y} and under the "no measurement error" assumption, \code{BudgetIV_scalar} 
#' returns an (asymptotically) valid confidence set for \eqn{\theta}.  
#' 
#' @return  
#' A dataframe consisting of .
#' 
#' @references  
#' Jordan Penn, Lee Gunderson, Gecia Bravo-Hermsdorff,
#' Ricardo Silva, and David Watson. (2024). BudgetIV: Optimal Partial Identification of Causal Effects with Mostly Invalid Instruments. \emph{arXiv}
#' preprint, 2411.06913.
#' 
#' Jack Bowden, Fabiola Del Greco M, Cosetta Minelli, George Davey Smith, Nuala A Sheehan, and John R Thompson. (2016). Assessing the suitability of summary data for 
#' two-sample Mendelian randomization analyses using MR-Egger regression: the role of the I^2 statistic. \emph{Int. J. Epidemiol.} 46.6, pp. 1985--1998.
#' 
#' Charles F Manski. (1990). Nonparametric bounds on treatment effects. \emph{Am. Econ. Rev.} 80.2, pp. 219--323.
#' 
#' Henri Theil. (1953). Repeated least-squares applied to complete equation systems. \emph{Centraal Planbureau Memorandum}.
#' 
#' @examples  
#' set.seed(123)
#' 
#' # Simple experiment with multidimensional exposure.
#' 
#' 
#' 
#' 
#' # Mock summary statistics
#' 
#' 
#' 
#' beta_y = matrix(c(1,2,3,4), nrow=1, ncol=4)
#' beta_phi = matrix(c(4, 3, 2, 1, -3, -1, 2, 1), nrow=2, ncol=4)
#' 
#' # 
#' 
#' 
#' 
#' @export 
#' @import data.table
#' @import arrangements
#' @import MASS
#' @import Rglpk

BudgetIV <- function(
    beta_y, 
    beta_phi, 
    phi_basis, 
    tau_vec, 
    b_vec, 
    ATE_search_domain, 
    X_baseline=0,
    delta_beta_y=NA,
    tol=1e-10,
    dummy_infinity=1e10
) {
  
  if(is.vector(beta_y)){beta_y <- matrix(beta_y, nrow=1)}
  
  d_X <- ncol(ATE_search_domain)
  d_Z <- ncol(beta_y)
  
  if (all(is.na(delta_beta_y))){delta_beta_y <- numeric(d_Z)}
  
  if(is.vector(delta_beta_y)){delta_beta_y <- matrix(delta_beta_y, nrow=1)}
    
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
  else if (ncol(beta_y) != ncol(beta_phi)) {
    stop('Cov(Y, Z) and Cov(Phi(X), Z) must have the same number of columns.')
  }
  else if (ncol(beta_y) < nrow(beta_phi)){
    stop('BudgetIV only supports partial identification in the regime d_{Phi} <= d_{Z}')
  }
  else if (ncol(beta_y) != ncol(delta_beta_y)){
    stop('If specifying half-width errors delta_beta_y, there must be as many errors as components of beta_y. Run get_covariance with confidence_threshold set to sum numeric in (0,1) to compute these')
  }
  else if (!all(delta_beta_y >= 0)){
    stop("Please ensure half-width errors are greater than or equal to zero")
  }
  else if (length(phi_basis) != nrow(beta_phi)){
    stop("Please ensure the number of basis features phi_i equals the number of columns (beta_phi)_i = cov(phi_i, Z)")
  }
  
  d_Z <- ncol(beta_y)
  d_Phi <- nrow(beta_phi)
  
  # If there are as many constraints as instruments, don't include tau_{K+1} = infinity
  if(b_vec[length(b_vec)] == d_Z){
    
    # (tau_1, ..., tau_K)
    taus <- tau_vec
    
    # Vector of differences (m_1, m_2 - m_1, ..., m_K - m_{K-1})
    b_deltas <- c(b_vec[1], diff(b_vec))
  }
  
  # Otherwise, include tau_{K+1} = infinity (problem is "under-constrained")
  else{
    
    # (tau_1, ..., tau_K, tau_{K+1})
    taus <- c(tau_vec, dummy_infinity)
    
    # Vector of differences (m_1, m_2 - m_1, ..., m_K - m_{K-1}, d_Z - m_K)
    b_deltas <- c(b_vec[1], diff(b_vec), d_Z - b_vec[length(b_vec)])
    
    }
  
  # List to contain all dose reponse ATE bounds and the corresponding decision variable U
  partial_identification_ATE <- data.table(
    "curve_index" = numeric(),
    "x" = list(),
    "lower_ATE_bound" = numeric(),
    "upper_ATE_bound" = numeric(),
    "U" = list()
  )
  
  curve_index <- 0
  
  phi_zero <- lapply(phi_basis, function(e) eval(e, X_baseline[1, ]))
  
  # Iterate through the values of the one-hot encoding U. The iterator maps tau_i to i, so we have to reverse this map. 
  U_perm_iter <- ipermutations(taus, freq=b_deltas) 
  
  curr_U_perm <- U_perm_iter$getnext()
  
  # For every possible perm for every unique way to satisfy the budget constraints. 
  # There are d_Z! /( b_1!(b_2 - b_1)!...(b_{K_red}-b_{K_red - 1})! ) perms. 
  while (!is.null(curr_U_perm)) {
    
    bounds <- taus[curr_U_perm] + delta_beta_y
    
    f.con <- rbind(t(beta_phi), t(beta_phi))  # Combine the upper and lower bound constraints
    f.dir <- c(rep("<=", d_Z), rep(">=", d_Z))  # Directions of the inequalities
    f.rhs <- c(beta_y + bounds, beta_y - bounds)  # Right-hand side for the inequalities
    
    f.obj <- rep(0, d_Phi)  # Objective function for theta (length d_Phi, the dimension of theta)
    
    search_bounds <- list(lower = list(ind = c(1L, 2L), val = c(-Inf, -Inf)),
                          upper = list(ind = c(1L, 2L), val = c(Inf, Inf)))

    constraint_satisfaction <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = TRUE, bounds = search_bounds)
    
    
    if (constraint_satisfaction$status == 0){
      
      curve_index <- curve_index + 1
      
      for (coord_index in 1:nrow(ATE_search_domain)) {
        
        phi_curr_coord <- lapply(phi_basis, function(e) eval(e, ATE_search_domain[coord_index, ]))
        
        f.obj <- unlist(phi_curr_coord) - unlist(phi_zero)
        
        ATE_min_curr_coord <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = FALSE, bounds = search_bounds)
        ATE_max_curr_coord <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs, max = TRUE, bounds = search_bounds)
        
        new_bound <- data.table(
          "curve_index" = curve_index,
          "x" = list(ATE_search_domain[coord_index, ]),
          "lower_ATE_bound" = ATE_min_curr_coord$optimum,
          "upper_ATE_bound" = ATE_max_curr_coord$optimum,
          "U" = list(curr_U_perm)
        )
        
        partial_identification_ATE <- rbind(partial_identification_ATE, new_bound)
         
      }
      
      curve_index <- curve_index + 1
      
    }
    
    
    curr_U_perm <- U_perm_iter$getnext()
    
  }
  
  if(is.data.frame(ATE_bounds_example$x[[1]])){
    
    partial_identification_ATE <- partial_identification_ATE[, cbind(.SD, rbindlist(x, use.names = TRUE, fill = TRUE)), .SDcols = !'x']
  
  }
  
  return(partial_identification_ATE)
  
}

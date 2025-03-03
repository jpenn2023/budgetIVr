#' BudgetIV: partial identification of causal effects with invalid instruments
#' 
#' Partial identification and coverage of a causal effect parameter using summary statistics and budget constraint assumptions.
#' 
#' @param beta_y Either \eqn{1 \times d_{Z}} matrix or a \eqn{d_{Z}}-dimensional vector representing the (estimated) cross covariance \eqn{\mathrm{Cov}(Y, Z)}.
#' @param beta_phi A \eqn{d_{\Phi} \times d_{Z}} matrix representing the (estimated) cross covariance \eqn{\mathrm{Cov}(\Phi (X), Z)}.
#' @param phi_basis A \eqn{d_{\Phi}}-dimensional expression (separated by commas) with each term representing a component of \eqn{\Phi (X)}.
#' The expression consists of \eqn{d_{X}} unique vars. 
#' @param tau_vec, A \eqn{K}-dimensional vector of strictly increasing, positive thresholds (representing degrees if IV invalidity). 
#' @param b_vec A \eqn{K}-dimensional vector of increasing positive integer budgets. 
#' @param delta_beta_y Either (a) a \eqn{d_{Z}}-dimensional vector of positive half-widths for box-shaped confidence bounds on \code{beta_y};
#' or (b) the value \code{NA} for partial identification or feasible region estimates without uncertainty quantification. 
#' @param ATE_search_domain A \eqn{d_{X}}-column data.table or data.frame with column names equal to the vars in \code{phi_basis}.
#' Rows correspond to values of the treatment \eqn{X}. 
#' @param X_baseline Either a data.table, data.frame or list. 
#' If entered as a data.table or data.frame, \code{X_baseline} must have \eqn{d_{X}} columns and one row, with column names equal to the
#' vars in \code{phi_basis}. 
#' If entered as a list, entries must be numeric or one-dimensional vectors, with entry names equal to the vars in \code{phi_basis}. 
#' Correspond to the baseline treatment \eqn{x_0}.
#' @param dummy_infinity Dummy value for \eqn{\tau_{K+1} := + \infty}. 
#' Necessary for linear programming approach and set to \code{1e10} by default. 
#' Increase if \code{beta_y}, \code{beta_Phi} have entries greater than or equal to \code{1e8}.
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
#' A data.table with each row corresponding to a set of bounds on the ATE at a given point in \code{ATE_search_domain}. 
#' Columns include: a non-unique identifier \code{curve_index} with a one-to-one mapping with \code{U}; \code{lower_ATE_bound} and \code{upper_ATE_bound} for the corresponding
#' bounds on the ATE; a list \code{U} for the corresponding budget assignment; and a column for each unique variable in \code{ATE_search_domain} to indicate the treatment value at 
#' which the bounds are being calculated.   
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
#' data(simulated_data_BudgetIV)
#'
#' beta_y <- simulated_data_BudgetIV$beta_y
#' 
#' beta_phi_1 <- simulated_data_BudgetIV$beta_phi_1
#' beta_phi_2 <- simulated_data_BudgetIV$beta_phi_2
#' 
#' beta_phi <- matrix(c(beta_phi_1, beta_phi_2), nrow = 2, byrow = TRUE)
#' 
#' delta_beta_y <- simulated_data_BudgetIV$delta_beta_y
#' 
#' tau_vec = c(0)
#' b_vec = c(3)
#' 
#' x_vals <- seq(from = 0, to = 1, length.out = 500)
#' 
#' ATE_search_domain <- expand.grid("x" = x_vals)
#' 
#' phi_basis <- expression(x, x^2)
#' 
#' X_baseline <- list("x" = c(0))
#' 
#' partial_identification_ATE <- BudgetIV(beta_y = beta_y, 
#'                                        beta_phi = beta_phi, 
#'                                        phi_basis = phi_basis, 
#'                                        tau_vec = tau_vec, 
#'                                        b_vec = b_vec, 
#'                                        ATE_search_domain = ATE_search_domain, 
#'                                        X_baseline = X_baseline,
#'                                        delta_beta_y = delta_beta_y)
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
    X_baseline,
    delta_beta_y=NA,
    dummy_infinity=1e10
) {
  
  if(is.vector(beta_y) && is.numeric(beta_y)){beta_y <- matrix(beta_y, nrow=1)}
  
  d_X <- ncol(ATE_search_domain)
  d_Z <- ncol(beta_y)
  
  if (all(is.na(delta_beta_y))){delta_beta_y <- numeric(d_Z)}
  
  if(is.vector(delta_beta_y)){delta_beta_y <- matrix(delta_beta_y, nrow=1)}
    
  # Error messages
  if(!is.matrix(beta_y)){
    stop("Argument 'beta_y' must be a vector or single-row matrix.")
  }
  else if(!is.numeric(beta_y)){
    stop("Argument 'beta_y' must have numeric entries.")
  }
  else if(nrow(beta_y) != 1){
    stop("Argument 'beta_y', if input as a matrix, must only have one row.")
  }
  else if(!is.matrix(beta_phi)){
    stop("Argument 'beta_phi' must be a matrix.")
  }
  else if(!is.numeric(beta_phi)){
    stop("Argument 'beta_phi' must have numeric entries.")
  }
  else if (ncol(beta_y) != ncol(beta_phi)) {
    stop("Arguments 'beta_y' and 'beta_phi' must have the same number of columns.")
  }
  else if (ncol(beta_phi) < nrow(beta_phi)){
    stop("Argument 'beta_phi' must have more columns than rows. BudgetIV only supports partial identification in the 'complete' in which the number of causal effect parameters 
         is no greater than the number of candidate instruments (d_{Phi} <= d_{Z}). See the package documentation or Penn et al. (2025) for further details.")
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
  
  if(is.data.frame(ATE_search_domain)){
    ATE_search_domain <- as.data.table(ATE_search_domain)
  }
  
  else if (is.list(ATE_search_domain)){
    if (all(sapply(ATE_search_domain, is.vector))){
      if (length(unique(sapply(lst, length))) == 1){ 
        ATE_search_domain <- as.data.table(ATE_search_domain)
      }
      else{
        stop("Argument 'ATE_search_domain', if entered as a list, must be compatible with a data.frame structure.")
      }
    }
    else{
      stop("Argument 'ATE_search_domain', if entered as a list, must consist of named vector entries.")
    }
  }
  
  if (is.list(X_baseline) || is.data.frame(X_baseline)){
    X_baseline <- as.data.table(X_baseline)
  }
  
  if(!is.data.table(ATE_search_domain) && !is.data.frame(ATE_search_domain)){
    stop("Argument 'ATE_search_domain' must be a data.table, a data.frame or a list with names corresponding to the variable names in 'phi_basis'.")
  }
  else if(!is.data.table(X_baseline) && !is.data.frame(X_baseline) && !is.list(ATE_search_domain)){
    stop("Argument 'X_baseline' must be a data.table or data.frame with names corresponding to the variable names in 'phi_basis'.")
  }
  
  if (any(is.na(delta_beta_y)) ){
    warning("No confidence bounds for agument 'beta_y' given: treating 'beta_y' as an oracle summary statistic.")
    delta_beta_y <- numeric(d_Z)
  }
  
  if (ncol(delta_beta_y) != ncol(beta_y)){
    stop("Argument 'delta_beta_y', if given, must be of the same length as beta_y.")
  }
  else if(!is.numeric(delta_beta_y)){
    stop("Argument 'delta_beta_y' must have numeric entries.")
  }
  else if(nrow(beta_y) != 1){
    stop("Argument 'delta_beta_y', if input as a matrix, must only have one row.")
  }
  else if(!is.numeric(delta_beta_y)){
    stop("Argument 'delta_beta_y' must have numeric entries.")
  }
  else if(any(delta_beta_y < 0)){
    stop("Argument 'delta_beta_y' must have positive entries.")
  }
  else if(is.list(phi_basis)){
    stop("Please input argument 'phi_basis' as a single expression (e.g., 'expression(x, x - y)') rather than a list (e.g., 'list(expression(x), expression(x - y))').")
  }
  else if(!is.expression(phi_basis)){
    stop("Argument 'phi_basis' must be an expression.")
  }
  
  if(ncol(ATE_search_domain) != ncol(X_baseline)){
    stop("Arguments 'ATE_search_domain' and 'X_baseline' must have the same number of columns.")
  }
  else if(any(all.vars(phi_basis) != names(ATE_search_domain))){
    if(length(all.vars(phi_basis)) > ncol(ATE_search_domain)){
      stop("The column names of argument 'ATE_search_domain' must match the variables in phi_basis. There are more variables in 'phi_basis' than in 'ATE_search_domain'.")
    }
    else if(length(all.vars(phi_basis)) < ncol(ATE_search_domain)){
      stop("The column names of argument 'ATE_search_domain' must match the variables in phi_basis. There are fewer variables in 'phi_basis' than in 'ATE_search_domain'.")
    }
    else{
      stop("The column names of argument 'ATE_search_domain' must match the variables in phi_basis.")
    }
  }
  else if(any(names(X_baseline) != names(ATE_search_domain))){
    stop("The column names of arguments 'ATE_search_domain' and 'X_baseline' must match.")
  }
  else if (length(phi_basis) != nrow(beta_phi)){
    stop("The length of the argument 'phi_basis' must be equal to 'nrow(beta_phi)'. Each row of 'beta_phi' should correspond to a unique basis function.")
  }
  else if (!all(sapply(ATE_search_domain, is.numeric))){
    stop("Argument 'ATE_search_domain' must contain only numeric data.")
  }
  else if (!all(sapply(X_baseline, is.numeric))){
    stop("Argument 'X_baseline' must contain only numeric data.")
  }
  
  if( nrow(beta_phi) == 1){
    warning("Since 'nrow(beta_phi) = 1', consider using BudgetIV_scalar. BudgetIV_scalar can partially identify the scalar 
            causal effect parameter of interest with superexponential improvement on time complexity.")
  }
  
  d_Z <- ncol(beta_y)
  d_Phi <- nrow(beta_phi)
  
  # If there are as many constraints as instruments, don't include tau_{K+1} = infinity
  if(b_vec[length(b_vec)] == d_Z){
    
    # (tau_1, ..., tau_K)
    taus <- tau_vec
    
    # Vector of differences (b_1, b_2 - b_1, ..., b_K - b_{K-1})
    b_deltas <- c(b_vec[1], diff(b_vec))
  }
  
  # Otherwise, include tau_{K+1} = infinity (problem is "under-constrained")
  else{
    
    # (tau_1, ..., tau_K, tau_{K+1})
    taus <- c(tau_vec, dummy_infinity)
    
    # Vector of differences (b_1, b_2 - b_1, ..., b_K - b_{K-1}, d_Z - b_K)
    b_deltas <- c(b_vec[1], diff(b_vec), d_Z - b_vec[length(b_vec)])
    
    }
  
  # List to contain bounded ATE curves and the corresponding decision variable U
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
      
    }
    
    curr_U_perm <- U_perm_iter$getnext()
    
  }
  
  if(nrow(partial_identification_ATE) == 0){
    
    warning("Identified/confidence set is empty in this case - the background assumptions are too strict.")
    return(partial_identification_ATE)
    
  }
  
  else{
  
  if(is.data.frame(partial_identification_ATE$x[[1]]) || is.data.table(partial_identification_ATE$x[[1]])){
    
    partial_identification_ATE <- partial_identification_ATE[, cbind(.SD, rbindlist(x, use.names = TRUE, fill = TRUE)), .SDcols = !'x']
  
  }
  
  return(partial_identification_ATE)
  
  }
}
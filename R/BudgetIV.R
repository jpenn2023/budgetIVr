#' Partially identify causal effects with invalid instruments
#' 
#' Computes the set of possible values of a causal parameter consistent with
#' observational data and given budget constraints.
#' See Penn et al. (2025) for technical definitions.
#' 
#' @param beta_y Either \eqn{1 \times d_{Z}} matrix or a \eqn{d_{Z}}-dimensional 
#' vector representing the (estimated) cross covariance \eqn{\mathrm{Cov}(Y, Z)}.
#' @param beta_phi A \eqn{d_{\Phi} \times d_{Z}} matrix representing the 
#' (estimated) cross covariance \eqn{\mathrm{Cov}(\Phi (X), Z)}.
#' @param phi_basis A \eqn{d_{\Phi}}-dimensional expression (separated by 
#' commas) with each term representing a component of \eqn{\Phi (X)}.
#' The expression consists of \eqn{d_{X}} unique vars. 
#' The default value \code{NULL} can be used for a \eqn{d_{X} = d_{\Phi}}-dimensional linear model. 
#' @param tau_vec A \eqn{K}-dimensional vector of increasing, positive 
#' thresholds representing degrees of IV invalidity. 
#' The default value \code{NULL} can be used for a single threshold at \eqn{0}.
#' @param b_vec A \eqn{K}-dimensional vector of increasing positive integers 
#' representing the maximum number of IVs that can surpass each threshold. 
#' The default value \code{NULL} can be used for a single threshold at \eqn{0}, with at least \eqn{50\%} of IVs assumed to be valid.
#' @param delta_beta_y A \eqn{d_{Z}}-dimensional vector of positive half-widths for box-shaped 
#' confidence bounds on \code{beta_y}. 
#' The default value \code{NULL} can be used to not include finite sample uncertainty.
#' @param ATE_search_domain A \eqn{d_{X}}-column data.frame with column names 
#' equal to the vars in \code{phi_basis}.
#' Rows correspond to values of the treatment \eqn{X}. 
#' The default value \code{NULL} can be used to generate a small \eqn{d_{X}}-dimensional grid.
#' @param X_baseline Either a data.frame or list representing a baseline 
#' treatment \eqn{x_0}, with names equal to the vars in \code{phi_basis}.
#' The default value \code{NULL} can be used for the baseline treatment \eqn{0} for each of of the \eqn{d_{X}} vars.
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: (A1) they are associated with the treatment; 
#' (A2) they are unconfounded with the outcome; and (A3) exclusively effect the 
#' outcome through the treatment. 
#' Of these, only (A1) can be tested without further assumptions. 
#' The \code{budgetIV} function allows for valid causal inference when some  
#' proportion (possibly a small minority) of candidate instruments satisfy 
#' both (A2) and (A3). Tuneable thresholds decided by the user also allow for 
#' bounds on the degree of invalidity for each instrument (i.e., bounds on the 
#' proportion of \eqn{\mathrm{Cov}(Y, Z)} not explained by the causal effect of 
#' \eqn{X} on \eqn{Z}). Full technical details are included in Penn et al. (2025).
#' 
#' \code{budgetIV} assumes that treatment effects are homogeneous, which implies 
#' a structural equation of the form \eqn{Y = \theta \cdot \Phi(X) + g_y(Z, \epsilon_x)}, 
#' where \eqn{\theta} and \eqn{\Phi(X)} are a \eqn{d_{\Phi}}-dimensional vector 
#' and vector-valued function respectively. A valid basis expansion \eqn{\Phi (X)} 
#' is assumed (e.g., linear, logistic, polynomial, RBF, neural embedding, PCA, UMAP etc.). 
#' It is also assumed that \eqn{d_{\Phi} <= d_{Z}}, which allows us to 
#' treat the basis functions as a complete linear model (see Theil (1953), or Sanderson et al. (2019) 
#' for a modern MR focused discussion).
#' The parameters \eqn{\theta} capture the unknown treatment effect. 
#' Violation of (A2) and/or (A3) will bias classical IV approaches through the statistical 
#' dependence between \eqn{Z} and \eqn{g_y(Z, \epsilon_x)}, summarized by the 
#' covariance parameter \eqn{\gamma := \mathrm{Cov} (g_y(Z, \epsilon_x), Z)}.
#' 
#' \code{budgetIV} constrains \eqn{\gamma} through a series of positive 
#' thresholds \eqn{0 \leq \tau_1 < \tau_2 < \ldots < \tau_K} and corresponding 
#' integer budgets \eqn{0 < b_1 < b_2 < \ldots < b_K \leq d_Z}. It is assumed 
#' for each \eqn{i \in \{ 1, \ldots, K\}} that no more than \eqn{b_i} components 
#' of \eqn{\gamma} are greater in magnitude than \eqn{\tau_i}. For instance, 
#' taking \eqn{d_Z = 100}, \eqn{K = 1}, \eqn{b_1 = 5} and \eqn{\tau_1 = 0} means 
#' assuming \eqn{5} of the \eqn{100} candidates are valid instrumental 
#' variables (in the sense that their ratio estimates \eqn{\theta_j := 
#' \mathrm{Cov}(Y, Z_j)/\mathrm{Cov}(\Phi(X), Z_j)} are unbiased).
#' 
#' With \code{delta_beta_y = NULL}, \code{budgetIV} returns the identified set 
#' of causal effects that agree with both the budget constraints described above 
#' and the values of \eqn{\mathrm{Cov}(Y, Z)} and \eqn{\mathrm{Cov}(Y, Z)}, 
#' assumed to be exactly precise. Unlike classical partial identification 
#' methods (see Manski (1990) for a canonical example), the non-convex 
#' mixed-integer budget constraints yield a possibly disconnected solution set. 
#' Each connected subset has a different interpretation as to which of the 
#' candidate instruments \eqn{Z} are valid up to each threshold. 
#' 
#' \code{delta_beta_y} represents box-constraints to 
#' quantify uncertainty in \code{beta_y}. In the examples, \code{delta_beta_y} 
#' is calculated through a Bonferroni correction and gives an (asymptotically) 
#' valid confidence set over \code{beta_y}. Under the so-called "no measurement 
#' error" assumption (see Bowden et al. (2016)), which is commonly applied in 
#' Mendelian randomization, it is assumed that the estimate of \code{beta_y} is 
#' the dominant source of finite-sample uncertainty, with uncertainty in 
#' \code{beta_x} considered negligible. With an (asymptotically) valid confidence 
#' set for \code{delta_beta_y}, and under the "no measurement error" assumption, 
#' \code{budgetIV} returns an (asymptotically) valid confidence set for 
#' \eqn{\theta} when using just a single exposure.  
#' 
#' @return  
#' A \code{data.table} with each row corresponding to a set of bounds on the ATE 
#' at a given point in \code{ATE_search_domain}. Columns include: a non-unique 
#' identifier \code{curve_index} with a one-to-one mapping with \code{U}; 
#' \code{lower_ATE_bound} and \code{upper_ATE_bound} for the corresponding
#' bounds on the ATE; a list \code{U} for the corresponding budget assignment; 
#' and a column for each unique variable in \code{ATE_search_domain} to indicate 
#' the treatment value at which the bounds are being calculated.   
#' 
#' @references  
#' Jordan Penn, Lee Gunderson, Gecia Bravo-Hermsdorff, Ricardo Silva, and David 
#' Watson. (2024). BudgetIV: Optimal Partial Identification of Causal Effects 
#' with Mostly Invalid Instruments. \emph{AISTATS} 2025.
#' 
#' Jack Bowden, Fabiola Del Greco M, Cosetta Minelli, George Davey Smith, 
#' Nuala A Sheehan, and John R Thompson. (2016). Assessing the suitability of 
#' summary data for two-sample Mendelian randomization analyses using MR-Egger 
#' regression: the role of the I^2 statistic. \emph{Int. J. Epidemiol.} 46.6, 
#' pp. 1985--1998.
#' 
#' Charles F Manski. (1990). Nonparametric bounds on treatment effects. 
#' \emph{Am. Econ. Rev.} 80.2, pp. 219--323.
#' 
#' Henri Theil. (1953). Repeated least-squares applied to complete equation 
#' systems. \emph{Centraal Planbureau Memorandum}.
#' 
#' Eleanor Sanderson, George Davey Smith, Frank Windmeijer and Jack Bowden. (2019). 
#' An examination of multivariable Mendelian randomization in the single-sample and 
#' two-sample summary data settings. \emph{Int. J. Epidemiol.} 48.3, pp. 713--727.
#' 
#' @examples  
#' data(simulated_data_budgetIV)
#'
#' beta_y <- simulated_data_budgetIV$beta_y
#' 
#' beta_phi_1 <- simulated_data_budgetIV$beta_phi_1
#' beta_phi_2 <- simulated_data_budgetIV$beta_phi_2
#' 
#' beta_phi <- matrix(c(beta_phi_1, beta_phi_2), nrow = 2, byrow = TRUE)
#' 
#' delta_beta_y <- simulated_data_budgetIV$delta_beta_y
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
#' solution_set <- budgetIV(beta_y = beta_y, 
#'                          beta_phi = beta_phi, 
#'                          phi_basis = phi_basis, 
#'                          tau_vec = tau_vec, 
#'                          b_vec = b_vec, 
#'                          ATE_search_domain = ATE_search_domain, 
#'                          X_baseline = X_baseline,
#'                          delta_beta_y = delta_beta_y)
#' 
#' @export 
#' @import data.table
#' @import arrangements
#' @import MASS
#' @import Rglpk

budgetIV <- function(
    beta_y, 
    beta_phi, 
    phi_basis=NULL, 
    tau_vec=NULL, 
    b_vec=NULL, 
    ATE_search_domain=NULL, 
    X_baseline=NULL, 
    delta_beta_y=NULL 
) {
  
  # Convert numeric vector beta_y into matrix 
  if(is.vector(beta_y) && is.numeric(beta_y)){beta_y <- matrix(beta_y, nrow=1)}
  
  # Convert numeric vector beta_y into matrix
  if(is.vector(beta_phi) && is.numeric(beta_phi)){beta_phi <- matrix(beta_phi, nrow=1)}
  
  # Errors for incorrect format of beta_y or beta_phi
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
    stop("Argument 'beta_phi' must have more columns than rows. budgetIV only supports partial identification in the 'complete' in which the number of causal effect parameters 
         is no greater than the number of candidate instruments (d_{Phi} <= d_{Z}). See the package documentation or Penn et al. (2025) for further details.")
  }
  
  d_Z <- ncol(beta_y)
  d_Phi <- nrow(beta_phi)
  
  # Dealing with NULL tau_vec or b_vec 
  if(is.null(tau_vec)){
    if(is.null(b_vec)){
      tau_vec <- c(0)
      b_vec <- c(d_Z %/% 2)
    }
    else if(length(b_vec) == 1){
      tau_vec <- c(0)
    }
    else{
      stop("Argument 'tau_vec' is NULL while argument 'b_vec' has length greater than 1. When specifying multiple budgets, please specify as many thresholds.")
    }
  }
  else if(is.null(b_vec)){
    if(length(tau_vec) == 1){
      b_vec <- c(d_Z %/% 2)
    }
    else{
      stop("Argument 'b_vec' is NULL while argument 'tau_vec' has length greater than 1. When specifying multiple thresholds, please specify as many budgets.")
    }
  }

  # Dealing NULL phi_basis
  if(is.null(phi_basis)){
    # Null phi_basis and ATE_search_domain: define d_Phi-dimensional linear model.   
    if(is.null(ATE_search_domain)){
      warning(paste0("Arguments 'phi_basis' and 'ATE_search_domain' are not given: assuming treatment effect is linear in d_Phi = ", d_Phi, " variables."))
      if(d_Phi == 1){
        phi_basis <- expression(x)
        
        ATE_search_domain <- expand.grid("x" = seq(from = 0, to = 1, length.out = 50))
        
        if(is.null(X_baseline)){
          X_baseline <- list("x" = 0)
        }
      }
      else{
        X_names <- lapply(1:d_Phi, function(i) as.symbol(paste0("x_", i)))
        
        phi_basis <- do.call(expression, X_names)
        
        ATE_search_domain <- expand.grid( rep( seq(0, 1, length.out = 2), d_Phi) )
        colnames(ATE_search_domain) <- X_names
      }
    }
    # Otherwise, if compatible with ATE_search_domain, use a d_Phi-dimensional linear model.
    else if(d_Phi == ncol(ATE_search_domain)){
      X_names <- colnames(ATE_search_domain)
      phi_basis <- do.call(expression, X_names)
      
      if(is.null(X_baseline)){
        X_baseline = setNames(list(rep(0, length(X_names))), X_names)
      }
    }
    else{
      stop("Argument 'phi_basis' is not given and argument ncol('ATE_search_domain') != nrow('beta_phi'). Cannot fit the default linear model: please specify a choice of 'phi_basis'.")
    }
  }
  # Dealing with NULL ATE_search_domain but well-defined phi_basis.
  else if(is.null(ATE_search_domain)){
    X_names <- all.vars(phi_basis)
    d_X <- length(X_names)
    if(d_X == 1){
      ATE_search_domain <- expand.grid(seq(from = 0, to = 1, length.out = 50))
      colnames(ATE_search_domain) <- X_names
    }
    else{
      ATE_search_domain <- expand.grid( rep( seq(from = 0, to = 1, length.out = 2), d_X ) )
      colnames(ATE_search_domain) <- X_names
    }
  }
  # Dealing with NULL X_baseline but well defined search domain and phi_basis.
  else if(is.null(X_baseline)){
    X_names <- colnames(ATE_search_domain)
    
    X_baseline <- setNames(list(rep(0, length(X_names))), X_names)
  }
  
  d_X <- ncol(ATE_search_domain)
  
  if (is.null(delta_beta_y)){
    warning("No confidence bounds for agument 'beta_y' given: treating 'beta_y' as an oracle summary statistic.")
    delta_beta_y <- numeric(d_Z)
    }
  
  if(is.vector(delta_beta_y)){delta_beta_y <- matrix(delta_beta_y, nrow=1)}
    
  # Error messages

  if(!is.vector(tau_vec)){
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
      if (length(unique(sapply(ATE_search_domain, length))) == 1){ 
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
  
  # if( nrow(beta_phi) == 1){
  #   warning("Since 'nrow(beta_phi) = 1', consider using budgetIV_scalar. budgetIV_scalar can partially identify the scalar 
  #           causal effect parameter of interest with super-exponential improvement on time complexity.")
  # }
  
  d_Z <- ncol(beta_y)
  d_Phi <- nrow(beta_phi)
  
  # List to contain bounded ATE curves and the corresponding decision variable U
  partial_identification_ATE <- data.table(
    "curve_index" = numeric(),
    "x" = list(),
    "lower_ATE_bound" = numeric(),
    "upper_ATE_bound" = numeric(),
    "U" = list()
  )
  
  phi_zero <- lapply(phi_basis, function(e) eval(e, X_baseline[1, ]))
  
  # Calling budgetIV_scalar for one-dimensional Phi.  
  #
  # General idea: (1) call to get bounds on theta for each U; (2) get ATE bounds using theta bounds
  #
  if(nrow(beta_phi) == 1){
    theta_bounds <- budgetIV_scalar(beta_y = beta_y, 
                                    beta_phi = beta_phi, 
                                    tau_vec = tau_vec,
                                    b_vec = b_vec, 
                                    delta_beta_y = delta_beta_y,
                                    bounds_only = FALSE
                                    )
    
    # 
    # For each X in ATE_search_domain and for each theta interval in theta_bounds optimize theta * phi_basis(X) (upper and lower bounds) 
    
    for(interval_index in 1:nrow(theta_bounds)){ 
    
      for (coord_index in 1:nrow(ATE_search_domain)) {
        
        phi_curr_coord <- lapply(phi_basis, function(e) eval(e, ATE_search_domain[coord_index, ]))
        
        delta_phi <- unlist(phi_curr_coord) - unlist(phi_zero)
        
        ATE_theta_plus <- theta_bounds$lower_bound[interval_index]
        ATE_theta_minus <- theta_bounds$upper_bound[interval_index]
          
        lower_ATE_bound <- min(ATE_theta_minus, ATE_theta_plus)
        upper_ATE_bound <- max(ATE_theta_minus, ATE_theta_plus)
        
        new_bound <- data.table(
          "curve_index" = interval_index,
          "x" = list(ATE_search_domain[coord_index, ]),
          "lower_ATE_bound" = lower_ATE_bound,
          "upper_ATE_bound" = upper_ATE_bound,
          "U" = theta_bounds$budget_assignment[interval_index]
        )
        
        partial_identification_ATE <- rbind(partial_identification_ATE, new_bound)
        
      }
      
    }
    
    return(partial_identification_ATE)
    
    
  }
  
  dummy_infinity <- signif(max(beta_y, beta_phi) * d_Z * 1e10, 1)
  
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
  
  curve_index <- 0
  
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
    
    f.obj <- rep(0, d_Phi)  # Objective function for theta (length d_phi, the dimension of theta)
    
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
    
    partial_identification_ATE <- partial_identification_ATE[, cbind(.SD, rbindlist(partial_identification_ATE$x, use.names = TRUE, fill = TRUE)), .SDcols = !'x']
    
  }
  
  return(partial_identification_ATE)
  
  }
}
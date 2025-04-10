#' Efficient partial identification of a scalar causal effect parameter with invalid instruments
#' 
#' Partial identification and coverage of a causal effect parameter using summary statistics and budget constraint assumptions.
#' See Penn et al. (2025) for technical definitions.
#' 
#' @param beta_y A \eqn{d_{Z}}-dimensional vector representing the (estimated) 
#' cross covariance \eqn{\mathrm{Cov}(Y, Z)}.
#' @param beta_phi A \eqn{d_{Z}}-dimensional vector representing the (estimated) 
#' cross covariance \eqn{\mathrm{Cov}(\Phi (X), Z)}.
#' @param tau_vec A \eqn{K}-dimensional vector of increasing, positive 
#' thresholds representing degrees of IV invalidity. 
#' The default value \code{NULL} can be used for a single threshold at \eqn{0}.
#' @param b_vec A \eqn{K}-dimensional vector of increasing positive integers 
#' representing the maximum number of IVs that can surpass each threshold. 
#' The default value \code{NULL} can be used for a single threshold at \eqn{0}, with at least \eqn{50\%} of IVs assumed to be valid.
#' @param delta_beta_y A \eqn{d_{Z}}-dimensional vector of positive half-widths for box-shaped 
#' confidence bounds on \code{beta_y}. 
#' The default value \code{NULL} can be used to not include finite sample uncertainty.
#' @param bounds_only A boolean \code{TRUE} or \code{FALSE}. \code{TRUE} will store overlapping intervals in the confidence set as a single interval, 
#' while \code{FALSE} will store different intervals for different values of \code{budget_assignment} (see return value of Penn et al. (2025) for further details).
#' The default is \code{TRUE}.
#' 
#' If TRUE (default), the output consists only of disjoint bounds. Otherwise, if FALSE, the output consists of bounds for 
#' possibly touching intervals (but never overlapping), as well as the budget assignment corresponding to each bound. 
#' 
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: (A1) they are associated with the treatment; 
#' (A2) they are unconfounded with the outcome; and (A3) they exclusively effect the outcome through the treatment. 
#' Assumption (A1) has a simple statistical test, whereas for many data generating processes (A2) and (A3) are 
#' unprovably false. 
#' The \code{budgetIV} and \code{budgetIV_scalar} algorithms allow for valid causal inference when some proportion, 
#' possibly a small minority, of candidate instruments satisfy both (A2) and (A3).
#' 
#' \code{budgetIV} & \code{budgetIV_scalar} assume a homogeneous treatment effect, which implies the separable structural 
#' equation \eqn{Y = \theta \Phi(X) + g_y(Z, \epsilon_x)}. 
#' The difference between the algorithms is that \code{budgetIV_scalar} assumes \eqn{\Phi(X)} and \eqn{\theta} take
#' scalar values, which is exploited for super-exponential computational speedup and allows for causal inference
#' with thousands of candidate instruments \eqn{Z}.
#' Both methods assume ground truth knowledge of the functional form of \eqn{\Phi (X)}, e.g., a linear, 
#' logistic, Cox hazard, principal component based or other model. 
#' The parameter \eqn{\theta} captures the unknown treatment effect.
#' Violation of (A2) and/or (A3) will bias classical IV approaches through the statistical dependence
#' between \eqn{Z} and \eqn{g_y(Z, \epsilon_x)}, summarized by the covariance parameter 
#' \eqn{\gamma := \mathrm{Cov} (g_y(Z, \epsilon_x), Z)}.
#' 
#' \code{budgetIV} & \code{budgetIV_scalar} constrain \eqn{\gamma} through a series of positive thresholds 
#' \eqn{0 \leq \tau_1 < \tau_2 < \ldots < \tau_K} and corresponding integer budgets \eqn{0 < b_1 < b_2 < \ldots < b_K \leq d_Z}. 
#' It is assumed for each \eqn{i \in \{ 1, \ldots, K\}} that no more than \eqn{b_i} components of \eqn{\gamma} are greater in 
#' magnitude than \eqn{\tau_i}.
#' For instance, taking \eqn{d_Z = 100}, \eqn{K = 1}, \eqn{b_1 = 5} and \eqn{\tau_1 = 0} means 
#' assuming \eqn{5} of the \eqn{100} candidates are valid instrumental variables (in the sense that their ratio 
#' estimates \eqn{\theta_j := \mathrm{Cov}(Y, Z_j)/\mathrm{Cov}(\Phi(X), Z_j)} are unbiased).
#' 
#' With \code{delta_beta_y = NA}, \code{budgetIV} & \code{budgetIV_scalar} return the identified set
#' of causal effects that agree with both the budget constraints described above and the values of
#' \eqn{\mathrm{Cov}(Y, Z)} and \eqn{\mathrm{Cov}(Y, Z)}, assumed to be exactly precise. 
#' Unlike classical partial identification methods (see Manski (1990) ofr a canonical example), the non-convex mixed-integer
#' budget constraints yield a possibly disconnected identified set. 
#' Each connected subset has a different interpretation as to which of the candidate instruments \eqn{Z} 
#' are valid up to each threshold.
#' \code{budgetIV_scalar} returns these interpretations alongside the corresponding bounds on \eqn{\theta}. 
#' 
#' When \code{delta_beta_y} is not null, it is used as box-constraints to quantify uncertainty in \code{beta_y}. 
#' In the examples, \code{delta_beta_y} is calculated through a Bonferroni correction and gives an (asymptotically) 
#' valid confidence set over \code{beta_y}. 
#' Under the so-called "no measurement error" (NOME) assumption (see Bowden et al. (2016)) which is commonly applied in Mendelian randomisation, it is
#' assumed that the estimate of \code{beta_y} is the dominant source of finite-sample uncertainty, with uncertainty in \code{beta_x}
#' entirely negligible. 
#' With an (asymptotically) valid confidence set over \code{delta_beta_y} and under the "no measurement error" assumption, \code{budgetIV_scalar} 
#' returns an (asymptotically) valid confidence set for \eqn{\theta}.  
#' 
#' @return A data.table with each row corresponding to bounds on the scalar causal effect parameter \eqn{\theta} corresponding to a particular budget assignment \eqn{U} (see Penn et al. (2025)). 
#' The return table has the following rows: a logical \code{is_point} determining whether the upper and lower bounds are equivalent; numerical \code{lower_bound}
#' and \code{upper_bound} giving the lower and upper bounds; and a list \code{budget_assignment} giving the value of \eqn{U} for each candidate instrument. 
#' \code{budget_assignment} will only be returned if \code{bounds_only == FALSE} as input by the user.
#' 
#' A list of two entries: \code{intervals}, which is a two-column matrix with rows corresponding to disjoint bounds containing plausible values of \eqn{\theta}; 
#' and \code{points}, which is a one-column matrix consisting of lone plausible values of \eqn{\theta}---relevant when using \eqn{\tau_1 = 0}.   
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
#' @examples  
#' data(Do_et_al_summary_statistics)
#' 
#' candidatesHDL = Do_et_al_summary_statistics[Do_et_al_summary_statistics$pHDL <= 1e-8, ]
#' 
#' candidate_labels <- candidatesHDL$rsID
#' d_Z <- length(candidate_labels)
#' 
#' beta_x <- candidatesHDL$betaHDL
#' 
#' beta_y <- candidatesHDL$betaCAD
#' 
#' SE_beta_y <- abs(beta_y) / qnorm(1-candidatesHDL$pCAD/2)
#' 
#' alpha = 0.05
#' delta_beta_y <- qnorm(1 - alpha/(2*d_Z))*SE_beta_y
#' 
#' feasible_region <- budgetIV_scalar(
#'                                    beta_y = candidatesHDL$betaCAD,
#'                                    beta_phi = beta_x,
#'                                    tau_vec = c(0),
#'                                    b_vec = c(30),
#'                                    delta_beta_y = delta_beta_y,
#'                                    bounds_only = FALSE
#'                                    )
#' 
#' @export 
#' @import data.table
#' 
#' 

budgetIV_scalar <- function(
    beta_y,
    beta_phi,
    tau_vec=NULL,
    b_vec=NULL,
    delta_beta_y=NULL,
    bounds_only=TRUE
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
    if(nrow(delta_beta_y) == 1){
      delta_beta_y <- as.vector(delta_beta_y)
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
  
  # Dealing with NULL tau_vec and/or b_vec 
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


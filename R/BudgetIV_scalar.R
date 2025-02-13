#' BudgetIV for scalar exposures
#' 
#' Partial identification and coverage of a causal effect parameter using summary statistics and budget constraint assumptions.
#' 
#' @param beta_y A \eqn{1 \times d_{Z}} vector representing the (estimated) 
#' cross covariance \eqn{\mathrm{Cov}(Y, Z)}.
#' @param beta_phi A \eqn{1 \times d_{Z}} vector representing the (estimated) 
#' cross covariance \eqn{\mathrm{Cov}(\Phi (X), Z)}.
#' @param tau_vec, A \eqn{1 \times K} vector of strictly increasing, positive budget thresholds. \eqn{K} is the number of budget groups.
#' @param b_vec A \eqn{1 \times K} vector of increasing positive integer budgets. 
#' Represents the constraint that at least \eqn{b_i} different values of \eqn{j}, 
#' the candidate instrument \eqn{Z_j} satisfies \eqn{\mathrm{Cov} (Y - \theta \Phi (X), Z_j) \leq \tau_j}.
#' @param delta_beta_y Either (a) a \eqn{1 \times d_{Z}} vector of positive half-widths for box-shaped confidence bounds on \code{beta_y};
#' or (b) the empty value \code{NA} for partial identification or feasible region estimates without uncertainty quantification. 
#' 
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: (A1) they are associated with the treatment; 
#' (A2) they are unconfounded with the outcome; and (A3) they exclusively effect the outcome through the treatment. 
#' Assumption (A1) has a simple statistical test, whereas for many data generating processes (A2) and (A3) are 
#' unprovably false. 
#' The \code{BudgetIV} and \code{BudgetIV_scalar} algorithms allow for valid causal inference when some proportion, 
#' possibly a small minority, of candidate instruments satisfy both (A2) and (A3).
#' 
#' \code{BudgetIV} & \code{BudgetIV_scalar} assume a homogeneous treatment effect, which implies the separable structural 
#' equation \eqn{Y = \theta \Phi(X) + g_y(Z, \epsilon_x)}. 
#' The difference between the algorithms is that \code{BudgetIV_scalar} assumes \eqn{\Phi(X)} and \eqn{\theta} take
#' scalar values, which is exploited for super-exponential computational speedup and allows for causal inference
#' with thousands of candidate instruments \eqn{Z}.
#' Both methods assume ground truth knowledge of the functional form of \eqn{\Phi (X)}, e.g., a linear, 
#' logistic, Cox hazard, principal component based or other model. 
#' The parameter \eqn{\theta} captures the unknown treatment effect.
#' Violation of (A2) and/or (A3) will bias classical IV approaches through the statistical dependence
#' between \eqn{Z} and \eqn{g_y(Z, \epsilon_x)}, summarized by the covariance parameter 
#' \eqn{\gamma := \mathrm{Cov} (g_y(Z, \epsilon_x), Z)}.
#' 
#' \code{BudgetIV} & \code{BudgetIV_scalar} constrain \eqn{\gamma} through a series of positive thresholds 
#' \eqn{0 \leq \tau_1 < \tau_2 < \ldots < \tau_K} and corresponding integer budgets \eqn{0 < b_1 < b_2 < \ldots < b_K \leq d_Z}. 
#' It is assumed for each \eqn{i \in \{ 1, \ldots, K\}} that no more than \eqn{b_i} components of \eqn{\gamma} are greater in 
#' magnitude than \eqn{\tau_i}.
#' For instance, taking \eqn{d_Z = 100}, \eqn{K = 1}, \eqn{b_1 = 5} and \eqn{\tau_1 = 0} means 
#' assuming \eqn{5%} of the 100 candidates are valid instrumental variables (in the sense that their ratio 
#' estimates \eqn{\theta_j := \mathrm{Cov}(Y, Z_j)/\mathrm{Cov}(\Phi(X), Z_j) are unbiased).
#' 
#' With \code{delta_beta_y = NA}, \code{BudgetIV} & \code{BudgetIV_scalar} return the identified set
#' of causal effects that agree with both the budget constraints described above and the values of
#' \eqn{\mathrm{Cov}(Y, Z)} and \eqn{\mathrm{Cov}(Y, Z)}, assumed to be exactly precise. 
#' Unlike classical partial identification methods (see Manski (1990) ofr a canonical example), the non-convex mixed-integer
#' budget constraints yield a possibly disconnected identified set. 
#' Each connected subset has a different interpretation as to which of the candidate instruments \eqn{Z} 
#' are valid up to each threshold.
#' \code{BudgetIV_scalar} returns these interpretations alongside the corresponding bounds on \eqn{\theta}. 
#' 
#' When \code{delta_beta_y} is not null, it is used as box-constraints to quantify uncertainty in \code{beta_y}. 
#' In the examples, \code{delta_beta_y} is calculated through a Bonferroni correction and gives an (asymptotically) 
#' valid confidence set over \code{beta_y}. 
#' Under the so-called "no measurement error" (NOME) assumption (see Bowden et al. (2016)) which is commonly applied in Mendelian randomisation, it is
#' assumed that the estimate of \code{beta_y} is the dominant source of finite-sample uncertainty, with uncertainty in \code{beta_x}
#' entirely negligible. 
#' With an (asymptotically) valid confidence set over \code{delta_beta_y} and under the "no measurement error" assumption, \code{BudgetIV_scalar} 
#' returns an (asymptotically) valid confidence set for \eqn{\theta}.  
#' 
#' @return  
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
#' set.seed(42)
#' 
#' # Simple experiment to compare BudgetIV with convex background constraint methods.
#' 
#' 
#' 
#' @export 
#' @import data.table
#'
#' 

BudgetIV_scalar_exposure <- function(
    beta_y,
    beta_phi,
    tau_vec,
    b_vec,
    delta_beta_y=NA
) {
  
  #
  if (is.unsorted(tau_vec)) {
    stop('Please input tau constraints in increasing order.')
  }
  else if (any(duplicated(tau_vec))){ 
    stop('The same tau constraint cannot be specified twice.')
  }
  else if (is.unsorted(b_vec) || any(duplicated(b_vec))) {
    stop('The vector m must be strictly increasing, please see the definition of boldface m in the manuscript.')
  }
  else if (length(beta_y) != length(beta_phi)) {
    stop('beta_y and beta_phi must be vectors of the same length for scalar Phi(X). Please call "BudgetIV" for treatment of vector Phi(X).')
  }
  else if (any(is.na(delta_beta_y)) ){
    warning('No confidence bounds for beta_y given: resorting to a point estimate.')
    delta_beta_y <- numeric(d_Z)
  }
  else if (length(delta_beta_y != length(beta_y))){
    stop('delta_beta_y, if given, must be of the same length as beta_y.')
  }
  else if (!is.numeric(delta_beta_y) || any(delta_beta_y < 0)){
    stop('delta_beta_y, if given, must consist only of positive real numbers corresponding to uncertainty half-widths for each element of beta_y.')
  }
  
  d_Z <- length(beta_y)
  
  tau_intervals_lower <- matrix(nrow = (length(tau_vec)), ncol = d_Z)
  
  tau_intervals_upper <- matrix(nrow = (length(tau_vec)), ncol = d_Z)
  
  for (k in 1:length(tau_vec)){
    for (j in 1:d_Z){
      
      tau_intervals_lower[k, j] <- beta_y[j] / beta_phi[j] - (tau_vec[k]+delta_beta_y[j]) / beta_phi[j]
      
      tau_intervals_upper[k, j] <- beta_y[j] / beta_phi[j] + (tau_vec[k]+delta_beta_y[j]) / beta_phi[j]
      
    }}
  
  feasible_points <- matrix(nrow = 0, ncol = 1)
  
  feasible_intervals <- matrix(nrow = 0, ncol = 2)
  
  possible_bounds <- sort(c(c(tau_intervals_lower), c(tau_intervals_upper)))
  
  in_feasible <- FALSE
  
  for (p in 1:(length(possible_bounds)-1)){
    
    curr_point <- possible_bounds[p]
    
    curr_point_interval <- (possible_bounds[p] + possible_bounds[p+1])/2
    
    curr_interval_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point_interval, b_vec, tau_vec, delta_beta_y)
    
    curr_point_feasible <- validPoint_scalar(beta_y, beta_phi, d_Z, curr_point, b_vec, tau_vec, delta_beta_y)
    
    if(curr_interval_feasible && !in_feasible){
      
      last_feasible_opening = curr_point
      in_feasible <- TRUE
      
    }
    
    else if(!curr_interval_feasible && in_feasible){
      
      in_feasible <- FALSE
      feasible_intervals <- rbind(feasible_intervals, c(last_feasible_opening, curr_point))
      
    }
    
    else if(!in_feasible && curr_point_feasible){
      
      feasible_points <- rbind(feasible_points, c(curr_point))
      
    }
    
  }
  
  if(in_feasible){ feasible_intervals <- rbind(feasible_intervals, c(last_feasible_opening,  possible_bounds[length(possible_bounds)])) }
  
  return(list("points" = feasible_points, "intervals" = feasible_intervals))
  
}

validPoint_scalar <- function(beta_y, beta_phi, d_Z, theta, b_vec, tau_vec, delta_beta_y){
  
  b_to_fill <- b_vec
  
  for(j in 1:d_Z){
    
    beta_theta_j <- beta_phi[j] * theta
    
    tol <- .Machine$double.eps * max(beta_y[j], beta_theta_j)
    
    gamma_j <- beta_y[j] - beta_phi[j] * theta
    
    for(k in 1:length(tau_vec)){
      
      if(abs(gamma_j) <= tau_vec[k]+delta_beta_y[j] + tol){ b_to_fill[k] <- b_to_fill[k] - 1 }
      
    }
  }
  
  #print(m_to_fill)
  
  return(all(b_to_fill <= 0))
  
}
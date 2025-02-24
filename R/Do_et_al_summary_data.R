#' Summary statistics from Do et al.(2013)
#'
#' Common variants associated with plasma triglycerides and risk for coronary artery disease. 
#' Preprocessed and harmonized summary statistics from a Mendelian randomization analysis, including
#' summary statistics for variants' association with plasma triglyceride levels, serum HDL levels, serum 
#' LDL levels and risk of coronary artery disease (log-odds ratio).  
#' Dataset previously applied in the mode-based estimate approach of Hartwig et al. (2017).
#'
#' @docType data
#'
#' @usage data(Do_et_al_summary_statistics)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references Ron Do et al. (2013). 
#' Common variants associated with plasma triglycerides and risk for coronary artery disease. 
#' \emph{Nat Genet.} 45.11, pp. 1345--52.
#' 
#' Fernando Pires Hartwig, George Davey Smith, and Jack Bowden. (2017). 
#' Robust inference in summary data Mendelian randomization via the zero modal pleiotropy assumption.
#' \emph{Int. J. Epidemiol.} 46.6, pp. 1985--1998.
#' 
#'
#' @examples
#' 
#' # Extracting relevant summary statistics to investigate the causal effect of HDL on CAD risk.
#' 
#' data(Do_et_al_summary_statistics)
#' 
#' MBE_data <- read.csv('MBE_data.csv')
#'
#' candidatesHDL = MBE_data[MBE_data$pHDL <= 1e-8, ]
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
#' # For confidence set in BudgetIV/BudgetIV_scalar.
#' alpha = 0.05
#' delta_beta_y <- qnorm(1 - alpha/(2*d_Z))*SE_beta_y
#' 
"Do_et_al_summary_statistics"
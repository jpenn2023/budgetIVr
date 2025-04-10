#' Summary statistics from Do et al. (2013)
#'
#' Common variants associated with plasma triglycerides and risk for coronary artery disease. 
#' Preprocessed and harmonized summary statistics from a Mendelian randomization analysis, including
#' summary statistics for variants' association with plasma triglyceride levels, serum HDL levels, serum 
#' LDL levels and risk of coronary artery disease (CAD).  
#' Dataset previously applied in the mode-based estimate approach of Hartwig et al. (2017).
#' Each row of the dataset corresponds to a single genetic variant (single nucleotide polymorphism) found to be associated with either 
#' the HDL, LDL, or triglyceride biomarkers across a population of 180,000 (HDL, LDL) or 86,000 (triglyceride) individuals.
#' Got further biological and statistical details, see Do et al. (2013).
#' 
#' @docType data
#'
#' @usage data(Do_et_al_summary_statistics)
#'
#' @format A data frame with 185 rows and 14 variables:
#' 
#' @details
#' \describe{
#'  \item{\code{X}}{A unique identifier from 1 to 185.}
#'  \item{\code{rsID}}{A unique string specifying each SNP using the rsID format.}
#'  \item{\code{chr}}{String specifying the chromosomal position of each SNP.}
#'  \item{\code{a1}}{Character specifying one allele of the SNP (all 185 SNPs are assumed to be biallelic).}
#'  \item{\code{a2}}{Character specifying the other allele of the SNP.}
#'  \item{\code{betaLDL}}{Effect size (linear regression) for association between SNP allele and LDL.}
#'  \item{\code{pLDL}}{p-value for testing association between SNP allele and LDL.}
#'  \item{\code{betaHDL}}{Effect size (linear regression) for association between SNP allele and HDL.}
#'  \item{\code{pHDL}}{p-value for testing association between SNP allele and HDL.}
#'  \item{\code{betaTri}}{Effect size (linear regression) for association between SNP allele and triglyceride.}
#'  \item{\code{pTri}}{p-value for testing association between SNP allele and triglyceride.}
#'  \item{\code{betaCAD}}{Effect size (logistic regression) for association between SNP allele and CAD.}
#'  \item{\code{pCAD}}{p-value for testing association between SNP allele and CAD.}
#' }
#'  
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
#' @examples
#' # Extracting relevant summary statistics to investigate the causal effect of HDL on CAD risk.
#' 
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
#' # For confidence set in budgetIV/budgetIV_scalar.
#' alpha = 0.05
#' delta_beta_y <- qnorm(1 - alpha/(2*d_Z))*SE_beta_y
#' 
"Do_et_al_summary_statistics"
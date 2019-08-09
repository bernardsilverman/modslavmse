#' \pkg{modslavmse}
#'
#' \pkg{modslavmse} is a package for Multiple Systems Estimation as applied particularly in the context of Modern Slavery
#' 
#'This package serves three purposes.   Firstly, it implements the method used in the original work by Bales, Hesketh and Silverman
#'  to approach the question of the prevalence of Modern Slavery in the UK by estimating the 'dark figure' underlying data
#'  collected in the 2013 National Referral Mechanism (NRM) by the National Crime Agency (NCA).  That work used the package
#'  \pkg{Rcapture} with a particular approach to the choice of two-factor interactions to fit in the model. This method is implemented
#' in the routine \code{\link{MSEfit}}.  The data set \code{\link{UKdat}} gives the data used in the original paper.
#'   The research also involved testing the robustness of the results by omitting some of the lists and/or combining some into single lists.
#'  The routines \code{\link{omitlists}} and \code{\link{mergelists}} allow these operations to be done.
#'
#' The other purpose of this package is to implement development versions of current research directions.  A current focus is on Monte
#' Carlo Markov Chain approaches which allow some sort of model averaging, rather than the focus on a particular model implicit
#' in \code{MSEfit}.   This makes use of the package MCMCpack.
#'
#' The third, more specific, purpose is to allow full reproducibility of the work presented in Silverman (2018).   
#' This is done through the scripts given in the Examples section below.
#'
#' @examples
#'
#' data(UKdat, UKdat_5, UKdat_4, Ned, Ned_5, NewOrl, NewOrl_5, Kosovo) # the datasets used in the paper
#' make_AIC_stepwise_table1()  # Table 5
#' make_AIC_stepwise_table2()  # Table 6
#' make_allmodels_plots_script() # Figures 1, 2, 3 and 4
#' make_MCMCfit_tables_script() # Tables 7, 8, 9, 10, 11 and 12
#' make_MCMCeffects_table_script() #   Table 13
#' make_madyor_table_script()     #  Table 14 and Figures 5 and 6, as well as some numbers in the text
#' make_LCMCR_table_script()      # Table 15
#'
#' @references 
#'
#' K. B. Bales,  O. Hesketh, and B. W. Silverman (2015). Modern Slavery in the UK: How many victims? Significance 12 (3), 16-21.
#'
#' B. W. Silverman (2018). Model fitting in Multiple Systems Analysis for the quantification of Modern Slavery: Classical and Bayesian approaches.
#' 
#'   
"_PACKAGE"
# > [1] '_PACKAGE'

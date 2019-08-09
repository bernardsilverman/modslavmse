#' Markov Chain Monte Carlo Multiple Systems Estimation
#' 
#' Use the \pkg{MCMCpack} routine \code{MCMCpoisson} to fit a multiple systems estimation model. 
#'
#' @details This routine uses the routine \code{\link{MCMCpoisson()}} from the \pkg{MCMCpack} package. 
#'   It works by casting the 
#'   MSE problem into the form required by that routine and passing the MCMC parameters to that routine. 
#'  It calls the routine \code{\link{cleanuplists}}.
#'
#'
#' @param zdata an incidence matrix of the form supplied to \code{MSEfit}
#' @param priorprec If positive, the reciprocal of the prior variance of the two factor interaction parameters.  
#'  If zero, an improper prior over \eqn{(-\infty, \infty)}. If negative, then only main effects are fitted.  
#' @param mainprec The reciprocal of the prior variance of the main effect parameters. Only used if \code{priorprec > 0}.
#' @param sigthresh If \code{sigthresh > 0} the analysis is repeated dropping out any effects whose ratio of posterior mean 
#'  to standard deviation does not reach \code{sigthresh}
#' @param totalpopest If \code{totalpopest=T} then routine outputs the MCMC of the dark figure plus the original data, 
#'   in other words an estimate of the total population
#' @param runzerothresh If supplied, the result of a previous run with \code{sigthresh=0} and \code{totalpopest=F}
#' @param ... control parameters such as \code{mcmc}, \code{thin} and \code{seed} passed to \code{MCMCpoisson()}
#' @export
# 


MCMCfit <- function(zdata, priorprec = 1, mainprec = 1e-04, sigthresh = 2, totalpopest = F, 
    runzerothresh = NULL, ...) {
    # given an incidence matrix of the Rcapture form, carry out a MCMC using the
    # MCMCpoisson prior Here priorprec=0 corresponds to an improper prior over (-inf,
    # inf) If priorprec is positive it is the reciprocal of the prior variance of the two
    # factor interactions priorprec0 is the reciprocal of the prior variance of the main
    # effects in that case If priorprec is negative then two-factor interactions are not
    # fitted If sigthresh > 0, then the analysis is repeated dropping out any effects
    # whose posterior mean/SD does not exceed sigthresh If totalpopest=T then output the
    # MCMC of the dark figure plus the original data observed
    nd2 = dim(zdata)[2]
    # set up model for MCMC poisson with interactions, distinguishing between priorprec=0
    # and priorprec > 0 because of a bug in the MCMCpoisson program, if priorprec>0 all
    # the diagonal elements of B0 have to be nonzero, so mainprec is used.
    varcount = dimnames(zdata)[[2]][nd2]
    varnames = dimnames(zdata)[[2]][-nd2]
    nobserved = sum(zdata[, nd2])
    ninteractions = (nd2 - 1) * (nd2 - 2)/2
    form = paste(varcount, "~ (", paste(varnames, collapse = "+"), ")^2", collapse = "")
    # find a start vector
    start1 = log(c(nobserved * 5, t(zdata[, -nd2]) %*% zdata[, nd2]/(nobserved * 5)))
    # 
    zfull = cleanuplists(zdata, includezerocounts = T)
    BB0 = mainprec
    if (priorprec < 0) 
        form = makeform(zfull, NULL, F)
    if (priorprec > 0) {
        form = makeform(zfull, NULL, T)
        BB0 = diag(c(rep(mainprec, nd2), rep(priorprec, ninteractions)))
        start1 = c(start1, rep(0, ninteractions))
    }
    if (priorprec == 0) {
        BB0 = 0
        zcfd = removeemptyoverlaps(zfull)
        zfull = zcfd$data
        rempairs = zcfd$removed
        if (!is.null(rempairs)) 
            nrem = length(rempairs)/2 else nrem = 0
        start1 = c(start1, rep(0, ninteractions - nrem))
        # Modify formula to remove effects which are minus inf
        form = makeform(zfull, rempairs, T)
    }
    # run MCMC poisson, unless runzerothresh is supplied, in which case that is the
    # result of a previous run with sigthresh=0 and totalpopest=F so should be read in
    if (is.null(runzerothresh)) 
        zMCMC = MCMCpoisson(formula = form, data = zfull, B0 = BB0, beta.start = start1, 
            ...) else zMCMC = runzerothresh
    # if sigthresh> 0 then run again, dropping out any non-significant effects,
    # thresholding at sigthresh but if sigthresh is zero or no effects are below
    # threshold, there is no refinement
    zzs = summary(zMCMC)$statistics
    zzsn = dimnames(zzs)[[1]]
    zznotsig = (abs(zzs[, 1]/zzs[, 2]) < sigthresh)
    zznotsig[1:nd2] = F
    npars = length(zzsn)
    nrem = sum(zznotsig)
    if (nrem > 0) {
        start1 = start1[1:(npars - nrem)]
        if (is.matrix(BB0)) {
            BB0 = diag(diag(BB0)[1:(npars - nrem)])
            if ((npars - nrem) == nd2) 
                BB0 = 0
        }
        remeffs = paste(zzsn[zznotsig], collapse = "-")
        form = paste(form, remeffs, sep = "-")
        zMCMC = MCMCpoisson(formula = form, data = zfull, B0 = BB0, beta.start = start1, 
            ...)
    }
    if (totalpopest) 
        zMCMC[, 1] = nobserved + exp(zMCMC[, 1])
    return(zMCMC)
    
}

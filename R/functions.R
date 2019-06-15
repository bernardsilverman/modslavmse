#' Postprocess results from MSEfit
#'
#' choose whether to return confidence levels or the full fit
#'
#' @param zdat Original data set
#' @param zfit The full fit as returned by \code{\link[Rcapture]{closedpCI.t}}
#' @param mX The particular model being fitted
#' @param retfit If T return the entire fit; if F return confidence intervals and the point estimate
#'
#' @export
MSEretconf <-
function (zdat,zfit,mX,retfit) 
{
# choose whether to return confidence levels or the full fit
    if (retfit) return(zfit)
# construct vector of estimate and confidence levels
  xx = rep(NA, 5)
  names(xx)= c("2.5%", "10%","point estimate", "90%", "97.5%")
  xx[c(1,3,5)]= as.numeric(zfit$CI)[c(2,1,3)]
  zfit1 = closedpCI.t(zdat, dfreq = T, mX = mX, alpha=0.2)
  xx[c(2,4)]= as.numeric(zfit1$CI)[c(2,3)]
  return(xx)
  }
#' Fit MSE model using the \code{Rcapture} routine \code{closedpCI.t()}
#'
#' Implements the method of Bales, Hesketh and Silverman.  The routine calls the routine \code{closedpCI.t()} from the package \code{Rcapture} and returns
#'   an object which has the structure described in the documentation of that package. 
#'
#' @param zdat Input data frame.  If there are \eqn{k} lists the first \eqn{k} columns are 1/0 to denote presence or absence on the relevant list; the last column is the count
#' @param mainonly If \code{mainonly = T} then the stepwise method of BHS is used to choose a model; if \code{F} then only main effects are fitted with no interactions
#' @param  pthresh The significance level that a new parameter has to reach to be included
#' @param returnfit Controls whether the entire fit is returned or just confidence intervals and a point estimate.
#'
#' @export
#' 
MSEfit <-
function(zdat = UKdat, mainonly = F, pthresh = 0.05,returnfit = T) 
    {
require(Rcapture)
    # if returnfit = T then return the fitted model in the form provided by closedpCI.t
    #   otherwise return 95% and 80% confidence intervals and the point estimate
    #
    # if mainonly= F carry out the stepwise model fitting procedure as set out in Bales,
    # Hesketh and Silverman, using pthresh as the p-value that new effects have to
    # achieve before being included.  
    # if mainonly= T then fit main effects only
m = dim(zdat)[2] - 1
    # create and fit a formula with just the main effects
fmain = "~."
mXform = as.formula(fmain)
zfit = closedpCI.t(zdat, dfreq = T, mX = mXform)
if (mainonly) return(MSEretconf(zdat,zfit,mXform,returnfit))
    #
aic0 = zfit$fit$aic
    # set up a list of interactions for the fitting of interactive models
    # let ints be the matrix of all possible interactions
ints = NULL
for (i in (1:(m - 1))) for (j in ((i + 1):m)) {
    ints = cbind(ints, c(i, j))
    }
nints = m * (m - 1)/2
intsincluded = rep(F, nints)
# add interactions one at a time trying each possible interaction not yet included
for (jcycle in (1:nints)) {
    # add each possible new interaction and check the improvement in aic
    aicnew = rep(aic0, nints)
    for (j in (1:nints)) {
        if (!intsincluded[j]) {
            form1 = paste(fmain, "+c", ints[1, j], "*c", 
                          ints[2, j], sep = "")
            mXf=as.formula(form1)
            zf = closedpCI.t(zdat, dfreq = T, mX = mXf)
            aicnew[j] = zf$fit$aic
        }
    }
    # now test if any interaction has made a difference
    if (min(aicnew) >= aic0) 
        return(MSEretconf(zdat,zfit,mXform,returnfit))
    aic0 = min(aicnew)
    jintmax = min((1:nints)[aicnew == aic0])
    fmain1 = paste(fmain, "+c", ints[1, jintmax], "*c", 
                   ints[2, jintmax], sep = "")
    mXf=as.formula(fmain1)
    zf1 = closedpCI.t(zdat, dfreq = T, mX = mXf)
    pval = rev(summary(zf1$fit)$coefficients[, 4])[1]
    if (pval > pthresh) 
        return(MSEretconf(zdat,zfit,mXform,returnfit))
    zfit = zf1
    mXform=mXf
    fmain = fmain1
    intsincluded[jintmax] = T
}
return(MSEretconf(zdat,zfit,mXform,returnfit))
}
#' Script to return tables for the UK data
#'
#' @export
make_AIC_stepwise_table1 <-
function () 
{
    # construct all the tables for the UK data
    z = matrix(nrow=9, ncol=5)
    zmain = c(T,F,F)
    zpthresh = c(NA, 0.001,0.05)
    zn1 = rep(c("UK six lists", "UK five lists", "UK four lists"),3)
    zn2 = rep(c("Main effects only", "AIC threshold 0.1%", "AIC threshold 5%"), each=3)
    dimnames(z)[[1]] = paste(zn1,zn2,sep=": ")
    dimnames(z)[[2]]=c("2.5%", "10%","point estimate", "90%", "97.5%")
    #
    #  fill in the matrix z
    #
    for (j in (1:3)) {
        z[3*j-2,] = MSEfit(UKdat,zmain[j],zpthresh[j],returnfit=F)
        z[3*j-1,] = MSEfit(UKdat_5,zmain[j],zpthresh[j],returnfit=F)
        z[3*j,] = MSEfit(UKdat_4,zmain[j],zpthresh[j],returnfit=F)
    }
    z = round(z/1000,1)
    
    return(z)
}
#' Script to return tables for the non-UK data sets
#'
#' @export
make_AIC_stepwise_table2 <-
function () 
{
    # construct all the tables for the UK data
    z = matrix(nrow=15, ncol=5)
    zmain = c(T,F,F)
    zpthresh = c(NA, 0.001,0.05)
    zn1 = rep(c("Netherlands", "NL five lists","New Orleans", "NO five lists", "Kosovo"),3)
    zn2 = rep(c("Main effects only", "AIC threshold 0.1%", "AIC threshold 5%"), each=5)
    dimnames(z)[[1]] = paste(zn1,zn2,sep=": ")
    dimnames(z)[[2]]=c("2.5%", "10%","point estimate", "90%", "97.5%")
    #
    #  fill in the matrix z
    #
    for (j in (1:3)) {
        z[5*j-4,] = MSEfit(Ned,zmain[j],zpthresh[j],F)
        z[5*j-3,] = MSEfit(Ned_5,zmain[j],zpthresh[j],F)
        z[5*j-2,] = MSEfit(NewOrl,zmain[j],zpthresh[j],F)
        z[5*j-1,] = MSEfit(NewOrl_5,zmain[j],zpthresh[j],F)
        z[5*j,] = MSEfit(Kosovo,zmain[j],zpthresh[j],F)
    }
    z = round(z/1000,1)
    
    return(z)
}
#'
#' Consolidate counts for identical capture histories
#'
#'Consider an incidence matrix \code{mse} with \eqn{k+1} columns, where the first \eqn{k} columns are 1/0 denoting whether
#'    the relevant list is in the capture history, and the last column is the count of cases with that capture history.  This 
#'    routine finds rows with the same capture history and consolidates them into a single row whose count is the sum of counts of
#'    the relevant rows.  If \code{includezerocounts = T} then it also includes all the capture histories with zero counts; otherwise
#'    these are all removed.  
#' 
#' @param includezerocounts  If F then remove zero rows.  If T then include all possible capture histories including those with zero count, excluding the all-zero row corresponding to the dark figure.
#'@export
cleanuplists <- function(zmse, includezerocounts = F) {
    
    zmse = as.matrix(zmse)
    m = dim(zmse)[2] - 1
    # construct full capture history matrix
    zm = NULL
    for (j in (1:m)) zm = rbind(cbind(1, zm), cbind(0, zm))
    ztot = apply(zm, 1, sum)
    zm = zm[order(ztot), ]
    zm = cbind(zm[-1, ], 0)
    dimnames(zm)[[2]] = dimnames(zmse)[[2]]
    # find row of zm corresponding to each row of zmse and update count
    bcode = zmse[, -(m + 1)] %*% (2^(1:m))
    bc = zm[, -(m + 1)] %*% (2^(1:m))
    for (j in (1:dim(zmse)[1])) {
        ij = (1:dim(zm)[1])[bc == bcode[j]]
        zm[ij, m + 1] = zm[ij, m + 1] + zmse[j, m + 1]
    }
    # remove rows with zero counts if includezerocounts is F
    if (!includezerocounts) 
        zm = zm[zm[, m + 1] > 0, ]
    zmse = as.data.frame(zm)
    return(zmse)
}

#' Plots of BIC against estimates for all pairwise models
#'
#' Produce all the plots in Section 3.3 of Silverman (2018)
#'
#' @export
make_allmodels_plots_script <-
function()
{
plot(closedpMS.t(UKdat_5,dfreq=T,maxorder=2),omitOutliers=T)
plot(closedpMS.t(UKdat_5,dfreq=T,maxorder=2),omitOutliers=F)
plot(closedpMS.t(UKdat_4,dfreq=T,maxorder=2),omitOutliers=F)
plot(closedpMS.t(Ned_5,dfreq=T,maxorder=2),omitOutliers=F)
return(NULL)
}
#'  A function to make individual Dirichlet entries
#'
#' @export
make_LCMCR_table <-
function (zdat=UKdat, thin=100, seed=12345,burnin=10000,...) 
{
    require(LCMCR)
    zdat1=createLCMCR(zdat)
zcr= lcmCR(zdat1, tabular=T, seed=seed,verbose=F, ...)
zsim= lcmCR_PostSampl(zcr, thinning=thin, burnin=burnin,output=F) 
quantiles= quantile(zsim, probs=c(0.025,0.1,0.5,0.9,0.975))
return(round(quantiles/1000, 1))
}
#' A script to make the whole table
#'
#' @export
make_LCMCR_table_script <-
function (burnin0 = 100000) 
{
    zout= matrix(NA,nrow=8,ncol=5)
    dimnames(zout)=list( c("UK six lists", "UK five lists","UK four lists", "Netherlands", "Netherlands five lists","New Orleans",
                           "NewOrleans five lists","Kosovo"), c("2.5%","10%","50%","90%","97.5%"))
    zout[1,]=make_LCMCR_table(UKdat, burnin=burnin0)
    zout[2,]=make_LCMCR_table(UKdat_5, , burnin=burnin0)
    zout[3,]=make_LCMCR_table(UKdat_4, burnin=burnin0)
    zout[4,]=make_LCMCR_table(Ned, burnin=burnin0)
    zout[5,]=make_LCMCR_table(Ned_5, , burnin=burnin0)
    zout[6,]=make_LCMCR_table(NewOrl, burnin=burnin0)
    zout[7,]=make_LCMCR_table(NewOrl_5, burnin=burnin0)
    zout[8,]=make_LCMCR_table(Kosovo, burnin=burnin0)
    return(zout)  
}
#'  Graphical models method
#'
#'  Calling routine for the methods in the DGA package implementing
#'   the approach of Madigan and York (1997)
#'
#' @references D. Madigan and J. C. York (1997). Bayesian methods for estimation of the size of a closed population. Biometrika 84, 19-31.
#'
#' @export
make_madyor_table <-
function (zdat,maxpop=NA, plotpost=F, returnmaxmodprob=T, ...) 
{
    #
    # find the quantiles of the posterior distribution of the total number in the population
    #
    require(dga)
    Y =createarraydata(zdat)
    p = dim(zdat)[2]-1
    nobserved = sum(zdat[,p+1])
    if (is.na(maxpop)) maxpop=11*nobserved
   if (p==3) graphs = graphs3
   if (p==4) graphs = graphs4
   if (p==5) graphs = graphs5
if(p>5 )  stop("Too many lists")
#
   Nmiss = seq(from=0, to=maxpop-nobserved, length=1000)
   delta = 2^{-p}
   postmax= bma.cr(Y,Nmiss, delta,graphs, ...)
ytot=nobserved+Nmiss
   if(plotpost) plotPosteriorN(postmax,ytot)
   zz= apply(postmax, 2,sum)
   if (returnmaxmodprob)  {
       model_posterior_probs = apply(postmax, 1, sum)
       cat("Posterior probability of most probable model for dataset",deparse(substitute(zdat)), "and maximum size", maxpop, ":", round(max(model_posterior_probs),2), "\n")
       } 
           py=cumsum(zz)
           quantiles = approx(x=py,y=ytot, xout=c(0.025,0.1,0.5,0.9,0.975))
           quant=round(quantiles$y/1000,1)
           names(quant)= paste(100*quantiles$x, "%", sep="")
           return(quant)       
}
#'  Graphical Models script
#' 
#'  Script to make all the tables and figures in the section (5.1) 
#'  in Silverman (2018)
#'   on graphical models using the method of Madigan and York (1997)
#'
#' @references D. Madigan and J. C. York (1997). Bayesian methods for estimation of the size of a closed population. Biometrika 84, 19-31.
#'
#' @export
make_madyor_table_script <-
function () 
{
    zout= matrix(NA,nrow=6,ncol=5)
    dimnames(zout)=list( c("UK five lists max 40K","UK four lists", "Netherlands five lists","Netherlands five lists max 200K",
                           "NewOrleans five lists","Kosovo"), c("2.5%","10%","50%","90%","97.5%"))
    zout[1,]=make_madyor_table(UKdat_5,maxpop=40000, plotpost=T)
    zout[2,]=make_madyor_table(UKdat_4, plotpost=T)
    zout[3,]=make_madyor_table(Ned_5, plotpost=T)
    zout[4,]=make_madyor_table(Ned_5,maxpop=250000, plotpost=T)
    zout[5,]=make_madyor_table(NewOrl_5)
    zout[6,]=make_madyor_table(Kosovo)
    
    return(zout)
    
}
#' Construct a model formula for use by MCMCfit()
#'
#' Given the data frame \code{zdata} in \code{Rcapture} form, 
#'    creates a formula with the interactions given. 
#'
#' @details  
#' This routine is called by the Markov Chain Monte Carlo routine \code{\link{MCMCfit}}.  The routine \code{\link{MSEfit}} uses a formula with 
#'  numbered variables \eqn{c1, c2, \ldots } rather than using the names of the variables. 
#' 
#' @param zdata input data matrix
#' @param interactions a \eqn{2 \times m} matrix of two-factor interactions
#' @param includetwofac If \code{T} include all two factor \emph{excluding} those in \code{interactions}.  If \code{F} then include all main effects but only those two factor interactions in \code{interactions}
#' @export
makeform <- function(zdata, interactions = NULL, includetwofac = T) {
    # The matrix interactions is a 2 x m matrix If includetwofac=T then include all two
    # factor interactions but then exclude those in interactions If includetwofac=F then
    # include all main effects and then include those in interactions the value colch
    # will be used to include or exclude in the formula
    znames = dimnames(zdata)[[2]]
    nz = length(znames)
    varcount = znames[nz]
    varnames = znames[-nz]
    nv = nz - 1
    colch = "+"
    # first construct the full effects part of the formula
    form = paste(varnames, collapse = "+")
    if (includetwofac) {
        form = paste("(", form, ")^2", sep = "")
        colch = "-"
    }
    form = (paste(varcount, "~", form, sep = ""))
    # now modify if interactions is not null first create a matrix of names rather than
    # numbers
    if (!is.null(interactions)) {
        internames = varnames[interactions]
        dim(internames) = dim(as.matrix(interactions))
        intervec = apply(internames, 2, paste, collapse = ":")
        form = paste(form, paste(intervec, collapse = colch), sep = colch)
    }
    return(form)
    
}
#' Markov Chain Monte Carlo Multiple Systems Estimation
#' 
#' Use the \pkg{MCMCpack} routine \code{MCMCpoisson} to fit a multiple systems estimation model. 
#'
#' @details This routine uses the routine \code{\link{MCMCpoisson()}} from the \pkg{MCMCpack} package. It works by casting the 
#'   MSE problem into the form required by that routine and passing the MCMC parameters to that routine. It calls the routine \code{\link{cleanuplists}}.
#'
#'
#' @param zdata an incidence matrix of the form supplied to \code{MSEfit}
#' @param priorprec If positive, the reciprocal of the prior variance of the two factor interaction parameters.  If zero, an improper prior over \eqn{(-\infty, \infty)}. If negative, then only main effects are fitted.  
#' @param mainprec The reciprocal of the prior variance of the main effect parameters. Only used if \code{priorprec > 0}.
#' @param sigthresh If \code{sigthresh > 0} the analysis is repeated dropping out any effects whose ratio of posterior mean to standard deviation does not reach \code{sigthresh}
#' @param totalpopest If \code{totalpopest=T} then routine outputs the MCMC of the dark figure plus the original data, in other words an estimate of the total population
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
#' Print out all the tables for the method of Silverman (2018)
#'
#' This gives all the tables for the Bayesian approach set out in Section 4 of Silverman (2018)
#'
#' @export
make_MCMCfit_tables_script <-
function () 
{
    cat(" \n UK data six lists \n")
    print(make_MCMC_table(UKdat))
    cat(" \n UK data five lists \n")
    print(make_MCMC_table(UKdat_5))
    cat(" \n UK data four lists \n")
    print(make_MCMC_table(UKdat_4))
    cat(" \n Netherlands data \n")
    print(make_MCMC_table(Ned))
    # make_MCMC_table(Ned_5)  [Not included in paper]
    cat(" \n New Orleans data five lists \n")
    print(make_MCMC_table(NewOrl_5))
    cat(" \n Kosovo data \n")
    print(make_MCMC_table(Kosovo))
    invisible()
}
#' Table for Bayesian approach of Silverman (2018)
#'
#' Find quantiles of the posterior distribution fo various 
#'  parameters of the prior and a specified range of values for
#' the thresholding step.
#'
#' @export
make_MCMC_table <-
function (zdata=UKdat, sigthreshseq=c(0,2,5), ... ) 
{
    require(MCMCpack)
nsig=length(sigthreshseq)
priorprec1= rep( c(0,0.1,1,10), times=nsig)
sigthresh1 = rep(sigthreshseq, each=4)
zr = MCMCfit(zdata, priorprec= -1, totalpopest=T)
zout= c( -1, 0 , quantile( zr[,1], probs=c(0.025,0.1,0.5,0.9,0.975)))
for (j in (1:(4*nsig))) {
 # print(".")
  zr = MCMCfit(zdata, priorprec = priorprec1[j], sigthresh=sigthresh1[j], totalpopest=T, ...)
  zout = rbind(zout, c( priorprec1[j], sigthresh1[j],  quantile( zr[,1], probs=c(0.025,0.1,0.5,0.9,0.975))))
}
zout= zout[,-1]
zout[,-1] = round(zout[,-1]/1000,1)
dimnames(zout)[[1]] = c("Main", rep(c("Uniform", "Variance 10", "Variance 1", "Variance 0.1"),nsig))
dimnames(zout)[[2]][1]=  "Threshold"
return(zout)
}
#'  Find effects that survive the threshold
#'
#'  In the Bayesian thresholding method of Silverman (2018), find the effects
#' that survive thresholding at various threshold levels for given \code{priorprec}
#'
#' @export
MCMCeffects <- 
function (zdata=UKdat,priorprec=1, ...) 
{
nv= dim(zdata)[2]
zout= rep(NA,2)
z = summary(MCMCfit(zdata, priorprec=priorprec,sigthresh=0, ...))$statistics[-(1:nv),]
zn = dimnames(z)[[1]]
sigrat= abs(z[,1]/z[,2])
#
i1= (sigrat > 2) & (sigrat <= 5)
i2 = (sigrat > 5)
#
zout= c( paste(zn[i1], collapse= ","),paste(zn[i2], collapse= ",") )
names(zout)=  c("threshold 2 only", "survives at threshold 5")
return(zout)
}
#' Script for finding all the effects that survive thresholding for the various datasets
#'
#' @export
make_MCMCeffects_table_script <-
function (priorprec=1) 
{
    zout= matrix(NA,nrow=6,ncol=2)
    dimnames(zout)=list( c("UK six lists", "UK five lists","UK four lists", "Netherlands", "New Orleans five lists",
                           "Kosovo"), c("2 to 5","5"))
    zout[1,]=MCMCeffects(UKdat, priorprec)
    zout[2,]=MCMCeffects(UKdat_5, priorprec)
    zout[3,]=MCMCeffects(UKdat_4, priorprec)
    zout[4,]=MCMCeffects(Ned, priorprec)
    zout[5,]=MCMCeffects(NewOrl_5, priorprec)
    zout[6,]=MCMCeffects(Kosovo, priorprec)
    return(zout) 
}
#' Merge specified lists into a single list
#'
#' In order to merge distinct sets of lists, the routine has to be called for each set, noting that the column numbers
#'  will change after previous calls
#'
#' @param i,j,k,l  The lists to be merged.  Only i and j have to be given. 
#' @param rowclean If rowclean=T then remove all redundant rows (either by combining rows with the
#'      same capture pattern after the lists have been merged or by removing rows with zero counts)
#' @export
#'
mergelists <- function(zmse, i, j = NULL, k = NULL, l = NULL, rowclean = T, includezerocounts = F) {
    zmse = as.matrix(zmse)
    im = sort(c(i, j, k, l))
    nmerge = length(im)
    if (nmerge > 2) {
        for (j in (2:nmerge)) zmse = mergelists(zmse, im[1], im[j] + 2 - j, rowclean = F)
    } else {
        zz = zmse
        zn = dimnames(zz)[[2]]
        ii = im[1]
        jj = im[2]
        newname = paste(zn[ii], zn[jj], sep = "", collapse = "")
        dimnames(zz)[[2]][ii] = newname
        zz[, ii] = pmax(zz[, ii], zz[, jj])
        zmse = zz[, -jj]
    }
    if (rowclean) 
        zmse = cleanuplists(zmse, includezerocounts = includezerocounts)
    zmse = as.data.frame(zmse)
    return(zmse)
}

# > [1] '_PACKAGE'
#'
#'  Omit one or more of the lists in an incidence matrix
#'
#'  Construct the incidence matrix obtained by omitting one or more of the original lists.  This is useful for diagnostic purposes
#'   and robustness testing, and also if some of the lists are less reliable.
#'  
#' @param i,j,k,l  The lists to be omitted.  Only i has to be given. 
#' @param rowclean If rowclean=T then combine rows with the same capture pattern
#' @param includezerocounts  If F remove all rows with zero counts.  If T then include all possible observed capture combinations.
#' @export
#'
omitlists <- function(zmse, i, j = NULL, k = NULL, l = NULL, rowclean = T, includezerocounts = F) {
    zmse = as.matrix(zmse)
    zmse = zmse[, -c(i, j, k, l)]
    if (rowclean) 
        zmse = cleanuplists(zmse, includezerocounts = includezerocounts)
    zmse = as.data.frame(zmse)
    return(zmse)
}
#'
#' Remove pairs of lists with zero overlap
#' 
#' All possible capture histories are included in the output, except for those which contain a pair of lists for which
#'  there is no overlap at all, whether or not in combination with other lists.
#'  So for a pair to be omitted, there must be zero counts in every capture history including both of those lists. 
#'  
#' 
#' @return data The incidence matrix omitting all capture histories which contain a pair with zero overlap
#' @return removed  A matrix whose columns give the omitted pairs
#' @export
#'
removeemptyoverlaps <- function(zdata) {
    
    zfull = cleanuplists(zdata, includezerocounts = T)
    nd2 = dim(zfull)[2]
    # for each pair, remove if all zeroes below it
    rempairs = NULL
    for (i in (1:(nd2 - 2))) for (j in ((i + 1):(nd2 - 1))) {
        countbelow = sum(zfull[, i] * zfull[, j] * zfull[, nd2])
        if (countbelow == 0) {
            rempairs = cbind(rempairs, c(i, j))
            zfull = zfull[zfull[, i] * zfull[, j] == 0, ]
        }
    }
    # 
    return(list(data = zfull, removed = rempairs))
}

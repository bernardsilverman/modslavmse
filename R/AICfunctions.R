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
        z[3*j-2,] = MSEfit(UKdat,zmain[j],zpthresh[j],F)
        z[3*j-1,] = MSEfit(UKdat_5,zmain[j],zpthresh[j],F)
        z[3*j,] = MSEfit(UKdat_4,zmain[j],zpthresh[j],F)
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

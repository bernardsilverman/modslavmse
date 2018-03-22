#' Fit MSE model using the \code{Rcapture} routine \code{closedpCI.t()}
#'
#' Implements the method of Bales, Hesketh and Silverman.  The routine calls the routine \code{closedpCI.t()} from the package \code{Rcapture} and returns
#'   an object which has the structure described in the documentation of that package. 
#'
#' @param zdat Input data frame.  If there are \eqn{k} lists the first \eqn{k} columns are 1/0 to denote presence or absence on the relevant list; the last column is the count
#' @param mainonly If \code{mainonly = T} then the stepwise method of BHS is used to choose a model; if \code{F} then only main effects are fitted with no interactions
#' @param  pthresh The significance level that a new parameter has to reach to be included
#'
#' @export
#' 
MSEfit <-
function (zdat=UKdat, mainonly=F, pthresh=0.05) 
{
#
#  if mainonly= F carry out the stepwise model fitting procedure as set out in Bales, Hesketh and Silverman,
#    using pthresh as the p-value that new effects have to achieve before being included.
#
#  if mainonly= T then fit main effects only
#
#  let m be the number of variables and ints be the matrix of all possible interactions
 m = dim(zdat)[2] - 1
 
#
# create and fit a formula with just the main effects
#
 fmain =  "~."
 zfit= closedpCI.t(zdat, dfreq=T, mX=as.formula(fmain))
 if (mainonly) return(zfit)
#
#   let aic0 be aic of result
#   set up a list of interactions for the fitting of interactive models
#
aic0  = zfit$fit$aic
ints = NULL
 for ( i in (1:(m-1))) for (j in ((i+1):m)) ints= cbind(ints, c(i,j))
 nints = m*(m-1)/2
 intsincluded = rep(F, nints)
#
#  add interactions one at a time trying each possible interaction not yet included
#
for (jcycle in (1:nints)) {
#
#  add each possible new interaction and check the improvement in aic
#
aicnew = rep(aic0,nints)
for ( j in (1:nints) ) {
  if (!intsincluded[j]) {
     form1 = paste( fmain, "+c", ints[1,j], "*c", ints[2,j],sep="")
     zf = closedpCI.t(zdat, dfreq=T, mX = as.formula(form1))
     aicnew[j] = zf$fit$aic    
  }}
# now test if any interaction has made a difference
#
   if ( min(aicnew) >= aic0) return(zfit)
   aic0 = min(aicnew)
   jintmax = min((1:nints)[aicnew==aic0])
   fmain1= paste( fmain, "+c", ints[1,jintmax], "*c", ints[2,jintmax],sep="")
   zf1 = closedpCI.t(zdat, dfreq=T, mX=as.formula(fmain1))
#
# check if latest addition of parameter is significant
#
   pval = rev(summary(zf1$fit)$coefficients[, 4])[1]  
   if (pval > pthresh) return(zfit) 
   zfit = zf1
   fmain = fmain1
   intsincluded[jintmax] =T
}
return(zfit)
 
}

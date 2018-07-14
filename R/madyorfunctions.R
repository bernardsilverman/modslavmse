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

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
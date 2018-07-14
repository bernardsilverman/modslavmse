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

#' Construct a model formula for use by MCMCfit()
#'
#'  Given the data frame \code{zdata} in \code{Rcapture} form, 
#'    creates a formula with the interactions given. 
#' 
#' @param zdata input data matrix
#' @param interactions a \eqn{2 \times m} matrix of two-factor interactions
#' @param includetwofac If \code{T} include all two factor \emph{excluding} those in \code{interactions}.  If \code{F} then include all main effects but only those two factor interactions in \code{interactions}
#' @export
makeform <-
function (zdata,interactions=NULL,includetwofac=T) 
{
#
#  The matrix interactions is a 2 x m matrix
#
#  If includetwofac=T then include all two factor interactions but then 
#    exclude those in interactions
#
# If includetwofac=F then include all main effects and then include those
#    in interactions
# the value colch will be used to include or exclude in the formula
#
znames=dimnames(zdata)[[2]]
nz=length(znames)
varcount=znames[nz]
varnames=znames[-nz]
nv=nz-1
colch= "+"
#
# first construct the full effects part of the formula
#
form=paste(varnames,collapse="+")
if (includetwofac) {
	form=paste("(",form,")^2", sep="")
	colch="-"}
form=(paste(varcount,"~",form, sep=""))
#
#  now modify if interactions is not null
#    first create a matrix of names rather than numbers
#
if (!is.null(interactions)) {
 internames=varnames[interactions]
 dim(internames)=dim(as.matrix(interactions))
 intervec = apply(internames,2,paste,collapse=":")
 form = paste(form,paste(intervec, collapse=colch),sep=colch)} 
return(form)

}

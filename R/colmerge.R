#'
#'  Preprocess incidence matrix
#'
#'  This function preprocesses an incidence matrix of the form used in Rcapture.
#'  It allows you to merge two or more capture lists into one, or to omit one or more lists altogether.
#'  In addition it can produce a cleaned-up list where rows corresponding to a particular
#'      combination of lists are consolidated into a single row, and can produce a full dataset where the combinations with zero count
#'       are also included.
#'
#'@param zmse An incidence matrix with k+1 columns, where the first k columns are 1/0 giving the pattern of lists in that row, and the final column is the count
#'@param i    If i is a positive number, merge lists i and j.  If i is a vector of positive numbers, merge all the lists in i. If i is a negative number or a vector of negative numbers, remove the corresponding lists. 
#'@param j    If i is a positive number, merge lists i and j. Otherwise j is ignored.  Note that to merge separate combinations of lists, or to both merge and omit, it is necessary to call the routine separately for each combination, noting that previous calls will alter the numbering of the lists
#'@param rowclean  If rowclean=T then remove all redundant rows (either by combining rows with the same capture pattern or by removing rows with zero counts)
#'@param removezerocounts  If this parameter is T then all rows with zero counts will be removed.  If F then return the full incidence matrix including combinations where the count is zero.
#' @export
colmerge <- function (zmse, i=0, j, rowclean=T, removezerocounts=T) 
{
#   First test whether to merge or omit
#     
zmse=as.matrix(zmse)
if (i[1] > 0 ) { if(length(i)>1) {
        i=sort(i)
        nmerge= length(i)
        for ( j in (2:nmerge)) zmse = colmerge(zmse,i[1],  i[j] + 2 - j, rowclean=F)
                } else {
        zz= zmse
        zn = dimnames(zz)[[2]]
        newname= paste(zn[i], zn[j], sep="", collapse="")
        dimnames(zz)[[2]][i] = newname
        zz[,i] = pmax(zz[,i], zz[,j])
        zmse = zz[,-j]} }
if (i[1] < 0) { zmse=zmse[, i] } 
if (!rowclean) return(zmse)
#
#  if rowclean=T clean up by merging rows with same capture patterns
#
#  first create matrix of all possible capture patterns and sort these into increasing order of interaction
#
m=dim(zmse)[2]-1
zm = NULL
for (j in (1:m)) zm = rbind( cbind(1, zm), cbind(0,zm))
ztot=apply(zm, 1, sum)
zm= zm[order(ztot),]
zm= cbind(zm[-1,],0)
dimnames(zm)[[2]]=dimnames(zmse)[[2]]
#
#  for each row in the original zmse update the corresponding count in zm 
#
bcode =zmse[, -(m+1)]%*%(2^(1:m))
bc =zm[, -(m+1)]%*%(2^(1:m))
for ( j in (1:dim(zmse)[1])) {
    ij= (1:dim(zm)[1])[bc==bcode[j]]
    zm[ij,m+1] = zm[ij,m+1] + zmse[j,m+1] }
#
# remove zero rows counts
#
if (removezerocounts) zm = zm[zm[,m+1]>0,]
zmse=as.data.frame(zm)
return(zmse)
}

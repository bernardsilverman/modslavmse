#' 
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

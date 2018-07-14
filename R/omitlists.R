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

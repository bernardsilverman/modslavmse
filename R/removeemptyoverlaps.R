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

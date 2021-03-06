#'
#' Consolidate counts for identical capture histories
#'
#'Consider an incidence matrix \code{mse} with \eqn{k+1} columns, where the first \eqn{k} columns are 1/0 denoting whether
#'    the relevant list is in the capture history, and the last column is the count of cases with that capture history.  This 
#'    routine finds rows with the same capture history and consolidates them into a single row whose count is the sum of counts of
#'    the relevant rows.  If \code{includezerocounts = T} then it also includes all the capture histories with zero counts; otherwise
#'    these are all removed.  
#' 
#' @param includezerocounts  If F then remove zero rows.  If T then include all possible capture histories including those with zero count, excluding the all-zero row corresponding to the dark figure.
#'@export
cleanuplists <- function(zmse, includezerocounts = F) {
    
    zmse = as.matrix(zmse)
    m = dim(zmse)[2] - 1
    # construct full capture history matrix
    zm = NULL
    for (j in (1:m)) zm = rbind(cbind(1, zm), cbind(0, zm))
    ztot = apply(zm, 1, sum)
    zm = zm[order(ztot), ]
    zm = cbind(zm[-1, ], 0)
    dimnames(zm)[[2]] = dimnames(zmse)[[2]]
    # find row of zm corresponding to each row of zmse and update count
    bcode = zmse[, -(m + 1)] %*% (2^(1:m))
    bc = zm[, -(m + 1)] %*% (2^(1:m))
    for (j in (1:dim(zmse)[1])) {
        ij = (1:dim(zm)[1])[bc == bcode[j]]
        zm[ij, m + 1] = zm[ij, m + 1] + zmse[j, m + 1]
    }
    # remove rows with zero counts if includezerocounts is F
    if (!includezerocounts) 
        zm = zm[zm[, m + 1] > 0, ]
    zmse = as.data.frame(zm)
    return(zmse)
}

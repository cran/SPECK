#' Clustered thresholding of a vector.
#'
#' Performs thresholding for a vector of length \eqn{m} from a corresponding \eqn{m x n} reduced rank
#' reconstructed (RRR) matrix. Thresholding of a vector is only performed
#' if more than one cluster is identified using the \code{\link[Ckmeans.1d.dp:Ckmeans.1d.dp]{Ckmeans.1d.dp::Ckmeans.1d.dp()}} function
#' based on a one-dimensional dynamic programming clustering algorithm, which functions
#' by minimizing the sum of squares of within-cluster distances from an element to its associated cluster mean. If more than
#' one cluster is present, then the RRR output corresponding to the nonzero elements of the least-valued cluster, as identified by
#' the cluster mean, is set to zero. All other values in the least and higher-valued clusters are retained.
#'
#' @param rrr.vector Vector of length \eqn{m} from the corresponding \eqn{m x n}
#' RRR matrix.
#' @param max.num.clusters Maximum number of clusters for computation.
#' @param seed.ckmeans Seed specified to ensure reproducibility of the clustered thresholding.
#' @return
#' \itemize{
#'   \item rrr.thresholded.vector - A thresholded vector of length \eqn{m}.
#'   \item num.centers - Number of identified clusters.
#'   \item max.clust.prop - Proportion of samples with the specified maximum number of clusters.
#' }
#'
#' @examples
#' set.seed(10)
#' data.mat <- matrix(data = rbinom(n = 18400, size = 230, prob = 0.01), nrow = 80)
#' rrr.object <- randomizedRRR(counts.matrix = data.mat, rank.range.end = 60,
#' min.consec.diff = 0.01, rep.consec.diff = 2,
#' manual.rank = NULL, seed.rsvd = 1)
#' thresh.full.output <- ckmeansThreshold(rrr.vector = rrr.object$rrr.mat[,1],
#' max.num.clusters = 4, seed.ckmeans = 2)
#' head(thresh.full.output$rrr.thresholded.vector)
#' print(thresh.full.output$num.centers)
#' print(thresh.full.output$max.clust.prop)
#'
#' @export
ckmeansThreshold <- function(rrr.vector, max.num.clusters = 4, seed.ckmeans = 2) {
    if (missing(rrr.vector)) {
        stop("RRR vector is missing.")
    }
    set.seed(seed.ckmeans)
    ckmeans.res <- suppressWarnings(Ckmeans.1d.dp(x = rrr.vector, k = c(1:max.num.clusters)))
    prop.max.clust <- as.numeric(table(ckmeans.res$cluster)[max.num.clusters]/sum(table(ckmeans.res$cluster))*100)
    prop.max.clust[is.na(prop.max.clust)] <- 0
    if (length(ckmeans.res$centers) > 1) {
        min.center <- which(ckmeans.res$centers == min(ckmeans.res$centers))
        min.center.ind <- which(ckmeans.res$cluster == min.center)
        rrr.vector[min.center.ind] <- 0
    }
    threshold.list <- list(rrr.thresholded.vector = rrr.vector, num.centers = length(ckmeans.res$centers), max.clust.prop = prop.max.clust)
    return(threshold.list)
}

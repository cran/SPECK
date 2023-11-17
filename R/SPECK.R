#' Abundance estimation for single cell RNA-sequencing (scRNA-seq) data.
#'
#' Performs normalization, reduced rank reconstruction (RRR) and thresholding for a \eqn{m x n} scRNA-seq matrix
#' with \eqn{m} samples and \eqn{n} genes. The \code{\link[=speck]{speck()}} function calls the
#' \code{\link[=randomizedRRR]{randomizedRRR()}} function on the scRNA-seq matrix. Thresholding is next
#' applied to each gene from the \eqn{m x n} RRR matrix using the \code{\link[=ckmeansThreshold]{ckmeansThreshold()}}
#' function, resulting in a \eqn{m x n} thresholded matrix. See documentation for the \code{\link[=randomizedRRR]{randomizedRRR()}} and
#' \code{\link[=ckmeansThreshold]{ckmeansThreshold()}} functions for individual implementation details.
#'
#' @param counts.matrix \eqn{m x n} scRNA-seq counts matrix with \eqn{m} samples
#' and \eqn{n} genes.
#' @param rank.range.end Upper value of the rank for RRR.
#' @param min.consec.diff Minimum difference in the rate of change between a pair of successive standard deviation estimate.
#' @param rep.consec.diff Frequency of the minimum difference in the rate of change between a pair of successive standard deviation estimate.
#' @param manual.rank Optional, user-specified upper value of the rank used
#' for RRR as an alternative to automatically computed rank.
#' @param max.num.clusters Maximum number of clusters for computation.
#' @param seed.rsvd Seed specified to ensure reproducibility of the RRR.
#' @param seed.ckmeans Seed specified to ensure reproducibility of the clustered thresholding.
#' @return
#' \itemize{
#'   \item thresholded.mat - A \eqn{m x n} thresholded RRR matrix with \eqn{m} samples and \eqn{n} genes.
#'   \item rrr.mat - A \eqn{m x n} RRR matrix with \eqn{m} samples and \eqn{n} genes.
#'   \item rrr.rank - Automatically computed rank.
#'   \item component.stdev - A vector corresponding to standard deviations of non-centered sample principal components.
#'   \item clust.num - A vector of length \eqn{n} indicating the number of clusters identified by the
#'   \code{\link[=Ckmeans.1d.dp]{Ckmeans.1d.dp()}} algorithm for each gene.
#'   \item clust.max.prop - A vector of length \eqn{n} indicating the proportion of samples with the
#'   specified maximum number of clusters for each gene.
#' }
#'
#' @examples
#' set.seed(10)
#' data.mat <- matrix(data = rbinom(n = 18400, size = 230, prob = 0.01), nrow = 80)
#' speck.full <- speck(counts.matrix = data.mat, rank.range.end = 60,
#' min.consec.diff = 0.01, rep.consec.diff = 2,
#' manual.rank = NULL, max.num.clusters = 4,
#' seed.rsvd = 1, seed.ckmeans = 2)
#' print(speck.full$component.stdev)
#' print(speck.full$rrr.rank)
#' head(speck.full$clust.num); table(speck.full$clust.num)
#' head(speck.full$clust.max.prop); table(speck.full$clust.max.prop)
#' speck.output <- speck.full$thresholded.mat
#' dim(speck.output); str(speck.output)
#'
#' @export
speck <- function(counts.matrix, rank.range.end = 100, min.consec.diff = 0.01, rep.consec.diff = 2, manual.rank = NULL, max.num.clusters = 4, seed.rsvd = 1, seed.ckmeans = 2) {
    if (missing(counts.matrix)) {
        stop("Gene expression matrix is missing.")
    }
    message("Performing normalization and reduced rank reconstruction\n")
    counts.matrix.rrr <- randomizedRRR(counts.matrix = counts.matrix, rank.range.end = rank.range.end, min.consec.diff = min.consec.diff, rep.consec.diff = rep.consec.diff, manual.rank = manual.rank, seed.rsvd = seed.rsvd)
    message("Performing thresholding\n")
    rrr.thresholded.output <- apply(X = counts.matrix.rrr$rrr.mat, MARGIN = 2, FUN = ckmeansThreshold, max.num.clusters = max.num.clusters,
        seed.ckmeans = seed.ckmeans)
    message("Structuring results\n")
    rrr.thresholded.reformat <- as.data.frame(do.call(rbind, rrr.thresholded.output))
    data.rrr.thresholded <- as.matrix(as.data.frame(rrr.thresholded.reformat$rrr.thresholded.vector, check.names = FALSE))
    clusters.rrr.thresholded <- as.matrix(t(as.data.frame(rrr.thresholded.reformat$num.centers, check.names = FALSE)))
    colnames(clusters.rrr.thresholded) <- "number.clusters"
    max.clust.prop <- as.matrix(t(as.data.frame(rrr.thresholded.reformat$max.clust.prop, check.names = FALSE)))
    colnames(max.clust.prop) <- "proportion.max.clust"
    thresholded.rrr.list <- list(thresholded.mat = data.rrr.thresholded, rrr.mat = counts.matrix.rrr$rrr.mat, rrr.rank = counts.matrix.rrr$rrr.rank, component.stdev = counts.matrix.rrr$component.stdev, clust.num = clusters.rrr.thresholded, clust.max.prop = max.clust.prop)
    return(thresholded.rrr.list)
}

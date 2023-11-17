#' Reduced rank reconstruction (RRR) of a matrix.
#'
#' Computes the rank and subsequent RRR of a \eqn{m x n} counts matrix. Log-normalization is first performed
#' using the \code{\link[Seurat:NormalizeData]{Seurat::NormalizeData()}} function. RRR is next performed on the normalized \eqn{m x n}
#' matrix using randomized Singular Value Decomposition with the \code{\link[rsvd:rsvd]{rsvd::rsvd()}} function. Estimated rank is selected via a
#' construction of the standard deviations of non-centered sample principal components, which are used in a subsequent rate of change computation
#' where each successive standard deviation value is compared to the previous to determine the rank at which the absolute value of the rate of change
#' between consecutive values is at least 0.01 for at least two value pairs.
#'
#' @param counts.matrix A \eqn{m x n} counts matrix.
#' @param rank.range.end Upper value of the rank for RRR.
#' @param min.consec.diff Minimum difference in the rate of change between a pair of successive standard deviation estimate.
#' @param rep.consec.diff Frequency of the minimum difference in the rate of change between a pair of successive standard deviation estimate.
#' @param manual.rank Optional, user-specified upper value of the rank used
#' for RRR as an alternative to automatically computed rank.
#' @param seed.rsvd Seed specified to ensure reproducibility of the RRR.
#' @return
#' \itemize{
#'   \item rrr.mat - A \eqn{m x n} RRR matrix.
#'   \item rrr.rank - Automatically computed rank.
#'   \item component.stdev - A vector corresponding to standard deviations of non-centered sample principal components.
#' }
#'
#' @examples
#' set.seed(10)
#' data.mat <- matrix(data = rbinom(n = 18400, size = 230, prob = 0.01), nrow = 80)
#' rrr.object <- randomizedRRR(counts.matrix = data.mat, rank.range.end = 60,
#' min.consec.diff = 0.01, rep.consec.diff = 2,
#' manual.rank = NULL, seed.rsvd = 1)
#' print(rrr.object$component.stdev)
#' print(rrr.object$rrr.rank)
#' dim(rrr.object$rrr.mat); str(rrr.object$rrr.mat)
#' @export
randomizedRRR <- function(counts.matrix, rank.range.end = 100, min.consec.diff = 0.01, rep.consec.diff = 2, manual.rank = NULL, seed.rsvd = 1) {
  if (missing(counts.matrix)) {
    stop("Counts matrix is missing.")
  }
  if (rank.range.end > min(dim(counts.matrix))) {
    stop("High value of the manually specified rank.range.end
         parameter must be smaller than the minimum dimension of the
         counts.matrix input.")
  }
  if (!is.null(manual.rank) && (manual.rank > rank.range.end)) {
    stop("Manually specified rank must be less than the high value of the specified rank.range.end parameter.")
  }
  data.object <- CreateSeuratObject(counts = Matrix(t(as.matrix(counts.matrix))))
  data.object <- NormalizeData(object = data.object, normalization.method = "LogNormalize", verbose = FALSE)
  val <- 0; initial <- 0; final <- 0;
  set.seed(seed.rsvd)
  rsvd.res <- rsvd(t(as.matrix(data.object@assays$RNA@layers$data)), k = rank.range.end)
  if (!is.null(manual.rank)) {
    xhat <- rsvd.res$u[, 1:manual.rank] %*% diag(rsvd.res$d[1:manual.rank]) %*% t(rsvd.res$v[, 1:manual.rank])
    rownames(xhat) <- colnames(data.object)
    colnames(xhat) <- rownames(data.object)
    rrr.list <- list(rrr.mat = xhat, component.stdev = NULL, rrr.rank = manual.rank)
  } else {
    rsvd.sdev <- sqrt(rsvd.res$d^2/(nrow(counts.matrix) - 1))
    rsvd.diff <- round(diff(rsvd.sdev, differences = 1, lag = 1), 2)
    rle.amt <- rle(rsvd.diff)
    rle.amt$values <- abs(rle.amt$values)
    min.val <- min(rle.amt$values[which(rle.amt$values >= min.consec.diff & rle.amt$lengths >= rep.consec.diff)])
    consec.inst <- cumsum(c(1, diff(which(abs(rsvd.diff) == min.val)) != 1))
    ind.min <- min(consec.inst[duplicated(consec.inst)])
    automated.rank <- which(abs(rsvd.diff) == min.val)[ind.min+1]
    if (length(automated.rank) == 0) {
      stop("No stable automated rank value identified. Specify a higher value for the rank.range.end parameter.")
    }
    xhat <- rsvd.res$u[, 1:automated.rank] %*% diag(rsvd.res$d[1:automated.rank]) %*% t(rsvd.res$v[, 1:automated.rank])
    rownames(xhat) <- colnames(data.object)
    colnames(xhat) <- rownames(data.object)
    rrr.list <- list(rrr.mat = xhat, component.stdev = rsvd.sdev, rrr.rank = automated.rank)
  }
  return(rrr.list)
}


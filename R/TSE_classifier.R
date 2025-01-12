#' Classify Samples Based on TSE Categories
#'
#' This function classifies samples into TSE categories ("TSE_positive", "TSE_neutral", "TSE_negative") based on correlation with predefined centroids.
#'
#' @param x A numeric matrix of normalized gene expression data (rows: genes, columns: samples).
#' @param minCor Minimum correlation required to assign a TSE category. Default is 0.2.
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{TSE_category}}{Assigned TSE category for each sample.}
#'   \item{\code{cor_pval}}{p-value of the correlation with the nearest centroid.}
#'   \item{\code{separationLevel}}{A measure of how distinct the sample is from other categories.}
#'   \item{\code{TSE_positive, TSE_neutral, TSE_negative}}{Correlation scores for each TSE category.}
#'   \item{\code{sampleId}}{The ID of the sample.}
#' }
#' @examples
#' # Example usage with mock data
#' load(system.file("extdata", "centroids_TSE.RData", package = "TSECalculatoR"))
#' count_matrix <- matrix(rnorm(nrow(centroids_TSE * 10), mean = 1, sd = 1), nrow = nrow(centroids_TSE), ncol = 10)
#' rownames(count_matrix) <- rownames(centroids_TSE)
#' colnames(count_matrix) <- paste0("Sample", 1:10)
#' TSE_classify(count_matrix)
#' @export
TSE_classify <- function (x, minCor = 0.2)
{
  load(system.file("extdata", "centroids_TSE.RData", package = "TSECalculatoR"))
  TSEclasses <- c("TSE_positive", "TSE_neutral", "TSE_negative")
  genesToKeep <- intersect(rownames(centroids_TSE), rownames(x))
  if (length(genesToKeep) == 0) {
    stop("Genes provided are not in the list of genes used for classification.\n Make sure that gene names are symbols") }
  if (length(genesToKeep) < 0.6 * nrow(centroids_TSE))  {
    warning("Less than 60% of the genes used for classification are present in the data provided. Results may not be relevant") }
  cor.dat <- as.data.frame(cor(x[genesToKeep, ], centroids_TSE[match(genesToKeep, rownames(centroids_TSE)), TSEclasses], use = "complete.obs"), row.names = colnames(x))
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y) {
    TSEclasses[which.max(y)]
  })
  cor.dat$corToNearest <- apply(cor.dat[, TSEclasses], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp) {
    cor.test(x[genesToKeep, smp], centroids_TSE[match(genesToKeep, rownames(centroids_TSE)), cor.dat[smp, "nearestCentroid"]])$p.value
  })
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest -
                                        cor.dat[, TSEclasses], 1, function(x) {
                                          sort(x)[2]
                                        })
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, TSEclasses], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed
  cor.dat$TSE_category <- cor.dat$nearestCentroid
  try(cor.dat[which(cor.dat$corToNearest < minCor), "TSE_category"] <- NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <- NA)
  cor.dat <- cor.dat[, c("TSE_category", "cor_pval", "separationLevel", TSEclasses)]
  cor.dat$sampleId <- rownames(cor.dat)
  return(cor.dat)
}

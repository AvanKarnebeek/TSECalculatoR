#' Calculate TSE Score
#'
#' Computes the TSE score using GSVA on predefined gene signatures and categorizes the result into TSE categories.
#' Additionally, calculates a proliferation signature.
#'
#' @param count_matrix A numeric matrix of normalized gene expression data (rows: genes, columns: samples).
#' Row names should represent genes, and column names should represent samples.
#' @return A data frame containing:
#' \describe{
#'   \item{\code{sampleId}}{Sample IDs.}
#'   \item{\code{TSE_Category}}{The category of the TSE score: "Positive", "Neutral", or "Negative".}
#'   \item{\code{TSE_Value}}{The calculated TSE score (T-cell signature - stromal signature).}
#'   \item{\code{Tcellsig}}{The aggregated T-cell signature score.}
#'   \item{\code{Stromalsig}}{The aggregated stromal signature score.}
#'   \item{\code{Proliferation}}{The calculated proliferation signature score.}
#' }
#' @examples
#' # Example data
#' count_matrix <- matrix(runif(1000), nrow = 100)
#' rownames(count_matrix) <- paste0("Gene", 1:100)
#' colnames(count_matrix) <- paste0("Sample", 1:10)
#' calculate_tsescore(count_matrix)
#' @export

calculate_tsescore <- function(count_matrix) {
  # Validate the input
  if (!is.matrix(count_matrix)) stop("The input count_matrix must be a numeric matrix.")
  if (is.null(rownames(count_matrix))) stop("The count_matrix must have row names representing genes.")

  # Load predefined gene signatures from the package
  file_path <- system.file("extdata", "GeneSignaturesTSEscore.xlsx", package = "TSECalculatoR")

  gene_signatures <- readxl::read_excel(file_path)

  # Split gene signatures by their IDs
  genesTargetTherapy <- split(gene_signatures$Gene, gene_signatures$signatureID)

  # Validate gene presence
  required_genes <- unique(unlist(genesTargetTherapy))
  genes_present <- intersect(required_genes, rownames(count_matrix))

  # Check if no genes are present
  if (length(genes_present) == 0) {
    stop("None of the required genes are present in the count matrix. Please ensure gene names are symbols.")
  }

  # Warn if less than 60% of required genes are present
  percentage_present <- length(genes_present) / length(required_genes)
  if (percentage_present < 0.6) {
    warning("Less than 60% of the required genes are present in the count matrix. Results may not be reliable.")
  }


  # Check for empty gene signatures and remove them
  for (signature in names(genesTargetTherapy)) {
    genes_in_signature <- genesTargetTherapy[[signature]]
    present_genes <- intersect(genes_in_signature, rownames(count_matrix))
    if (length(present_genes) == 0) {
      message(signature, " signature doesn't have any genes and will be removed from the analysis. Results may not be reliable.")
      genesTargetTherapy[[signature]] <- NULL
    }
  }

  # Remove NULL entries
  genesTargetTherapy <- genesTargetTherapy[!sapply(genesTargetTherapy, is.null)]

  # Run GSVA to calculate scores
  gsvapar <- GSVA::gsvaParam(count_matrix, genesTargetTherapy)
  gsva_results <- GSVA::gsva(gsvapar)

  # Convert GSVA results to a data frame
  results_targetTherapyMarkers <- as.data.frame(gsva_results)
  results_targetTherapyMarkers$targetGroup <- rownames(results_targetTherapyMarkers)

  # Reshape the GSVA results: wide -> long -> wide
  results_targetTherapyMarkers <- reshape2::melt(results_targetTherapyMarkers, id.vars = "targetGroup")
  ResultsGeneSignatures <- reshape2::dcast(results_targetTherapyMarkers, variable ~ targetGroup, value.var = "value")

  # Dynamically calculate T-cell and stromal signatures based on available columns
  tcell_signatures <- intersect(c("IFN gamma", "T cell inflamed GEP", "T cell signature",
                                  "Immune gene signature", "Chemoattractants", "tGE8",
                                  "Cytotoxic CD8 T cell"), colnames(ResultsGeneSignatures))

  stromal_signatures <- intersect(c("Stromal signature", "Fibroblasts", "CAF", "TBRS", "EMT/stroma core genes"),
                                  colnames(ResultsGeneSignatures))

  ResultsGeneSignatures <- ResultsGeneSignatures %>%
    dplyr::mutate(
      Tcellsig = rowMeans(dplyr::select(., all_of(tcell_signatures)), na.rm = TRUE),
      Stromalsig = rowMeans(dplyr::select(., all_of(stromal_signatures)), na.rm = TRUE),
      TSE_Value = Tcellsig - Stromalsig,
      sampleId = variable
    ) %>%
    dplyr::select(sampleId, TSE_Value, Tcellsig, Stromalsig)

  return(ResultsGeneSignatures)
}

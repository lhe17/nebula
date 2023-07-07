#' Retrieve data from Seurat or Single Cell Experiment object to prepare for use in Nebula
#'
#' @param obj \code{Seurat} or \code{SingleCellExperiment} object to create data set for Nebula.
#' @param assay Assay to retrieve counts from the corresponding \code{Seurat} count matrix.
#' @param sid Sample_id to use metadata object i.e. \code{obj$sid}.
#' @param pred Character vector of predictors from metadata in \code{Seurat} or \code{SingleCellExperiment} objects.
#' @param offset Metadata column corresponding to per-cell scaling factor e.g. TMM.
#' @return data_neb: A list usable for Nebula.
#' @export
#' @examples
#' library(Seurat)
#' library(nebula)
#'
#' data("sample_seurat")
#' re <- scToNeb(obj = seu_obj, assay = "RNA", sid = "replicate", pred = c("celltype", "tech"))
#'


scToNeb <- function(obj, assay = NULL, sid = NULL, pred = NULL, offset = NULL, verbose = TRUE)
{
  p_df <- list()
  if ("SingleCellExperiment" %in% class(obj)) {
    counts <- counts(obj)
    for (k in pred){
      p_df[[k]] <- obj[[k]]
    }
    p_df <- data.frame(p_df, row.names = colnames(obj))
    data_neb <- list(counts = counts, pred = p_df)
    if (!is.null(sid)) {
      data_neb$sid <- as.character(obj[[sid]])
    }
    if (!is.null(offset)) {
      data_neb$offset <- as.double(obj[[offset]])
    }
  } else if ("Seurat" %in% class(obj)) {
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(obj)
      if (verbose)
        {Biobase::note("No assay provided. Using the default assay:", assay)}
    }
    if (is.null(sid)) {
      warning("No sample id provided, please add it prior to running nebula.")
    }
    Seurat::DefaultAssay(obj) <- assay
    for (k in pred){
      p_df[k] <- obj[[k]]
    }
    p_df <- data.frame(p_df, row.names = colnames(obj))
    data_neb <- list(
      counts = Seurat::GetAssayData(obj, slot = "counts"),
      pred = p_df
    )
    if (!is.null(sid)) {
      data_neb$sid <- as.character(obj[[sid]][[sid]])
    }
    if (!is.null(offset)) {
      data_neb$offset <- as.double(obj[[offset]][[offset]])
    }
  } else {stop("Please provide either a SingleCellExperiment or Seurat object")}
  return(data_neb)
}
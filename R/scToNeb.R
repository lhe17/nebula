#' Retrieve data from Seurat or SingleCellExperiment object to prepare for use in nebula
#'
#' @param obj \code{Seurat} or \code{SingleCellExperiment} object to create data set for Nebula.
#' @param assay Assay to retrieve counts from the corresponding \code{Seurat} count matrix.
#' @param id Sample ID to use metadata object i.e. \code{obj$id}.
#' @param pred Character vector of predictors from metadata in \code{Seurat} or \code{SingleCellExperiment} objects.
#' @param offset Metadata column corresponding to per-cell scaling factor e.g. TMM.
#' @param verbose Indicating whether to print additional messages.
#' @return data_neb: A list usable for nebula.
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(nebula)
#'
#' data("sample_seurat")
#' re <- scToNeb(obj = sample_seurat, assay = "RNA", id = "replicate", pred = c("celltype", "tech"))
#' }


scToNeb <- function(obj, assay = NULL, id = NULL, pred = NULL, offset = NULL, verbose = TRUE)
{
  p_df <- list()
  if ("SingleCellExperiment" %in% class(obj)) {
    covs = colnames(SingleCellExperiment::colData(obj))
    count <- SingleCellExperiment::counts(obj)
    for (k in pred){
      if(k %in% covs)
      {p_df[[k]] <- obj[[k]]}else{
        stop("A specified variable is not available in the meta data of the SingleCellExperiment object.")
      }
    }
    p_df <- data.frame(p_df, row.names = colnames(obj))
    data_neb <- list(count = count, pred = p_df)
    if (!is.null(id)) {
      if(id %in% covs)
      {data_neb$id <- as.character(obj[[id]])}else{
        stop("The sample ID is not available in the meta data of the SingleCellExperiment object.")
      }
    }else{
      warning("No sample ID provided, please specify it when running nebula.")
    }
    if (!is.null(offset)) {
      if(offset %in% covs)
      {data_neb$offset <- as.double(obj[[offset]])}else{
        stop("The offset is not available in the meta data of the SingleCellExperiment object.")
      }
    }
  } else if ("Seurat" %in% class(obj)) {
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(obj)
      if (verbose)
        {cat("No assay provided. Using the default assay:", assay)}
    }
    if (is.null(id)) {
      warning("No sample ID provided, please specify it when running nebula.")
    }
    Seurat::DefaultAssay(obj) <- assay
    covs = colnames(obj@meta.data)
    for (k in pred){
      if(k %in% covs)
      {p_df[k] <- obj[[k]]}else{
        stop("A specified variable is not available in the meta data of the Seurat object.")
      }
    }
    p_df <- data.frame(p_df, row.names = colnames(obj))
    data_neb <- list(
      count = Seurat::GetAssayData(obj, slot = "counts"),
      pred = p_df
    )
    if (!is.null(id)) {
      if(id %in% covs)
      {data_neb$id <- as.character(obj[[id]][[id]])}else{
        stop("The sample ID is not available in the meta data of the Seurat object.")
      }
    }
    if (!is.null(offset)) {
      if(offset %in% covs)
      {data_neb$offset <- as.double(obj[[offset]][[offset]])}else{
        stop("The offset is not available in the meta data of the Seurat object.")
      }
    }
  } else {stop("Please provide either a SingleCellExperiment or Seurat object")}
  return(data_neb)
}